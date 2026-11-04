"""Four confirmations of the root cause found by root_cause.py.

  1 open loop     projection and relaxed estimator at very large rho on one recorded log, window by window, without feedback chaos
  2 acquisition   which test blocks the first k_good in the scenes that never acquire it, against the tuned projection
  3 all cells     first gain supplied from truth for the full candidate and the projection, seeds of ridge_final.json
  4 tuning        30 gate evaluations for the candidate, the projection's budget; same selection seeds, then the same evaluation seeds
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import common
import physics_plant as PP
from dynamic_gain_value import Truth, stats
from failure_replay import Measured, cell_scene
from ridge_estimator import RidgeProjection
from ridge_final import Candidate, BASE as FINAL_BASE
from root_cause import Full, Plain, BASE as RC_BASE, estimator

HERE = Path(__file__).resolve().parent
LN12 = math.log(1.2)
FULL = dict(nw=24, order=4, fc=14., memory=True, sigma_scale=.796, chi_scale=.551)


def open_loop(controls):
    sc = cell_scene('breath', 'mid', 0, FINAL_BASE); seed = FINAL_BASE
    p = dict(S.FINAL_C); p.update(controls)
    pol = Plain(estimator('P3_16'), single_model=False, **p)
    log = PP.episode(sc, Truth(pol), seed)
    rel = RidgeProjection(nw=16, order=4, fc=1., memory=False)
    rel.reset(pol)

    class O:
        pass
    k, s, c = [], [], []
    for n in range(40, len(log)):
        f = common.window_fit(log, n, 16, 3, pol.params['noise'], pol.ns, pol.eb, pol.afc, pol.ndel, pol.params['enc_y'])
        o = O(); o.log = log; o.n = n
        a = rel.estimate(o, n, pol)
        if np.isfinite(f['k']) and np.isfinite(a[0]) and f['k'] > 0 and a[0] > 0:
            k.append(math.log(a[0]/f['k'])); s.append(a[1]*math.sqrt(f['info'])); c.append(a[2]/max(f['chi2'], 1e-12))
    return dict(windows=len(k), ln_k_ratio_p99=float(np.percentile(np.abs(k), 99)), sigma_ratio=[float(np.min(s)), float(np.max(s))], chi_ratio_p50=float(np.median(c)), chi_ratio_p1=float(np.percentile(c, 1)), chi_ratio_p99=float(np.percentile(c, 99)))


def acquisition(a):
    name, c, i, controls = a
    sc = cell_scene(c, 'hard', i, RC_BASE); seed = RC_BASE+sum(ord(x) for x in c+'hard')*100+i
    p = dict(S.FINAL_C); p.update(controls)
    if name == 'Full':
        pol = Full(RidgeProjection(**FULL), single_model=False, **p); pol.sc = sc; pol.centre = True
    else:
        pol = Plain(estimator('P5_16'), single_model=False, **p)
    w = Truth(pol); tr = []

    class Trace:
        def reset(self, platform):
            w.reset(platform)

        def __call__(self, o):
            u = w(o)
            chi = pol.fit.get('chi2', np.nan) if pol.fit else np.nan
            dof = getattr(pol.estimator, 'dof', 9)
            rel = pol.rel if hasattr(pol, 'rel') else np.nan
            tr.append((o.log[o.n-1, 0], rel, chi, dof, pol.k_seek, pol.k_good, float(o.log[o.n-1, 5]), pol.k_hat))
            return u
    log = PP.episode(sc, Trace(), seed)
    x = np.array(tr, float); t = x[:, 0]
    fin = np.isfinite(x[:, 5]); first = t[np.argmax(fin)] if fin.any() else -1.
    win = (t < (first if first > 0 else 1.0)) & np.isfinite(x[:, 1]) & np.isfinite(x[:, 2])
    if win.sum() == 0:
        return dict(name=name, cond=c, i=i, first=float(first), windows=0)
    rel, chi, dof = x[win, 1], x[win, 2], x[win, 3]
    precise = rel <= p['rel_max']; consistent = chi <= p['chi2_max']
    bound = rel*(1+np.sqrt(dof*np.maximum(0., chi-1.))); within = bound <= LN12
    kt = log[:len(tr), 8][win]; kh = x[win, 7]
    return dict(name=name, cond=c, i=i, first=float(first), windows=int(win.sum()), precise=float(precise.mean()), consistent=float(consistent.mean()), within=float(within.mean()),
                all_three=float((precise & consistent & within).mean()), seek_end=float(x[win, 4][-1]), rel_min=float(np.nanmin(rel)), chi_med=float(np.nanmedian(chi)), bound_min=float(np.nanmin(bound)),
                dof_med=float(np.median(dof)), ln_khat=float(np.nanmedian(np.abs(np.log(kh/np.maximum(kt, 1.))))), F_max=float(np.nanmax(x[win, 6])))


def full_cell(a):
    name, c, t, i, controls, spec, gift = a
    sc = cell_scene(c, t, i, FINAL_BASE); seed = FINAL_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(controls)
    if name.startswith('P5'):
        pol = Plain(estimator('P5_16'), single_model=False, **p)
    else:
        pol = Full(RidgeProjection(**spec), single_model=False, **p); pol.sc = sc; pol.centre = True
    pol.intervene = ('gift',) if gift else ()
    m = S._metric(PP.episode(sc, Truth(pol), seed), pol, sc, 20.8)
    return dict(name=name, cond=c, tissue=t, i=i, metrics=m)


def tune_one(a):
    spec, controls, ci, c, t, i = a
    sc = cell_scene(c, t, i, 340000); seed = 340000+ci*100+i
    p = dict(S.FINAL_C); p.update(controls)
    pol = Full(RidgeProjection(**spec), single_model=False, **p); pol.sc = sc; pol.centre = True
    return math.log1p(min(S._metric(PP.episode(sc, Truth(pol), seed), pol, sc, 20.8)['J'], 1000.))


def report(rows, names, ref):
    L = lambda m: math.log1p(min(m['J'], 1000)); key = lambda r: (r['cond'], r['tissue'], r['i'])
    base = {key(r): r['metrics'] for r in rows if r['name'] == ref}
    boot = lambda d: np.percentile(d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1), [2.5, 97.5])
    groups = (('soft', lambda r: r['tissue'] == 'soft'), ('mid', lambda r: r['tissue'] == 'mid'), ('hard', lambda r: r['tissue'] == 'hard'),
              ('event', lambda r: r['cond'] in ('flinch', 'cough')), ('fast', lambda r: r['cond'] == 'shift_fast'), ('smooth', lambda r: r['cond'] in ('breath', 'shift_slow')),
              ('rib+combo', lambda r: r['cond'] in ('rib', 'combo')))
    out = {}
    for n in names:
        rr = [r for r in rows if r['name'] == n]; s = stats(rr)
        d = np.array([L(r['metrics'])-L(base[key(r)]) for r in rr]); c = boot(d)
        g = {}
        for gn, f in groups:
            q = [r for r in rr if f(r)]; dg = np.array([L(r['metrics'])-L(base[key(r)]) for r in q]); cg = boot(dg)
            g[gn] = dict(d=float(dg.mean()), ci=[float(cg[0]), float(cg[1])], rmse=stats(q)['rmse_med'], fail=stats(q)['fail'])
        out[n] = dict(s, vs_ref=[float(d.mean()), float(c[0]), float(c[1])], groups=g)
        print(f"{n:18s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.3f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs {ref} {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]")
        print('    '+' '.join(f"{gn}:{v['d']:+.2f}[{v['ci'][0]:+.2f},{v['ci'][1]:+.2f}]r{v['rmse']:.3f}f{v['fail']}" for gn, v in g.items()))
    return out


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']
    res = {}
    res['open_loop'] = open_loop(controls); print('open loop', res['open_loop'], flush=True)
    rc = json.loads((HERE/'results'/'root_cause.json').read_text(encoding='utf-8'))['rows']
    never = [(r['cond'], r['i']) for r in rc if r['name'] in ('Full', 'R16_limit', 'R24+mem') and r['cause'] and r['cause']['kind'] == 'never_acquired']
    never = sorted(set(never))
    print('never acquired scenes', never, flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        acq = list(ex.map(acquisition, [(n, c, i, controls) for n in ('Full', 'P5_16') for c, i in never], chunksize=1))
        res['acquisition'] = acq
        for q in acq:
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in q.items()})
        rows = list(ex.map(full_cell, [(n, c, ti, i, controls, FULL, g) for n, g in (('P5_16', False), ('P5_16 gift', True), ('Full', False), ('Full gift', True)) for c, ti in S.CELLS for i in range(10)], chunksize=1))
        res['all_cells'] = report(rows, ['P5_16', 'P5_16 gift', 'Full', 'Full gift'], 'P5_16')
        rng = np.random.default_rng(340000)
        cands = [dict(FULL)]+[dict(FULL, sigma_scale=float(np.exp(rng.uniform(np.log(.2), np.log(5.)))), chi_scale=float(np.exp(rng.uniform(np.log(.2), np.log(5.))))) for _ in range(29)]
        scored = []
        for spec in cands:
            scored.append((float(np.mean(list(ex.map(tune_one, [(spec, controls, ci, c, ti, i) for ci, (c, ti) in enumerate(S.CELLS) for i in range(3)], chunksize=1)))), spec))
        v, best = min(scored, key=lambda x: x[0])
        print('tuned gate', round(v, 3), {k: round(x, 3) for k, x in best.items() if 'scale' in k}, 'score range', round(min(s for s, _ in scored), 3), round(float(np.median([s for s, _ in scored])), 3), flush=True)
        res['tuning'] = dict(best=best, scores=[s for s, _ in scored])
        rows2 = rows+list(ex.map(full_cell, [('Full tuned', c, ti, i, controls, best, False) for c, ti in S.CELLS for i in range(10)], chunksize=1))
        res['tuned'] = report(rows2, ['Full tuned'], 'P5_16')
    (HERE/'results'/'root_cause_confirm.json').write_text(json.dumps(res, indent=1, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
