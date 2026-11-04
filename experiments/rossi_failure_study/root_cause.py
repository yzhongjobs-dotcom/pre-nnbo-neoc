"""Why the relaxed projection loses on hard tissue: separate implementation, setting, integration and the method itself.

Hard tissue only, same scenes and noise for every variant (seeds 403000, not used before). Each comparison answers one question:
  implementation   relaxed estimator at a very large rho against the cubic projection of the same window: must nearly coincide
  estimator form   16 and 24 ticks, projection and relaxed, all with the same gate scales
  integration      the first gain supplied from truth (removes acquisition); true k_good or true k_bar throughout (removes one signal)
Every runaway is classified from the trace: never acquired, after contact loss, fast event, slow divergence.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP
from dynamic_gain_value import Truth, stats
from estimator_swap import SwapPolicy
from failure_replay import cell_scene
from ridge_estimator import RidgeProjection
from ridge_final import Candidate

HERE = Path(__file__).resolve().parent
BASE = 403000
N_SCENE = 15
HARD = [(c, 'hard') for c in ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor')]
GATE = dict(sigma_scale=.796, chi_scale=.551)


class Mixin:
    intervene = ()

    def _stiffness_est(self, o, n):
        k8 = float(o.log[n-1, 8]) if n >= 1 else float('nan')
        if 'gift' in self.intervene and self.first and np.isfinite(k8) and k8 > 1 and o.log[n-1, 5] > .5:
            self.first = False; self.k_good = k8; self.kbar = k8*1.2; self.f_set = float(o.log[n-1, 5])
        out = super()._stiffness_est(o, n)
        if 'true_good' in self.intervene and np.isfinite(self.k_good) and np.isfinite(k8) and k8 > 1:
            self.k_good = k8
        return out


class Plain(Mixin, SwapPolicy):
    pass


class Full(Mixin, Candidate):
    pass


def estimator(kind):
    if kind == 'P3_16':
        return S.estimator('Projection', dict(nw=16, degree=3, sigma_scale=1., chi_scale=1.))
    if kind == 'R16_limit':
        return RidgeProjection(nw=16, order=4, fc=1., memory=False, sigma_scale=1., chi_scale=1.)
    if kind == 'P5_16':
        return S.estimator('Projection', dict(nw=16, degree=5, **GATE))
    if kind == 'P5_24':
        return S.estimator('Projection', dict(nw=24, degree=5, **GATE))
    if kind == 'P3_24':
        return S.estimator('Projection', dict(nw=24, degree=3, **GATE))
    if kind == 'R16':
        return RidgeProjection(nw=16, order=4, fc=14., memory=False, **GATE)
    if kind == 'R24':
        return RidgeProjection(nw=24, order=4, fc=14., memory=False, **GATE)
    if kind == 'R24+mem':
        return RidgeProjection(nw=24, order=4, fc=14., memory=True, **GATE)
    raise ValueError(kind)


VARIANTS = {
    'P3_16': ('P3_16', ()), 'R16_limit': ('R16_limit', ()),
    'P5_16': ('P5_16', ()), 'P5_24': ('P5_24', ()), 'P3_24': ('P3_24', ()), 'R16': ('R16', ()), 'R24': ('R24', ()), 'R24+mem': ('R24+mem', ()),
    'P5_16 gift': ('P5_16', ('gift',)), 'R24 gift': ('R24', ('gift',)),
    'P5_16 true_good': ('P5_16', ('true_good',)), 'R24 true_good': ('R24', ('true_good',)),
    'P5_16 true_bar': ('P5_16', ('true_bar',)), 'R24 true_bar': ('R24', ('true_bar',)),
    'Full': ('R24+mem', ('full',)), 'Full gift': ('R24+mem', ('full', 'gift')),
}


def classify(log, trace, sc):
    F = log[:len(trace), 3]; t = log[:len(trace), 0]; g = trace[:, 1]
    if F.max() <= 50:
        return None
    pk = int(np.argmax(F > 50)); on = pk
    while on > 0 and F[on] > 6:
        on -= 1
    rise = float(t[pk]-t[on])
    if not np.isfinite(g[:pk+1]).any():
        return dict(kind='never_acquired', t_on=float(t[on]), rise=rise)
    if not np.isfinite(g[on]):
        return dict(kind='acquired_during_runaway', t_on=float(t[on]), rise=rise)
    if (t[on] > .5) and (F[(t > t[on]-.5) & (t <= t[on])] < .5).any():
        return dict(kind='after_contact_loss', t_on=float(t[on]), rise=rise)
    kt = log[on, 8]
    return dict(kind='fast_event' if rise <= .15 else 'slow_divergence', t_on=float(t[on]), rise=rise,
                ln_good=float(np.log(g[on]/max(kt, 1.))), ln_bar=float(np.log(trace[on, 0]/max(kt, 1.))))


def job(a):
    name, c, t, i, controls = a
    kind, flags = VARIANTS[name]
    sc = cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(controls)
    if 'true_bar' in flags:
        p['oracle_fb'] = True
    if 'full' in flags:
        pol = Full(estimator(kind), single_model=False, **p); pol.sc = sc; pol.centre = True
    else:
        pol = Plain(estimator(kind), single_model=False, **p)
    pol.intervene = flags
    w = Truth(pol); tr = []

    class Trace:
        def reset(self, platform):
            w.reset(platform)

        def __call__(self, o):
            u = w(o); tr.append((pol.kbar, pol.k_good)); return u
    log = PP.episode(sc, Trace(), seed)
    m = S._metric(log, pol, sc, 20.8)
    out = dict(name=name, cond=c, tissue=t, i=i, metrics=m, cause=None)
    if log is None:
        out['cause'] = dict(kind='numeric', t_on=None, rise=None)
    else:
        out['cause'] = classify(log, np.array(tr, float), sc)
        kt = log[:len(tr), 8]; g = np.array(tr, float)[:, 1]; ok = (log[:len(tr), 0] > 1.) & (kt > 0) & np.isfinite(g) & (log[:len(tr), 3] < 20.)
        out['good_err'] = [float(np.percentile(np.log(g[ok]/kt[ok]), q)) for q in (10, 50, 90)] if ok.sum() > 20 else None
    return out


def equivalence(controls):
    """Commands of the projection and the relaxed estimator at very large rho, one scene and the same noise."""
    sc = cell_scene('shift_fast', 'mid', 2, BASE); seed = BASE+2
    p = dict(S.FINAL_C); p.update(controls)
    logs = []
    for kind in ('P3_16', 'R16_limit'):
        pol = Plain(estimator(kind), single_model=False, **p)
        logs.append(PP.episode(sc, Truth(pol), seed))
    n = min(len(logs[0]), len(logs[1]))
    return float(np.max(np.abs(logs[0][:n, 11]-logs[1][:n, 11]))), float(np.std(logs[0][:n, 11]))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']
    d, s = equivalence(controls)
    print(f'equivalence: max command difference {d*1e6:.3f} um, command std {s*1e6:.1f} um')
    args = [(n, c, ti, i, controls) for n in VARIANTS for c, ti in HARD for i in range(N_SCENE)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(job, args, chunksize=1))
    L = lambda m: math.log1p(min(m['J'], 1000)); key = lambda r: (r['cond'], r['i'])
    ref = {key(r): r['metrics'] for r in rows if r['name'] == 'P5_16'}
    boot = lambda x: np.percentile(x[np.random.default_rng(1).integers(0, len(x), (10000, len(x)))].mean(1), [2.5, 97.5])
    summary = {}
    for n in VARIANTS:
        rr = [r for r in rows if r['name'] == n]; st = stats(rr)
        dl = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(dl)
        kinds = {}
        for r in rr:
            if r['cause']:
                kinds[r['cause']['kind']] = kinds.get(r['cause']['kind'], 0)+1
        ge = [r['good_err'] for r in rr if r.get('good_err')]
        med = np.median(np.array(ge), axis=0) if ge else [np.nan]*3
        summary[n] = dict(st, vs_p5=[float(dl.mean()), float(c[0]), float(c[1])], causes=kinds, good_err=[float(x) for x in med])
        print(f"{n:16s} L {st['L']:.3f} rmse {st['rmse_med']:.3f}/{st['rmse_p90']:.2f} fail {st['fail']:2d} long {st['long_loss']:2d} | vs P5 {dl.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] | k_good P10/50/90 {med[0]:+.2f} {med[1]:+.2f} {med[2]:+.2f} | {kinds}")
    (HERE/'results'/'root_cause.json').write_text(json.dumps(dict(equivalence=[d, s], summary=summary, rows=rows), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
