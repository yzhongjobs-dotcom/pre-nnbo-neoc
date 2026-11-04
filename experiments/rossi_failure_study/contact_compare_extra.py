"""Three follow-ups to contact_compare.py.

fair hc     Measured force in the regressor was the harmful part of gamma, and the Hunt-Crossley form carried it. Repeat that form without it,
            on the matrix seeds, and compare with the stored elastic and weak-power rows.
minimal     y memory, z memory and centring were each removable while removing all three was not: which single one is enough?
lock        The hard-tissue scenes where the first gain was never accepted, without the acquisition rule, after the memory discretisation fix.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP
import contact_compare as C
from contact_compare import ModelPolicy, RIDGE, CTL, L, key, boot, scene, MATRIX_BASE, ABLATE_BASE, N_MATRIX, N_ABLATE, PLANTS
from dynamic_gain_value import Truth, stats
from failure_replay import cell_scene
from ridge_estimator import RidgeProjection, POP_ALPHA, POP_TAU
from root_cause import BASE as RC_BASE, HARD, N_SCENE as RC_N, classify

HERE = Path(__file__).resolve().parent
FAIR = {'hc rate only': dict(lam=.25), 'hc rate + conv gamma': dict(lam=.25, conv=.206), 'weak power + conv gamma': dict(memory=True, conv=.206)}
MINIMAL = {'relaxed + centring': dict(centre=True), 'relaxed + y memory': dict(ymem=True), 'relaxed + z memory': dict(zmem=True)}


def fair_job(a):
    plant, model, c, t, i = a
    f = FAIR[model]
    sc = scene(c, t, i, MATRIX_BASE, plant); seed = MATRIX_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    mem = f.get('memory', False)
    est = RidgeProjection(memory=mem, lam=f.get('lam', 0.), alpha=POP_ALPHA if mem else None, tau=POP_TAU if mem else None, **RIDGE)
    pol = ModelPolicy(est, single_model=False, **p)
    pol.lam = f.get('lam', 0.); pol.gamma = f.get('conv', 0.)
    if mem:
        pol.zalpha = POP_ALPHA; pol.ztau = POP_TAU
    return dict(plant=plant, model=model, cond=c, tissue=t, i=i, metrics=S._metric(PP.episode(sc, Truth(pol), seed), pol, sc, 20.8))


def minimal_job(a):
    name, c, t, i = a
    f = MINIMAL[name]
    sc = cell_scene(c, t, i, ABLATE_BASE); seed = ABLATE_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    pol = ModelPolicy(RidgeProjection(memory=f.get('ymem', False), **RIDGE), single_model=False, **p)
    if f.get('zmem'):
        pol.zalpha = POP_ALPHA; pol.ztau = POP_TAU
    pol.centre = f.get('centre', False); pol.acquire = True
    return dict(name=name, cond=c, tissue=t, i=i, metrics=S._metric(PP.episode(sc, Truth(pol), seed), pol, sc, 20.8))


def lock_job(a):
    c, t, i = a
    sc = cell_scene(c, t, i, RC_BASE); seed = RC_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    pol = ModelPolicy(RidgeProjection(memory=True, **RIDGE), single_model=False, **p)
    pol.zalpha = POP_ALPHA; pol.ztau = POP_TAU; pol.acquire = False
    w = Truth(pol); tr = []

    class Trace:
        def reset(self, platform):
            w.reset(platform)

        def __call__(self, o):
            u = w(o); tr.append((pol.kbar, pol.k_good)); return u
    log = PP.episode(sc, Trace(), seed)
    return dict(cond=c, i=i, metrics=S._metric(log, pol, sc, 20.8), cause=classify(log, np.array(tr, float), sc) if log is not None else dict(kind='numeric'),
                acquired=bool(np.isfinite(pol.k_good)))


def main():
    stored = json.loads((HERE/'results'/'contact_matrix.json').read_text(encoding='utf-8'))['rows']
    ablation = json.loads((HERE/'results'/'ablation.json').read_text(encoding='utf-8'))['rows']
    with ProcessPoolExecutor(max_workers=15) as ex:
        lock = list(ex.map(lock_job, [(c, t, i) for c, t in HARD for i in range(RC_N)], chunksize=1))
        kinds = {}
        for r in lock:
            if r['cause']:
                kinds[r['cause']['kind']] = kinds.get(r['cause']['kind'], 0)+1
        s = stats(lock)
        print(f"lock check, hard tissue seeds {RC_BASE}, no acquisition rule: L {s['L']:.3f} fail {s['fail']} long {s['long_loss']} causes {kinds} never had k_good at the end {sum(not r['acquired'] for r in lock)}", flush=True)
        mini = list(ex.map(minimal_job, [(n, c, t, i) for n in MINIMAL for c, t in S.CELLS for i in range(N_ABLATE)], chunksize=1))
        full = {key(r): r['metrics'] for r in ablation if r['name'] == 'Full'}; base = {key(r): r['metrics'] for r in ablation if r['name'] == 'P5_16 only'}
        rel = {key(r): r['metrics'] for r in ablation if r['name'] == 'relaxed prior only'}
        for n in MINIMAL:
            rr = [r for r in mini if r['name'] == n]; s = stats(rr)
            d = np.array([L(r['metrics'])-L(full[key(r)]) for r in rr]); c = boot(d)
            e = np.array([L(r['metrics'])-L(base[key(r)]) for r in rr]); f = boot(e)
            print(f"{n:20s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.2f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs Full {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] vs P5 {e.mean():+.3f} [{f[0]:+.3f},{f[1]:+.3f}]", flush=True)
        for n, b in (('Full', full), ('relaxed prior only', rel)):
            rr = [r for r in ablation if r['name'] == n]; s = stats(rr)
            e = np.array([L(r['metrics'])-L(base[key(r)]) for r in rr]); f = boot(e)
            print(f"{n:20s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.2f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs P5 {e.mean():+.3f} [{f[0]:+.3f},{f[1]:+.3f}]", flush=True)
        fair = list(ex.map(fair_job, [(pl, m, c, t, i) for pl in PLANTS for m in FAIR for c, t in S.CELLS for i in range(N_MATRIX)], chunksize=1))
    mean = {}
    for pl in PLANTS:
        el = {key(r): r['metrics'] for r in stored if r['plant'] == pl and r['model'] == 'elastic'}
        wp = {key(r): r['metrics'] for r in stored if r['plant'] == pl and r['model'] == 'weak power'}
        print(f"[{pl}] elastic L {stats([r for r in stored if r['plant'] == pl and r['model'] == 'elastic'])['L']:.3f}, weak power L {stats([r for r in stored if r['plant'] == pl and r['model'] == 'weak power'])['L']:.3f}")
        for m in FAIR:
            rr = [r for r in fair if r['plant'] == pl and r['model'] == m]; s = stats(rr); mean.setdefault(m, []).append(s['L'])
            d = np.array([L(r['metrics'])-L(el[key(r)]) for r in rr]); c = boot(d)
            g = np.array([L(r['metrics'])-L(wp[key(r)]) for r in rr]); h = boot(g)
            print(f"  {m:24s} L {s['L']:.3f} rmse {s['rmse_med']:.3f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs elastic {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] vs weak power {g.mean():+.3f} [{h[0]:+.3f},{h[1]:+.3f}]", flush=True)
    print('mean L over plants', {m: round(float(np.mean(v)), 3) for m, v in mean.items()})
    (HERE/'results'/'contact_compare_extra.json').write_text(json.dumps(dict(lock=lock, minimal=mini, fair=fair), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
