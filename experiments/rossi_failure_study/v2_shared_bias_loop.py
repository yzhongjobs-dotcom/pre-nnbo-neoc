"""Controlled closed loops for the bias shared by short and long history: remove one factor, let the controller run, compare the same rib crossings.

The replay counterfactual (v2_shared_bias_source.py) is confounded: commands that tracked breathing remain, so removing breathing moves the depth by mm.
Here each variant is a fresh closed loop (Huber-D4 gated, as huber_windows.job, with its noise-free replay of the same commands), one factor off:
  base; no_breath (scene breathing component removed, cardiac kept); no_friction (mu_k = mu_s = k_d = 0); no_inertia (tool mass and its estimate 0);
  no_memory (viscoelastic branches off). Diagnostic counterfactuals, not a change of the settled setting.
For change windows (rising, falling, reversal outside events) of each run, noise-free local-linear endpoint errors over 24 and 64 ticks against that
run's own k_eff: the share of windows where the long span is harmful, median e24 and e64 there, and the share of harmful windows with common bias
(|e64 - e24| < .5 |e64|). Ribs sit at the same scan positions in every variant.
Writes results/v2_shared_bias_loop.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import huber_prior_check as H
import huber_windows as W
import v2_detection_cause as DC
import v2_history_ceiling as C
import v2_verify as V

HERE = Path(__file__).resolve().parent
OUT = Path(r'D:/pip_env/temp/claude/8fbcb2674095e6be/af85b341-794d-4ff3-a383-3afb99750be5/scratchpad/shared_bias_loop')
CASES = (0, 1, 2, 3, 4, 5, 36, 37, 38, 39, 40, 41)
VARIANTS = ('base', 'no_breath', 'no_friction', 'no_inertia', 'no_memory')


def modify(sc, variant):
    sc = dict(sc)
    if variant == 'no_breath':
        sc['u'] = sc['u']-np.asarray(sc.get('motion_breath', np.zeros_like(sc['u'])))
    elif variant == 'no_friction':
        sc.update(mu_k=0., mu_s=0., k_d=0.)
    elif variant == 'no_inertia':
        sc.update(tool_mass=0., tool_mass_hat=0.)
    if variant in ('no_memory', 'no_memory_no_pose'):
        sc['alpha'] = np.zeros_like(np.asarray(sc['alpha'], float)); sc.pop('memory_start', None)
    if variant in ('no_pose', 'no_memory_no_pose'):
        sc['theta0'] = 0.
    if variant == 'cover18':
        f = .018/float(sc['cover']); ribs = np.asarray(sc['ribs'], float).copy()
        if len(ribs):
            ribs[:, 2] *= f; ribs[:, 3] = .018
        sc['ribs'] = ribs; sc['cover'] = .018
    if variant == 'no_surface':
        hs = np.asarray(sc['h_s'], float); grid = sc['field_x0']+np.arange(len(hs))*sc['dx']
        radius = sc['surface_radius']; c = sc['surface_center']
        sc['h_s'] = -(grid-c)**2/(2.*radius) if np.isfinite(radius) else np.zeros_like(hs)
    return sc


def capture(task):
    variant, idx = task
    original = W.SE.MC.scene_for
    W.SE.MC.scene_for = lambda cell, base, var: modify(original(cell, base, var), variant)
    W.OUT = OUT/variant
    W.OUT.mkdir(parents=True, exist_ok=True)
    try:
        W.job(CC.cases()[idx])
    finally:
        W.SE.MC.scene_for = original
    return variant, idx


def analyse(task):
    variant, idx = task
    path = OUT/variant/f'case_{idx:02d}.npz'
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    dy0, z0 = C.rows(data, True)
    rows = []
    for j in range(100, n):
        if k[j] <= 200. or not any(seg[q][j] for q in ('rising', 'falling', 'reversal')):
            continue
        q24, q64 = DC.weights(dy0, t, j, 24), DC.weights(dy0, t, j, 64)
        r64 = np.arange(j-64, j)
        if q24 is None or q64 is None or not np.isfinite(z0[r64]).all():
            continue
        e24 = float(q24@z0[r64[-24:]])/k[j]-1; e64 = float(q64@z0[r64])/k[j]-1
        rows.append((e24, e64))
    return dict(variant=variant, index=idx, group=meta['group'], rows=rows, breath_pp=float(np.ptp(meta.get('breath_pp', [0])) if False else 0.))


def main():
    variants = VARIANTS if len(sys.argv) < 2 else tuple(sys.argv[1].split(','))
    cases = CASES if len(sys.argv) < 3 else tuple(int(x) for x in sys.argv[2].split(','))
    tasks = [(v, i) for v in variants for i in cases]
    with ProcessPoolExecutor(8) as ex:
        list(ex.map(capture, tasks, chunksize=1))
        res = list(ex.map(analyse, tasks, chunksize=1))
    rep = {}
    print('group variant | change windows | harmful share | among harmful: median e24, median e64, common-bias share | median |e24| all windows')
    for g in ('A', 'A_large'):
        for v in variants:
            a = np.array([x for r in res if r['variant'] == v and r['group'] == g for x in r['rows']])
            if not len(a):
                continue
            e24, e64 = a[:, 0], a[:, 1]
            harm = np.abs(e64)-np.abs(e24) > .05
            common = harm & (np.abs(e64-e24) < .5*np.abs(e64))
            cell = dict(n=len(a), harmful=float(harm.mean()), e24_harm=float(np.median(e24[harm])), e64_harm=float(np.median(e64[harm])), common=float(common.sum()/max(harm.sum(), 1)),
                        abs24=float(np.median(np.abs(e24))), abs64=float(np.median(np.abs(e64))))
            rep[f'{g}|{v}'] = cell
            print(f"  {g:8s} {v:12s} | {cell['n']:5d} | {cell['harmful']:.2f} | {cell['e24_harm']:+.3f} {cell['e64_harm']:+.3f} {cell['common']:.2f} | |e24| {cell['abs24']:.3f} |e64| {cell['abs64']:.3f}")
    path = HERE/'results/v2_shared_bias_loop.json'
    old = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    old.update(rep); path.write_text(json.dumps(old, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
