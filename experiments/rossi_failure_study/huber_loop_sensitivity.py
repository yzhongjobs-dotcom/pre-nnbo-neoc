"""Are the large closed-loop differences between Huber-D4 and same-config D4 caused by the prior, or by a loop that amplifies anything?

Takes the cases with the largest |dL| in results/huber_same_config_loop.json and reruns each of the two arms with the stiffness estimate
multiplied by (1 + eps), eps = 1e-6 and -1e-6: far below any estimator difference, a pure test of sensitivity.
If L moves by as much as the Huber-D4 difference, that difference says nothing about the prior.
Writes results/huber_loop_sensitivity.json.
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
import huber_same_config_loop as L
import overlay_compare as O
import ridge_estimator as RE
import d_robust as DR

HERE = Path(__file__).resolve().parent
EPS = (0., 1e-6, -1e-6)
LARGER = (1e-3, -1e-3, 5e-3, -5e-3)


def run(task):
    arm, case, eps = task
    originals = RE.RidgeProjection.estimate, DR.RobustRidge.estimate

    def wrap(fn):
        def inner(self, o, n, policy):
            k, sigma, chi = fn(self, o, n, policy)
            return k*(1.+eps), sigma, chi
        return inner

    RE.RidgeProjection.estimate = wrap(originals[0]); DR.RobustRidge.estimate = wrap(originals[1])
    try:
        row = L.worker((arm, case))
    finally:
        RE.RidgeProjection.estimate, DR.RobustRidge.estimate = originals
    return dict(arm=arm, index=case['index'], group=case['group'], eps=eps, L=row['L'], peak=row['peak'])


def main():
    data = json.loads((HERE/'results/huber_same_config_loop.json').read_text(encoding='utf-8'))
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    now = {(r['method'], r['index']): r for r in data['rows']}
    gaps = []
    for c in CC.cases():
        i = c['index']
        gaps.append((abs(stored[('Huber-D4', i)]['L']-now[('D4', i)]['L']), 'gated', c))
        gaps.append((abs(stored[('Huber-D4_nogate', i)]['L']-now[('D4_nogate', i)]['L']), 'nogate', c))
    gaps.sort(key=lambda g: -g[0])
    picked = gaps[:8]
    tasks = []
    for gap, mode, c in picked:
        for eps in EPS:
            tasks.append(('Huber-D4' if mode == 'gated' else 'Huber-D4_nogate', c, eps))
            tasks.append(('D4' if mode == 'gated' else 'D4_nogate', c, eps))
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(run, tasks, chunksize=1))
    out = []
    print('group idx mode | prior gap dL | Huber L at 0,+,- | D4 L at 0,+,- | Huber spread, D4 spread')
    for gap, mode, c in picked:
        pick = lambda prefix: [r['L'] for r in rows if r['index'] == c['index'] and r['arm'].startswith(prefix) and r['arm'].endswith('_nogate') == (mode == 'nogate')]
        h, d = pick('Huber'), pick('D4')
        s = dict(group=c['group'], index=c['index'], mode=mode, gap=gap, huber=h, d4=d, spread_huber=float(np.ptp(h)), spread_d4=float(np.ptp(d)))
        out.append(s)
        print(f"{c['group']:12s} {c['index']:2d} {mode:6s} | {gap:.3f} | {h[0]:.3f} {h[1]:.3f} {h[2]:.3f} | {d[0]:.3f} {d[1]:.3f} {d[2]:.3f} | {s['spread_huber']:.3f} {s['spread_d4']:.3f}")
    # The three differences that survive 1e-6: do they survive a perturbation still far below the noise?
    survivors = [g for g in picked if g[0] > .2 and max(next(o for o in out if o['index'] == g[2]['index'] and o['mode'] == g[1])[k] for k in ('spread_huber', 'spread_d4')) < .3*g[0]]
    tasks = []
    for gap, mode, c in survivors:
        for eps in LARGER:
            tasks.append(('Huber-D4' if mode == 'gated' else 'Huber-D4_nogate', c, eps))
            tasks.append(('D4' if mode == 'gated' else 'D4_nogate', c, eps))
    with ProcessPoolExecutor(8) as ex:
        more = list(ex.map(run, tasks, chunksize=1))
    print('larger perturbations, eps 1e-3, -1e-3, 5e-3, -5e-3')
    larger = []
    for gap, mode, c in survivors:
        pick = lambda prefix: [r['L'] for r in more if r['index'] == c['index'] and r['arm'].startswith(prefix) and r['arm'].endswith('_nogate') == (mode == 'nogate')]
        h, d = pick('Huber'), pick('D4')
        larger.append(dict(group=c['group'], index=c['index'], mode=mode, gap=gap, huber=h, d4=d))
        print(f"{c['group']:12s} {c['index']:2d} {mode:6s} gap {gap:.3f} | Huber "+' '.join(f'{x:.3f}' for x in h)+' | D4 '+' '.join(f'{x:.3f}' for x in d))
    (HERE/'results/huber_loop_sensitivity.json').write_text(json.dumps(dict(eps=EPS, cases=out, larger_eps=LARGER, larger=larger), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
