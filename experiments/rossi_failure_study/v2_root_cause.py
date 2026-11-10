"""Step 2: which property of the stiffness the controller uses explains force control, on the original method's own closed loops.

Arms: Huber-D4 ungated (24 ticks) and Huber-D4 ungated + IMM, unchanged (CC.job definitions). Per tick: k_good, k_eff, force.
Per case: plateau spread of ln(k_good/k_eff), sustained under-estimate (worst 0.1 s mean below truth), sustained over-estimate, and their timing
relative to the force peak. Question: does L follow plateau noise, under-estimates or over-estimates? No tuning, no new estimator.
Writes results/v2_root_cause.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.stats import spearmanr

import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O
import rib_transition_eval as R

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent


def job(task):
    arm, case = task
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    if arm == 'imm':
        R.make_filter(pol, 'imm')
    trace = []
    step = pol.step

    def capture(o):
        u = step(o); trace.append(pol.k_good); return u
    pol.step = capture
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    n = min(len(log), len(trace))
    kg = np.asarray(trace[:n], float); t = log[:n, 0]; k = log[:n, 26]; F = log[:n, 18]
    _, masks = SE.segment_masks(log[:n], sc)
    ok = (k > 200.) & np.isfinite(kg) & (kg > 0) & (t >= 1.)
    e = np.full(n, np.nan); e[ok] = np.log(kg[ok]/k[ok])
    sm = uniform_filter1d(np.nan_to_num(e), 13, mode='nearest')        # 0.1 s mean
    sm[~ok] = np.nan
    pk = int(np.nanargmax(np.abs(F-sc['FD'])))
    win = (t >= t[pk]-.3) & (t <= t[pk])
    pe = e[masks['plateau'][:n] & ok]
    return dict(arm=arm, index=case['index'], group=case['group'], L=float(np.log1p(m['J'])), peak=float(m['peak']), rmse=float(m['rmse']),
                plateau_sd=float(np.std(pe)) if len(pe) > 20 else np.nan, under=float(max(0., -np.nanmin(sm))), over=float(max(0., np.nanmax(sm))),
                under_pre=float(max(0., -np.nanmin(sm[win]))) if np.isfinite(sm[win]).any() else np.nan, over_pre=float(max(0., np.nanmax(sm[win]))) if np.isfinite(sm[win]).any() else np.nan,
                frac_under=float(np.nanmean(e[ok] < math.log(.7))), frac_over=float(np.nanmean(e[ok] > math.log(1.4))))


def main():
    cases = CC.cases()
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, c) for a in ('nogate', 'imm') for c in cases], chunksize=1))
    stored = {(r['method'], r['index']): r['L'] for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    worst = max(abs(r['L']-stored[('Huber-D4_nogate' if r['arm'] == 'nogate' else 'Huber-D4_imm', r['index'])]) for r in rows)
    print('regression to stored L, max |dL|', worst)
    rep = dict(regression=worst, rows=rows, corr={})
    keys = ('plateau_sd', 'under', 'over', 'under_pre', 'over_pre', 'frac_under', 'frac_over')
    print('Spearman with L across 84 runs (and excluding combo):')
    for kk in keys:
        a = np.array([r[kk] for r in rows]); L = np.array([r['L'] for r in rows]); ok = np.isfinite(a)
        nc = np.array([r['group'] != 'combo' for r in rows]) & ok
        c1 = spearmanr(a[ok], L[ok]).correlation; c2 = spearmanr(a[nc], L[nc]).correlation
        rep['corr'][kk] = dict(all=float(c1), no_combo=float(c2))
        print(f'  {kk:11s} {c1:+.2f}  {c2:+.2f}')
    print('worst cases by L: arm idx group L peak | under over (0.1 s) | before the peak under over | plateau sd')
    for r in sorted(rows, key=lambda r: -r['L'])[:16]:
        print(f"  {r['arm']:6s} {r['index']:2d} {r['group']:12s} {r['L']:.3f} {r['peak']:6.1f} | {r['under']:.2f} {r['over']:.2f} | {r['under_pre']:.2f} {r['over_pre']:.2f} | {r['plateau_sd']:.3f}")
    (HERE/'results/v2_root_cause.json').write_text(json.dumps(rep, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
