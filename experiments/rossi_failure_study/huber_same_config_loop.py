"""Closed loop: D4 with Huber-D4's window and gate scales, so only the prior differs.

Table 1's D4 row uses nw 32 and sigma_scale 1.5; Huber-D4 uses nw 24 and sigma_scale 0.796.
Predictions fixed before running, from the prior-level replay: scene B identical (no Huber weight ever
activates); A, A_slow_thin, A+B near zero; only cough and combo separate, more with gates than without,
because Huber also changes chi. Huber rows come from completion_closed_loop.json; the Huber-D4 row is
rerun here as a regression. Writes results/huber_same_config_loop.json.
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
import overlay_compare as O

HERE = Path(__file__).resolve().parent
SAME = dict(nw=24, fc=14., gate=.551, sigma_scale=.796)
GROUPS = ('A', 'A_slow_thin', 'B', 'A+B', 'cough', 'combo', 'A_large')


def worker(task):
    method, case = task
    original = O.chosen
    configs = original()
    configs['D4'] = dict(SAME)
    O.chosen = lambda: configs
    try:
        row = CC.job((method, case))
    finally:
        O.chosen = original
    return dict(method=method, index=case['index'], group=case['group'], L=row['L'], peak=row['peak'], force_rmse=row['force_rmse'])


def main():
    assert O.chosen()['Huber-D4']['nw'] == 24 and O.chosen()['Huber-D4']['sigma_scale'] == SAME['sigma_scale'] and O.chosen()['Huber-D4']['gate'] == SAME['gate']
    tasks = [(m, c) for m in ('D4', 'D4_nogate', 'Huber-D4') for c in CC.cases()]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(worker, tasks, chunksize=1))
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    worst = max(abs(r['L']-stored[('Huber-D4', r['index'])]['L']) for r in rows if r['method'] == 'Huber-D4')
    assert worst == 0., worst
    now = {(r['method'], r['index']): r for r in rows}
    summary = {}
    print('regression Huber-D4 rerun vs stored, max |dL|', worst)
    print('group: gated Huber L, D4_24 L, D4_32 L (table) | no-gate Huber L, D4_24 L, D4_32 L (table) | identical cases gated, no-gate; max |dL| | peaks')
    for g in GROUPS+('all',):
        cs = [c['index'] for c in CC.cases() if g == 'all' or c['group'] == g]
        val = lambda m, src: np.array([src[(m, i)]['L'] for i in cs])
        hg, dg, tg = val('Huber-D4', stored), val('D4', now), val('D4', stored)
        hn, dn, tn = val('Huber-D4_nogate', stored), val('D4_nogate', now), val('D4_nogate', stored)
        summary[g] = dict(huber_gated=float(hg.mean()), d4_24_gated=float(dg.mean()), d4_32_gated=float(tg.mean()),
                          huber_nogate=float(hn.mean()), d4_24_nogate=float(dn.mean()), d4_32_nogate=float(tn.mean()),
                          same_gated=int(np.sum(np.abs(hg-dg) < 1e-9)), same_nogate=int(np.sum(np.abs(hn-dn) < 1e-9)), n=len(cs),
                          max_gated=float(np.max(np.abs(hg-dg))), max_nogate=float(np.max(np.abs(hn-dn))),
                          peak_huber_gated=max(stored[('Huber-D4', i)]['peak'] for i in cs), peak_d4_24_gated=max(now[('D4', i)]['peak'] for i in cs),
                          peak_huber_nogate=max(stored[('Huber-D4_nogate', i)]['peak'] for i in cs), peak_d4_24_nogate=max(now[('D4_nogate', i)]['peak'] for i in cs))
        s = summary[g]
        print(f"{g:12s} {s['huber_gated']:.3f} {s['d4_24_gated']:.3f} {s['d4_32_gated']:.3f} | {s['huber_nogate']:.3f} {s['d4_24_nogate']:.3f} {s['d4_32_nogate']:.3f} | same {s['same_gated']}/{s['n']} {s['same_nogate']}/{s['n']} max {s['max_gated']:.3f} {s['max_nogate']:.3f} | peak {s['peak_huber_gated']:.1f} {s['peak_d4_24_gated']:.1f} / {s['peak_huber_nogate']:.1f} {s['peak_d4_24_nogate']:.1f}")
    (HERE/'results/huber_same_config_loop.json').write_text(json.dumps(dict(config=SAME, regression=worst, summary=summary, rows=rows), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
