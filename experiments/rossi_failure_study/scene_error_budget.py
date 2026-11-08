"""Where the force error comes from in the two core scenes, from existing closed-loop logs only.

Reads results/baseline_logs (overlay protocol, 8 N, 10 s, scored from t>=1 s) and splits the axial force error by
frequency band, by J term, and by whether the true contact gain is changing. No simulation is run.
Usage: python -X utf8 scene_error_budget.py
"""
import glob
import json
import re
from pathlib import Path

import numpy as np

import method3 as M

HERE = Path(__file__).resolve().parent
LOGS = HERE/'results'/'baseline_logs'
FD = 8.
BANDS = (('slow <2 Hz', 0., 2.), ('mid 2-15 Hz', 2., 15.), ('probe 15-27 Hz', 15., 27.), ('high >27 Hz', 27., 1e9))
CHANGING = .5    # |d ln k/dt| above this per second counts as a stiffness transition


def one(path):
    with np.load(path) as z:
        saved = z['log']; cols = list(z['cols'])
    g = lambda c: saved[:, cols.index(c)]
    t, f, k, depth, y = g(0), g(18), g(21), g(6), g(2)
    w = t >= 1.
    e = f[w]-FD
    dt = float(np.median(np.diff(t)))
    a, b, c = np.mean((e/M.TOL_F)**2), 20.*np.mean(f[w] < M.FMIN), np.mean(np.maximum(0., e-M.F_OVER)**2)
    spec = np.abs(np.fft.rfft(e-e.mean()))**2
    freq = np.fft.rfftfreq(len(e), dt)
    total = spec.sum()
    bands = {n: float(spec[(freq >= lo) & (freq < hi)].sum()/total) if total > 0 else 0. for n, lo, hi in BANDS}
    # smooth the true gain over one probe window before differentiating
    n = 24
    ls = np.convolve(np.log(np.maximum(k, 1.)), np.ones(n)/n, mode='same')
    rate = np.abs(np.gradient(ls, dt))
    rate[:n] = 0.; rate[-n:] = 0.    # the moving average falls off both ends
    moving = rate[w] > CHANGING
    return dict(file=Path(path).name, L=float(np.log1p(a+b+c)), J=float(a+b+c), J_err=float(a), J_loss=float(b), J_over=float(c),
                rmse=float(np.sqrt(np.mean(e**2))), mean_err=float(e.mean()), peak=float(f[w].max()), fmin=float(f[w].min()),
                bands=bands, moving_time=float(moving.mean()),
                moving_err_share=float(np.sum(e[moving]**2)/np.sum(e**2)),
                k_ratio=float(np.percentile(k[w], 99)/max(np.percentile(k[w], 1), 1.)),
                k_med=float(np.median(k[w])), rate_p99=float(np.percentile(rate[w], 99)),
                depth_range_mm=float(1e3*(np.percentile(depth[w], 99)-np.percentile(depth[w], 1))),
                y_range_mm=float(1e3*(np.percentile(y[w], 99)-np.percentile(y[w], 1))))


def summary(rows):
    arr = lambda key: np.array([r[key] for r in rows])
    J = arr('J')
    order = np.sort(J)[::-1]
    out = dict(n=len(rows), L_mean=float(arr('L').mean()), L_median=float(np.median(arr('L'))),
               rmse_median=float(np.median(arr('rmse'))), rmse_p90=float(np.percentile(arr('rmse'), 90)),
               peak_median=float(np.median(arr('peak'))), peak_max=float(arr('peak').max()),
               top1_J_share=float(order[0]/J.sum()), top3_J_share=float(order[:3].sum()/J.sum()),
               J_err_share=float(arr('J_err').sum()/J.sum()), J_loss_share=float(arr('J_loss').sum()/J.sum()), J_over_share=float(arr('J_over').sum()/J.sum()),
               moving_time=float(arr('moving_time').mean()), moving_err_share_median=float(np.median(arr('moving_err_share'))),
               k_ratio_median=float(np.median(arr('k_ratio'))), k_med_median=float(np.median(arr('k_med'))), rate_p99_median=float(np.median(arr('rate_p99'))),
               depth_range_mm_median=float(np.median(arr('depth_range_mm'))), y_range_mm_median=float(np.median(arr('y_range_mm'))))
    for n, _, _ in BANDS:
        v = np.array([r['bands'][n] for r in rows])
        out['band '+n+' median'] = float(np.median(v))
    return out


def main():
    res = {}
    for method in ('Huber-D4', 'Xing-DAHIC'):
        for scene, pat in (('A breath_rib', 'breath_rib_rib_*'), ('B breath_fast', 'breath_fast_relaxed_*'), ('A+B breath_rib_fast', 'breath_rib_fast_rib_*')):
            files = sorted(glob.glob(str(LOGS/f'{method}_{pat}.npz')))
            rows = [one(p) for p in files]
            res[f'{method} | {scene}'] = dict(summary(rows), episodes=rows)
            s = res[f'{method} | {scene}']
            print(f'{method} | {scene}:', {k: round(v, 3) for k, v in s.items() if k != 'episodes'})
    (HERE/'results'/'scene_error_budget.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
