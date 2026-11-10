"""Can the slow force channel time the start of a monotone stiffness run? Causal band-passed psi = (dF - k_hat dy)/k_hat (v2_slope_sources.py).
For each true run (|d ln k| >= .15 over the run, 0.2 s smoothed rate beyond .2/s), the nearest causal psi sign change of the same direction:
lag after the true start. And false boundaries: psi sign changes per second on plateaus outside events.
Writes results/v2_cue_timing.json.
"""
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import huber_prior_check as H
import v2_slope_sources as S
import v2_verify as V

HERE = Path(__file__).resolve().parent


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]; F = log[:, 5]; y = log[:, 2]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    raw = np.asarray(data['raw'], float); good = np.isfinite(raw) & (raw > 0)
    kh = uniform_filter1d(np.interp(np.arange(n), np.flatnonzero(good), raw[good]), 12, mode='nearest')
    dF = np.gradient(S.band(F, .05, 5.), t); dy = np.gradient(S.band(y, .05, 5.), t)
    psi = S.band((dF-kh*dy)/kh, .5, 5., causal=True)
    ell = np.log(np.maximum(k, 1.)); s = uniform_filter1d(ell, 25, mode='nearest'); r = np.gradient(s, t)
    sign = np.where(r > .2, 1, np.where(r < -.2, -1, 0))
    starts = [j for j in range(1, n) if sign[j] != 0 and sign[j] != sign[j-1] and t[j] > 1.2 and k[j] > 200.]
    runs = []
    for j in starts:
        e = j
        while e < n-1 and sign[e+1] == sign[j]:
            e += 1
        if abs(s[e]-s[j]) >= .15:
            runs.append((j, e, sign[j]))
    sp = np.sign(psi)
    flips = [j for j in range(1, n) if sp[j] != sp[j-1]]
    lags = []
    for j, e, sg in runs:
        cand = [f for f in flips if t[j]-.2 <= t[f] <= t[e] and sp[f] == sg]
        lags.append(float(t[cand[0]]-t[j]) if cand else None)
    pl = np.asarray(seg['plateau'])
    false = sum(1 for f in flips if pl[f])/max(pl.sum()*.008, 1e-9)
    return dict(group=meta['group'], lags=lags, false_per_s=false, plateau_s=float(pl.sum()*.008))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rep = {}
    for g in ('A', 'A_slow_thin', 'A_large', 'all'):
        rr = [r for r in res if g == 'all' or r['group'] == g]
        lags = [x for r in rr for x in r['lags']]
        found = [x for x in lags if x is not None]
        fp = [r['false_per_s'] for r in rr if r['plateau_s'] > 1]
        rep[g] = dict(runs=len(lags), found=len(found), lag_median=float(np.median(found)) if found else None, lag_p90=float(np.percentile(found, 90)) if found else None,
                      within_50ms=float(np.mean(np.abs(found) <= .05)) if found else None, false_flips_per_s_plateau=float(np.median(fp)) if fp else None)
        print(g, rep[g])
    (HERE/'results/v2_cue_timing.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
