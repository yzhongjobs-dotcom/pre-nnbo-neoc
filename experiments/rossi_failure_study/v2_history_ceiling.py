"""Information ceiling for the slope: what if history were cut exactly at the true start of the current monotone run?

Slope information grows as T^3 with the data span T, so the only lever inside the stiffness channel is a longer history within one run.
Estimators on the recorded rows (noisy and noise-free), W = I (Huber weights are inactive outside events), D4 with the same rho:
  C24      constant k, last 24 ticks                       (Huber-D4 outside events)
  C_run    constant k over the oracle run span              (longer history alone)
  L_run    k and free slope over the oracle run span        (perfect segmentation, V2's model)
  R_run    k and slope with the oracle ridge sv2 = v^2 over the run span
  L48      k and free slope over a fixed 48 ticks           (history without segmentation)
Run span: from the first tick of the current sign run of the 0.2 s smoothed ln k rate (|rate| > .2), at least 24 and at most 96 ticks.
Predictions fixed before running: for |v/k| near 1 /s, L_run and R_run do not beat C24 in rms; a gain appears only above about 2.5 /s; L48 is worse at reversals.
Writes results/v2_history_ceiling.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import huber_prior_check as H
import v2_verify as V

HERE = Path(__file__).resolve().parent
TC, SIG, RHO = .008, H.SIGMA, H.RHO
_M = {}


def mat(m):
    if m not in _M:
        D = np.diff(np.eye(m), 4, axis=0)
        _M[m] = np.eye(m)-np.linalg.inv(np.eye(m)+RHO*D.T@D)
    return _M[m]


def rows(data, noise_free):
    Y = data['Y0' if noise_free else 'Y']; Z = data['z0' if noise_free else 'z']
    n = len(Y); dy = np.full(n, np.nan); z = np.full(n, np.nan)
    for j in range(n):
        if np.isfinite(Y[j, -1]) and np.isfinite(Y[j, -2]):
            dy[j-1] = Y[j, -1]-Y[j, -2]; z[j-1] = Z[j, -1]
    return dy, z


def fit(dy, z, t, j, L, mode, v2=None):
    """Window of rows j-L .. j-1 (row r spans log r to r+1), endpoint t[j-1]."""
    r = np.arange(j-L, j)
    if r[0] < 1 or not (np.isfinite(dy[r]).all() and np.isfinite(z[r]).all()):
        return np.nan
    y = np.cumsum(dy[r]); h0 = y-y[-1]; zz = z[r]
    tt = t[r]+.0035-t[j-1]
    h1 = tt*h0
    M = mat(L); mh0, mh1 = M@h0, M@h1
    a, b, c1 = h0@mh0, h0@mh1, h1@mh1
    p0, p1 = h0@(M@zz), h1@(M@zz)
    if mode == 'const':
        return p0/a
    if mode == 'free':
        return (p0*c1-b*p1)/(a*c1-b*b)
    cp = c1+SIG**2/max(v2, 1e-9)
    return (p0-(b/cp)*p1)/(a-b*b/cp)


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    ell = np.log(np.maximum(k, 1.)); s = uniform_filter1d(ell, 25, mode='nearest'); rate = np.gradient(s, t)
    sign = np.where(rate > .2, 1, np.where(rate < -.2, -1, 0))
    start = np.zeros(n, int); cur = 0
    for j in range(1, n):
        if sign[j] != sign[j-1]:
            cur = j
        start[j] = cur
    out = {}
    for nf in (False, True):
        dy, z = rows(data, nf)
        est = {q: np.full(n, np.nan) for q in ('C24', 'C_run', 'L_run', 'R_run', 'L48')}
        vrel = np.full(n, np.nan); span = np.full(n, np.nan)
        for j in range(100, n):
            if k[j] <= 200. or not any(seg[q][j] for q in ('rising', 'falling', 'reversal', 'plateau')):
                continue
            L = int(np.clip(j-start[j], 24, 96))
            rr = np.arange(j-L, j)
            if np.any(k[rr] <= 200.):
                continue
            kr = .5*(k[rr]+k[rr+1]); v = float(np.polyfit(t[rr]+.0035-t[j-1], kr, 1)[0])
            est['C24'][j] = fit(dy, z, t, j, 24, 'const'); est['C_run'][j] = fit(dy, z, t, j, L, 'const')
            est['L_run'][j] = fit(dy, z, t, j, L, 'free'); est['R_run'][j] = fit(dy, z, t, j, L, 'ridge', v*v)
            est['L48'][j] = fit(dy, z, t, j, 48, 'free')
            vrel[j] = v/k[j]; span[j] = L*TC
        out['nf' if nf else 'noisy'] = {q: np.log(np.maximum(v, 200.)/k).tolist() for q, v in est.items()}
        out['vrel'] = vrel.tolist(); out['span'] = span.tolist()
    return dict(index=meta['index'], group=meta['group'], seg={q: m.tolist() for q, m in seg.items()}, **out)


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    names = ('C24', 'C_run', 'L_run', 'R_run', 'L48')
    report = {}
    print('ln error |e| median / rms / bias by segment and |v/k| bin; noisy, then noise-free C24 and L_run medians; median span s')
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        for lo, hi in ((0, .5), (.5, 1.5), (1.5, 2.5), (2.5, 99)):
            cell = {}
            for src in ('noisy', 'nf'):
                for q in names:
                    x = np.concatenate([np.asarray(r[src][q])[np.asarray(r['seg'][s]) & (np.abs(np.asarray(r['vrel'])) >= lo) & (np.abs(np.asarray(r['vrel'])) < hi)] for r in res])
                    x = x[np.isfinite(x)]
                    if len(x) < 15:
                        continue
                    cell[f'{src}_{q}'] = dict(n=len(x), med=float(np.median(np.abs(x))), rms=float(np.sqrt(np.mean(x*x))), bias=float(np.median(x)))
            if not cell:
                continue
            sp = np.concatenate([np.asarray(r['span'])[np.asarray(r['seg'][s]) & (np.abs(np.asarray(r['vrel'])) >= lo) & (np.abs(np.asarray(r['vrel'])) < hi)] for r in res])
            cell['span'] = float(np.nanmedian(sp))
            report[f'{s}|{lo}-{hi}'] = cell
            f = lambda q: f"{cell[q]['med']:.3f}/{cell[q]['rms']:.3f}/{cell[q]['bias']:+.3f}" if q in cell else '-'
            print(f"{s:9s} |v/k| {lo}-{hi:<4} n {cell.get('noisy_C24', {}).get('n', 0):5d} span {cell['span']:.2f} | "+' '.join(f"{q} {f('noisy_'+q)}" for q in names)+f" | nf C24 {f('nf_C24')} L_run {f('nf_L_run')}")
    (HERE/'results/v2_history_ceiling.json').write_text(json.dumps(report, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
