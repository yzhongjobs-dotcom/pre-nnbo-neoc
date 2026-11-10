"""Adaptive history for the slope without knowing where the run starts: intersection of confidence intervals (ICI) over local-linear spans.

Root cause established in v2_change_decomposition.py and v2_history_ceiling.py: one 24-tick window has slope SNR below 1; slope information grows as T^3,
so the slope needs history, but history is only valid inside the current monotone run, whose start is unknown. ICI (Goldenshluger-Nemirovski,
Katkovnik) picks, for the endpoint value, the largest span whose confidence interval still intersects those of all shorter spans; its risk is
within a log factor of the best span chosen with hindsight. Spans 24, 32, 48, 64, 96 ticks; local-linear D4 fit per span (W = I);
endpoint sd from the exact linear weights, sigma = measured tick noise. Gamma chosen on even case indices, reported on odd ones.
Predictions fixed before running: for |v/k| >= 1.5 rises and reversals, error between C24 and the oracle-run ceiling; plateau equal to C24;
falls still gain little, since the excitation deficit on falls is untouched.
Writes results/v2_ici_replay.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import v2_history_ceiling as C
import v2_verify as V

HERE = Path(__file__).resolve().parent
TC, SIG = .008, H.SIGMA
SPANS = (24, 32, 48, 64, 96)
GAMMAS = (1.0, 1.5, 2.0, 2.5, 3.0)


def local_linear(dy, z, t, j, L):
    r = np.arange(j-L, j)
    if r[0] < 1 or not (np.isfinite(dy[r]).all() and np.isfinite(z[r]).all()):
        return None
    y = np.cumsum(dy[r]); h0 = y-y[-1]; zz = z[r]; h1 = (t[r]+.0035-t[j-1])*h0
    M = C.mat(L); mh0, mh1 = M@h0, M@h1
    a, b, c1 = h0@mh0, h0@mh1, h1@mh1
    det = a*c1-b*b
    if det <= 0:
        return None
    q = (c1*mh0-b*mh1)/det
    k = float(q@zz)
    return k, SIG*math.sqrt(float(q@q)), float(mh0@zz/a)


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    dy, z = C.rows(data, False)
    fits = np.full((n, len(SPANS), 2), np.nan); c24 = np.full(n, np.nan); vrel = np.full(n, np.nan)
    for j in range(100, n):
        if k[j] <= 200.:
            continue
        for i, L in enumerate(SPANS):
            f = local_linear(dy, z, t, j, L)
            if f is not None:
                fits[j, i] = f[:2]
                if L == 24:
                    c24[j] = f[2]
        rr = np.arange(j-24, j); kr = .5*(k[rr]+k[rr+1])
        vrel[j] = float(np.polyfit(t[rr]-t[j-1], kr, 1)[0])/k[j]
    ici = {}
    for g in GAMMAS:
        out = np.full(n, np.nan)
        for j in range(n):
            lo, hi = -np.inf, np.inf; best = np.nan
            for i in range(len(SPANS)):
                kk, s = fits[j, i]
                if not np.isfinite(kk):
                    break
                lo, hi = max(lo, kk-g*s), min(hi, kk+g*s)
                if lo > hi:
                    break
                best = kk
            out[j] = best
        ici[g] = out
    lk = np.log(np.where(k > 200., k, np.nan))
    err = lambda x: np.log(np.clip(x, 200., 3e5))-lk
    return dict(index=meta['index'], group=meta['group'], seg={q: m.tolist() for q, m in seg.items()}, vrel=vrel.tolist(),
                C24=err(c24).tolist(), L24=err(fits[:, 0, 0]).tolist(), **{f'ICI{g}': err(v).tolist() for g, v in ici.items()})


def stats(x):
    x = x[np.isfinite(x)]
    return dict(n=len(x), med=float(np.median(np.abs(x))), rms=float(np.sqrt(np.mean(x*x))), bias=float(np.median(x))) if len(x) > 10 else None


def collect(res, key, s, lo=0., hi=99., parity=None):
    xs = []
    for r in res:
        if parity is not None and r['index'] % 2 != parity:
            continue
        m = np.asarray(r['seg'][s]) & (np.abs(np.asarray(r['vrel'])) >= lo) & (np.abs(np.asarray(r['vrel'])) < hi)
        xs.append(np.asarray(r[key])[m])
    return np.concatenate(xs)


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    # Development objective: mean rms over plateau, rising, falling, reversal, even cases.
    obj = {}
    for g in GAMMAS:
        obj[g] = float(np.mean([stats(collect(res, f'ICI{g}', s, parity=0))['rms'] for s in ('plateau', 'rising', 'falling', 'reversal')]))
    gstar = min(obj, key=obj.get)
    print('development objective by gamma', {g: round(v, 4) for g, v in obj.items()}, 'chosen', gstar)
    rep = dict(objective=obj, gamma=gstar, test={})
    print('test half (odd cases): |e| median / rms / bias')
    for s in ('plateau', 'rising', 'falling', 'reversal', 'event'):
        for lo, hi in ((0, 99), (0, 1.5), (1.5, 99)):
            cell = {q: stats(collect(res, q, s, lo, hi, parity=1)) for q in ('C24', 'L24', f'ICI{gstar}')}
            if cell['C24'] is None:
                continue
            rep['test'][f'{s}|{lo}-{hi}'] = cell
            print(f"  {s:9s} |v/k| {lo}-{hi:<3}  n {cell['C24']['n']:5d} | "+' '.join(f"{q} {c['med']:.3f}/{c['rms']:.3f}/{c['bias']:+.3f}" for q, c in cell.items() if c))
    (HERE/'results/v2_ici_replay.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
