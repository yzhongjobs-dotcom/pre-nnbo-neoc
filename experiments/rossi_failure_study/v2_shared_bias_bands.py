"""Where in frequency is the shared bias, and is the long-span harm in A_large the curvature of k itself?

On the noise-free windows of the controlled loops (v2_shared_bias_loop.py, variants base and no_memory), local-linear endpoint estimator q over 24 or 64 ticks:
  e = q'(z0 - k_T h0 - v h1)/k_T split by DFT bins of the window (q'x = sum over bins), into low (<= 10.4 Hz), probe (15.6-26 Hz), high (>= 31 Hz), and
  curvature part q' cumsum((k - k_T - v tau) dY0)/k_T, v the window's own least-squares slope of k.
Prediction: in A the shared bias sits mainly outside the probe band and shrinks without memory; in A_large e64 is mostly curvature.
Writes results/v2_shared_bias_bands.json.
"""
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import v2_detection_cause as DC
import v2_history_ceiling as C
import v2_shared_bias_loop as SL
import v2_verify as V

HERE = Path(__file__).resolve().parent
BANDS = dict(low=lambda f: f <= 10.5, probe=lambda f: (f > 15) & (f < 27), high=lambda f: f >= 31)


def split(q, x):
    Q, X = np.fft.rfft(q), np.fft.rfft(x)
    part = np.real(np.conj(Q)*X)*2/len(q); part[0] /= 2
    if len(q) % 2 == 0:
        part[-1] /= 2
    f = np.fft.rfftfreq(len(q), .008)
    return {b: float(part[m(f)].sum()) for b, m in BANDS.items()}


def job(task):
    variant, idx = task
    data = H.load(SL.OUT/variant/f'case_{idx:02d}.npz'); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    dy0, z0 = C.rows(data, True)
    out = []
    for j in range(100, n):
        if k[j] <= 200. or not any(seg[q][j] for q in ('rising', 'falling', 'reversal')):
            continue
        rec = {}
        for L in (24, 64):
            q = DC.weights(dy0, t, j, L); r = np.arange(j-L, j)
            if q is None or not np.isfinite(z0[r]).all():
                rec = None; break
            y = np.cumsum(dy0[r]); h0 = y-y[-1]; tau = t[r]+.0035-t[j-1]
            kr = .5*(k[r]+k[r+1]); v = np.polyfit(tau, kr, 1)[0]
            curv = np.r_[0., np.cumsum((kr-k[j]-v*tau)[1:]*np.diff(y))]
            resid = z0[r]-k[j]*h0-v*tau*h0
            parts = split(q, resid-curv)
            rec[L] = dict(e=float(q@resid)/k[j], curv=float(q@curv)/k[j], **{b: p/k[j] for b, p in parts.items()})
        if rec:
            out.append(rec)
    return dict(variant=variant, index=idx, group=meta['group'], rows=out)


def main():
    tasks = [(v, i) for v in ('base', 'no_memory') for i in SL.CASES]
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(job, tasks, chunksize=1))
    rep = {}
    print('mean over change windows of each part of the noise-free endpoint error (relative), and rms; harmful = |e64| > |e24| + .05')
    for g in ('A', 'A_large'):
        for v in ('base', 'no_memory'):
            rows = [x for r in res if r['variant'] == v and r['group'] == g for x in r['rows']]
            harm = [x for x in rows if abs(x[64]['e'])-abs(x[24]['e']) > .05]
            for label, sel in (('all', rows), ('harmful', harm)):
                for L in (24, 64):
                    cell = {p: dict(mean=float(np.mean([x[L][p] for x in sel])), rms=float(np.sqrt(np.mean([x[L][p]**2 for x in sel])))) for p in ('e', 'curv', 'low', 'probe', 'high')}
                    rep[f'{g}|{v}|{label}|{L}'] = dict(n=len(sel), **cell)
                    print(f"  {g:8s} {v:10s} {label:8s} L{L} n {len(sel):5d} | " + ' '.join(f"{p} {cell[p]['mean']:+.3f}({cell[p]['rms']:.3f})" for p in ('e', 'curv', 'low', 'probe', 'high')))
    (HERE/'results/v2_shared_bias_bands.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
