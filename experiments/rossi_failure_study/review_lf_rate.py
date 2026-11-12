"""Verification on the stored closed loops (results/huber_windows): is the low-frequency contact disturbance the stiffness rate?

At constant force a multiplicative stiffness field gives  dF/dt = k_LF (dy/dt + du/dt) + n F dln(k)/dt.
Column 1 uses the true tissue motion u (physics only, sign of u found from tracking); column 2 uses the encoder alone with the true k.
Run from this folder: python -X utf8 review_lf_rate.py
"""
import json

import numpy as np

TC = .008


def ma(v, w=6):     # six ticks are one probe cycle, so the 20.8 Hz ripple is removed
    out = np.convolve(v, np.ones(w)/w)[:len(v)]; out[:w-1] = np.nan
    return out


print('case group shape phase k_med v | corr with true motion, gain | corr encoder only | body-motion share of the disturbance variance | rms dlnk/dt')
for i in list(range(0, 12))+list(range(36, 42)):
    d = np.load(f'results/huber_windows/case_{i:02d}.npz'); meta = json.loads(str(d['meta'])); log = d['log']
    t = log[:, 0]; y = ma(log[:, 2]); u = ma(log[:, 6]); F = np.maximum(ma(log[:, 7]), 1.); k = ma(np.maximum(log[:, 9], 1.))
    dF, dy, du, dl = (ma(np.gradient(v, TC), 5) for v in (F, y, u, np.log(k)))
    g = (t > 1.5) & (log[:, 9] > 0) & (log[:, 7] > .5) & np.isfinite(dl) & np.isfinite(dy) & np.isfinite(du)
    s = -np.sign(np.sum(dy[g]*du[g]))
    R = (dF-k*(dy+s*du))/F; D = (dF-k*dy)/F; mot = k*s*du/F
    print(f"{i:2d} {meta['group']:11s} {meta['cell'][2]:12s} {meta['cell'][3]:4.0f} {np.nanmedian(k[g]):6.0f} {np.median(np.abs(np.diff(log[:, 1])))/TC*1e3:5.1f} mm/s | "
          f"{np.corrcoef(R[g], dl[g])[0, 1]:+.2f} {np.polyfit(dl[g], R[g], 1)[0]:+.2f} | {np.corrcoef(D[g], dl[g])[0, 1]:+.2f} | {np.var(mot[g])/np.var(D[g]):.2f} | {np.sqrt(np.mean(dl[g]**2)):.2f}")
