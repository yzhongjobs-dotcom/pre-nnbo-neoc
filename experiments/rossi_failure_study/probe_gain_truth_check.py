"""Measured y->F gain at 20.8 Hz against the logged truth |K(20.8 Hz)|, on saved 8 N closed-loop logs. No simulation.

Lock-in over 0.5 s (about ten probe cycles) with a cubic removed from each channel. Uses the plant's true axial force (col 18)
and true encoder (col 2), so the ratio is free of sensor noise. It tests whether the logged truth is the gain the probe sees.
Usage: python -X utf8 probe_gain_truth_check.py
"""
import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOGS = HERE/'results'/'baseline_logs'
F0 = 20.8; NW = 62


def lockin(t, y, f, k):
    out = []
    for i in range(int(1./.008), len(t)-NW, NW//2):
        tt = t[i:i+NW]; c = tt.mean(); s = tt-c
        A = np.c_[np.cos(2*np.pi*F0*tt), np.sin(2*np.pi*F0*tt), np.ones(NW), s, s**2, s**3]
        a = np.linalg.lstsq(A, y[i:i+NW], rcond=None)[0]; b = np.linalg.lstsq(A, f[i:i+NW], rcond=None)[0]
        ya = complex(a[0], -a[1]); fa = complex(b[0], -b[1])
        if abs(ya) < 2e-7 or f[i:i+NW].min() < .5:
            continue
        g = fa/ya
        out.append((c, abs(g), np.degrees(np.angle(g)), float(np.median(k[i:i+NW])), abs(ya)*1e6, abs(fa)))
    return np.array(out)


def main():
    res = {}
    print('scene | logged k | measured |G| | ratio med p10 p90 | phase deg | probe amp um | force amp N')
    for pat in ('breath_rib_rib', 'breath_fast_relaxed', 'breath_relaxed', 'rib_rib'):
        rows = []
        for p in sorted(glob.glob(str(LOGS/f'Huber-D4_{pat}_*.npz'))):
            z = np.load(p); L = z['log']; c = list(z['cols']); g = lambda i: L[:, c.index(i)]
            r = lockin(g(0), g(2), g(18), g(21))
            if len(r):
                rows.append(dict(file=Path(p).name, k=float(np.median(r[:, 3])), G=float(np.median(r[:, 1])), ratio=float(np.median(r[:, 1]/r[:, 3])), ratio_p10=float(np.percentile(r[:, 1]/r[:, 3], 10)), ratio_p90=float(np.percentile(r[:, 1]/r[:, 3], 90)), phase=float(np.median(r[:, 2])), amp_um=float(np.median(r[:, 4])), famp=float(np.median(r[:, 5]))))
        res[pat] = rows
        rows.sort(key=lambda x: x['k'])
        print(pat, 'n', len(rows))
        for x in rows[:: max(1, len(rows)//8)]+rows[-1:]:
            print(f"  {x['file'][9+len(pat):-4]:20s} k {x['k']:7.0f} G {x['G']:7.0f} ratio {x['ratio']:.2f} [{x['ratio_p10']:.2f},{x['ratio_p90']:.2f}] phase {x['phase']:6.1f} amp {x['amp_um']:6.1f} um  F {x['famp']:.3f} N")
    (HERE/'results'/'probe_gain_truth_check.json').write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
