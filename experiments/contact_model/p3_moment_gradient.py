"""Premise check: does the moment channel measure the path gradient of log stiffness?

For any symmetric contact pressure w(xi) over a footprint of half-width a, a modulus field that is locally
log-linear, E(x+xi) ~ E(x) exp(g xi), gives My/F ~ g <xi^2>_w and d ln k / dx = g (shift invariance), with
<xi^2> = a^2/4 for the Hertz elliptic pressure. So g_m = (My/F) / (a^2/4) should track the true d ln k/dx, with a
from the Hertz relation a = sqrt(R d) and d from the force/stiffness ratio. Here: closed-loop logs of the final
method (diagnosis S2 and T3m, 10 each); measured moment (column 13, sensor noise/bias/filter of the Nano25-grade
default) and a 5x noisier RM65-6F-grade copy; true d ln k/dx from column 8 along the path. Reported: correlation and
regression slope of g_m on the true gradient over ticks with a large true gradient, and the false-gradient level in
homogeneous stretches. Run: python -X utf8 p3_moment_gradient.py
"""
import math
import os
import sys
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE

R_PROBE = .03


def _ep(a):
    fam, sid, p, T = a
    log, _ = common.run(variants.build('Proposed-IMM', dict(CE.FINAL)), p, T, sid)
    t, x, F, My, kt = log[:, 0], log[:, 1], log[:, 5], log[:, 13], np.maximum(log[:, 8], 1.)
    sel = (t > 1.) & (F > 1.)
    # depth from force/stiffness (window estimate would be used online; here the true tangent stiffness for the
    # premise check), Hertz half-width
    d = np.clip(F/kt, 1e-5, 2e-2); a = np.sqrt(R_PROBE*d)
    rng = np.random.default_rng(sid)
    out = []
    for grade, extra in (('Nano25', 0.), ('RM65-6F', 4*7.5e-4)):
        Mym = My+extra*rng.standard_normal(len(My))*.3     # filtered-noise scale of the 30 Hz sensor filter
        gm = (Mym/np.maximum(F, .5))/(a*a/4)
        lk = np.log(kt); dx = np.gradient(x); gt = np.gradient(lk)/np.where(np.abs(dx) > 1e-7, dx, np.nan)
        # 40 ms moving average for both (same smoothing)
        w = 5; ker = np.ones(w)/w
        gm_s = np.convolve(gm, ker, 'same'); gt_s = np.convolve(np.nan_to_num(gt), ker, 'same')
        out.append((grade, gm_s[sel], gt_s[sel]))
    return fam, out


if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 10) if s[0] in ('S2', 'T3m')]
    with Pool(16) as po:
        rows = po.map(_ep, sc, chunksize=1)
    for fam in ('S2', 'T3m'):
        for gi, grade in enumerate(('Nano25', 'RM65-6F')):
            gm = np.concatenate([r[1][gi][1] for r in rows if r[0] == fam])
            gt = np.concatenate([r[1][gi][2] for r in rows if r[0] == fam])
            big = np.abs(gt) > 50.; flat = np.abs(gt) < 2.
            slope = float(np.sum(gm[big]*gt[big])/np.sum(gt[big]**2))
            corr = float(np.corrcoef(gm[big], gt[big])[0, 1])
            print('%-4s %-8s ticks with |true g|>50/m: %5d corr %.2f slope %.2f | homogeneous: |g_m| median %.1f p90 %.1f /m (true |g| p90 there %.1f)'
                  % (fam, grade, big.sum(), corr, slope, np.median(np.abs(gm[flat])), np.percentile(np.abs(gm[flat]), 90),
                     np.percentile(np.abs(gt[flat]), 90)), flush=True)
