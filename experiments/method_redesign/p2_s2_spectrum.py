"""Where (in frequency) is the stiffness-scan error of the final method versus original C1? Error PSD (Welch, true force)
summed into bands, S2 and T3m, diagnosis cohort 20 per family. Configs: final, TrueMeas (true measurement gain), C1-GN,
final with smoother motion model (q_smooth/10, r_scale x3). Run: python -X utf8 p2_s2_spectrum.py
"""
import os, sys
import numpy as np
from scipy.signal import welch
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
from common import P
from contact_update import ContactUpdate
F = dict(CE.FINAL)
CFG = {'Proposed': lambda: variants.build('Proposed-IMM', F),
       'TrueMeas': lambda: variants.build('Proposed-IMM', dict(F, oracle_meas=True)),
       'C1-GN': lambda: ContactUpdate(),
       'Smooth': lambda: variants.build('Proposed-IMM', dict(F, q_smooth=F['q_smooth']/10, r_scale=F['r_scale']*3))}
B = ((0, 2), (2, 5), (5, 10), (10, 19), (19, 23), (23, 40), (40, 62.5))
def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(CFG[cfg](), p, T, sid)
    e = log[log[:, 0] >= 1., 3]-P.FD
    f, S = welch(e-e.mean(), fs=125., nperseg=256)
    df = f[1]-f[0]
    return cfg, fam, [float(np.sqrt(np.sum(S[(f >= lo) & (f < hi)])*df)) for lo, hi in B], float(np.sqrt(np.mean(e**2)))
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('S2', 'T3m')]
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    print('band RMS (median over scenes), N; bands Hz:', B)
    for fam in ('S2', 'T3m'):
        for c in CFG:
            v = np.array([r[2] for r in rows if r[0] == c and r[1] == fam]); t = [r[3] for r in rows if r[0] == c and r[1] == fam]
            print('%-4s %-9s rmse %.3f/%.3f |' % (fam, c, np.median(t), np.mean(t)), ' '.join('%.3f' % x for x in np.median(v, 0)), flush=True)
