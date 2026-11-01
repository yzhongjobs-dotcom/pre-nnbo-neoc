"""Decompose the stiffness-scan RMSE gap between original C1 (GN, dither 0.1 N) and the final method (dither 0.17 N):
excitation ripple at 20.8 Hz, slow error (<5 Hz moving average) and remaining band. Diagnosis cohort, 20 per family.
Run: python -X utf8 p2_c1_decompose.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
from common import P
from contact_update import ContactUpdate
def parts(log, Tc=.008, warm=1.):
    sel = log[:, 0] >= warm; e = log[sel, 3]-P.FD; t = log[sel, 0]
    w = int(round(.2/Tc)); slow = np.convolve(e, np.ones(w)/w, mode='same')[w:-w]; hp = (e[w:-w]-slow); tt = t[w:-w]
    A = np.vstack([np.cos(2*np.pi*20.8*tt), np.sin(2*np.pi*20.8*tt)]).T; c, *_ = np.linalg.lstsq(A, hp, rcond=None)
    rip = A@c; rest = hp-rip
    return [float(np.sqrt(np.mean(x**2))) for x in (e[w:-w], slow, rip, rest)]
CFG = {'Proposed': lambda: variants.build('Proposed-IMM', dict(CE.FINAL)),
       'Proposed-dith0.10': lambda: variants.build('Proposed-IMM', dict(CE.FINAL, a_force=.10)),
       'C1-GN': lambda: ContactUpdate(), 'C1-GN-dith0.17': lambda: ContactUpdate(Fdith=.17)}
def _ep(a):
    cfg, (fam, sid, p, T) = a
    log = P.episode(p, CFG[cfg](), T=T, seed=sid)
    return cfg, fam, parts(log) if log is not None else [np.nan]*4
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('S2', 'T3m')]
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    for fam in ('S2', 'T3m'):
        for c in CFG:
            v = np.array([r[2] for r in rows if r[0] == c and r[1] == fam])
            print('%-4s %-18s total %.3f/%.3f | slow %.3f/%.3f | ripple %.3f | band %.3f' % (fam, c, np.median(v[:, 0]), np.mean(v[:, 0]),
                  np.median(v[:, 1]), np.mean(v[:, 1]), np.median(v[:, 2]), np.median(v[:, 3])), flush=True)
