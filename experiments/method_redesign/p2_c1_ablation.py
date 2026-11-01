"""Which part of the original C1 gives its advantage on stiffness-only scans?
C1 variants: default; iters=0 (no GN optimisation, gains from the stiffness-scaled prior rule); reference=False
(no u_ff reference shaping); ki=0 (no integral force action). Diagnosis cohort, 20 scenes per family.
Run: python -X utf8 p2_c1_ablation.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, cohort_eval as CE
from common import P
from contact_update import ContactUpdate
CFG = {'C1': {}, 'C1 no-GN': dict(iters=0), 'C1 no-ref': dict(reference=False), 'C1 ki=0': dict(ki=0.)}
def _ep(a):
    cfg, (fam, sid, p, T) = a
    try:
        log = P.episode(p, ContactUpdate(**CFG[cfg]), T=T, seed=sid)
        if log is None: return cfg, fam, np.nan, np.inf
        m = common.task_metrics(log, 20.8); return cfg, fam, m['rmse'], m['peak']
    except Exception as ex:
        return cfg, fam, np.nan, repr(ex)
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('S2', 'T3m', 'T2m-soft')]
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    errs = [r for r in rows if isinstance(r[3], str)]
    if errs: print('errors:', errs[:2])
    for fam in ('S2', 'T3m', 'T2m-soft'):
        print('%-9s' % fam, ' | '.join('%s %.3f/%.3f' % (c, np.nanmedian([r[2] for r in rows if r[0] == c and r[1] == fam]), np.nanmean([r[2] for r in rows if r[0] == c and r[1] == fam])) for c in CFG))
