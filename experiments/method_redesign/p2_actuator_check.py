"""P2 fact check: is the mismatch degradation caused by the controller's wrong actuator/sensor model?
Same perturbed arm (servo lag x0.6, or 20 Hz force filter); controller either keeps the nominal model (as in P1)
or is told the true values. Diagnosis cohort, 20 scenes per family. Run: python -X utf8 p2_actuator_check.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE, robustness as RB
from common import P
COND = {'lag x0.6': dict(tau_a=0.019*0.6), 'filter 20 Hz': dict(fc=20.)}
def _ep(a):
    cond, informed, (fam, sid, p, T) = a
    plat = dict(RB.NOMINAL); plat.update(COND[cond])
    inner = variants.build('Proposed-IMM', dict(CE.FINAL))
    pol = inner if informed else RB.NominalModel(inner, seed=sid)
    log = P.episode(p, pol, T=T, platform=plat, seed=sid)
    if log is None: return cond, informed, fam, np.nan, np.inf
    m = common.task_metrics(log, 20.8)
    return cond, informed, fam, m['rmse'], m['peak']
if __name__ == '__main__':
    sc = CE.cohort('diag', 20)
    jobs = [(c, inf, s) for c in COND for inf in (False, True) for s in sc]
    with Pool(16) as po: rows = po.map(_ep, jobs, chunksize=1)
    for c in COND:
        for fam in CE.FAMS:
            out = []
            for inf in (False, True):
                r = [x for x in rows if x[0] == c and x[1] == inf and x[2] == fam]
                v = [x[3] for x in r if np.isfinite(x[3])]; bad = sum((not np.isfinite(x[3])) or x[4] > 50 for x in r)
                out.append('%s %.3f/%.3f r%d' % ('told' if inf else 'nominal', np.median(v), np.mean(v), bad))
            print('%-12s %-10s %s' % (c, fam, ' | '.join(out)), flush=True)
