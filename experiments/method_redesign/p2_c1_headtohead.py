"""P2 fact check: original C1 (per-cycle Gauss-Newton receding-horizon optimisation of D, K_PF and u_ff, repaired
projection estimator; cobot_scan_benchmark/contact_update.py, default parameters tuned on the earlier three-task
protocol) versus the final method, same diagnosis scenes (20 per family). Run: python -X utf8 p2_c1_headtohead.py
"""
import os, sys, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
from common import P
from contact_update import ContactUpdate
def _ep(a):
    cfg, (fam, sid, p, T) = a
    pol = variants.build('Proposed-IMM', dict(CE.FINAL)) if cfg == 'Proposed' else ContactUpdate()
    t0 = time.perf_counter()
    log = P.episode(p, pol, T=T, seed=sid)
    if log is None: return cfg, fam, np.nan, np.inf, 0
    m = common.task_metrics(log, 20.8)
    return cfg, fam, m['rmse'], m['peak'], m['loss_max']
if __name__ == '__main__':
    sc = CE.cohort('diag', 20)
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in ('Proposed', 'C1-GN') for s in sc], chunksize=1)
    for fam in CE.FAMS:
        out = []
        for c in ('Proposed', 'C1-GN'):
            r = [x for x in rows if x[0] == c and x[1] == fam]; v = [x[2] for x in r if np.isfinite(x[2])]
            out.append('%s %.3f/%.3f r%d L%d' % (c, np.median(v), np.mean(v), sum((not np.isfinite(x[2])) or x[3] > 50 for x in r), sum(x[4] > .5 for x in r)))
        print('%-10s %s' % (fam, ' | '.join(out)), flush=True)
