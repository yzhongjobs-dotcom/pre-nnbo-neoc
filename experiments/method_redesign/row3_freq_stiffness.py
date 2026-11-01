"""Hypothesis check: does the force-to-motion conversion need the low-frequency (relaxed) stiffness rather than the
instantaneous one measured by the 20.8 Hz excitation? Measurement gain = factor x true instantaneous tangent stiffness.
Diagnosis cohort, stiff push families, 20 scenes each. Run: python -X utf8 row3_freq_stiffness.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE

class ScaledTrueMeas(variants.MotionAwareIMM):
    def __init__(self, factor, **kw):
        super().__init__(**kw); self.factor = factor
    def step(self, o):
        # true instantaneous tangent stiffness (log col 8) scaled; diagnosis only
        self.k_good = max(float(o.log[o.n-1, 8]), 1.)*self.factor
        return super().step(o)

def _ep(a):
    fac, (fam, sid, p, T) = a
    pol = ScaledTrueMeas(fac, single_model=False, **dict(CE.FINAL, oracle_meas=False))
    log, _ = common.run(pol, p, T, sid)
    m = common.task_metrics(log, 20.8)
    return fac, fam, m['rmse'], m['peak']

if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('T2m-stiff', 'T2f-stiff')]
    facs = (1.0, 0.85, 1/1.4)
    with Pool(16) as po: rows = po.map(_ep, [(f, s) for f in facs for s in sc], chunksize=1)
    for fam in ('T2m-stiff', 'T2f-stiff'):
        print(fam, ' | '.join('x%.2f %.3f/%.3f' % (f, np.median([r[2] for r in rows if r[0] == f and r[1] == fam]), np.mean([r[2] for r in rows if r[0] == f and r[1] == fam])) for f in facs))
