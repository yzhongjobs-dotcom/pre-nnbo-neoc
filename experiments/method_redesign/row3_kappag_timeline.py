"""Fast pushes: when is the measurement gain out of the feed-forward band, and why?
Logs kappa_g = k_true / k_good per tick for the final method on diagnosis T2f-stiff and T2m-stiff scenes.
Run: python -X utf8 row3_kappag_timeline.py
"""
import os, sys
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE

class Rec(variants.MotionAwareIMM):
    def step(self, o):
        out = super().step(o)
        n = o.n
        if not hasattr(self, 'rec'): self.rec = []; self._last_upd = 0.; self._kg_prev = np.nan
        t = o.log[n-1, 0]
        if self.k_good != self._kg_prev: self._last_upd = t; self._kg_prev = self.k_good
        self.rec.append((t, o.log[n-1, 8], self.k_good, o.log[n-1, 5], t-self._last_upd, self.kbar))
        return out

def _ep(sc):
    fam, sid, p, T = sc
    pol = Rec(single_model=False, **CE.FINAL)
    log, pol = common.run(pol, p, T, sid)
    R = np.array(pol.rec); R = R[R[:, 0] > 1.]
    return fam, R

if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('T2f-stiff', 'T2m-stiff')]
    with Pool(16) as po: res = po.map(_ep, sc)
    for fam in ('T2m-stiff', 'T2f-stiff'):
        R = np.vstack([r for f, r in res if f == fam])
        ok = (R[:, 1] > 0) & np.isfinite(R[:, 2])
        kg = R[ok, 1]/R[ok, 2]; F = R[ok, 3]; age = R[ok, 4]
        out = (kg > 1.2) | (kg < 1/1.2)
        print(fam, 'ticks %d, out of band %.1f%% (kappa_g>1.2: %.1f%%, <0.83: %.1f%%)' % (len(kg), 100*out.mean(), 100*(kg > 1.2).mean(), 100*(kg < 1/1.2).mean()))
        for lo, hi in ((0, .1), (.1, .3), (.3, 1.), (1., 99)):
            m = (age >= lo) & (age < hi)
            if m.sum(): print('   age since last accepted update in [%.1f,%.1f) s: share %.1f%%, out of band %.1f%%, median |ln kg| %.3f' % (lo, hi, 100*m.mean(), 100*out[m].mean(), np.median(np.abs(np.log(kg[m])))))
        for lo, hi in ((0, 3), (3, 6), (6, 12), (12, 1e9)):
            m = (F >= lo) & (F < hi)
            if m.sum(): print('   force in [%g,%g) N: share %.1f%%, out of band %.1f%%' % (lo, hi, 100*m.mean(), 100*out[m].mean()))
