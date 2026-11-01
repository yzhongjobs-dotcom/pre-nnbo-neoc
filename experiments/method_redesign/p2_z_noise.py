"""Is the motion-measurement noise model of the IMM calibrated? Records Z and R per tick in the final method
(diagnosis cohort, 10 per family), reports (a) empirical variance of the high-passed Z (Z minus its 5-tick centred
mean, scaled to white-noise variance) over R/r_scale (=(sigma_F/k)^2) and over the tracker's R, (b) the share of Z
high-band power at 20.8 Hz (excitation leak), (c) the normalised innovation squared of the IMM mixture prediction.
Run: python -X utf8 p2_z_noise.py
"""
import os, sys
import numpy as np
from scipy.signal import welch
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE
F = dict(CE.FINAL)
class Rec(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform); self.rec = []
    def drift_pred(self, o, n):
        pre = (self.tracker.x.copy(), self.tracker.P.copy()) if self.tracker.initialized else None
        r = super().drift_pred(o, n)
        if self.tracker.initialized and pre is not None and self._fin_prev is not None:
            A = self.tracker.A; xp = A@pre[0]; Sp = (A@pre[1]@A.T)[0, 0]+self.tracker.R
            self.rec.append((n, self.Z, self.tracker.R, (self.Z-xp[0])**2/Sp))
        return r
def _ep(a):
    fam, sid, p, T = a
    pol = variants.build('Proposed-IMM', F); pol.__class__ = Rec
    log, _ = common.run(pol, p, T, sid)
    R = np.array(pol.rec)
    if len(R) < 300: return fam, None
    R = R[R[:, 0] > 125]; z = R[:, 1]; w = 5
    # 5-tick centred mean removes slow motion; residual variance of white noise is (1-1/w) sigma^2
    hp = z[w//2:-(w//2)]-np.convolve(z, np.ones(w)/w, mode='valid')
    v = np.var(hp)/(1-1/w); Rm = np.median(R[:, 2])
    f, S = welch(z-np.convolve(z, np.ones(25)/25, mode='same'), fs=125., nperseg=256)
    hb = (f > 10); line = (f > 19.5) & (f < 22.)
    return fam, (v/(Rm/F['r_scale']), v/Rm, float(S[line].sum()/S[hb].sum()), float(np.median(R[:, 3])), float(np.mean(R[:, 3])))
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 10)]
    with Pool(16) as po: rows = po.map(_ep, sc, chunksize=1)
    print('family | var(HP Z)/(sigma_F/k)^2 | var/R_tracker | 20.8 Hz share of >10 Hz power | NIS median/mean')
    for fam in CE.FAMS:
        v = np.array([r[1] for r in rows if r[0] == fam and r[1] is not None])
        print('%-9s %.2f | %.1f | %.2f | %.2f/%.1f' % (fam, *np.median(v, 0)), flush=True)
