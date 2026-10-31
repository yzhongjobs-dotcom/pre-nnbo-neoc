"""How accurate is the measurement stiffness used by the motion channel? (diagnosis for property 4)

For each configuration, logs the gain actually used to convert force into motion (kref inside drift_pred) and the true
tangent stiffness, and reports per family: the fraction of in-contact ticks with kappa_g = k/kref outside the frozen
band [0.75, 1.2], the 95th percentile of |ln kappa_g|, and the same restricted to push/retract ramps.
Run: CONFIGS=final,F1,F8 COHORT=diag N=30 FAMS=T2m-stiff,T2f-stiff python -X utf8 diag_kappag.py
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import cohort_eval as CE
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']


class Logged(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.kref_log = {}

    def drift_pred(self, o, n):
        r = super().drift_pred(o, n)
        if self._fin_prev is not None and self.tracker.initialized:
            self.kref_log[n] = self._kref_used
        return r


_orig = variants.MotionAwareIMM.drift_pred


def _patched(self, o, n):
    # record the conversion gain: recompute it exactly as drift_pred does
    j = n-1
    F = o.log[j, 5]
    kref = self.k_good
    g = self.params.get('kg_gamma', 0.)
    if g and np.isfinite(kref) and np.isfinite(getattr(self, 'f_set', np.nan)):
        kref = kref*(max(F, .3)/self.f_set)**g
    if self.params.get('oracle_meas'):
        kref = max(float(o.log[j, 8]), 1.)
    if not np.isfinite(kref) or kref <= 0:
        kref = self.kbar
    self._kref_used = kref
    return _orig(self, o, n)


Logged.drift_pred.__wrapped__ = None
variants.MotionAwareIMM.drift_pred = _patched


def ramp_mask(t, p):
    m = np.zeros(len(t), bool)
    for k in range(int(p[15])):
        t0, tr, hold = p[16+4*k], p[18+4*k], p[19+4*k]
        for a in (t0, t0+tr+hold):
            m |= (t >= a) & (t <= a+tr+0.1)
    return m


def _ep(a):
    cfg, (fam, sid, p, T) = a
    extra = dict(CE.extra_configs()).get(cfg, {})
    pol = Logged(single_model=False, **dict(BASE, **extra))
    log, pol = common.run(pol, p, T, sid)
    t, F, k = log[:, 0], log[:, 3], log[:, 8]
    idx = np.array(sorted(pol.kref_log)); idx = idx[idx < len(t)]
    kr = np.array([pol.kref_log[i] for i in idx])
    valid = (t[idx] > 1) & (F[idx] > 0.5) & (k[idx] > 0)
    lk = np.log(k[idx][valid]/kr[valid]); rm = ramp_mask(t[idx], p)[valid]
    out = lambda x: float(np.mean((x > math.log(1.2)) | (x < math.log(0.75)))) if len(x) else np.nan
    q95 = lambda x: float(np.percentile(np.abs(x), 95)) if len(x) else np.nan
    return dict(cfg=cfg, family=fam, out_all=out(lk), out_ramp=out(lk[rm]), q95_all=q95(lk), q95_ramp=q95(lk[rm]),
                bias_ramp=float(np.median(lk[rm])) if rm.any() else np.nan, peak=float(F[t > 1].max()))


if __name__ == '__main__':
    cfgs = os.environ.get('CONFIGS', 'final').split(',')
    fams = os.environ.get('FAMS', 'T2m-stiff,T2f-stiff').split(',')
    sc = [s for s in CE.cohort(os.environ.get('COHORT', 'diag'), int(os.environ.get('N', '30'))) if s[0] in fams]
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in cfgs for s in sc])
    for fam in fams:
        for c in cfgs:
            R = [r for r in rows if r['cfg'] == c and r['family'] == fam]
            print('%-9s %-6s outside band: all %.3f ramps %.3f | q95 |ln kappa_g| all %.2f ramps %.2f | median ln kappa_g on ramps %+.2f | peak>50N %d'
                  % (fam, c, np.nanmean([r['out_all'] for r in R]), np.nanmean([r['out_ramp'] for r in R]),
                     np.nanmean([r['q95_all'] for r in R]), np.nanmean([r['q95_ramp'] for r in R]),
                     np.nanmean([r['bias_ramp'] for r in R]), sum(r['peak'] > 50 for r in R)), flush=True)
