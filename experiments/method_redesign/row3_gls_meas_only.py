"""Clean ablation of E4 (GLS drift-prior estimator): change only the measurement-gain channel.
Feedback channel (conservative k_bar), excitation, IMM, supervision: unchanged FINAL (projection).
  A: acceptance by the projection window's risk score b (unchanged rule); accepted value taken from E4.
  B: acceptance by the same risk score b computed on the E4 window (chi2, sigma from GLS); value from E4.
Diagnosis cohort (110000+), 20 scenes per family. Run: python -X utf8 row3_gls_meas_only.py
"""
import math, os, sys, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, HERE)
import common, variants, cohort_eval as CE
import gls_estimator as G
CFG_GLS = (0.12, 3.0)              # selected by row3_gls.py on the tuning seeds (ell = 120 ms, alpha = 3)

class MeasOnly(variants.MotionAwareIMM):
    def __init__(self, mode, **kw):
        super().__init__(**kw); self.mode_ = mode
    def reset(self, platform):
        super().reset(platform); self.gls = G.GLS(CFG_GLS)
    def _stiffness_est(self, o, n):
        kg_before = self.k_good
        out = super()._stiffness_est(o, n)           # projection: feedback k_bar and (possibly) k_good
        accepted_proj = self.k_good is not kg_before and self.k_good != kg_before
        f = self.gls.fit(o.log, n)
        if not (np.isfinite(f['k']) and f['info'] > 0):
            if self.mode_ == 'B': self.k_good = kg_before
            return out
        k = float(np.clip(f['k'], 200., 3e5)); sig = 1/math.sqrt(f['info'])
        if self.mode_ == 'A':
            if accepted_proj: self.k_good = k
        else:
            self.k_good = kg_before
            p = self.params; dof = max(p['nw']-p['degree']-2, 1)
            b = (sig/k)*(1.+math.sqrt(dof*max(0., f['chi2']-1.)))
            if sig <= p['rel_max']*k and f['chi2'] <= p['chi2_max'] and b <= math.log(1.2):
                self.k_good = k; self.f_set = max(float(o.log[n-1, 5]), .3); self.kg_rel = sig/k
        return out

def build(cfg):
    if cfg == 'Proposed': return variants.build('Proposed-IMM', dict(CE.FINAL))
    return MeasOnly(cfg[-1], single_model=False, **CE.FINAL)

def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(build(cfg), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(cfg=cfg, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'], loss_max=m['loss_max'])

if __name__ == '__main__':
    sc = CE.cohort('diag', 20); CFG = ('Proposed', 'E4-A', 'E4-B')
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    for fam in CE.FAMS:
        print('%-10s' % fam, ' | '.join('%s %.3f/%.3f (%d)' % (c, np.median([r['rmse'] for r in rows if r['cfg'] == c and r['family'] == fam and r['ok']]),
              np.mean([r['rmse'] for r in rows if r['cfg'] == c and r['family'] == fam and r['ok']]),
              sum((not r['ok']) or r['peak'] > 50 for r in rows if r['cfg'] == c and r['family'] == fam)) for c in CFG), flush=True)
    json.dump(rows, open(os.path.join(HERE, 'results', 'row3_gls_meas_only.json'), 'w'), indent=1)
