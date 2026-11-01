"""V3: temporal estimators in the measurement channel only (projection acceptance and feedback channel unchanged).
  E2-A        : accepted value from the E2 MHE (projection + k arrival prior + drift continuity), tuned in step 2
  E4-A        : accepted value from the E4 GLS drift-prior estimator (row3)
  E4arr-A     : E4 value fused with a k random-walk arrival prior (E2's q_k), scalar update each tick
Diagnosis cohort, 20 scenes per family. Run: python -X utf8 v3_temporal_estimators.py
"""
import math, os, sys, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, HERE)
import common, variants, cohort_eval as CE
import gls_estimator as G
import step2_estimators as S2
E2_CFG = (0.4919356373319487, 0.024402086755545314, 3.2928832250993345)
GLS_CFG = (0.12, 3.0)

class MeasChannel(variants.MotionAwareIMM):
    def __init__(self, mode, **kw):
        super().__init__(**kw); self.mode_ = mode
    def reset(self, platform):
        super().reset(platform)
        self.e2 = S2.MHE(E2_CFG); self.gls = G.GLS(GLS_CFG); self.ka = np.nan; self.va = np.nan
    def _value(self, o, n):
        if self.mode_ == 'E2-A':
            k, v = self.e2.step(S2.sanitize(o.log[:n]) if False else o.log, n); return float(k[0])
        f = self.gls.fit(o.log, n)
        if not (np.isfinite(f['k']) and f['info'] > 0): return np.nan
        k = float(np.clip(f['k'], 200., 3e5)); v = 1./f['info']
        if self.mode_ == 'E4-A': return k
        qk = E2_CFG[0]
        if not np.isfinite(self.ka): self.ka, self.va = k, v
        else:
            vp = self.va+(qk*self.ka)**2/S2.NW
            g = vp/(vp+v); self.ka = self.ka+g*(k-self.ka); self.va = (1-g)*vp
        return self.ka
    def _stiffness_est(self, o, n):
        kg_before = self.k_good
        out = super()._stiffness_est(o, n)
        val = self._value(o, n)
        if self.k_good != kg_before and np.isfinite(val): self.k_good = val
        return out

def build(cfg):
    if cfg == 'Proposed': return variants.build('Proposed-IMM', dict(CE.FINAL))
    return MeasChannel(cfg, single_model=False, **CE.FINAL)

def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(build(cfg), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(cfg=cfg, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'], loss_max=m['loss_max'])

if __name__ == '__main__':
    sc = CE.cohort('diag', 20); CFG = ('Proposed', 'E2-A', 'E4-A', 'E4arr-A')
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    json.dump(rows, open(os.path.join(HERE, 'results', 'v3_temporal_estimators.json'), 'w'), indent=1)
    for fam in CE.FAMS:
        print('%-10s' % fam, ' | '.join('%s %.3f/%.3f r%d L%d' % (c, np.median([r['rmse'] for r in rows if r['cfg'] == c and r['family'] == fam and r['ok']]),
              np.mean([r['rmse'] for r in rows if r['cfg'] == c and r['family'] == fam and r['ok']]),
              sum((not r['ok']) or r['peak'] > 50 for r in rows if r['cfg'] == c and r['family'] == fam),
              sum(r['loss_max'] > .5 for r in rows if r['cfg'] == c and r['family'] == fam)) for c in CFG), flush=True)
