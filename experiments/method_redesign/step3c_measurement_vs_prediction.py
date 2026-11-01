"""Step 3 diagnostic 2: is the remaining gap prediction (horizon) or measurement content?
All rows: true stiffness, excitation on, same admittance + feed-forward structure.
  TrueK        : feed-forward from the IMM on the force-derived motion measurement (as deployed)
  TrueS-now    : feed-forward from the true surface at the current time (no prediction)
  TrueS-12     : true surface 12 ms ahead, no servo-lag term
  TrueS-pred   : true surface 12 ms ahead + tau_s * velocity (= Oracle+exc)
Diagnosis cohort (110000+), 20 scenes per motion family. Run: python -X utf8 step3c_measurement_vs_prediction.py
"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan')
sys.path.insert(0, H); sys.path.insert(0, os.path.join(H, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
import proposed as Pm
import plant as P

FAMS = ('T2m-soft', 'T2m-stiff', 'T2f-stiff')
CFG = ('TrueK', 'TrueS-now', 'TrueS-12', 'TrueS-pred')

class TrueSurface(Pm.Oracle):
    def __init__(self, p_scene, lead=0., use_vel=False, **kw):
        super().__init__(p_scene, **kw); self.lead = lead; self.use_vel = use_vel
    def drift_pred(self, o, n):
        x, t = o.log[n-1, 1], o.log[n-1, 0]
        hh, ht, _ = P.surface(x, t+self.lead, self.p_scene)
        q = hh+(self.tau*ht if self.use_vel else 0.)
        return self.kbar*q, q

def build(cfg, p):
    F = dict(CE.FINAL)
    if cfg == 'TrueK':
        return variants.build('Proposed-IMM', dict(F, oracle_fb=True, oracle_meas=True))
    kw = {k: F[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
    kw.update(excitation=True, a_force=F['a_force'], f_exc=20.8)
    lead = {'TrueS-now': 0., 'TrueS-12': 0.012, 'TrueS-pred': 0.012}[cfg]
    return TrueSurface(p, lead=lead, use_vel=(cfg == 'TrueS-pred'), **kw)

def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(build(cfg, p), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(cfg=cfg, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'])

if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in FAMS]
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    for fam in FAMS:
        line = '%-10s' % fam
        for c in CFG:
            r = [x['rmse'] for x in rows if x['cfg'] == c and x['family'] == fam and x['ok']]
            line += ' | %-10s %.3f/%.3f' % (c, np.median(r), np.mean(r))
        print(line, flush=True)
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'step3c_measurement_vs_prediction.json'), 'w'), indent=1)
