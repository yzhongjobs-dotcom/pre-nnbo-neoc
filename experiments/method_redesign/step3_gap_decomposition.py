"""Step 3 diagnostic: where does the gap between the proposed method and the motion oracle come from?
Same control structure in every row (stiffness-normalised admittance + motion feed-forward):
  Proposed                : estimated stiffness, IMM motion prediction, excitation on
  TrueK                   : true stiffness for feedback and measurement gain, IMM prediction, excitation on
  Oracle+exc              : true stiffness and true short-horizon surface prediction, excitation on
  Oracle                  : as above, excitation off (existing reference controller)
Diagnosis cohort (seeds 110000+), motion families, 20 scenes each. Run: python -X utf8 step3_gap_decomposition.py
"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan')
sys.path.insert(0, H)
import common, variants, cohort_eval as CE
import proposed as Pm

FAMS = ('T2m-soft', 'T2m-stiff', 'T2f-stiff', 'T3m')
CFG = ('Proposed', 'TrueK', 'Oracle+exc', 'Oracle')

def build(cfg, p):
    F = dict(CE.FINAL)
    if cfg == 'Proposed':
        return variants.build('Proposed-IMM', F)
    if cfg == 'TrueK':
        return variants.build('Proposed-IMM', dict(F, oracle_fb=True, oracle_meas=True))
    kw = {k: F[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
    if cfg == 'Oracle+exc':
        return Pm.Oracle(p, excitation=True, a_force=F.get('a_force', .17), f_exc=F.get('f_exc', 20.8), **kw)
    return Pm.Oracle(p, **kw)

def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(build(cfg, p), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(cfg=cfg, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'])

if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in FAMS]
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    out = {}
    for fam in FAMS:
        line = '%-10s' % fam
        for c in CFG:
            r = [x['rmse'] for x in rows if x['cfg'] == c and x['family'] == fam and x['ok']]
            bad = sum((not x['ok']) or x['peak'] > 50 for x in rows if x['cfg'] == c and x['family'] == fam)
            out[(fam, c)] = (float(np.median(r)), float(np.mean(r)), bad)
            line += ' | %-10s %.3f/%.3f (%d)' % (c, np.median(r), np.mean(r), bad)
        print(line, flush=True)
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'step3_gap_decomposition.json'), 'w'), indent=1)
