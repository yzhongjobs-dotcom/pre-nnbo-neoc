"""Where does the true-stiffness gain in pushes come from? Split TrueK into measurement gain vs feedback stiffness,
and test removing feedback conservatism. Diagnosis cohort (110000+), 20 scenes per family.
Run: python -X utf8 row3_fastpush_split.py
"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE
FAMS = ('T2m-stiff', 'T2f-stiff')
CFG = {'Proposed': {}, 'TrueMeas': dict(oracle_meas=True), 'TrueFb': dict(oracle_fb=True),
       'TrueK': dict(oracle_meas=True, oracle_fb=True), 'NoConservative': dict(conservative=False)}
def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(variants.build('Proposed-IMM', dict(CE.FINAL, **CFG[cfg])), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(cfg=cfg, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'])
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in FAMS]
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    for fam in FAMS:
        print(fam, ' | '.join('%s %.3f/%.3f (%d)' % (c, np.median([x['rmse'] for x in rows if x['cfg'] == c and x['family'] == fam and x['ok']]),
              np.mean([x['rmse'] for x in rows if x['cfg'] == c and x['family'] == fam and x['ok']]),
              sum((not x['ok']) or x['peak'] > 50 for x in rows if x['cfg'] == c and x['family'] == fam)) for c in CFG), flush=True)
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'row3_fastpush_split.json'), 'w'), indent=1)
