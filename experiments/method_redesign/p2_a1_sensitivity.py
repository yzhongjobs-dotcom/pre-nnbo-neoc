"""Is the stiffness-scan performance sensitive to the choice of the unit-stiffness admittance A1 inside the certified set?
Fixed A1 designs from motion_aware_scan/A1_feasible.json (535 designs passing frozen stability on [0.5, 2.14] and the
circle criterion): the final one and designs spanning the bandwidth range; everything else = final method.
Diagnosis cohort, 20 scenes per family (S2, T3m, T2m-soft, T2m-stiff). Run: python -X utf8 p2_a1_sensitivity.py
"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
MA = os.path.join(HERE, '..', 'motion_aware_scan'); sys.path.insert(0, MA)
import common, variants, cohort_eval as CE
FEAS = json.load(open(os.path.join(MA, 'A1_feasible.json')))
bw = np.array([d['bw'] for d in FEAS]); order = np.argsort(bw)
pick = [FEAS[order[int(q*(len(order)-1))]] for q in (0., .25, .5, .75, .9, 1.)]
DES = {'final': {k: CE.FINAL[k] for k in ('M1', 'D1', 'KPF1', 'ki')}}
for d in pick: DES['bw%.0f' % d['bw']] = {k: d[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
def _ep(a):
    name, (fam, sid, p, T) = a
    log, _ = common.run(variants.build('Proposed-IMM', dict(CE.FINAL, **DES[name])), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return name, fam, m['rmse'], m['peak'], not m['diverged']
if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 20) if s[0] in ('S2', 'T3m', 'T2m-soft', 'T2m-stiff')]
    with Pool(16) as po: rows = po.map(_ep, [(n, s) for n in DES for s in sc], chunksize=1)
    print({n: DES[n] for n in DES})
    for fam in ('S2', 'T3m', 'T2m-soft', 'T2m-stiff'):
        print('%-9s' % fam, ' | '.join('%s %.3f/%.3f r%d' % (n, np.median([r[2] for r in rows if r[0] == n and r[1] == fam and r[4]]),
              np.mean([r[2] for r in rows if r[0] == n and r[1] == fam and r[4]]), sum((not r[4]) or r[3] > 50 for r in rows if r[0] == n and r[1] == fam)) for n in DES))
