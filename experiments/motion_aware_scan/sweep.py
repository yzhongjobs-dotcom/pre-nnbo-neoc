"""Dose-response sweeps over the variables that create the difficulty (paper draft §五).

(a) Push speed on stiff tissue: ramp time 0.6 -> 0.1 s (amplitude 3-8 mm as in the benchmark), breathing on.
(b) Stiffness span on the static scan: ridge/base stiffness ratio 3 -> 60.
20 fresh scenes per dose (seeds 160000+ / 170000+), never used for tuning, selection or diagnosis.
Run: python -X utf8 sweep.py -> results_v2/sweep.json
"""
import json, os
import numpy as np
from multiprocessing import Pool
import common
import cohort_eval as CE
import scenes_v2 as SV

METHODS = ('Proposed', 'RossiType', 'WijType', 'DOMPAC', 'YuenType')
RAMPS = (0.6, 0.45, 0.3, 0.2, 0.15, 0.1)
RATIOS = (3., 10., 20., 30., 45., 60.)
N = 20


def scenes():
    out = []
    for i, r in enumerate(RAMPS):
        rng = np.random.default_rng(160000+100*i)
        for j in range(N):
            p, T = SV.S.make('T2', rng, dict(stiff=True, ramp=r)); out.append(('ramp', r, 160000+100*i+j, p, T))
    for i, q in enumerate(RATIOS):
        rng = np.random.default_rng(170000+100*i)
        for j in range(N):
            p, T = SV.S.make('T1', rng, dict(A=0.0, ratio=q)); out.append(('ratio', q, 170000+100*i+j, p, T))
    return out


def _ep(a):
    m, (kind, dose, sid, p, T) = a
    try:
        log, _ = common.run(CE.build(m, p), p, T, sid)
        mt = common.task_metrics(log, 20.8)
        return dict(cfg=m, kind=kind, dose=dose, sid=sid, ok=not mt['diverged'], rmse=mt['rmse'], peak=mt['peak'],
                    loss=mt['loss'], loss_max=mt['loss_max'])
    except Exception as ex:
        return dict(cfg=m, kind=kind, dose=dose, sid=sid, ok=False, err=str(ex))


if __name__ == '__main__':
    sc = scenes()
    with Pool(16) as po:
        rows = po.map(_ep, [(m, s) for m in METHODS for s in sc], chunksize=1)
    json.dump(rows, open(os.path.join(common.HERE, 'results_v2', 'sweep.json'), 'w'), indent=1)
    for kind, doses in (('ramp', RAMPS), ('ratio', RATIOS)):
        for d in doses:
            line = '%-5s %5.2f' % (kind, d)
            for m in METHODS:
                R = [r for r in rows if r['cfg'] == m and r['kind'] == kind and r['dose'] == d]
                ok = [r['rmse'] for r in R if r['ok']]
                line += ' | %s %.3f/%.3f (%d)' % (m[:4], np.median(ok), np.mean(ok), sum((not r['ok']) or r.get('peak', 0) > 50 for r in R))
            print(line, flush=True)
