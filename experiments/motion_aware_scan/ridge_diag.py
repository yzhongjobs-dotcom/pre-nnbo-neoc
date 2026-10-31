"""Diagnosis set for high ridge ratios (30-60x; seeds 180000+, separate from the sweep and the confirmation cohort)."""
import json, os, sys
import numpy as np
from multiprocessing import Pool
import common, cohort_eval as CE, scenes_v2 as SV
def scenes(n=30):
    rng = np.random.default_rng(180000); out = []
    for i in range(n):
        p, T = SV.S.make('T1', rng, dict(A=0.0, ratio=float(rng.uniform(30, 60)))); out.append(('ridge-hi', 180000+i, p, T))
    return out
def _ep(a):
    m, (fam, sid, p, T) = a
    log, _ = common.run(CE.build(m, p), p, T, sid); mt = common.task_metrics(log, 20.8)
    return dict(cfg=m, family=fam, sid=sid, ok=not mt['diverged'], rmse=mt['rmse'], peak=mt['peak'])
if __name__ == '__main__':
    cfgs = sys.argv[1].split(',')
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in cfgs for s in scenes()])
    for c in cfgs:
        R = [r for r in rows if r['cfg'] == c]; ok = [r['rmse'] for r in R if r['ok']]
        print('ridge-hi', c, 'median %.3f mean %.3f bad %d' % (np.median(ok), np.mean(ok), sum((not r['ok']) or r['peak'] > 50 for r in R)))
