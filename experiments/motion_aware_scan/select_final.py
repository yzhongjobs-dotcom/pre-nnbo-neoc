"""Pre-declared final-method selection on an independent selection set (10 scenes x 5 families,
seeds 80000+, never used for tuning or testing). Criterion declared before running: lowest mean
log(1+J); if within 5% of the best, prefer the simpler candidate in the order KF < KFg < IMM < IMMg.
Run: python -X utf8 select_final.py -> results_v2/selection.json"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common, variants, scenes_v2 as SV

CANDS = ['Proposed-KF', 'Proposed-KFg', 'Proposed-IMM', 'Proposed-IMMg']


def tuned(m):
    return json.load(open(os.path.join(common.HERE, 'tuned_v2', m+'.json')))['best']['params']


def ev(a):
    name, (fam, sid, p, T) = a
    log, _ = common.run(variants.build(name, tuned(name)), p, T, sid)
    m = common.task_metrics(log, 20.8)
    J = m['J'] if np.isfinite(m['J']) else 1e3
    return dict(name=name, fam=fam, sid=sid, L=math.log1p(min(J, 1e3)), rmse=m['rmse'], peak=m['peak'])


if __name__ == '__main__':
    sc = SV.scenes('select')
    with Pool(16) as pool:
        res = pool.map(ev, [(n, s) for n in CANDS for s in sc], chunksize=1)
    out = {}
    for n in CANDS:
        r = [x for x in res if x['name'] == n]
        out[n] = dict(mean_L=float(np.mean([x['L'] for x in r])), max_L=float(np.max([x['L'] for x in r])),
                      n_peak_over_50N=int(sum(x['peak'] > 50 for x in r)),
                      by_family={f: float(np.mean([x['L'] for x in r if x['fam'] == f])) for f in SV.FAMILIES})
    best = min(out, key=lambda k: out[k]['mean_L'])
    final = next(c for c in CANDS if out[c]['mean_L'] <= 1.05*out[best]['mean_L'])
    json.dump(dict(criterion='min mean log(1+J); within 5% prefer simpler (KF<KFg<IMM<IMMg)', results=out,
                   best=best, final=final, rows=res), open(os.path.join(common.HERE, 'results_v2', 'selection.json'), 'w'), indent=1)
    for k, v in out.items():
        print(k, round(v['mean_L'], 3), 'max', round(v['max_L'], 2), 'peaks>50N', v['n_peak_over_50N'], {f: round(x, 2) for f, x in v['by_family'].items()})
    print('best', best, 'final', final)
