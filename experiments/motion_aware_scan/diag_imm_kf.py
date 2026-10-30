"""Dev-set robustness of IMM vs KF tracker and the failing test T3m scenes."""
import json, math, numpy as np, common, variants, scenes_v2 as SV
from multiprocessing import Pool

def tuned(m): return json.load(open('tuned_v2/%s.json' % m))['best']['params']

def ev(a):
    name, (fam, sid, p, T) = a
    log, _ = common.run(variants.build(name, tuned(name)), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return (name, fam, sid, round(math.log1p(min(m['J'], 1e3)), 2), round(m['rmse'], 3), round(m['peak'], 1))

if __name__ == '__main__':
    sc = SV.scenes('dev')
    with Pool(16) as pool:
        res = pool.map(ev, [(n, s) for n in ('Proposed-IMM', 'Proposed-KF') for s in sc])
    for n in ('Proposed-IMM', 'Proposed-KF'):
        r = [x for x in res if x[0] == n]
        L = [x[3] for x in r]
        print(n, 'mean L', round(np.mean(L), 3), 'max L', max(L), 'worst', sorted(r, key=lambda x: -x[3])[:3])
    t = json.load(open('results_v2/test/main.json'))
    print('IMM T3m bad test scenes', [(r['sid'], round(r['rmse'], 1), round(r['peak'])) for r in t
                                      if r['method'] == 'Proposed' and r['family'] == 'T3m' and r['rmse'] > 2])
