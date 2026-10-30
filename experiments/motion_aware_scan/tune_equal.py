"""Equal-budget tuning on development scenes (methods §10): every method gets the same
random-search budget over its own free parameters (plus its default), scored by the
benchmark's physical composite J (force error /0.5 N, contact loss, peak excess) averaged
over the development cohort. Test cohorts use different seeds and are never touched here.
Run: python -X utf8 tune_equal.py [method ...]  -> tuned/<method>.json
"""
import json, os, sys, math, time
import numpy as np
from multiprocessing import Pool
import common
from common import S

N_RANDOM = 40
DEV = [('T1', 6, 9100), ('T2', 6, 9200), ('T3', 6, 9300)]
OUT = os.path.join(common.HERE, 'tuned'); os.makedirs(OUT, exist_ok=True)
FEAS = json.load(open(os.path.join(common.HERE, 'A1_feasible.json')))


def logu(r, lo, hi): return float(math.exp(r.uniform(math.log(lo), math.log(hi))))


def space(method, r):
    if method == 'Proposed':
        d = FEAS[int(r.integers(0, min(len(FEAS), 120)))]          # circle-certified designs only
        return dict(M1=d['M1'], D1=d['D1'], KPF1=d['KPF1'], ki=0.0, a_force=logu(r, .08, .6),
                    c=float(r.uniform(1., 3.)), rho_d=logu(r, .005, .06), nw=int(r.choice([12, 16, 24])),
                    chi2_max=logu(r, 2., 10.), rel_max=float(r.uniform(.15, .4)))
    if method == 'RossiType':
        return dict(K_I=logu(r, 2., 120.), lam=float(r.uniform(.95, .999)), Q=logu(r, 1., 8.),
                    A_E=logu(r, 2e-6, 1e-4))
    if method == 'WijType':
        return dict(K_P=float(r.uniform(0., 1.5)), K_I=logu(r, 2., 120.), innovation_scale=logu(r, .2, 5.))
    if method == 'DOMPAC':
        return dict(m=logu(r, .3, 10.), b=logu(r, 50., 5000.), k_h0=logu(r, 500., 5e4),
                    eso_bandwidth=logu(r, 2., 60.), Q=logu(r, .1, 100.), R=logu(r, .01, 10.))
    if method == 'RossiDOMPAC':
        return dict(m=logu(r, .3, 10.), b=logu(r, 50., 5000.), eso_bandwidth=logu(r, 2., 60.),
                    Q=logu(r, .1, 100.), R=logu(r, .01, 10.), lam=float(r.uniform(.95, .999)),
                    bp_Q=logu(r, 1., 8.), A_E=logu(r, 2e-6, 1e-4))
    raise ValueError(method)


def make(method, params):
    import comparators as C, proposed as Pm
    if method == 'Proposed':
        return Pm.MotionAware(**params)
    return getattr(C, method)(**params)


def scenes(spec):
    out = []
    for task, n, seed in spec:
        for i, (p, T) in enumerate(S.cohort(task, n, seed)):
            out.append((task, seed+i, p, T))
    return out


def _eval(args):
    method, params, (task, sid, p, T) = args
    try:
        log, _ = common.run(make(method, params), p, T, sid)
        m = common.task_metrics(log, 20.8)
    except Exception as ex:                      # invalid parameter combination
        return dict(J=1e3, rmse=np.nan, err=str(ex))
    return dict(J=min(m['J'], 1e3), rmse=m['rmse'])


def tune(method, pool):
    r = np.random.default_rng(12345)
    cands = [{}]+[space(method, r) for _ in range(N_RANDOM)]
    sc = scenes(DEV)
    jobs = [(method, c, s) for c in cands for s in sc]
    res = pool.map(_eval, jobs, chunksize=1)
    k = len(sc); rows = []
    for i, c in enumerate(cands):
        Js = [x['J'] for x in res[i*k:(i+1)*k]]
        rows.append(dict(params=c, J=float(np.mean(Js)), rmse=float(np.nanmean([x['rmse'] for x in res[i*k:(i+1)*k]]))))
    best = min(rows, key=lambda x: x['J'])
    json.dump(dict(method=method, budget=len(cands), dev=DEV, best=best, all=rows), open(os.path.join(OUT, method+'.json'), 'w'), indent=1)
    return best


if __name__ == '__main__' and os.environ.get('TUNE_STAGE') != '2':
    methods = sys.argv[1:] or ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC']
    with Pool(16) as pool:
        for mth in methods:
            t0 = time.time(); b = tune(mth, pool)
            print(f"{mth}: J={b['J']:.3f} rmse={b['rmse']:.3f} params={b['params']} ({time.time()-t0:.0f}s)", flush=True)


N_LOCAL = 20


def perturb(method, params, r):
    """stage 2: local search around the stage-1 best (same budget for every method)."""
    q = dict(params)
    if method == 'Proposed':
        d0 = q['D1']
        near = [d for d in FEAS[:120] if 0.6 <= d['D1']/d0 <= 1.7]
        d = near[int(r.integers(0, len(near)))] if near else dict(M1=q['M1'], D1=q['D1'], KPF1=q['KPF1'])
        q.update(M1=d['M1'], D1=d['D1'], KPF1=d['KPF1'])
        for k in ('a_force', 'c', 'rho_d', 'chi2_max', 'rel_max'):
            q[k] = float(q[k]*math.exp(r.normal(0, .25)))
        q['nw'] = int(r.choice([16, 24, 32])) if r.random() < .3 else q['nw']
        return q
    for k, v in params.items():
        if isinstance(v, float):
            q[k] = float(v*math.exp(r.normal(0, .3)))
    if 'lam' in q: q['lam'] = float(min(q['lam'], .9995))
    return q


def refine(method, pool):
    path = os.path.join(OUT, method+'.json')
    rec = json.load(open(path))
    base = rec['best']['params'] or space(method, np.random.default_rng(0)) if rec['best']['params'] else {}
    if not rec['best']['params']:
        # default won stage 1: perturb around the class defaults
        import comparators as C, proposed as Pm
        base = dict(getattr(C, method).params) if method != 'Proposed' else Pm.MotionAware().params
    r = np.random.default_rng(54321)
    cands = [perturb(method, base, r) for _ in range(N_LOCAL)]
    sc = scenes(DEV)
    res = pool.map(_eval, [(method, c, s) for c in cands for s in sc], chunksize=1)
    k = len(sc); rows = rec['all']
    for i, c in enumerate(cands):
        Js = [x['J'] for x in res[i*k:(i+1)*k]]
        rows.append(dict(params=c, J=float(np.mean(Js)), rmse=float(np.nanmean([x['rmse'] for x in res[i*k:(i+1)*k]]))))
    best = min(rows, key=lambda x: x['J'])
    rec.update(best=best, all=rows, budget=len(rows))
    json.dump(rec, open(path, 'w'), indent=1)
    return best


if __name__ == '__main__' and os.environ.get('TUNE_STAGE') == '2':
    methods = [m for m in sys.argv[1:]] or ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC']
    with Pool(16) as pool:
        for mth in methods:
            b = refine(mth, pool)
            print(f"{mth} (stage 2): J={b['J']:.3f} rmse={b['rmse']:.3f} params={b['params']}", flush=True)
