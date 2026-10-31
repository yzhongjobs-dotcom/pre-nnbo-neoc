"""Equal-budget tuning after the setup audit.

Same budget for every method (41 random + 20 local evaluations), development scenes from all
five families (scenes_v2, 6 each), robust objective: mean over scenes of log(1 + J), where J is
the benchmark's physical composite (force error/0.5 N, contact loss, peak excess). The log keeps
one extreme stress scene from deciding the parameters. Test and replication seeds are not used.
Run: python -X utf8 tune_v2.py [methods]  ->  tuned_v2/<method>.json
"""
import json, math, os, sys, time
import numpy as np
from multiprocessing import Pool
import common
import tune_equal as TE
import scenes_v2 as SV

OUT = os.path.join(common.HERE, 'tuned_v2'); os.makedirs(OUT, exist_ok=True)
N_RANDOM, N_LOCAL = 40, 20


VARIANTS = ('Proposed-IMM', 'Proposed-KF', 'Proposed-IMMg', 'Proposed-KFg')


def make(method, params, p_scene=None):
    if method == 'YuenType':
        from comparators import YuenType
        return YuenType(p_scene, **params)
    if method in VARIANTS:
        import variants
        return variants.build(method, dict(params))
    return TE.make(method, params)


def space(method, r):
    if method in VARIANTS:
        import variants
        q = TE.space('Proposed', r)
        q['kgood_rule'] = str(r.choice(['precise', 'consistent']))
        q.update(variants.space(method, r))
        return q
    if method == 'Proposed':
        q = TE.space('Proposed', r)
        q['kgood_rule'] = str(r.choice(['precise', 'consistent']))
        return q
    import space_patch as SP
    return SP.space(method, r)


def perturb(method, params, r):
    if method == 'YuenType':
        import space_patch as SP
        return SP.yuen_perturb(params, r)
    if method in VARIANTS:
        import variants
        base = {k: v for k, v in params.items() if k not in variants.DEFAULTS[method] and k != 'kgood_rule'}
        q = TE.perturb('Proposed', base, r)
        q['kgood_rule'] = params.get('kgood_rule', 'precise')
        q.update(variants.perturb(method, {k: params[k] for k in variants.DEFAULTS[method] if k in params}, r))
        return q
    q = TE.perturb(method, {k: v for k, v in params.items() if k != 'kgood_rule'}, r)
    if method == 'Proposed':
        q['kgood_rule'] = params.get('kgood_rule', 'precise') if r.random() > .2 else str(r.choice(['precise', 'consistent']))
    else:
        for k, v in params.items():
            if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
                q[k] = v                       # keep discrete choices (e.g. Np) unless re-drawn below
        if 'Np' in params and r.random() < .3:
            q['Np'] = int(r.choice([8, 15, 25, 40]))
        for k, v in params.items():
            if isinstance(v, str):
                q[k] = v
        q = clamp(method, q, params)
    return q


BOUNDS = dict(K_I=(2., 120.), lam=(.8, .999), Q=(1., 8.), A_E=(2e-6, 1e-4), f_exc=(2., 30.),
              K_P=(0., 1.5), probe_hz=(20., 125.), innovation_scale=(.2, 5.), m=(.3, 10.), b=(50., 5000.),
              k_h0=(500., 5e4), eso_bandwidth=(2., 60.), q_over_r=(1e3, 1e8), bp_Q=(1., 8.))


def clamp(method, q, orig):
    """keep stage-2 proposals inside the stage-1 ranges and the DO-MPAC nominal pole screen."""
    for k, (lo, hi) in BOUNDS.items():
        if k in q and isinstance(q[k], float):
            if k == 'lam' and method != 'WijType':
                lo = .95
            q[k] = float(min(max(q[k], lo), hi))
    if 'q_over_r' in q:
        import space_patch as SP
        if not SP.nominal_ok(q['q_over_r'], q['Np']):
            q['q_over_r'], q['Np'] = orig['q_over_r'], orig['Np']
    return q


def _eval(args):
    method, params, (fam, sid, p, T) = args
    try:
        log, _ = common.run(make(method, params, p), p, T, sid)
        m = common.task_metrics(log, 20.8)
        J = m['J'] if np.isfinite(m['J']) else 1e3
    except Exception as ex:
        return dict(L=math.log1p(1e3), rmse=np.nan, err=str(ex))
    return dict(L=math.log1p(min(J, 1e3)), rmse=m['rmse'])


def evaluate(method, cands, pool, sc):
    res = pool.map(_eval, [(method, c, s) for c in cands for s in sc], chunksize=1)
    k = len(sc); rows = []
    for i, c in enumerate(cands):
        rr = res[i*k:(i+1)*k]
        rows.append(dict(params=c, L=float(np.mean([x['L'] for x in rr])),
                         rmse=float(np.nanmean([x['rmse'] for x in rr]))))
        if method == 'YuenType':
            rows[-1]['n_errors'] = sum('err' in x for x in rr)
    return rows


def tune(method, pool):
    sc = SV.scenes('dev')
    r = np.random.default_rng(777)
    import comparators as C, proposed as Pm
    default = {} if method != 'Proposed' else {}
    rows = evaluate(method, [default]+[space(method, r) for _ in range(N_RANDOM)], pool, sc)
    best = min(rows, key=lambda x: x['L'])
    if not best['params'] and method in VARIANTS:
        import variants
        best = min([x for x in rows if x['params']], key=lambda x: x['L'])
    base = best['params'] or (dict(getattr(C, method).params) if method != 'Proposed' else
                              {k: v for k, v in Pm.MotionAware().params.items() if k in TE.space('Proposed', r)})
    r2 = np.random.default_rng(888)
    rows += evaluate(method, [perturb(method, base, r2) for _ in range(N_LOCAL)], pool, sc)
    best = min(rows, key=lambda x: x['L'])
    rec = dict(method=method, budget=len(rows), objective='mean log(1+J)', best=best, all=rows)
    if method == 'YuenType':
        rec.update(objective='mean log(1+min(J,1e3))', split='dev',
                   scenes=[dict(family=f, sid=sid, T=T) for f, sid, _, T in sc],
                   random_seed=777, local_seed=888, sensor=dict(f_s=30., tau_m=.03, sigma_m=.0001, seed=0),
                   defaults=dict(C.YuenType.params))
    json.dump(rec, open(os.path.join(OUT, method+'.json'), 'w'), indent=1, default=str)
    return best


if __name__ == '__main__':
    methods = sys.argv[1:] or ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC']
    with Pool(16) as pool:
        for m in methods:
            t0 = time.time(); b = tune(m, pool)
            print(f"{m}: L={b['L']:.3f} rmse={b['rmse']:.3f} params={b['params']} ({time.time()-t0:.0f}s)", flush=True)
