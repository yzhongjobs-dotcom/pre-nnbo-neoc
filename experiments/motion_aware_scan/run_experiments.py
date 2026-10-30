"""Claim-driven simulation study (methods §10; paper draft §四).

S1 mechanism   : estimators on identical closed-loop data, homogeneous tissue at three
                 stiffness levels x {static, breathing, push}, with/without excitation;
                 plus the combination hypothesis (Rossi estimate inside its own loop vs
                 inside the DO-MPAC disturbance-compensation loop).
S2 stiffness   : static heterogeneous scans (T1: ridges 10-30x) - consistency of the force loop.
S3 main task   : spine scan with breathing and push events (T3) and push/retract on soft and
                 stiff tissue (T2).
Ablation       : S3 scenes, proposed method minus one part at a time.
Timing         : per-call controller time (P50/P99), all methods.
Every method uses the parameters selected by tune_equal.py on development seeds (equal
budget); test seeds below were never used for development.
Run: python -X utf8 run_experiments.py  -> results/*.json, results/tables.md
"""
import json, math, os, time
import numpy as np
from multiprocessing import Pool
import common
from common import S, P
import tune_equal as TU

RES = os.path.join(common.HERE, 'results'); os.makedirs(RES, exist_ok=True)
F_EXC = 20.8
TEST = dict(T1=('T1', 12, 9500), T2=('T2', 12, 9600), T3=('T3', 12, 9700))
if os.environ.get('REPLICATE') == '1':     # fresh replication cohort
    TEST = dict(T1=('T1', 12, 10500), T2=('T2', 12, 10600), T3=('T3', 12, 10700))
    RES = os.path.join(common.HERE, 'results_replication'); os.makedirs(RES, exist_ok=True)
METHODS = ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC']


def tuned(method):
    return json.load(open(os.path.join(TU.OUT, method+'.json')))['best']['params']


def build(name, p_scene=None):
    import proposed as Pm, comparators as C
    if name in METHODS:
        return TU.make(name, tuned(name))
    base = tuned('Proposed')
    if name == 'Oracle':
        kw = {k: base[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
        return Pm.Oracle(p_scene, **kw)
    if name.startswith('Fixed'):
        kw = {k: base[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
        return Pm.FixedEnd(float(name.split(':')[1]), **kw)
    abl = dict(**base)
    if name == 'Abl-noProjection': abl['degree'] = 0
    elif name == 'Abl-noConservative': abl['conservative'] = False
    elif name == 'Abl-noFeedforward': abl['feedforward'] = False
    elif name == 'Abl-noExcitation': abl['excitation'] = False
    else: raise ValueError(name)
    pol = Pm.MotionAware(**abl); pol.name = name
    return pol


def _episode(args):
    name, task, sid, p, T = args
    pol = build(name, p)
    t0 = time.time()
    log, _ = common.run(pol, p, T, sid)
    m = common.task_metrics(log, F_EXC)
    tc = np.array(pol.t_call[10:])*1e3 if len(pol.t_call) > 10 else np.array([np.nan])
    kap = np.nan; cover = np.nan
    if log is not None and pol.trace and len(pol.trace[0]) >= 2:
        tr = np.array([r[:2] for r in pol.trace], float)
        n = min(len(tr), len(log)); kt = log[1:n+1, 8] if len(log) > n else log[:n, 8]
        kt = log[:n, 8]; sel = (tr[:n, 0] > 1.) & (kt > 0) & np.isfinite(tr[:n, 1])
        if sel.any():
            r = kt[sel]/tr[:n][sel, 1]
            kap = [float(np.percentile(r, q)) for q in (5, 50, 95)]
            cover = float(np.mean(r <= 1.0))
    return dict(method=name, task=task, sid=sid, wall=time.time()-t0,
                rmse=m['rmse'], p95=m['p95'], peak=m['peak'], n_loss=m['n_loss'], loss=m['loss'],
                osc=m['osc_ex'], ripple=m['ripple'], J=m['J'], diverged=m['diverged'],
                t_p50=float(np.nanpercentile(tc, 50)), t_p99=float(np.nanpercentile(tc, 99)),
                kappa=kap, k_le_kbar=cover)


def scene_list(spec, split_t2=False):
    task, n, seed = spec
    out = []
    for i, (p, T) in enumerate(S.cohort(task, n, seed)):
        tag = task
        if split_t2:
            tag = 'T2-stiff' if p[0] > 2e5 else 'T2-soft'
        out.append((tag, seed+i, p, T))
    return out


# ---------------- S1: estimation mechanism on identical data ----------------
def s1_scenes():
    out = []
    r = np.random.default_rng(9800)
    for E, lvl in ((1.5e4, 'soft'), (1.5e5, 'mid'), (1.0e6, 'stiff')):
        for motion in ('static', 'breathing', 'push'):
            for i in range(4):
                p, T = S.make('T2', r)
                p = np.array(p); p[0] = E*float(r.uniform(.8, 1.25))
                if motion != 'breathing': p[2] = 0.
                if motion != 'push': p[15] = 0.
                if motion == 'push' and p[15] == 0:
                    p[15] = 1; p[16:20] = [2.5, 0.006, 0.25, 1.0]
                out.append((lvl, motion, 9800+len(out), p, 6.0))
    return out


def _s1(args):
    import comparators as C
    lvl, motion, sid, p, T, exc = args
    import proposed as Pm
    base = dict(Pm.MotionAware(**tuned('Proposed')).params); base['excitation'] = exc
    pol = Pm.MotionAware(**base)
    log, _ = common.run(pol, p, T, sid)
    if log is None:
        return dict(level=lvl, motion=motion, exc=exc, diverged=1)
    pl = P.PLATFORMS['PC125']; d = common.DirectPolicy(); d.reset(pl)
    rest = C._RossiEstimator(d.Tc, F_EXC, tuned('RossiType')['Q'], tuned('RossiType')['lam'], 100., 3., .01)
    out = dict(level=lvl, motion=motion, exc=exc, diverged=0, static=[], proj=[], rossi=[], cover=[])
    c = base['c']; nw = base['nw']
    for n in range(nw+3, len(log)):
        kt = log[n-1, 8]
        kr = rest.update(log[n-1, 5], log[n-1, 2])
        if log[n-1, 0] < 1.0 or kt <= 0:
            continue
        for key, deg in (('static', 0), ('proj', 3)):
            f = common.window_fit(log, n, nw, deg, base['noise'], d.ns, d.eb, d.afc, d.ndel)
            if np.isfinite(f['k']) and f['info'] > 0:
                out[key].append(f['k']/kt)
                if key == 'proj':
                    out['cover'].append(float(f['k']+c/math.sqrt(f['info']) >= kt))
        out['rossi'].append(kr/kt)
    m = common.task_metrics(log, F_EXC)
    out['ripple'] = m['ripple']
    for key in ('static', 'proj', 'rossi'):
        a = np.array(out[key])
        out[key] = [float(np.percentile(a, q)) for q in (10, 50, 90)] if len(a) else [np.nan]*3
    out['cover'] = float(np.mean(out['cover'])) if out['cover'] else np.nan
    return out


def _combo(args):
    """Rossi estimate quality inside its own loop vs inside the DO-MPAC loop (hypothesis 3)."""
    name, task, sid, p, T = args
    pol = build(name)
    log, _ = common.run(pol, p, T, sid)
    if log is None or not pol.trace:
        return dict(method=name, task=task, ratio=[np.nan]*3)
    tr = np.array(pol.trace, float); n = min(len(tr), len(log)); kt = log[:n, 8]
    sel = (tr[:n, 0] > 1.) & (kt > 0)
    r = tr[:n][sel, 1]/kt[sel]
    return dict(method=name, task=task, ratio=[float(np.percentile(r, q)) for q in (10, 50, 90)])


def summarize(rows, key_task='task'):
    out = {}
    for r in rows:
        out.setdefault((r[key_task], r['method']), []).append(r)
    tab = {}
    for (task, mth), rs in out.items():
        def agg(k):
            a = np.array([x[k] for x in rs], float)
            return [float(np.nanmean(a)), float(np.nanstd(a))]
        tab[f'{task}|{mth}'] = dict(n=len(rs), rmse=agg('rmse'), p95=agg('p95'), peak_max=float(np.nanmax([x['peak'] for x in rs])),
                                     peak=agg('peak'), n_loss=float(np.nansum([x['n_loss'] for x in rs])), osc=agg('osc'),
                                     ripple=agg('ripple'), J=agg('J'), diverged=int(sum(x['diverged'] for x in rs)),
                                     t_p99=float(np.nanmedian([x['t_p99'] for x in rs])), t_p50=float(np.nanmedian([x['t_p50'] for x in rs])))
    return tab


def paired(rows, task, a, b, n_boot=4000, seed=0):
    """paired bootstrap CI of mean RMSE difference (b - a) over common scenes."""
    ra = {r['sid']: r['rmse'] for r in rows if r['task'] == task and r['method'] == a}
    rb = {r['sid']: r['rmse'] for r in rows if r['task'] == task and r['method'] == b}
    ks = sorted(set(ra) & set(rb)); d = np.array([rb[k]-ra[k] for k in ks], float)
    d = d[np.isfinite(d)]
    if len(d) == 0: return [np.nan]*3
    g = np.random.default_rng(seed); bs = [np.mean(g.choice(d, len(d))) for _ in range(n_boot)]
    return [float(np.mean(d)), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def k_bounds():
    ks = []
    for task, n, seed in TU.DEV:
        for p, T in S.cohort(task, n, seed):
            import proposed as Pm
            log, _ = common.run(Pm.Oracle(p, **{k: tuned('Proposed')[k] for k in ('M1', 'D1', 'KPF1', 'ki')}), p, T, 1)
            if log is not None:
                kk = log[log[:, 8] > 0, 8]; ks.extend(kk.tolist())
    return float(np.percentile(ks, 5)), float(np.percentile(ks, 95))


if __name__ == '__main__':
    t_start = time.time()
    klo, khi = k_bounds()
    fixed = [f'Fixed:{klo:.0f}', f'Fixed:{khi:.0f}']
    main_methods = METHODS+['Oracle']+fixed
    abl = ['Abl-noProjection', 'Abl-noConservative', 'Abl-noFeedforward', 'Abl-noExcitation']
    with Pool(16) as pool:
        s1 = pool.map(_s1, [(*sc, ex) for sc in s1_scenes() for ex in (True, False)], chunksize=1)
        json.dump(s1, open(os.path.join(RES, 'S1_estimation.json'), 'w'), indent=1)
        s2sc = scene_list(TEST['T1'])
        s2 = pool.map(_episode, [(m, *sc) for m in main_methods for sc in s2sc], chunksize=1)
        json.dump(s2, open(os.path.join(RES, 'S2_stiffness_span.json'), 'w'), indent=1)
        s3sc = scene_list(TEST['T3'])+scene_list(TEST['T2'], split_t2=True)
        s3 = pool.map(_episode, [(m, *sc) for m in main_methods for sc in s3sc], chunksize=1)
        json.dump(s3, open(os.path.join(RES, 'S3_main.json'), 'w'), indent=1)
        ab = pool.map(_episode, [(m, *sc) for m in abl for sc in s3sc], chunksize=1)
        json.dump(ab, open(os.path.join(RES, 'S3_ablation.json'), 'w'), indent=1)
        combo = pool.map(_combo, [(m, *sc) for m in ('RossiType', 'RossiDOMPAC') for sc in s3sc], chunksize=1)
        json.dump(combo, open(os.path.join(RES, 'S1_combination.json'), 'w'), indent=1)
    summary = dict(k_bounds=[klo, khi], S2=summarize(s2), S3=summarize(s3), ablation=summarize(ab + [r for r in s3 if r['method'] == 'Proposed']),
                   paired={t: {m: paired(s3+s2, t, 'Proposed', m) for m in main_methods[1:]+abl}
                           for t in ('T1', 'T3', 'T2-soft', 'T2-stiff')},
                   wall_s=time.time()-t_start)
    for t in ('T1', 'T3', 'T2-soft', 'T2-stiff'):
        for m in abl:
            summary['paired'][t][m] = paired(s3+ab, t, 'Proposed', m)
    json.dump(summary, open(os.path.join(RES, 'summary.json'), 'w'), indent=1)
    print(json.dumps(dict(k_bounds=[klo, khi], wall_s=summary['wall_s']), indent=1))
