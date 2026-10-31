"""Claim-driven simulation study after the setup audit (paper draft §五).

S1  estimation on identical closed-loop data (the observations the controller consumed), full
    push events, phase-wise statistics (breathing / ramp / hold / static); independent seeds for
    the main and the replication cohort.
S2  static heterogeneous scan (stiffness span only).
S3  motion: T3m, T2m-soft, T2m-stiff (moderate pushes); T2f-stiff reported as a stress test.
Ablation on the S3 families. Timing in a single process (no pool).
Failures are counted separately; RMSE statistics are conditional on success and say so; contact
retention is reported as loss fraction, longest loss and episode count.
Env: SPLIT=test|replication (default test); METHOD_SET=v2 (tuned_v2) ; EXTRA=comma list of extra
proposed variants registered in variants.py.
Set METHODS=comma-separated method names to select comparators, including YuenType.
Run: python -X utf8 run_v2.py
"""
import json, math, os, time
import numpy as np
from multiprocessing import Pool
import common
from common import S, P
import tune_equal as TE
import scenes_v2 as SV

SPLIT = os.environ.get('SPLIT', 'test')
RES = os.path.join(common.HERE, 'results_v2', SPLIT); os.makedirs(RES, exist_ok=True)
TUNED = os.path.join(common.HERE, 'tuned_v2')
F_EXC = 20.8
METHODS = ['Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC']
if os.environ.get('METHODS'):
    METHODS = [x.strip() for x in os.environ['METHODS'].split(',') if x.strip()]
FINAL = os.environ.get('FINAL', 'Proposed-IMM')      # final method; 'Proposed' below refers to it
ALT = 'Abl-IMM' if FINAL == 'Proposed-KF' else 'Abl-KF'
ABL = ['Abl-polyFF', ALT, 'Abl-noProjection', 'Abl-noConservative', 'Abl-noFeedforward', 'Abl-noExcitation']
EXTRA = [x for x in os.environ.get('EXTRA', '').split(',') if x]


def tuned(m):
    return json.load(open(os.path.join(TUNED, m+'.json')))['best']['params']


def final_params():
    return dict(tuned(FINAL))


def build(name, p_scene=None):
    import proposed as Pm, variants
    if name == 'YuenType':
        from comparators import YuenType
        return YuenType(p_scene, **tuned(name))
    if name == 'Proposed':
        pol = variants.build(FINAL, final_params()); pol.name = 'Proposed'
        return pol
    if name in METHODS:
        return TE.make(name, tuned(name))
    if name in EXTRA:
        return variants.build(name, tuned(name))
    base = final_params()
    kw = {k: base[k] for k in ('M1', 'D1', 'KPF1', 'ki')}
    if name == 'Oracle':
        return Pm.Oracle(p_scene, **kw)
    if name.startswith('Fixed:'):
        return Pm.FixedEnd(float(name.split(':')[1]), **kw)
    if name == 'Abl-polyFF':                     # earlier window-polynomial feed-forward, tuned separately
        pol = Pm.MotionAware(**tuned('Proposed')); pol.name = name; return pol
    if name == 'Abl-KF':                         # single-model Kalman tracker, tuned separately
        pol = variants.build('Proposed-KF', tuned('Proposed-KF')); pol.name = name; return pol
    if name == 'Abl-IMM':                        # two-model IMM tracker, tuned separately
        pol = variants.build('Proposed-IMM', tuned('Proposed-IMM')); pol.name = name; return pol
    flag = {'Abl-noProjection': ('degree', 0), 'Abl-noConservative': ('conservative', False),
            'Abl-noFeedforward': ('feedforward', False), 'Abl-noExcitation': ('excitation', False)}[name]
    q = dict(base); q[flag[0]] = flag[1]
    pol = variants.build(FINAL, q); pol.name = name
    return pol


def _episode(args):
    name, fam, sid, p, T = args
    pol = build(name, p)
    log, _ = common.run(pol, p, T, sid)
    m = common.task_metrics(log, F_EXC)
    ok = not m['diverged']
    rip = m['ripple'] if ok else np.nan
    return dict(method=name, family=fam, sid=sid, ok=bool(ok), rmse=m['rmse'], p95=m['p95'], peak=m['peak'],
                n_loss=m['n_loss'], loss_frac=m['loss'], loss_max=m['loss_max'], osc=m['osc_ex'], ripple=rip,
                rmse_noexc=float(math.sqrt(max(m['rmse']**2-rip**2, 0.))) if ok else np.nan, J=m['J'])


def summarize(rows):
    tab = {}
    for key in sorted({(r['family'], r['method']) for r in rows}):
        rs = [r for r in rows if (r['family'], r['method']) == key]
        ok = [r for r in rs if r['ok']]
        def mean(k): return float(np.nanmean([r[k] for r in ok])) if ok else np.nan
        def sd(k): return float(np.nanstd([r[k] for r in ok])) if ok else np.nan
        tab['%s|%s' % key] = dict(n=len(rs), n_fail=len(rs)-len(ok), rmse=[mean('rmse'), sd('rmse')],
                                  rmse_noexc=mean('rmse_noexc'), p95=mean('p95'),
                                  peak_max=float(np.nanmax([r['peak'] for r in ok])) if ok else np.nan,
                                  loss_frac=mean('loss_frac'), loss_max=float(np.nanmax([r['loss_max'] for r in ok])) if ok else np.nan,
                                  n_loss=float(np.nansum([r['n_loss'] for r in ok])), osc=mean('osc'), ripple=mean('ripple'))
    return tab


def paired(rows, fam, a, b, key='rmse', n_boot=4000):
    ra = {r['sid']: r for r in rows if r['family'] == fam and r['method'] == a}
    rb = {r['sid']: r for r in rows if r['family'] == fam and r['method'] == b}
    ks = [k for k in sorted(set(ra) & set(rb)) if ra[k]['ok'] and rb[k]['ok']]
    d = np.array([rb[k][key]-ra[k][key] for k in ks], float)
    if len(d) == 0:
        return dict(n=0, mean=np.nan, ci=[np.nan, np.nan])
    g = np.random.default_rng(0); bs = [np.mean(g.choice(d, len(d))) for _ in range(n_boot)]
    return dict(n=int(len(d)), mean=float(d.mean()), ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                wins=int(np.sum(d > 0)), per_scene=[float(x) for x in d])


# ---------------- S1 ----------------
def s1_scenes():
    base = {'test': 50000, 'replication': 60000, 'confirm': 70000}.get(SPLIT, 75000)
    r = np.random.default_rng(base); out = []
    for E, lvl in ((1.5e4, 'soft'), (1.5e5, 'mid'), (1.0e6, 'stiff')):
        for motion in ('static', 'breathing', 'push'):
            for i in range(4):
                p, T = S.make('T2', r, dict(ramp=float(r.uniform(.2, .5))))
                p = np.array(p); p[0] = E*float(r.uniform(.8, 1.25))
                if motion != 'breathing': p[2] = 0.
                if motion != 'push': p[15] = 0.
                if motion == 'push':
                    p[15] = 1; p[16:20] = [2.5, float(r.choice([-1, 1]))*float(r.uniform(.003, .008)), float(r.uniform(.2, .5)), 1.0]
                out.append((lvl, motion, base+len(out), p, 8.0))
    return out


def _phase(t, p, motion):
    lab = np.full(len(t), motion, object)
    if motion == 'push':
        t0, amp, tr, hold = p[16:20]
        lab[:] = 'push-rest'
        lab[((t >= t0) & (t < t0+tr)) | ((t >= t0+tr+hold) & (t < t0+2*tr+hold))] = 'push-ramp'
        lab[(t >= t0+tr) & (t < t0+tr+hold)] = 'push-hold'
    return lab


def _s1(args):
    import comparators as C, variants
    lvl, motion, sid, p, T, exc = args
    q = final_params(); q['excitation'] = exc
    pol = variants.build(FINAL, q); base = pol.params
    log, _ = common.run(pol, p, T, sid)
    if log is None:
        return dict(level=lvl, motion=motion, exc=exc, diverged=1)
    d = common.DirectPolicy(); d.reset(P.PLATFORMS['PC125'])
    rp = tuned('RossiType')
    # Rossi estimator with the matched (force-domain) displacement, as in the tuned comparator;
    # it is fed the same closed-loop data (whose excitation is the proposed method's 20.8 Hz)
    rest = C._RossiEstimator(d.Tc, F_EXC, rp.get('Q', 3.), rp.get('lam', .99), 100., 3., .01)
    md = C._MatchedDisplacement(d)
    c, nw = base['c'], base['nw']
    labs = _phase(log[:, 0], p, motion)
    acc = {}
    for n in range(nw+3, len(log)):
        kt = log[n-1, 8]
        kr = rest.update(log[n-1, 5], md.update(P.Obs(log, n)))
        if log[n-1, 0] < 1.0 or kt <= 0:
            continue
        ph = labs[n-1]
        a = acc.setdefault(ph, dict(static=[], proj=[], rossi=[], cover=[]))
        for key, deg in (('static', 0), ('proj', 3)):
            f = common.window_fit(log, n, nw, deg, base['noise'], d.ns, d.eb, d.afc, d.ndel)
            if np.isfinite(f['k']) and f['info'] > 0:
                a[key].append(f['k']/kt)
                if key == 'proj':
                    a['cover'].append(float(f['k']+c/math.sqrt(f['info']) >= kt))
        a['rossi'].append(kr/kt)
    out = dict(level=lvl, motion=motion, exc=exc, diverged=0, ripple=common.task_metrics(log, F_EXC)['ripple'], phases={})
    for ph, a in acc.items():
        out['phases'][ph] = {k: ([float(np.nanpercentile(v, q)) for q in (10, 50, 90)] if k != 'cover' else float(np.mean(v)))
                             if len(v) else None for k, v in a.items()}
    return out


def timing(names, n_scenes=3):
    out = {}
    sc = SV.scenes('dev', ('T3m', 'T2m-stiff'))[:n_scenes]+SV.scenes('dev', ('T2m-stiff',))[:n_scenes]
    for name in names:
        calls = []; worst = 0.
        for fam, sid, p, T in sc:
            pol = build(name, p); common.run(pol, p, T, sid)
            tc = np.array(pol.t_call[10:])*1e3; calls.extend(tc.tolist()); worst = max(worst, float(np.percentile(tc, 99)))
        calls = np.array(calls)
        out[name] = dict(p50=float(np.median(calls)), p99=float(np.percentile(calls, 99)), worst_scene_p99=worst,
                         frac_over_8ms=float(np.mean(calls > 8.)), n=int(len(calls)))
    return out


if __name__ == '__main__':
    t0 = time.time()
    sc = SV.scenes(SPLIT)
    klo, khi = 1109., 28083.          # development stiffness range (5th/95th pct), as before
    main = METHODS+['Oracle', 'Fixed:%d' % klo, 'Fixed:%d' % khi]+EXTRA
    with Pool(16) as pool:
        s1 = pool.map(_s1, [(*x, ex) for x in s1_scenes() for ex in (True, False)], chunksize=1)
        json.dump(s1, open(os.path.join(RES, 'S1.json'), 'w'), indent=1)
        rows = pool.map(_episode, [(m, *x) for m in main for x in sc], chunksize=1)
        json.dump(rows, open(os.path.join(RES, 'main.json'), 'w'), indent=1)
        ab = pool.map(_episode, [(m, *x) for m in ABL for x in sc if x[0] != 'S2'], chunksize=1)
        json.dump(ab, open(os.path.join(RES, 'ablation.json'), 'w'), indent=1)
    tim = timing(METHODS+['Abl-polyFF']+EXTRA)
    allrows = rows+ab
    summ = dict(split=SPLIT, summary=summarize(allrows), timing=tim,
                paired={f: {m: paired(allrows, f, 'Proposed', m) for m in main[1:]+ABL if m != 'Proposed'}
                        for f in SV.FAMILIES},
                wall_s=time.time()-t0)
    for x in EXTRA:
        summ['paired_extra_'+x] = {f: {m: paired(allrows, f, x, m) for m in main if m != x} for f in SV.FAMILIES}
    json.dump(summ, open(os.path.join(RES, 'summary.json'), 'w'), indent=1)
    print('done', SPLIT, round(time.time()-t0), 's')
