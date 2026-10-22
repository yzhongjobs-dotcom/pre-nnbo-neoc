"""Reviewer-gap experiments (research/paper_experiments_draft.md, Sec. V-G, gaps 2, 3, 5, 10).

No existing file is changed; this module imports the frozen methods and writes gaps_results.json.

1. c2fresh   Candidate 2 (full variant, the three trained seeds, no retraining) on the seed-5000 fresh cohort
             (40 scenes) of all six cells of candidate1_v3.py. B1, EtR-A, B3b 0.1 s and C1v3 are re-run on the
             same scenes and noise seeds and checked against the per-scene J stored in candidate1_v3_results.json.
2. plat NAME Platform sensitivity on the two main cells (T3 0.5 s / 1 mm, T1 r30 / 1 mm), test seed 6000
             (never used before), 40 scenes. Platforms (PLATS): PC125 (nominal, reference row), PC125dead,
             PC500 (negative control) and a servo-lag sweep on PC125 (tau_a 0.5x, 2x nominal; optional 2 cycles
             of dead time). Per platform:
               B1        retuned with tune.best_fixed (the gate.py procedure: grid + Nelder-Mead, 20 train scenes,
                         peak rule), i.e. the platform's own best fixed controller;
               B1@PC125  the PC125 B1 gains transferred (reference only);
               B3b 0.1 s rule retuned on train (gate.py procedure, 0.1 s only);
               EtR-A     (a) PC125 tuning transferred, (b) retuned on train (baselines_adaptive.tune_etr procedure,
                         72 configurations);
               C1v3      (a) PC125 configuration transferred, robot model constants of the true platform;
                         (a') same, but the policy's internal robot model keeps the PC125 latency (ndel 1,
                              tau_a 19 ms) -- the unidentified-latency case (only where Tc = 8 ms);
                         (b) retuned: the candidate1_model.tune_variant search (49 configurations, 20 train
                             scenes, v1 policy) on the platform, wrapped exactly like v3 (v2 supervisor, k0 3e4,
                             fallback = the platform's retuned B1);
               C2        the three PC125-trained seeds transferred unchanged (no retraining).
3. t2        T2 push task (homogeneous soft or stiff tissue with patient pushes; false-adaptation test), PC125,
             dose ramp 0.5 s / breathing 1 mm, test seed 6000: T2-tuned B1 and B3b, the T3-cell settings
             transferred (B1, EtR-A, C1v3) and retuned on T2 train (EtR-A, C1v3 with T2 B1 fallback), C2.

Selection uses train (1000+) only; seeds 3000/4000/5000 are seen; 6000 is the new test seed.
Run (<= 2 worker processes, single-threaded BLAS): python -B gaps_eval.py c2fresh | plat PC125dead | t2 | report
"""
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ[_k] = '1'
import json, sys, time, math
import numpy as np
import plant as P
import scenarios as S
import tune as U
import candidate1_model as C
import candidate1_v2 as V2
import candidate1_v3 as V3
import candidate2_learned as C2
import baselines_adaptive as BA
from baselines_ref import TruthK

OUT = 'gaps_results.json'
NPROC = 2
TR0, FRESH5, FRESH6 = 1000, 5000, 6000
NTR, NTE = 20, 40
BASE = dict(P.PLATFORMS['PC125'])
NOMINAL_MODEL = dict(ndel=BASE['ndel'], tau_a=BASE['tau_a'])
PLATS = {
    'PC125': dict(BASE),
    'PC125dead': dict(P.PLATFORMS['PC125dead']),
    'PC500': dict(P.PLATFORMS['PC500']),
    'PC125_lag0.5': dict(BASE, tau_a=0.5 * BASE['tau_a']),
    'PC125_lag2': dict(BASE, tau_a=2.0 * BASE['tau_a']),
    'PC125_dead2': dict(BASE, ndel=2),
}
C2_SEEDS = (0, 1, 2)
MAIN = ['T3', 'T1']


# ============================================================================ worker side
class C1Belief(V3.Candidate1v3):
    """Candidate1v3 whose internal robot model (ndel, tau_a) may differ from the plant it runs on."""
    def __init__(self, belief=None, **kw):
        super().__init__(**kw); self.belief = belief

    def reset(self, pl):
        super().reset(dict(pl, **self.belief) if self.belief else pl)


_NETS = {}


def _net(variant, seed):
    if (variant, seed) not in _NETS: _NETS[(variant, seed)] = C2.load_net(variant, seed)
    return _NETS[(variant, seed)]


def _warm():
    V3.warmup()
    p, T = S.make('N1', np.random.default_rng(0), dict(E=2e4, A=0.001))
    net, meta = _net('full', 0)
    P.episode(p, C2.Cand2Policy(net), T=0.3, seed=0)


def _job(a):
    p, T, spec, plat, seed = a
    kind, kw = spec
    pl = PLATS[plat] if isinstance(plat, str) else plat
    upd = None
    if kind == 'fixed': ctrl = tuple(kw)
    elif kind == 'B3b': ctrl = TruthK(*kw)
    elif kind == 'EtR': ctrl = BA.EtR(*kw[0], **kw[1])
    elif kind == 'C1': ctrl = C.Candidate1(**kw); upd = ctrl.update
    elif kind == 'C1v3':
        kw = dict(kw); b = kw.pop('belief', None); ctrl = C1Belief(belief=b, **kw); upd = ctrl.update
    elif kind == 'C2':
        net, meta = _net(*kw); cf = meta['cfg']
        ctrl = C2.Cand2Policy(net, every=cf['every'], gains_only=cf['gains_only'], noforce=cf['noforce'])
    else: raise ValueError(kind)
    log = P.episode(p, ctrl, T=T, platform=pl, seed=seed, update=upd)
    m = S.metrics(log, Tc=pl['Tc'])
    if kind == 'C1v3': m['n_loss_entries'] = ctrl.n_loss_entries
    return m


_POOL = None


def pool():
    global _POOL
    if _POOL is None:
        from multiprocessing import Pool
        _POOL = Pool(NPROC, initializer=_warm)
        U._pool = _POOL                      # tune.best_fixed / tune.evaluate use the same 2 workers
    return _POOL


def run(scenes, specs, plat, seed0):
    jobs = [(p, T, s, plat, seed0 + i) for i, (p, T) in enumerate(scenes) for s in specs]
    res = pool().map(_job, jobs, chunksize=max(1, len(jobs) // (8 * NPROC)))
    ns = len(specs)
    return [res[i * ns:(i + 1) * ns] for i in range(len(scenes))]


# ============================================================================ summaries and CIs
def summ(rows):
    out = U.summary(rows)
    Js = np.array([r['J'] for r in rows])
    out['J_median'] = float(np.median(Js)); out['n_J_over_20'] = int((Js > 20).sum())
    out['frac_scenes_with_loss'] = float(np.mean([r['loss'] > 0 for r in rows]))
    return out


def per_scene(rows):
    return dict(J=[float(r['J']) for r in rows], rmse=[float(r['rmse']) for r in rows],
                over7=[float(not (r['peak'] <= 7.0)) for r in rows], lossscene=[float(r['loss'] > 0) for r in rows],
                loss=[float(r['loss']) for r in rows], peak=[float(r['peak']) for r in rows])


def seedavg(list_of_rows):
    """per-scene averages over seeds (C2)."""
    ps = [per_scene(r) for r in list_of_rows]
    return {k: list(np.mean([q[k] for q in ps], 0)) for k in ps[0]}


def cis(a, b):
    return dict(dJ=U.paired_ci(a['J'], b['J']), drmse=U.paired_ci(a['rmse'], b['rmse']),
                dover7=U.paired_ci(a['over7'], b['over7']), dlossscene=U.paired_ci(a['lossscene'], b['lossscene']),
                wins=int(sum(x < y for x, y in zip(a['J'], b['J']))))


def entry(rowsets):
    """rowsets: list of row lists (one per seed, or one)."""
    if len(rowsets) == 1:
        return dict(summary=summ(rowsets[0]), per_scene=per_scene(rowsets[0]))
    ss = [summ(r) for r in rowsets]
    keys = [k for k in ss[0] if isinstance(ss[0][k], (int, float))]
    return dict(summary={k: float(np.mean([s[k] for s in ss])) for k in keys},
                summary_sd={k: float(np.std([s[k] for s in ss])) for k in keys},
                per_seed=[s for s in ss], per_scene=seedavg(rowsets))


def add_cis(out, refs):
    for nm, e in out.items():
        e['ci'] = {}
        for rf in refs:
            if rf in out and rf != nm:
                e['ci']['vs_' + rf] = cis(e['per_scene'], out[rf]['per_scene'])


# ============================================================================ results file
def load():
    try: return json.load(open(OUT))
    except Exception: return {}


def save(R):
    json.dump(R, open(OUT, 'w'), indent=1, default=float)


LOGF = None


def log(s):
    print(s, flush=True)
    if LOGF: open(LOGF, 'a').write(s + '\n')


def line(nm, s):
    return (f"  {nm:22s} J {s['J']:8.3f} (med {s['J_median']:.3f}) rmse {s['rmse']:.3f} >7N {s['frac_peak_over']:.3f} "
            f"loss {s['loss']:.4f} lossscenes {s['frac_scenes_with_loss']:.3f} J>20 {s['n_J_over_20']:.1f}")


# ============================================================================ 1. C2 on the seed-5000 cohort
def c2fresh(R):
    v3 = json.load(open('candidate1_v3_results.json'))['final_seed5000']
    res = R.setdefault('c2_seed5000', {})
    for cell in V3.ORDER:
        t0 = time.time()
        task, dose = V2.CELLS[cell]; scenes = S.cohort(task, NTE, FRESH5, dose)
        bs, bn, src = V2.baselines(cell)
        specs = [('C1v3', V3.v3cfg(cell))] + bs + [('C2', ('full', s)) for s in C2_SEEDS]
        names = ['C1v3'] + bn + ['C2_s%d' % s for s in C2_SEEDS]
        r = run(scenes, specs, 'PC125', FRESH5)
        per = {nm: [x[j] for x in r] for j, nm in enumerate(names)}
        out = {nm: entry([per[nm]]) for nm in ['C1v3'] + bn}
        out['C2'] = entry([per['C2_s%d' % s] for s in C2_SEEDS])
        # same scenes as candidate1_v3_results.json?
        ver = {}
        for nm in ['C1v3'] + bn:
            old = np.array(v3[cell]['methods'][nm]['J_per_scene']); new = np.array(out[nm]['per_scene']['J'])
            ver[nm] = dict(max_abs_dJ=float(np.max(np.abs(old - new))), n_differ_1e6=int(np.sum(np.abs(old - new) > 1e-6)),
                           J_old=float(old.mean()), J_new=float(new.mean()))
        add_cis(out, ['B1', 'EtR-A', 'B3b_0.1', 'C1v3'])
        for s, nm in zip(C2_SEEDS, ['C2_s%d' % s for s in C2_SEEDS]):
            ps = per_scene(per[nm])
            out['C2'].setdefault('per_seed_ci', {})[str(s)] = {rf: cis(ps, out[rf]['per_scene']) for rf in ['B1', 'EtR-A']}
        res[cell] = dict(task=task, dose=dose, seed0=FRESH5, n=NTE, etr_source=src, verify_vs_v3_json=ver, methods=out,
                         seconds=time.time() - t0)
        save(R)
        log(f"[c2fresh {cell}] verify {json.dumps({k: (round(v['max_abs_dJ'], 6), v['n_differ_1e6']) for k, v in ver.items()})}")
        for nm in out:
            ex = ''
            if 'vs_B1' in out[nm]['ci']:
                a = out[nm]['ci']['vs_B1']['dJ']; b = out[nm]['ci'].get('vs_EtR-A', {}).get('dJ', (np.nan,) * 3)
                ex = f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f}) vs EtR-A {b[0]:+.3f} ({b[1]:+.3f},{b[2]:+.3f})"
            log(line(nm, out[nm]['summary']) + ex)


# ============================================================================ tuning procedures on a platform
def tune_b1(train, plat):
    b = U.best_fixed(train, PLATS[plat], seed0=TR0)
    return dict(gains=list(b['gains']), train_J=b['train_J'], grid_best=list(b['grid_best']),
                any_feasible_on_grid=b['any_feasible_on_grid'])


def select(r):
    return C.select(r)


def tune_b3b(train, plat, lat=0.1):
    import gate as G
    ctr = [('B3b', (lat,) + tuple(q)) for q in G.RULE]
    r = run(train, ctr, plat, TR0); j, Jm, f = select(r)
    return dict(rule=list(G.RULE[j]), train_J=float(Jm[j]), train_feasible=bool(f[j]))


def tune_etr(train, plat):
    """baselines_adaptive.tune_etr (as re-implemented in candidate1_v2.etr_tune) on the given platform."""
    E1 = dict(BA.EST0); Q = BA.RULE
    r = run(train, [('EtR', (tuple(q), dict(E1, rule='truthk'))) for q in Q], plat, TR0); j, Jm, f = select(r)
    order = np.argsort(np.where(f, Jm, Jm + 1e6)); tops = [Q[i] for i in order[:3]]
    cand = [(q, e) for q in tops for e in BA.EST_VAR]
    r2 = run(train, [('EtR', (tuple(q), dict(e, rule='truthk'))) for q, e in cand], plat, TR0); j2, Jm2, f2 = select(r2)
    if (not f2[j2] and f[j]) or (f[j] == f2[j2] and Jm[j] <= Jm2[j2]): best = (Q[j], E1, float(Jm[j]), bool(f[j]))
    else: best = (cand[j2][0], cand[j2][1], float(Jm2[j2]), bool(f2[j2]))
    return dict(rule_params=list(best[0]), est_params=best[1], train_J=best[2], train_feasible=best[3],
                n_configs=len(Q) + len(cand))


def etr_spec_of(t):
    return ('EtR', (tuple(t['rule_params']), dict(t['est_params'], rule='truthk')))


def tune_c1(train, plat):
    """candidate1_model.tune_variant('C1') on the given platform: hand-set default + 32 random (rng 12345) +
    16 one-parameter perturbations (rng 777) of the stage-A best = 49 configurations, v1 policy."""
    fixed = C.VARIANTS['C1']
    rng = np.random.default_rng(12345)
    A = [dict()] + C.stage_a(rng)
    rA = run(train, [('C1', dict(c, **fixed)) for c in A], plat, TR0); jA, JmA, fA = select(rA)
    B = C.stage_b(A[jA], np.random.default_rng(777))
    rB = run(train, [('C1', dict(c, **fixed)) for c in B], plat, TR0); _, JmB, fB = select(rB)
    Jall = np.concatenate([JmA, JmB]); fall = np.concatenate([fA, fB]); cand = A + B
    j = int(np.argmin(np.where(fall, Jall, Jall + 1e6)))
    return dict(config=dict(cand[j], **fixed), train_J=float(Jall[j]), train_feasible=bool(fall[j]), n_configs=len(cand),
                n_feasible=int(fall.sum()))


def v3wrap(base, b1gains):
    c = dict(C.DEFAULTS); c.update(base); c.update(k0=3e4, **V3.SUP)
    g = b1gains; c['fb'] = (g[0], g[1], g[3], g[4]); c['solver'] = 'ord'; c['iters'] = 1
    return c


# ============================================================================ 2. platform sensitivity
def platform(R, plat, cells=MAIN):
    res = R.setdefault('platform_seed6000', {}).setdefault(plat, {})
    pl = PLATS[plat]
    for cell in cells:
        if 'methods' in res.get(cell, {}): log(f"[{plat} {cell}] already complete, skipped"); continue
        t0 = time.time()
        task, dose = V2.CELLS[cell]
        train = S.cohort(task, NTR, TR0, dose); test = S.cohort(task, NTE, FRESH6, dose)
        g125 = V2.gate_ref(cell)
        etr125, _ = V2.etr_spec(cell)
        tun = res.get(cell, {}).get('tuning', {})
        if plat == 'PC125':
            tun = dict(B1=dict(gains=list(g125['B1']['gains']), source='gate_PC125.json'),
                       B3b=dict(rule=list(g125['B3b']['0.1']['rule']), source='gate_PC125.json'))
        else:
            if 'B1' not in tun:
                tun['B1'] = tune_b1(train, plat); log(f"[{plat} {cell}] B1 {tun['B1']} ({time.time() - t0:.0f}s)"); res[cell] = dict(task=task, dose=dose, platform=pl, tuning=tun); save(R)
            if 'B3b' not in tun:
                tun['B3b'] = tune_b3b(train, plat); log(f"[{plat} {cell}] B3b {tun['B3b']} ({time.time() - t0:.0f}s)"); res[cell] = dict(task=task, dose=dose, platform=pl, tuning=tun); save(R)
            if 'EtR-A' not in tun:
                tun['EtR-A'] = tune_etr(train, plat); log(f"[{plat} {cell}] EtR-A {tun['EtR-A']} ({time.time() - t0:.0f}s)"); res[cell] = dict(task=task, dose=dose, platform=pl, tuning=tun); save(R)
            if 'C1' not in tun:
                tun['C1'] = tune_c1(train, plat); log(f"[{plat} {cell}] C1 {tun['C1']} ({time.time() - t0:.0f}s)"); res[cell] = dict(task=task, dose=dose, platform=pl, tuning=tun); save(R)
        res[cell] = dict(task=task, dose=dose, platform=pl, tuning=tun); save(R)
        b1 = tun['B1']['gains']
        specs = {'B1': ('fixed', tuple(b1)), 'B3b_0.1': ('B3b', (0.1,) + tuple(tun['B3b']['rule'])),
                 'EtR-A_tr': etr125, 'C1v3_tr': ('C1v3', V3.v3cfg(cell))}
        if plat != 'PC125':
            specs['B1@PC125'] = ('fixed', tuple(g125['B1']['gains']))
            specs['EtR-A_rt'] = etr_spec_of(tun['EtR-A'])
            specs['C1v3_rt'] = ('C1v3', v3wrap(tun['C1']['config'], b1))
            if pl['Tc'] == BASE['Tc']:
                specs['C1v3_tr_nominalmodel'] = ('C1v3', dict(V3.v3cfg(cell), belief=NOMINAL_MODEL))
        names = list(specs)
        sp = [specs[k] for k in names] + [('C2', ('full', s)) for s in C2_SEEDS]
        r = run(test, sp, plat, FRESH6)
        per = {nm: [x[j] for x in r] for j, nm in enumerate(names)}
        out = {nm: entry([per[nm]]) for nm in names}
        out['C2_tr'] = entry([[x[len(names) + q] for x in r] for q in range(len(C2_SEEDS))])
        add_cis(out, ['B1', 'EtR-A_tr', 'EtR-A_rt', 'B3b_0.1'])
        res[cell].update(methods=out, seed0=FRESH6, n=NTE, seconds=time.time() - t0); save(R)
        log(f"[{plat} {cell} seed {FRESH6}] ({time.time() - t0:.0f}s)")
        for nm in out:
            ex = ''
            if 'vs_B1' in out[nm]['ci']:
                a = out[nm]['ci']['vs_B1']['dJ']; ex = f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f})"
            log(line(nm, out[nm]['summary']) + ex)



# ============================================================================ 2b. zero-latency oracle (headroom check)
def b3b0(R, plats):
    """B3b at 0 s information delay (rule retuned on train, gate.py procedure), added to finished platform cells."""
    for plat in plats:
        for cell in MAIN:
            c = R.get('platform_seed6000', {}).get(plat, {}).get(cell)
            if not c or 'methods' not in c or 'B3b_0.0' in c['methods']: continue
            t0 = time.time(); task, dose = V2.CELLS[cell]
            train = S.cohort(task, NTR, TR0, dose); test = S.cohort(task, NTE, FRESH6, dose)
            if plat == 'PC125': c['tuning']['B3b0'] = dict(rule=list(V2.gate_ref(cell)['B3b']['0.0']['rule']), source='gate_PC125.json')
            else: c['tuning']['B3b0'] = tune_b3b(train, plat, lat=0.0)
            r = run(test, [('B3b', (0.0,) + tuple(c['tuning']['B3b0']['rule']))], plat, FRESH6)
            c['methods']['B3b_0.0'] = entry([[x[0] for x in r]])
            add_cis(c['methods'], ['B1', 'EtR-A_tr', 'EtR-A_rt', 'B3b_0.1'])
            save(R)
            a = c['methods']['B3b_0.0']['ci']['vs_B1']['dJ']
            log(f"[b3b0 {plat} {cell}] rule {c['tuning']['B3b0']['rule']} " + line('B3b_0.0', c['methods']['B3b_0.0']['summary']) +
                f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f}) ({time.time() - t0:.0f}s)")

# ============================================================================ 3. T2 push task
T2DOSE = dict(ramp=0.5, A=0.001)


def t2(R):
    t0 = time.time(); plat = 'PC125'
    train = S.cohort('T2', NTR, TR0, T2DOSE); test = S.cohort('T2', NTE, FRESH6, T2DOSE)
    res = R.setdefault('T2_seed6000', dict(task='T2', dose=T2DOSE, platform=plat))
    tun = res.setdefault('tuning', {})
    if 'B1' not in tun: tun['B1'] = tune_b1(train, plat); log(f"[T2] B1 {tun['B1']} ({time.time() - t0:.0f}s)"); save(R)
    if 'B3b' not in tun: tun['B3b'] = tune_b3b(train, plat); log(f"[T2] B3b {tun['B3b']}"); save(R)
    if 'EtR-A' not in tun: tun['EtR-A'] = tune_etr(train, plat); log(f"[T2] EtR-A {tun['EtR-A']}"); save(R)
    if 'C1' not in tun: tun['C1'] = tune_c1(train, plat); log(f"[T2] C1 {tun['C1']} ({time.time() - t0:.0f}s)"); save(R)
    g3 = V2.gate_ref('T3'); etr3, _ = V2.etr_spec('T3')
    specs = {'B1': ('fixed', tuple(tun['B1']['gains'])), 'B1@T3': ('fixed', tuple(g3['B1']['gains'])),
             'B3b_0.1': ('B3b', (0.1,) + tuple(tun['B3b']['rule'])),
             'EtR-A_tr': etr3, 'EtR-A_rt': etr_spec_of(tun['EtR-A']),
             'C1v3_tr': ('C1v3', V3.v3cfg('T3')), 'C1v3_rt': ('C1v3', v3wrap(tun['C1']['config'], tun['B1']['gains']))}
    names = list(specs)
    sp = [specs[k] for k in names] + [('C2', ('full', s)) for s in C2_SEEDS]
    r = run(test, sp, plat, FRESH6)
    per = {nm: [x[j] for x in r] for j, nm in enumerate(names)}
    out = {nm: entry([per[nm]]) for nm in names}
    out['C2_tr'] = entry([[x[len(names) + q] for x in r] for q in range(len(C2_SEEDS))])
    add_cis(out, ['B1', 'B1@T3', 'EtR-A_tr', 'EtR-A_rt', 'B3b_0.1'])
    stiff = [bool(p[0] > 1e5) for p, T in test]
    sub = {}
    for lab, msk in [('soft', [not s for s in stiff]), ('stiff', stiff)]:
        sub[lab] = {nm: {k: float(np.mean([v for v, m in zip(out[nm]['per_scene'][k], msk) if m]))
                         for k in ['J', 'rmse', 'over7', 'lossscene', 'loss']} for nm in out}
        sub[lab]['n'] = int(sum(msk))
    res.update(methods=out, by_tissue=sub, seed0=FRESH6, n=NTE, stiff_per_scene=stiff, seconds=time.time() - t0); save(R)
    log(f"[T2 seed {FRESH6}] ({time.time() - t0:.0f}s), stiff scenes {sum(stiff)}/{len(stiff)}")
    for nm in out:
        ex = ''
        if 'vs_B1' in out[nm]['ci']:
            a = out[nm]['ci']['vs_B1']['dJ']; ex = f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f})"
        log(line(nm, out[nm]['summary']) + ex)


# ============================================================================ report (markdown tables)
def fci(c, d=3):
    return f"{c[0]:+.{d}f} ({c[1]:+.{d}f}, {c[2]:+.{d}f})"


def row(nm, e, refs):
    s = e['summary']; sd = e.get('summary_sd')
    J = f"{s['J']:.3f}" + (f" ± {sd['J']:.3f}" if sd else '')
    cells = [nm, J, f"{s['rmse']:.3f}", f"{s['frac_peak_over']:.3f}", f"{s['loss']:.4f} / {s['frac_scenes_with_loss']:.2f}",
             (f"{s['n_J_over_20']:.1f}" if sd else f"{int(s['n_J_over_20'])}")]
    for rf in refs:
        c = e['ci'].get('vs_' + rf)
        cells.append(fci(c['dJ']) if c else '–')
    return '| ' + ' | '.join(cells) + ' |'


def report(R):
    L = []
    hdr = lambda refs: ['| Method | J | RMSE N | frac peak > 7 N | loss time / loss scenes | J > 20 | ' +
                        ' | '.join('ΔJ vs ' + r for r in refs) + ' |', '|' + '---|' * (6 + len(refs))]
    if 'c2_seed5000' in R:
        L.append('## 1. C2 on the seed-5000 cohort\n')
        for cell, c in R['c2_seed5000'].items():
            refs = ['B1', 'EtR-A', 'C1v3']
            L.append(f"### {cell} {c['task']} {c['dose']} (EtR-A from {c['etr_source']})\n")
            L.append('verify vs candidate1_v3_results.json (max |ΔJ| per scene, scenes differing): ' +
                     ', '.join(f"{k} {v['max_abs_dJ']:.2g}/{v['n_differ_1e6']}" for k, v in c['verify_vs_v3_json'].items()) + '\n')
            L += hdr(refs)
            for nm in ['B1', 'EtR-A', 'B3b_0.1', 'C1v3', 'C2']:
                L.append(row(nm, c['methods'][nm], refs))
            ps = c['methods']['C2']['per_seed_ci']
            L.append('\nC2 per seed ΔJ vs B1: ' + '; '.join(f"s{s} {fci(v['B1']['dJ'])}" for s, v in ps.items()) +
                     '. vs EtR-A: ' + '; '.join(f"s{s} {fci(v['EtR-A']['dJ'])}" for s, v in ps.items()) + '\n')
    if 'platform_seed6000' in R:
        L.append('## 2. Platform sensitivity (test seed 6000, 40 scenes)\n')
        for plat, pc in R['platform_seed6000'].items():
            for cell, c in pc.items():
                if 'methods' not in c: continue
                refs = ['B1', 'EtR-A_tr', 'EtR-A_rt'] if plat != 'PC125' else ['B1', 'EtR-A_tr']
                L.append(f"### {plat} {json.dumps(c['platform'])} — {cell}\n")
                tu = c['tuning']
                L.append(f"B1 gains {np.round(tu['B1']['gains'], 6).tolist()}" +
                         (f", train J {tu['B1']['train_J']:.3f}" if 'train_J' in tu['B1'] else '') +
                         (f"; EtR-A rt {tu['EtR-A']['rule_params']} train J {tu['EtR-A']['train_J']:.3f} feas {tu['EtR-A']['train_feasible']}" if 'EtR-A' in tu else '') +
                         (f"; C1 rt train J {tu['C1']['train_J']:.3f} feas {tu['C1']['train_feasible']} ({tu['C1']['n_feasible']}/49 feasible)" if 'C1' in tu else '') + '\n')
                L += hdr(refs)
                for nm in c['methods']:
                    L.append(row(nm, c['methods'][nm], refs))
                L.append('')
    if 'T2_seed6000' in R and 'methods' in R['T2_seed6000']:
        c = R['T2_seed6000']; refs = ['B1', 'B1@T3', 'EtR-A_rt']
        L.append(f"## 3. T2 push task {c['dose']} PC125 seed 6000\n")
        L += hdr(refs)
        for nm in c['methods']: L.append(row(nm, c['methods'][nm], refs))
        for lab in ['soft', 'stiff']:
            b = c['by_tissue'][lab]
            L.append(f"\n{lab} (n {b['n']}): " + '; '.join(f"{nm} J {v['J']:.3f} rmse {v['rmse']:.3f} >7N {v['over7']:.2f} lossscenes {v['lossscene']:.2f}"
                                                          for nm, v in b.items() if nm != 'n'))
    print('\n'.join(L))


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'report'
    LOGF = '/tmp/claude-0/-home-user/f9a8790b-fc73-5507-98e3-ab5273024fda/scratchpad/gaps_console.log'
    R = load()
    t0 = time.time()
    if what == 'c2fresh': c2fresh(R)
    elif what == 'plat':
        for pn in sys.argv[2].split(','):
            platform(R, pn, sys.argv[3].split(',') if len(sys.argv) > 3 else MAIN)
    elif what == 't2': t2(R)
    elif what == 'b3b0': b3b0(R, sys.argv[2].split(','))
    elif what == 'report': report(R)
    if what != 'report':
        R.setdefault('runs', []).append(dict(what=' '.join(sys.argv[1:]), seconds=time.time() - t0, finished=time.strftime('%Y-%m-%d %H:%M:%S')))
        save(R)
    if _POOL is not None:
        _POOL.close(); _POOL.join()
