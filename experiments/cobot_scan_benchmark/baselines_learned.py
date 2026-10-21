"""E27-style learned variable impedance control (comparison baseline, causal).

Reference: Li et al., "Impedance learning-based adaptive force tracking for robot on unknown
terrains", T-RO 2025 (E27).  E27 pre-trains, in simulation on randomly generated terrain, an
MLP that tunes impedance parameters (dK, dB, dx_r) from observations o=[x, e] (robot state and
force error) with an actor-critic DRL agent, reward r = r_const - e' Re e, and runs it at a
0.05 s control/update period (E27:579, 719-720).

Our transfer to the shared benchmark interface:
  * policy = one-hidden-layer tanh MLP (as E27 eq. 15), numpy only;
  * inputs (12, causal only): force error averaged over the last four 48 ms windows, force
    slope, encoder velocity, command-encoder lag, force std (oscillation), 200 ms encoder
    displacement and force change, and the previously applied D and KPF (log columns 9/10);
    never k_true or log columns 3,4,6,7,8;
  * outputs: log D and log KPF as bounded offsets (+-2 e-folds) around learned base values;
    M and ki are learned constants, K = 0 (as in the B3b rule family); u_ff is not used;
  * update period 0.05 s (plant.episode(..., update=0.05)); ticks in between hold the gains;
  * trained with separable CMA-ES (Ros & Hansen 2008) on the episode cost instead of DRL
    (no autodiff available; ~130 parameters; episode-level composite J with peak/loss
    terms; common random numbers across a generation).  Fitness = minibatch mean of J / J_ref, with
    J_ref the starting fixed controller BASE on the same scene and noise seed, so T1 (small J)
    and T3 (large J) scenes count by relative change while rare failures stay linear in J
    (a first attempt with mean log J improved train fitness but blew up on dev: aborted).
  * training scenes: mixed T1/T3 with randomised doses, scene rng seeds 1000 + generation
    (train split, < 2000); selection of the checkpoint by a dev cohort (seed 2000) at the two
    reported doses, subject to the B1 feasibility rule (peak <= 7 N in >= 95 % of scenes).

Usage:  python baselines_learned.py train SEED [GENS]     -> learned_vic_seed{SEED}.npz
        python baselines_learned.py test SEED [SEED ...]  -> learned_results.json
"""
import json, sys, time, numpy as np
import plant as P, scenarios as S

UPDATE = 0.05
NF, NH, NO = 12, 8, 2
NSCAL = 4                                   # logD0, logKPF0, logM, ki
NPAR = NF * NH + NH + NH * NO + NO + NSCAL
LOGD0, LOGKPF0, LOGM0, KI0 = np.log(120.), np.log(1.15e-4), np.log(1.5), 10.   # start: a dev-feasible fixed point
SPAN = 2.0                                  # +- e-folds of the variable part
LD_LO, LD_HI = np.log(20.), np.log(2500.)
LK_LO, LK_HI = np.log(1e-6), np.log(1e-3)
BASE = (1.5, 120., 0., 10., 1.15e-4)      # the policy's starting point (network output zero)
CELLS = [('T3', dict(ramp=0.5, A=0.001)), ('T1', dict(ratio=30., A=0.001))]


def unpack(th):
    i = 0
    W1 = th[i:i + NF * NH].reshape(NH, NF); i += NF * NH
    b1 = th[i:i + NH]; i += NH
    W2 = th[i:i + NH * NO].reshape(NO, NH); i += NH * NO
    b2 = th[i:i + NO]; i += NO
    s = th[i:i + NSCAL]
    return W1, b1, W2, b2, s


class LearnedVIC:
    """Causal policy: reads only obs.t/force/y/y_cmd and its own applied D, KPF (log 9/10)."""
    oracle = False

    def __init__(self, theta, zero_net=False):
        self.W1, self.b1, self.W2, self.b2, s = unpack(np.asarray(theta, float))
        self.zero = zero_net
        self.lD0 = float(np.clip(LOGD0 + s[0], LD_LO, LD_HI))
        self.lK0 = float(np.clip(LOGKPF0 + s[1], LK_LO, LK_HI))
        self.M = float(np.clip(np.exp(LOGM0 + 0.7 * s[2]), 0.25, 3.0))
        self.ki = float(np.clip(KI0 * np.exp(0.7 * s[3]), 0., 30.))

    def reset(self, pl):
        self.Tc = pl['Tc']; self.w = max(int(round(UPDATE / self.Tc)), 1)

    def features(self, o):
        n = o.n; w = self.w; Tc = self.Tc
        F = o.log[:n, 5]; y = o.log[:n, 2]; yc = o.log[:n, 11]
        x = np.empty(NF)
        for k in range(4):                                   # force error, 4 windows
            hi = max(n - k * w, 1); lo = min(max(n - (k + 1) * w, 0), hi - 1)
            x[k] = P.FD - F[lo:hi].mean()
        j = max(n - 1 - w, 0); j4 = max(n - 1 - 4 * w, 0)
        x[4] = (F[n - 1] - F[j]) / (w * Tc) / 20.
        x[5] = (y[n - 1] - y[j]) / (w * Tc) / 0.01
        x[6] = (yc[n - 1] - y[n - 1]) / 1e-3
        x[7] = F[max(n - 2 * w, 0):n].std() / 0.3
        x[8] = (np.log(max(o.log[n - 1, 9], 1e-9)) - self.lD0) / SPAN
        x[9] = (np.log(max(o.log[n - 1, 10], 1e-9)) - self.lK0) / SPAN
        x[10] = (y[n - 1] - y[j4]) / 1e-3
        x[11] = F[n - 1] - F[j4]
        return np.clip(x, -5., 5.)

    def __call__(self, o):
        if o.n == 0 or self.zero:
            lD, lK = self.lD0, self.lK0
        else:
            z = self.W2 @ np.tanh(self.W1 @ self.features(o) + self.b1) + self.b2
            lD = self.lD0 + SPAN * np.tanh(z[0]); lK = self.lK0 + SPAN * np.tanh(z[1])
        lD = min(max(lD, LD_LO), LD_HI); lK = min(max(lK, LK_LO), LK_HI)
        return (self.M, float(np.exp(lD)), 0., self.ki, float(np.exp(lK)))


def init_theta(seed):
    rng = np.random.default_rng(seed)
    th = np.zeros(NPAR)
    th[:NF * NH] = rng.standard_normal(NF * NH) / np.sqrt(NF)   # W1 random, W2 = 0 -> fixed start
    return th


# ---------------- evaluation (worker side) ----------------
def _run(job):
    th, zero, p, T, nseed, plat, full = job
    c0 = time.process_time()
    pol = LearnedVIC(th, zero) if th is not None else None
    log = P.episode(p, pol if pol is not None else full, T=T, platform=plat, seed=nseed,
                    update=UPDATE if pol is not None else None)
    m = S.metrics(log, Tc=P.PLATFORMS[plat]['Tc'])
    if log is not None and pol is not None:
        D = log[:, 9]; K = log[:, 10]
        m['D_med'] = float(np.median(D)); m['D_p5'] = float(np.percentile(D, 5)); m['D_p95'] = float(np.percentile(D, 95))
        m['KPF_med'] = float(np.median(K)); m['KPF_p5'] = float(np.percentile(K, 5)); m['KPF_p95'] = float(np.percentile(K, 95))
    m['cpu'] = time.process_time() - c0
    return m


def _run_ctrl(job):
    """evaluate an arbitrary picklable controller (B1 gains tuple or TruthK) for paired CIs."""
    ctrl, p, T, nseed, plat = job
    c0 = time.process_time()
    log = P.episode(p, ctrl, T=T, platform=plat, seed=nseed)
    m = S.metrics(log, Tc=P.PLATFORMS[plat]['Tc']); m['cpu'] = time.process_time() - c0
    return m


# ---------------- training ----------------
def train_batch(g, n_each=4):
    """mixed T1/T3 scenes with randomised doses; scene rng seed 1000+g (train split)."""
    seed = S.SPLIT['train'] + g; assert seed < S.SPLIT['dev']
    rng = np.random.default_rng(seed); out = []
    for task in ['T3'] * n_each + ['T1'] * n_each:
        dose = dict(A=rng.uniform(0.0005, 0.005), ratio=rng.uniform(10, 30))
        if task == 'T3': dose['ramp'] = rng.uniform(0.2, 0.6)
        out.append(S.make(task, rng, dose))
    return out


def dev_set(n=20):
    return [(t, S.cohort(t, n, S.SPLIT['dev'], d)) for t, d in CELLS]


class SepCMA:
    def __init__(self, m, sigma, lam, seed):
        n = len(m); self.n = n; self.m = m.copy(); self.sigma = sigma; self.lam = lam
        mu = lam // 2; w = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1)); self.w = w / w.sum()
        self.mu = mu; self.mueff = 1. / np.sum(self.w ** 2); me = self.mueff
        self.cc = (4 + me / n) / (n + 4 + 2 * me / n); self.cs = (me + 2) / (n + me + 5)
        c1 = 2 / ((n + 1.3) ** 2 + me); cmu = min(1 - c1, 2 * (me - 2 + 1 / me) / ((n + 2) ** 2 + me))
        self.c1 = c1 * (n + 2) / 3; self.cmu = min(1 - self.c1, cmu * (n + 2) / 3)
        self.damps = 1 + 2 * max(0, np.sqrt((me - 1) / (n + 1)) - 1) + self.cs
        self.chiN = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
        self.pc = np.zeros(n); self.ps = np.zeros(n); self.C = np.ones(n); self.g = 0
        self.rng = np.random.default_rng(seed)

    def ask(self):
        z = self.rng.standard_normal((self.lam, self.n))
        return self.m + self.sigma * np.sqrt(self.C) * z

    def tell(self, X, f):
        idx = np.argsort(f); Y = (X[idx[:self.mu]] - self.m) / self.sigma
        ym = self.w @ Y; self.m = self.m + self.sigma * ym
        self.ps = (1 - self.cs) * self.ps + np.sqrt(self.cs * (2 - self.cs) * self.mueff) * ym / np.sqrt(self.C)
        hs = np.linalg.norm(self.ps) / np.sqrt(1 - (1 - self.cs) ** (2 * (self.g + 1))) / self.chiN < 1.4 + 2 / (self.n + 1)
        self.pc = (1 - self.cc) * self.pc + hs * np.sqrt(self.cc * (2 - self.cc) * self.mueff) * ym
        self.C = ((1 - self.c1 - self.cmu) * self.C + self.c1 * (self.pc ** 2 + (1 - hs) * self.cc * (2 - self.cc) * self.C)
                  + self.cmu * (self.w @ (Y ** 2)))
        self.sigma *= np.exp(self.cs / self.damps * (np.linalg.norm(self.ps) / self.chiN - 1))
        self.g += 1


def dev_eval(pool, th, dev, plat='PC125'):
    out = {}; cpu = 0.; ne = 0
    for task, sc in dev:
        jobs = [(th, False, p, T, S.SPLIT['dev'] + i, plat, None) for i, (p, T) in enumerate(sc)]
        r = pool.map(_run, jobs, chunksize=5); cpu += sum(m['cpu'] for m in r); ne += len(r)
        pk = np.array([m['peak'] if np.isfinite(m['peak']) else 1e9 for m in r])
        out[task] = dict(J=float(np.mean([m['J'] for m in r])), rmse=float(np.nanmean([m['rmse'] for m in r])),
                         feas=bool(np.mean(pk <= P.FD + 3.0) >= 0.95))
    return out, cpu, ne


def train(seed, gens=150, lam=16, sigma0=0.08, dev_every=10, plat='PC125'):
    from multiprocessing import Pool
    t0 = time.time(); cpu = 0.; neps = 0
    es = SepCMA(init_theta(seed), sigma0, lam, seed)
    dev = dev_set(); hist = []
    with Pool(2) as pool:
        d0, c, ne = dev_eval(pool, es.m, dev, plat); cpu += c; neps += ne
        ref = {t: d0[t]['J'] for t in d0}

        def score(d):
            return float(np.mean([d[t]['J'] / ref[t] for t in d])) + (0 if all(d[t]['feas'] for t in d) else 1e3)
        best = (score(d0), es.m.copy(), 0, d0)
        print(f'seed {seed} gen 0 dev {d0} score {best[0]:.3f}', flush=True)
        for g in range(gens):
            X = es.ask(); sc = train_batch(g)
            jobs = [(X[k], False, p, T, 10 ** 6 + 100 * g + i, plat, None) for k in range(lam) for i, (p, T) in enumerate(sc)]
            jobs += [(None, False, p, T, 10 ** 6 + 100 * g + i, plat, BASE) for i, (p, T) in enumerate(sc)]
            r = pool.map(_run, jobs, chunksize=8); cpu += sum(m['cpu'] for m in r); neps += len(r)
            Jm = np.array([m['J'] for m in r[:lam * len(sc)]]).reshape(lam, len(sc))
            Jref = np.array([m['J'] for m in r[lam * len(sc):]])
            f = (Jm / (Jref + 0.01)).mean(1); es.tell(X, f)
            rec = dict(gen=g + 1, f_best=float(f.min()), f_mean=float(f.mean()), sigma=float(es.sigma))
            if (g + 1) % dev_every == 0 or g + 1 == gens:
                d, c, ne = dev_eval(pool, es.m, dev, plat); cpu += c; neps += ne
                s = score(d); rec['dev'] = d; rec['dev_score'] = s
                if s < best[0]: best = (s, es.m.copy(), g + 1, d)
                print(f'seed {seed} gen {g+1} f {f.mean():.3f} sigma {es.sigma:.3f} dev T3 {d["T3"]["J"]:.3f} T1 {d["T1"]["J"]:.3f} '
                      f'score {s:.3f} best@{best[2]} cpu {cpu:.0f}s wall {time.time()-t0:.0f}s', flush=True)
            hist.append(rec)
    np.savez(f'learned_vic_seed{seed}.npz', theta=best[1], theta_final=es.m, best_gen=best[2], dev_score=best[0],
             dev_ref=json.dumps(ref), cfg=json.dumps(dict(NF=NF, NH=NH, NO=NO, UPDATE=UPDATE, lam=lam, sigma0=sigma0, gens=gens)))
    json.dump(dict(seed=seed, gens=gens, lam=lam, sigma0=sigma0, episodes=neps, cpu_s=cpu, wall_s=time.time() - t0,
                   best_gen=best[2], dev_best=best[3], dev_ref=ref, hist=hist),
              open(f'learned_vic_seed{seed}_train.json', 'w'), indent=1)
    return best


# ---------------- test ----------------
def test(seeds, plat='PC125', NTE=40):
    from multiprocessing import Pool
    from baselines_ref import TruthK
    import tune as U
    gate = {(r['task'], json.dumps(r['dose'], sort_keys=True)): r for r in json.load(open('gate_PC125.json'))}
    out = dict(method='E27-style learned VIC (numpy MLP, sep-CMA-ES)', platform=plat, update_s=UPDATE,
               arch=dict(inputs=NF, hidden=NH, outputs=['logD', 'logKPF'], learned_constants=['D0', 'KPF0', 'M', 'ki'], K=0.,
                         n_params=NPAR), cells=[], training={})
    cpu = 0.
    for s in seeds:
        tr = json.load(open(f'learned_vic_seed{s}_train.json'))
        out['training'][str(s)] = {k: tr[k] for k in ['gens', 'lam', 'episodes', 'cpu_s', 'wall_s', 'best_gen', 'dev_best', 'dev_ref']}
    with Pool(2) as pool:
        for task, dose in CELLS:
            test_sc = S.cohort(task, NTE, S.SPLIT['test'], dose); g = gate[(task, json.dumps(dose, sort_keys=True))]
            cell = dict(task=task, dose=dose, gate_B1=g['B1']['test'], gate_B1_gains=g['B1']['gains'],
                        gate_B3b={k: v['test'] for k, v in g['B3b'].items()}, learned={}, learned_zero_net={})
            # re-run B1 and B3b@0.1 on the same test cohort for paired CIs
            b3 = TruthK(0.1, *g['B3b']['0.1']['rule'])
            rb1 = pool.map(_run_ctrl, [(tuple(g['B1']['gains']), p, T, S.SPLIT['test'] + i, plat) for i, (p, T) in enumerate(test_sc)])
            rb3 = pool.map(_run_ctrl, [(b3, p, T, S.SPLIT['test'] + i, plat) for i, (p, T) in enumerate(test_sc)])
            cpu += sum(m['cpu'] for m in rb1 + rb3)
            cell['B1_rerun'] = U.summary(rb1); cell['B3b01_rerun'] = U.summary(rb3)
            jB1 = [m['J'] for m in rb1]; jB3 = [m['J'] for m in rb3]
            allJ = []
            for s in seeds:
                th = np.load(f'learned_vic_seed{s}.npz')['theta']
                for zero, key in [(False, 'learned'), (True, 'learned_zero_net')]:
                    r = pool.map(_run, [(th, zero, p, T, S.SPLIT['test'] + i, plat, None) for i, (p, T) in enumerate(test_sc)])
                    cpu += sum(m['cpu'] for m in r); jl = [m['J'] for m in r]
                    e = dict(test=U.summary(r), dJ_vs_B1=U.paired_ci(jl, jB1), dJ_vs_B3b01=U.paired_ci(jl, jB3),
                             drmse_vs_B1=U.paired_ci([m['rmse'] for m in r], [m['rmse'] for m in rb1]))
                    if not zero:
                        allJ.append(jl)
                        e['gains'] = {k: float(np.mean([m[k] for m in r if k in m])) for k in
                                      ['D_med', 'D_p5', 'D_p95', 'KPF_med', 'KPF_p5', 'KPF_p95']}
                        pol = LearnedVIC(th); e['constants'] = dict(M=pol.M, ki=pol.ki, D0=float(np.exp(pol.lD0)), KPF0=float(np.exp(pol.lK0)))
                    cell[key][str(s)] = e
            cell['learned_seed_mean'] = {k: float(np.mean([cell['learned'][str(s)]['test'][k] for s in seeds]))
                                         for k in ['J', 'rmse', 'peak', 'frac_peak_over', 'inband', 'loss']}
            cell['learned_seed_std'] = {k: float(np.std([cell['learned'][str(s)]['test'][k] for s in seeds]))
                                        for k in ['J', 'rmse', 'peak', 'frac_peak_over', 'inband', 'loss']}
            out['cells'].append(cell)
    out['test_cpu_s'] = cpu
    out['total_train_cpu_s'] = float(sum(out['training'][str(s)]['cpu_s'] for s in seeds))
    out['total_train_episodes'] = int(sum(out['training'][str(s)]['episodes'] for s in seeds))
    json.dump(out, open('learned_results.json', 'w'), indent=1)
    for c in out['cells']:
        print(c['task'], c['dose'], 'B1', {k: round(c['B1_rerun'][k], 3) for k in ['J', 'rmse', 'peak', 'frac_peak_over', 'inband']},
              'B3b0.1', {k: round(c['B3b01_rerun'][k], 3) for k in ['J', 'rmse', 'peak', 'frac_peak_over', 'inband']})
        for s in seeds:
            e = c['learned'][str(s)]; z = c['learned_zero_net'][str(s)]
            print('  seed', s, {k: round(e['test'][k], 3) for k in ['J', 'rmse', 'peak', 'frac_peak_over', 'inband']},
                  'dJ vs B1', tuple(round(x, 3) for x in e['dJ_vs_B1']), 'zero-net J', round(z['test']['J'], 3), e['gains'])
    return out


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'test'
    if cmd == 'train':
        train(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 150)
    elif cmd == 'test':
        test([int(x) for x in sys.argv[2:]] or [0])
