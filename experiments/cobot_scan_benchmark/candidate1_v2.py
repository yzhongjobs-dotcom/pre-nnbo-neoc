"""Candidate 1, v2: contact-loss supervisor for the moving-window estimator + OCP of candidate1_model.py.

Everything in candidate1_model.Candidate1 is kept (estimator, OCP, ordinary warm-started GN solve, tuned
configurations, k0 = 3e4 start prior of the post-test finding). v2 adds one supervisory layer, motivated by the
diverging fast-push traces (candidate1_v2_log.md, 'Failure mechanism'):

  Failure mechanism (T3, 0.25 s push ramp): a 7 mm push or retraction at up to 42-44 mm/s, often on bone, makes
  the window fit trade stiffness against drift (k_hat 26 -> 3 kN/m, vs at its 0.2 m/s bound) and the OCP
  over-retracts (command 15 mm back for a 7 mm push) -> contact loss. In loss, (i) the estimator keeps fitting
  the few pre-loss samples left in its window (vs at +-0.2 m/s), (ii) the local contact model is anchored at the
  current robot position with force level ~0, so the predictor believes the surface is touching the probe while
  it is 10-40 mm away, and (iii) the force integral winds up (ef = FD for the whole loss). The OCP then uses the
  feedforward to fight the integrator-driven approach it believes will hit the surface: u_ff runs to its -20 mm
  bound and stays there, the loss is prolonged (0.6-1.2 s), and the re-contact happens with a wound integrator,
  a saturated feedforward and an estimator window that mixes pre- and post-loss data -> 20-100 N impacts and a
  bounce limit cycle.

  Fix (v2 supervisor; the estimator and OCP are unchanged while the supervisor is in 'normal' mode). Final setting
  rec_D = 'fb' (R['sup_final'] in candidate1_v2_results.json); the other recovery laws are kept for the record:
   1. Simplex-style fallback on contact loss. Entered when the filtered force F < F_off (0.3 N final). While in it
      the model-based layer is suspended (none of its local-contact information is valid) and the cell's best fixed
      admittance B1 (gate_PC125.json, tuned on the same train split) is applied with K = 0: g = (M_B1, D_B1, 0,
      ki_app, KPF_B1, u_ff held). Left after n_back consecutive ticks with F in [band_lo, band_hi] (12 ticks in
      [1, 6] N final).
   2. Bumpless integrator (conditional integration through the gain interface, aw = 1). The plant force integral I
      keeps integrating ef ~ FD in loss; the supervisor replicates I exactly and accumulates the loss-time part in
      I_off; it then applies ki_app = clip(ki (I - I_off) / I, 0, ki) so that the applied integral action equals the
      virtual integral that was frozen during loss, and the OCP predictor is started from the applied integral.
      (In loss, ki_app = 0 when I <= 0, so the integral never drives the probe away.)
   3. Estimator freeze and reset: k held at the (conservatively released) k_used at loss entry, drift zeroed; after
      the fallback the window only uses samples taken after the re-contact tick, and until n_settle (12) such
      samples exist the stiffness stays at k_hold with zero drift and the measured force as level; the OCP restarts
      warm from the fallback gains (bumpless) and the held u_ff.
   4. Physical bound on the estimated surface speed |vs| <= vs_max (0.1 m/s final, 0.2 in v1).
  Also implemented and rejected on the fast-push dev split (candidate1_v2_log.md): rec_D='vapp' (damping set so the
  approach speed is v_app, with optional chase of the pre-loss drift estimate, stiffness-scaled speed c_app/k), rec_D=
  'rule' (rule gains at k_hold), a feedforward rate clamp du_max, a feedforward force floor F_lo, fallback entry on
  F > F_fb_hi, and settle-phase variants. Defaults of the class are not the final setting; sup = 0 reproduces
  Candidate1 exactly (checked).

Protocol (candidate1_v2_log.md): developed on train (1000+) / dev (2000+) splits only; fresh final test seed 4000.
Run (<= 2 processes): python -B candidate1_v2.py dev|tune|fresh|oldtest|known|etrtune|c1tune [cells]
"""
import math, time, json, sys, itertools
import numpy as np
import plant as P
import scenarios as S
import tune as U
import candidate1_model as C
from candidate1_model import _sync, _est, _solve_ord, XLO, XHI, LDMIN, LDMAX, FD

V2_DEFAULTS = dict(
    sup=1, F_off=0.5, F_on=1.0, rec_D='vapp', v_app=0.02, u_rate=0.0, n_settle=3,
    F_lo=0.0, vs_max=0.2, aw=1, du_max=0.0, chase=0.0, vs_lim=0.05, t_chase=0.15, v_min=0.002, c_app=0.0, hand='rule', settle_rec=0,
    fb=None, F_fb_hi=1e9, band_lo=1.0, band_hi=7.0, n_back=12,
)


class Candidate1v2(C.Candidate1):
    def __init__(self, **kw):
        v2 = {k: kw.pop(k) for k in list(kw) if k in V2_DEFAULTS}
        if kw.get('solver', 'ord') != 'ord': raise ValueError('v2 supports the final ordinary solver only')
        kw['solver'] = 'ord'
        super().__init__(**kw)
        self.v2 = dict(V2_DEFAULTS); self.v2.update(v2)
        self.__dict__.update(self.v2)

    def reset(self, pl):
        super().reset(pl)
        self.est_lo[2] = -self.vs_max; self.est_hi[2] = self.vs_max
        self.mode = 0; self.n_rc = 0; self.k_hold = self.kuse
        self.n_loss_entries = 0; self.n_loss_ticks = 0
        self.I_off = 0.; self.I_prev = 0.

    # ---------------------------------------------------------------- bumpless integrator (conditional integration)
    def _ki_sched(self, ef, loss_next):
        """the plant integral I keeps integrating ef ~ FD during contact loss. The supervisor keeps the replicated
        increments of I made while in loss in I_off (conditional integration of the virtual integral Iv = I - I_off,
        frozen in loss) and applies Iv through the gain: ki_app I_next = ki Iv_next, ki_app clipped to [0, ki] (the
        applied integral action is never larger than the virtual one). Identical to v1 until the first loss."""
        In = min(max(self.rs[2] + ef * self.Tc, -20.), 20.)
        if loss_next and In <= 0.: return 0.                 # no retracting integral drive during re-approach
        if not self.aw or (self.I_off == 0. and not loss_next): return self.ki
        Ivn = (self.rs[2] if loss_next else In) - self.I_off
        if abs(In) < 1e-9: return 0.
        return min(max(self.ki * Ivn / In, 0.), self.ki)

    # ---------------------------------------------------------------- recovery law (contact-loss mode)
    def _recover(self, n, log, F):
        M_, D_, K_, ki_, KP_, u_ = self.g
        ef = FD - F
        if self.mode == 0:                       # settle phase after re-contact (settle_rec): hold the gains
            kia = self._ki_sched(ef, False)
            self.g = (self.M, D_, self.K, kia, KP_, self.u_hold)
            return self.g
        kia = self._ki_sched(ef, True)
        if self.rec_D == 'fb':                   # Simplex fallback: the cell's best fixed admittance (K = 0)
            Mf, Df, kif, KPf = self.fb
            kia = min(kia, kif)
            uf = self.u_hold
            if self.u_rate > 0:
                uf = uf - math.copysign(min(self.u_rate * self.Tc, abs(uf)), uf)
            self.u_hold = uf
            self.g = (Mf, Df, 0., kia, KPf, uf)
            self.n_loss_ticks += 1
            if n % 4 == 0:
                self.hist.append((log[n - 1, 0], self.k_hold, self.k_hold, 0., Df, KPf, uf, float('nan')))
            return self.g
        if self.rec_D == 'vapp':
            In = min(max(self.rs[2] + ef * self.Tc, -20.), 20.)
            u = ef + kia * In - self.K * (self.rs[0] - self.rs[3])
            va = self.v_app
            if self.c_app > 0:                   # impact-limited approach speed: k_hold v T_react <= dF
                va = max(self.v_min, min(va, self.c_app / self.k_hold))
            if self.chase > 0 and (n - self.n_entry) * self.Tc < self.t_chase:
                # approach relative to the surface motion extrapolated from the pre-loss drift estimate
                va = max(self.v_min, va - self.chase * self.vs_pre)
            D = min(max(u / va, P.GMIN[1]), P.GMAX[1]) if u > 0 else P.GMAX[1]
        else:
            D = self.rule(self.k_hold)[0]
        uf = self.u_hold
        if self.u_rate > 0:
            step = self.u_rate * self.Tc
            uf = uf - math.copysign(min(step, abs(uf)), uf)
        self.u_hold = uf
        self.g = (self.M, D, self.K, kia, KP_, uf)
        self.n_loss_ticks += 1
        if n % 4 == 0:
            self.hist.append((log[n - 1, 0], self.k_hold, self.k_hold, 0., D, KP_, uf, float('nan')))
        return self.g

    # ---------------------------------------------------------------- per update
    def __call__(self, o):
        t0 = time.perf_counter()
        n = o.n; log = o.log
        if n == 0:
            D, KP = self.rule(self.kuse)
            self.x = np.array([math.log(D), KP / 1e-4, 0., 0.])
            self.dith_a = min(max(self.Fdith / self.kuse, self.amin), self.amax) if self.Fdith > 0 else 0.
            self.g = (self.M, D, self.K, self.ki, KP, self._dither(n))
            return self.g
        if self.rs is None:
            self.rs = np.array([log[0, 2], 0., 0., log[0, 2]]); self.m = 0
        M_, D_, K_, ki_, KP_, u_ = self.g
        _sync(log, self.m, n, self.rs, M_, D_, K_, ki_, self.Tc, self.VMAX)
        self.m = n
        if self.sup and self.aw and self.mode == 1: self.I_off += self.rs[2] - self.I_prev   # integral made in loss
        self.I_prev = self.rs[2]
        F = log[n - 1, 5]
        dtu = self.per * self.Tc
        # ---------------- supervisor
        if self.sup and self.rec_D == 'fb':
            if self.mode == 0 and (F < self.F_off or F > self.F_fb_hi):
                self.mode = 1; self.n_loss_entries += 1; self.n_entry = n; self.n_inband = 0
                self.k_hold = self.kuse; self.th[2] = 0.; self.th[3] = 0.
                self.u_hold = float(self.x[2]) * 1e-3
            elif self.mode == 1:
                self.n_inband = self.n_inband + 1 if self.band_lo <= F <= self.band_hi else 0
                if self.n_inband == 1: self.n_rc = n - 1
                if self.n_inband >= self.n_back:
                    self.mode = 0
                    self.th[0] = math.log(self.k_hold); self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
                    self.x = np.array([math.log(self.g[1]), self.g[4] / 1e-4, min(max(self.u_hold / 1e-3, XLO[2]), XHI[2]), 0.])
                    self.kuse = self.k_hold
            if self.mode == 1:
                g = self._recover(n, log, F)
                self.t_call.append(time.perf_counter() - t0)
                return g
        elif self.sup:
            if self.mode == 0 and F < self.F_off:
                self.mode = 1; self.n_loss_entries += 1; self.n_entry = n
                self.vs_pre = min(max(self.th[2], -self.vs_lim), self.vs_lim)
                self.k_hold = self.kuse
                self.th[2] = 0.; self.th[3] = 0.
                self.u_hold = float(self.x[2]) * 1e-3          # OCP part of u_ff (dither removed)
            elif self.mode == 1 and F > self.F_on:
                self.mode = 0; self.n_rc = n - 1
                self.th[0] = math.log(self.k_hold); self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
                Dr, KPr = self.rule(self.k_hold)
                if self.hand == 'cur': Dr, KPr = self.g[1], self.g[4]       # bumpless in the gains too
                self.x = np.array([math.log(Dr), KPr / 1e-4, min(max(self.u_hold / 1e-3, XLO[2]), XHI[2]), 0.])
                self.kuse = self.k_hold
            if self.mode == 1 or (self.settle_rec and self.n_rc > 0 and n - 1 - self.n_rc < self.n_settle):
                g = self._recover(n, log, F)
                self.t_call.append(time.perf_counter() - t0)
                return g
        kia = self._ki_sched(FD - F, False)
        # ---------------- estimator (window restricted to samples after the last re-contact)
        t1 = time.perf_counter()
        mu = self.th.copy()
        mu[2] = self.th[2] + (self.th[3] * dtu if self.nfree == 4 else 0.)
        pw = np.array([1.0 / (self.sk ** 2 * self.per), 0.0, 1.0 / (self.sv ** 2 * self.per), 1.0 / (self.sa ** 2 * self.per)])
        info = np.zeros(2)
        nw = self.Nw
        if self.sup and self.n_rc > 0: nw = min(self.Nw, n - 1 - self.n_rc)
        if self.sup and self.n_rc > 0 and nw < self.n_settle:
            nv = -1
            self.th[0] = mu[0]; self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
        else:
            nv = _est(log, n, nw, self.lam, self.ndel, self.eb, self.afc, self.ns, self.th, mu, pw, self.nfree,
                      self.est_lo, self.est_hi, self.iters_est, self.sz, info)
            if nv == 0:
                self.th[0] = mu[0]; self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
        kraw = math.exp(self.th[0])
        if self.klat > 0:
            self.kq.append((o.t[n - 1], kraw))
            while len(self.kq) > 1 and self.kq[1][0] <= o.t[n - 1] - self.klat + 1e-9: self.kq.pop(0)
            kraw = self.kq[0][1] if self.kq[0][0] <= o.t[n - 1] - self.klat + 1e-9 else self.k0
        if self.trel > 0:
            kraw = max(kraw, self.kuse * math.exp(-dtu / self.trel))
        self.kuse = kraw
        t2 = time.perf_counter()
        # ---------------- OCP (unchanged)
        nd = self.ndel
        tgt0 = log[n - 1 - nd, 11] if n - 1 - nd >= 0 else log[0, 2]
        q = np.array([log[j, 11] if j >= 0 else log[0, 2] for j in range(n - nd, n)], dtype=np.float64)
        kp = self.kuse
        # predictor integral state: the one actually applied at the next tick, ki * I0' = kia * I_next
        In_ = min(max(self.rs[2] + (FD - F) * self.Tc, -20.), 20.)
        I0 = (kia / self.ki * In_ - (FD - F) * self.Tc) if (self.ki > 0 and kia != self.ki) else self.rs[2]
        st = np.array([self.rs[0], self.rs[1], I0, self.rs[3], log[n - 1, 2], F, tgt0,
                       self.th[1], kp, self.th[2], self.th[3], log[n - 1, 2], 0., 0., 0., 0.])
        self.dith_a = min(max(self.Fdith / kp, self.amin), self.amax) if self.Fdith > 0 else 0.
        dith = np.array([self._dither(n + (i // self.per) * self.per) for i in range(self.H)])
        Dr, KPr = self.rule(kp)
        prior = np.array([math.log(Dr), KPr / 1e-4, 0., 0.])
        xprev = self.x.copy(); xprev[2] = xprev[2] + xprev[3] * 1e-2 * min(dtu, self.tu) / 1e-3
        xprev[2] = min(max(xprev[2], XLO[2]), XHI[2])
        lo = XLO.copy(); hi = XHI.copy()
        lo[0] = max(LDMIN, xprev[0] - self.rate_slow * self.per); hi[0] = min(LDMAX, xprev[0] + self.rate_fast * self.per)
        lo[1] = max(0., xprev[1] - self.rateK_fast * self.per); hi[1] = min(10., xprev[1] + self.rateK_slow * self.per)
        if self.ffmode == 0:
            lo[2] = hi[2] = 0.; lo[3] = hi[3] = 0.; xprev[2] = 0.; xprev[3] = 0.
        if self.F_lo > 0 and F < self.F_lo:      # feedforward force floor: no further retraction by u_ff
            lo[2] = max(lo[2], xprev[2]); lo[3] = max(lo[3], 0.)
        if self.du_max > 0:                      # feedforward rate clamp (|du_ff/dt| <= du_max)
            du = self.du_max * dtu / 1e-3
            lo[2] = max(lo[2], xprev[2] - du); hi[2] = min(hi[2], xprev[2] + du)
            if lo[2] > hi[2]: lo[2] = hi[2] = min(max(xprev[2], lo[2]), hi[2])
        x = np.minimum(np.maximum(xprev, lo), hi)
        nrinfo = np.zeros(4)
        c = _solve_ord(x, self.iters, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                       lo, hi, self.Jc, self.rprev, nrinfo)
        nroll = nrinfo[0]
        t3 = time.perf_counter()
        self.x = x
        self.n_upd += 1
        D = math.exp(x[0]); KP = x[1] * 1e-4; uff = x[2] * 1e-3
        self.g = (self.M, D, self.K, kia, KP, uff + dith[0])
        tend = time.perf_counter()
        self.t_est.append(t2 - t1); self.t_ocp.append(t3 - t2); self.t_call.append(tend - t0); self.nroll.append(nroll)
        if n % 4 == 0:
            self.hist.append((o.t[n - 1], math.exp(self.th[0]), kp, self.th[2], D, KP, uff, c))
        return self.g


# ============================================================================ evaluation (<= 2 processes)
NPROC = 2
_POOL = None


def warmup():
    C.warmup()
    p, T = S.make('N1', np.random.default_rng(0), dict(E=2e4, A=0.001))
    P.episode(p, Candidate1v2(), T=0.3, seed=0, update=0.008)


def pool():
    global _POOL
    if _POOL is None:
        from multiprocessing import Pool
        _POOL = Pool(NPROC, initializer=warmup)
    return _POOL


def close_pool():
    global _POOL
    if _POOL is not None:
        _POOL.close(); _POOL.join(); _POOL = None


def _job(args):
    p, T, spec, seed, keep = args
    kind, kw = spec
    if kind != 'C1v2':
        return C._job(args)
    pol = Candidate1v2(**kw)
    log = P.episode(p, pol, T=T, platform='PC125', seed=seed, update=pol.update)
    m = S.metrics(log)
    tc = np.array(pol.t_call) if pol.t_call else np.zeros(1)
    m['t_call_mean'] = float(tc.mean()); m['t_call_p99'] = float(np.percentile(tc, 99)); m['t_call_max'] = float(tc.max())
    m['n_loss_entries'] = pol.n_loss_entries; m['n_loss_ticks'] = pol.n_loss_ticks
    if keep:
        m['hist'] = np.array(pol.hist); m['log'] = log
    return m


def run(scenes, specs, seed0, keep=False, procs=NPROC):
    jobs = [(p, T, s, seed0 + i, keep) for i, (p, T) in enumerate(scenes) for s in specs]
    if procs > 1 and len(jobs) > 1:
        res = pool().map(_job, jobs, chunksize=max(1, len(jobs) // (8 * procs)))
    else:
        res = [_job(j) for j in jobs]
    ns = len(specs)
    return [res[i * ns:(i + 1) * ns] for i in range(len(scenes))]


def summ(rows):
    out = U.summary(rows)
    Js = np.array([r['J'] for r in rows])
    out['J_median'] = float(np.median(Js)); out['n_J_over_20'] = int((Js > 20).sum())
    out['n_runaway'] = int(sum(1 for r in rows if (r['J'] > 20) or r.get('diverged', 0)))
    out['frac_scenes_with_loss'] = float(np.mean([r['loss'] > 0 for r in rows]))
    for k in ['n_loss_entries', 't_call_mean']:
        if k in rows[0]: out[k] = float(np.mean([r[k] for r in rows]))
    return out


def line(name, s):
    return (f"{name:28s} J {s['J']:7.3f} (med {s['J_median']:.3f}) rmse {s['rmse']:.3f} >7N {s['frac_peak_over']:.3f} "
            f"peak {s['peak']:.2f} loss {s['loss']:.4f} J>20 {s['n_J_over_20']} div {s['diverged']:.0f}")


# ============================================================================ protocol
R1 = json.load(open('candidate1_results.json'))
CFG1 = {c: dict(R1['tuning'][c]['C1']['config']) for c in ['T3', 'T1', 'T3fast']}     # tuned v1 configs (k0 = 5000)
CELLS = {'T3': ('T3', dict(ramp=0.5, A=0.001)), 'T1': ('T1', dict(ratio=30., A=0.001)),
         'T3fast': ('T3', dict(ramp=0.25, A=0.001)),
         'T3fast4': ('T3', dict(ramp=0.25, A=0.004)), 'T3_4': ('T3', dict(ramp=0.5, A=0.004)),
         'T1_4': ('T1', dict(ratio=30., A=0.004))}
TR0, DEV0, TE0, FRESH0 = 1000, 2000, 3000, 4000
NFRESH = 40


def v2cfg(cell, fix):
    """the v1 configuration tuned for the cell (or its 1 mm counterpart / own v2 tuning when present), with the v2
    start prior k0 = 3e4 and the supervisor setting `fix`."""
    base = CFG1.get(cell)
    if base is None:
        base = R.get('c1_tuning', {}).get(cell, {}).get('config')
    if base is None:
        base = CFG1[{'T3fast4': 'T3fast', 'T3_4': 'T3', 'T1_4': 'T1'}[cell]]
    out = dict(base, k0=3e4, **fix)
    if out.get('rec_D') == 'fb' and out.get('fb') is None:
        g = gate_ref(cell)['B1']['gains']; out['fb'] = (g[0], g[1], g[3], g[4])
    return out


def gate_ref(cell):
    task, dose = CELLS[cell]
    for c in json.load(open('gate_PC125.json')):
        if c['task'] == task and c['dose'] == dose: return c


def etr_spec(cell):
    """EtR-A: the tuning of adaptive_results.json for the two original cells, else this module's own tuning
    (same procedure, baselines_adaptive.tune_etr)."""
    task, dose = CELLS[cell]
    for c in json.load(open('adaptive_results.json'))['cells']:
        if c['task'] == task and c['dose'] == dose:
            t = c['methods']['EtR-A']['tuning']
            return ('EtR', (tuple(t['rule_params']), dict(t['est_params'], rule='truthk'))), 'adaptive_results.json'
    t = R.get('etr_tuning', {}).get(cell)
    if t is None: return None, None
    return ('EtR', (tuple(t['rule_params']), dict(t['est_params'], rule='truthk'))), 'candidate1_v2 etrtune (train)'


def etr_tune(cell, log=print):
    """baselines_adaptive.tune_etr re-implemented on this module's 2-process pool: stage 1 rule grid (gate.RULE) with
    the default lock-in settings, stage 2 the 8 estimator variants on the 3 best rules; train split, 20 scenes."""
    import baselines_adaptive as BA
    task, dose = CELLS[cell]; train = S.cohort(task, 20, TR0, dose)
    E1 = dict(BA.EST0)
    Q = BA.RULE
    sp = [('EtR', (tuple(q), dict(E1, rule='truthk'))) for q in Q]
    r = run(train, sp, TR0); j, Jm, f = C.select(r)
    order = np.argsort(np.where(f, Jm, Jm + 1e6))
    tops = [Q[i] for i in order[:3]]
    cand = [(q, e) for q in tops for e in BA.EST_VAR]
    r2 = run(train, [('EtR', (tuple(q), dict(e, rule='truthk'))) for q, e in cand], TR0); j2, Jm2, f2 = C.select(r2)
    if (not f2[j2] and f[j]) or (f[j] == f2[j2] and Jm[j] <= Jm2[j2]): best = (Q[j], E1, float(Jm[j]), bool(f[j]))
    else: best = (cand[j2][0], cand[j2][1], float(Jm2[j2]), bool(f2[j2]))
    out = dict(rule_params=list(best[0]), est_params=best[1], train_J=best[2], train_feasible=best[3], n_configs=len(Q) + len(cand))
    log(f"[etrtune] {cell}: {out}")
    return out


def compare(scenes, seed0, specs, names, ref='B1'):
    res = run(scenes, specs, seed0)
    out = {}
    per = {nm: [r[j] for r in res] for j, nm in enumerate(names)}
    for nm in names:
        rows = per[nm]
        d = dict(summary=summ(rows), J_per_scene=[float(r['J']) for r in rows], rmse_per_scene=[float(r['rmse']) for r in rows],
                 peak_per_scene=[float(r['peak']) for r in rows], loss_per_scene=[float(r['loss']) for r in rows])
        out[nm] = d
    return out, per


def cis(out, per, cand_names, refs):
    for nm in cand_names:
        for rf in refs:
            if rf in per and rf != nm:
                out[nm]['dJ_vs_' + rf] = U.paired_ci([r['J'] for r in per[nm]], [r['J'] for r in per[rf]])
                out[nm]['drmse_vs_' + rf] = U.paired_ci([r['rmse'] for r in per[nm]], [r['rmse'] for r in per[rf]])
                out[nm]['dover7_vs_' + rf] = U.paired_ci([float(r['peak'] > 7.0) for r in per[nm]],
                                                         [float(r['peak'] > 7.0) for r in per[rf]])
                out[nm]['wins_vs_' + rf] = int(sum(a['J'] < b['J'] for a, b in zip(per[nm], per[rf])))


def baselines(cell):
    g = gate_ref(cell)
    sp = [('fixed', tuple(g['B1']['gains'])), ('B3b', (0.1,) + tuple(g['B3b']['0.1']['rule']))]
    nm = ['B1', 'B3b_0.1']
    e, src = etr_spec(cell)
    if e is not None: sp.append(e); nm.append('EtR-A')
    return sp, nm, src


def evaluate_cell(cell, seed0, n, cfgs, log=print, with_base=True):
    task, dose = CELLS[cell]; scenes = S.cohort(task, n, seed0, dose)
    names = list(cfgs); specs = [('C1v2', cfgs[k]) for k in names]
    src = None
    if with_base:
        bs, bn, src = baselines(cell); specs += bs; names += bn
    out, per = compare(scenes, seed0, specs, names)
    cis(out, per, list(cfgs), ['B1', 'EtR-A', 'B3b_0.1'] + [k for k in cfgs if k.startswith('C1v1')])
    for nm in names:
        s = out[nm]['summary']
        extra = ''
        if 'dJ_vs_B1' in out[nm]:
            a = out[nm]['dJ_vs_B1']; extra = f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f}) wins {out[nm]['wins_vs_B1']}/{n}"
        log(f"[{cell} seed {seed0} n {n}] " + line(nm, s) + extra)
    return dict(cell=cell, task=task, dose=dose, seed0=seed0, n=n, configs=cfgs, etr_source=src, methods=out)


def _save():
    json.dump(R, open('candidate1_v2_results.json', 'w'), indent=1, default=float)


try:
    R = json.load(open('candidate1_v2_results.json'))
except Exception:
    R = {}

# --------------------------------------------------------------- declared settings
FIX0 = dict(sup=0)                                                   # = Candidate1 (v1) with k0 = 3e4
FIX_FINAL = None                                                     # set by 'tune' (stored in the JSON)

# declared search over the supervisor parameters (49 configurations: default + 32 random + 16 one-parameter
# perturbations of the stage-A best, same style as candidate1_model.SPACE), fast-push train split, v1 T3fast config
SUP_SPACE = dict(F_off=[0.3, 0.5, 1.0], band_lo=[0.5, 1.0, 2.0], band_hi=[6.0, 7.0, 8.0], n_back=[6, 12, 24, 48],
                 n_settle=[3, 6, 12], u_rate=[0.0, 0.01, 0.03], aw=[0, 1], vs_max=[0.06, 0.1, 0.2], F_fb_hi=[10.0, 1e9])
SUP_DEFAULT = dict(sup=1, rec_D='fb', F_off=0.5, band_lo=1.0, band_hi=7.0, n_back=12, n_settle=6, u_rate=0.0, aw=1,
                   vs_max=0.2, F_fb_hi=1e9)


def sup_draws():
    rng = np.random.default_rng(4242)
    A = [dict(SUP_DEFAULT)]
    for _ in range(32):
        c = dict(SUP_DEFAULT); c.update({k: v[int(rng.integers(len(v)))] for k, v in SUP_SPACE.items()})
        A.append(c)
    return A


def sup_perturb(best):
    rng = np.random.default_rng(777); out = []; keys = list(SUP_SPACE)
    while len(out) < 16:
        k = keys[int(rng.integers(len(keys)))]
        alt = [v for v in SUP_SPACE[k] if v != best.get(k)]
        c = dict(best); c[k] = alt[int(rng.integers(len(alt)))]
        if c not in out and c != best: out.append(c)
    return out


def tune_sup(log=print, ntr=40):
    cell = 'T3fast'; task, dose = CELLS[cell]; train = S.cohort(task, ntr, TR0, dose)
    A = sup_draws(); t0 = time.time()
    rA = run(train, [('C1v2', v2cfg(cell, c)) for c in A + [FIX0]], TR0)
    jv1 = float(np.mean([r[-1]['J'] for r in rA])); rA = [r[:-1] for r in rA]
    jA, JmA, fA = C.select(rA)
    B = sup_perturb(A[jA])
    rB = run(train, [('C1v2', v2cfg(cell, c)) for c in B], TR0); _, JmB, fB = C.select(rB)
    cand = A + B; Jall = np.concatenate([JmA, JmB]); fall = np.concatenate([fA, fB])
    # selection: the declared C1 rule (min mean train J among peak-feasible settings, plain min J if none is feasible)
    j = int(np.argmin(np.where(fall, Jall, Jall + 1e6)))
    order = np.argsort(np.where(fall, Jall, Jall + 1e6))
    log(f"[tune_sup] v1 (sup=0, k0=3e4) train J {jv1:.3f}")
    log(f"[tune_sup] T3fast train n={ntr}: best J {Jall[j]:.3f} feasible {bool(fall[j])} ({len(cand)} configs, {int(fall.sum())} feasible, "
        f"{time.time() - t0:.0f}s) {cand[j]}")
    for i in order[:6]: log(f"   {Jall[i]:.3f} {cand[i]}")
    return dict(best=cand[j], train_J=float(Jall[j]), feasible=bool(fall[j]), n_configs=len(cand), n_train=ntr,
                top6=[(float(Jall[i]), cand[i]) for i in order[:6]], all=[(float(Jall[i]), cand[i]) for i in range(len(cand))],
                v1_train_J=jv1, seconds=time.time() - t0)


def c1_tune(cell, fix, log=print):
    """the declared candidate1_model 49-configuration tuning (SPACE, same draws) for a cell without a v1 tuning
    (4 mm cells), run with the v2 supervisor setting `fix` and k0 = 3e4."""
    task, dose = CELLS[cell]; train = S.cohort(task, C.NTR, TR0, dose)
    rng = np.random.default_rng(12345)
    A = [dict()] + C.stage_a(rng); fixed = dict(solver='ord', iters=1, k0=3e4, **fix)
    if fixed.get('rec_D') == 'fb':
        g = gate_ref(cell)['B1']['gains']; fixed['fb'] = (g[0], g[1], g[3], g[4])
    t0 = time.time()
    rA = run(train, [('C1v2', dict(c, **fixed)) for c in A], TR0); jA, JmA, fA = C.select(rA)
    B = C.stage_b(A[jA], np.random.default_rng(777))
    rB = run(train, [('C1v2', dict(c, **fixed)) for c in B], TR0); _, JmB, fB = C.select(rB)
    cand = A + B; Jall = np.concatenate([JmA, JmB]); fall = np.concatenate([fA, fB])
    j = int(np.argmin(np.where(fall, Jall, Jall + 1e6)))
    cfg = dict(cand[j], solver='ord', iters=1)          # v2cfg adds k0, the supervisor setting and fb
    log(f"[c1tune] {cell}: best train J {Jall[j]:.3f} feasible {bool(fall[j])} ({len(cand)} configs, {int(fall.sum())} feasible, "
        f"{time.time() - t0:.0f}s) {cand[j]}")
    return dict(config=cfg, train_J=float(Jall[j]), train_feasible=bool(fall[j]), n_configs=len(cand), seconds=time.time() - t0)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'dev'
    arg = sys.argv[2] if len(sys.argv) > 2 else ''
    L = []
    def log(s):
        print(s, flush=True); L.append(s)
    fixf = R.get('sup_final')
    if what == 'dev':
        # named settings on train/dev of the given cells; spec: cells|name=json;name=json
        cells = arg.split(',') if arg else ['T3fast']
        sets = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {'v1': FIX0, 'v2default': dict(V2_DEFAULTS)}
        split = sys.argv[4] if len(sys.argv) > 4 else 'dev'
        nsc = int(sys.argv[5]) if len(sys.argv) > 5 else 40
        s0 = {'train': TR0, 'dev': DEV0}[split]
        for c in cells:
            r = evaluate_cell(c, s0, nsc, {k: v2cfg(c, v) for k, v in sets.items()}, log, with_base=len(sys.argv) > 6)
            R.setdefault('development', {}).setdefault(c, {})[f'{split}{nsc}:' + ','.join(sets)] = r; _save()
    elif what == 'tune':
        R['sup_tuning'] = tune_sup(log); _save()
    elif what == 'etrtune':
        for c in arg.split(','):
            R.setdefault('etr_tuning', {})[c] = etr_tune(c, log); _save()
    elif what == 'c1tune':
        for c in arg.split(','):
            R.setdefault('c1_tuning', {})[c] = c1_tune(c, fixf, log); _save()
    elif what in ('fresh', 'oldtest'):
        s0 = FRESH0 if what == 'fresh' else TE0
        for c in arg.split(','):
            cf = {'C1v2': v2cfg(c, fixf)}
            if c in CFG1: cf['C1v1'] = dict(CFG1[c], sup=0)             # v1 as tested (k0 = 5000, no supervisor)
            r = evaluate_cell(c, s0, NFRESH, cf, log)
            R.setdefault(what, {})[c] = r; _save()
    elif what == 'known':
        c = arg or 'T3'
        cfg = v2cfg(c, fixf); b1 = gate_ref(c)['B1']['gains']
        rng = np.random.default_rng(7); scenes = []
        for E in [1.5e4, 1.0e6]:
            for k in range(3): scenes.append(S.make('N1', rng, dict(E=E, A=0.001)))
        res = run(scenes, [('C1v2', cfg), ('fixed', tuple(b1))], 500, keep=True)
        out = {'soft_15kPa': [], 'stiff_1MPa': []}
        for (p, T), r in zip(scenes, res):
            m, mb = r[0], r[1]; h = m['hist']; lg = m['log']
            kt = lg[::4, 8][:len(h)]; ok = (h[:, 0] > 1.0) & (kt > 0)
            d = dict(E=float(p[0]), k_true_median=float(np.median(kt[ok])), k_hat_median=float(np.median(h[ok, 1])),
                     k_hat_over_true_median=float(np.median(h[ok, 1] / kt[ok])), D_median=float(np.median(h[ok, 4])),
                     KPF_median=float(np.median(h[ok, 5])), J=m['J'], rmse=m['rmse'], peak=m['peak'], loss=m['loss'],
                     n_loss_entries=m['n_loss_entries'], J_B1=mb['J'], rmse_B1=mb['rmse'], peak_B1=mb['peak'])
            out['soft_15kPa' if p[0] < 1e5 else 'stiff_1MPa'].append(d)
        for k, Ls in out.items():
            log(f"[known {c}] {k}: k_true {np.mean([d['k_true_median'] for d in Ls]):.0f} k_hat {np.mean([d['k_hat_median'] for d in Ls]):.0f} "
                f"D {np.mean([d['D_median'] for d in Ls]):.0f} KPF {np.mean([d['KPF_median'] for d in Ls]):.2e} "
                f"J {np.mean([d['J'] for d in Ls]):.4f} vs B1 {np.mean([d['J_B1'] for d in Ls]):.4f} "
                f"rmse {np.mean([d['rmse'] for d in Ls]):.4f} vs {np.mean([d['rmse_B1'] for d in Ls]):.4f} "
                f"peak max {np.max([d['peak'] for d in Ls]):.2f} loss entries {np.sum([d['n_loss_entries'] for d in Ls])}")
        R.setdefault('known_behaviour', {})[c] = out; _save()
    R.setdefault('console', []).extend(L); _save()
    close_pool()
