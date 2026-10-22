"""Candidate 1, v3: configuration fix for the T3 0.5 s / 4 mm negative, plus a single shared configuration, for the
estimator + OCP + supervisor of candidate1_model.py / candidate1_v2.py (neither file is changed; v3 imports them).

Mechanism of the v2 negative on T3 0.5 s / 4 mm (candidate1_v3_log.md, 'Mechanism'):
  The 4 mm cell had no v1 tuning, so v2 ran the declared 49-configuration C1 search on 20 train scenes for it. Only 8
  of 49 settings met the peak constraint (<= 1 of 20 scenes above 7 N), so the choice was decided by one or two
  scenes, and the search space did not contain the transferable 1 mm configuration (it differs from the hand-set
  default in probe_v only). The selected setting has no excitation dither (Fdith = 0) and a 2x faster stiffness
  random walk (sk = 0.3); with that, the window fit cannot separate stiffness from surface motion while the probe
  follows the breathing surface (robot motion is nearly collinear with the drift basis {1, tau, tau^2}), and k_hat
  drifts high (k_hat/k_true median 2.1-2.5 on soft tissue and 1.4 on bone in quiet periods, p90 4-5, on train, dev
  and the seed-4000 traces; 1.1 / 1.05 with the 1 mm configuration). The rule prior then gives too-stiff gains; at a
  push onset the stiffness spike plus conservative release sends D to ~1500 N s/m and KPF to ~1e-5, the probe cannot
  follow, the force overshoots to ~10 N, and the long feedforward ramp (tu 0.08) then over-retracts into contact
  loss. The same setting is as bad on the 1 mm cell (dev J 0.58 vs 0.12), so the negative is a selection failure,
  not an effect of 4 mm breathing: the constant-velocity + acceleration drift model fits 0.2-0.33 Hz breathing to
  ~1e-8 m over the 0.1 s window (drift-velocity error on bone +-0.3 mm/s at both amplitudes); no periodic drift
  term is needed.

v3 final design (developed on train 1000+ / dev 2000+ only; fresh test seed 5000, run once):
  C1v3: no per-cell re-tuning for breathing amplitude. Each 4 mm cell uses the v1 configuration tuned on its 1 mm
     counterpart (T3_4 <- T3, T1_4 <- T1, T3fast4 <- T3fast, the last one unchanged from v2); estimator, OCP, v2
     supervisor, k0 = 3e4 and the cell's B1 fallback are unchanged. Class defaults below reproduce v2 exactly (checked).
  C1v3-shared (reported separately): one estimator/OCP configuration for all six cells from a declared pooled search
     (shared_tune: 3 anchors + 32 random + 16 perturbations of C.SPACE, 40 train scenes per cell, score = mean over
     cells of log(J / J_B1), peak-feasible in the four non-fast cells). Fallback = the cell's B1 gains.

Implemented, tested on train/dev and rejected (defaults off; numbers in candidate1_v3_log.md):
  dmode = 1/2  excitation-aware dither: amplitude scaled by s = clip((r - r_lo) / (r_hi - r_lo), 0, 1), r = std of log k of
               the window fit (dmode 1: posterior incl. the arrival prior, relative to sk sqrt(per); dmode 2: data only).
               r stays at 0.92-0.99 (dmode 1) and the data-only std is 0.34-0.47 even with the dither, so the window never
               becomes confident enough to switch the excitation off (dmode 2 was measured open loop only, not tuned).
  fbk = 1..4   stiffness-aware fallback: D_fb = max(D_B1, k_hold / wI), KPF_fb = min(KPF_B1, cP / k_hold) (1: always;
               2: KPF only; 3: only on ticks with F >= band_lo; 4: latched from the first re-contact tick), k_hold = max
               k_used over the last kmem s. All slow the in-loss approach or remove the fast re-contact path and raise
               fast-push J.
  KPFcap       cap of the rule prior and OCP box on KPF (soft-side aggressiveness): <= 0.004 N RMSE on soft N1.

Run (<= 2 processes): python -B candidate1_v3.py shared_tune | dev [cells] | final [cells] | known
"""
import os
for _k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ.setdefault(_k, '1')
import math, time, json, sys
import numpy as np
from numba import njit
import plant as P
import scenarios as S
import tune as U
import candidate1_model as C
import candidate1_v2 as V
from candidate1_model import _sync, _est, _solve_ord, XLO, XHI, LDMIN, LDMAX, FD

V3_DEFAULTS = dict(dmode=0, r_lo=0.3, r_hi=0.7, KPFcap=1e-3, fbk=0, kmem=0.0, rec_sig=0)


@njit(cache=False)
def _kvar(log, nend, Nw, lam, ndel, eb, afc, ns, th, pw, nfree, sz, kprior):
    """posterior variance of kappa = log k of the window fit of candidate1_model._est at the estimate th: inverse of
    the Gauss-Newton normal matrix (same data, weights and arrival-prior weights pw; parameters >= nfree fixed)."""
    j0 = max(ndel + 1, nend - 1 - Nw); j1 = nend - 2
    m = j1 - j0 + 1
    if m < 3: return -1.0
    te = log[nend - 1, 0]; ye = log[nend - 1, 2]
    a1 = 1.0 - afc; A8 = a1 ** ns; Ws = 1.0 - A8
    A = np.zeros((4, 4)); Jr = np.zeros(4)
    k = math.exp(th[0]); nv = 0
    for jj in range(m):
        j = j0 + jj
        Fj = log[j, 5]; Fj1 = log[j + 1, 5]
        if not ((Fj > 0.5) and (Fj1 > 0.5)): continue
        tgt = log[j - ndel, 11]; yy = log[j, 2]; tj = log[j, 0]
        sY = 0.; sT = 0.; sT2 = 0.
        for mm in range(1, ns + 1):
            yy = tgt + (yy - tgt) * eb
            w = afc * a1 ** (ns - mm)
            tau = tj + mm * 1e-3 - te
            sY += w * (yy - ye); sT += w * tau; sT2 += w * 0.5 * tau * tau
        Jr[0] = k * (sY + th[2] * sT + th[3] * sT2) / sz
        Jr[1] = Ws / sz; Jr[2] = k * sT / sz; Jr[3] = k * sT2 / sz
        w = lam ** (j1 - j); nv += 1
        for p in range(4):
            for q in range(4): A[p, q] += w * Jr[p] * Jr[q]
    if nv < 3: return -1.0
    for p in range(4): A[p, p] += (pw[p] if (p > 0 or kprior) else 0.) + 1e-9
    Af = A[:nfree, :nfree].copy()
    Ai = np.linalg.inv(Af)
    return Ai[0, 0]


class Candidate1v3(V.Candidate1v2):
    def __init__(self, **kw):
        v3 = {k: kw.pop(k) for k in list(kw) if k in V3_DEFAULTS}
        super().__init__(**kw)
        self.v3 = dict(V3_DEFAULTS); self.v3.update(v3)
        self.__dict__.update(self.v3)

    def reset(self, pl):
        super().reset(pl)
        self.fb0 = self.fb
        self.kmem_q = []
        self.r_hist = []; self.s_hist = []
        self.dscale = 1.0

    def rule(self, k):
        D, KP = super().rule(k)
        return D, min(KP, self.KPFcap)

    def _enter_fallback(self, n):
        k_hold = self.kuse
        if self.kmem > 0 and self.kmem_q:
            k_hold = max(k_hold, max(q[1] for q in self.kmem_q))
        self.k_hold_fb = k_hold
        self.fb = self.fb0
        if self.fbk and self.fb0 is not None:
            Mf, Df, kif, KPf = self.fb0
            Dk = max(Df, min(k_hold / self.wI, P.GMAX[1])); KPk = min(KPf, self.cP / k_hold)
            if self.fbk == 1: self.fb_c = (Mf, Dk, kif, KPk)          # D and KPF, always
            elif self.fbk == 2: self.fb_c = (Mf, Df, kif, KPk)        # KPF only (proportional loop gain), always
            else: self.fb_c = (Mf, Dk, kif, KPk)                      # 3: D and KPF only while in contact
            if self.fbk in (1, 2): self.fb = self.fb_c

    def _recover(self, n, log, F):
        if self.fbk == 3 and self.mode == 1:
            self.fb = self.fb_c if F >= self.band_lo else self.fb0
        if self.fbk == 4 and self.mode == 1 and F >= self.band_lo:
            self.fb = self.fb_c                              # 4: latched once contact is regained (approach on B1)
        return super()._recover(n, log, F)

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
        if self.sup and self.aw and self.mode == 1: self.I_off += self.rs[2] - self.I_prev
        self.I_prev = self.rs[2]
        F = log[n - 1, 5]
        dtu = self.per * self.Tc
        # ---------------- supervisor (v2, final 'fb' law; v3: stiffness-aware fallback gains at entry)
        if self.sup:
            if self.rec_D != 'fb': raise ValueError('v3 supports the v2 final fallback law only')
            if self.mode == 0 and (F < self.F_off or F > self.F_fb_hi):
                self.mode = 1; self.n_loss_entries += 1; self.n_entry = n; self.n_inband = 0
                self.k_hold = self.kuse; self.th[2] = 0.; self.th[3] = 0.
                self.u_hold = float(self.x[2]) * 1e-3
                self._enter_fallback(n)
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
        kia = self._ki_sched(FD - F, False)
        # ---------------- estimator (window restricted to samples after the last re-contact)
        t1 = time.perf_counter()
        mu = self.th.copy()
        mu[2] = self.th[2] + (self.th[3] * dtu if self.nfree == 4 else 0.)
        pw = np.array([1.0 / (self.sk ** 2 * self.per), 0.0, 1.0 / (self.sv ** 2 * self.per), 1.0 / (self.sa ** 2 * self.per)])
        info = np.zeros(2)
        nw = self.Nw
        if self.sup and self.n_rc > 0: nw = min(self.Nw, n - 1 - self.n_rc)
        rr = 1.0
        if self.sup and self.n_rc > 0 and nw < self.n_settle:
            nv = -1
            self.th[0] = mu[0]; self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
        else:
            nv = _est(log, n, nw, self.lam, self.ndel, self.eb, self.afc, self.ns, self.th, mu, pw, self.nfree,
                      self.est_lo, self.est_hi, self.iters_est, self.sz, info)
            if nv == 0:
                self.th[0] = mu[0]; self.th[1] = F; self.th[2] = 0.; self.th[3] = 0.
            elif self.dmode or self.rec_sig:
                vk = _kvar(log, n, nw, self.lam, self.ndel, self.eb, self.afc, self.ns, self.th, pw, self.nfree, self.sz,
                           1 if self.dmode == 1 else 0)
                if vk > 0: rr = math.sqrt(vk) / ((self.sk * math.sqrt(self.per)) if self.dmode == 1 else 1.0)
                else: rr = 1e3
        kraw = math.exp(self.th[0])
        if self.klat > 0:
            self.kq.append((o.t[n - 1], kraw))
            while len(self.kq) > 1 and self.kq[1][0] <= o.t[n - 1] - self.klat + 1e-9: self.kq.pop(0)
            kraw = self.kq[0][1] if self.kq[0][0] <= o.t[n - 1] - self.klat + 1e-9 else self.k0
        if self.trel > 0:
            kraw = max(kraw, self.kuse * math.exp(-dtu / self.trel))
        self.kuse = kraw
        if self.kmem > 0:
            self.kmem_q.append((o.t[n - 1], kraw))
            while self.kmem_q and self.kmem_q[0][0] < o.t[n - 1] - self.kmem: self.kmem_q.pop(0)
        t2 = time.perf_counter()
        # ---------------- OCP (unchanged except the optional KPF cap)
        nd = self.ndel
        tgt0 = log[n - 1 - nd, 11] if n - 1 - nd >= 0 else log[0, 2]
        q = np.array([log[j, 11] if j >= 0 else log[0, 2] for j in range(n - nd, n)], dtype=np.float64)
        kp = self.kuse
        In_ = min(max(self.rs[2] + (FD - F) * self.Tc, -20.), 20.)
        I0 = (kia / self.ki * In_ - (FD - F) * self.Tc) if (self.ki > 0 and kia != self.ki) else self.rs[2]
        st = np.array([self.rs[0], self.rs[1], I0, self.rs[3], log[n - 1, 2], F, tgt0,
                       self.th[1], kp, self.th[2], self.th[3], log[n - 1, 2], 0., 0., 0., 0.])
        if self.dmode:
            self.dscale = min(max((rr - self.r_lo) / max(self.r_hi - self.r_lo, 1e-9), 0.), 1.)
        if self.rec_sig and n % 4 == 0:
            self.r_hist.append((o.t[n - 1], rr, self.dscale))
        self.dith_a = (self.dscale * min(max(self.Fdith / kp, self.amin), self.amax)) if self.Fdith > 0 else 0.
        dith = np.array([self._dither(n + (i // self.per) * self.per) for i in range(self.H)])
        Dr, KPr = self.rule(kp)
        prior = np.array([math.log(Dr), KPr / 1e-4, 0., 0.])
        xprev = self.x.copy(); xprev[2] = xprev[2] + xprev[3] * 1e-2 * min(dtu, self.tu) / 1e-3
        xprev[2] = min(max(xprev[2], XLO[2]), XHI[2])
        lo = XLO.copy(); hi = XHI.copy()
        lo[0] = max(LDMIN, xprev[0] - self.rate_slow * self.per); hi[0] = min(LDMAX, xprev[0] + self.rate_fast * self.per)
        lo[1] = max(0., xprev[1] - self.rateK_fast * self.per); hi[1] = min(10., xprev[1] + self.rateK_slow * self.per)
        hi[1] = min(hi[1], self.KPFcap / 1e-4)
        if lo[1] > hi[1]: lo[1] = hi[1]
        if self.ffmode == 0:
            lo[2] = hi[2] = 0.; lo[3] = hi[3] = 0.; xprev[2] = 0.; xprev[3] = 0.
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
    V.warmup()
    p, T = S.make('N1', np.random.default_rng(0), dict(E=2e4, A=0.001))
    P.episode(p, Candidate1v3(dmode=1, fbk=1, kmem=0.2, rec_sig=1, sup=1, rec_D='fb', fb=(1., 150., 30., 1e-4)),
              T=0.3, seed=0, update=0.008)


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
    if kind != 'C1v3':
        return V._job(args)
    pol = Candidate1v3(**kw)
    log = P.episode(p, pol, T=T, platform='PC125', seed=seed, update=pol.update)
    m = S.metrics(log)
    tc = np.array(pol.t_call) if pol.t_call else np.zeros(1)
    m['t_call_mean'] = float(tc.mean()); m['t_call_p99'] = float(np.percentile(tc, 99)); m['t_call_max'] = float(tc.max())
    m['n_loss_entries'] = pol.n_loss_entries; m['n_loss_ticks'] = pol.n_loss_ticks
    if pol.r_hist:
        rh = np.array(pol.r_hist); m['r_median'] = float(np.median(rh[:, 1])); m['dscale_mean'] = float(np.mean(rh[:, 2]))
    if keep:
        m['hist'] = np.array(pol.hist); m['log'] = log
        if pol.r_hist: m['r_hist'] = np.array(pol.r_hist)
    return m


def run(scenes, specs, seed0, keep=False, procs=NPROC):
    jobs = [(p, T, s, seed0 + i, keep) for i, (p, T) in enumerate(scenes) for s in specs]
    if procs > 1 and len(jobs) > 1:
        res = pool().map(_job, jobs, chunksize=max(1, len(jobs) // (8 * procs)))
    else:
        res = [_job(j) for j in jobs]
    ns = len(specs)
    return [res[i * ns:(i + 1) * ns] for i in range(len(scenes))]


# ============================================================================ protocol
CELLS = V.CELLS                                   # T3, T1, T3fast, T3fast4, T3_4, T1_4
ORDER = ['T3', 'T3_4', 'T3fast', 'T3fast4', 'T1', 'T1_4']
COUNTERPART = {'T3_4': 'T3', 'T1_4': 'T1', 'T3fast4': 'T3fast', 'T3': 'T3', 'T1': 'T1', 'T3fast': 'T3fast'}
TR0, DEV0, FRESH5 = 1000, 2000, 5000
NFINAL = 40
SUP = V.R['sup_final']                            # v2 supervisor, unchanged


def fb_of(cell):
    g = V.gate_ref(cell)['B1']['gains']; return (g[0], g[1], g[3], g[4])


def with_sup(base, cell):
    c = dict(C.DEFAULTS); c.update(base); c.update(k0=3e4, **SUP); c['fb'] = fb_of(cell)
    c['solver'] = 'ord'; c['iters'] = 1
    return c


def v3cfg(cell):
    """C1v3: the v1 configuration tuned on the 1 mm counterpart (no per-cell re-tuning for breathing amplitude),
    v2 supervisor and start prior, the cell's B1 gains as fallback."""
    return with_sup(V.CFG1[COUNTERPART[cell]], cell)


def v2cfg(cell):
    return V.v2cfg(cell, SUP)


def shared_cfg(cell):
    s = R.get('shared_tuning', {}).get('best')
    if s is None: raise RuntimeError('run shared_tune first')
    return with_sup(s, cell)


# declared shared-configuration search (one estimator/OCP configuration for all six cells)
NSH_TR = 40


def shared_candidates_A():
    rng = np.random.default_rng(12345)
    anchors = [dict(V.CFG1['T3']), dict(V.CFG1['T1']), dict(V.CFG1['T3fast'])]
    for a in anchors:
        a.pop('solver', None); a.pop('iters', None)
    return anchors + C.stage_a(rng)


def shared_score(Jc, Jb1, feas_slow):
    """mean over cells of log(mean J / mean J_B1) (geometric-mean J ratio to the cell's best fixed controller);
    infeasible (peak > 7 N in > 5 % of train scenes in a non-fast cell) gets +1e6."""
    return float(np.mean([math.log(Jc[c] / Jb1[c]) for c in Jc])) + (0. if feas_slow else 1e6)


def shared_eval(cands, train, log=print, tag=''):
    out = []
    specs_by_cell = {}
    for c in ORDER:
        specs_by_cell[c] = [('C1v3', with_sup(k, c)) for k in cands]
    Jc = {c: [] for c in ORDER}; feas = {c: [] for c in ORDER}
    for c in ORDER:
        r = run(train[c], specs_by_cell[c], TR0)
        for j in range(len(cands)):
            rows = [x[j] for x in r]
            Jc[c].append(float(np.mean([q['J'] for q in rows])))
            feas[c].append(bool(np.mean([q['peak'] <= FD + 3.0 for q in rows]) >= 0.95))
        log(f"[shared {tag}] {c} done")
    return Jc, feas


def shared_tune(log=print):
    t0 = time.time()
    train = {c: S.cohort(CELLS[c][0], NSH_TR, TR0, CELLS[c][1]) for c in ORDER}
    Jb1 = {}
    for c in ORDER:
        r = run(train[c], [('fixed', tuple(V.gate_ref(c)['B1']['gains']))], TR0)
        Jb1[c] = float(np.mean([x[0]['J'] for x in r]))
    A = shared_candidates_A()
    JA, fA = shared_eval(A, train, log, 'A')
    slow = [c for c in ORDER if 'fast' not in c]
    sc = lambda J, F, j: shared_score({c: J[c][j] for c in ORDER}, Jb1, all(F[c][j] for c in slow))
    sA = [sc(JA, fA, j) for j in range(len(A))]
    jA = int(np.argmin(sA))
    B = C.stage_b(A[jA], np.random.default_rng(777))
    JB, fB = shared_eval(B, train, log, 'B')
    sB = [sc(JB, fB, j) for j in range(len(B))]
    cand = A + B; sall = sA + sB
    J = {c: JA[c] + JB[c] for c in ORDER}; Fe = {c: fA[c] + fB[c] for c in ORDER}
    j = int(np.argmin(sall))
    order = np.argsort(sall)
    log(f"[shared] B1 train J {Jb1}")
    for i in order[:8]:
        log(f"   score {sall[i]:.3f} " + ' '.join(f"{c} {J[c][i]:.3f}" for c in ORDER) + f" cfg {cand[i]}")
    log(f"[shared] selected #{j} score {sall[j]:.3f} ({len(cand)} configs, {time.time() - t0:.0f}s)")
    return dict(best=cand[j], score=sall[j], n_configs=len(cand), n_train_per_cell=NSH_TR, B1_train_J=Jb1,
                train_J={c: J[c][j] for c in ORDER}, all=[dict(cfg=cand[i], score=sall[i], J={c: J[c][i] for c in ORDER},
                                                             feas={c: Fe[c][i] for c in ORDER}) for i in range(len(cand))],
                anchor_scores=dict(T3cfg=sall[0], T1cfg=sall[1], T3fastcfg=sall[2]), seconds=time.time() - t0)


def summ(rows):
    out = V.summ(rows)
    out['frac_scenes_with_loss'] = float(np.mean([r['loss'] > 0 for r in rows]))
    return out


def evaluate(cell, seed0, n, cfgs, log=print, with_base=True):
    task, dose = CELLS[cell]; scenes = S.cohort(task, n, seed0, dose)
    names = list(cfgs); specs = [('C1v3', cfgs[k]) for k in names]
    src = None
    if with_base:
        bs, bn, src = V.baselines(cell); specs += bs; names += bn
    res = run(scenes, specs, seed0)
    per = {nm: [r[j] for r in res] for j, nm in enumerate(names)}
    out = {}
    for nm in names:
        rows = per[nm]
        out[nm] = dict(summary=summ(rows), J_per_scene=[float(r['J']) for r in rows], rmse_per_scene=[float(r['rmse']) for r in rows],
                       peak_per_scene=[float(r['peak']) for r in rows], loss_per_scene=[float(r['loss']) for r in rows])
    for nm in cfgs:
        for rf in names:
            if rf == nm: continue
            out[nm]['dJ_vs_' + rf] = U.paired_ci([r['J'] for r in per[nm]], [r['J'] for r in per[rf]])
            out[nm]['dover7_vs_' + rf] = U.paired_ci([float(r['peak'] > 7.0) for r in per[nm]], [float(r['peak'] > 7.0) for r in per[rf]])
            out[nm]['drmse_vs_' + rf] = U.paired_ci([r['rmse'] for r in per[nm]], [r['rmse'] for r in per[rf]])
            out[nm]['dlossscene_vs_' + rf] = U.paired_ci([float(r['loss'] > 0) for r in per[nm]], [float(r['loss'] > 0) for r in per[rf]])
            out[nm]['wins_vs_' + rf] = int(sum(a['J'] < b['J'] for a, b in zip(per[nm], per[rf])))
    for nm in names:
        s = out[nm]['summary']
        ex = ''
        if 'dJ_vs_B1' in out[nm]:
            a = out[nm]['dJ_vs_B1']; ex = f" dJ vs B1 {a[0]:+.3f} ({a[1]:+.3f},{a[2]:+.3f})"
        log(f"[{cell} seed {seed0} n {n}] {nm:12s} J {s['J']:7.3f} (med {s['J_median']:.3f}) rmse {s['rmse']:.3f} >7N {s['frac_peak_over']:.3f} "
            f"loss {s['loss']:.4f} lossscenes {s['frac_scenes_with_loss']:.3f} J>20 {s['n_J_over_20']}" + ex)
    return dict(cell=cell, task=task, dose=dose, seed0=seed0, n=n, configs=cfgs, etr_source=src, methods=out)


def known(cfg, b1, log=print, rng_seed=7, nper=3, seed0=500):
    """v2 known-behaviour check (N1, 1 mm breathing, 15 kPa and 1 MPa)."""
    rng = np.random.default_rng(rng_seed); scenes = []
    for E in [1.5e4, 1.0e6]:
        for k in range(nper): scenes.append(S.make('N1', rng, dict(E=E, A=0.001)))
    res = run(scenes, [('C1v3', cfg), ('fixed', tuple(b1))], seed0, keep=True)
    out = {'soft_15kPa': [], 'stiff_1MPa': []}
    for (p, T), r in zip(scenes, res):
        m, mb = r[0], r[1]; h = m['hist']; lg = m['log']
        kt = lg[::4, 8][:len(h)]; ok = (h[:, 0] > 1.0) & (kt > 0)
        d = dict(E=float(p[0]), k_true_median=float(np.median(kt[ok])), k_hat_median=float(np.median(h[ok, 1])),
                 k_hat_over_true_median=float(np.median(h[ok, 1] / kt[ok])), D_median=float(np.median(h[ok, 4])),
                 KPF_median=float(np.median(h[ok, 5])), J=m['J'], rmse=m['rmse'], peak=m['peak'], loss=m['loss'],
                 osc=m['osc'], n_loss_entries=m['n_loss_entries'], J_B1=mb['J'], rmse_B1=mb['rmse'], peak_B1=mb['peak'], osc_B1=mb['osc'])
        out['soft_15kPa' if p[0] < 1e5 else 'stiff_1MPa'].append(d)
    for k, Ls in out.items():
        log(f"[known] {k}: k_true {np.mean([d['k_true_median'] for d in Ls]):.0f} k_hat {np.mean([d['k_hat_median'] for d in Ls]):.0f} "
            f"ratio {np.mean([d['k_hat_over_true_median'] for d in Ls]):.2f} D {np.mean([d['D_median'] for d in Ls]):.0f} KPF {np.mean([d['KPF_median'] for d in Ls]):.2e} "
            f"J {np.mean([d['J'] for d in Ls]):.4f} vs B1 {np.mean([d['J_B1'] for d in Ls]):.4f} rmse {np.mean([d['rmse'] for d in Ls]):.4f} vs "
            f"{np.mean([d['rmse_B1'] for d in Ls]):.4f} peak max {np.max([d['peak'] for d in Ls]):.2f} loss entries {np.sum([d['n_loss_entries'] for d in Ls])}")
    return out


def _save():
    json.dump(R, open('candidate1_v3_results.json', 'w'), indent=1, default=float)


try:
    R = json.load(open('candidate1_v3_results.json'))
except Exception:
    R = {}


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'dev'
    arg = sys.argv[2] if len(sys.argv) > 2 else ''
    L = []
    def log(s):
        print(s, flush=True); L.append(s)
    if what == 'shared_tune':
        R['shared_tuning'] = shared_tune(log); _save()
    elif what == 'dev':
        # non-degradation check on train / dev: C1v2 vs C1v3 (and the shared configuration when tuned)
        cells = arg.split(',') if arg else ORDER
        for c in cells:
            cf = {'C1v2': v2cfg(c), 'C1v3': v3cfg(c)}
            if 'shared_tuning' in R: cf['C1v3-shared'] = shared_cfg(c)
            ntr, ndv = (40, 60) if 'fast' in c else (40, 40)
            for split, s0, n in [('train', TR0, ntr), ('dev', DEV0, ndv)]:
                R.setdefault('nondegradation', {}).setdefault(c, {})[split] = evaluate(c, s0, n, cf, log); _save()
    elif what == 'final':
        cells = arg.split(',') if arg else ORDER
        for c in cells:
            cf = {'C1v3': v3cfg(c), 'C1v2': v2cfg(c)}
            if 'shared_tuning' in R: cf['C1v3-shared'] = shared_cfg(c)
            R.setdefault('final_seed5000', {})[c] = evaluate(c, FRESH5, NFINAL, cf, log); _save()
    elif what == 'known':
        out = {}
        for nm, cfgf in [('C1v3_T3', lambda: v3cfg('T3')), ('C1v3_T1', lambda: v3cfg('T1')), ('C1v3_T3fast', lambda: v3cfg('T3fast'))] + \
                ([('C1v3-shared_T3', lambda: shared_cfg('T3'))] if 'shared_tuning' in R else []):
            cell = nm.split('_', 1)[1]
            log(f"[known] configuration {nm}")
            out[nm] = known(cfgf(), V.gate_ref(cell)['B1']['gains'], log)
        R['known_behaviour'] = out; _save()
    R.setdefault('console', []).extend(L); _save()
    close_pool()
