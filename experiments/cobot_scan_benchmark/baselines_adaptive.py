"""Causal adaptive comparison baselines (README 'To be added before candidates').

All policies follow the plant.episode interface: reset(platform_dict), __call__(obs) -> gains
(M, D, K, ki, KPF[, u_ff]) once per control period. They read only causal signals: obs.t,
obs.force (filtered measured normal force), obs.y (encoder position), the tick index and their
own state. No policy here reads obs.k_true or log columns 3, 4, 6, 7, 8 (those are used only in
the offline verification `check_estimator`, to score the estimate).

Plant convention (plant.py): penetration coordinate y (= -z), ef = FD - F_filt,
    M yc'' + D yc' + K (yc - y_ref) = ef + ki * int(ef)        y_cmd = yc + KPF * ef + u_ff

1. EtR, estimate-then-rule
   Estimator EtR-A (lock-in, primary). A sinusoidal command dither of amplitude a is added via
   u_ff at f_d = 1/(N_c Tc) with N_c = round(1/(10 Hz * Tc)) samples per cycle (10.42 Hz on PC125,
   an integer number of samples per cycle, so a boxcar of whole cycles rejects DC and harmonics).
   Both the measured force F and the encoder position y are first-differenced (removes ramps
   from breathing / scan drift; scales F and y by the same factor at f_d) and demodulated against
   the same reference exp(-j 2 pi i / N_c) over a sliding window of n_cyc cycles:
       Z_F = sum_w dF_i e^{-j theta_i},  Z_y = sum_w dy_i e^{-j theta_i},
       k_hat = |Z_F| / (|Z_y| |H_ft(f_d)|)
   H_ft is the known first-order F/T low-pass (platform fc, sampled at 1 kHz). k_hat is the
   dynamic contact stiffness at f_d (SLS instantaneous branch, which is what plant.k_eff
   reports). The window is frozen while F < 0.5 N in it (contact loss). The dither amplitude is
   force-targeted, a = clip(F_dith / k_hat, a_min, a_max) (10-30 um on soft tissue, smaller on
   bone so the dither stays a few hundredths of a newton).
   Estimator EtR-R (RLS, optional): scalar recursive least squares with forgetting lambda on the
   increments, dF = k dy:  S_xy <- lam S_xy + dy dF,  S_yy <- lam S_yy + dy^2,  k_hat = S_xy/S_yy.
   The surface moves (breathing, patient events), so dF = k (dh + dy) and the estimate is biased
   by the unmeasured dh; kept for comparison only.
   Rule 'truthk' (same as oracle B3b, baselines_ref.TruthK): D = max(Dmin, k_hat/wI),
   KPF = cP/k_hat, ki fixed, M = 0.3 kg, K = 0.
   Rule 'haninger' (Haninger, Hegeler, Peternel RAS 2023, well-damped contact constraint,
   literature/extracted_text/haninger.txt:452-480, Eq. (20) at line 476):
       D_i >= 2 zeta sqrt(M_i K_e,i(x_t))
   applied with equality (the least damping that satisfies the constraint): D = max(Dmin,
   2 zeta sqrt(M k_hat)); in Haninger the stiffness comes from the GP force-model gradient
   (Eq. (19), line 469); here it comes from the lock-in estimate. KPF = cP/k_hat, ki fixed.

2. Duan RAS 2018 adaptive variable impedance (literature/extracted_text/field_supplements/
   F02.txt). Contact impedance with K = 0 (Eq. (9), line 271 and lines 266-286):
       f_e - f_d = m e'' + b e',  e = x_c - x_e
   Adaptive damping (Eqs. (13)-(14), lines 231-247):
       f_e - f_d = m e'' + (b + Delta_b(t)) e',  Delta_b(t) = b Phi(t) / e'(t),
       Phi(t) = Phi(t - lambda) + sigma (f_d(t - lambda) - f_e(t - lambda)) / b,
   lambda = controller sampling period, sigma = update rate. Stability bound (Eq. (30), line
   346): 0 < sigma < lambda b / (m + lambda b). Paper values: m = 1, b = 100, sigma = 0.5 with
   5 ms cycle (lines 465-468, 556, 624-625).
   Mapping. Delta_b e' = b Phi, and unrolling the recursion b Phi(t) = sigma sum_{j<t} (f_d - f_e)
   (t_j) = (sigma/lambda) int (f_d - f_e) dt (one-sample delay). With e = x_c - x_e and the
   penetration coordinate y = -x_c (+const), f_d - f_e = ef:
       m y'' + b y' = ef + (sigma/lambda) int ef dt
   so Duan is EXACTLY our admittance with M = m, D = b, K = 0, KPF = 0, ki = sigma/lambda
   (lambda = Tc = 8 ms). The 'variable damping' Delta_b = b Phi / e' is singular at e' = 0 and is
   only the paper's way of writing an integral force term; implementing it literally as a time-
   varying D would be numerically ill-posed and is algebraically identical to the integral.
   x_e is replaced by its initial value (no surface model is available, F02 Eq. (7), line
   236), so the e_e'' / e_e' feedforward terms of Eq. (15) vanish. ki <= 30 (plant bound) caps
   sigma at 0.24 at Tc = 8 ms (the paper's sigma = 0.5 at 5 ms is ki = 100, outside the bound).
   Duan is therefore a member of the B1 fixed family (K = 0, KPF = 0), tuned with its own budget.

3. Song, Fang, Qian RA-L 2025 NAIPC (literature/extracted_text/field_papers/E04.txt).
   Impedance with compensator (Eqs. (14)-(15), lines 228-262), f~ = f_d - f_e:
       m_d x_c'' + b_d x_c' + f~ = u,  u = -d + f~ - alpha f~ - gamma f~',
       d' = beta (-u + f~ - d)
   The paper's own reduction (Eqs. (21)-(25), lines 343-365): u = (1 - alpha - gamma beta -
   gamma s - alpha beta / s) f~. Hence, in the penetration coordinate (y = -x_c):
       m_d y'' + b_d y' = (alpha + gamma beta) ef + gamma ef' + alpha beta int ef + const
   With c = alpha + gamma beta (const = d(0) - beta gamma f~(t0) = 0 here: every scene starts in
   contact at FD, so d(0) = 0 is the consistent initial state; d(0) is a contact-transition
   device in the paper):
       M = m_d/c, D = b_d/c, K = 0, ki = alpha beta / c,  plus  gamma ef'/(m_d s^2 + b_d s)
       = (gamma/b_d) ef / (tau s + 1),  tau = m_d/b_d = M/D,
   i.e. a direct force path KPF = gamma/b_d through a first-order lag tau. Implemented exactly:
   plant KPF = gamma/b_d on the current sample, and u_ff = KPF (lp(ef) - ef) with lp the ZOH
   lag of tau (on the latest available sample; one sample older than the plant's ef, the only
   approximation). Constraints implied by the paper's alpha, beta, gamma > 0 (Eq. (30), line
   388-394): with g = gamma/c = KPF D and q = alpha beta / c = ki,  g q = alpha beta gamma / c^2
   <= 1/4 (AM-GM), i.e. KPF * D * ki <= 1/4; and Eq. (30) m_d alpha beta < b_d c, i.e. M ki < D.
   Scale invariance (m_d, b_d, alpha, gamma) -> s (...) leaves the controller unchanged; we tune
   in the mapped space (M, D, ki, KPF) and report one native representative with c = 1.
"""
import math, json, sys, time, itertools
import numpy as np
from multiprocessing import Pool
import plant as P, scenarios as S, tune as U

# ------------------------------------------------------------------ estimators
class LockInK:
    """incremental lock-in estimate of contact stiffness from causal F and encoder y.
    update(i, F, y, a): tick index i, filtered force, encoder position, dither amplitude that
    was commanded for that tick (used only for the excitation check rho)."""
    def __init__(self, Tc, fc, ncyc=2, f_dither=10., k0=2e4, kmin=300., kmax=3e5, proj=False, rho=0., ndiff=1):
        self.Nc = max(4, int(round(1.0 / (Tc * f_dither)))); self.f = 1.0 / (self.Nc * Tc)
        self.W = int(ncyc * self.Nc)
        self.ref = np.exp(-2j * np.pi * np.arange(self.Nc) / self.Nc)
        afc = 1.0 - math.exp(-2 * math.pi * fc * 1e-3); z = np.exp(-2j * np.pi * self.f * 1e-3)
        self.Hft = afc / (1.0 - (1.0 - afc) * z)        # complex response of the known F/T low-pass
        self.ndiff = ndiff                                # order of the differencing (removes trends)
        self.dgain = (2 * math.sin(math.pi / self.Nc)) ** ndiff   # |1 - e^{-j w Tc}|^ndiff
        self.proj = proj; self.rho = rho
        self.bF = np.zeros(self.W, complex); self.bY = np.zeros(self.W, complex); self.bok = np.ones(self.W, bool)
        self.bA = np.zeros(self.W)
        self.SF = 0j; self.SY = 0j; self.SA = 0.; self.nbad = 0; self.cnt = 0; self.prev = None
        self.k = k0; self.kmin = kmin; self.kmax = kmax
    def update(self, i, F, y, a=0.):
        if self.prev is None:
            self.prev = [(F, y)] * (self.ndiff + 1); return self.k
        self.prev = self.prev[1:] + [(F, y)]
        hF = np.array([q[0] for q in self.prev]); hY = np.array([q[1] for q in self.prev])
        for _ in range(self.ndiff): hF = np.diff(hF); hY = np.diff(hY)
        dF = hF[0]; dy = hY[0]
        r = self.ref[i % self.Nc]; j = self.cnt % self.W
        zF = dF * r; zY = dy * r; ok = F > 0.5
        self.SF += zF - self.bF[j]; self.SY += zY - self.bY[j]; self.SA += a - self.bA[j]
        self.nbad += int(not ok) - int(not self.bok[j])
        self.bF[j] = zF; self.bY[j] = zY; self.bok[j] = ok; self.bA[j] = a; self.cnt += 1
        if self.cnt % (4 * self.W) == 0:          # re-sum to bound round-off
            self.SF = self.bF.sum(); self.SY = self.bY.sum(); self.SA = self.bA.sum()
        # excitation check: encoder response at f_d at least rho x the commanded dither
        thr = max(1e-9 * self.W, self.rho * 0.5 * self.dgain * self.SA)
        if self.cnt >= self.W and self.nbad == 0 and abs(self.SY) > thr:
            ZF = self.SF / self.Hft
            kk = (ZF * self.SY.conjugate()).real / abs(self.SY) ** 2 if self.proj else abs(ZF) / abs(self.SY)
            self.k = min(max(kk, self.kmin), self.kmax)
        return self.k

class RLSK:
    """scalar RLS with forgetting on increments dF = k dy (biased by surface motion)."""
    def __init__(self, Tc, fc, lam=0.98, k0=2e4, kmin=300., kmax=3e5, f_dither=10.):
        self.lam = lam; self.Sxy = 0.; self.Syy = 0.; self.prev = None; self.k = k0
        self.kmin = kmin; self.kmax = kmax; self.Nc = max(4, int(round(1.0 / (Tc * f_dither))))
    def update(self, i, F, y, a=0.):
        if self.prev is None:
            self.prev = (F, y); return self.k
        dF = F - self.prev[0]; dy = y - self.prev[1]; self.prev = (F, y)
        if F > 0.5:
            self.Sxy = self.lam * self.Sxy + dy * dF; self.Syy = self.lam * self.Syy + dy * dy
            if self.Syy > 1e-14 and self.Sxy > 0:
                self.k = min(max(self.Sxy / self.Syy, self.kmin), self.kmax)
        return self.k

# ------------------------------------------------------------------ 1. estimate-then-rule
class EtR:
    """estimate-then-rule. rule 'truthk': D = max(Dmin, k/wI); 'haninger': D = max(Dmin,
    2 zeta sqrt(M k)); both KPF = cP/k, fixed ki, K = 0; 'fixed': given gains (estimator check).
    est 'lockin' (EtR-A) or 'rls' (EtR-R). trel > 0: peak hold of the estimate with release time
    constant trel (fast attack, slow release; conservative at bone entry)."""
    oracle = False
    def __init__(self, wI=80., cP=2., ki=30., Dmin=20., M=0.3, rule='truthk', zeta=1., est='lockin',
                 ncyc=2, f_dither=10., Fdith=0.1, amin=2e-6, amax=30e-6, rho=0., trel=0., ndiff=1, lam=0.98,
                 k0=2e4, dither=True, proj=False, fixed=None):
        self.wI = wI; self.cP = cP; self.ki = ki; self.Dmin = Dmin; self.M = M; self.rule = rule
        self.zeta = zeta; self.est = est; self.ncyc = ncyc; self.fd = f_dither; self.Fdith = Fdith
        self.amin = amin; self.amax = amax; self.rho = rho; self.trel = trel; self.ndiff = ndiff; self.lam = lam
        self.k0 = k0; self.dither = dither; self.proj = proj; self.fixed = fixed
    def reset(self, pl):
        self.Tc = pl['Tc']
        if self.est == 'lockin':
            self.E = LockInK(pl['Tc'], pl['fc'], self.ncyc, self.fd, k0=self.k0, proj=self.proj, rho=self.rho,
                            ndiff=self.ndiff)
        else: self.E = RLSK(pl['Tc'], pl['fc'], self.lam, k0=self.k0, f_dither=self.fd)
        self.k = self.k0; self.khist = []; self.a_prev = 0.
    def gains(self, k):
        if self.rule == 'fixed': return tuple(self.fixed)
        if self.rule == 'truthk': D = max(self.Dmin, k / self.wI)
        else: D = max(self.Dmin, 2 * self.zeta * math.sqrt(self.M * k))
        return (self.M, D, 0., self.ki, self.cP / k)
    def __call__(self, o):
        n = o.n
        if n > 0:
            ke = self.E.update(n - 1, o.log[n - 1, 5], o.log[n - 1, 2], self.a_prev)  # causal: F_filt, encoder y
            if self.trel > 0: ke = max(ke, self.k * math.exp(-self.Tc / self.trel))
            self.k = ke
        self.khist.append(self.k)
        u = 0.; a = 0.
        if self.dither:
            a = min(max(self.Fdith / self.k, self.amin), self.amax)
            u = a * math.sin(2 * math.pi * (n % self.E.Nc) / self.E.Nc)
        self.a_prev = a
        return self.gains(self.k) + (u,)

# ------------------------------------------------------------------ 2. Duan RAS2018
def duan_gains(m, b, sigma, Tc=0.008):
    """F02 Eqs. (13)-(14) mapped onto the plant admittance: (M, D, K, ki, KPF)."""
    return (m, b, 0., sigma / Tc, 0.)

def duan_sigma_max(m, b, Tc=0.008):
    """F02 Eq. (30): 0 < sigma < lambda b / (m + lambda b)."""
    return Tc * b / (m + Tc * b)

class Duan:
    """policy form of Duan (fixed gains; the adaptation lives in the plant's force integral)."""
    oracle = False
    def __init__(self, m, b, sigma): self.m = m; self.b = b; self.sigma = sigma
    def reset(self, pl): self.g = duan_gains(self.m, self.b, self.sigma, pl['Tc'])
    def __call__(self, o): return self.g

# ------------------------------------------------------------------ 3. Song RA-L2025 NAIPC
def song_native(M, D, ki, KPF):
    """one native (m_d, b_d, alpha, beta, gamma) representative (c = alpha + beta gamma = 1) of
    the mapped gains; None if outside alpha, beta, gamma > 0 (needs KPF D ki <= 1/4)."""
    g = KPF * D
    if g <= 0 or ki <= 0:
        return None if (g < 0 or ki < 0) else dict(md=M, bd=D, alpha=1.0 if ki > 0 else 0., beta=ki, gamma=g)
    disc = 1 - 4 * g * ki
    if disc < 0: return None
    beta = (1 - math.sqrt(disc)) / (2 * g); alpha = 1 - beta * g
    return dict(md=M, bd=D, alpha=alpha, beta=beta, gamma=g)

def song_feasible(M, D, ki, KPF):
    return KPF * D * ki <= 0.25 + 1e-12 and M * ki < D

class Song:
    """NAIPC (E04 Eqs. (14)-(15)) in mapped form: PI admittance (K = 0) + force path
    KPF/(tau s + 1), tau = M/D. Parameters are the mapped gains (see module docstring)."""
    oracle = False
    def __init__(self, M, D, ki, KPF): self.M = M; self.D = D; self.ki = ki; self.KPF = KPF
    def reset(self, pl):
        tau = self.M / self.D; self.a = math.exp(-pl['Tc'] / tau); self.lp = 0.
    def __call__(self, o):
        corr = 0.
        if o.n > 0:
            ef = P.FD - o.log[o.n - 1, 5]
            self.lp = self.a * self.lp + (1 - self.a) * ef
            corr = self.KPF * (self.lp - ef)
        return (self.M, self.D, 0., self.ki, self.KPF, corr)

# ------------------------------------------------------------------ evaluation (<= 2 processes)
_POOL = None
def pool():
    global _POOL
    if _POOL is None: _POOL = Pool(2)
    return _POOL

def evaluate(scenes, ctrls, plat='PC125', seed0=0):
    jobs = [(p, T, c, plat, seed0 + i, None) for i, (p, T) in enumerate(scenes) for c in ctrls]
    res = pool().map(U._eval, jobs, chunksize=max(1, len(jobs) // 32))
    nc = len(ctrls); return [res[i * nc:(i + 1) * nc] for i in range(len(scenes))]

def select(res):
    Jm = U.J(res).mean(0); feas = U.feasible(res)
    j = int(np.argmin(np.where(feas, Jm, Jm + 1e6)))
    return j, float(Jm[j]), bool(feas[j]), np.argsort(np.where(feas, Jm, Jm + 1e6))

def refine(make, z0, lo, hi, train, seed0, maxfev):
    """Nelder-Mead on a box-clipped parameter vector z (objective: mean train J, +1e3 if the
    peak constraint fails), same selection rule as tune.best_fixed."""
    from scipy.optimize import minimize
    def obj(z):
        zc = np.clip(z, lo, hi); r = evaluate(train, [make(zc)], 'PC125', seed0)
        return float(np.mean([row[0]['J'] for row in r])) + (0 if U.feasible(r)[0] else 1e3)
    r = minimize(obj, z0, method='Nelder-Mead', options=dict(maxfev=maxfev, xatol=1e-3, fatol=1e-4))
    zc = np.clip(r.x, lo, hi); return zc, obj(zc)

# ------------------------------------------------------------------ tuning (train cohort)
CELLS = [('T3', dict(ramp=0.5, A=0.001)), ('T1', dict(ratio=30., A=0.001))]
NTR, NTE = 20, 40
RULE = list(itertools.product([20., 40., 80., 160.], [0.5, 1., 2., 4.], [3., 10., 30.], [20., 100.]))  # = gate.RULE
HAN = list(itertools.product([0.5, 1., 2., 4.], [0.3, 1.], [0.5, 1., 2., 4.], [3., 10., 30.]))     # zeta, M, cP, ki
EST0 = dict(est='lockin', ncyc=2, f_dither=20., ndiff=2, rho=0.3, trel=0.1, Fdith=0.1)
EST_VAR = [dict(EST0, ncyc=nc, trel=tr, ndiff=nd) for nc in [1, 2] for tr in [0., 0.1] for nd in [2, 3]]
TR0 = S.SPLIT['train']; TE0 = S.SPLIT['test']

def _grid_select(train, ctrls):
    r = evaluate(train, ctrls, 'PC125', TR0); j, Jm, f, order = select(r)
    return j, Jm, f, order

def tune_etr(train, rule):
    """stage 1: rule grid with the default lock-in settings; stage 2: estimator variants on the 3
    best rules (no local refine, same budget style as the B3b oracle in gate.py)."""
    if rule == 'truthk':
        mk = lambda q, e: EtR(q[0], q[1], q[2], q[3], rule='truthk', **e)
        Q = RULE
    elif rule == 'haninger':
        mk = lambda q, e: EtR(cP=q[2], ki=q[3], M=q[1], zeta=q[0], Dmin=20., rule='haninger', **e)
        Q = HAN
    else:   # RLS estimator with the truthk rule, reduced rule grid
        mk = lambda q, e: EtR(q[0], q[1], q[2], q[3], rule='truthk', **e)
        Q = list(itertools.product([40., 80., 160.], [1., 2., 4.], [10., 30.], [20., 100.]))
    E1 = EST0 if rule != 'rls' else dict(EST0, est='rls', lam=0.98)
    j, Jm, f, order = _grid_select(train, [mk(q, E1) for q in Q])
    tops = [Q[i] for i in order[:3]]
    EV = EST_VAR if rule != 'rls' else [dict(E1, lam=l) for l in [0.9, 0.95, 0.99]]
    cand = [(q, e) for q in tops for e in EV]
    j2, Jm2, f2, _ = _grid_select(train, [mk(q, e) for q, e in cand])
    if (not f2 and f) or (f == f2 and Jm <= Jm2): best = (Q[j], E1, Jm, f)
    else: best = (cand[j2][0], cand[j2][1], Jm2, f2)
    return dict(rule_params=best[0], est_params=best[1], train_J=best[2], train_feasible=best[3],
                n_configs=len(Q) + len(cand)), mk(best[0], best[1])

def _logz(g, off): return np.log(np.asarray(g) + off)

def tune_duan(train, Tc=0.008, maxfev=80):
    """grid M x D x ki (ki = sigma/Tc), feasible under F02 Eq. (30), then Nelder-Mead in log
    space from the 3 best points (fixed gains: exact, fast numba path)."""
    grid = [(m, b, ki) for m in U.GRID_M for b in U.GRID_D for ki in U.GRID_KI if ki * Tc < duan_sigma_max(m, b, Tc)]
    ctr = [duan_gains(m, b, ki * Tc, Tc) for (m, b, ki) in grid]
    j, Jm, f, order = _grid_select(train, ctr)
    off = np.array([0., 0., 0.1]); lo = np.log(np.array([0.25, 20., 0.]) + off); hi = np.log(np.array([3., 2500., 30.]) + off)
    def make(z):
        m, b, ki = np.exp(z) - off; ki = min(ki, 0.999 * duan_sigma_max(m, b, Tc) / Tc)
        return duan_gains(m, b, ki * Tc, Tc)
    best = (Jm + (0 if f else 1e3), ctr[j])
    for i in order[:3]:
        z, v = refine(make, _logz(grid[i], off), lo, hi, train, TR0, maxfev)
        if v < best[0]: best = (v, make(z))
    g = best[1]
    return dict(m=g[0], b=g[1], sigma=g[3] * Tc, ki=g[3], sigma_max_eq30=duan_sigma_max(g[0], g[1], Tc),
                train_J=best[0], grid_best=grid[j], n_grid=len(grid)), tuple(float(x) for x in g)

def tune_song(train, maxfev=60):
    """grid in the mapped space (M, D, ki, KPF) restricted to alpha, beta, gamma > 0
    (KPF D ki <= 1/4) and E04 Eq. (30) (M ki < D); then Nelder-Mead from the 2 best points."""
    GD = list(np.round(np.geomspace(20, 2500, 6), 1))
    grid = [q for q in itertools.product(U.GRID_M, GD, U.GRID_KI, [0., 5e-5, 1e-4, 2e-4, 4e-4]) if song_feasible(*q)]
    j, Jm, f, order = _grid_select(train, [Song(*q) for q in grid])
    off = np.array([0., 0., 0.1, 1e-6]); lo = np.log(np.array([0.25, 20., 0., 0.]) + off)
    hi = np.log(np.array([3., 2500., 30., 1e-3]) + off)
    def proj(z):
        M, D, ki, K = np.exp(np.clip(z, lo, hi)) - off
        ki = min(ki, 0.25 / max(K * D, 1e-12), 0.999 * D / M); return (M, D, max(ki, 0.), max(K, 0.))
    best = (Jm + (0 if f else 1e3), grid[j])
    for i in order[:2]:
        z, v = refine(lambda z: Song(*proj(z)), _logz(grid[i], off), lo, hi, train, TR0, maxfev)
        if v < best[0]: best = (v, proj(z))
    q = tuple(float(x) for x in best[1])
    return dict(mapped=dict(M=q[0], D=q[1], ki=q[2], KPF=q[3], tau=q[0] / q[1]), native_c1=song_native(*q),
                train_J=best[0], grid_best=grid[j], n_grid=len(grid)), Song(*q)

# ------------------------------------------------------------------ known-behaviour checks
SAFE = (1.0, 600., 0., 10., 2e-5)   # conservative fixed gains, stable from 10 kPa to 1 MPa

def check_estimator(est_params, rule_policy=None):
    """offline: lock-in estimate vs obs.k_true (truth used ONLY for scoring) on homogeneous soft /
    stiff tissue (N1, breathing 1 mm) under a conservative fixed controller ('fixed' mode, so the
    estimate does not act on the plant) and in closed loop with the tuned EtR-A policy; plus the
    T1 ladder step response (delay of the estimate at ridge entry)."""
    out = {}
    for E in [1e4, 3e4, 1e5, 3e5, 1e6]:
        sc = S.cohort('N1', 6, 5000, dict(E=E, A=0.001)); row = {}
        for mode in ['fixed', 'closed']:
            rs = []
            for i, (p, T) in enumerate(sc):
                if mode == 'fixed': pol = EtR(rule='fixed', fixed=SAFE, **est_params)
                else:
                    if rule_policy is None: continue
                    import copy; pol = copy.deepcopy(rule_policy)
                log = P.episode(p, pol, T=T, seed=5000 + i)
                if log is None: continue
                kh = np.array(pol.khist[:len(log)]); kt = log[:, 8]; m = (log[:, 0] > 1.) & (kt > 0)
                r = kh[m] / kt[m]; rs.append((np.median(r), np.percentile(r, 10), np.percentile(r, 90), kt[m].mean(), np.median(kh[m]),
                                              S.metrics(log)['rmse']))
            if rs:
                a = np.array(rs).mean(0)
                row[mode] = dict(ratio_median=a[0], ratio_p10=a[1], ratio_p90=a[2], k_true_mean=a[3], k_hat_median=a[4], rmse=a[5])
        out['N1_E%.0e' % E] = row
    # T1 ladder: delay of the estimate at ridge entry (closed loop, tuned EtR-A)
    if rule_policy is not None:
        import copy
        sc = S.cohort('T1', 10, 5000, dict(ratio=30., A=0.001)); dl = []; soft = []; bone = []
        for i, (p, T) in enumerate(sc):
            pol = copy.deepcopy(rule_policy); log = P.episode(p, pol, T=T, seed=5000 + i)
            if log is None: continue
            kh = np.array(pol.khist[:len(log)]); kt = log[:, 8]; t = log[:, 0]
            ks = np.median(kt[(t > 1) & (t < 1.3)]); kb = kt.max(); thr = math.sqrt(ks * kb)
            it = np.argmax(kt > thr); ih = it + np.argmax(kh[it:] > thr) if (kh[it:] > thr).any() else None
            if ih is not None: dl.append(t[ih] - t[it])
            ms = (t > 1) & (kt > 0) & (kt < 1.5 * ks); mb = kt > 0.8 * kb
            soft.append(np.median(kh[ms] / kt[ms])); bone.append(np.median(kh[mb] / kt[mb]))
        out['T1_ladder'] = dict(entry_delay_s_median=float(np.median(dl)), entry_delay_s_p90=float(np.percentile(dl, 90)),
                                n_detected=len(dl), ratio_soft_median=float(np.median(soft)), ratio_bone_median=float(np.median(bone)))
    return out

def check_slope(duan_g, song_pol, Tc=0.008):
    """slope-like drift: homogeneous 30 kPa, no breathing, one slow push (12 mm over a 6 s
    smoothstep, near-constant 3 mm/s in its middle). Constant impedance (Duan with sigma = 0, i.e.
    F02 Eq. (11)) has the bias f_e - f_d = b x_e' (F02 Sec. 5.1.2); Duan (sigma > 0) and Song
    should remove it. Reports mean F - FD (true force) over the middle 2 s of the ramp."""
    rng = np.random.default_rng(0); p = np.zeros(P.NP)
    p[0] = 3e4; p[1] = 0.2; p[2] = 0.; p[3] = 0.25; p[5] = 0.015; p[13] = 1.; p[14] = 1.
    p[15] = 1; p[16:20] = [1.0, 0.012, 6.0, 1.0]; T = 9.0
    res = {}
    def run(name, pol):
        log = P.episode(p, pol, T=T, seed=1); t = log[:, 0]; m = (t > 3.0) & (t < 5.0)
        res[name] = dict(mean_err_mid_ramp=float(np.mean(log[m, 3] - P.FD)), rmse_all=S.metrics(log)['rmse'])
    m_, b_ = 1.0, 300.
    run('constant_impedance_m1_b300_sigma0', duan_gains(m_, b_, 0., Tc))
    run('duan_m1_b300_sigma0.1', duan_gains(m_, b_, 0.1, Tc))
    run('duan_tuned', duan_g)
    run('song_tuned', song_pol)
    res['predicted_bias_constant_impedance'] = b_ * 1.5 * 0.012 / 6.0
    return res

def check_song_mapping():
    """numerical check of the E04 reduction: native NAIPC (Eqs. (14)-(15), forward Euler, 0.1 ms)
    vs the mapped PI admittance + lagged force path, on a spring environment with a moving
    surface (no delay). Returns max |F_native - F_mapped|."""
    md, bd, al, be, ga = 2., 200., 1., 5., 0.01; ke = 5000.; fd = 4.; dt = 1e-4; N = 40000
    def surf(t): return 0.004 * math.sin(2 * math.pi * 0.5 * t)
    # native
    y = fd / ke; v = 0.; d = 0.; fprev = 0.; Fn = []
    for i in range(N):
        t = i * dt; fe = ke * (y + surf(t)); ft = fd - fe; ftd = (ft - fprev) / dt if i else 0.; fprev = ft
        u = -d + ft - al * ft - ga * ftd; d += dt * be * (-u + ft - d)
        acc = (u - ft + bd * v) / md          # md xc'' + bd xc' + ft = u with xc = -y, xc' = -v
        v -= dt * acc; y += dt * v; Fn.append(fe)
    # mapped: M y'' + D y' = ef + g ef' + q int ef  (c = al + be ga)
    c = al + be * ga; M, D, q, g = md / c, bd / c, al * be / c, ga / c
    y = fd / ke; v = 0.; I = 0.; fprev = 0.; Fm = []
    for i in range(N):
        t = i * dt; fe = ke * (y + surf(t)); ef = fd - fe; efd = (ef - fprev) / dt if i else 0.; fprev = ef; I += ef * dt
        acc = (ef + g * efd + q * I - D * v) / M; v += dt * acc; y += dt * v; Fm.append(fe)
    return float(np.max(np.abs(np.array(Fn) - np.array(Fm))))

# ------------------------------------------------------------------ main
def _sum(rows): return U.summary(rows)

def main():
    t_all = time.time()
    gate = json.load(open('gate_PC125.json'))
    def gate_cell(task, dose):
        for c in gate:
            if c['task'] == task and c['dose'] == dose: return c
    out = {'platform': 'PC125', 'n_train': NTR, 'n_test': NTE, 'cells': []}
    pols_for_checks = {}
    for task, dose in CELLS:
        t0 = time.time()
        train = S.cohort(task, NTR, TR0, dose); test = S.cohort(task, NTE, TE0, dose)
        cell = {'task': task, 'dose': dose, 'methods': {}}
        tuned = {}
        for name, rule in [('EtR-A', 'truthk'), ('EtR-A-Haninger', 'haninger'), ('EtR-R', 'rls')]:
            info, pol = tune_etr(train, rule); tuned[name] = (info, pol)
            print(task, name, info, '%.0fs' % (time.time() - t0), flush=True)
        info, g = tune_duan(train); tuned['Duan2018'] = (info, g); print(task, 'Duan', info, '%.0fs' % (time.time() - t0), flush=True)
        info, pol = tune_song(train); tuned['Song2025'] = (info, pol); print(task, 'Song', info, '%.0fs' % (time.time() - t0), flush=True)
        names = list(tuned); res = evaluate(test, [tuned[n][1] for n in names], 'PC125', TE0)
        gc = gate_cell(task, dose)
        # B1 and B3b re-evaluated on the same test seeds for paired CIs
        from baselines_ref import TruthK
        rb = evaluate(test, [tuple(gc['B1']['gains']), TruthK(0.1, *gc['B3b']['0.1']['rule']), TruthK(0., *gc['B3b']['0.0']['rule'])], 'PC125', TE0)
        jB1 = [r[0]['J'] for r in rb]; jB3 = [r[1]['J'] for r in rb]
        for k, n in enumerate(names):
            rows = [r[k] for r in res]; jj = [m['J'] for m in rows]
            cell['methods'][n] = dict(tuning=tuned[n][0], test=_sum(rows), dJ_vs_B1=U.paired_ci(jj, jB1), dJ_vs_B3b01=U.paired_ci(jj, jB3),
                                      drmse_vs_B1=U.paired_ci([m['rmse'] for m in rows], [r[0]['rmse'] for r in rb]))
        cell['B1'] = dict(gains=gc['B1']['gains'], test_gate_json=gc['B1']['test'], test_rerun=_sum([r[0] for r in rb]))
        cell['B3b_0.1'] = dict(rule=gc['B3b']['0.1']['rule'], test_gate_json=gc['B3b']['0.1']['test'], test_rerun=_sum([r[1] for r in rb]))
        cell['B3b_0.0'] = dict(rule=gc['B3b']['0.0']['rule'], test_gate_json=gc['B3b']['0.0']['test'], test_rerun=_sum([r[2] for r in rb]))
        cell['seconds'] = time.time() - t0
        out['cells'].append(cell); pols_for_checks[task] = tuned
        json.dump(out, open('adaptive_results.json', 'w'), indent=1, default=float)
        for n in names + ['B1', 'B3b_0.1']:
            s = cell['methods'][n]['test'] if n in cell['methods'] else cell[n]['test_rerun']
            print(f"{task} {n:16s} J {s['J']:.3f} rmse {s['rmse']:.3f} peak {s['peak']:.2f} over {s['frac_peak_over']:.2f} inband {s['inband']:.3f}", flush=True)
    # known behaviour
    kT1 = 'T1' if 'T1' in pols_for_checks else list(pols_for_checks)[0]
    kT3 = 'T3' if 'T3' in pols_for_checks else list(pols_for_checks)[0]
    tk = pols_for_checks[kT1]['EtR-A']
    out['checks'] = {'estimator_T1tuned': check_estimator(tk[0]['est_params'], tk[1]),
                     'slope_drift': check_slope(pols_for_checks[kT3]['Duan2018'][1], pols_for_checks[kT3]['Song2025'][1]),
                     'song_mapping_max_abs_force_diff_N': check_song_mapping()}
    out['seconds_total'] = time.time() - t_all
    json.dump(out, open('adaptive_results.json', 'w'), indent=1, default=float)
    print(json.dumps(out['checks'], indent=1, default=float))
    pool().close(); pool().join()

if __name__ == '__main__':
    main()
