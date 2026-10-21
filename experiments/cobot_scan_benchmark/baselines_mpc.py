"""MPC-based variable-impedance baselines adapted to the benchmark's normal-axis
admittance + force structure (plant.py):  M yc'' + D yc' + K (yc - y_ref) = ef + ki I,
y_cmd = yc + KPF ef (+ u_ff),  ef = FD - F_filt.

Line references: E41 = literature/extracted_text/field_papers/E41.txt (Jin et al., TMECH 2023);
XUE = literature/extracted_text/xue.txt (Xue et al., RAS 2025).

1. JinMPC  (E41)  'os' one-step and 'ms' multi-step MPC-VIC.
   Jin's 1-D impedance m e'' + d e' + k e = f (E41:161) is our admittance with e = yc - y_ref and
   f = ef + ki I.  Jin's "novel impedance model" (E41:270-281, eq. 13-14) makes the impedance
   parameters the MPC input; here the inputs are u = (log D, KPF) with M, K, ki fixed
   (free-damping impedance constraints, E41:406-414 eq. 21 / E41:435-439 eq. 28, with the box mapped
   to plant.GMIN/GMAX: D in [20, 2500] N s/m, KPF in [0, 1e-3] m/N).
   Cost mapping (E41:329, eq. 15 / E41:409, eq. 26):  Jin's "precision" term s'Qs on the pose error
   becomes the predicted contact-force error ((F_pred - FD)/0.5 N)^2 plus a velocity term
   qv (yc'/0.01 m/s)^2 (Jin's Q also weights e', E41:647);  Jin's "compliance" term u'Ru (penalises
   impedance magnitude) becomes rc * normalised (log D - log Dmin)^2 (low admittance damping = more
   compliant to contact force);  a gain-change term rd * (Delta log D)^2 + rd * (Delta KPF/1e-4)^2
   regularises the receding-horizon sequence.
   Nominal anchor (deviation, documented): Jin's u'Ru is centred at zero impedance.  Centred there
   (family 'Jin_MS_literal': Dref = Dmin, rc = 0.01, M 0.3, K 0, ki 30) the MPC, fed by a causal
   stiffness estimate biased low on bone, drives KPF / D bang-bang and fails (train J >> B1).  The
   main families centre the compliance cost at the cell's B1 gains (Dref, KPref) and take M, K, ki
   from B1 (tuned on the same train cohort), i.e. MPC adaptation on top of the best fixed controller;
   rc large -> B1.  A per-update rate box |dlog D| <= 0.15, |dKPF| <= 1e-5 is part of Omega_t (Jin
   allows state-dependent impedance constraint sets, E41:373-380).
   Estimator gain kgain (tuned, {1, 1.7}): offline the RLS estimate is biased low by x0.4-0.6
   against the true instantaneous stiffness (closed-loop confound: scanning onto bone raises F while
   the robot retracts); kgain is a train-calibrated correction.  The state constraint of the MS version (E41:416-418)
   is the force ceiling F <= FLIM, softened with Jin's barrier eta = alpha (1 - exp(-phi/delta))
   (E41:488, eq. 32).
   Environment model.  Jin assumes the model is known (E41:845-850 leaves model learning to future
   work) and holds the external force constant over the step (E41:192 'f' = 0 assumption).  Force
   regulation needs a force prediction, so the predictor here is the known robot side (exact ZOH
   admittance, KPF path, ndel-cycle transport delay and first-order servo lag from the platform
   dict, 30 Hz sensor filter) closed through a linear contact model F = F0 + k_hat (y - y0)
   [+ k_hat * v_s * t when dist=1], with k_hat from causal recursive least squares on band-passed
   increments of the filtered force and the (same-filtered) encoder position, and v_s an EMA
   estimate of surface velocity.  dist=0 is the literal analogue of Jin's constant-force assumption;
   it fails under patient motion (the MPC brakes when the surface retracts, predicting an overshoot
   against a static surface), so the anchored families use dist=1.
   OS (E41:318-372, eq. 15-16, 23-25): one decision (log D, KPF) held for one update; predicted
   state s_{t+1} taken one update after the transport delay (ndel + per ticks); the residuals are
   linearised by finite differences (exact in KPF) and the resulting 2-variable box QP is solved
   exactly by active-set enumeration; 2 SQP passes.
   MS (E41:404-440, eq. 26-35): H = 10 updates of T = 16 ms (E41:746 uses H = 10, T = 15 ms) with a
   (log D, KPF) pair per update, tanh reparameterisation of the box (E41:454-464, eq. 29, applied
   to log D), barrier for the force ceiling; solved by Levenberg-Marquardt on the stacked residuals
   (finite-difference Jacobian, warm-started by shifting) instead of iLQR (same optimum class,
   different iteration).
   Update period: plant.episode(update=0.015) -> 2 ticks = 16 ms on PC125 (E41:646-648: T = 15 ms).
   Energy tank (E41:497-559, eq. 39-42): Gamma' = sigma1 (D yc'^2 - 0.5 ef^2 KPF'), sigma1 per (40);
   the only element that can inject energy is the KPF path (a virtual spring of compliance KPF
   between yc and y_cmd; Jin's variable stiffness k_v); a KPF increase is applied only when
   sigma2 = 1 per (42), otherwise KPF is held.  D variation is always dissipative in Jin's storage
   function (E41:490-492), so the tank never gates D.  Gamma0 = 0.1, Gamma_max = 1,
   Gamma_min = 0.01 J (E41:648).

2. XueMPVIC (XUE).  Their impedance m_i e'' + d_i(t) e' + k_i(t) e = e_f (XUE:220) is our admittance:
   m_i -> M (constant), d_i -> D, k_i -> K (admittance spring to the initial position y_ref,
   XUE:215-216 z_r' = 0), e_f -> ef.  Not in Xue: ki and KPF; they are kept as fixed tuned constants
   (ki adds ki * int(measured ef) to the model-based drive; KPF acts through plant's direct path).
   - Adaptive robust HC estimator (XUE:247-295, eq. 12-16): theta = [log k, n, c/k],
     phi = [1, log z, z'], y = log F; dead zone a = 1 if |e_y| > 2 Delta (eq. 15, Delta = 0.03,
     XUE:489); normalised-gradient update (16).  z is penetration in mm (units not stated in the
     paper; mm keeps phi well scaled).  Penetration needs the initial environment position:
     Xue's Assumption 1 (XUE:159-162) - known initial environment position, robot starts in contact.
     ys0_mode='A1' passes the scene's surface position at t = 0 (scene prior used at t = 0 only, the
     paper's assumption; NOT a causal signal); ys0_mode='nominal' uses y(0) - z_nom with z_nom the
     4 N indentation of a nominal 20 kPa tissue (fully causal).  z <= zmin (environment pushed past
     the initial surface, where the HC model is undefined; Xue stop the robot there, XUE:563-566)
     freezes the estimator and clamps z.
   - Modes (XUE:286-367, Table 1, eq. 24-26): mu(ef, ef') smooth-step indicator; safety mode
     K = 0, D = d_max; performance K = k_max, D = d_min; omega1 = 0.3 N (XUE:486), omega2 = efb/2
     (paper: efb = 1, omega2 = 0.5 N/s; efb is re-tuned because our filtered-force derivative noise
     is several N/s).
   - Passivity (XUE:308-331, eq. 27-30): constrained stiffness K_bar tracks mu*k_max with rate
     lambda = 5, energy tank V_tk with p1, p2 = 1e-3, eps = 1e-4 m; f_r -> FD + ki I.
     V_tk0 = 0.1 J (scaled to our impedance energies; paper's 5 is in unstated units).
   - MPC (XUE:368-412, eq. 32-36): linearised impedance + HC environment residual
     r = c1 e'' + c2 e' + c3 e + c4 with c1 = M, c2 = D + c z^n, c3 = K_bar + n (k + c z') z^(n-1),
     c4 = F_hat(z0) - FD - ki I - (environment slope terms at the measured state);  cost
     sum_{k=1..hp} r_k^2 over the Euler triple-integrator model (36) with jerk input, hp = 5, hc = 3
     (XUE:515-517), |e''| <= 5 m/s^2 (XUE:485).  The explicit-Euler model (36) is only valid for
     (d/m) Ts < 2 (paper: m 0.01, d <= 0.4, Ts 0.01 -> 0.4), so (M, d_min, d_max) pairs are restricted
     to (1, 50, 200) and (3, 150, 600).  The plant's u_ff realises the model's own step x1 = A x0 + B u.
     The estimated HC damping is projected to c >= 0 and (k + c z') clamped at 0 (eq. 9).
     f0src='model' linearises at the estimator output (literal); 'measured' at the measured force
     (the Delta -> 0 limit at the current state; the HC model still gives slope and prediction).
     ki (not in Xue) is restricted to {0, 3, 10}: the model-force loop has a slow coupled pole
     ~ kz / c2, and ki = 30 through it is unstable.  The MPC state is the commanded admittance state;
     the environment is linearised at the measured encoder state (so the predicted force accounts
     for commands still in flight).  Physical sign convention used (the extracted eq. 10-11 has a
     sign inconsistency that would give positive force feedback).  Realisation on the position
     interface (XUE:420-421 Remark 5: acceleration command = integral of the jerk): the optimal
     acceleration is integrated over one tick into the commanded position, which is sent through
     plant's u_ff offset while plant's own admittance is neutralised (M 0.25, D 2500, K 0, ki 0;
     its slow drift is replicated exactly and cancelled).  Runs every control tick (125 Hz; paper
     100 Hz, XUE:475).

Tuning (tune.best_fixed-comparable budget): grids of a few dozen settings on the train cohort
(seed SPLIT['train'], 20 scenes), selection = min mean J s.t. peak <= 7 N in >= 95 % of scenes.

Run:  python baselines_mpc.py [check|tune|all]   (single or two processes; module work guarded).
"""
import math, time, json, sys, itertools
import numpy as np
from numba import njit
import plant as P
import scenarios as S

FD = P.FD
DLO, DHI = math.log(P.GMIN[1]), math.log(P.GMAX[1])
KPMAX = float(P.GMAX[4])
NPROC = 2

# ----------------------------------------------------------------------------------------------
# numba core: closed-loop predictor of the plant's admittance + actuator + linear contact model
# st: [yc, v, I, y_prev, tgt_prev, Ff_prev, F0, kh, vs]   (state after tick n-1)
# q : in-flight actuator targets for ticks n .. n+ndel-1 (= y_cmd[n-ndel .. n-1])
# pr: [M, K, ki, yref, beta, af, Tc, VMAX, wf, wv, v0, flim, alpha_b, delta_b]
# ----------------------------------------------------------------------------------------------
@njit(cache=True)
def _predict(nt, per, logD, KP, st, q, pr, last_only, res):
    """simulate nt ticks; gains block b = min(i//per, nb-1). Writes 3 residuals per tick into res
    (or 3 residuals of the last tick only if last_only) and returns the number written."""
    M = pr[0]; K = pr[1]; ki = pr[2]; yref = pr[3]; beta = pr[4]; af = pr[5]; Tc = pr[6]
    VMAX = pr[7]; wf = pr[8]; wv = pr[9]; v0 = pr[10]; flim = pr[11]; ab = pr[12]; db = pr[13]
    yc = st[0]; v = st[1]; I = st[2]; y = st[3]; tgt = st[4]; Ff = st[5]; F0 = st[6]; kh = st[7]; vs = st[8]
    y0 = y; Fprev = F0
    nd = q.shape[0]
    qq = np.empty(nd + 1)
    for i in range(nd): qq[i] = q[i]
    nb = logD.shape[0]
    sf = math.sqrt(wf); sv = math.sqrt(wv)
    m = 0
    for i in range(nt):
        b = min(i // per, nb - 1)
        D = math.exp(logD[b]); kp = KP[b]
        y = tgt + (y - tgt) * beta                 # lag over the previous tick
        F = F0 + kh * (y - y0) + kh * vs * (i + 1) * Tc
        Ff = Ff + af * (0.5 * (F + Fprev) - Ff); Fprev = F
        ef = FD - Ff
        I = min(max(I + ef * Tc, -20.), 20.)
        u = ef + ki * I - K * (yc - yref)
        e = math.exp(-D * Tc / M)
        vn = v * e + (1 - e) * u / D
        vn = min(max(vn, -VMAX), VMAX)
        yc += 0.5 * (v + vn) * Tc; v = vn
        ycmd = yc + kp * ef
        if nd > 0:
            tgt = qq[0]
            for k in range(nd - 1): qq[k] = qq[k + 1]
            qq[nd - 1] = ycmd
        else:
            tgt = ycmd
        if (not last_only) or i == nt - 1:
            ov = F - flim
            bar = ab * (1.0 - math.exp(-ov / db)) if ov > 0 else 0.0
            res[m] = sf * (F - FD) / 0.5; res[m + 1] = sv * v / v0; res[m + 2] = math.sqrt(bar)
            m += 3
    return m


@njit(cache=True)
def _ms_resid(z, nb, per, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, out):
    logD = np.empty(nb); KP = np.empty(nb)
    for b in range(nb):
        logD[b] = DLO + (DHI - DLO) * 0.5 * (math.tanh(z[b]) + 1.0)          # E41 eq. (29), log D
        KP[b] = KPMAX * 0.5 * (math.tanh(z[nb + b]) + 1.0) if use_kpf else kpfix
    m = _predict(nb * per, per, logD, KP, st, q, pr, False, out)
    srd = math.sqrt(rd); src = math.sqrt(rc)
    for b in range(nb):
        pd = ldp if b == 0 else logD[b - 1]; pk = kpp if b == 0 else KP[b - 1]
        out[m] = srd * (logD[b] - pd); out[m + 1] = srd * (KP[b] - pk) / 1e-4
        out[m + 2] = src * (logD[b] - ldr); out[m + 3] = src * (KP[b] - kpr) / 1e-4
        m += 4
    return m


@njit(cache=True)
def _ms_solve(z, iters, nb, per, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr):
    nz = z.shape[0]
    nr = 3 * nb * per + 4 * nb
    r = np.empty(nr); r2 = np.empty(nr); J = np.empty((nr, nz))
    m = _ms_resid(z, nb, per, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r)
    c = np.sum(r[:m] ** 2)
    lam = 1e-2; h = 1e-4
    for it in range(iters):
        for j in range(nz):
            zj = z[j]; z[j] = zj + h
            _ms_resid(z, nb, per, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r2)
            z[j] = zj
            for i in range(m): J[i, j] = (r2[i] - r[i]) / h
        JJ = J[:m].T @ J[:m]; g = J[:m].T @ r[:m]
        improved = False
        for _ in range(6):
            A = JJ + lam * (np.diag(np.diag(JJ)) + 1e-6 * np.eye(nz))
            dz = -np.linalg.solve(A, g)
            zn = z + dz
            for j in range(nz): zn[j] = min(max(zn[j], -4.), 4.)
            _ms_resid(zn, nb, per, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r2)
            cn = np.sum(r2[:m] ** 2)
            if cn < c:
                z[:] = zn; r[:m] = r2[:m]; c = cn; lam = max(lam * 0.3, 1e-6); improved = True
                break
            lam *= 10.
        if not improved: break
    return c


@njit(cache=True)
def _boxqp2(H, g, lo, hi):
    """exact minimiser of 0.5 x'Hx + g'x on a 2-D box (convex) by active-set enumeration."""
    best = 1e300; bx0 = lo[0]; bx1 = lo[1]
    for s0 in range(3):
        for s1 in range(3):
            x = np.zeros(2); fr = np.zeros(2, dtype=np.bool_)
            if s0 == 0: fr[0] = True
            else: x[0] = lo[0] if s0 == 1 else hi[0]
            if s1 == 0: fr[1] = True
            else: x[1] = lo[1] if s1 == 1 else hi[1]
            if fr[0] and fr[1]:
                det = H[0, 0] * H[1, 1] - H[0, 1] * H[1, 0]
                if abs(det) < 1e-300: continue
                x[0] = (-g[0] * H[1, 1] + g[1] * H[0, 1]) / det
                x[1] = (-g[1] * H[0, 0] + g[0] * H[1, 0]) / det
            elif fr[0]:
                if H[0, 0] <= 0: continue
                x[0] = -(g[0] + H[0, 1] * x[1]) / H[0, 0]
            elif fr[1]:
                if H[1, 1] <= 0: continue
                x[1] = -(g[1] + H[1, 0] * x[0]) / H[1, 1]
            if x[0] < lo[0] - 1e-12 or x[0] > hi[0] + 1e-12 or x[1] < lo[1] - 1e-12 or x[1] > hi[1] + 1e-12:
                continue
            f = 0.5 * (x[0] * (H[0, 0] * x[0] + H[0, 1] * x[1]) + x[1] * (H[1, 0] * x[0] + H[1, 1] * x[1])) + g[0] * x[0] + g[1] * x[1]
            if f < best: best = f; bx0 = x[0]; bx1 = x[1]
    return bx0, bx1


@njit(cache=True)
def _os_resid(x, nt, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, out):
    logD = np.array([x[0]]); KP = np.array([x[1] if use_kpf else kpfix])
    m = _predict(nt, nt, logD, KP, st, q, pr, True, out)
    out[m] = math.sqrt(rd) * (x[0] - ldp); out[m + 1] = math.sqrt(rd) * (KP[0] - kpp) / 1e-4
    out[m + 2] = math.sqrt(rc) * (x[0] - ldr); out[m + 3] = math.sqrt(rc) * (KP[0] - kpr) / 1e-4
    return m + 4


@njit(cache=True)
def _os_solve(x, nt, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, passes, rateD, rateK, ldr, kpr):
    """Jin OS-MPC: residuals linearised at x (Gauss-Newton => QP of E41 eq. 23-25), exact 2-D box QP.
    Box = global impedance bounds intersected with a per-update rate box (state-dependent Omega_t)."""
    lo = np.array([max(DLO, ldp - rateD), max(0.0, kpp - rateK)])
    hi = np.array([min(DHI, ldp + rateD), min(KPMAX, kpp + rateK) if use_kpf else 0.0])
    if not use_kpf: x[1] = 0.0
    r = np.empty(7); r2 = np.empty(7); J = np.zeros((7, 2)); hs = np.array([1e-3, 1e-6])
    for _ in range(passes):
        _os_resid(x, nt, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r)
        for j in range(2):
            xj = x[j]; x[j] = xj + hs[j]
            _os_resid(x, nt, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r2)
            x[j] = xj
            for i in range(7): J[i, j] = (r2[i] - r[i]) / hs[j]
        # scale KPF to O(1)
        Js = J.copy(); Js[:, 1] *= 1e-4
        H = Js.T @ Js + 1e-9 * np.eye(2); g = Js.T @ r
        xs0 = np.array([x[0], x[1] / 1e-4])
        los = np.array([lo[0], lo[1] / 1e-4]); his = np.array([hi[0], hi[1] / 1e-4])
        # min 0.5 d'Hd + g'd over d in [lo - x0, hi - x0]
        d0, d1 = _boxqp2(H, g, los - xs0, his - xs0)
        x[0] = xs0[0] + d0; x[1] = (xs0[1] + d1) * 1e-4
    _os_resid(x, nt, st, q, pr, rd, rc, ldp, kpp, use_kpf, kpfix, ldr, kpr, r)
    return np.sum(r ** 2)


@njit(cache=True)
def _xue_mpc(c1, c2, c3, c4, e0, v0, a0, Ts, hp, hc, amax):
    """Xue eq. (35)-(36): min sum_{k=1..hp} (c1 a_k + c2 v_k + c3 e_k + c4)^2 over the Euler
    triple integrator, decisions a_1..a_hc (a_k = a_hc afterwards), |a_k| <= amax; bounded LS by
    active-set enumeration (3^hc). Returns a_1."""
    A = np.zeros((hp, hc)); b = np.zeros(hp)
    # e_k, v_k are affine in (a_0 known, a_1..a_hc); track coefficient vectors
    ce = np.zeros(hc); cv = np.zeros(hc); ca = np.zeros(hc)
    ke = e0; kv = v0; ka = a0
    for k in range(1, hp + 1):
        # step: e += Ts v ; v += Ts a ; a = a_k
        ke = ke + Ts * kv; ce = ce + Ts * cv
        kv = kv + Ts * ka; cv = cv + Ts * ca
        ca = np.zeros(hc); ca[min(k, hc) - 1] = 1.0; ka = 0.0
        for j in range(hc): A[k - 1, j] = c1 * ca[j] + c2 * cv[j] + c3 * ce[j]
        b[k - 1] = c1 * ka + c2 * kv + c3 * ke + c4
    best = 1e300; bestx = np.zeros(hc)
    ncomb = 3 ** hc
    for code in range(ncomb):
        cc = code; fixed = np.zeros(hc); isf = np.zeros(hc, dtype=np.bool_)
        for j in range(hc):
            s = cc % 3; cc //= 3
            if s == 0: isf[j] = True
            else: fixed[j] = -amax if s == 1 else amax
        nf = 0
        for j in range(hc):
            if isf[j]: nf += 1
        rhs = b.copy()
        for j in range(hc):
            if not isf[j]: rhs += A[:, j] * fixed[j]
        x = fixed.copy()
        if nf > 0:
            Af = np.zeros((hp, nf)); t = 0
            for j in range(hc):
                if isf[j]: Af[:, t] = A[:, j]; t += 1
            G = Af.T @ Af + 1e-12 * np.eye(nf)
            sol = np.linalg.solve(G, -(Af.T @ rhs))
            t = 0; ok = True
            for j in range(hc):
                if isf[j]:
                    x[j] = sol[t]; t += 1
                    if abs(x[j]) > amax + 1e-9: ok = False
            if not ok: continue
        rr = A @ x + b
        f = np.sum(rr ** 2)
        if f < best: best = f; bestx = x
    return bestx[0]


# ----------------------------------------------------------------------------------------------
# shared causal machinery
# ----------------------------------------------------------------------------------------------
class _Causal:
    """replicates the plant's admittance state from causal signals and applied gains; never reads
    obs.k_true or truth log columns (3, 4, 6, 7, 8)."""
    oracle = False

    def reset(self, pl):
        self.Tc = pl['Tc']; self.ndel = int(pl['ndel']); self.tau = pl['tau_a']; self.VMAX = pl['VMAX']
        dt = 5e-4; nc = int(round(self.Tc / dt))
        self.beta = (1 - dt / self.tau) ** nc if self.tau > 0 else 0.0
        afc = 1.0 - math.exp(-2 * math.pi * pl['fc'] * 1e-3)
        self.af = 1.0 - (1.0 - afc) ** int(round(self.Tc / 1e-3))
        self.m = 0; self.g = None; self.yc = None; self.v = 0.; self.I = 0.; self.yref = None
        self.tsum = 0.; self.tmax = 0.; self.ncall = 0; self.maxrep = 0.
        self._reset()

    def _sync(self, o):
        n = o.n
        if n == 0: return
        f = o.force; yc_ = o.y_cmd
        if self.yc is None:
            self.yref = float(o.y[0]); self.yc = self.yref; self.v = 0.; self.I = 0.
        M, D, K, ki, KPF, uff = self.g
        Tc = self.Tc
        for j in range(self.m, n):
            ef = FD - f[j]
            self.I = min(max(self.I + ef * Tc, -20.), 20.)
            u = ef + ki * self.I - K * (self.yc - self.yref)
            e = math.exp(-D * Tc / M)
            vn = self.v * e + (1 - e) * u / D
            vn = min(max(vn, -self.VMAX), self.VMAX)
            self.yc += 0.5 * (self.v + vn) * Tc; self.v = vn
            self.maxrep = max(self.maxrep, abs(self.yc + KPF * ef + uff - yc_[j]))
            self._tick(o, j, ef)
        self.m = n

    def _pred_state(self, o, F0, kh, vs):
        n = o.n; nd = self.ndel; yc_ = o.y_cmd
        tp = yc_[n - 1 - nd] if n - 1 - nd >= 0 else self.yref
        q = np.array([yc_[j] if j >= 0 else self.yref for j in range(n - nd, n)], float)
        st = np.array([self.yc, self.v, self.I, o.y[n - 1], tp, o.force[n - 1], F0, kh, vs])
        return st, q

    def __call__(self, o):
        t0 = time.perf_counter()
        if o.n > 0 and self.g is not None: self._sync(o)
        g = tuple(float(x) for x in self._decide(o))
        if len(g) == 5: g = g + (0.0,)
        self.g = tuple(np.clip(np.array(g), P.GMIN, P.GMAX))
        dt = time.perf_counter() - t0
        self.tsum += dt; self.ncall += 1; self.tmax = max(self.tmax, dt)
        return self.g


class StiffRLS:
    """causal linear contact-stiffness estimate from band-passed increments of the filtered force
    and the identically filtered encoder position (k in N/mm), RLS with exponential forgetting."""
    def __init__(self, Tc, af, tau_k, tau_hp=0.01, k0=1.5, kmin=0.3, kmax=40.):
        self.lam = math.exp(-Tc / tau_k); self.ahp = 1 - math.exp(-Tc / tau_hp); self.af = af
        self.k = k0; self.P = 1e4; self.kmin = kmin; self.kmax = kmax
        self.yf = None; self.Fp = None; self.my = 0.; self.mF = 0.
    def update(self, y_mm, F):
        if self.yf is None:
            self.yf = y_mm; self.Fp = F; return self.k
        yfn = self.yf + self.af * (y_mm - self.yf)
        dy = yfn - self.yf; dF = F - self.Fp; self.yf = yfn; self.Fp = F
        self.my += self.ahp * (dy - self.my); self.mF += self.ahp * (dF - self.mF)
        x = dy - self.my; z = dF - self.mF
        Px = self.P * x; den = self.lam + x * Px
        gk = Px / den
        self.k = min(max(self.k + gk * (z - self.k * x), self.kmin), self.kmax)
        self.P = min((self.P - gk * Px) / self.lam, 1e8)
        return self.k


class JinMPC(_Causal):
    """Jin et al. TMECH 2023 MPC-VIC adapted (see module docstring). mode 'os' | 'ms'."""
    def __init__(self, mode='ms', M=0.3, K=0., ki=30., qv=0.1, rd=0.1, rc=0.01, H=10, iters=4,
                 tau_k=0.2, dist=0, use_kpf=True, kpf0=1e-4, D0=100., flim=6.5, alpha_b=50.,
                 delta_b=0.5, tank=True, G0=0.1, Gmax=1.0, Gmin=0.01, rateD=0.15, rateK=1e-5, Dref=None,
                 KPref=None, kgain=1.0):
        # nominal impedance of the compliance cost (Jin's u'Ru is centred at zero impedance; here at
        # (Dref, KPref), default = the initial gains; Dref = Dmin, rc -> Jin's literal form)
        Dref = D0 if Dref is None else Dref; KPref = kpf0 if KPref is None else KPref
        self.__dict__.update(dict(rateD=rateD, rateK=rateK, ldr=math.log(Dref), kpr=KPref, kgain=kgain))
        self.__dict__.update(dict(mode=mode, M=M, K=K, ki=ki, qv=qv, rd=rd, rc=rc, H=H, iters=iters,
                                  tau_k=tau_k, dist=dist, use_kpf=use_kpf, kpf0=kpf0, D0=D0, flim=flim,
                                  alpha_b=alpha_b, delta_b=delta_b, tank=tank, G0=G0, Gmax=Gmax, Gmin=Gmin))

    def _reset(self):
        self.per = None
        self.est = StiffRLS(self.Tc, self.af, self.tau_k)
        self.vs = 0.; self.Ffp = None; self.ylast = None
        self.ld = math.log(self.D0); self.kp = self.kpf0 if self.use_kpf else 0.0
        self.z = None; self.G = self.G0; self.ngate = 0; self.nkup = 0
        self.hist = []
        self.pr = None

    def _tick(self, o, j, ef):
        y = o.y[j]; F = o.force[j]
        k = self.est.update(y * 1e3, F)          # N/mm
        if self.Ffp is not None:
            vsi = ((F - self.Ffp) / (k * 1e3) - (y - self.ylast)) / self.Tc
            self.vs += (1 - math.exp(-self.Tc / 0.1)) * (vsi - self.vs)
        self.Ffp = F; self.ylast = y
        # energy tank (E41 eq. 39-40) integrated per tick with the applied gains
        if self.tank:
            D = self.g[1]
            rate = D * self.v ** 2
            if self.G > self.Gmax and rate > 0: pass
            else: self.G += rate * self.Tc

    def _decide(self, o):
        n = o.n
        if n < 3 or self.g is None:
            if self.per is None and n > 0: pass
            self._last_n = n
            return (self.M, math.exp(self.ld), self.K, self.ki, self.kp)
        if self.per is None: self.per = max(n - self._last_n, 1)
        per = self.per
        if self.pr is None:
            self.pr = np.array([self.M, self.K, self.ki, self.yref, self.beta, self.af, self.Tc, self.VMAX,
                                1.0, self.qv, 0.01, self.flim, self.alpha_b, self.delta_b])
        f = o.force
        Fbar = f[n - 2] + (f[n - 1] - f[n - 2]) / self.af
        kh = self.est.k * 1e3 * self.kgain          # train-calibrated estimator gain (closed-loop bias)
        st, q = self._pred_state(o, Fbar, kh, self.vs if self.dist else 0.0)
        ldp, kpp = self.ld, self.kp
        if self.mode == 'os':
            x = np.array([ldp, kpp])
            _os_solve(x, self.ndel + per, st, q, self.pr, self.rd, self.rc, ldp, kpp, self.use_kpf, kpp, 2,
                      self.rateD, self.rateK, self.ldr, self.kpr)
            ld, kp = float(x[0]), float(x[1])
        else:
            H = self.H
            if self.z is None:
                zd = math.atanh(min(max(2 * (ldp - DLO) / (DHI - DLO) - 1, -0.999), 0.999))
                zk = math.atanh(min(max(2 * kpp / KPMAX - 1, -0.999), 0.999))
                self.z = np.array([zd] * H + [zk] * H)
            else:
                self.z[:H - 1] = self.z[1:H].copy(); self.z[H:2 * H - 1] = self.z[H + 1:].copy()
            _ms_solve(self.z, self.iters, H, per, st, q, self.pr, self.rd, self.rc, ldp, kpp,
                      self.use_kpf, kpp, self.ldr, self.kpr)
            ld = DLO + (DHI - DLO) * 0.5 * (math.tanh(self.z[0]) + 1)
            kp = KPMAX * 0.5 * (math.tanh(self.z[H]) + 1) if self.use_kpf else 0.0
            # rate box of Omega_t on the applied move (the planned sequence is penalised by rd)
            ld = min(max(ld, ldp - self.rateD), ldp + self.rateD)
            kp = min(max(kp, kpp - self.rateK), kpp + self.rateK)
        # tank gating of energy-injecting KPF increases (E41 eq. 41-42)
        if self.tank and kp > kpp:
            ef = FD - f[n - 1]; inj = 0.5 * ef * ef * (kp - kpp)
            self.nkup += 1
            if self.G - inj < self.Gmin:
                kp = kpp; self.ngate += 1
            else:
                self.G -= inj
        self.ld, self.kp = ld, kp
        self.hist.append((o.t[n - 1], math.exp(ld), kp, self.est.k))
        return (self.M, math.exp(ld), self.K, self.ki, kp)


def _sig(x, a, b):
    """Xue eq. (24)."""
    if x < a: return 0.0
    if x > b: return 1.0
    w = b - a
    return 0.5 + 0.5 * math.sin(math.pi / w * (x - b) + math.pi / 2)


Z_NOM_MM = None
def _znom_mm():
    """4 N static indentation of a nominal homogeneous 20 kPa tissue (causal prior for ys0)."""
    global Z_NOM_MM
    if Z_NOM_MM is None:
        p = np.zeros(P.NP); p[0] = 2e4; p[13] = 1.0
        lo, hi = 1e-7, 0.03
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if P.footprint(P.X0, mid, p)[0] > FD: hi = mid
            else: lo = mid
        Z_NOM_MM = 0.5 * (lo + hi) * 1e3
    return Z_NOM_MM


class XueMPVIC(_Causal):
    """Xue et al. RAS 2025 MPVIC adapted (see module docstring). Runs every control tick."""
    BUILTIN = (0.25, 2500., 0., 0.)          # neutralised plant admittance (M, D, K, ki)

    def __init__(self, ys0_mode='A1', M=1.0, dmin=100., dmax=1000., kmax=0., ki=30., KPF=1e-4,
                 w1=0.3, efb=20., w2=None, Delta=0.03, lam=5., p2=1e-3, Vtk0=0.1, eps=1e-4,
                 hp=5, hc=3, amax=5., zmin_mm=0.3, fdot_tau=0.02, n0=1.3, f0src='model'):
        self.__dict__.update(dict(ys0_mode=ys0_mode, M=M, dmin=dmin, dmax=dmax, kmax=kmax, ki=ki, KPF=KPF,
                                  w1=w1, efb=efb, w2=(efb / 2 if w2 is None else w2), Delta=Delta, lam=lam,
                                  p2=p2, Vtk0=Vtk0, eps=eps, hp=hp, hc=hc, amax=amax, zmin_mm=zmin_mm,
                                  fdot_tau=fdot_tau, n0=n0, f0src=f0src))
        self.ys0_scene = None

    def bind(self, p):
        """Xue Assumption 1 (XUE:159-162): the initial environment position is known. Scene prior at
        t = 0 only (surface height at x0, t = 0 in the penetration coordinate)."""
        self.ys0_scene = -float(P.surface(P.X0, 0.0, np.asarray(p, float))[0])

    def _reset(self):
        self.th = None; self.phi_prev = None; self.ys0 = None
        self.e = 0.; self.ed = 0.; self.edd = 0.            # MPC (commanded admittance) state, e = yc - yref
        self.Ie = 0.; self.efd = 0.; self.efp = None
        self.Kbar = 0.; self.kun_prev = 0.; self.Vtk = self.Vtk0
        self.ycb = None; self.vb = 0.; self.zprev = None
        self.hist = []; self.nupd = 0; self.nest = 0; self.nfreeze = 0; self.mu = 0.; self.phi_valid = False

    def _sync(self, o):
        # replicate the neutralised built-in admittance (yc_b) instead of the generic replica
        n = o.n
        if n == 0: return
        f = o.force
        if self.ycb is None:
            self.yref = float(o.y[0]); self.ycb = self.yref; self.vb = 0.; self.yc = self.yref
        M, D, K, ki, KPF, uff = self.g
        Tc = self.Tc
        for j in range(self.m, n):
            ef = FD - f[j]
            u = ef - K * (self.ycb - self.yref)
            e = math.exp(-D * Tc / M)
            vn = self.vb * e + (1 - e) * u / D
            vn = min(max(vn, -self.VMAX), self.VMAX)
            self.ycb += 0.5 * (self.vb + vn) * Tc; self.vb = vn
            self.maxrep = max(self.maxrep, abs(self.ycb + KPF * ef + uff - o.y_cmd[j]))
        self.m = n

    def _decide(self, o):
        n = o.n; Tc = self.Tc
        if n == 0:
            return self.BUILTIN + (self.KPF, 0.0)
        f = o.force; y = o.y
        if self.ys0 is None:
            if self.ys0_mode == 'A1' and self.ys0_scene is not None: self.ys0 = self.ys0_scene
            else: self.ys0 = float(y[0]) - _znom_mm() * 1e-3
        F = f[n - 1]
        z = (y[n - 1] - self.ys0) * 1e3                                  # mm
        zd = (y[n - 1] - y[n - 2]) * 1e3 / Tc if n > 1 else 0.0          # mm/s
        valid = z > self.zmin_mm
        zc = max(z, self.zmin_mm)
        phi = np.array([1.0, math.log(zc), zd])
        if self.th is None:
            self.th = np.array([math.log(max(F, 0.1)) - self.n0 * math.log(zc), self.n0, 0.0])
        # ---- adaptive robust estimator, XUE eq. (12)-(16) ----
        if self.phi_prev is not None and F > 0.05 and self.phi_valid:
            ey = math.log(F) - float(self.th @ self.phi_prev)
            if abs(ey) > 2 * self.Delta:
                self.th = self.th + self.phi_prev * ey / (1.0 + float(self.phi_prev @ self.phi_prev))
                self.th[2] = max(self.th[2], 0.0)       # projection: HC damping c >= 0 (physical)
                self.nupd += 1
            self.nest += 1
        if not valid: self.nfreeze += 1
        self.phi_prev = phi; self.phi_valid = valid
        logk, nn, ck = float(self.th[0]), float(self.th[1]), float(self.th[2])
        nn = min(max(nn, 0.3), 4.0)
        kk = math.exp(logk); cc = ck * kk                                # F = (k + c z') z^n, z in mm
        # ---- modes, XUE Table 1 / eq. (24)-(26) ----
        ef = FD - F
        if self.efp is not None:
            self.efd += (1 - math.exp(-Tc / self.fdot_tau)) * ((ef - self.efp) / Tc - self.efd)
        self.efp = ef
        mu = (_sig(ef, -self.w1, 0.) * (1 - _sig(ef, FD - self.w1, FD)) *
              _sig(self.efd, -self.efb, -self.efb + self.w2) * (1 - _sig(self.efd, self.efb - self.w2, self.efb)))
        self.mu = mu
        d = self.dmax + mu * (self.dmin - self.dmax)
        kun = mu * self.kmax
        # ---- constrained stiffness and tank, XUE eq. (27)-(30); e_z -> yc - yref ----
        self.Ie = min(max(self.Ie + ef * Tc, -20.), 20.)
        fr = FD + self.ki * self.Ie
        ez = self.e; ezd = self.ed
        ek = kun - self.Kbar
        kund = (kun - self.kun_prev) / Tc; self.kun_prev = kun
        if abs(ez) <= self.eps:
            kbd = kund + self.lam * ek
        else:
            p1 = self.lam * ek * ez * ez / 2 if self.Vtk > 0 else 0.0
            kbd = min(kund + self.lam * ek, 2 * (d * ezd * ezd - fr * ezd + p1) / (ez * ez))
        self.Kbar = max(self.Kbar + kbd * Tc, 0.0)
        self.Vtk = max(self.Vtk + ((d - self.p2) * ezd * ezd - fr * ezd - 0.5 * kbd * ez * ez) * Tc, 0.0)
        # ---- linearised MPC, XUE eq. (32)-(36) ----
        zm = zc * 1e-3
        kc = max(kk + cc * zd, 0.0)                                        # HC eq. (9): max(., 0)
        Fhat = kc * zc ** nn                                               # estimator output at current state
        F0 = Fhat if self.f0src == 'model' else F   # 'measured': linearise at the measured force (Delta->0 limit)
        kz = nn * kc * zc ** (nn - 1) * 1e3                                # N/m
        cz = cc * zc ** nn * 1e3                                           # N s/m
        eact = y[n - 1] - self.yref
        c1 = self.M; c2 = d + cz; c3 = self.Kbar + kz
        c4 = F0 - fr - kz * eact - cz * zd * 1e-3
        a1 = _xue_mpc(c1, c2, c3, c4, self.e, self.ed, self.edd, Tc, self.hp, self.hc, self.amax)
        # realise one step of the MPC model (36): x1 = A x0 + B u with jerk u = (a1 - a0)/Ts, i.e. the
        # acceleration command is the integral of the jerk (Remark 5) and acts from the next step
        self.e += Tc * self.ed
        self.ed = min(max(self.ed + Tc * self.edd, -self.VMAX), self.VMAX)
        self.edd = a1
        # command through u_ff, cancelling the neutralised built-in admittance's predicted state
        Mb, Db, Kb, kib = self.BUILTIN
        efn = ef
        ee = math.exp(-Db * Tc / Mb)
        vbn = self.vb * ee + (1 - ee) * (efn - Kb * (self.ycb - self.yref)) / Db
        ycb_next = self.ycb + 0.5 * (self.vb + vbn) * Tc
        uff = self.yref + self.e - ycb_next
        if n % 4 == 0: self.hist.append((o.t[n - 1], d, self.Kbar, mu, kk, nn, ck, Fhat, z))
        return (Mb, Db, Kb, kib, self.KPF, uff)

    def __call__(self, o):
        t0 = time.perf_counter()
        if o.n > 0 and self.g is not None: self._sync(o)
        g = self._decide(o)
        self.g = tuple(float(x) for x in np.clip(np.array(g, float), P.GMIN, P.GMAX))
        dt = time.perf_counter() - t0
        self.tsum += dt; self.ncall += 1; self.tmax = max(self.tmax, dt)
        return self.g


# ----------------------------------------------------------------------------------------------
# evaluation (own small pool: <= 2 processes)
# ----------------------------------------------------------------------------------------------
CLASSES = {'JinMPC': JinMPC, 'XueMPVIC': XueMPVIC}
UPDATE = {'JinMPC': 0.015, 'XueMPVIC': None}


def _job(args):
    p, T, spec, seed, keep = args
    name, kw = spec
    pol = CLASSES[name](**kw)
    if hasattr(pol, 'bind'): pol.bind(p)
    log = P.episode(p, pol, T=T, platform='PC125', seed=seed, update=UPDATE[name])
    m = S.metrics(log, Tc=P.PLATFORMS['PC125']['Tc'])
    m['t_call_mean'] = pol.tsum / max(pol.ncall, 1); m['t_call_max'] = pol.tmax; m['replica_err'] = pol.maxrep
    if name == 'JinMPC':
        h = np.array(pol.hist) if pol.hist else np.zeros((1, 4))
        m['D_med'] = float(np.median(h[:, 1])); m['KPF_med'] = float(np.median(h[:, 2]))
        m['khat_med'] = float(np.median(h[:, 3])); m['tank_gate'] = pol.ngate; m['kpf_up'] = pol.nkup
    else:
        h = np.array(pol.hist) if pol.hist else np.zeros((1, 9))
        m['D_med'] = float(np.median(h[:, 1])); m['mu_mean'] = float(np.mean(h[:, 3]))
        m['est_upd_frac'] = pol.nupd / max(pol.nest, 1); m['freeze_frac'] = pol.nfreeze / max(pol.ncall, 1)
    if keep:
        m['hist'] = pol.hist
        m['log'] = log
    return m


_pool = None
def run(scenes, specs, seed0, procs=NPROC, keep=False):
    """metrics[scene][spec]"""
    global _pool
    jobs = [(p, T, s, seed0 + i, keep) for i, (p, T) in enumerate(scenes) for s in specs]
    if procs > 1:
        if _pool is None:
            from multiprocessing import Pool
            _pool = Pool(procs)
        res = _pool.map(_job, jobs, chunksize=max(1, len(jobs) // (8 * procs)))
    else:
        res = [_job(j) for j in jobs]
    ns = len(specs)
    return [res[i * ns:(i + 1) * ns] for i in range(len(scenes))]


def summary(rows):
    keys = ['J', 'rmse', 'p95', 'maxerr', 'inband', 'peak', 't_over', 'loss', 'loss_max', 'osc', 'shear', 'diverged']
    out = {k: float(np.nanmean([r[k] for r in rows])) for k in keys}
    pk = np.array([r['peak'] for r in rows]); out['frac_peak_over'] = float(np.mean(~(pk <= FD + 3.0)))
    for k in ['t_call_mean', 't_call_max', 'D_med', 'KPF_med', 'khat_med', 'mu_mean', 'est_upd_frac',
              'freeze_frac', 'tank_gate', 'replica_err']:
        if k in rows[0]: out[k] = float(np.nanmean([r[k] for r in rows])) if k != 't_call_max' else float(np.max([r[k] for r in rows]))
    return out


def select(res):
    Jm = np.array([[m['J'] for m in row] for row in res]).mean(0)
    pk = np.array([[m['peak'] if np.isfinite(m['peak']) else 1e9 for m in row] for row in res])
    feas = (pk <= FD + 3.0).mean(0) >= 0.95
    j = int(np.argmin(np.where(feas, Jm, Jm + 1e6)))
    return j, Jm, feas


def paired_ci(a, b, n=10000, seed=0):
    d = np.asarray(a) - np.asarray(b); rng = np.random.default_rng(seed)
    bs = d[rng.integers(0, len(d), (n, len(d)))].mean(1)
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


# tuning grids (few dozen settings each)
def _b1(task, dose):
    g = gate_ref(task, dose)['B1']['gains']
    return dict(M=g[0], K=g[2], ki=g[3], D0=g[1], kpf0=g[4])

def grid_os(task, dose):
    b = _b1(task, dose)
    return [('JinMPC', dict(mode='os', qv=qv, rd=rd, rc=rc, kgain=kg, dist=1, **b))
            for qv, rd, rc, kg in itertools.product([0.0, 0.1], [0.1, 1.0], [0.3, 3.0], [1.0, 1.7])]

def grid_ms(task, dose):
    b = _b1(task, dose)
    return [('JinMPC', dict(mode='ms', qv=0.1, rd=rd, rc=rc, kgain=kg, tau_k=tk, dist=1, **b))
            for rd, rc, kg, tk in itertools.product([0.1, 1.0], [1.0, 3.0, 10.0], [1.0, 1.7], [0.1, 0.3])]

def grid_ms_literal(task, dose):
    """no nominal anchor: Jin's compliance cost toward the minimum impedance (Dref = Dmin, rc small),
    M = 0.3, K = 0, ki = 30 (the B3b rule family's fixed part)."""
    return [('JinMPC', dict(mode='ms', M=0.3, K=0., ki=30., D0=100., kpf0=1e-4, Dref=20., KPref=1e-4,
                            rc=0.01, qv=0.1, rd=rd, dist=ds, kgain=kg))
            for rd, ds, kg in itertools.product([0.1, 1.0], [0, 1], [1.0, 1.7])]

def grid_xue(mode='A1'):
    return [('XueMPVIC', dict(ys0_mode=mode, M=M, dmin=dm, dmax=dx, ki=ki, efb=efb, f0src=fs))
            for (M, dm, dx), ki, efb, fs in itertools.product([(1., 50., 200.), (3., 150., 600.)], [0., 3., 10.],
                                                              [10., 40.], ['model', 'measured'])]


CELLS = [('T3', dict(ramp=0.5, A=0.001)), ('T1', dict(ratio=30., A=0.001))]


def gate_ref(task, dose):
    g = json.load(open('gate_PC125.json'))
    for r in g:
        if r['task'] == task and r['dose'] == dose:
            return r
    return None


def tune_cell(task, dose, methods, ntr=20, nte=40, procs=NPROC, log=print):
    train = S.cohort(task, ntr, S.SPLIT['train'], dose); test = S.cohort(task, nte, S.SPLIT['test'], dose)
    out = {'task': task, 'dose': dose}
    ref = gate_ref(task, dose)
    for meth, grid in methods:
        t0 = time.time()
        if isinstance(grid, tuple):          # ('from', other_method, overrides): test-only variant
            src = out[grid[1]]; spec = (src['spec'][0], dict(src['spec'][1], **grid[2]))
            grid = [spec]; j = 0; Jm = np.array([np.nan]); feas = np.array([False])
        else:
            rt = run(train, grid, S.SPLIT['train'], procs)
            j, Jm, feas = select(rt)
        rs = run(test, [grid[j]], S.SPLIT['test'], procs)
        rows = [r[0] for r in rs]
        o = {'spec': grid[j], 'params': grid[j][1], 'train_J': float(Jm[j]), 'train_feasible': bool(feas[j]),
             'n_feasible': int(feas.sum()), 'n_grid': len(grid), 'test': summary(rows),
             'train_grid_J': [float(x) for x in Jm], 'seconds': time.time() - t0,
             'test_J': [r['J'] for r in rows]}
        if ref is not None:
            o['dJ_vs_B1'] = None  # filled below when B1 per-scene J available
        out[meth] = o
        s = o['test']
        log(f"{task} {dose} {meth}: {grid[j][1]} trainJ {Jm[j]:.3f} feas {feas[j]} | test J {s['J']:.3f} rmse {s['rmse']:.3f} "
            f"peak {s['peak']:.2f} over {s['frac_peak_over']:.2f} inband {s['inband']:.3f} t/call {s['t_call_mean']*1e3:.3f} ms "
            f"({o['seconds']:.0f}s)")
    if ref is not None:
        # recompute B1 and B3b(0.1 s) per-scene test J on the same seeds for paired CIs
        from baselines_ref import TruthK
        g1 = tuple(ref['B1']['gains']); r3 = ref['B3b']['0.1']['rule']
        jb1 = []; jb3 = []
        for i, (p, T) in enumerate(test):
            l1 = P.episode(p, g1, T=T, platform='PC125', seed=S.SPLIT['test'] + i)
            jb1.append(S.metrics(l1)['J'])
            l3 = P.episode(p, TruthK(0.1, *r3), T=T, platform='PC125', seed=S.SPLIT['test'] + i)
            jb3.append(S.metrics(l3)['J'])
        out['B1'] = {'gains': ref['B1']['gains'], 'test': ref['B1']['test'], 'test_J_recomputed_mean': float(np.mean(jb1))}
        out['B3b_0.1'] = {'rule': r3, 'test': ref['B3b']['0.1']['test'], 'test_J_recomputed_mean': float(np.mean(jb3))}
        for meth, _ in methods:
            out[meth]['dJ_vs_B1'] = paired_ci(out[meth]['test_J'], jb1)
            out[meth]['dJ_vs_B3b_0.1'] = paired_ci(out[meth]['test_J'], jb3)
    return out


def checks(procs=1, log=print, specs=None):
    """known-behaviour checks on homogeneous tissue (N1-type scenes, breathing 1 mm)."""
    out = {}
    rng = np.random.default_rng(7)
    scenes = []
    for E in [1.5e4, 1.0e6]:
        for k in range(3):
            p, T = S.make('N1', rng, dict(E=E, A=0.001)); scenes.append((p, T))
    if specs is None:
        specs = [('JinMPC', dict(mode='ms')), ('JinMPC', dict(mode='os')), ('JinMPC', dict(mode='ms')),
                 ('XueMPVIC', dict(ys0_mode='A1')), ('XueMPVIC', dict(ys0_mode='nominal'))]
    specs = [(a, dict(b)) for a, b in specs]
    res = run(scenes, specs, 500, procs=procs, keep=True)
    names = ['jin_ms', 'jin_os', 'jin_ms_literal', 'xue_A1', 'xue_nominal']
    for si, (p, T) in enumerate(scenes):
        tag = 'soft_15kPa' if p[0] < 1e5 else 'stiff_1MPa'
        for ci, nm in enumerate(names):
            m = res[si][ci]; log_ = m['log']
            kt = float(np.median(log_[log_[:, 0] > 1, 8]))       # truth used for the CHECK only
            d = {'J': m['J'], 'rmse': m['rmse'], 'peak': m['peak'], 'D_med': m['D_med'], 'k_true_med_Npm': kt,
                 'replica_err': m['replica_err']}
            if nm.startswith('jin'):
                h = np.array(m['hist']); d['khat_med_Npm'] = float(np.median(h[h[:, 0] > 1, 3]) * 1e3)
                d['KPF_med'] = m['KPF_med']
            else:
                h = np.array(m['hist'])     # t, d, Kbar, mu, k, n, c/k, Fhat, z
                late = h[:, 0] > T - 2.0; mid = (h[:, 0] > 1.0) & (h[:, 0] < 2.5)
                # HC local stiffness n F/z vs truth k_eff (instantaneous incl. SLS factor)
                ks = h[:, 5] * h[:, 7] / np.maximum(h[:, 8], 1e-3) * 1e3
                d.update({'theta_late_mean': [float(np.mean(h[late, 4])), float(np.mean(h[late, 5])), float(np.mean(h[late, 6]))],
                          'theta_late_std': [float(np.std(h[late, 4])), float(np.std(h[late, 5])), float(np.std(h[late, 6]))],
                          'logk_drift_mid_to_late': float(np.mean(np.log(h[late, 4])) - np.mean(np.log(h[mid, 4]))),
                          'Fhat_vs_F_rms_late': float(np.sqrt(np.mean((h[late, 7] - np.interp(h[late, 0], log_[:, 0], log_[:, 5])) ** 2))),
                          'local_k_HC_med_Npm': float(np.median(ks[h[:, 0] > 1])),
                          'est_upd_frac': m['est_upd_frac'], 'mu_mean': m['mu_mean']})
            out.setdefault(nm, {}).setdefault(tag, []).append(d)
    # concise printout
    for nm in names:
        for tag in ['soft_15kPa', 'stiff_1MPa']:
            L = out[nm][tag]
            log(f"check {nm} {tag}: D_med {np.mean([x['D_med'] for x in L]):.0f} rmse {np.mean([x['rmse'] for x in L]):.3f} "
                f"peak {np.mean([x['peak'] for x in L]):.2f} k_true {np.mean([x['k_true_med_Npm'] for x in L]):.0f} "
                + (f"khat {np.mean([x['khat_med_Npm'] for x in L]):.0f}" if nm.startswith('jin') else
                   f"kHC {np.mean([x['local_k_HC_med_Npm'] for x in L]):.0f} Fhat-F {np.mean([x['Fhat_vs_F_rms_late'] for x in L]):.3f} "
                   f"theta_std_late {np.round(np.mean([x['theta_late_std'] for x in L],0),4).tolist()} drift {np.mean([x['logk_drift_mid_to_late'] for x in L]):.3f}")
                + f" repl {max(x['replica_err'] for x in L):.1e}")
    return out


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    try:
        results = json.load(open('mpc_results.json'))
    except Exception:
        results = {}
    if what in ('tune', 'all'):
        cells = []
        for task, dose in CELLS:
            methods = [('Jin_OS', grid_os(task, dose)), ('Jin_MS', grid_ms(task, dose)),
                       ('Jin_MS_literal', grid_ms_literal(task, dose)), ('Xue_A1', grid_xue('A1')),
                       ('Xue_nominal_ys0', ('from', 'Xue_A1', dict(ys0_mode='nominal')))]
            cells.append(tune_cell(task, dose, methods))
            results['cells'] = cells
            json.dump(results, open('mpc_results.json', 'w'), indent=1, default=float)
    if what in ('check', 'all'):
        specs = None
        if 'cells' in results:
            c0 = results['cells'][0]
            specs = [tuple(c0['Jin_MS']['spec']), tuple(c0['Jin_OS']['spec']), tuple(c0['Jin_MS_literal']['spec']),
                     tuple(c0['Xue_A1']['spec']), tuple(c0['Xue_nominal_ys0']['spec'])]
        results['checks'] = checks(procs=NPROC, specs=specs)
        json.dump(results, open('mpc_results.json', 'w'), indent=1, default=float)
