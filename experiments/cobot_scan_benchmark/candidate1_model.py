"""Candidate 1 (model / optimisation route): moving-window contact estimation that separates surface
motion from stiffness, a finite-horizon OCP over impedance and feedforward through the exact robot-side
predictor, and a warm-started / neighbouring update. Causal policy for plant.episode.

Interface. reset(platform_dict); __call__(obs) -> g = (M, D, K, ki, KPF, u_ff) every update period
(final: 8 ms = every control tick on PC125). Reads only obs.t, obs.force (filtered F/T), obs.y (encoder),
obs.y_cmd and the tick index (obs.log columns 0, 2, 5, 11). Never reads k_true, log columns 3, 4, 6, 7, 8 or
the scene vector p. (The attributes _vs/_as are an oracle diagnostic hook used only in offline scratch
tests, never in any reported run.) Robot-side constants used: Tc, ndel, tau_a, fc, VMAX of plant.PLATFORMS
(calibratable) and the shared admittance law; the 1 kHz F/T sampling.

Plant convention (plant.py): penetration y (= -z), ef = FD - F_filt,
    M yc'' + D yc' + K (yc - y_ref) = ef + ki I,  I = int ef dt,   y_cmd = yc + KPF ef + u_ff.

1. Estimator (moving-horizon, arrival-cost regularised, constrained Gauss-Newton; numba _est)
   Robot position between ticks is reconstructed from the encoder sample and the known in-flight target
   (servo first-order lag, ndel transport delay):  y(t_j + m ms) = y*_j + (y_j - y*_j) e^{-m ms/tau_a},
   y*_j = y_cmd[j - ndel].  The F/T filter (first-order, a_f = 1 - e^{-2 pi fc 1ms}, sampled at 1 kHz) is
   inverted one step:  z_j = F_{j+1} - (1-a_f)^ns F_j = sum_{m=1..ns} a_f (1-a_f)^{ns-m} f(t_j + m ms) + noise,
   ns = Tc / 1 ms.  Local contact model with window-end reference (t_e, y_e = last logged tick):
       f(t) = g0 + k (y(t) - y_e) + k [ vs tau + as tau^2 / 2 ],   tau = t - t_e,   k = exp(kappa)
   i.e. stiffness k (slope with respect to robot position) and the equivalent surface trajectory
   h(tau) = vs tau + as tau^2/2 (patient motion, plus any force change at fixed robot position such as the
   footprint moving onto bone, which is equivalent for force prediction) plus the force level g0 (absorbs
   the sensor bias and the SLS relaxation level). theta = (kappa, g0, vs, as) minimises
       sum_j lam^{j1-j} ((z_j - zhat_j(theta)) / sigma_z)^2 + sum_{i in kappa,vs,as} (theta_i - thetabar_i)^2 / (s_i^2 per)
   over the last Nw intervals (samples with F < 0.5 N excluded), arrival prior thetabar = previous estimate
   (vs propagated with as), random-walk scales s_k (log k), s_v, s_a per update; no prior on g0.  Bounds:
   k in [200, 3e5] N/m, |vs| <= 0.2 m/s, |as| <= 5 m/s^2, g0 in [-5, 30] N; 3 projected GN iterations
   warm-started at the prior.  Optional conservative release (from EtR-A): k_used = max(k_hat,
   k_used e^{-dt/trel}).  Optional excitation dither (from EtR-A) added to u_ff:
   a sin(2 pi n / 6) with a = clip(Fdith / k_used, 2 um, 60 um) (20.8 Hz on PC125); it is a plant command,
   so it is inside the force metric; it is known to the predictor.  No viscous (Hunt-Crossley) or explicit SLS
   term: both were negligible/absorbed in the window (SLS relaxation shows as drift).
   Ablation (a) 'drift=0' fixes vs = as = 0 (stiffness + level only).

2. Decision: finite-horizon OCP (numba _roll/_resid), solved every update
   Decision x = (log D, KPF, u0, u1); M, K, ki held at tuned values.  Feedforward over the horizon:
   u_ff(t_{n+i}) = u0 + u1 min(i Tc, t_u) (ramp then hold) + known dither.
   Predictor = exact robot side: ZOH admittance with velocity cap and anti-windup, KPF path, ndel delay
   buffer, first-order servo lag at 1 ms, 1 kHz F/T filter, started from the replicated admittance state
   (yc, yc', I) and the last logged (y, F_filt, in-flight commands), closed through the estimated contact
       f = max(0, g0 + k_used (y - y_e + hhat(tau))),  hhat = vs min(tau, t_v) + g_a as min(tau, t_v)^2 / 2.
   Cost J(x; c) = 1/2 ||r(x; c)||^2, r stacking, for each scenario s with weight w_s and each predicted tick
   i = 1..H:  (f_i - FD)/0.5/sqrt(H),  sqrt(w_c) max(0, f_i - F_c)/sqrt(H) (F_c = 6.5 N, 7 N ceiling margin),
   sqrt(w_l) max(0, F_l - f_i)/0.5/sqrt(H) (F_l = 1 N, contact loss); a terminal admittance-velocity term
   sqrt(w_v)(yc'_H + hhat'_H)/0.01 (w_v = 0 in all tuned runs); and, once: gain-rate
   sqrt(r_D)(log D - log D_prev), sqrt(r_K)(KPF - KPF_prev)/1e-4, feedforward-rate sqrt(r_u)(u0 - u_prev)/1 mm,
   gain prior sqrt(p_D)(log D - log D_r), sqrt(p_K)(KPF - KPF_r)/1e-4 with the B3b-form rule evaluated at the
   estimate D_r = max(Dmin, k_used/wI), KPF_r = min(cP/k_used, 1e-3), feedforward prior sqrt(p_u) u0/1 mm,
   sqrt(p_u1) u1/(10 mm/s).  Scenarios: s = 0 nominal (current state, estimated drift extrapolated);
   s = 1, 2 robustness probes from equilibrium (F = FD at rest, no feedforward) against an unpredicted surface
   ramp of speed v_d on contacts of stiffness k_used and rho k_used (weight 'probe').  The probes give the gain
   choice a disturbance-rejection objective that the feedforward cannot cancel.
   Bounds: plant.GMIN/GMAX intersected with per-update rate limits
   log D in [log D_prev - 0.15 per, log D_prev + 2 per], KPF in [KPF_prev - 1e-3 per, KPF_prev + 1e-4 per]
   (fast toward stiff/safe gains, slow toward soft), |u0| <= 20 mm, |u1| <= 0.2 m/s.

3. Update rules
   ORD (final design): warm start = previous decision shifted along its feedforward ramp; projected
   Gauss-Newton / Levenberg-Marquardt: forward-difference Jacobian of r (4 rollout sets), box QP
   min 1/2 d'(J'J + mu diag)d + (J'r)'d on the shifted box (projected Newton active set + Gauss-Seidel
   polish), cost acceptance (1 rollout set); 1 iteration per update = 6 rollout sets.
   NU (neighbouring / quasi-Newton update of the earlier candidate, ablation c), per update:
     dF = F(a_k; c_{k+1}) - F(a_k; c_k) ~ Jc'(r(a_k; c_{k+1}) - r(a_k; c_k))   (one rollout set at the new
          context; Jc cached Jacobian)
     a_pred = Pi(a_k - B dF), componentwise step limits, Pi = box projection;
     cost rejection (keep a_k and re-anchor next update with an ordinary GN iteration + new Jc, B);
     on acceptance BFGS correction of the cached inverse curvature B with s = a_pred - a_k,
     y = Jc'(r(a_pred) - r(a_k)) (positivity check s'y > 1e-8 |s||y|); then ncorr projected quasi-Newton
     corrections a <- Pi(a - B Jc' r(a)) with cost acceptance and BFGS; Jc and B = (Jc'Jc + mu diag)^-1
     refreshed every nref updates.  ncorr = 1, nref = 2 gives ~5.5 rollout sets per update, matched to ORD.
   Compute per update (estimator + OCP + Python glue) is measured with time.perf_counter inside the policy.

Departures from the earlier candidate design (research/separate_candidate_design_and_iteration.md Secs. 3, 10;
experiments/dual_candidate_task_simulation):
 - No HJB critic / trajectory-local iLQR value: a direct finite-horizon OCP with 4 decision variables (one
   gain block, ramp-then-hold feedforward) instead of 10-dimensional two-segment actions with gain-transition
   dynamics and an energy state; no energy tank / stored-energy term (not part of this benchmark's interface).
 - 1-D normal axis with this benchmark's admittance + KPF + u_ff structure instead of the B0 2-D force/moment
   interface with four gains.
 - Estimator: local log-stiffness + quadratic equivalent-surface trajectory + force level, one-step filter
   inversion with the reconstructed 1 kHz robot position (equation error), instead of the two-region Hertz
   output-error fit with moments; acceleration term added; bias absorbed in g0; no SLS/viscous term.
 - Directional residual with a cached Jacobian (one rollout set, not a fresh gradient); BFGS secant pairs use
   the same cached-Jacobian projection; re-anchoring by an ordinary iteration after any rejection.
 - Equilibrium robustness probes added to the cost (not in the earlier design): without them the OCP's gain
   freedom raised train J (see candidate1_log.md).
 - Ramp-then-hold feedforward (a horizon-long ramp produced a slow feedforward/integrator null-mode
   oscillation); rule-form gain prior; asymmetric rate limits; conservative release and dither from EtR-A.
 - Final solver is the ordinary warm-started solve: the neighbouring update did not help on train (it is
   kept and reported as ablation c).

Post-test findings (no test re-run; the reported test numbers are for the design above, k0 = 5000 N/m):
 - Start-up prior. On homogeneous 1 MPa tissue (known-behaviour check) the soft start prior k0 = 5000 N/m gives
   soft gains at t = 0 on a 26-31 kN/m contact; one of three scenes oscillates, loses contact and saturates u_ff
   (J 33). v2: k0 = 3e4 N/m (start on stiff, safe gains). It is J-identical on T3/T1 train and dev (all those
   episodes start on soft tissue) and removes the failure (stiff N1 J 0.014-0.025). Recommended setting.
 - Fast push (T3, 0.25 s ramp): retractions of up to 44 mm/s outrun the robot; after contact loss the clipped
   model force gives the OCP no gradient and the drift estimate degrades, and 2/40 test scenes diverge.
   The flag lossfix = 1 (loss barrier on the unclipped force, predictor-side u_ff clip, drift reset when the
   window is mostly in loss) made dev worse (J 92 vs 5.8) and is not used (default 0). Unresolved.

Run (<= 2 processes): python -B candidate1_model.py tune|test|solver|known [cells] [variants]
"""
import math, time, json, sys, itertools, copy
import numpy as np
from numba import njit
import plant as P
import scenarios as S
import tune as U

FD = P.FD
LDMIN, LDMAX = math.log(P.GMIN[1]), math.log(P.GMAX[1])
NX = 4                                         # decision x = [log D, KPF/1e-4, u0/1e-3 m, u1/1e-2 m/s]
XLO = np.array([LDMIN, 0.0, -20.0, -20.0]); XHI = np.array([LDMAX, 10.0, 20.0, 20.0])
FDH = np.array([2e-3, 2e-3, 2e-3, 2e-3])       # finite-difference steps (scaled units)


# ============================================================================ numba kernels
@njit(cache=False)
def _sync(log, m0, n, rs, M, D, K, ki, Tc, VMAX):
    """replicate the plant admittance (yc, v, I) over logged ticks m0..n-1 with the applied gains."""
    yc = rs[0]; v = rs[1]; I = rs[2]; yref = rs[3]
    e = math.exp(-D * Tc / M)
    for j in range(m0, n):
        ef = FD - log[j, 5]
        I = min(max(I + ef * Tc, -20.), 20.)
        u = ef + ki * I - K * (yc - yref)
        vn = v * e + (1 - e) * u / D
        vn = min(max(vn, -VMAX), VMAX)
        yc += 0.5 * (v + vn) * Tc; v = vn
    rs[0] = yc; rs[1] = v; rs[2] = I


@njit(cache=False)
def _est(log, nend, Nw, lam, ndel, eb, afc, ns, th, mu, pw, nfree, lo, hi, iters, sz, info):
    """moving-window constrained Gauss-Newton fit (MHE with arrival cost) of
         f(t) = g0 + k (y(t) - y_e) + k (vs tau + as tau^2 / 2),   tau = t - t_e,  k = exp(kappa)
    to the one-step filter-inverted force data z_j = F_{j+1} - a^ns F_j = sum_m w_m f(t_j + m ms),
    with y(t) reconstructed at 1 kHz from the encoder sample and the known servo target/lag.
    th = [kappa, g0, vs, as] (in/out); mu, pw: arrival prior mean and weight; nfree: 2 (no drift),
    3 (+vs), 4 (+as). Window: intervals j -> j+1 with rows j+1 <= nend-1."""
    j0 = max(ndel + 1, nend - 1 - Nw); j1 = nend - 2
    m = j1 - j0 + 1
    info[0] = 0.; info[1] = 0.
    if m < 3: return 0
    te = log[nend - 1, 0]; ye = log[nend - 1, 2]
    a1 = 1.0 - afc; A8 = a1 ** ns; Ws = 1.0 - A8
    Z = np.empty(m); Y = np.empty(m); T1 = np.empty(m); T2 = np.empty(m); Wt = np.empty(m)
    nv = 0
    for jj in range(m):
        j = j0 + jj
        Fj = log[j, 5]; Fj1 = log[j + 1, 5]
        Z[jj] = Fj1 - A8 * Fj
        tgt = log[j - ndel, 11]; yy = log[j, 2]; tj = log[j, 0]
        sY = 0.; sT = 0.; sT2 = 0.
        for mm in range(1, ns + 1):
            yy = tgt + (yy - tgt) * eb
            w = afc * a1 ** (ns - mm)
            tau = tj + mm * 1e-3 - te
            sY += w * (yy - ye); sT += w * tau; sT2 += w * 0.5 * tau * tau
        Y[jj] = sY; T1[jj] = sT; T2[jj] = sT2
        ok = (Fj > 0.5) and (Fj1 > 0.5)
        if ok:
            Wt[jj] = lam ** (j1 - j); nv += 1
        else:
            Wt[jj] = 0.
    if nv < 3: return 0
    A = np.zeros((4, 4)); b = np.zeros(4); Jr = np.zeros(4)
    ssr = 0.
    for it in range(iters):
        k = math.exp(th[0])
        A[:, :] = 0.; b[:] = 0.; ssr = 0.
        for jj in range(m):
            if Wt[jj] == 0.: continue
            pred = Ws * th[1] + k * (Y[jj] + th[2] * T1[jj] + th[3] * T2[jj])
            r = (Z[jj] - pred) / sz
            Jr[0] = k * (Y[jj] + th[2] * T1[jj] + th[3] * T2[jj]) / sz
            Jr[1] = Ws / sz; Jr[2] = k * T1[jj] / sz; Jr[3] = k * T2[jj] / sz
            w = Wt[jj]; ssr += w * r * r
            for p in range(4):
                b[p] += w * Jr[p] * r
                for q in range(4): A[p, q] += w * Jr[p] * Jr[q]
        for p in range(4):
            A[p, p] += pw[p] + 1e-9; b[p] += pw[p] * (mu[p] - th[p])
        for p in range(nfree, 4):
            for q in range(4):
                A[p, q] = 0.; A[q, p] = 0.
            A[p, p] = 1.; b[p] = 0.
        d = np.linalg.solve(A, b)
        for p in range(nfree):
            th[p] = min(max(th[p] + d[p], lo[p]), hi[p])
    info[0] = nv; info[1] = math.sqrt(ssr / max(nv, 1))
    return nv


@njit(cache=False)
def _roll(x, H, per, st, q, dith, pr, kmult, dh, wsc, out, o0, mode, vd):
    """closed-loop predictor of the exact robot side (admittance ZOH, KPF path, u_ff, ndel delay,
    first-order servo lag at 1 kHz, 1 kHz first-order F/T filter) through the estimated contact
    model with the surface motion extrapolated. Starts at the last logged tick t_{n-1}; interval 0
    is already committed. Writes 3 residuals per predicted tick t_{n+1}..t_{n+H} and one terminal term.
    mode 0: from the current (replicated) state with the estimated drift; dh = surface step [m].
    mode 1: robustness probe from equilibrium (F = FD, at rest, no feedforward/dither) against an
            unpredicted surface ramp of velocity vd on a contact of stiffness kmult * k_hat."""
    M = pr[0]; K = pr[1]; ki = pr[2]; Tc = pr[3]; VMAX = pr[4]; eb = pr[5]; afc = pr[6]; ns = int(pr[7])
    ga = pr[8]; tv = pr[9]; Fc = pr[10]; wc = math.sqrt(pr[11]); Fl = pr[12]; wl = math.sqrt(pr[13])
    D = math.exp(x[0]); KP = x[1] * 1e-4; u0 = x[2] * 1e-3; u1 = x[3] * 1e-2
    yc = st[0]; v = st[1]; I = st[2]; yref = st[3]; y = st[4]; Ff = st[5]; tgt = st[6]
    g0 = st[7]; k = st[8] * kmult; vs = st[9]; acc = st[10]; ye = st[11]
    nd = q.shape[0]
    buf = np.empty(max(nd, 1))
    for i in range(nd): buf[i] = q[i]
    if mode == 1:
        yc = ye; v = 0.; I = 0.; yref = ye; y = ye; Ff = FD; tgt = ye; g0 = FD; acc = 0.; vs = vd; u0 = 0.; u1 = 0.
        tv = 1e9
        for i in range(nd): buf[i] = ye
    e = math.exp(-D * Tc / M)
    sH = wsc / math.sqrt(H)
    m = o0; f = 0.
    for i in range(H + 1):
        for mm in range(1, ns + 1):
            y = tgt + (y - tgt) * eb
            tau = (i * ns + mm) * 1e-3
            tt = min(tau, tv)
            hh = vs * tt + 0.5 * ga * acc * tt * tt + dh
            f = g0 + k * (y - ye + hh)
            flin = f
            if f < 0.: f = 0.
            Ff += afc * (f - Ff)
        if i >= 1:
            out[m] = sH * (f - FD) / 0.5
            out[m + 1] = sH * wc * (f - Fc if f > Fc else 0.)
            fl_ = flin if pr[16] > 0 else f      # lossfix: barrier on the unclipped model force (gradient in loss)
            out[m + 2] = sH * wl * ((Fl - fl_) / 0.5 if fl_ < Fl else 0.)
            m += 3
        if i == H: break
        ef = FD - Ff
        I = min(max(I + ef * Tc, -20.), 20.)
        uu = ef + ki * I - K * (yc - yref)
        vn = v * e + (1 - e) * uu / D
        vn = min(max(vn, -VMAX), VMAX)
        yc += 0.5 * (v + vn) * Tc; v = vn
        blk = (i // per) * per
        uf = u0 + u1 * min(blk * Tc, pr[15]) + (dith[i] if mode == 0 else 0.)
        if pr[16] > 0: uf = min(max(uf, -0.02), 0.02)     # lossfix: plant clip of u_ff known to the predictor
        ycmd = yc + KP * ef + uf
        if nd > 0:
            tgt = buf[0]
            for kk in range(nd - 1): buf[kk] = buf[kk + 1]
            buf[nd - 1] = ycmd
        else:
            tgt = ycmd
    # terminal velocity-consistency term: the admittance (not the feedforward offset) should carry the
    # extrapolated surface velocity at the horizon end (prevents integrator / feedforward cancellation)
    tend = (H + 1) * ns * 1e-3
    vse = (vs + ga * acc * min(tend, tv)) if tend < tv else 0.0
    out[m] = wsc * math.sqrt(pr[14]) * (v + vse) / 0.01
    return m + 1


@njit(cache=False)
def _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, out):
    """stacked residual vector: per scenario (kmult, surface step dh, weight) 3H force/barrier terms,
    then 7 regularisation terms (gain/feedforward rate to the previously applied decision, gain prior,
    feedforward prior)."""
    m = 0
    for s in range(scen.shape[0]):
        m = _roll(x, H, per, st, q, dith, pr, scen[s, 0], scen[s, 1], math.sqrt(scen[s, 2]), out, m,
                  int(scen[s, 3]), scen[s, 4])
    out[m] = rw[0] * (x[0] - prev[0]); out[m + 1] = rw[1] * (x[1] - prev[1]); out[m + 2] = rw[2] * (x[2] - prev[2])
    out[m + 3] = rw[3] * (x[0] - prior[0]); out[m + 4] = rw[4] * (x[1] - prior[1])
    out[m + 5] = rw[5] * x[2]; out[m + 6] = rw[6] * x[3]
    return m + 7


@njit(cache=False)
def _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, J, r2):
    nr = r.shape[0]
    for j in range(NX):
        h = FDH[j]
        if x[j] + h > hi[j]: h = -h
        xj = x[j]; x[j] = xj + h
        _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r2)
        x[j] = xj
        for i in range(nr): J[i, j] = (r2[i] - r[i]) / h


@njit(cache=False)
def _boxqp(Hm, g, lo, hi):
    """min 0.5 d'Hd + g'd on lo <= d <= hi (convex, small n): projected-Newton active set, then
    projected Gauss-Seidel polish."""
    n = g.shape[0]
    d = np.zeros(n)
    for i in range(n): d[i] = min(max(0., lo[i]), hi[i])
    free = np.ones(n, dtype=np.bool_)
    for it in range(8):
        gr = Hm @ d + g
        nf = 0
        for i in range(n):
            free[i] = not ((d[i] <= lo[i] + 1e-12 and gr[i] > 0) or (d[i] >= hi[i] - 1e-12 and gr[i] < 0))
            if free[i]: nf += 1
        if nf == 0: break
        idx = np.empty(nf, dtype=np.int64); t = 0
        for i in range(n):
            if free[i]: idx[t] = i; t += 1
        Hf = np.empty((nf, nf)); rf = np.empty(nf)
        for a in range(nf):
            rf[a] = -gr[idx[a]]
            for c in range(nf): Hf[a, c] = Hm[idx[a], idx[c]]
        step = np.linalg.solve(Hf, rf)
        moved = 0.
        for a in range(nf):
            i = idx[a]; dn = min(max(d[i] + step[a], lo[i]), hi[i]); moved += abs(dn - d[i]); d[i] = dn
        if moved < 1e-12: break
    for sweep in range(20):
        for i in range(n):
            gi = g[i]
            for c in range(n): gi += Hm[i, c] * d[c]
            if Hm[i, i] > 0:
                d[i] = min(max(d[i] - gi / Hm[i, i], lo[i]), hi[i])
    return d


@njit(cache=False)
def _cost(r):
    return 0.5 * np.sum(r * r)


@njit(cache=False)
def _solve_ord(x, iters, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, Jc, rc, info):
    """ordinary re-solve: warm-started projected Gauss-Newton / Levenberg-Marquardt with cost
    acceptance. Returns final cost; info[0] = rollouts used. Leaves J and r at the final point in Jc, rc
    (the last Jacobian evaluated, i.e. at the start point of the last iteration)."""
    nr = rc.shape[0]; nsc = scen.shape[0]
    r = np.empty(nr); r2 = np.empty(nr); rn = np.empty(nr); J = np.empty((nr, NX))
    _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r); nroll = nsc
    c = _cost(r); lamb = 1e-3
    for it in range(iters):
        _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, J, r2); nroll += NX * nsc
        Hm = J.T @ J; g = J.T @ r
        for i in range(NX): Hm[i, i] += lamb * (Hm[i, i] + 1e-6)
        d = _boxqp(Hm, g, lo - x, hi - x)
        xn = x + d
        _resid(xn, H, per, st, q, dith, pr, scen, prev, prior, rw, rn); nroll += nsc
        cn = _cost(rn)
        Jc[:, :] = J
        if cn < c:
            x[:] = xn; r[:] = rn; c = cn; lamb = max(lamb * 0.3, 1e-6)
        else:
            lamb *= 10.
    rc[:] = r
    info[0] = nroll
    return c


@njit(cache=False)
def _bfgs(B, s, yv):
    sy = 0.; ss = 0.; yy = 0.
    for i in range(NX): sy += s[i] * yv[i]; ss += s[i] * s[i]; yy += yv[i] * yv[i]
    if sy <= 1e-8 * math.sqrt(ss * yy) or sy <= 1e-14: return False    # positivity check
    rho = 1.0 / sy
    I = np.eye(NX)
    A = I - rho * np.outer(s, yv)
    B[:, :] = A @ B @ A.T + rho * np.outer(s, s)
    return True


@njit(cache=False)
def _solve_nu(x, ncorr, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, Jc, B, rprev, smax, info):
    """neighbouring / quasi-Newton update (earlier candidate, Sec. 3.2):
         dF = F(a_k; c_{k+1}) - F(a_k; c_k),  F = grad_a J ~ Jc' r   (directional residual: one rollout
              of the new context at the old decision, cached Jacobian Jc)
         a_pred = Pi(a_k - B dF)            (B cached inverse curvature, step limit smax, box Pi)
       cost rejection, BFGS correction of B on accepted steps (positivity check), then ncorr projected
       quasi-Newton corrections with the fresh residual. info: [rollouts, accepted_pred, n_bfgs, rejected].
       On return rprev holds the residual at the returned decision (the anchor for the next update)."""
    nr = rprev.shape[0]; nsc = scen.shape[0]
    r1 = np.empty(nr); rp = np.empty(nr)
    _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r1); nroll = nsc
    c1 = _cost(r1)
    dF = Jc.T @ (r1 - rprev)
    step = -(B @ dF)
    xp = np.empty(NX)
    for i in range(NX):
        s_ = min(max(step[i], -smax[i]), smax[i])
        xp[i] = min(max(x[i] + s_, lo[i]), hi[i])
    _resid(xp, H, per, st, q, dith, pr, scen, prev, prior, rw, rp); nroll += nsc
    cp = _cost(rp)
    nb = 0; acc = 0; rej = 0
    if cp < c1:
        acc = 1
        s = xp - x; yv = Jc.T @ (rp - r1)
        if _bfgs(B, s, yv): nb += 1
        x[:] = xp; cur = rp.copy(); cc = cp
    else:
        rej = 1; cur = r1.copy(); cc = c1
    # projected quasi-Newton corrections
    for it in range(ncorr):
        g = Jc.T @ cur
        step = -(B @ g)
        for i in range(NX):
            s_ = min(max(step[i], -smax[i]), smax[i])
            xp[i] = min(max(x[i] + s_, lo[i]), hi[i])
        _resid(xp, H, per, st, q, dith, pr, scen, prev, prior, rw, rp); nroll += nsc
        cp = _cost(rp)
        if cp < cc:
            s = xp - x; yv = Jc.T @ (rp - cur)
            if _bfgs(B, s, yv): nb += 1
            x[:] = xp; cur[:] = rp; cc = cp
        else:
            break
    rprev[:] = cur
    info[0] = nroll; info[1] = acc; info[2] = nb; info[3] = rej
    return cc


@njit(cache=False)
def _refresh(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, Jc, B, lamb):
    """new cached Jacobian at x and inverse curvature B = (J'J + lamb diag)^-1."""
    r2 = np.empty(r.shape[0])
    _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, Jc, r2)
    Hm = Jc.T @ Jc
    for i in range(NX): Hm[i, i] += lamb * (Hm[i, i] + 1e-6)
    B[:, :] = np.linalg.inv(Hm)
    return NX * scen.shape[0]


# ============================================================================ policy
DEFAULTS = dict(
    # admittance parts held fixed (tuned)
    M=0.3, K=0., ki=30.,
    # gain prior rule (B3b form evaluated at the estimate): D_r = max(Dmin, k/wI), KPF_r = min(cP/k, KPFmax)
    wI=80., cP=2., Dmin=20.,
    # estimator
    Nw=16, lam=0.85, sz=0.03, sk=0.15, sv=0.01, sa=0.5, k0=5000., iters_est=3, drift=2, trel=0.1,
    klat=0.0,                                  # extra latency (s) of the stiffness information used by the OCP
    # dither (excitation through u_ff; inside the plant and the force metric)
    Fdith=0.1, amin=2e-6, amax=60e-6, Ncyc=6,
    # OCP
    H=16, ga=0.0, tv=0.1, tu=0.05, Fc=6.5, wc=10., Fl=1.0, wl=20., wv=0.0,
    rD=0.3, rK=0.3, ru=0.0, pD=1.0, pK=1.0, pu=0.1, pu1=0.01,
    rate_fast=2.0, rate_slow=0.15, rateK_fast=10., rateK_slow=1.0, ffmode=1,
    probe=10.0, probe_k=2.0, probe_v=0.01,
    # solver
    solver='nu', iters=1, ncorr=1, nref=2, lamb=1e-3,
    update=0.008,
    lossfix=0,     # post-test fix (v2, not run on test): see docstring 'Post-test findings'
)


class Candidate1:
    oracle = False

    def __init__(self, **kw):
        bad = set(kw) - set(DEFAULTS)
        if bad: raise ValueError(bad)
        self.hp = dict(DEFAULTS); self.hp.update(kw)
        self.__dict__.update(self.hp)
        self.diag = False; self.rec = None

    # ---------------------------------------------------------------- setup
    def reset(self, pl):
        self.Tc = pl['Tc']; self.ndel = int(pl['ndel']); self.tau = pl['tau_a']; self.VMAX = pl['VMAX']
        self.eb = math.exp(-1e-3 / self.tau) if self.tau > 0 else 0.0
        self.afc = 1.0 - math.exp(-2 * math.pi * pl['fc'] * 1e-3)
        self.ns = int(round(self.Tc / 1e-3))
        self.per = max(int(round(self.update / self.Tc)), 1)
        self.rs = None; self.m = 0; self.g = None
        self.th = np.array([math.log(self.k0), FD, 0., 0.])
        self.kuse = self.k0; self.khist_raw = []
        self.est_lo = np.array([math.log(200.), -5., -0.2, -5.]); self.est_hi = np.array([math.log(3e5), 30., 0.2, 5.])
        self.nfree = 2 if self.drift == 0 else (3 if self.drift == 1 else 4)
        self.pr = np.array([self.M, self.K, self.ki, self.Tc, self.VMAX, self.eb, self.afc, self.ns, self.ga, self.tv,
                            self.Fc, self.wc, self.Fl, self.wl, self.wv, self.tu, self.lossfix])
        sc = [[1.0, 0.0, 1.0, 0, 0.]]
        if self.probe > 0:      # equilibrium robustness probes: unpredicted ramp vd at k_hat and probe_k * k_hat
            sc.append([1.0, 0.0, self.probe, 1, self.probe_v])
            if self.probe_k > 0: sc.append([self.probe_k, 0.0, self.probe, 1, self.probe_v])
        self.scen = np.array(sc)
        self.nr = (3 * self.H + 1) * len(sc) + 7
        self.rw = np.sqrt(np.array([self.rD, self.rK, self.ru, self.pD, self.pK, self.pu, self.pu1]))
        self.x = None
        self.Jc = np.zeros((self.nr, NX)); self.B = np.eye(NX); self.rprev = np.zeros(self.nr)
        self.age = 10 ** 9; self.need_ref = True
        self.smax = np.array([0.5, 3.0, 2.0, 5.0])
        self.dith_a = 0.
        # bookkeeping
        self.t_est = []; self.t_ocp = []; self.t_call = []; self.nroll = []; self.hist = []
        self.n_acc = 0; self.n_rej = 0; self.n_bfgs = 0; self.n_ref = 0; self.n_upd = 0; self.maxrep = 0.
        self.kq = []                                 # queue for added stiffness-information latency
        self.diag_gap = []

    def rule(self, k):
        D = max(self.Dmin, min(k / self.wI, P.GMAX[1])); KP = min(self.cP / k, P.GMAX[4])
        return D, KP

    def _dither(self, n):
        if self.Fdith <= 0: return 0.
        return self.dith_a * math.sin(2 * math.pi * (n % self.Ncyc) / self.Ncyc)

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
        if self.rs is None:     # plant init: yc = y_ref = encoder y(0), v = I = 0
            self.rs = np.array([log[0, 2], 0., 0., log[0, 2]]); self.m = 0
        M_, D_, K_, ki_, KP_, u_ = self.g
        _sync(log, self.m, n, self.rs, M_, D_, K_, ki_, self.Tc, self.VMAX)
        ef = FD - log[n - 1, 5]
        self.maxrep = max(self.maxrep, abs(self.rs[0] + KP_ * ef + u_ - log[n - 1, 11]))
        self.m = n
        # ---------------- estimator
        t1 = time.perf_counter()
        mu = self.th.copy(); dtu = self.per * self.Tc
        mu[2] = self.th[2] + (self.th[3] * dtu if self.nfree == 4 else 0.)
        pw = np.array([1.0 / (self.sk ** 2 * self.per), 0.0, 1.0 / (self.sv ** 2 * self.per), 1.0 / (self.sa ** 2 * self.per)])
        info = np.zeros(2)
        nv = _est(log, n, self.Nw, self.lam, self.ndel, self.eb, self.afc, self.ns, self.th, mu, pw, self.nfree,
                  self.est_lo, self.est_hi, self.iters_est, self.sz, info)
        if self.lossfix and nv < self.Nw // 2:      # v2: window mostly in contact loss -> no drift extrapolation
            self.th[2] = 0.; self.th[3] = 0.
        if nv == 0:          # no valid data (contact loss): keep stiffness, zero drift, current force
            self.th[0] = mu[0]; self.th[1] = log[n - 1, 5]; self.th[2] = 0.; self.th[3] = 0.
        kraw = math.exp(self.th[0])
        if self.klat > 0:
            self.kq.append((o.t[n - 1], kraw))
            while len(self.kq) > 1 and self.kq[1][0] <= o.t[n - 1] - self.klat + 1e-9: self.kq.pop(0)
            kraw = self.kq[0][1] if self.kq[0][0] <= o.t[n - 1] - self.klat + 1e-9 else self.k0
        if self.trel > 0:
            kraw = max(kraw, self.kuse * math.exp(-dtu / self.trel))
        self.kuse = kraw
        t2 = time.perf_counter()
        # ---------------- OCP context
        nd = self.ndel
        tgt0 = log[n - 1 - nd, 11] if n - 1 - nd >= 0 else log[0, 2]
        q = np.array([log[j, 11] if j >= 0 else log[0, 2] for j in range(n - nd, n)], dtype=np.float64)
        kp = self.kuse
        st = np.array([self.rs[0], self.rs[1], self.rs[2], self.rs[3], log[n - 1, 2], log[n - 1, 5], tgt0,
                       self.th[1], kp, getattr(self, '_vs', self.th[2]), getattr(self, '_as', self.th[3]),   # _vs/_as: oracle diagnostic hook only
                       log[n - 1, 2], 0., 0., 0., 0.])
        self.dith_a = min(max(self.Fdith / kp, self.amin), self.amax) if self.Fdith > 0 else 0.
        dith = np.array([self._dither(n + (i // self.per) * self.per) for i in range(self.H)])
        Dr, KPr = self.rule(kp)
        prior = np.array([math.log(Dr), KPr / 1e-4, 0., 0.])
        # previous applied decision (u0 shifted along its ramp)
        xprev = self.x.copy(); xprev[2] = xprev[2] + xprev[3] * 1e-2 * min(dtu, self.tu) / 1e-3
        xprev[2] = min(max(xprev[2], XLO[2]), XHI[2])
        lo = XLO.copy(); hi = XHI.copy()
        # rate limits per update on the gains (fast toward stiff / safe, slow toward soft)
        lo[0] = max(LDMIN, xprev[0] - self.rate_slow * self.per); hi[0] = min(LDMAX, xprev[0] + self.rate_fast * self.per)
        lo[1] = max(0., xprev[1] - self.rateK_fast * self.per); hi[1] = min(10., xprev[1] + self.rateK_slow * self.per)
        if self.ffmode == 0:        # gains only: u_ff (OCP part) fixed at 0 (drift still used in the prediction)
            lo[2] = hi[2] = 0.; lo[3] = hi[3] = 0.; xprev[2] = 0.; xprev[3] = 0.
        x = np.minimum(np.maximum(xprev, lo), hi)
        nrinfo = np.zeros(4)
        if self.solver == 'nu':
            if self.need_ref or self.age >= self.nref:
                r0 = np.empty(self.nr)
                _resid(x, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw, r0)
                nroll = self.scen.shape[0]
                if self.x is not None and self.need_ref:
                    # first / after rejection: ordinary solve to re-anchor
                    c = _solve_ord(x, max(self.iters, 1), self.H, self.per, st, q, dith, self.pr, self.scen, xprev,
                                   prior, self.rw, lo, hi, self.Jc, self.rprev, nrinfo)
                    nroll = nrinfo[0]
                    rr = self.rprev.copy()
                    nroll += _refresh(x, rr, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                                      lo, hi, self.Jc, self.B, self.lamb)
                    self.need_ref = False
                else:
                    # periodic refresh of the cached Jacobian / curvature at the warm start, then the update
                    nroll += _refresh(x, r0, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                                      lo, hi, self.Jc, self.B, self.lamb)
                    c = _solve_nu(x, self.ncorr, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                                  lo, hi, self.Jc, self.B, self.rprev, self.smax, nrinfo)
                    nroll += nrinfo[0]
                    self.n_acc += nrinfo[1]; self.n_bfgs += nrinfo[2]; self.n_rej += nrinfo[3]
                    if nrinfo[3]: self.need_ref = True
                self.age = 0; self.n_ref += 1
            else:
                c = _solve_nu(x, self.ncorr, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                              lo, hi, self.Jc, self.B, self.rprev, self.smax, nrinfo)
                nroll = nrinfo[0]
                self.n_acc += nrinfo[1]; self.n_bfgs += nrinfo[2]; self.n_rej += nrinfo[3]
                if nrinfo[3]: self.need_ref = True
                self.age += 1
        else:
            c = _solve_ord(x, self.iters, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw,
                           lo, hi, self.Jc, self.rprev, nrinfo)
            nroll = nrinfo[0]
        t3 = time.perf_counter()
        if self.diag:
            xs = x.copy(); Jt = np.zeros_like(self.Jc); rt = np.zeros(self.nr); it_ = np.zeros(4)
            xw = np.minimum(np.maximum(xprev, lo), hi)
            cs = _solve_ord(xw, 25, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw, lo, hi, Jt, rt, it_)
            cs2 = _solve_ord(xs, 25, self.H, self.per, st, q, dith, self.pr, self.scen, xprev, prior, self.rw, lo, hi, Jt, rt, it_)
            self.diag_gap.append((c, min(cs, cs2)))
        self.x = x
        if self.rec is not None: self.rec.append((n, st.copy(), q.copy(), dith.copy(), x.copy(), self.scen.copy()))
        self.n_upd += 1
        D = math.exp(x[0]); KP = x[1] * 1e-4; uff = x[2] * 1e-3
        self.g = (self.M, D, self.K, self.ki, KP, uff + dith[0])
        tend = time.perf_counter()
        self.t_est.append(t2 - t1); self.t_ocp.append(t3 - t2); self.t_call.append(tend - t0); self.nroll.append(nroll)
        if n % 4 == 0:
            self.hist.append((o.t[n - 1], math.exp(self.th[0]), kp, self.th[2], D, KP, uff, c))
        return self.g


# ============================================================================ evaluation (<= 2 processes)
NPROC = 2
_POOL = None
def warmup():
    """compile the numba kernels (so that measured compute per update excludes JIT time)."""
    p, T = S.make('N1', np.random.default_rng(0), dict(E=2e4, A=0.001))
    for kw in [dict(), dict(solver='ord'), dict(drift=0, probe=0.)]:
        P.episode(p, Candidate1(**kw), T=0.3, seed=0, update=0.008)


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
    plat = 'PC125'
    if kind in ('C1', 'C1diag'):
        pol = Candidate1(**kw); pol.diag = kind == 'C1diag'
        log = P.episode(p, pol, T=T, platform=plat, seed=seed, update=pol.update)
        m = S.metrics(log)
        tc = np.array(pol.t_call) if pol.t_call else np.zeros(1)
        m['t_call_mean'] = float(tc.mean()); m['t_call_p99'] = float(np.percentile(tc, 99)); m['t_call_max'] = float(tc.max())
        m['t_est_mean'] = float(np.mean(pol.t_est)) if pol.t_est else 0.; m['t_ocp_mean'] = float(np.mean(pol.t_ocp)) if pol.t_ocp else 0.
        m['t_ocp_p99'] = float(np.percentile(pol.t_ocp, 99)) if pol.t_ocp else 0.
        m['rollouts_mean'] = float(np.mean(pol.nroll)) if pol.nroll else 0.
        m['n_upd'] = pol.n_upd; m['nu_acc'] = pol.n_acc; m['nu_rej'] = pol.n_rej; m['nu_bfgs'] = pol.n_bfgs; m['n_ref'] = pol.n_ref
        m['replica_err'] = pol.maxrep
        if pol.diag_gap:
            gp = np.array(pol.diag_gap); m['gap_mean'] = float(np.mean(gp[:, 0] - gp[:, 1]))
            m['gap_rel_mean'] = float(np.mean((gp[:, 0] - gp[:, 1]) / np.maximum(gp[:, 1], 1e-9)))
            m['cost_mean'] = float(np.mean(gp[:, 0])); m['cost_star_mean'] = float(np.mean(gp[:, 1]))
        if keep:
            m['hist'] = np.array(pol.hist); m['log'] = log
    else:
        if kind == 'fixed': ctrl = tuple(kw)
        elif kind == 'B3b':
            from baselines_ref import TruthK; ctrl = TruthK(*kw)
        elif kind == 'EtR':
            from baselines_adaptive import EtR; ctrl = EtR(*kw[0], **kw[1])
        log = P.episode(p, ctrl, T=T, platform=plat, seed=seed)
        m = S.metrics(log)
        if keep: m['log'] = log
    return m


def run(scenes, specs, seed0, keep=False, procs=NPROC):
    jobs = [(p, T, s, seed0 + i, keep) for i, (p, T) in enumerate(scenes) for s in specs]
    if procs > 1 and len(jobs) > 1:
        res = pool().map(_job, jobs, chunksize=max(1, len(jobs) // (8 * procs)))
    else:
        res = [_job(j) for j in jobs]
    ns = len(specs)
    return [res[i * ns:(i + 1) * ns] for i in range(len(scenes))]


def summary(rows):
    out = U.summary(rows)
    for k in ['t_call_mean', 't_est_mean', 't_ocp_mean', 'rollouts_mean', 'replica_err', 'gap_mean', 'gap_rel_mean',
              'cost_mean', 'cost_star_mean']:
        if k in rows[0]: out[k] = float(np.mean([r[k] for r in rows]))
    for k in ['t_call_p99', 't_ocp_p99']:
        if k in rows[0]: out[k] = float(np.mean([r[k] for r in rows]))
    if 't_call_max' in rows[0]: out['t_call_max'] = float(np.max([r['t_call_max'] for r in rows]))
    for k in ['nu_acc', 'nu_rej', 'nu_bfgs', 'n_ref', 'n_upd']:
        if k in rows[0]: out[k] = float(np.sum([r[k] for r in rows]))
    return out


def select(res):
    Jm = np.array([[m['J'] for m in row] for row in res]).mean(0)
    pk = np.array([[m['peak'] if np.isfinite(m['peak']) else 1e9 for m in row] for row in res])
    feas = (pk <= FD + 3.0).mean(0) >= 0.95
    j = int(np.argmin(np.where(feas, Jm, Jm + 1e6)))
    return j, Jm, feas


# ============================================================================ protocol
TR0 = S.SPLIT['train']; DEV0 = S.SPLIT['dev']; TE0 = S.SPLIT['test']
NTR, NTE = 20, 40
CELLS = {'T3': ('T3', dict(ramp=0.5, A=0.001)), 'T1': ('T1', dict(ratio=30., A=0.001)),
         'T3fast': ('T3', dict(ramp=0.25, A=0.001))}
# declared search space (stage A: 32 random draws, stage B: 16 one-parameter perturbations of the
# stage-A best, plus the hand-set default = 49 configurations per method variant and cell)
SPACE = dict(
    ki=[10., 30.], M=[0.3, 1.0], wI=[40., 80., 160.], cP=[1., 2., 4.],
    pD=[0.1, 1., 10.], probe=[3., 10., 30.], probe_k=[1.5, 2., 3.], probe_v=[0.003, 0.01],
    pu=[0.03, 0.1, 0.3], tu=[0.03, 0.05, 0.08], ga=[0., 1.], tv=[0.05, 0.1, 0.2], H=[16, 24],
    lam=[0.8, 0.85, 0.9], sk=[0.08, 0.15, 0.3], Nw=[12, 16, 24], trel=[0., 0.1], Fdith=[0., 0.05, 0.1],
    sv=[0.005, 0.01, 0.02])
NA, NB = 32, 16
VARIANTS = {
    # final design: ordinary warm-started projected Gauss-Newton (1 iteration / update). Chosen on TRAIN before any
    # test run because the neighbouring update (as specified) had higher train J on T3 (0.210 vs 0.144) and equal on T1.
    'C1': dict(solver='ord', iters=1),
    'abl_a_no_surface_motion': dict(solver='ord', iters=1, drift=0),
    'abl_b_gains_only': dict(solver='ord', iters=1, ffmode=0),
    'abl_c_neighbouring_update': dict(solver='nu', ncorr=1, nref=2),   # ~same rollout budget per update
    'abl_d1_klat_0.1s': dict(solver='ord', iters=1, klat=0.1),
    'abl_d2_update_40ms': dict(solver='ord', iters=1, update=0.04),
}


def _cfg(d):
    d = dict(d)
    if 'pD' in d: d['pK'] = d['pD']
    if 'H' in d: d['H'] = int(d['H'])
    if 'Nw' in d: d['Nw'] = int(d['Nw'])
    return d


def stage_a(rng):
    return [_cfg({k: v[int(rng.integers(len(v)))] for k, v in SPACE.items()}) for _ in range(NA)]


def stage_b(best, rng):
    out = []; keys = list(SPACE)
    while len(out) < NB:
        k = keys[int(rng.integers(len(keys)))]
        alt = [v for v in SPACE[k] if v != best.get(k, DEFAULTS[k])]
        if not alt: continue
        c = dict(best); c[k] = alt[int(rng.integers(len(alt)))]; c = _cfg(c)
        if c not in out: out.append(c)
    return out


def tune_variant(cell, variant, log=print):
    task, dose = CELLS[cell]; fixed = VARIANTS[variant]
    train = S.cohort(task, NTR, TR0, dose)
    rng = np.random.default_rng(12345)                    # same draws for every variant (paired budget)
    A = [dict()] + stage_a(rng)
    t0 = time.time()
    rA = run(train, [('C1', dict(c, **fixed)) for c in A], TR0)
    jA, JmA, fA = select(rA)
    best = A[jA]
    B = stage_b(best, np.random.default_rng(777))
    rB = run(train, [('C1', dict(c, **fixed)) for c in B], TR0)
    Jall = np.concatenate([JmA, select(rB)[1]]); fall = np.concatenate([fA, select(rB)[2]])
    cand = A + B
    j = int(np.argmin(np.where(fall, Jall, Jall + 1e6)))
    cfg = dict(cand[j], **fixed)
    log(f"[tune] {cell} {variant}: best train J {Jall[j]:.3f} feasible {bool(fall[j])} ({len(cand)} configs, "
        f"{int(fall.sum())} feasible, {time.time() - t0:.0f}s) cfg {cand[j]}")
    return dict(config=cfg, searched=cand[j], train_J=float(Jall[j]), train_feasible=bool(fall[j]), n_configs=len(cand),
                n_feasible=int(fall.sum()), train_J_all=[float(x) for x in Jall], stageA_best=A[jA], seconds=time.time() - t0)


def gate_ref(task, dose):
    for c in json.load(open('gate_PC125.json')):
        if c['task'] == task and c['dose'] == dose: return c


def etr_spec(cell):
    """EtR-A as tuned in adaptive_results.json (only the two main cells)."""
    task, dose = CELLS[cell]
    try:
        for c in json.load(open('adaptive_results.json'))['cells']:
            if c['task'] == task and c['dose'] == dose:
                t = c['methods']['EtR-A']['tuning']; q = t['rule_params']
                return ('EtR', (tuple(q), dict(t['est_params'], rule='truthk')))
    except Exception:
        return None


def test_cell(cell, cfgs, log=print):
    """cfgs: {name: config}; recomputes B1, B3b(0.1 s, 0 s) and EtR-A on the same test seeds."""
    task, dose = CELLS[cell]; test = S.cohort(task, NTE, TE0, dose); g = gate_ref(task, dose)
    names = list(cfgs); specs = [('C1', cfgs[k]) for k in names]
    base = [('B1', ('fixed', tuple(g['B1']['gains']))), ('B3b_0.1', ('B3b', (0.1,) + tuple(g['B3b']['0.1']['rule']))),
            ('B3b_0.0', ('B3b', (0.0,) + tuple(g['B3b']['0.0']['rule'])))]
    e = etr_spec(cell)
    if e is not None: base.append(('EtR-A', e))
    res = run(test, specs + [b[1] for b in base], TE0)
    out = {'cell': cell, 'task': task, 'dose': dose, 'n_test': NTE, 'methods': {}}
    per = {}
    for j, nm in enumerate(names + [b[0] for b in base]):
        rows = [r[j] for r in res]; per[nm] = rows
        out['methods'][nm] = dict(test=summary(rows), J_per_scene=[float(r['J']) for r in rows],
                                  rmse_per_scene=[float(r['rmse']) for r in rows], peak_per_scene=[float(r['peak']) for r in rows])
    for nm in names:
        jj = [r['J'] for r in per[nm]]
        for ref in [b[0] for b in base]:
            out['methods'][nm]['dJ_vs_' + ref] = U.paired_ci(jj, [r['J'] for r in per[ref]])
        out['methods'][nm]['drmse_vs_B1'] = U.paired_ci([r['rmse'] for r in per[nm]], [r['rmse'] for r in per['B1']])
    for nm in names + [b[0] for b in base]:
        s_ = out['methods'][nm]['test']
        log(f"[test] {cell} {nm:26s} J {s_['J']:.3f} rmse {s_['rmse']:.3f} over {s_['frac_peak_over']:.2f} peak {s_['peak']:.2f} "
            f"inband {s_['inband']:.3f} loss {s_['loss']:.4f}" + (f"  dJ vs B1 {out['methods'][nm]['dJ_vs_B1']}" if nm in names else ''))
    return out


def solver_study(cell, cfg, nsc=6, log=print):
    """decision quality and compute of the neighbouring update vs ordinary re-solve at matched and larger
    budgets: same hyperparameters, solver swapped. Dev scenes; diag mode adds a converged (25 GN iterations)
    reference solve at every update (not applied) to measure the cost gap. Timing is measured in a separate
    single-process pass without diag."""
    task, dose = CELLS[cell]; dev = S.cohort(task, nsc, DEV0, dose)
    sol = {'NU(ncorr1,nref2)': dict(solver='nu', ncorr=1, nref=2), 'NU(ncorr2,nref4)': dict(solver='nu', ncorr=2, nref=4),
           'NU(ncorr0,nref1)': dict(solver='nu', ncorr=0, nref=1),
           'ORD(1 iter)': dict(solver='ord', iters=1), 'ORD(2 iter)': dict(solver='ord', iters=2),
           'ORD(5 iter)': dict(solver='ord', iters=5)}
    names = list(sol)
    rd = run(dev, [('C1diag', dict(cfg, **sol[k])) for k in names], DEV0)
    warmup()
    rt = run(dev[:3], [('C1', dict(cfg, **sol[k])) for k in names], DEV0, procs=1)
    out = {}
    for j, k in enumerate(names):
        a = [r[j] for r in rd]; b = [r[j] for r in rt]
        out[k] = dict(closed_loop_J=float(np.mean([m['J'] for m in a])), gap_mean=float(np.mean([m['gap_mean'] for m in a])),
                      gap_rel_mean=float(np.mean([m['gap_rel_mean'] for m in a])), cost_mean=float(np.mean([m['cost_mean'] for m in a])),
                      rollout_sets_per_update=float(np.mean([m['rollouts_mean'] for m in a])),
                      t_update_mean_us=1e6 * float(np.mean([m['t_call_mean'] for m in b])),
                      t_ocp_mean_us=1e6 * float(np.mean([m['t_ocp_mean'] for m in b])),
                      t_ocp_p99_us=1e6 * float(np.mean([m['t_ocp_p99'] for m in b])),
                      t_est_mean_us=1e6 * float(np.mean([m['t_est_mean'] for m in b])),
                      t_update_max_us=1e6 * float(np.max([m['t_call_max'] for m in b])),
                      nu_accept=float(np.sum([m['nu_acc'] for m in a])), nu_reject=float(np.sum([m['nu_rej'] for m in a])),
                      nu_bfgs=float(np.sum([m['nu_bfgs'] for m in a])), n_upd=float(np.sum([m['n_upd'] for m in a])))
        o = out[k]
        log(f"[solver] {cell} {k:18s} J {o['closed_loop_J']:.3f} gap {o['gap_mean']:.4f} (rel {o['gap_rel_mean']:.3f}) cost {o['cost_mean']:.3f} "
            f"rollout-sets {o['rollout_sets_per_update']:.1f} t_upd {o['t_update_mean_us']:.0f}us t_ocp {o['t_ocp_mean_us']:.0f}us (p99 {o['t_ocp_p99_us']:.0f}) "
            f"t_est {o['t_est_mean_us']:.0f}us acc {o['nu_accept']:.0f} rej {o['nu_reject']:.0f} bfgs {o['nu_bfgs']:.0f}/{o['n_upd']:.0f}")
    return out


def known_behaviour(cfg, b1, log=print):
    """homogeneous N1 tissue (breathing 1 mm), soft 15 kPa and stiff 1 MPa, 3 scenes each: estimated vs true
    stiffness (truth used only for scoring), chosen D / KPF, and J vs B1."""
    rng = np.random.default_rng(7); scenes = []
    for E in [1.5e4, 1.0e6]:
        for k in range(3):
            scenes.append(S.make('N1', rng, dict(E=E, A=0.001)))
    res = run(scenes, [('C1', dict(cfg)), ('fixed', tuple(b1))], 500, keep=True)
    out = {'soft_15kPa': [], 'stiff_1MPa': []}
    for (p, T), r in zip(scenes, res):
        m, mb = r[0], r[1]; h = m['hist']; lg = m['log']
        kt = lg[::4, 8][:len(h)]; ok = (h[:, 0] > 1.0) & (kt > 0)
        d = dict(E=float(p[0]), k_true_median=float(np.median(kt[ok])), k_hat_median=float(np.median(h[ok, 1])),
                 k_hat_over_true_median=float(np.median(h[ok, 1] / kt[ok])), k_hat_over_true_p10=float(np.percentile(h[ok, 1] / kt[ok], 10)),
                 k_hat_over_true_p90=float(np.percentile(h[ok, 1] / kt[ok], 90)),
                 D_median=float(np.median(h[ok, 4])), KPF_median=float(np.median(h[ok, 5])), uff_rms_mm=float(1e3 * np.sqrt(np.mean(h[ok, 6] ** 2))),
                 J=m['J'], rmse=m['rmse'], peak=m['peak'], J_B1=mb['J'], rmse_B1=mb['rmse'], peak_B1=mb['peak'])
        out['soft_15kPa' if p[0] < 1e5 else 'stiff_1MPa'].append(d)
    for k, L in out.items():
        log(f"[known] {k}: k_true {np.mean([d['k_true_median'] for d in L]):.0f} k_hat {np.mean([d['k_hat_median'] for d in L]):.0f} "
            f"ratio {np.mean([d['k_hat_over_true_median'] for d in L]):.2f} D {np.mean([d['D_median'] for d in L]):.0f} "
            f"KPF {np.mean([d['KPF_median'] for d in L]):.2e} J {np.mean([d['J'] for d in L]):.3f} vs B1 {np.mean([d['J_B1'] for d in L]):.3f} "
            f"rmse {np.mean([d['rmse'] for d in L]):.3f} vs {np.mean([d['rmse_B1'] for d in L]):.3f}")
    return out


def _save(R, path='candidate1_results.json'):
    json.dump(R, open(path, 'w'), indent=1, default=float)


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    cells = sys.argv[2].split(',') if len(sys.argv) > 2 else ['T3', 'T1']
    variants = sys.argv[3].split(',') if len(sys.argv) > 3 else list(VARIANTS)
    try:
        R = json.load(open('candidate1_results.json'))
    except Exception:
        R = {}
    R.setdefault('protocol', dict(n_train=NTR, n_test=NTE, train_seed=TR0, test_seed=TE0, dev_seed=DEV0, platform='PC125',
                                  search_space=SPACE, n_configs_per_variant=1 + NA + NB, selection='min mean train J s.t. peak <= 7 N in >= 95 % of train scenes',
                                  update_period_s=DEFAULTS['update'], defaults=DEFAULTS, variants=VARIANTS))
    R['protocol']['variants'] = VARIANTS
    if what == 'tune':
        for c in cells:
            for v in variants:
                R.setdefault('tuning', {}).setdefault(c, {})[v] = tune_variant(c, v); _save(R)
    elif what == 'test':
        for c in cells:
            cf = {v: R['tuning'][c][v]['config'] for v in variants if v in R['tuning'].get(c, {})}
            R.setdefault('test', {})[c] = test_cell(c, cf); _save(R)
    elif what == 'solver':
        for c in cells:
            R.setdefault('solver_study', {})[c] = solver_study(c, R['tuning'][c]['C1']['config']); _save(R)
    elif what == 'known':
        c = cells[0]; task, dose = CELLS[c]
        R.setdefault('known_behaviour', {})[c + '_tuned'] = known_behaviour(R['tuning'][c]['C1']['config'], gate_ref(task, dose)['B1']['gains']); _save(R)
    close_pool()
