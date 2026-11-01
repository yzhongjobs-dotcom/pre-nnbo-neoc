"""C1 GN kernels: same cost/solver, known current FF + direct-command rate cap.
Source: candidate1_model.py. mode1 robustness probes are unchanged.
A 16-state input is the exact original kernel; 17th state enables V4b input/cap.
"""
import math
import numpy as np
from numba import njit
from candidate1_model import FD,FDH,NX

@njit(cache=True)
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
    lastcmd = q[nd-1] if nd > 0 else st[6]
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
        if mode == 0 and st.shape[0] > 16:
            # Allocate the physical slew budget including the known current FF.
            # Future FF increments are zero. Return to the SAME residual frame
            # before propagating the servo: r=s+a_ff already includes FF motion.
            # Adding FF to both this servo and r would double-count it.
            physical = ycmd + st[16]
            physical = min(max(physical, lastcmd-VMAX*Tc), lastcmd+VMAX*Tc)
            lastcmd = physical
            ycmd = physical - st[16]
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

@njit(cache=True)
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

@njit(cache=True)
def _jac(x, r, H, per, st, q, dith, pr, scen, prev, prior, rw, lo, hi, J, r2):
    nr = r.shape[0]
    for j in range(NX):
        h = FDH[j]
        if x[j] + h > hi[j]: h = -h
        xj = x[j]; x[j] = xj + h
        _resid(x, H, per, st, q, dith, pr, scen, prev, prior, rw, r2)
        x[j] = xj
        for i in range(nr): J[i, j] = (r2[i] - r[i]) / h

@njit(cache=True)
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

@njit(cache=True)
def _cost(r):
    return 0.5 * np.sum(r * r)

@njit(cache=True)
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
