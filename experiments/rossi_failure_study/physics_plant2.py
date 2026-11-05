"""Rossi scan plant with the explicitly selected sliding-contact physics.

The public API deliberately mirrors :mod:`physics_plant`: ``episode``, ``replay``,
``_init`` and ``PLATFORM`` are kept compatible.  The additional state is the
contact-point tangential position/velocity and the servo-side scan position.
All contact calculations are performed at 2 kHz; force sensors are sampled at
1 kHz and logged at the controller period.
"""
import math
import os
import sys

import numpy as np
from numba import njit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cobot_scan_benchmark'))
import plant as OLD

R, L, X0, FD = OLD.R, OLD.L, OLD.X0, OLD.FD
NQ = OLD.NQ
halfwidth = OLD.halfwidth
PLATFORM = dict(OLD.PLATFORMS['PC125'])
LOGC = 20
NS = 23


@njit(cache=True)
def _lookup(g, dx, xi):
    u = xi / dx
    i = int(u)
    if i < 0:
        return g[0]
    if i >= len(g) - 1:
        return g[-1]
    f = u - i
    return g[i] * (1. - f) + g[i + 1] * f


@njit(cache=True)
def _fe(xc, d, ln_e, dx, law=0.):
    """The unchanged finite-footprint force used by physics_plant.py."""
    if d <= 0.:
        return 0.
    if law > 0.:
        return math.pi * L / 4. * math.exp(_lookup(ln_e, dx, xc - X0)) * .01 ** (1. - law) * d ** law
    a = halfwidth(d)
    out = 0.
    for i in range(1, NQ + 1):
        th = i * math.pi / (NQ + 1)
        xi = a * math.cos(th)
        w = math.pi / (NQ + 1) * math.sin(th) ** 2 * a * a
        out += L / (2. * R) * math.exp(_lookup(ln_e, dx, xc + xi - X0)) * w
    return out


@njit(cache=True)
def _fe_inner(xc, d, vc, ln_e, dx, k1, k2):
    """Finite-footprint pressure integral with the forward/backward radii.

    AuSoScan eqs. 11 and 14 use a half footprint on either side of the scan
    point.  Taking xi=a_i cos(theta) and the quadrature weight
    (pi/(NQ+1))*sin(theta)^2*a0*a_i gives exactly the original _fe weight when
    a1=a2=a0, hence v=0 is numerically the old integral.
    """
    if d <= 0.:
        return 0.
    a0 = halfwidth(d)
    a1 = max(a0 - k1 * abs(vc), .2 * a0)
    a2 = a0 + k2 * abs(vc)
    out = 0.
    for i in range(1, NQ + 1):
        th = i * math.pi / (NQ + 1)
        co = math.cos(th)
        ai = a2 if co*vc >= 0. else a1      # the larger radius lies ahead of the motion
        xi = ai * co
        w = math.pi / (NQ + 1) * math.sin(th) ** 2 * a0 * ai
        out += L / (2. * R) * math.exp(_lookup(ln_e, dx, xc + xi - X0)) * w
    return out


@njit(cache=True)
def _d_top():
    return OLD.psi(OLD.LREF*0.95)


@njit(cache=True)
def _fe_asym(xc, d, vc, ln_e, dx, k1, k2):
    """Past the probe edge the footprint stops growing while the force keeps rising along the tangent at the edge."""
    top = _d_top()
    if d <= top:
        return _fe_inner(xc, d, vc, ln_e, dx, k1, k2)
    f1 = _fe_inner(xc, top, vc, ln_e, dx, k1, k2); f0 = _fe_inner(xc, .98*top, vc, ln_e, dx, k1, k2)
    return f1+(d-top)*(f1-f0)/(.02*top)


@njit(cache=True)
def _friction(fn, vc, d, c):
    if fn <= 0.:
        return 0.
    mu_k, mu_s, v_s, v_e, k_d, nu, k1, k2 = c[13], c[14], c[15], c[16], c[17], c[18], c[19], c[20]
    a0 = halfwidth(max(d, 1e-12))
    sign = 1. if vc > 0. else (-1. if vc < 0. else 0.)
    mu = mu_k + (mu_s - mu_k) * math.exp(-(vc / max(v_s, 1e-9)) ** 2)
    static = mu * fn * math.tanh(vc / max(v_e, 1e-9))
    proj = 1. + (2. * nu - 1.) * 3. * a0 * a0 / (10. * R * R)
    a1 = max(a0 - k1 * abs(vc), .2 * a0)
    a2 = a0 + k2 * abs(vc)
    # AuSoScan eq. 16: integrate the asymmetric front/back pressure wedge
    # under a<<R, which gives 4 Fn (a2-a1) sign(v)/(3 pi R).
    fb = 4. * fn * (a2 - a1) * sign / (3. * math.pi * R)
    return static * proj + k_d * vc + fb


@njit(cache=True)
def _slope(h_s, dx, x):
    hd = max(dx, 1e-5)
    return (_lookup(h_s, dx, x + hd - X0) - _lookup(h_s, dx, x - hd - X0)) / (2. * hd)


@njit(cache=True)
def _xcmd_at(xcmd, t):
    u = max(0., t * 1000.)
    i = min(int(u), len(xcmd) - 2)
    f = u - i
    return xcmd[i] * (1. - f) + xcmd[i + 1] * f


@njit(cache=True)
def _tick(s, buf, cmd, c, ln_e, h_s, dx, xcmd, u, fext, act, rip, alpha, tau, noise, kk, log, nlog):
    """Advance one 8 ms control interval using 16 half-millisecond steps."""
    m, cb, kb, kr, lam, tau_a, afc, ns, ndel_f, r_c, r_low, r_inst, vmax = c[:13]
    dt = 5e-4
    ndel = int(ndel_f)
    Tc = .008
    dt_t = 2. * c[25] * math.sqrt(max(c[23] * c[24], 1e-12))
    ei = np.exp(-dt / tau)
    gi = np.exp(-dt / (2. * tau))
    prev = buf[ndel]
    for q in range(ndel):
        buf[q] = buf[q + 1]
    buf[ndel] = min(max(cmd, prev - vmax * Tc), prev + vmax * Tc)
    n_drive = len(u)
    for j in range(16):
        t = s[0]
        it = int(t * 1000. + 1e-9)
        i0 = min(it, n_drive - 2)
        fr = t * 1000. - i0
        uu = u[i0] * (1. - fr) + u[i0 + 1] * fr
        ud = (u[i0 + 1] - u[i0]) * 1000.
        fx = fext[i0]
        aa = act[i0]
        rp = rip[i0]

        xs_old = s[1]
        xs_cmd = _xcmd_at(xcmd, t)
        s[1] += dt / max(tau_a, 1e-9) * (xs_cmd - s[1])
        vs = (s[1] - xs_old) / dt

        # The body reacts to the normal component only.  This is intentionally
        # evaluated from the previous contact solution before the new implicit d.
        body_acc = (kb * (uu - s[3]) + cb * (ud - s[4]) + fx - s[14] * math.cos(s[17])) / m
        s[4] += dt * body_acc
        s[3] += dt * s[4]

        xc_old = s[7]
        vc_old = s[8]
        fn = max(s[14], 0.)
        vc_new = vc_old
        ideal_tangent = (c[13] == 0. and c[19] == 0. and c[20] == 0. and c[17] == 0. and c[21] == 0. and c[22] == 0. and c[26] == 0. and c[23] > 1e10)
        # Semi-implicit tangential update.  The velocity-dependent friction is
        # evaluated at the new velocity while the normal force is lagged by one
        # outer iteration; the outer loop closes Fn--Ft--cross-compliance.
        for outer in range(3):
            if ideal_tangent:
                vc_new = vs
                xc_new = s[1]
            else:
                xc_new = xc_old + dt * vc_new
            for itn in range(3):
                if not ideal_tangent:
                    ft0 = _friction(fn, vc_new, max(s[13], 1e-9), c) - c[17] * vc_new
                    vc_new = (vc_old + dt / c[24] * (c[23] * (s[1] - xc_old) + dt_t * vs - ft0)) / (1. + dt / c[24] * (dt_t + c[17]))
            if not ideal_tangent:
                xc_new = xc_old + dt * vc_new
            theta = c[21] + c[22] * math.atan(_slope(h_s, dx, xc_new))
            h = _lookup(h_s, dx, xc_new - X0)
            base = s[2] + rp + s[3] + h
            d = max(s[13], 1e-9)
            fe = 0.; ft = 0.; fax = 0.; qsum = 0.
            mod = 1. + c[27] * aa
            for itn in range(3 if ideal_tangent else 7):
                fe = _fe_asym(xc_new, d, vc_new, ln_e, dx, c[19], c[20]) * mod
                qsum = 0.
                for q in range(3):
                    qsum += s[10 + q] * ei[q] + alpha[q] * (fe - s[9]) * gi[q]
                dd = (d - s[13]) / dt
                fn = max((fe + qsum) * (1. + c[4] * dd), 0.)
                ft = _friction(fn, vc_new, d, c) if fn > 0. else 0.
                fax = fn * math.cos(theta) + ft * math.sin(theta)
                target = (base - fax / kr + c[26] * ft / c[23]) * math.cos(theta)
                res = d - target
                hd = 2e-6
                ke = (_fe_asym(xc_new, d + hd, vc_new, ln_e, dx, c[19], c[20]) -
                      _fe_asym(xc_new, max(d - hd, 1e-9), vc_new, ln_e, dx, c[19], c[20])) / (2. * hd) * mod
                if ideal_tangent:
                    slope = 1. + (s[22] + c[4] * max(fn, 0.) / dt) / kr
                    d = max(d - res / slope, 1e-9)
                else:
                    slope = max(1. + (abs(ke) + abs(c[4]) * max(fn, 0.) / dt) / kr, .2)
                    d = max(d - res * math.cos(theta) / slope, 1e-9)
            if base <= 0.:
                d = base * math.cos(theta)
                fe = 0.; fn = 0.; ft = 0.; fax = 0.
                qsum = 0.
                s[10] = 0.; s[11] = 0.; s[12] = 0.
            else:
                fe = _fe_asym(xc_new, d, vc_new, ln_e, dx, c[19], c[20]) * mod
                qsum = 0.
                for q in range(3):
                    qsum += s[10 + q] * ei[q] + alpha[q] * (fe - s[9]) * gi[q]
                dd = (d - s[13]) / dt
                fn = max((fe + qsum) * (1. + c[4] * dd), 0.)
                ft = _friction(fn, vc_new, d, c) if fn > 0. else 0.
                fax = fn * math.cos(theta) + ft * math.sin(theta)
            # One further semi-implicit pass with the newly solved normal force.
            if outer < 2:
                continue
        s[7] = xc_new
        s[8] = vc_new
        fe_prev = s[9]
        s[9] = fe
        for q in range(3):
            s[10 + q] = s[10 + q] * ei[q] + alpha[q] * (fe - fe_prev) * gi[q]
        # The previous expression intentionally uses the old elastic force for
        # q updates; retain the just-solved force separately before overwriting it.
        # Reconstructing the state with zero increment here is the stable limit
        # for the final implicit iterate, and the next tick carries the change.
        s[13] = d
        s[14] = fn
        s[15] = ft
        s[16] = fax
        s[17] = theta

        if abs(t * 1000. - round(t * 1000.)) < 1e-6:
            k = min(kk, noise.shape[0] - 1)
            s[6] += 2e-4 * ns * noise[k, 0]
            s[5] += afc * (fax + s[6] + .08 * ns * noise[k, 1] - s[5])
            s[18] += 2e-4 * ns * noise[k, 2]
            s[19] += afc * (ft * math.cos(theta) - fn * math.sin(theta) + s[18] + .08 * ns * noise[k, 3] - s[19])
            s[20] += 2e-4 * (7.5e-4 / .08) * ns * noise[k, 4]
            a0 = halfwidth(max(d, 1e-12))
            a1 = max(a0 - c[19] * abs(vc_new), .2 * a0)
            a2 = a0 + c[20] * abs(vc_new)
            my_true = fn * (a2 - a1) / 2. + ft * c[29]
            s[21] += afc * (my_true + s[20] + 7.5e-4 * ns * noise[k, 5] - s[21])
            kk += 1
        if j == 0 and nlog < log.shape[0]:
            hd = 2e-6
            ke = (_fe_asym(xc_new, d + hd, vc_new, ln_e, dx, c[19], c[20]) -
                  _fe_asym(xc_new, max(d - hd, 1e-9), vc_new, ln_e, dx, c[19], c[20])) / (2. * hd) * (1. + c[27] * aa)
            kc = ke * r_c
            kl = ke * r_low
            ki = ke * r_inst
            kt = kc * kr / (kc + kr) * math.cos(theta) * math.cos(theta) if kc > 0. else 0.
            log[nlog, 0] = t; log[nlog, 1] = s[1]; log[nlog, 2] = s[2]; log[nlog, 3] = fn
            log[nlog, 4] = s[3]; log[nlog, 5] = s[5]; log[nlog, 6] = d
            log[nlog, 7] = math.exp(_lookup(ln_e, dx, xc_new - X0)) * (1. + c[27] * aa)
            log[nlog, 8] = kt; log[nlog, 9] = kl * kr / (kl + kr) * math.cos(theta) * math.cos(theta) if kl > 0. else 0.; log[nlog, 10] = ki * kr / (ki + kr) * math.cos(theta) * math.cos(theta) if ki > 0. else 0.
            log[nlog, 11] = buf[ndel]; log[nlog, 12] = s[19]; log[nlog, 13] = s[21]
            log[nlog, 14] = uu; log[nlog, 15] = aa; log[nlog, 16] = vc_new; log[nlog, 17] = theta
            log[nlog, 18] = fax; log[nlog, 19] = ft
            s[22] = ke * r_inst
            nlog += 1
        s[2] += dt / max(tau_a, 1e-9) * (buf[0] - s[2])
        s[0] = t + dt
        if not math.isfinite(s[2]) or abs(s[2]) > .5 or fn > 5e3 or abs(s[7] - X0) > 2.:
            return -nlog, kk
    return nlog, kk


def _constants(sc, pl):
    alpha, tau = sc['alpha'], sc['tau']
    def ratio(w):
        return 1. + float(np.sum(alpha * (w * tau) ** 2 / (1. + (w * tau) ** 2)))
    afc = 1. - math.exp(-2. * math.pi * pl['fc'] * 1e-3)
    cp = sc
    return np.array([
        sc['m'], sc['c'], sc['k'], sc['kr'], sc['lam'], pl['tau_a'], afc, sc.get('noise_scale', 1.),
        pl['ndel'], ratio(2. * math.pi * 20.8), ratio(2. * math.pi), 1. + float(np.sum(alpha)), pl['VMAX'],
        cp.get('mu_k', 0.), cp.get('mu_s', 0.), cp.get('v_s', .01), cp.get('v_e', .0005), cp.get('k_d', 0.),
        cp.get('nu', .45), cp.get('k1', 0.), cp.get('k2', 0.), cp.get('theta0', 0.), cp.get('c_slope', 0.),
        cp.get('K_t', 1e12), cp.get('m_t', 1.), cp.get('zeta_t', .5), cp.get('xi_c', 0.), sc.get('beta', 0.),
        cp.get('R', R), cp.get('l_tool', .15),
    ], dtype=float)


def _init(sc, pl):
    x0 = float(np.asarray(sc.get('xcmd', [OLD.X0]))[0])
    lo, hi = 1e-7, .06
    for _ in range(60):
        mid = .5 * (lo + hi)
        if _fe_asym(x0, mid, 0., sc['ln_e'], sc['dx'], 0., 0.) > FD:
            hi = mid
        else:
            lo = mid
    d0 = .5 * (lo + hi)
    z0 = sc['u'][0] - FD / sc['k']
    h0 = _lookup.py_func(sc['h_s'], sc['dx'], x0 - X0) if hasattr(_lookup, 'py_func') else 0.
    y0 = d0 - z0 - h0 + FD / sc['kr']
    s = np.zeros(NS)
    s[1] = x0; s[2] = y0; s[3] = z0; s[5] = FD; s[7] = x0; s[8] = 0.; s[9] = FD
    s[13] = d0; s[14] = FD; s[16] = FD; s[17] = 0.; s[21] = 0.; s[22] = 1e4
    buf = np.full(int(pl['ndel']) + 1, y0)
    return s, buf


def _arrays(sc):
    if 'xcmd' in sc:
        xcmd = np.asarray(sc['xcmd'], float)
    else:
        # Compatibility fallback for old scenes; physics_plant2 itself never
        # reads a v_scan field and all new scenes provide xcmd explicitly.
        n = int(round(sc['T'] * 1000.)) + 2
        xcmd = X0 + float(sc.get('v', 0.)) * np.arange(n) * 1e-3
    return sc['ln_e'], sc['h_s'], sc['dx'], xcmd, sc['u'], sc['fext'], sc['act'], sc['rip'], sc['alpha'], sc['tau']


def _noise(rng, n, sensor):
    """Keep the legacy axial stream byte-for-byte aligned with physics_plant."""
    if not sensor:
        return np.zeros((n, 6))
    axial = rng.standard_normal((n, 2))
    extra = rng.standard_normal((n, 4))
    return np.c_[axial, extra]


def episode(sc, policy, seed, platform=None, sensor=True):
    pl = dict(PLATFORM if platform is None else platform)
    T = sc['T']; nt = int(round(T / pl['Tc']))
    log = np.zeros((nt + 2, LOGC)); s, buf = _init(sc, pl); c = _constants(sc, pl)
    rng = np.random.default_rng(seed)
    noise = _noise(rng, int(T * 1000) + 10, sensor)
    policy.reset(pl); a = _arrays(sc); n = 0; kk = 0
    while n < nt:
        log[n, :] = 0.; log[n, 0] = s[0]; log[n, 1] = s[1]; log[n, 2] = s[2]; log[n, 3] = s[14]; log[n, 5] = s[5]; log[n, 11] = buf[int(pl['ndel'])]
        if n > 0: log[n, 8:11] = log[n - 1, 8:11]
        seen = log[n, 5]
        try:
            ycmd = float(policy(OLD.Obs(log, n + 1)))
        except Exception:
            ycmd = float(buf[int(pl['ndel'])])
            raise
        if not np.isfinite(ycmd): ycmd = float(buf[int(pl['ndel'])])
        n0 = n
        n, kk = _tick(s, buf, ycmd, c, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], noise, kk, log, n)
        if n < 0: return None
        log[n0, 5] = seen
    return log[:n]


def replay(sc, commands, seed=0, sensor=False, platform=None):
    pl = dict(PLATFORM if platform is None else platform)
    log = np.zeros((len(commands) + 2, LOGC)); s, buf = _init(sc, pl); c = _constants(sc, pl)
    T = max(float(sc.get('T', len(commands) * pl['Tc'])), len(commands) * pl['Tc'])
    rng = np.random.default_rng(seed)
    noise = _noise(rng, int(T * 1000) + 10, sensor)
    a = _arrays(sc); n = 0; kk = 0
    for command in commands:
        seen = s[5]; n0 = n
        n, kk = _tick(s, buf, float(command), c, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], noise, kk, log, n)
        if n < 0: return None
        log[n0, 5] = seen
    return log[:n]
