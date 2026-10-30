"""Minimal linear check for the problem-stage mechanism claims (analysis, not a method benchmark).

Plant (C1 benchmark settings): Tc=8 ms, one-cycle command delay, first-order servo 19 ms,
force f = k (y - s) + noise, first-order 30 Hz force filter. Admittance (K=0):
M yc'' + D yc' = e + ki I, I' = e, y_cmd = yc + KPF e + u.

(1) Loop gain scales with kappa = k / k_design: report the stable kappa interval and the
    0-dB crossover of the stiffness-normalized design, and the fixed-design range over 30x.
(2) Windowed stiffness regression while the surface moves and force is regulated:
    static-surface model (intercept only) vs cubic-drift projection, with/without excitation.
"""
import numpy as np

Tc, tau_s, fc = 0.008, 0.019, 30.0
a_s = np.exp(-Tc / tau_s)
a_f = np.exp(-2 * np.pi * fc * Tc)


def gains(kd):
    # stiffness-normalized design: whole admittance scaled by kd (M, D ~ kd; KPF ~ 1/kd)
    # reference values follow the C1 rule prior at kd = 2.4e4 (D = kd/80, KPF = 2/kd), M scaled with D
    return dict(M=0.3 * kd / 2.4e4, D=kd / 80, ki=30.0, KPF=2.0 / kd)


def simulate(k_true, kd, T=8.0, surf=None, exc=0.0, fexc=20.8, Fd=4.0, noise=0.08, seed=0, ff=None):
    g = gains(kd)
    n = int(T / Tc)
    rng = np.random.default_rng(seed)
    yc = vc = I = 0.0
    y = Fd / k_true  # start in equilibrium against s=0
    yc = y
    cmd_prev = y
    Ff = Fd
    out = np.zeros((n, 4))
    for i in range(n):
        t = i * Tc
        s = surf(t) if surf else 0.0
        f = max(k_true * (y - s), 0.0)
        y_meas = y  # position that produced f (same sample)
        Ff = a_f * Ff + (1 - a_f) * (f + noise * rng.standard_normal())
        e = Fd - Ff
        I += e * Tc
        # semi-implicit admittance step (stiff D/M); explicit Euler is unstable at 8 ms
        vc = (vc + (e + g['ki'] * I) * Tc / g['M']) / (1 + g['D'] * Tc / g['M'])
        yc += vc * Tc
        u = exc * np.sin(2 * np.pi * fexc * t)
        if ff is not None:
            u += ff(t)
        cmd = yc + g['KPF'] * e + u
        # one-cycle delay then first-order servo
        y = a_s * y + (1 - a_s) * cmd_prev
        cmd_prev = cmd
        out[i] = (t, y_meas, f, Ff)
    return out


def stable_kappa(kd=2.4e4):
    # closed-loop stability of the linearized loop vs kappa (k_true = kappa*kd), surface static
    res = []
    for kappa in [0.03, 0.1, 0.3, 1, 2, 3, 4, 6, 8, 12]:
        o = simulate(kappa * kd, kd, T=6.0, noise=0.0, surf=lambda t: 0.5e-3 * (t > 1.0))
        tail = o[int(3.0 / Tc):, 2]
        osc = np.ptp(tail)
        # 10-90 % rise of the force after the 0.5 mm surface step
        seg = o[int(1.0 / Tc):int(3.0 / Tc), 2]
        res.append((kappa, osc, seg.max()))
    return res


def regress(o, k_true, window, deg, t0):
    # windowed regression of measured force on probe position over `window` samples ending at index
    t, y, f, Ff = o.T
    rows = []
    idx = np.arange(int(t0 / Tc), len(t), window)
    for j in idx:
        sl = slice(j - window, j)
        Y = y[sl]; Z = f[sl] + 0.08 * np.random.default_rng(j).standard_normal(window)
        tt = (t[sl] - t[sl].mean()) / (window * Tc)
        B = np.vstack([tt ** p for p in range(deg + 1)]).T
        P = np.eye(window) - B @ np.linalg.pinv(B)
        Yp, Zp = P @ Y, P @ Z
        if Yp @ Yp < 1e-16:
            continue
        rows.append((Yp @ Zp) / (Yp @ Yp) / k_true)
    r = np.array(rows)
    return np.median(r), np.percentile(r, 10), np.percentile(r, 90)


if __name__ == '__main__':
    print('(1) stiffness-normalized loop, kappa = k_true / k_design (surface step 0.5 mm, no noise)')
    for kappa, osc, pk in stable_kappa():
        print(f'   kappa={kappa:5.2f}  steady ptp={osc:9.3g} N  peak after step={pk:9.3g} N')

    k = 2.4e4
    surf = lambda t: 1e-3 * np.sin(2 * np.pi * 0.25 * t) + 5e-3 * 0.5 * (1 - np.cos(np.pi * np.clip((t - 2) / 1.5, 0, 1)))
    print('\n(2) window regression k_hat/k (median [p10,p90]); surface: 1 mm 0.25 Hz breathing + 5 mm smooth push')
    for label, exc in (('no excitation', 0.0), ('20.8 Hz excitation 0.02 mm', 2e-5)):
        o = simulate(k, k, T=8.0, surf=surf, exc=exc)
        for W in (16, 64):
            for deg, name in ((0, 'static surface (intercept)'), (3, 'cubic drift projection')):
                m, lo, hi = regress(o, k, W, deg, 1.0)
                print(f'   {label:28s} window={W*Tc*1e3:5.0f} ms  {name:28s} {m:7.3f} [{lo:7.3f},{hi:7.3f}]')
