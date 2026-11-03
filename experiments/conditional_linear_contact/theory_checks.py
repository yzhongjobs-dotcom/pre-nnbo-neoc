"""理论与合成检验；所有数字均由本文件实际计算得到。"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
from conditional_filter import FourModeIMM  # noqa: E402
import cohort_eval as CE  # noqa: E402

OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)


def weights(ns=8, fc=30.):
    afc = 1. - math.exp(-2. * math.pi * fc * .001)
    a = 1. - afc
    m = np.arange(1, ns + 1, dtype=float)
    w = afc * a ** (ns - m)
    tau = (m - 1.) * 1e-3
    return afc, a, w, np.array([w.sum(), w @ tau, w @ (tau * tau / 2.)])


def synth(seed, feedback=True, N=650, k=2.5e4, q=3e-5, noise=.04):
    """与增量观测方程一致的闭环合成数据，包含简单一阶伺服。"""
    rng = np.random.default_rng(seed)
    dt, ns, tau = .008, 8, .019
    afc, a_f, w, W = weights(ns)
    rho = math.exp(-dt / tau)
    qF = k * k * q
    Q3 = qF * np.array([[dt ** 5 / 20., dt ** 4 / 8., dt ** 3 / 6.],
                        [dt ** 4 / 8., dt ** 3 / 3., dt ** 2 / 2.],
                        [dt ** 3 / 6., dt ** 2 / 2., dt]])
    x = np.array([k, 4., 0., 0.], dtype=float)
    y = 0.; command = 0.; F = 4.
    Hs, zs, Fs, xs, dys = [], [], [], [], []
    for _ in range(N):
        e = 4. - F
        if feedback:
            command += .35 * e / 7e4
        # A known small displacement excitation supplies the information that
        # force feedback alone cannot create near a perfectly regulated force.
        command += 2.5e-4 * math.sin(2. * math.pi * 20.8 * (_ * dt + dt))
        y1 = y + (1. - rho) * (command - y)
        dy = y1 - y
        fr = -np.expm1(-np.arange(ns) * .001 / tau) / (1. - math.exp(-dt / tau))
        Y = float(w @ (fr * dy))
        H = np.array([Y, W[0], W[1], W[2]])
        z = float(H @ x + noise * rng.standard_normal())
        F1 = a_f ** ns * F + z
        Fmat = np.array([[1., 0., 0., 0.],
                         [dy, 1., dt, dt * dt / 2.],
                         [0., 0., 1., dt], [0., 0., 0., 1.]])
        eta = np.r_[0., rng.multivariate_normal(np.zeros(3), Q3)]
        x = Fmat @ x + eta
        x[0] = k  # known-constant-k data-generating process
        Hs.append(H); zs.append(z); Fs.append(Fmat); xs.append(x.copy()); dys.append(dy)
        y, F = y1, F1
    return np.asarray(Hs), np.asarray(zs), np.asarray(Fs), np.asarray(xs)


def run_synth(Hs, zs, Fs, k, q=3e-5, noise=.04):
    filt = FourModeIMM(dt=.008, q_smooth=q, q_maneuver=q * 8., p_stay_motion=.99,
                       q_lnk0=1e-8, q_lnk1=1e-6, p_stay_stiff=.99, R=noise ** 2,
                       r_scale=1.)
    filt.initialize_if_needed(k, 4.)
    est, sig = [], []
    for H, z, F in zip(Hs, zs, Fs):
        filt.update_interval(float(z), H, F, k)
        est.append(filt.x[0]); sig.append(math.sqrt(max(filt.P[0, 0], 0.)))
    return np.asarray(est), np.asarray(sig)


def synthetic_check():
    k = 2.5e4
    out = {}
    for feedback in (False, True):
        rows = []
        for seed in range(40):
            H, z, F, x = synth(700000 + seed, feedback=feedback, k=k, N=450)
            est, sig = run_synth(H, z, F, k)
            sl = slice(100, None)
            err = est[sl] - k
            ss = np.maximum(sig[sl], 1e-12)
            rows.append(dict(mean_error=float(np.mean(err)), rmse=float(np.sqrt(np.mean(err ** 2))),
                             normalized_error_variance=float(np.mean((err / ss) ** 2)),
                             coverage95=float(np.mean(np.abs(err) <= 1.96 * ss))))
        out['feedback' if feedback else 'no_feedback'] = {
            k: float(np.mean([r[k] for r in rows])) for k in rows[0]
        }
    return dict(k=k, repetitions=40, intervals=450, burn_intervals=100, results=out,
                interpretation='闭环与无反馈均按同一增量模型生成；这是条件后验校准检验，不是实际对象证明。')


def frequency_check():
    afc, a_f, w, W = weights()
    dt = .008
    H = W.copy()
    q = float(CE.FINAL['q_smooth'])
    R = .04 ** 2
    rows = []
    for k in (1e3, 1e4, 1e5):
        qF = k * k * q
        Q = qF * np.array([[dt ** 5 / 20., dt ** 4 / 8., dt ** 3 / 6.],
                            [dt ** 4 / 8., dt ** 3 / 3., dt ** 2 / 2.],
                            [dt ** 3 / 6., dt ** 2 / 2., dt]])
        A = np.array([[1., dt, dt ** 2 / 2.], [0., 1., dt], [0., 0., 1.]])
        P = np.eye(3)
        for _ in range(20000):
            Pp = A @ P @ A.T + Q
            S = float(H @ Pp @ H + R)
            K = Pp @ H / S
            I = np.eye(3) - np.outer(K, H)
            Pn = (I @ Pp @ I.T + R * np.outer(K, K) + (I @ Pp @ I.T + R * np.outer(K, K)).T) / 2.
            if np.max(abs(Pn - P)) < 1e-18:
                P = Pn
                break
            P = Pn
        Pp = A @ P @ A.T + Q
        S = float(H @ Pp @ H + R)
        K = Pp @ H / S
        freqs = [0.2, 1., 5., 10., 20., 40., 55.]
        for f in freqs:
            omega = 2. * math.pi * f
            z = np.sin(omega * dt * np.arange(1600))
            xp = np.zeros(3); innov = []
            for zi in z:
                v = float(zi - H @ xp)
                innov.append(v)
                xp = A @ (xp + K * v)
            tt = np.arange(400, 1600) * dt
            Afit = np.column_stack([np.sin(omega * tt), np.cos(omega * tt)])
            coef = np.linalg.lstsq(Afit, np.asarray(innov)[400:], rcond=None)[0]
            mag = float(np.linalg.norm(coef))
            wc = (k * k * q / (R * .001)) ** (1. / 6.)
            butter = float(omega ** 3 / math.sqrt(omega ** 6 + wc ** 6))
            rows.append(dict(k=k, frequency_hz=f, omega_c=wc, butterworth=butter,
                             innovation_magnitude=mag, absolute_error=mag - butter,
                             relative_error=(mag - butter) / max(butter, 1e-12)))
    rng = np.random.default_rng(3)
    Y = rng.normal(size=200); innovation = rng.normal(size=200)
    alpha = .013
    k0 = 1000.
    k_a = k0 + alpha * np.cumsum(Y * innovation)
    k_b = k0 + alpha * np.cumsum(Y * innovation)
    return dict(weights=dict(afc=afc, a=a_f, W=W.tolist()), rows=rows,
                k_channel_innovation_integral_max_abs_difference=float(np.max(abs(k_a - k_b))))


def prior_limit_check():
    rng = np.random.default_rng(5)
    n = 40; t = np.linspace(-1., 1., n)
    B = np.column_stack([t ** d for d in range(4)])
    Y = np.sin(2.7 * t) + .2 * t ** 4
    z = 4.2 * Y + B @ np.array([.7, -.3, .2, .1]) + .03 * rng.standard_normal(n)
    Q, *_ = np.linalg.qr(B, mode='reduced')
    Yp = Y - Q @ (Q.T @ Y); zp = z - Q @ (Q.T @ z)
    k_proj = float(Yp @ zp / (Yp @ Yp))
    rows = []
    for lam in (1e-3, 1e-1, 1., 1e2, 1e4, 1e8, 1e12):
        A = np.column_stack([Y, B])
        penalty = np.diag([0., 1. / lam ** 2, 1. / lam ** 2, 1. / lam ** 2, 1. / lam ** 2])
        theta = np.linalg.solve(A.T @ A + penalty, A.T @ z)
        rows.append(dict(prior_sd=lam, k=float(theta[0]), abs_difference=float(abs(theta[0] - k_proj))))
    return dict(k_projection=k_proj, rows=rows,
                conditions=['B has full column rank', 'Y is not in span(B)',
                            'Gaussian/weighted least-squares likelihood and independent nuisance prior',
                            'prior variance of every polynomial coefficient tends to infinity'])


def transition_table():
    q_s, q_m = float(CE.FINAL['q_smooth']), float(CE.FINAL['q_maneuver'])
    R, dt = .04 ** 2, .001
    rows = []
    for qname, q in [('q_smooth', q_s), ('q_maneuver', q_m)]:
        for k in (1e3, 1e4, 1e5):
            omega = (k * k * q / (R * dt)) ** (1. / 6.)
            rows.append(dict(q=qname, q_value=q, k=k, omega_c=omega, frequency_hz=omega / (2 * math.pi)))
    return dict(noise_variance=R, dt=dt, rows=rows)


def main():
    result = dict(synthetic=synthetic_check(), frequency=frequency_check(),
                  infinite_prior=prior_limit_check(), transition_frequencies=transition_table())
    path = OUT / 'theory_checks.json'
    path.write_text(json.dumps(result, indent=2, default=float), encoding='utf8')
    print(path)
    print(json.dumps(result['synthetic'], indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
