"""Exact linear stability boundary and bandwidth of the stiffness-normalized loop (companion to gain_and_bias.py).

Discrete closed loop of gain_and_bias.simulate without noise, saturation or contact loss: states
(y, cmd_prev, Ff, I, vc, yc). Design gains at k_design; plant stiffness k = kappa * k_design.
Reports the spectral radius over kappa, the stable kappa interval, and the -3 dB bandwidth of the
force-reference-to-force response (a proxy for following ability).
"""
import numpy as np
from gain_and_bias import Tc, a_s, a_f, gains


def closed_loop(kappa, kd=2.4e4):
    g = gains(kd); k = kappa * kd
    M, D, ki, KPF = g['M'], g['D'], g['ki'], g['KPF']
    den = 1 + D * Tc / M
    # x = [y, cmd_prev, Ff, I, vc, yc]; input Fd, disturbance s enter via f = k (y - s)
    A = np.zeros((6, 6)); Bf = np.zeros(6)
    # Ff+ = a_f Ff + (1-a_f) k y
    A[2, 0] = (1 - a_f) * k; A[2, 2] = a_f
    # e = Fd - Ff+
    def e_row():
        r = -A[2].copy(); return r, 1.0
    er, eb = e_row()
    A[3] = np.eye(6)[3] + Tc * er; Bf[3] = Tc * eb                        # I+
    A[4] = (np.eye(6)[4] + (er + ki * A[3]) * Tc / M) / den               # vc+
    Bf[4] = (eb + ki * Bf[3]) * Tc / M / den
    A[5] = np.eye(6)[5] + Tc * A[4]; Bf[5] = Tc * Bf[4]                   # yc+
    cmd = A[5] + KPF * er; cmdb = Bf[5] + KPF * eb
    A[0] = a_s * np.eye(6)[0] + (1 - a_s) * np.eye(6)[1]                  # y+ uses cmd_prev
    A[1] = cmd; Bf[1] = cmdb                                              # cmd_prev+ = cmd
    return A, Bf


def bandwidth(kappa, kd=2.4e4):
    A, Bf = closed_loop(kappa, kd)
    C = np.zeros(6); C[0] = kappa * kd                                     # actual force k*y
    w = np.logspace(-1, np.log10(np.pi / Tc), 4000)
    H = np.array([C @ np.linalg.solve(np.exp(1j * wi * Tc) * np.eye(6) - A, Bf) for wi in w])
    H0 = abs(H[0]); idx = np.where(abs(H) < H0 / np.sqrt(2))[0]
    return w[idx[0]] / (2 * np.pi) if len(idx) else np.nan


if __name__ == '__main__':
    ks = np.logspace(np.log10(0.01), np.log10(10), 400)
    rho = np.array([max(abs(np.linalg.eigvals(closed_loop(x)[0]))) for x in ks])
    stable = ks[rho < 1]
    print(f'stable kappa interval (grid 0.01-10): [{stable.min():.3f}, {stable.max():.3f}]')
    for x in (0.033, 0.1, 0.3, 0.5, 1.0, 1.5, 1.8):
        r = max(abs(np.linalg.eigvals(closed_loop(x)[0])))
        print(f'kappa={x:5.3f}  spectral radius={r:.4f}  closed-loop force bandwidth={bandwidth(x):6.2f} Hz')
