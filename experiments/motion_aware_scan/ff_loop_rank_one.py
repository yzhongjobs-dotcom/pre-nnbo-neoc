"""Exact admissible measurement-gain interval of the force-closed feed-forward loop as a rank-one perturbation (method_reassessment_20261108.md §2.3).

In ff_loop.py's frozen linearised model the closed-loop matrix is affine in kappa_g, A(kg) = A1 + (kg-1) b c^T with A1 = A(1).
Hence det(zI - A(kg)) = det(zI - A1) (1 - (kg-1) L1(z)), L1(z) = c^T (zI - A1)^-1 b, and the stability interval containing kg = 1
is bounded by the reciprocals of the real-axis crossings of L1(e^{j theta}). The script checks this prediction against ff_loop's
spectral-radius scan (tracker gain frozen at its kg = 1 value so the perturbation is exactly rank-one) and reports where |L1| peaks.
Run: python -X utf8 ff_loop_rank_one.py -> results_v2/ff_loop_rank_one.json
"""
import json, os
import numpy as np
import ff_loop as F


def amat(kappa, kg, K):
    A = np.zeros((F.NX, F.NX))
    for i in range(F.NX):
        x = np.zeros(F.NX); x[i] = 1.
        A[:, i] = F.tick(x, 0., kappa, kg, K)
    keep = [i for i in range(F.NX) if i != 3]
    return A[np.ix_(keep, keep)]


def analyse(k_phys, q, kappa):
    K = F.kalman_gain(q, k_phys)
    A1 = amat(kappa, 1.0, K); Bc = amat(kappa, 2.0, K)-A1
    U, S, Vt = np.linalg.svd(Bc); b = U[:, 0]*S[0]; c = Vt[0]
    assert S[1] < 1e-9*S[0]
    th = np.linspace(1e-4, np.pi, 40000); I = np.eye(len(A1))
    L = np.array([c@np.linalg.solve(np.exp(1j*t)*I-A1, b) for t in th])
    idx = np.where(np.sign(L.imag[:-1]) != np.sign(L.imag[1:]))[0]
    re = L.real[idx]; fz = th[idx]/(2*np.pi*F.Tc)
    pos, neg = re > 0, re < 0
    d_hi = 1/re[pos].max() if pos.any() else np.inf
    d_lo = 1/re[neg].min() if neg.any() else -np.inf
    grid = np.arange(0.2, 3.0, 0.005)
    r = [max(abs(e) for e in np.linalg.eigvals(amat(kappa, g, K)) if abs(e-1) > 1e-7) for g in grid]
    i1 = int(np.argmin(abs(grid-1))); j = h = i1
    while j > 0 and r[j-1] < 1: j -= 1
    while h < len(r)-1 and r[h+1] < 1: h += 1
    f = th/(2*np.pi*F.Tc)
    return dict(pred=(round(1+d_lo, 3), round(1+d_hi, 3)), scan=(round(float(grid[j]), 3), round(float(grid[h]), 3)),
                peak_absL=round(float(abs(L).max()), 3), f_peak_Hz=round(float(f[np.argmax(abs(L))]), 2),
                f_hi_cross_Hz=round(float(fz[pos][np.argmax(re[pos])]), 2) if pos.any() else None,
                f_lo_cross_Hz=round(float(fz[neg][np.argmin(re[neg])]), 2) if neg.any() else None)


if __name__ == '__main__':
    out = {}
    for k_phys in (1e3, 1e4, 1e5):
        for mode, q in (('smooth', F.P['q_smooth']), ('maneuver', F.P['q_maneuver'])):
            for kappa in (0.5, 1.0, 2.0):
                key = f'{k_phys:g}|{mode}|{kappa}'
                out[key] = analyse(k_phys, q, kappa); print(key, out[key])
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_v2', 'ff_loop_rank_one.json'), 'w'), indent=1)
