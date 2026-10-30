"""Frozen-gain analysis of the force-closed motion feed-forward loop (paper draft §三.5).

Structure: the motion measurement is Z = F_in/k_ref - y_bar. With f = k (y - s), Z = (k/k_ref - 1) y - (k/k_ref) s,
so for k_ref = k the measurement does not depend on the probe command and the feed-forward loop is open
(the closed-loop poles are those of the force loop and of the tracker). For k_ref != k it closes with gain
(kappa_g - 1), kappa_g = k/k_ref. This script builds the linearised 8 ms platform model of design_A1.py
(one-period delay, 19 ms servo and 30 Hz force filter at 1 ms sub-steps, unit-stiffness admittance A1 at
feedback ratio kappa = k/k_bar) plus the exact one-step filter inversion, filter-weighted encoder
interpolation, a steady-state constant-acceleration Kalman filter (one IMM mode frozen) and the 12 ms + tau_s
predictor, and reports the closed-loop spectral radius over (kappa, kappa_g).
Run: python -X utf8 ff_loop.py -> results_v2/ff_loop.json
"""
import json, math, os
import numpy as np
import design_A1 as D

Tc, tau_s, NS = D.Tc, D.tau_s, D.NSUB
BETA, ALPHA = D.BETA, D.ALPHA
P = json.load(open(os.path.join(D.HERE, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
M1, D1, KPF1 = P['M1'], P['D1'], P['KPF1']
NOISE = 0.04
DEC = (1-ALPHA)**NS
_m = np.arange(NS)
_w = ALPHA*(1-ALPHA)**(NS-1-_m)
FRAC = float(_w @ (-np.expm1(-_m*1e-3/tau_s)/(-math.expm1(-Tc/tau_s)))/_w.sum())
L = 1.5*Tc
At = np.array([[1, Tc, Tc*Tc/2], [0, 1, Tc], [0, 0, 1.]])
Q1 = np.array([[Tc**5/20, Tc**4/8, Tc**3/6], [Tc**4/8, Tc**3/3, Tc**2/2], [Tc**3/6, Tc**2/2, Tc]])
CP = np.array([1, L+tau_s, L*L/2+tau_s*L])          # s(t+L) + tau_s * sdot(t+L)
NX = 7+2+3+1                                          # A1 loop | F_prev, y_prev | tracker | sp_prev


def kalman_gain(q, k_phys):
    R = (NOISE/k_phys)**2*P['r_scale']
    Pm = np.eye(3)*1e-6
    for _ in range(20000):
        Pp = At@Pm@At.T+q*Q1
        K = Pp[:, 0]/(Pp[0, 0]+R)
        Pn = (np.eye(3)-np.outer(K, [1, 0, 0]))@Pp
        if np.max(np.abs(Pn-Pm)) < 1e-30:
            break
        Pm = Pn
    return K


def tick(x, s, kappa, kappa_g, K):
    y, cp, Ff, I, wd, cmd, mp = x[:7]
    Fp, yp = x[7:9]
    xh = x[9:12]; spp = x[12]
    e, ep = -Ff, -mp                                  # force error (target removed: linear deviation)
    eD = math.exp(-D1*Tc/M1)
    wd2 = wd*eD+(1-eD)*e/D1
    fin = (Ff-DEC*Fp)/(1-DEC)
    ybar = yp+FRAC*(y-yp)
    Z = kappa_g/kappa*fin-ybar                         # F_in / k_ref in units of k_bar = 1
    xm = At@xh
    xh2 = xm+K*(Z-xm[0])
    sp = CP@xh2
    cmd2 = cmd+0.5*Tc*(wd+wd2)+KPF1*(e-ep)-(sp-spp)
    cp2 = cmd
    for _ in range(NS):
        Ff = Ff+ALPHA*(kappa*(y-s)-Ff)
        y = y+BETA*(cp2-y)
    return np.concatenate([[y, cp2, Ff, I, wd2, cmd2, -e], [x[2], x[0]], xh2, [sp]])


def radius(kappa, kappa_g, K):
    A = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.
        A[:, i] = tick(x, 0., kappa, kappa_g, K)
    keep = [i for i in range(NX) if i != 3]           # unused integral state
    ev = np.linalg.eigvals(A[np.ix_(keep, keep)])
    # the accumulated command and tracker position have a structural unit pole (absolute level);
    # remove eigenvalues at exactly 1 that are uncontrollable level modes
    ev = ev[np.abs(ev-1) > 1e-7]
    return float(np.max(np.abs(ev)))


if __name__ == '__main__':
    kg_grid = np.round(np.concatenate([np.arange(0.2, 1.0, 0.05), np.arange(1.0, 4.01, 0.1)]), 3)
    out = {}
    for k_phys in (1e3, 3e3, 1e4, 3e4, 1e5):
        for mode, q in (('smooth', P['q_smooth']), ('maneuver', P['q_maneuver'])):
            K = kalman_gain(q, k_phys)
            for kappa in (0.5, 1.0, 2.0):
                r = [radius(kappa, kg, K) for kg in kg_grid]
                st = [float(kg) for kg, v in zip(kg_grid, r) if v < 1-1e-9]
                rng = (min(st), max(st)) if st else None
                # stable interval containing kappa_g = 1
                lo = hi = None
                if radius(kappa, 1.0, K) < 1:
                    i1 = list(kg_grid).index(1.0)
                    j = i1
                    while j > 0 and r[j-1] < 1: j -= 1
                    h = i1
                    while h < len(r)-1 and r[h+1] < 1: h += 1
                    lo, hi = float(kg_grid[j]), float(kg_grid[h])
                key = '%g|%s|%.1f' % (k_phys, mode, kappa)
                out[key] = dict(rho_at_1=radius(kappa, 1.0, K), interval=(lo, hi), K=K.tolist())
                print(key, 'rho(kg=1)=%.4f' % out[key]['rho_at_1'], 'stable k/k_ref in', (lo, hi))
    json.dump(out, open(os.path.join(D.HERE, 'results_v2', 'ff_loop.json'), 'w'), indent=1)
