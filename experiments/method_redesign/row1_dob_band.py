"""Row 1 (DOB view of the coupling): stable kappa_g band of the force-closed feed-forward loop versus a
unit-DC first-order Q filter on the motion feed-forward. Same linearised platform model as
motion_aware_scan/ff_loop.py (one-period delay, 19 ms servo, 30 Hz filter at 1 ms sub-steps, A1 at
feedback ratio kappa, exact one-step filter inversion, steady-state CA Kalman, 12 ms + tau_s predictor),
with one extra state q: q_j = a q_{j-1} + (1-a) sp_j, a = exp(-2 pi f_Q Tc); the command uses -(q_j - q_{j-1}).
f_Q = inf reproduces ff_loop.py. Run: python -X utf8 row1_dob_band.py -> results/row1_dob_band.json
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan'))
import ff_loop as FL

NX = FL.NX+1

def tick(x, kappa, kappa_g, K, a):
    y, cp, Ff, I, wd, cmd, mp = x[:7]
    Fp, yp = x[7:9]; xh = x[9:12]; spp = x[12]; q = x[13]
    e, ep = -Ff, -mp
    eD = math.exp(-FL.D1*FL.Tc/FL.M1)
    wd2 = wd*eD+(1-eD)*e/FL.D1
    fin = (Ff-FL.DEC*Fp)/(1-FL.DEC)
    ybar = yp+FL.FRAC*(y-yp)
    Z = kappa_g/kappa*fin-ybar
    xm = FL.At@xh; xh2 = xm+K*(Z-xm[0]); sp = FL.CP@xh2
    q2 = a*q+(1-a)*sp
    cmd2 = cmd+0.5*FL.Tc*(wd+wd2)+FL.KPF1*(e-ep)-(q2-q)
    cp2 = cmd
    for _ in range(FL.NS):
        Ff = Ff+FL.ALPHA*(kappa*(y-0.)-Ff)
        y = y+FL.BETA*(cp2-y)
    return np.concatenate([[y, cp2, Ff, I, wd2, cmd2, -e], [x[2], x[0]], xh2, [sp, q2]])

def radius(kappa, kappa_g, K, a):
    A = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.; A[:, i] = tick(x, kappa, kappa_g, K, a)
    keep = [i for i in range(NX) if i != 3]
    ev = np.linalg.eigvals(A[np.ix_(keep, keep)]); ev = ev[np.abs(ev-1) > 1e-7]
    return float(np.max(np.abs(ev)))

def band(kappa, k_phys, q, a, grid):
    r = [radius(kappa, kg, FL.kalman_gain(q, k_phys/kg), a) for kg in grid]
    i1 = int(np.argmin(np.abs(grid-1.)))
    if r[i1] >= 1: return None
    j = h = i1
    while j > 0 and r[j-1] < 1: j -= 1
    while h < len(r)-1 and r[h+1] < 1: h += 1
    return float(grid[j]), float(grid[h])

if __name__ == '__main__':
    grid = np.round(np.concatenate([np.arange(0.2, 1.0, 0.025), np.arange(1.0, 3.01, 0.025)]), 3)
    out = {}
    for fQ in (np.inf, 20., 10., 5., 3., 2., 1.):
        a = 0. if not np.isfinite(fQ) else math.exp(-2*math.pi*fQ*FL.Tc)
        worst_lo, worst_hi = 0., 9.
        for k_phys in (1e3, 1e4, 1e5):
            for mode in ('q_smooth', 'q_maneuver'):
                for kappa in (0.5, 1.0):
                    b = band(kappa, k_phys, FL.P[mode], a, grid)
                    out['%g|%g|%s|%.1f' % (fQ, k_phys, mode, kappa)] = b
                    if b is None: worst_lo, worst_hi = 1., 1.
                    else: worst_lo, worst_hi = max(worst_lo, b[0]), min(worst_hi, b[1])
        out['common|%g' % fQ] = (worst_lo, worst_hi)
        print('f_Q = %5s Hz  common stable kappa_g band over k in [1e3,1e5], both modes, kappa in {0.5,1}: [%.3f, %.3f]' % (fQ, worst_lo, worst_hi), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'row1_dob_band.json'), 'w'), indent=1)
