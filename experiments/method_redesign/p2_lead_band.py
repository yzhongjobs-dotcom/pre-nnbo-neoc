"""Stable kappa_g band of the force-closed feed-forward loop when only the servo-lag lead is low-passed (the
candidate set of the supervisory selection), same linearised platform as motion_aware_scan/ff_loop.py and
row1_dob_band.py: feed-forward = s(t+L) + tau_s * v_f, v_f_j = a v_f_{j-1} + (1-a) sdot(t+L), a = exp(-2 pi f Tc).
f = inf reproduces ff_loop.py. Common band over k in {1e3,1e4,1e5}, both IMM modes, kappa in {0.5, 1}.
Run: python -X utf8 p2_lead_band.py -> results/p2_lead_band.txt
"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'motion_aware_scan'))
import ff_loop as FL
NX = FL.NX+1
CPOS = np.array([1, FL.L, FL.L*FL.L/2]); CVEL = np.array([0, 1, FL.L])
def tick(x, kappa, kappa_g, K, a):
    y, cp, Ff, I, wd, cmd, mp = x[:7]
    Fp, yp = x[7:9]; xh = x[9:12]; spp = x[12]; vf = x[13]
    e, ep = -Ff, -mp
    eD = math.exp(-FL.D1*FL.Tc/FL.M1)
    wd2 = wd*eD+(1-eD)*e/FL.D1
    fin = (Ff-FL.DEC*Fp)/(1-FL.DEC)
    ybar = yp+FL.FRAC*(y-yp)
    Z = kappa_g/kappa*fin-ybar
    xm = FL.At@xh; xh2 = xm+K*(Z-xm[0])
    vf2 = a*vf+(1-a)*(CVEL@xh2)
    sp = CPOS@xh2+FL.tau_s*vf2
    cmd2 = cmd+0.5*FL.Tc*(wd+wd2)+FL.KPF1*(e-ep)-(sp-spp)
    cp2 = cmd
    for _ in range(FL.NS):
        Ff = Ff+FL.ALPHA*(kappa*(y-0.)-Ff)
        y = y+FL.BETA*(cp2-y)
    return np.concatenate([[y, cp2, Ff, I, wd2, cmd2, -e], [x[2], x[0]], xh2, [sp, vf2]])
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
    out = []
    for f in [float(x) for x in sys.argv[1:]] or (np.inf, 10., 6., 3.):
        a = 0. if not np.isfinite(f) else math.exp(-2*math.pi*f*FL.Tc)
        lo, hi = 0., 9.
        for k_phys in (1e3, 1e4, 1e5):
            for mode in ('q_smooth', 'q_maneuver'):
                for kappa in (0.5, 1.0):
                    b = band(kappa, k_phys, FL.P[mode], a, grid)
                    lo, hi = (1., 1.) if b is None else (max(lo, b[0]), min(hi, b[1]))
        out.append('lead low-pass %5s Hz: common stable kappa_g band [%.3f, %.3f]' % (f, lo, hi)); print(out[-1], flush=True)
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'p2_lead_band.txt'), 'a').write('\n'.join(out)+'\n')
