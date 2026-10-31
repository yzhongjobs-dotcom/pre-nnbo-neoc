"""Sufficient condition for the time-varying loop: frozen stability plus bounded parameter variation (paper draft §三).

Model: the linearised 8 ms platform of design_A1/ff_loop (one-period delay, 19 ms servo, 30 Hz force filter at 1 ms
sub-steps), unit-stiffness admittance A1 at feedback ratio kappa = k/k_bar, exact one-step filter inversion,
filter-weighted encoder interpolation, and the motion channel with the actual two-model IMM linearised at frozen mode
probability mu (mixing, model-conditioned steady-state Kalman gains, output combination) and the 12 ms + tau_s predictor;
measurement ratio kappa_g = k/k_ref; the tracker gains use the controller's R = r(sigma/k_ref)^2.
Parameters p = (ln kappa, ln kappa_g, mu, ln k). The closed loop x+ = A(p) x has structural unit eigenvalues (absolute
level of the accumulated command and of the tracker position, which do not affect the force); we verify that their
eigen-subspace does not depend on p and remove it exactly. On the reduced system, with P(p) solving
A'PA - P = -I and a = sup ||A||, V = x'P(p_t)x decreases whenever ||P(p_{t+1}) - P(p_t)|| < 1/a^2, which holds if
sum_i L_i |dp_i| < 1/a^2 with L_i = sup ||dP/dp_i|| (finite differences over the grid). Jumps (immediate rises of k_bar)
are covered by an average dwell time N_d >= ln(c2/c1) / (-ln(1 - 1/c2)) ticks between jumps.
Run: python -X utf8 theory_ltv.py -> results_v2/theory_ltv.json
"""
import json, math, os
import numpy as np
from scipy.linalg import solve_discrete_lyapunov, null_space
import design_A1 as D

Tc, tau_s, NS = D.Tc, D.tau_s, D.NSUB
BETA, ALPHA = D.BETA, D.ALPHA
FIN = json.load(open(os.path.join(D.HERE, 'tuned_v2', 'Proposed-IMMb.json')))['best']['params']
M1, D1, KPF1 = FIN['M1'], FIN['D1'], FIN['KPF1']
QS, QM, PST, RS = FIN['q_smooth'], FIN['q_maneuver'], FIN['p_stay'], FIN['r_scale']
NOISE = 0.04
DEC = (1-ALPHA)**NS
_m = np.arange(NS); _w = ALPHA*(1-ALPHA)**(NS-1-_m)
FRAC = float(_w @ (-np.expm1(-_m*1e-3/tau_s)/(-math.expm1(-Tc/tau_s)))/_w.sum())
L = 1.5*Tc
At = np.array([[1, Tc, Tc*Tc/2], [0, 1, Tc], [0, 0, 1.]])
Q1 = np.array([[Tc**5/20, Tc**4/8, Tc**3/6], [Tc**4/8, Tc**3/3, Tc**2/2], [Tc**3/6, Tc**2/2, Tc]])
CP = np.array([1, L+tau_s, L*L/2+tau_s*L])
NX = 7+2+6+1


def kgain(q, k_ref):
    R = (NOISE/k_ref)**2*RS; Pm = np.eye(3)*1e-6
    for _ in range(20000):
        Pp = At@Pm@At.T+q*Q1; K = Pp[:, 0]/(Pp[0, 0]+R)
        Pn = (np.eye(3)-np.outer(K, [1, 0, 0]))@Pp
        if np.max(np.abs(Pn-Pm)) < 1e-30: break
        Pm = Pn
    return K


def tick(x, kappa, kappa_g, mu, Ks, Km):
    y, cp, Ff, I, wd, cmd, mp = x[:7]; Fp, yp = x[7:9]
    xs, xm = x[9:12], x[12:15]; spp = x[15]
    e, ep = -Ff, -mp
    eD = math.exp(-D1*Tc/M1); wd2 = wd*eD+(1-eD)*e/D1
    fin = (Ff-DEC*Fp)/(1-DEC); ybar = yp+FRAC*(y-yp)
    Z = kappa_g/kappa*fin-ybar
    ms, mm = 1-mu, mu
    cs = PST*ms+(1-PST)*mm; cm = (1-PST)*ms+PST*mm
    mix_s = (PST*ms*xs+(1-PST)*mm*xm)/cs; mix_m = ((1-PST)*ms*xs+PST*mm*xm)/cm
    ps, pm = At@mix_s, At@mix_m
    xs2 = ps+Ks*(Z-ps[0]); xm2 = pm+Km*(Z-pm[0])
    sp = CP@(ms*xs2+mm*xm2)
    cmd2 = cmd+0.5*Tc*(wd+wd2)+KPF1*(e-ep)-(sp-spp)
    cp2 = cmd
    for _ in range(NS):
        Ff = Ff+ALPHA*(kappa*y-Ff); y = y+BETA*(cp2-y)
    return np.concatenate([[y, cp2, Ff, I, wd2, cmd2, -e], [x[2], x[0]], xs2, xm2, [sp]])


def Amat(lkap, lkg, mu, lk):
    kappa, kappa_g, k = math.exp(lkap), math.exp(lkg), math.exp(lk)
    k_ref = k/kappa_g
    Ks, Km = kgain(QS, k_ref), kgain(QM, k_ref)
    A = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.; A[:, i] = tick(x, kappa, kappa_g, mu, Ks, Km)
    keep = [i for i in range(NX) if i != 3]            # the integral state is unused (k_i = 0)
    return A[np.ix_(keep, keep)]


if __name__ == '__main__':
    # design region (paper property 7): kappa in [0.5, 1.4]; the wider kappa = 2 row is reported separately below
    grid = dict(lkap=np.log([0.5, 0.7, 1.0, 1.2, 1.4]), lkg=np.log([0.85, 0.92, 1.0, 1.08, 1.15]),
                mu=np.array([0., .25, .5, .75, 1.]), lk=np.log([1e3, 3e3, 1e4, 3e4, 1e5]))
    pts = [(a, b, c, d) for a in grid['lkap'] for b in grid['lkg'] for c in grid['mu'] for d in grid['lk']]
    As = {pt: Amat(*pt) for pt in pts}
    # unit-eigenvalue subspace: invariance check across the grid
    def unit_space(A):
        w, V = np.linalg.eig(A); U = V[:, np.abs(w-1) < 1e-6]
        return np.linalg.qr(np.real(U))[0] if U.shape[1] else np.zeros((A.shape[0], 0))
    U0 = unit_space(As[pts[0]])
    dev = max(np.linalg.norm(U - U0@(U0.T@U)) for U in (unit_space(A) for A in As.values()))
    # left eigenvectors define the invariant complement; reduce with a fixed basis
    W = null_space(U0.T) if U0.shape[1] else np.eye(U0.shape[0])   # orthonormal complement (identity if none)
    # projected dynamics are exact iff U0 is invariant for all A and the complement map is A-independent of the
    # unit modes; use the quotient system on R^n / span(U0): A_r = W' A W (valid since A U0 = U0)
    red = {pt: W.T@A@W for pt, A in As.items()}
    rho = {pt: float(np.max(np.abs(np.linalg.eigvals(Ar)))) for pt, Ar in red.items()}
    stable = [pt for pt in pts if rho[pt] < 1]
    Ps = {pt: solve_discrete_lyapunov(red[pt].T, np.eye(red[pt].shape[0])) for pt in stable}
    a = max(np.linalg.norm(red[pt], 2) for pt in stable)
    c2 = max(np.linalg.eigvalsh(Ps[pt]).max() for pt in stable); c1 = min(np.linalg.eigvalsh(Ps[pt]).min() for pt in stable)
    # Lipschitz constants along each axis (finite differences between neighbouring grid points)
    names = ('lkap', 'lkg', 'mu', 'lk'); Lc = {}
    for ax, nm in enumerate(names):
        g = grid[nm]; best = 0.
        for pt in stable:
            i = int(np.argmin(np.abs(g-pt[ax])))
            if i+1 < len(g):
                q = list(pt); q[ax] = g[i+1]; q = tuple(q)
                if q in Ps:
                    best = max(best, np.linalg.norm(Ps[q]-Ps[pt], 2)/abs(g[i+1]-g[i]))
        Lc[nm] = best
    budget = 1/a**2
    allow = {nm: budget/Lc[nm] if Lc[nm] > 0 else np.inf for nm in names}
    Nd = math.log(c2/c1)/(-math.log(1-1/c2))
    out = dict(unit_modes=int(U0.shape[1]), unit_space_deviation=float(dev), n_points=len(pts), n_stable=len(stable),
               unstable=[[float(math.exp(p[0])), float(math.exp(p[1])), float(p[2]), float(math.exp(p[3]))] for p in pts if rho[p] >= 1],
               rho_max=max(rho[pt] for pt in stable), a=a, c1=c1, c2=c2, lipschitz=Lc, budget=budget,
               allowable_rate_per_tick=allow, average_dwell_ticks=Nd)
    wide = [(math.log(2.0), b_, c_, d_) for b_ in grid['lkg'] for c_ in grid['mu'] for d_ in grid['lk']]
    out['kappa2_unstable'] = int(sum(np.max(np.abs(np.linalg.eigvals(Amat(*w_)))) >= 1 for w_ in wide))
    out['kappa2_points'] = len(wide)
    out['note'] = ('balanced-coordinate slow-variation certificate (scratchpad check) allows ~1e-7 per tick; '
                   'the raw-coordinate numbers below are not meaningful as rates (state units unbalanced)')
    json.dump(out, open(os.path.join(D.HERE, 'results_v2', 'theory_ltv.json'), 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != 'unstable'}, indent=1, default=float))
    print('unstable grid points (kappa, kappa_g, mu, k):', out['unstable'][:20], '... total', len(out['unstable']))
