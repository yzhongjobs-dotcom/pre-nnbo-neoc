"""Frozen IMM + unconstrained MPC stability screen before task tuning.

States, in kbar-scaled displacement units except physical force:
[y, last_sent, F, Fprev, yprev, x_s(3), x_m(3), model_filtered_F].
The frozen measurement gain uses the exact same R as FINAL. Model-conditioned
gains are frozen separately; IMM mixing and merged output retain frozen mu.
The accumulated measurement level is represented algebraically as in theory_ltv.
"""
import itertools
from functools import lru_cache
import numpy as np
from scipy.linalg import solve_discrete_are, cho_factor, cho_solve
from mpc_controller import (TC, TAU, AF, AS, AM, CY, CU, DEC, WEIGHTS,
                            prediction_matrices, CE)

GRID = dict(kappa=[.5, .7, 1., 1.4, 1.8, 2., 2.14],
            kappa_g=[.82, .91, 1., 1.09, 1.18],
            k=[1e3, 3e3, 1e4, 3e4, 1e5], mu=[0., .5, 1.])
AT = np.array([[1, TC, TC*TC/2], [0, 1, TC], [0, 0, 1.]])
Q1 = np.array([[TC**5/20, TC**4/8, TC**3/6],
               [TC**4/8, TC**3/3, TC**2/2], [TC**3/6, TC**2/2, TC]])
FRAC = float(WEIGHTS@((1-AM)/(1-AS))/WEIGHTS.sum())
LAG = (1-FRAC)*TC
SHIFT = np.array([[1, LAG, .5*LAG*LAG], [0, 1, LAG], [0, 0, 1.]])
NX = 12


@lru_cache(None)
def kgain(q, kref):
    R = CE.FINAL['r_scale']*(.04/kref)**2
    scale = np.diag([1., TC, TC*TC])
    a = scale@AT@np.linalg.inv(scale)
    cov = solve_discrete_are(a.T, np.array([[1.], [0.], [0.]]), scale@(q*Q1)@scale/R, np.ones((1, 1)))
    return np.linalg.solve(scale, cov[:, 0]/(cov[0, 0]+1))


@lru_cache(None)
def base_map(kap, kg, k, mu):
    Ks = kgain(CE.FINAL['q_smooth'], k/kg)
    Km = kgain(CE.FINAL['q_maneuver'], k/kg)
    pst = CE.FINAL['p_stay']
    def before(x):
        y, cmd, F, Fp, yp = x[:5]
        xs, xm = x[5:8], x[8:11]
        fin = (F-DEC*Fp)/(1-DEC)
        Z = kg/kap*fin-(yp+FRAC*(y-yp))
        ms, mm = 1-mu, mu
        cs, cm = pst*ms+(1-pst)*mm, (1-pst)*ms+pst*mm
        ps = AT@((pst*ms*xs+(1-pst)*mm*xm)/cs)
        pm = AT@(((1-pst)*ms*xs+pst*mm*xm)/cm)
        xs2, xm2 = ps+Ks*(Z-ps[0]), pm+Km*(Z-pm[0])
        eta = SHIFT@(ms*xs2+mm*xm2)
        r = kap/kg  # khat/kbar
        state = np.array([r*(cmd-y), x[11], r*(y+eta[0]), r*eta[1], r*eta[2]])
        bias = F-x[11]
        nx = np.r_[AS*y+(1-AS)*cmd, cmd, DEC*F+kap*(CY*y+CU*cmd), F, y, xs2, xm2, 0.]
        return nx, state, bias
    vals = [before(np.eye(NX)[j]) for j in range(NX)]
    A = np.array([v[0] for v in vals]).T
    S = np.array([v[1] for v in vals]).T
    B = np.array([v[2] for v in vals])
    return A, S, B


def screen(params, retain=False):
    H, qr, terminal = params['H'], params['qr'], params['terminal']
    BF, BR, GF, GR = prediction_matrices(H)
    weights = np.ones(H); weights[-1] = terminal
    gtg = GF.T@(weights[:, None]*GF)
    points = list(itertools.product(GRID['kappa'], GRID['kappa_g'], GRID['k'], GRID['mu']))
    mats, units = [], []
    for kap, kg, k, mu in points:
        r = kap/kg
        gain = cho_solve(cho_factor(gtg+np.eye(H)/(qr*r*r), check_finite=False), GF.T*weights, check_finite=False)[0]
        a, s, b = base_map(kap, kg, k, mu)
        a = a.copy()
        a[1] -= (gain@(BF@s+np.ones((H, 1))*b))/r
        a[11] = BF[0]@s
        mats.append(a)
    ev = np.linalg.eigvals(np.array(mats))
    # The algebraic measurement fixes its arbitrary integration origin, as in
    # ff_loop.py. Keep ALL eigenvalues of this represented dynamic system;
    # no near-unit pole may be discarded to make a candidate pass.
    unit = np.abs(ev-1.) < 1e-7
    rho = np.max(np.abs(ev), axis=1)
    worst = int(np.argmax(rho))
    result = dict(params=params, points=len(points), n_stable=int(np.sum(rho < 1-1e-8)),
                  rho_max=float(rho[worst]), worst=dict(zip(('kappa', 'kappa_g', 'k', 'mu'), points[worst])),
                  unit_count_range=[int(unit.sum(axis=1).min()), int(unit.sum(axis=1).max())],
                  accepted=bool(np.all(rho < 1-1e-8)))
    if retain:
        result['rows'] = [dict(kappa=p[0], kappa_g=p[1], k=p[2], mu=p[3], rho=float(r)) for p, r in zip(points, rho)]
    return result


if __name__ == '__main__':
    import json
    for qr in (.001, .01, .1, 1., 10., 100.):
        print(json.dumps(screen(dict(H=15, qr=qr, terminal=1.))), flush=True)
