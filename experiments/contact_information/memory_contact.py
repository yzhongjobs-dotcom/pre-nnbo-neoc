"""Local SLS transfer, with unknown moving surface and memory initial state.

F=(ke_end+ke_slope*t/T)*(y-s)+km*r, rdot=(ydot-sdot)-r/tau_m.
For cubic s(t), its contribution lies in B4 plus one exponential. This
profiles surface/initial memory without pretending absolute indentation is
measured. km is locally constant; ke may vary linearly in this short window.
The relaxation times below are explicit hypotheses, NOT identified truths.
For constrained fits, the covariance is the unconstrained design reference
inflated by residual error; it is NOT a calibrated post-constraint posterior.
"""
import math
import numpy as np
from local_contact import residualize

TIMES = dict(memory_040=.04, memory_200=.2, memory_600=.6,
             memoryc_040=.04, memoryc_200=.2, drift_040=.04)


def features(t, y0, y1, te, ye, tw, tau, tau_m, origin, q0=0.):
    dt = .008
    m = np.arange(8)*.001
    af = -math.expm1(-2*math.pi*30*.001)
    w = af*(1-af)**(7-np.arange(8))
    frac = -np.expm1(-m/tau)/(-math.expm1(-dt/tau))
    u = y0[:, None]+(y1-y0)[:, None]*frac-ye
    v = (t[:, None]+m-te)/tw
    B = np.column_stack([(v**p)@w for p in range(5)]+
                        [np.exp(-(t[:, None]+m-origin)/tau_m)@w])
    rows = []
    q = float(q0)
    den = 1/tau_m-1/tau
    for delta in y1-y0:
        speed = delta/(tau*(-math.expm1(-dt/tau)))
        if abs(den) < 1e-10:
            response = speed*m*np.exp(-m/tau)
            increment = speed*dt*math.exp(-dt/tau)
        else:
            response = speed*(np.exp(-m/tau)-np.exp(-m/tau_m))/den
            increment = speed*(math.exp(-dt/tau)-math.exp(-dt/tau_m))/den
        rows.append(q*np.exp(-m/tau_m)+response)
        q = q*math.exp(-dt/tau_m)+increment
    U = np.column_stack([u@w, (v*u)@w, np.asarray(rows)@w])
    return B, U, q


def fit(measured, model='memory_200', nw=32, tau=.019, noise=.04):
    if measured.ndim != 2 or measured.shape[1] != 7:
        raise ValueError('seven measured channels required')
    tau_m = TIMES[model]
    n = len(measured)
    j = np.arange(max(1, n-1-nw), n-1)
    if len(j) < 16 or np.any(measured[j, 3] <= .5) or np.any(measured[j+1, 3] <= .5):
        return None  # this local memory law assumes maintained contact
    te, ye, tw, origin = measured[-1, 0], measured[-1, 2], nw*.008, measured[j[0], 0]
    B, U, qend = features(measured[j, 0], measured[j, 2], measured[j+1, 2], te, ye, tw, tau, tau_m, origin)
    drift_only = model.startswith('drift_')
    if drift_only:
        U = U[:, :2]
    z = measured[j+1, 3]-math.exp(-2*math.pi*30*.008)*measured[j, 3]
    ur, zr = residualize(B, U), residualize(B, z)
    norm = np.linalg.norm(ur, axis=0)
    if np.min(norm) < 1e-13:
        return None
    _, s, vt = np.linalg.svd(ur/norm, full_matrices=False)
    if s[-1] < s[0]*1e-8:
        return None
    coef = np.linalg.lstsq(ur/norm, zr, rcond=1e-8)[0]/norm
    constrained = model.startswith('memoryc_')
    if constrained:
        from scipy.optimize import lsq_linear
        # Positive endpoint relaxed gains and nonnegative memory modulus.
        # [a,b,c] -> [ke_end,ke_slope,km]=[a,a-b,c].
        transform = np.array([[1., 0., 0.], [1., -1., 0.], [0., 0., 1.]])
        Ac = ur@transform
        cn = np.linalg.norm(Ac, axis=0)
        solution = lsq_linear(Ac/cn, zr, bounds=(0., np.inf), tol=1e-12)
        if not solution.success:
            return None
        coef = transform@(solution.x/cn)
    residual = zr-ur@coef
    dof = len(j)-B.shape[1]-U.shape[1]
    chi = float(residual@residual/(dof*noise**2))
    variance = noise**2*max(1., chi)
    cov = variance*((vt.T/s**2)@vt)/norm[:, None]/norm[None, :]
    beta = np.linalg.lstsq(B, z-U@coef, rcond=1e-12)[0]
    X = np.column_stack([B, U]); xn = np.linalg.norm(X, axis=0)
    _, ss, vv = np.linalg.svd(X/xn, full_matrices=False)
    full_cov = variance*((vv.T/ss**2)@vv)/xn[:, None]/xn[None, :]
    w = 2*np.pi*20.8
    functional = np.array([1., 0.] if drift_only else [1., 0., (w*tau_m)**2/(1+(w*tau_m)**2)])
    k = float(functional@coef)
    sd = float(np.sqrt(functional@cov@functional))
    physical = bool(coef[0] >= 0 and coef[0]-coef[1] >= -1e-8 and (drift_only or coef[2] >= 0))
    return dict(k=k, sigma=sd, coefficients=coef, cov=cov, beta=beta, full_cov=full_cov, variance=variance,
                te=float(te), ye=float(ye), tw=tw, origin=float(origin), qend=qend, tau_m=tau_m,
                risk=sd/max(k, 1.)*(1+math.sqrt(dof*max(chi-1., 0.))), chi2=chi, physical=physical,
                min_gain=float(min(coef[0], coef[0]-coef[1])) if physical else -1.,
                condition=float(s[0]/s[-1]), residual_rms=float(np.sqrt(np.mean(residual**2))),
                low_frequency_gain=float(coef[0]), low_frequency_sd=float(np.sqrt(cov[0, 0])),
                drift_only=drift_only,
                quadrature_gain=0. if drift_only else float(coef[2]*w*tau_m/(1+(w*tau_m)**2)))


def predict(r, measured, end, horizon, tau, model):
    j = np.arange(end-1, end-1+horizon)
    B, U, _ = features(measured[j, 0], measured[j, 2], measured[j+1, 2],
                       r['te'], r['ye'], r['tw'], tau, r['tau_m'], r['origin'], r['qend'])
    if r['drift_only']:
        U = U[:, :2]
    decay = math.exp(-2*math.pi*30*.008)
    weights = decay**np.arange(horizon-1, -1, -1)
    a = weights@np.column_stack([B, U])
    value = decay**horizon*measured[end-1, 3]+a@np.r_[r['beta'], r['coefficients']]
    var = r['variance']*float(weights@weights)+float(a@r['full_cov']@a)
    return float(value), float(max(var, 1e-15))
