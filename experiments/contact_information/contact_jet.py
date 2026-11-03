"""Measured-coordinate contact expansions. No plant, material or truth access.

These are local input/output derivatives, not identified elastic constants.
Products are formed at the 1-ms quadrature points BEFORE force filtering.
"""
import math
import numpy as np
from local_contact import residualize

MODELS = ('lean', 'motion3', 'normal_curvature', 'gain_curvature', 'centered_curvature')
SCALE = .001  # metres, numerical units only


def features(t, y0, y1, te, ye, tw, tau=.019, fc=30., model='lean', center=None):
    if model not in MODELS:
        raise ValueError(model)
    m = np.arange(8)*.001
    af = -math.expm1(-2*math.pi*fc*.001)
    w = af*(1-af)**(7-np.arange(8))
    frac = -np.expm1(-m/tau)/(-math.expm1(-.008/tau))
    u = y0[:, None]+(y1-y0)[:, None]*frac-ye
    v = (t[:, None]+m-te)/tw
    if center is not None:
        u = u-np.polynomial.polynomial.polyval(v, center)
    degree = 3 if model == 'motion3' else (5 if model == 'gain_curvature' else 4)
    B = np.column_stack([(v**p)@w for p in range(degree+1)])
    cols = [u@w, (v*u)@w]
    if model in ('normal_curvature', 'centered_curvature'):
        cols.append((.5*u*u/SCALE)@w)
    elif model == 'gain_curvature':
        cols.append((.5*v*v*u)@w)
    return B, np.column_stack(cols)


def fit(measured, nw=32, tau=.019, fc=30., noise=.04, model='lean'):
    if measured.ndim != 2 or measured.shape[1] != 7:
        raise ValueError('seven measured channels required')
    n = len(measured)
    j = np.arange(max(1, n-1-nw), n-1)
    j = j[(measured[j, 3] > .5) & (measured[j+1, 3] > .5)]
    if len(j) < 16:
        return None
    te, ye, tw = float(measured[-1, 0]), float(measured[-1, 2]), nw*.008
    center = None
    if model == 'centered_curvature':
        # Encoder-only moving reference. This is NOT an estimated tissue surface.
        ix = np.r_[j, j[-1]+1]
        tt = (measured[ix, 0]-te)/tw
        center = np.polynomial.polynomial.polyfit(tt, measured[ix, 2]-ye, 3)
    B, U = features(measured[j, 0], measured[j, 2], measured[j+1, 2], te, ye, tw, tau, fc, model, center)
    z = measured[j+1, 3]-math.exp(-2*math.pi*fc*.008)*measured[j, 3]
    ur, zr = residualize(B, U), residualize(B, z)
    norms = np.linalg.norm(ur, axis=0)
    if min(norms) < 1e-13:
        return None
    A = ur/norms
    _, s, vt = np.linalg.svd(A, full_matrices=False)
    if s[-1] < s[0]*1e-8:
        return None
    coef = np.linalg.lstsq(A, zr, rcond=1e-8)[0]/norms
    dof = max(len(z)-B.shape[1]-U.shape[1], 1)
    residual = zr-ur@coef
    chi = float(residual@residual/(dof*noise**2))
    variance = noise**2*max(1., chi)
    cov = variance*((vt.T/s**2)@vt)/norms[:, None]/norms[None, :]
    beta = np.linalg.lstsq(B, z-U@coef, rcond=1e-10)[0]
    # Full response covariance for held-out prediction includes drift uncertainty.
    X = np.column_stack([B, U])
    xn = np.linalg.norm(X, axis=0)
    _, ss, vv = np.linalg.svd(X/xn, full_matrices=False)
    full_cov = variance*((vv.T/ss**2)@vv)/xn[:, None]/xn[None, :]
    functional = np.zeros(len(coef)); functional[0] = 1.
    if center is not None:
        functional[2] = -center[0]/SCALE
    k, sd = float(functional@coef), float(np.sqrt(functional@cov@functional))
    risk = sd/max(k, 1.)*(1+math.sqrt(dof*max(chi-1., 0.)))
    # Positivity only on fitted coordinate box, NOT a passivity guarantee.
    v = np.linspace(-1., 0., 9)
    gain = coef[0]+coef[1]*v
    if model == 'gain_curvature':
        gain += .5*coef[2]*v*v
    if model in ('normal_curvature', 'centered_curvature'):
        displacement = measured[j, 2]-ye
        if center is not None:
            displacement = displacement-np.polynomial.polynomial.polyval((measured[j, 0]-te)/tw, center)
        min_gain = float(np.min(gain)+min(coef[2]*np.min(displacement)/SCALE,
                                         coef[2]*np.max(displacement)/SCALE))
    else:
        min_gain = float(np.min(gain))
    return dict(k=k, sigma=sd, relative_sd=sd/max(k, 1.), risk=risk, chi2=chi,
                coefficients=coef, cov=cov, beta=beta, full_cov=full_cov,
                variance=variance, te=te, ye=ye, tw=tw, dof=dof, min_gain=min_gain,
                center=center, functional=functional,
                condition=float(s[0]/s[-1]), residual_rms=float(np.sqrt(np.mean(residual**2))))


def predict(r, measured, end, horizon, tau, model):
    """Conditional future force: future encoders visible, future forces hidden.

    Multi-step propagation starts from F[end-1] only; no intermediate future
    force resets. This evaluates contact representation, not motion forecasts.
    """
    j = np.arange(end-1, end-1+horizon)
    B, U = features(measured[j, 0], measured[j, 2], measured[j+1, 2],
                    r['te'], r['ye'], r['tw'], tau, model=model, center=r['center'])
    decay = math.exp(-2*math.pi*30*.008)
    weights = decay**np.arange(horizon-1, -1, -1)
    a = weights@np.column_stack([B, U])
    value = decay**horizon*measured[end-1, 3]+a@np.r_[r['beta'], r['coefficients']]
    # Conditional fitted-model noise, not calibrated closed-loop confidence.
    var = r['variance']*float(weights@weights)+float(a@r['full_cov']@a)
    return float(value), float(max(var, 1e-15))
