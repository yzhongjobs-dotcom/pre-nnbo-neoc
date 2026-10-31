"""Shared pieces for the motion-aware contact study (paper draft §三–四, methods §10).

All policies here run their own controller and return an absolute position command
every control period (plant direct mode). They read only what a causal controller
may read through plant.Obs: time, filtered force, encoder position and sent commands.
Oracle policies additionally read truth columns and the surface function; they are
reference controllers, not deployable methods.
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import plant as P          # noqa: E402
import scenarios as S      # noqa: E402

FD = P.FD
K_MIN, K_MAX = 200., 3e5          # stiffness search range [N/m], as in C1


class DirectPolicy:
    """Base class: subclasses implement step(obs) -> absolute command [m]."""
    direct = True
    name = 'direct'

    def reset(self, platform):
        self.Tc = platform['Tc']; self.ndel = int(platform['ndel'])
        self.tau = platform['tau_a']; self.fc = platform['fc']; self.VMAX = platform['VMAX']
        self.ns = int(round(self.Tc/1e-3))                     # 1 kHz force samples per tick
        self.afc = 1.0-math.exp(-2*math.pi*self.fc*1e-3)       # force filter coefficient
        self.eb = math.exp(-1e-3/self.tau) if self.tau > 0 else 0.
        self.t_call = []; self.trace = []
        self.y_cmd = None

    def __call__(self, o):
        t0 = time.perf_counter()
        if o.n == 0:
            self.t_call.append(time.perf_counter()-t0)
            return float('nan')                              # plant holds its initial command
        if self.y_cmd is None:
            self.y_cmd = float(o.log[0, 11])
        y = self.step(o)
        self.y_cmd = y
        self.t_call.append(time.perf_counter()-t0)
        return y

    def step(self, o):
        raise NotImplementedError


from numba import njit


@njit(cache=True)
def _build(log, j0, m, degree, ns, eb, afc, ndel, te, ye, tw, enc):
    a = 1.-afc
    X = np.zeros((m, degree+1)); Y = np.zeros(m); Z = np.zeros(m)
    keep = np.zeros(m, np.bool_)
    for jj in range(m):
        j = j0+jj
        if min(log[j, 5], log[j+1, 5]) <= .5:
            continue
        keep[jj] = True
        yy = log[j, 2]
        tgt = log[j-ndel, 11]
        y1 = log[j+1, 2]
        ebn = eb**ns
        for mm in range(1, ns+1):
            # plant order per 1 ms sample: the force is sampled with the current position, then the
            # servo advances; sample mm-1 therefore sees the position after mm-1 ms of servo motion
            w = afc*a**(ns-mm)
            tau = (log[j, 0]+(mm-1)*1e-3-te)/tw
            if enc:
                # measured endpoints y_j, y_{j+1}; only the intra-period shape uses the servo time constant, so a
                # servo-lag or delay mismatch distorts Y far less than propagating the sent command through the model
                fr = (1.-eb**(mm-1))/(1.-ebn) if ebn < 1. else (mm-1)/ns
                Y[jj] += w*(log[j, 2]+fr*(y1-log[j, 2])-ye)
            else:
                Y[jj] += w*(yy-ye)
            for d in range(degree+1):
                X[jj, d] += w*tau**d
            yy = tgt+(yy-tgt)*eb
        Z[jj] = log[j+1, 5]-a**ns*log[j, 5]
    return X, Y, Z, keep


def window_fit(log, n, nw, degree, noise, ns, eb, afc, ndel, enc=False):
    """Equation-error contact regression on the last nw control intervals.

    z_j = F[j+1]-a^ns F[j] (one-step inverse of the 30 Hz force filter) is regressed on
    the filtered, servo-interpolated probe displacement Y_j and a degree-p time basis
    B_j filtered the same way (same processing domain), i.e. z = k Y + B beta + eps.
    The basis is eliminated first (variable projection); beta is then recovered from
    z - k Y. Returns dict(k, beta, info=||Y_perp||^2/noise^2, chi2, valid, tw, te).
    degree=0 is the static-surface estimator (constant offset only).
    """
    j0 = max(ndel+1, n-1-nw)
    m = n-2-j0+1
    out = dict(k=np.nan, beta=np.zeros(degree+1), info=0., chi2=np.inf, valid=0, tw=nw*ns*1e-3,
               te=log[n-1, 0] if n > 0 else 0.)
    if m < degree+3:
        return out
    te, ye = log[n-1, 0], log[n-1, 2]
    tw = nw*ns*1e-3
    X, Y, Z, keep = _build(log, j0, m, degree, ns, eb, afc, ndel, te, ye, tw, enc)
    X, Y, Z = X[keep], Y[keep], Z[keep]
    valid = int(keep.sum())
    out['valid'] = valid
    if valid < degree+3:
        return out
    gram = X.T@X+1e-12*np.eye(degree+1)
    cy = np.linalg.solve(gram, X.T@Y); cz = np.linalg.solve(gram, X.T@Z)
    yr = Y-X@cy; zr = Z-X@cz
    syy = float(yr@yr)
    if syy <= 0:
        return out
    k = float(yr@zr)/syy
    beta = cz-k*cy
    res = Z-k*Y-X@beta
    dof = max(valid-degree-2, 1)
    out.update(k=k, beta=beta, info=syy/noise**2, chi2=float(res@res)/noise**2/dof)
    return out


def drift_poly(fit, h):
    """equivalent force drift and its time derivative at te+h from a window fit."""
    b, tw = fit['beta'], fit['tw']
    x = h/tw
    g = sum(b[d]*x**d for d in range(len(b)))
    dg = sum(d*b[d]*x**(d-1) for d in range(1, len(b)))/tw
    return g, dg


def run(policy, p, T, seed, platform='PC125'):
    """one episode; returns (log, policy)."""
    log = P.episode(p, policy, T=T, platform=platform, seed=seed)
    return log, policy


def task_metrics(log, f_exc=None, warm=1.0, Tc=0.008):
    """physical metrics on the true contact force; oscillation excludes the
    excitation line so that deliberate dither is reported separately as ripple."""
    m = S.metrics(log, warm=warm, Tc=Tc)
    if log is None or log.shape[0] < 10:
        m.update(osc_ex=np.nan, ripple=np.nan, n_loss=np.nan)
        return m
    sel = log[:, 0] >= warm
    f = log[sel, 3]; t = log[sel, 0]
    w = max(int(round(0.1/Tc)), 1)
    hp = f-np.convolve(f, np.ones(w)/w, mode='same')
    hp = hp[w:-w]; tt = t[w:-w]
    ripple = 0.
    if f_exc:
        c, s_ = np.cos(2*np.pi*f_exc*tt), np.sin(2*np.pi*f_exc*tt)
        A = np.vstack([c, s_]).T
        coef, *_ = np.linalg.lstsq(A, hp, rcond=None)
        hp = hp-A@coef
        ripple = float(np.hypot(*coef)/math.sqrt(2))
    loss = f < S.FMIN
    n_loss = int(np.sum(loss[1:] & ~loss[:-1])+(1 if loss[0] else 0))
    m.update(osc_ex=float(np.sqrt(np.mean(hp**2))), ripple=ripple, n_loss=n_loss)
    return m
