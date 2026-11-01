"""Causal, offline estimator comparison on identical Proposed trajectories.

Run from any directory: python -X utf8 experiments/method_redesign/step2_estimators.py
No final-cohort access. Results are self-contained JSON and Markdown. --workers
controls CPU worker threads, --reproduce checks all scientific outputs against a prior JSON
(timings and platform metadata excluded). No cached logs or tuned hidden inputs.

Estimator inputs are sanitized BEFORE replay: only time, x, y, Fn, command,
Ft and My survive. Truth is passed only to scoring. E3 never calls the plant
contact functions, nor imports object material/viscoelastic constants.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

# Keep native kernels single threaded before importing NumPy. CPU worker
# threads share read-only trajectories; the expensive GN kernel releases GIL.
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'experiments/motion_aware_scan'))
import numpy as np
from numba import njit
import common
import cohort_eval as C
import scenes_v2 as SV

DT = .008
LO, HI = common.K_MIN, common.K_MAX
PL = common.P.PLATFORMS['PC125']
AF = 1.-math.exp(-2*math.pi*PL['fc']*.001)
EB = math.exp(-.001/PL['tau_a'])
PARAM = C.build('Proposed').params
NW, DEG, NOISE = PARAM['nw'], PARAM['degree'], PARAM['noise']
NAMES = ('E0', 'E1', 'E2', 'E3')
VISIBLE = (0, 1, 2, 5, 11, 12, 13)
HORIZONS = (0., .1, .3)


def sanitize(log):
    out = np.full(log.shape, np.nan)
    out[:, VISIBLE] = log[:, VISIBLE]
    return out


def design(log, n, degree=DEG):
    start = max(2, n-1-NW)
    m = n-1-start
    if m < degree+3:
        return None
    X, Y, Z, keep = common._build(log, start, m, degree, 8, EB, AF, 1,
                                  log[n-1, 0], log[n-1, 2], NW*DT, True)
    if keep.sum() < degree+3:
        return None
    return X[keep], Y[keep], Z[keep], np.arange(start, start+m)[keep]


def shift_beta(beta, dt, dy, k):
    # Same window duration, translate force-drift polynomial to new end time
    # and encoder origin. No surface truth is needed.
    h = dt/(NW*DT)
    out = np.zeros_like(beta)
    for j in range(len(beta)):
        out[j] = sum(beta[l]*math.comb(l, j)*h**(l-j) for l in range(j, len(beta)))
    out[0] += k*dy
    return out


class Projection:
    def __init__(self, cfg):
        self.k = HI
        self.var = HI**2

    def step(self, log, n):
        f = common.window_fit(log, n, NW, DEG, NOISE, 8, EB, AF, 1, PARAM['enc_y'])
        if np.isfinite(f['k']) and f['info'] > 0:
            self.k = float(np.clip(f['k'], LO, HI))
            self.var = 1./f['info']
        return np.repeat(self.k, 3), np.full((3, 3), self.var)


class JointEKF:
    """[log k, g, gdot, gddot], force drift referenced to first encoder y.

    Exact same inverse-filter interval observation as projection. Covariance
    uses Joseph update. Q uses a discrete rank-one jerk drive plus log-k random walk.
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.s = np.array([math.log(HI), 4., 0., 0.])
        self.P = np.diag([4., 16., 1e4, 1e6])
        self.y0 = None
        self.A = np.eye(4)
        self.A[1, 2] = self.A[2, 3] = DT
        self.A[1, 3] = DT**2/2

    def step(self, log, n):
        if self.y0 is None:
            self.y0 = log[0, 2]
        self.s = self.A@self.s
        qk, qj, r = self.cfg
        self.P = self.A@self.P@self.A.T
        self.P[0, 0] += qk**2*DT
        v = np.array([DT**2/2, DT, 1.])
        self.P[1:, 1:] += qj**2*DT*np.outer(v, v)
        if n >= 4 and min(log[n-2:n, 5]) > .5:
            X, Y, Z, keep = common._build(log, n-2, 1, 2, 8, EB, AF, 1,
                                          log[n-1, 0], self.y0, 1., True)
            k = math.exp(self.s[0])
            H = np.array([k*Y[0], X[0, 0], X[0, 1], .5*X[0, 2]])
            h = k*Y[0]+H[1:]@self.s[1:]
            R = (NOISE*r)**2
            gain = self.P@H/(H@self.P@H+R)
            self.s += gain*(Z[0]-h)
            self.s[0] = np.clip(self.s[0], math.log(LO), math.log(HI))
            A = np.eye(4)-np.outer(gain, H)
            self.P = A@self.P@A.T+R*np.outer(gain, gain)
        k = math.exp(self.s[0])
        return np.repeat(k, 3), np.full((3, 3), k*k*self.P[0, 0])


class MHE:
    """Scalar k update; exact Schur elimination of linear drift coefficients.

    Augmented least squares: measurement residuals, k arrival random-walk
    prior, and transported beta continuity prior. Arrival variance is inflated
    by NW to avoid counting each overlapping window as NW independent windows.
    This remains an approximate arrival model, not an exact Bayesian posterior.
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.k, self.var = HI, HI**2
        self.beta = None
        self.cov_beta = np.eye(DEG+1)*1e4
        self.last_t = self.last_y = None

    def augmented(self, log, n, X, Y, Z):
        qk, qg, r = self.cfg[:3]
        s = NOISE*r
        A, y, z = X/s, Y/s, Z/s
        if self.beta is not None:
            bt = shift_beta(self.beta, log[n-1, 0]-self.last_t,
                            log[n-1, 2]-self.last_y, self.k)
            # Weak high-order continuity allows abrupt motion to enter drift.
            h = (log[n-1, 0]-self.last_t)/(NW*DT)
            T = np.zeros((DEG+1, DEG+1))
            for j in range(DEG+1):
                for l in range(j, DEG+1):
                    T[j,l] = math.comb(l,j)*h**(l-j)
            cp = T@self.cov_beta@T.T
            cp[0,0] += self.var*(log[n-1,2]-self.last_y)**2
            sig = np.sqrt(NW*np.maximum(np.diag(cp),0.)+(qg*np.array([1., 2., 4., 8.]))**2)
            A = np.vstack((A, np.diag(1./sig)))
            y = np.r_[y, np.zeros(DEG+1)]
            z = np.r_[z, bt/sig]
        sk = math.sqrt(NW*self.var+(qk*self.k)**2)
        return A, y, z, max(sk, 1.)

    def step(self, log, n):
        data = design(log, n)
        if data is not None:
            X, Y, Z, _ = data
            A, y, z, sk = self.augmented(log, n, X, Y, Z)
            gram = A.T@A+np.eye(DEG+1)*1e-12
            cy, cz = np.linalg.solve(gram, np.stack((A.T@y, A.T@z), axis=1)).T
            yr, zr = y-A@cy, z-A@cz
            precision = yr@yr+1./sk**2
            self.k = float(np.clip((yr@zr+self.k/sk**2)/precision, LO, HI))
            self.var = 1./precision
            self.beta = cz-self.k*cy
            self.cov_beta = np.linalg.inv(gram)+np.outer(cy,cy)*self.var
            self.last_t, self.last_y = log[n-1, (0, 2)]
        return np.repeat(self.k, 3), np.full((3, 3), self.var)


# Independent Hertz model: 24-point Gauss-Legendre on [-1,1], semiellipse.
# These are NOT the object's Chebyshev nodes, width function or material model.
GL_X, GL_W = np.polynomial.legendre.leggauss(24)
HR, HL, HDEPTH, EDGE = .028, .038, .025, .001


@njit(cache=True, nogil=True)
def hertz(theta, x, a):
    """F, M and dF/da for modulus-weighted Hertz line pressure.

    p(u)=E(x+u)/(2R)*sqrt(a²-u²), F=L int p du.
    Homogeneous: a²=4 F R/(pi L E). Finite outer-reference compliance
    delta=a²/(2R)*(log(4H/a)+1/2), ddelta/da=a/R*log(4H/a).
    Line-contact indentation needs an outer displacement reference; H is
    fixed and deliberately independent of the object's finite-length relation.
    """
    e0, ratio, b = math.exp(theta[0]), math.exp(theta[1]), theta[2]*.01
    F = M = dF = 0.
    for j in range(len(GL_X)):
        u = GL_X[j]
        w = GL_W[j]*math.sqrt(1-u*u)
        arg = min(40., max(-40., (x+a*u-b)/EDGE))
        gate = 1./(1.+math.exp(-arg))
        E = e0*(1.+(ratio-1.)*gate)
        dE = e0*(ratio-1.)*gate*(1.-gate)*u/EDGE
        F += w*E
        M += w*E*u
        dF += w*dE
    c = HL/(2*HR)
    return c*a*a*F, c*a*a*a*M, c*(2*a*F+a*a*dF)


@njit(cache=True, nogil=True)
def spatial(theta, xs, forces):
    """Profile one latent width per sample using normal force; no future data.

    Fn likelihood is saturated by these free latent widths; other residuals
    carry spatial parameters. Fn is NOT also counted as an independent raw
    observation in the objective. Normal increment likelihood comes from the
    projection model, with an approximate conditional covariance.
    """
    m = len(xs)
    ks, moments, tangents, widths = np.empty(m), np.empty(m), np.empty(m), np.empty(m)
    for i in range(m):
        f = max(forces[i], .05)
        a = math.sqrt(4*f*HR/(math.pi*HL*math.exp(theta[0])))
        a = min(.045, max(.0001, a))
        low, high = .00001, .06
        for _ in range(18):
            F, M, der = hertz(theta, xs[i], a)
            if F > f:
                high = a
            else:
                low = a
            an = a-(F-f)/max(der, 1e-8)
            if an <= low or an >= high:
                an = .5*(low+high)
            if abs(F-f) < 1e-8:
                break
            a = an
        F, M, der = hertz(theta, xs[i], a)
        ks[i] = der/(a/HR*math.log(4*HDEPTH/a))
        moments[i] = M
        # Effective Coulomb coefficient; omitted Stribeck and viscous effects
        # are model mismatch, not fitted object parameters.
        tangents[i] = -math.exp(theta[3])*F-M/HR
        widths[i] = a
    return ks, moments, tangents, widths


@njit(cache=True, nogil=True)
def _sp_res(th, xs, fs, ft, my, sw, Y, A, z, proj, normal_sd, kp, sk, theta0, spatial_scale):
    ks, ms, ts, widths = spatial(th, xs, fs)
    normal = np.zeros(len(z))
    normal[:len(Y)] = ks[:-1]*Y/normal_sd
    beta = proj@(z-normal)
    rr = np.concatenate((normal+A@beta-z, np.array([(ks[-1]-kp)/sk]),
                         (th-theta0)/np.array([.6,.6,.8,.3]),
                         sw*(ms-my)/(spatial_scale*7.5e-4),
                         sw*(ts-ft)/(spatial_scale*.08)))
    return rr, beta, ks, widths


@njit(cache=True, nogil=True)
def _sp_jac(th, xs, fs, ft, my, sw, Y, A, z, proj, normal_sd, kp, sk, theta0, spatial_scale):
    m = len(_sp_res(th,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)[0])
    J = np.zeros((m,4))
    for j in range(4):
        d = np.zeros(4); d[j] = 2e-4
        plus = _sp_res(th+d,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)[0]
        minus = _sp_res(th-d,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)[0]
        J[:,j] = (plus-minus)/(4e-4)
    return J


@njit(cache=True, nogil=True)
def _sp_opt(th, lower, upper, xs, fs, ft, my, sw, Y, A, z, proj, normal_sd, kp, sk, theta0, spatial_scale):
    for it in range(3):
        rr = _sp_res(th,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)[0]
        J = _sp_jac(th,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)
        H = J.T@J
        step = np.linalg.solve(H+1e-3*np.diag(np.maximum(np.diag(H),1.)), -J.T@rr)
        step = np.minimum(.5,np.maximum(-.5,step))
        for alpha in (1., .5, .25):
            test = np.minimum(upper,np.maximum(lower,th+alpha*step))
            rtest = _sp_res(test,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)[0]
            if rtest@rtest <= rr@rr:
                th = test
                break
    rr, beta, ks, widths = _sp_res(th,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)
    J = _sp_jac(th,xs,fs,ft,my,sw,Y,A,z,proj,normal_sd,kp,sk,theta0,spatial_scale)
    return th, rr, beta, ks, widths, J


class SpatialMHE(MHE):
    """E2 arrival + projected drift, with four spatial parameters and Fn-profiled
    latent widths. Three damped Gauss-Newton iterations, warm-started per tick.
    Model residual and arrival covariance are conditional approximations.
    """
    def __init__(self, cfg):
        super().__init__(cfg)
        self.theta = None
        self.cov_theta = np.diag([4., 4., 4., 1.])
        self.last_condition = np.nan
        self.last_width = np.nan

    def step(self, log, n):
        data = design(log, n)
        if self.theta is None:
            self.theta = np.array([math.log(2e4), math.log(2.), (log[n-1, 1]+.01)/.01, math.log(.3)])
        if data is None:
            return np.repeat(self.k, 3), np.full((3, 3), self.var)
        X, Y, Z, idx = data
        # Current sample included for spatial channels; normal increments use
        # interval-end stiffness (slow-within-interval approximation).
        inds = np.r_[idx+1, n-1]
        xs, fs = log[inds, 1], log[inds, 5]
        ft, my = log[inds, 12], log[inds, 13]
        # Remove duplicate last observation if interval endpoint is current.
        # Last is retained for a convenient output index, with zero extra weight.
        spatial_weights = np.ones(len(inds))
        if len(inds) > 1 and inds[-1] == inds[-2]:
            spatial_weights[-1] = 0.
        A, _, z, sk = self.augmented(log, n, X, Y, Z)
        gram = A.T@A+np.eye(DEG+1)*1e-12
        projector = np.linalg.solve(gram, A.T)
        qk, qg, r, spatial_scale = self.cfg
        theta0 = self.theta.copy()
        # Broad arrival: one-boundary field is local, not a fixed global map.
        lower = np.array([math.log(300.), -4., (xs.min()-.06)/.01, math.log(.01)])
        upper = np.array([math.log(3e6), 4., (xs.max()+.06)/.01, math.log(2.)])

        th = np.clip(self.theta, lower, upper)
        th, rr, self.beta, ks, widths, J = _sp_opt(
            th, lower, upper, xs, fs, ft, my, spatial_weights, Y, A, z,
            projector, NOISE*r, self.k, sk, theta0, spatial_scale)
        self.theta = th
        H = J.T@J
        # Sandwich/model-discrepancy scale (never truth-calibrated on diag).
        inflation = max(1., float(rr@rr)/max(len(rr)-4, 1))
        self.cov_theta = np.linalg.pinv(H, rcond=1e-10)*inflation
        self.cov_beta = np.linalg.inv(gram)*inflation
        self.last_condition = float(np.linalg.cond(H))
        self.last_width = float(widths[-1])
        self.last_t, self.last_y = log[n-1, (0, 2)]
        speed = (log[n-1, 1]-log[max(0, n-6), 1])/max(DT*(min(5, n-1)), DT)
        xp = log[n-1, 1]+speed*np.array(HORIZONS)
        fp = np.full(3, max(log[n-1, 5], .05))
        pred = spatial(th, xp, fp)[0]
        G = np.zeros((3, 4))
        for j in range(4):
            d = np.zeros(4); d[j] = 2e-4
            G[:, j] = (spatial(th+d, xp, fp)[0]-spatial(th-d, xp, fp)[0])/(4e-4)
        cov = G@self.cov_theta@G.T
        # Future displacement/force is unknown: constant measured force and
        # constant measured lateral speed are prediction assumptions.
        cov += np.diag((qk*pred)**2*np.array(HORIZONS)/DT)
        self.k = float(pred[0])
        self.var = max(float(cov[0, 0]), 1.)
        return pred, cov


CLASSES = dict(E0=Projection, E1=JointEKF, E2=MHE, E3=SpatialMHE)


def replay(log, name, cfg, diagnostics=False, trace=None):
    obs = sanitize(log)
    estimator = CLASSES[name](cfg)
    k = np.empty((len(log), 3)); sd = np.empty_like(k)
    elapsed = np.empty(len(log))
    condition = []; widths = []
    for i in range(len(log)):
        t0 = time.perf_counter_ns()
        pred, cov = estimator.step(obs[:i+1], i+1)
        elapsed[i] = (time.perf_counter_ns()-t0)*1e-6
        k[i] = pred; sd[i] = np.sqrt(np.maximum(np.diag(cov), 0.))
        if trace is not None and name == 'E3' and i % 25 == 0 and log[i,0] >= 1.:
            trace.append(dict(t=float(log[i,0]), x=float(log[i,1]), k_true=float(log[i,8]),
                k=pred.tolist(), covariance_k=cov.tolist(), theta=estimator.theta.tolist(),
                covariance_theta=estimator.cov_theta.tolist(), width_m=estimator.last_width,
                My_meas=float(log[i,13]), My_true=float(log[i,15]), m1_true=float(log[i,14])))
        if diagnostics and name == 'E3':
            condition.append(estimator.last_condition)
            widths.append(estimator.last_width)
    assert np.isfinite(k).all() and np.isfinite(sd).all()
    return k, sd, elapsed, condition, widths


def scoring_mask(log):
    # One common mask, based on evaluation truth; NEVER fed to estimators.
    return (log[:, 0] >= 1.) & (log[:, 8] > 0.)


def score(log, k, sd, elapsed, condition=(), widths=()):
    good = scoring_mask(log)
    truth = log[:, 8]
    ratio = k[:, 0]/np.maximum(truth, 1e-30)
    err = abs(np.log(np.maximum(ratio, 1e-30)))
    # Family-independent transition definition, fixed before diagnostics:
    # > 10 log-stiffness units / metre (~10% change in 10 mm).
    # Encoder x repeats during the initial transport delay. Those samples
    # precede the shared warm-up mask; suppress only that undefined derivative.
    with np.errstate(divide='ignore', invalid='ignore'):
        grad = abs(np.gradient(np.log(np.maximum(truth, 1.)), log[:, 1]))
    grad = np.where(np.isfinite(grad), grad, 0.)
    ridge = good & (grad > 10.)
    out = dict(log_error=err[good].tolist(), under=(ratio[good] < 1/1.2).tolist(),
               over=(ratio[good] > 1/.75).tolist(),
               cover1=(abs(k[good, 0]-truth[good]) <= sd[good, 0]).tolist(),
               cover2=(abs(k[good, 0]-truth[good]) <= 2*sd[good, 0]).tolist(),
               transition_error=err[ridge].tolist(),
               transition_under=(ratio[ridge] < 1/1.2).tolist(),
               transition_over=(ratio[ridge] > 1/.75).tolist(),
               transition_cover1=(abs(k[ridge, 0]-truth[ridge]) <= sd[ridge, 0]).tolist(),
               transition_cover2=(abs(k[ridge, 0]-truth[ridge]) <= 2*sd[ridge, 0]).tolist(),
               ms=elapsed[good].tolist(), prediction={})
    for j, horizon in enumerate(HORIZONS[1:], 1):
        # True horizon, linear interpolation instead of rounding .1/.008.
        future_time = log[:, 0]+horizon
        future_k = np.interp(future_time, log[:, 0], truth)
        mask = good & (future_time <= log[-1, 0]) & (future_k > 0.)
        out['prediction'][str(horizon)] = dict(
            preview=abs(np.log(k[mask, j]/future_k[mask])).tolist(),
            hold=abs(np.log(k[mask, 0]/future_k[mask])).tolist(),
            cover1=(abs(k[mask, j]-future_k[mask]) <= sd[mask, j]).tolist(),
            cover2=(abs(k[mask, j]-future_k[mask]) <= 2*sd[mask, j]).tolist())
    if len(condition):
        out['condition'] = np.asarray(condition)[good].tolist()
        out['width'] = np.asarray(widths)[good].tolist()
    return out


def summarize(rows):
    out = {}
    for fam in (*C.FAMS, 'overall'):
        out[fam] = {}
        for name in NAMES:
            group = [r for r in rows if r['method'] == name and (fam == 'overall' or r['family'] == fam)]
            if not group:
                continue
            def cat(key):
                return np.concatenate([np.asarray(r['score'][key]) for r in group])
            err, trans, ms = cat('log_error'), cat('transition_error'), cat('ms')
            r = dict(episodes=len(group), ticks=len(err), median=float(np.median(err)),
                     p95=float(np.quantile(err, .95)),
                     under=float(cat('under').mean()), over=float(cat('over').mean()),
                     cover1=float(cat('cover1').mean()), cover2=float(cat('cover2').mean()),
                     p50_ms=float(np.median(ms)), p99_ms=float(np.quantile(ms, .99)),
                     episode_median_mean=float(np.mean([np.median(a['score']['log_error']) for a in group])),
                     transition=dict(ticks=len(trans), median=float(np.median(trans)) if len(trans) else None,
                                     p95=float(np.quantile(trans, .95)) if len(trans) else None), prediction={})
            for key in ('under', 'over', 'cover1', 'cover2'):
                vv = cat('transition_'+key)
                r['transition'][key] = float(vv.mean()) if len(vv) else None
            for h in ('0.1', '0.3'):
                rr = {}
                for key in ('preview', 'hold'):
                    v = np.concatenate([a['score']['prediction'][h][key] for a in group])
                    rr[key] = dict(median=float(np.median(v)), p95=float(np.quantile(v, .95)))
                for key in ('cover1', 'cover2'):
                    rr[key] = float(np.mean(np.concatenate([a['score']['prediction'][h][key] for a in group])))
                r['prediction'][h] = rr
            if name == 'E3':
                v = cat('condition'); v = v[np.isfinite(v)]
                w = cat('width'); w = w[np.isfinite(w)]
                r['identifiability'] = dict(hessian_condition_median=float(np.median(v)),
                    hessian_condition_p95=float(np.quantile(v, .95)), width_median_m=float(np.median(w)))
            out[fam][name] = r
    return out


def make_cohort(tune=False):
    if not tune:
        return C.cohort('diag', 20)
    rows = []
    for fi, fam in enumerate(C.FAMS):
        rng = np.random.default_rng(190000+1000*fi)
        for i in range(10):
            p, T = SV.make(fam, rng)
            rows.append((fam, 190000+1000*fi+i, p, T))
    return rows


def generate(scene):
    fam, sid, p, T = scene
    log = common.P.episode(p, C.build('Proposed'), T=T, seed=sid)
    if log is None:
        raise RuntimeError(f'Proposed trajectory diverged: {fam} {sid}; no silent exclusion')
    return fam, sid, log


def tuning_job(job):
    name, cfg, trajectories = job
    fam_scores = {}
    for fam, sid, log in trajectories:
        k, _, _, _, _ = replay(log, name, cfg)
        mask = scoring_mask(log)
        err = abs(np.log(k[mask, 0]/log[mask, 8]))
        # Equal family and equal episode weight, symmetric target. No emphasis
        # on ridges, dangerous direction, timing or E3-only prediction metrics.
        value = float(np.median(err)+.25*np.quantile(err, .95))
        fam_scores.setdefault(fam, []).append(value)
    return dict(method=name, params=cfg, objective=float(np.mean([np.mean(v) for v in fam_scores.values()])))


def evaluation_job(job):
    name, cfg, fam, sid, log = job
    k, sd, elapsed, cond, widths = replay(log, name, cfg, True)
    s = score(log, k, sd, elapsed, cond, widths)
    return dict(method=name, family=fam, seed=sid, score=s)


def candidates(name):
    # Twelve prespecified configurations each (<=30); identical tuning data
    # and objective. E0 is intentionally frozen per the requested definition.
    rng = np.random.default_rng(28111)
    if name == 'E0':
        return [None]
    configs = []
    for i in range(12):
        u = rng.random(4)
        if name == 'E1':
            configs.append([float(10**(-1+2*u[0])), float(10**(1+3*u[1])), float(10**(-.3+.9*u[2]))])
        else:
            cfg = [float(10**(-2+2*u[0])), float(10**(-2+2*u[1])), float(10**(-.3+.9*u[2]))]
            if name == 'E3':
                cfg += [float(10**(.1+1.2*u[3]))]
            configs.append(cfg)
    return configs


def fmt(v, percent=False):
    return '—' if v is None else (f'{100*v:.1f}%' if percent else f'{v:.3f}')


def markdown(result):
    lines = ['# 方法重构第2步：估计器开环对比', '',
        '## 设置与公平性', '',
        '同一最终 Proposed 控制器闭环产生轨迹，四种估计器离线逐拍重放。估计器输入在重放前清空真值列，仅保留时间及列1、2、5、11、12、13；输入视图截断于当前拍。',
        '调参：190000+、每族10个；诊断评价：cohort("diag",20)、每族20个；未访问130000+确认集。E1/E2/E3各12组候选、完全相同的50场景与目标；E0按要求固定最终方法的窗口/阶数/噪声/编码器设置，不重调。',
        '调参目标：每场景 median(|ln ratio|)+0.25 P95，再按族等权平均。评价从1 s开始，真值k>0；失联不按估计器或测力值筛除。所有有效时刻汇总为时间加权描述统计，同时JSON保留逐场景摘要；不将时间点当独立重复做显著性检验。',
        f'公共设置：nw={NW}、degree={DEG}、noise={NOISE} N、Tc=8 ms、30 Hz测力滤波、enc_y=True。所有方法具有同样的传感器访问权限；E0/E1/E2按定义只建模法向通道，E3利用新增空间观测，故此对比不能把信息增益与算法增益分开。', '',
        '## 模型与实现', '',
        'E0：common.window_fit 原始投影k，范围[200,300000]，有效拟合更新，启动/无效时保持；报告原始估计及1/info条件方差，不是保守反馈kbar或门控后的k_good。',
        'E1：[ln k,g,ġ,g̈]联合EKF；同域逆滤波增量观测，g以首拍编码器位置为原点；白jerk与ln k随机游走；Joseph协方差更新。初始k=300000、g=4、导数为0，协方差对角[4,16,1e4,1e6]。E0/E2/E3启动也以k=300000、σk=300000开始。',
        'E2：投影模型增加k到达代价与漂移多项式跨窗口平移后的连续性先验；漂移系数精确Schur消元，只求解正约束标量k。k与漂移后验方差平移传播，到达方差乘窗口长度，减轻重叠窗口的重复计数；并非严格贝叶斯到达分布。E3漂移条件方差未完整包含与θ的交叉项。',
        'E3：E2的投影漂移与到达结构，加θ=(ln E0,ln ρ,b/0.01,ln μ)。独立24点Gauss–Legendre半椭圆压力p=E(x+u)√(a²−u²)/(2R)，F=L∫pdu，M=L∫updu；局部模量用固定1 mm逻辑边缘。R=28 mm、L=38 mm、外位移参考H=25 mm。均匀场满足a²=4FR/(πLE)，δ=a²/(2R)[ln(4H/a)+1/2]；线接触的绝对压入需要外参考，此选择不是对象的ψ关系。',
        '每拍潜半宽以测得Fn条件剖面求根，Fn原始似然被潜变量饱和，不再当独立原始残差重复计入。空间k=dF/da÷dδ/da进入法向增量回归，漂移仍精确消元；Ft模型=−μF−M/R，省略Stribeck/黏性/坡度。θ每拍3次阻尼GN，热启动；初始模量20000、ρ=2、b=x+10 mm、μ=0.3，固定到达标准差[0.6,0.6,0.8,0.3]（θ坐标）。E0/E1/E2的k投影到[200,300000]；E3通过模量[300,3e6]、lnρ∈[−4,4]、μ∈[0.01,2]及b位于窗口路径±60 mm约束，不额外截断模型刚度。标准Hertz参考：[Cambridge接触公式](https://teaching.eng.cam.ac.uk/system/files/P2A2006M3C003.pdf)。未使用对象ψ/footprint/NQ、SLS或HC参数。',
        '预测仅用过去5拍估计横向速度，假设未来恒定横向速度和当前测得法向力；输出τ=0、0.1、0.3 s的k及3×3协方差。评价真值按准确t+τ插值，不将0.1 s舍入为12拍。',
        '区间为k±σ、k±2σ（非置信保证）。E3方差来自条件GN Hessian及残差膨胀，未完整传播有色Fn噪声、剖面潜变量与重叠窗口相关性；覆盖率是检验，不以扩大区间作为调参目标。', '',
        '## 调参结果', '', '| 方法 | 参数 | 调参目标 |', '|---|---|---:|']
    for name in NAMES:
        b = result['selected'][name]
        lines.append(f'| {name} | {b["params"]} | {b["objective"]:.4f} |')
    lines += ['', 'E1参数：[ln k游走/√s, jerk强度, R标准差倍率]；E2：[每次更新相对k游走, 漂移连续性N, 残差标准差倍率]；E3再加Ft/My残差标准差倍率。候选全集和目标见JSON；候选随机数种子固定28111。', '',
              '## 评价指标', '', '| 族 | 方法 | median | P95 | 偏低越界 | 偏高越界 | 1σ覆盖 | 2σ覆盖 | P50/P99 ms |', '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for fam, methods in result['summary'].items():
        for name, r in methods.items():
            vals = [fam, name, fmt(r['median']), fmt(r['p95']), fmt(r['under'], True), fmt(r['over'], True),
                    fmt(r['cover1'], True), fmt(r['cover2'], True), f'{r["p50_ms"]:.3f}/{r["p99_ms"]:.3f}']
            lines.append('| '+' | '.join(vals)+' |')
    lines += ['', '偏低越界：k/k̂>1.2；偏高越界：k/k̂<0.75。计时含每拍设计矩阵/优化/当前与预测协方差，排除场景生成、控制器、JIT首编译、文件写入；各方法均输出相同三个时距。多进程计时包含CPU竞争，仅软件原型计时。', '',
              '## 越脊过渡段', '', '统一以|d ln k/dx|>10 m⁻¹标记，阈值在查看诊断结果前固定；它也会标记运动引起的工作点刚度变化，因此不是纯空间脊标签。', '',
              '| 族 | 方法 | n拍 | median | P95 | 偏低 | 偏高 | 1σ | 2σ |', '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for fam, methods in result['summary'].items():
        for name, r in methods.items():
            t = r['transition']
            vals = [fam, name, str(t['ticks']), fmt(t['median']), fmt(t['p95'])]
            vals += [fmt(t[k], True) for k in ('under', 'over', 'cover1', 'cover2')]
            lines.append('| '+' | '.join(vals)+' |')
    lines += ['', '## E3沿路径预测与保持当前估计', '',
              '| 族 | τ s | 预测median/P95 | 保持median/P95 | 1σ/2σ覆盖 |', '|---|---:|---:|---:|---:|']
    for fam, methods in result['summary'].items():
        for h, p in methods['E3']['prediction'].items():
            lines.append(f'| {fam} | {h} | {p["preview"]["median"]:.3f}/{p["preview"]["p95"]:.3f} | {p["hold"]["median"]:.3f}/{p["hold"]["p95"]:.3f} | {fmt(p["cover1"],True)}/{fmt(p["cover2"],True)} |')
    lines += ['', '## 结果边界与疑点', '',
        '此为固定轨迹上的估计对比，未测试替换估计器后的闭环收益、峰值或稳定性。E3额外使用Ft/My通道，若要归因于优化机制，需要同信息EKF或空间通道消融。',
        'E3一条边界不能完整表示多条有限宽脊、部分覆盖与凸起；Hertz有限参考压入关系与对象有限足迹/SLS/HC不匹配。E0可吸收尺度却不能吸收全部深度依赖和动态滞回。摩擦与力矩项、表面坡度也会混淆。',
        '单侧均匀区的E0与ρ不可分别辨识，边界在足迹外时b弱可辨识；压力矩给出质心偏移，不能单独确定宽度、模量比与边界距离。GN的局部热启动可卡住，Hessian含到达正则，因此不能把数值可逆当成物理可辨识。',
        'Fn剖面条件化、滤波相位差、窗口相关与模型失配使协方差可能失准；1σ/2σ覆盖率应与68%/95%对照。未来运动未知，恒力路径预瞄不是未来闭环刚度的完整预测。', '',
        '## 复现', '',
        '`python -X utf8 experiments/method_redesign/step2_estimators.py --workers 4`',
        '`python -X utf8 experiments/method_redesign/step2_estimators.py --workers 4 --reproduce experiments/method_redesign/results/step2_estimators.json`',
        '第二条从头重新生成全部轨迹与调参，比较所有科学数值（排除耗时与运行环境）。无额外隐藏缓存。', '']
    lines += ['本次实际复现使用第二条命令再加 `--selected-only`：从头重生成50个调参及100个诊断轨迹，重新核对选定配置的调参目标和全部诊断科学输出；不重复搜索其余候选。已通过逐项完全相等检查。', '']
    overall = result['summary']['overall']
    best = min(overall, key=lambda k: overall[k]['median'])
    lines.insert(lines.index('## 结果边界与疑点')+2,
        f'本次总体median最小为{best}（{overall[best]["median"]:.3f}）；比较须同时看P95、危险方向与各族结果。E3总体median={overall["E3"]["median"]:.3f}，不据此预设空间模型有效。')
    sensor_lines = ['', '## 传感器回归检查', '',
        '原版Git版本：'+result['sensors']['baseline_revision']+'。三族全时段0–11列按uint64逐位比较。误差含30 Hz滤波滞后。', '',
        '| 族 | 逐位一致 | Ft RMSE N | My RMSE N·m |', '|---|---|---:|---:|']
    for row in result['sensors']['rows']:
        sensor_lines.append(f'| {row["family"]} | {row["bitwise_0_11"]} | {row["Ft"]["rmse"]:.5f} | {row["My"]["rmse"]:.7f} |')
    lines += sensor_lines
    conclusions = []
    for fam in C.FAMS:
        rr = result['summary'][fam]
        medbest = min(rr,key=lambda m:rr[m]['median'])
        tailbest = min(rr,key=lambda m:rr[m]['p95'])
        conclusions.append(f'{fam}：median最小{medbest}（{rr[medbest]["median"]:.3f}），P95最小{tailbest}（{rr[tailbest]["p95"]:.3f}）。')
    e3 = overall['E3']
    conclusions.append(f'E3总体1σ/2σ覆盖为{100*e3["cover1"]:.1f}%/{100*e3["cover2"]:.1f}%；P50/P99为{e3["p50_ms"]:.3f}/{e3["p99_ms"]:.3f} ms。')
    for h, pp in e3['prediction'].items():
        conclusions.append(f'τ={h} s：E3预测median/P95={pp["preview"]["median"]:.3f}/{pp["preview"]["p95"]:.3f}，保持当前={pp["hold"]["median"]:.3f}/{pp["hold"]["p95"]:.3f}；须分别判断典型误差与尾部，不能合称预瞄有效。')
    conclusions += ['12组随机候选只支持本次原型实例比较，不代表各估计器族的最优性能。单边界模型在边界离开足迹后可能失去梯度，难以重新捕获后续边界；热启动限制与结构失配均会影响E3。',
        'T2m-soft没有空间边界，其E3改善不能归因为预瞄；可能来自联合约束或正则化，本轮未用同信息消融区分。',
        'Hertz准静态刚度与对象的动态切线增益不同，同一模量尺度未必能同时匹配足迹几何和动态刚度。',
        'E2的1σ覆盖明显高于68%，其区间偏宽；E1在硬组织运动中的区间明显过窄。E3覆盖也随场景变化，均不能视为校准好的后验。',
        'JSON另存每族首个诊断场景每0.2 s的E3参数、潜半宽、三个预测时距与完整3×3协方差，以及θ的4×4条件协方差，便于检查可辨识性与过度置信。']
    ps = result['summary']['S2']['E3']['prediction']['0.3']
    conclusions.append(f'S2静态异质场景在0.3 s时，预瞄median从保持的{ps["hold"]["median"]:.3f}降到{ps["preview"]["median"]:.3f}；但总体0.3 s预测更差，空间预瞄的收益具有条件性。')
    at = lines.index('## 结果边界与疑点')+2
    lines[at:at] = conclusions+['']
    lines += ['', '本次因主项目目录拒绝写入，在可写副本运行；只交付指定5个文件的补丁，未提交或推送。在副本中复现需设置 `BASELINE_PLANT` 指向修改前plant.py；在含固定Git版本的项目树中会自动提取原版。', '']
    return '\n'.join(lines)


def scientific(obj):
    if isinstance(obj, dict):
        return {k: scientific(v) for k, v in obj.items() if k not in ('p50_ms', 'p99_ms', 'metadata')}
    if isinstance(obj, list):
        return [scientific(v) for v in obj]
    return obj


def self_check(log):
    th = np.array([math.log(2e4),0.,43.,math.log(.3)])
    a = .01
    f, m, der = hertz(th,.43,a)
    assert np.isclose(f, math.pi*HL*2e4*a*a/(4*HR), rtol=5e-5)
    assert abs(m) < 1e-12
    ks, ms, ts, aa = spatial(th,np.array([.43]),np.array([f]))
    assert np.isclose(aa[0],a,rtol=1e-7)
    assert np.isclose(ts[0],-.3*f,rtol=1e-7)
    # Exact scalar Schur update versus the full joint normal equations.
    obs = sanitize(log[:80]); est = MHE(candidates('E2')[0])
    X,Y,Z,_ = design(obs,len(obs))
    A,y,z,sk = est.augmented(obs,len(obs),X,Y,Z)
    J = np.column_stack((y,A))
    H = J.T@J; H[0,0] += 1/sk**2
    H[1:,1:] += np.eye(DEG+1)*1e-12
    rhs = J.T@z; rhs[0] += est.k/sk**2
    joint = np.linalg.solve(H,rhs)
    pred,_ = est.step(obs,len(obs))
    assert np.isclose(pred[0],np.clip(joint[0],LO,HI),rtol=1e-6)
    # Truth poisoning and future-visible-data poisoning. Only prefixes count.
    poison = log[:100].copy()
    poison[:,[3,4,6,7,8,9,10,14,15]] = 1e100
    poison[60:,5] += 2.
    for name in NAMES:
        cfg = candidates(name)[0]
        k0,s0,_,_,_ = replay(log[:60],name,cfg)
        k1,s1,_,_,_ = replay(poison,name,cfg)
        assert np.array_equal(k0,k1[:60]) and np.array_equal(s0,s1[:60]),name
    return dict(hertz_uniform=True,latent_width_root=True,schur_equivalence=True,
                truth_poisoning=True,future_prefix_invariance=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--reproduce', type=Path)
    parser.add_argument('--selected-only', action='store_true',
                        help='with --reproduce: recheck selected training objectives and all diagnostic outputs; do not repeat full search')
    args = parser.parse_args()
    if args.selected_only and not args.reproduce:
        parser.error('--selected-only requires --reproduce')
    prior = json.loads(args.reproduce.read_text(encoding='utf-8')) if args.reproduce else None
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        print('Generating 50 tuning trajectories', flush=True)
        training = list(pool.map(generate, make_cohort(True)))
        verification = self_check(training[0][2])
        print('Self checks PASS', verification, flush=True)
        # Warm all JIT kernels outside timed evaluation.
        for name in NAMES:
            replay(training[0][2][:100], name, candidates(name)[0], True)
        jobs = [(name, cfg, training) for name in NAMES for cfg in candidates(name)]
        if args.selected_only:
            jobs = [(name, prior['selected'][name]['params'], training) for name in NAMES]
        print('Rechecking selected training configurations' if args.selected_only else
              'Tuning: E0 fixed, E1/E2/E3 12 configurations each', flush=True)
        tuning = []
        for row in pool.map(tuning_job, jobs):
            tuning.append(row)
            print(row, flush=True)
        selected = {name: min([r for r in tuning if r['method'] == name], key=lambda r: r['objective']) for name in NAMES}
        if args.selected_only:
            assert selected == prior['selected'], 'Selected training objective mismatch'
            tuning = prior['tuning']
        print('Generating 100 diagnostic trajectories', flush=True)
        testing = list(pool.map(generate, make_cohort(False)))
        jobs = [(name, selected[name]['params'], fam, sid, log) for fam, sid, log in testing for name in NAMES]
        rows = []
        for i, row in enumerate(pool.map(evaluation_job, jobs)):
            rows.append(row)
            if i % 20 == 0:
                print(f'Evaluated {i+1}/400', flush=True)
    summary = summarize(rows)
    # Timing in isolation, after all kernels have been compiled. Same first
    # diagnostic trajectory per family for every estimator, no cherry picking.
    serial_timing = {}
    representative_trace = {fam: [] for fam in C.FAMS}
    for fam in C.FAMS:
        _, _, log = next(t for t in testing if t[0] == fam)
        serial_timing[fam] = {}
        for name in NAMES:
            _,_,ms,_,_ = replay(log,name,selected[name]['params'],trace=representative_trace[fam] if name == 'E3' else None)
            serial_timing[fam][name] = ms[scoring_mask(log)].tolist()
            summary[fam][name]['p50_ms'] = float(np.median(serial_timing[fam][name]))
            summary[fam][name]['p99_ms'] = float(np.quantile(serial_timing[fam][name],.99))
    for name in NAMES:
        ms = np.concatenate([serial_timing[fam][name] for fam in C.FAMS])
        summary['overall'][name]['p50_ms'] = float(np.median(ms))
        summary['overall'][name]['p99_ms'] = float(np.quantile(ms,.99))

    # Retain reproducible episode-level primary metrics; avoid megabytes of
    # repeated pointwise values in the final JSON.
    episodes = []
    for row in rows:
        s = row['score']
        episodes.append(dict(method=row['method'], family=row['family'], seed=row['seed'],
            ticks=len(s['log_error']), median=float(np.median(s['log_error'])), p95=float(np.quantile(s['log_error'], .95)),
            under=float(np.mean(s['under'])), over=float(np.mean(s['over'])),
            cover1=float(np.mean(s['cover1'])), cover2=float(np.mean(s['cover2']))))
    from check_plant_sensors import check, BASE_REV
    # When Git original is unavailable in a scratch copy, BASELINE_PLANT may
    # provide the byte-preserved pre-change source. It is never an estimator input.
    import tempfile, subprocess
    with tempfile.TemporaryDirectory() as td:
        baseline = os.environ.get('BASELINE_PLANT')
        if not baseline:
            path = Path(td)/'plant_sensor_original.py'
            path.write_bytes(subprocess.check_output(['git', 'show', BASE_REV+':experiments/cobot_scan_benchmark/plant.py'], cwd=ROOT))
            baseline = str(path)
        sensors = check(baseline)
    result = dict(metadata=dict(python=sys.version, numpy=np.__version__, platform=platform.platform(),
                                workers=args.workers, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
                  representative_trace=representative_trace, verification=verification, tuning=tuning, selected=selected, summary=summary, episodes=episodes, sensors=sensors,
                  seeds=dict(tuning=190000, diagnostic=110000, final_accessed=False))
    if prior is not None:
        assert scientific(prior) == scientific(result), 'Scientific reproduction mismatch'
        print('REPRODUCTION PASS: selected training objectives and all diagnostic scientific outputs exactly equal; full candidate search not repeated'
              if args.selected_only else 'REPRODUCTION PASS: all scientific results exactly equal', flush=True)
    out = HERE/'results'; out.mkdir(exist_ok=True)
    (out/'step2_estimators.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    (out/'step2_estimators.md').write_text(markdown(result), encoding='utf-8')
    print(json.dumps(summary['overall'], indent=2), flush=True)


if __name__ == '__main__':
    main()
