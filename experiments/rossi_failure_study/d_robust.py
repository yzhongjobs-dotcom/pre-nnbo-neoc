"""D⁴ with a heavy-tailed prior: one convex objective, no detection or switching.

    min_{k,d} ||z - kY - d||^2 + 2 rho s^2 sum_i h_c( (D⁴d)_i / s ),   s = noise / sqrt(rho),   h_c the standard Huber function (u^2/2 for |u| <= c)

For small differences the penalty is rho ||D⁴d||^2, the D⁴ prior exactly.

Quadratic priors weight d by frequency, so a short burst and the 20.8 Hz probe response, which share a band, cannot be told apart.
The fourth difference of a burst is non-zero at a few samples, that of the probe at every sample: a prior that is quadratic for small
innovations and linear for large ones (Huber) separates them by structure. For small innovations it is the D⁴ prior itself, so smooth
windows give the D⁴ estimate. c is Huber's constant (1.345 by default). Solved by iteratively reweighted least squares, each step one D⁴
solve with weights on the differences; the stiffness statistics use the final weighted metric.
"""
import math

import numpy as np
from scipy.linalg import solveh_banded

import common
import contact_candidates as C
from contact_candidates import RIDGE, CandidatePolicy, Preprocess
from estimators import _clip_k, _fit_record
from ridge_estimator import MIN_TICKS, RidgeProjection


class RobustRidge(RidgeProjection):
    def __init__(self, c=1.345, iters=6, **kw):
        super().__init__(**kw)
        self.c = float(c); self.iters = int(iters)

    def reset(self, policy):
        super().reset(policy)
        self._band_cache = {}

    def _difference_bands(self, m):
        """Map weights to lower bands: ab[q, j] = (rho D' W D)[j+q, j]."""
        if m not in self._band_cache:
            D = np.diff(np.eye(m), self.order, axis=0)
            basis = np.zeros((self.order+1, m, len(D)))
            for q in range(self.order+1):
                basis[q, :m-q] = self.rho*(D[:, q:]*D[:, :m-q]).T
            self._band_cache[m] = D, basis.reshape((-1, len(D)))
        return self._band_cache[m]

    def estimate(self, o, n, policy):
        """Same cold-start IRLS steps as estimate_reference; two RHS, no dense inverse.

        Mathematical audit (standard Huber h(u) = u^2/2 for |u| <= c):
        the existing steps target ||z-kY-d||^2 + 2*rho*s^2*sum h(Dd/s),
        rather than the objective written in the module header. Six steps are
        an approximation to its minimizer. Sigma/chi freeze the final weights;
        they are surrogate GLS quantities, not exact robust sampling statistics.
        This optimization deliberately preserves those existing choices.
        """
        p = policy.params; log = self._log(o, n); nw = self.nw
        self.fit = _fit_record(degree=0, tw=nw*policy.Tc)
        j0 = max(policy.ndel+1, n-1-nw); m = n-2-j0+1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        _, Y, z, keep = common._build(log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel, te, ye, nw*policy.ns*1e-3, p['enc_y'])
        broken = np.flatnonzero(~keep)
        if len(broken):
            Y, z = Y[broken[-1]+1:], z[broken[-1]+1:]
        m = len(Y)
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        D, basis = self._difference_bands(m)
        rhs = np.column_stack((Y, z))
        scale = p['noise']/math.sqrt(self.rho)
        w = np.ones(m-self.order)
        for it in range(self.iters):
            ab = (basis@w).reshape((self.order+1, m))
            ab[0] += 1.
            smooth = solveh_banded(ab, rhs, lower=True, overwrite_ab=True, check_finite=False)
            sy, sz = smooth[:, 0], smooth[:, 1]
            my = Y-sy; info = float(Y@my)
            if info <= 0:
                return np.nan, np.inf, np.inf
            k = float(my@z)/info
            if it == self.iters-1:
                break
            u = np.abs(D@(sz-k*sy))/scale
            w = np.minimum(1., self.c/np.maximum(u, 1e-12))
        dof = max(m-self.order-1, 1)
        res = z-k*Y
        sigma = self.sigma_scale*p['noise']/math.sqrt(info)
        chi = self.chi_scale*float(res@((z-sz)-k*my))/p['noise']**2/dof
        self.dof = dof; self.weights = w
        k = _clip_k(k)
        if np.isfinite(k) and sigma <= .25*k:
            self.k_ref = k
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); self.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi

    def estimate_reference(self, o, n, policy):
        p = policy.params; log = self._log(o, n); nw = self.nw
        self.fit = _fit_record(degree=0, tw=nw*policy.Tc)
        j0 = max(policy.ndel+1, n-1-nw); m = n-2-j0+1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        _, Y, z, keep = common._build(log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel, te, ye, nw*policy.ns*1e-3, p['enc_y'])
        broken = np.flatnonzero(~keep)
        if len(broken):
            Y, z = Y[broken[-1]+1:], z[broken[-1]+1:]
        m = len(Y)
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        D = np.diff(np.eye(m), self.order, axis=0); DtD = D.T@D
        scale = p['noise']/math.sqrt(self.rho)          # prior scale of one fourth difference
        w = np.ones(m-self.order)
        for it in range(self.iters):
            S = np.linalg.inv(np.eye(m)+self.rho*(D.T*w)@D)
            M = np.eye(m)-S
            my = M@Y; info = float(Y@my)
            if info <= 0:
                return np.nan, np.inf, np.inf
            k = float(my@z)/info
            if it == self.iters-1:
                break
            u = np.abs(D@(S@(z-k*Y)))/scale
            w = np.minimum(1., self.c/np.maximum(u, 1e-12))
        dof = max(m-self.order-1, 1)
        res = z-k*Y
        sigma = self.sigma_scale*p['noise']/math.sqrt(info)
        chi = self.chi_scale*float(res@(M@res))/p['noise']**2/dof
        self.dof = dof; self.weights = w
        k = _clip_k(k)
        if np.isfinite(k) and sigma <= .25*k:
            self.k_ref = k
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); self.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi


class ContactRobust(Preprocess, RobustRidge):
    def __init__(self, kind, c=1.345, iters=6, **ridge):
        RobustRidge.__init__(self, c=c, iters=iters, memory=False, **dict(RIDGE, **ridge))
        self.kind = kind

    def reset(self, policy):
        RobustRidge.reset(self, policy)
        self._setup(self.kind, policy)


def build(sc, c=1.345, iters=6, overrides=None, ridge=None):
    """The selected Method (H on the D⁴ host) with the heavy-tailed prior, through method3's interfaces."""
    import method3 as M
    import sweeps as S
    p = dict(S.FINAL_C); p.update(C.CTL); p.update(dict(M.TRACKER, **(overrides or {})))
    return M.attach(CandidatePolicy('H', ContactRobust('H', c, iters, **(ridge or {})), single_model=False, **p), sc)
