"""E4: flat cubic nuisance + force-domain squared-exponential drift prior.

No truth columns, future observations, plant constants, or recursive stiffness
prior. Covariance is conditional on a causal force scale, not calibrated truth.
"""
import math
import time
import numpy as np
import step2_estimators as S
from variants import MotionAwareIMM


class GLS:
    def __init__(self, cfg):
        self.ell, self.alpha = map(float, cfg)
        self.k, self.var = S.HI, S.HI**2
        self.cache = {}

    def geometry(self, X, relative_ticks):
        key = tuple(relative_ticks)
        if key not in self.cache:
            # Complete QR gives Q^T B = 0 exactly (up to floating point).
            basis, R = np.linalg.qr(X, mode='complete')
            rank = np.linalg.matrix_rank(X)
            Q = basis[:, rank:]
            t = np.asarray(relative_ticks)*S.DT
            kernel = np.exp(-.5*((t[:, None]-t[None, :])/self.ell)**2)
            eig, U = np.linalg.eigh(Q.T@kernel@Q)
            self.cache[key] = (Q@U, np.maximum(eig, 0.), np.linalg.pinv(X))
        return self.cache[key]

    def fit(self, log, n):
        out = dict(k=np.nan, beta=np.zeros(S.DEG+1), info=0., chi2=np.inf,
                   valid=0, tw=S.NW*S.DT, te=log[n-1, 0] if n else 0.)
        data = S.design(log, n)
        if data is None:
            return out
        X, Y, Z, idx = data
        # Same causal normal-force observations as E0. Units are N in the
        # equation-error domain z, not raw unfiltered displacement units.
        scale = max(float(np.median(np.abs(log[idx+1, 5]))), .5)
        V, eig, inverse_B = self.geometry(X, idx-(n-1))
        yr, zr = V.T@Y, V.T@Z
        omega = S.NOISE**2+(self.alpha*scale)**2*eig
        info = float(np.sum(yr*yr/omega))
        if info <= 0 or not np.isfinite(info):
            return out
        k = float(np.sum(yr*zr/omega)/info)
        residual = zr-k*yr
        out.update(k=k, beta=inverse_B@(Z-k*Y), info=info,
                   chi2=float(np.sum(residual**2/omega))/max(len(yr)-1, 1),
                   valid=len(Y), force_scale=scale,
                   sigma=1/math.sqrt(info), gate=bool(S.LO <= k <= S.HI and
                       2/math.sqrt(info)/k <= math.log(1.2)))
        return out

    def step(self, log, n):
        fit = self.fit(log, n)
        self.last = fit
        if np.isfinite(fit['k']) and fit['info'] > 0:
            self.k = float(np.clip(fit['k'], S.LO, S.HI))
            self.var = 1/fit['info']
        return np.repeat(self.k, 3), np.full((3, 3), self.var)


class GLSController(MotionAwareIMM):
    """FINAL controller; override only stiffness fit and measurement bound.

    The original feedback precision/consistency/rate/raise logic is retained;
    chi2 and sigma now use the GLS covariance. Measurement acceptance replaces
    b by 2 sigma/k, retaining FINAL's precision and consistency prerequisites.
    """
    def __init__(self, cfg):
        super().__init__(**dict(S.C.FINAL))
        self.gls_cfg = cfg

    def reset(self, platform):
        super().reset(platform)
        assert self.ns == 8 and self.ndel == 1
        self.gls = GLS(self.gls_cfg)
        self.accepted = 0
        self.attempted = 0
        self.fit_ms = []

    def _stiffness_est(self, o, n):
        begin = time.perf_counter_ns()
        fit = self.gls.fit(o.log, n)
        self.fit_ms.append((time.perf_counter_ns()-begin)*1e-6)
        self.fit = fit
        if not np.isfinite(fit['k']) or fit['info'] <= 0:
            return self.kbar, False
        self.attempted += 1
        p = self.params
        k = float(np.clip(fit['k'], S.LO, S.HI))
        sig = 1/math.sqrt(fit['info'])
        self.k_hat, self.rel = k, sig/k
        precise = sig <= p['rel_max']*k
        consistent = fit['chi2'] <= p['chi2_max']
        if precise:
            self.k_check = k
        if precise and consistent and fit['gate']:
            self.k_good = k
            self.kg_rel = sig/k
            self.f_set = max(float(o.log[n-1, 5]), .3)
            self.accepted += 1
        if not p['conservative']:
            return k, True
        cand = float(np.clip(k+p['c']*sig, S.LO, S.HI))
        if self.first:
            if precise and consistent:
                self.first = False
                return cand, True
            return self.kbar, False
        if precise and consistent:
            return max(cand, (1-p['rho_d'])*self.kbar), True
        if not consistent and p['edge_rule'] == 'raise':
            return max(cand, self.kbar), False
        if precise:
            return max(k, self.kbar), False
        return self.kbar, False
