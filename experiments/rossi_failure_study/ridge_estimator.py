"""Projection relaxed by a prior on d: M = rho Dn'Dn (I + rho Dn'Dn)^-1 replaces I - B(B'B)^-1B'.

k = Y'Mz / Y'MY. As rho grows, Dn recovers a degree n-1 polynomial d of the same window, so n = 4 contains the cubic projection.
With memory, y first passes through a weak power-law response normalised at the probe frequency, so k stays the gain there.
"""
import math

import numpy as np

import common
from estimators import _clip_k, _fit_record

# Population weak power law, beta about .1; the branches differ from the plant's.
POP_TAU = np.array([.003, .03, .3])
POP_ALPHA = np.array([.41, .33, .26])
MIN_TICKS = 8     # half a projection window; a shorter run has too few probe cycles


class RidgeProjection:
    name = 'Ridge'

    def __init__(self, nw=32, order=4, fc=14., memory=False, sigma_scale=1., chi_scale=1., alpha=None, tau=None, gamma=0., lam=0., increment=True):
        self.nw = int(nw); self.order = int(order); self.fc = float(fc); self.memory = bool(memory)
        self.sigma_scale = float(sigma_scale); self.chi_scale = float(chi_scale)
        self.alpha = POP_ALPHA if alpha is None else np.asarray(alpha, float); self.tau = POP_TAU if tau is None else np.asarray(tau, float)
        self.gamma = float(gamma); self.lam = float(lam); self.k_ref = None
        self.increment = bool(increment)  # False reproduces historical levels.

    def reset(self, policy):
        # The window grows from MIN_TICKS after contact or start, as the projection's does: no estimate must wait for a full window.
        self.rho = (2*math.sin(math.pi*self.fc*policy.Tc))**(-2*self.order)
        self.matrices = {}
        self.fit = _fit_record(degree=0, tw=self.nw*policy.Tc)
        self.private = None; self.count = 0; self.dof = max(self.nw-self.order-2, 1); self.k_ref = None
        self.Tc = policy.Tc; self.w0 = 2*math.pi*float(policy.params.get('f_exc', 20.8))
        if self.memory:
            z = np.exp(-1j*self.w0*policy.Tc); self.e = np.exp(-policy.Tc/self.tau)
            self.b = self.tau/policy.Tc*(1-self.e)      # exact for motion linear within a tick
            self.mag = abs(1+np.sum(self.alpha*self.b*(1-z)/(1-self.e*z)))
            self.state = np.zeros(len(self.tau))
        self.y_prev = None; self.y_ref = None; self.y_out = None
        self.rate_out = None; self.base_prev = None; self.rate_dy_prev = 0.

    def _matrix(self, m):
        if m not in self.matrices:
            D = np.diff(np.eye(m), self.order, axis=0)
            M = np.eye(m)-np.linalg.inv(np.eye(m)+self.rho*D.T@D)
            self.matrices[m] = (M, max(m-self.order-1, 1))     # residual dof of k and the prior's null space
        return self.matrices[m]

    def _log(self, o, n):
        if not self.memory and not self.lam:
            return o.log
        if self.private is None:
            self.private = np.zeros((4096, o.log.shape[1]))
        for j in range(self.count, n):
            row = o.log[j].copy(); y = float(row[2]); dy = 0.
            if self.y_prev is None:
                self.y_prev = y; self.y_ref = y; self.y_out = y
            else:
                dy = y-self.y_prev; self.y_prev = y
                if self.memory:
                    previous = self.state
                    self.state = self.e*self.state+self.alpha*self.b*dy
                    if self.increment:
                        self.y_out += (dy+float(np.sum(self.state-previous)))/self.mag
                    else:
                        self.y_out = self.y_ref+((y-self.y_ref)+float(np.sum(self.state)))/self.mag
                else:
                    self.y_out = y
            out = self.y_out
            tau = 0.
            if self.lam and self.k_ref:
                # Hunt-Crossley rate term, F = k d (1 + lam d'): the response leads motion by tau = lam F/k, normalised at the probe frequency
                kt = self.k_ref
                for _ in range(2):      # k_ref is the magnitude at the probe frequency, tau needs the elastic tangent
                    tau = self.lam*max(float(row[5]), .3)/kt; kt = self.k_ref/math.sqrt(1+(self.w0*tau)**2)
                if not self.increment:
                    out = out+tau*dy/self.Tc
                    out = self.y_ref+(out-self.y_ref)/math.sqrt(1+(self.w0*tau)**2)
            if self.lam and self.increment:
                if self.rate_out is None:
                    self.rate_out = self.y_out
                else:
                    self.rate_out += (self.y_out-self.base_prev+
                        tau*(dy-self.rate_dy_prev)/self.Tc)/math.sqrt(1+(self.w0*tau)**2)
                self.base_prev = self.y_out; self.rate_dy_prev = dy
                out = self.rate_out
            row[2] = out; self.private[j] = row
        self.count = n
        return self.private[:n]

    def estimate(self, o, n, policy):
        p = policy.params; log = self._log(o, n); nw = self.nw
        self.fit = _fit_record(degree=0, tw=nw*policy.Tc)
        j0 = max(policy.ndel+1, n-1-nw); m = n-2-j0+1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        _, Y, z, keep = common._build(log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel, te, ye, nw*policy.ns*1e-3, p['enc_y'])
        if self.gamma:
            # static stiffening, k = k_end (F/F_end)^gamma, so the estimate is the gain at the window end
            F = np.maximum(log[j0:j0+m, 5], .3); Y = Y*(F/F[-1])**self.gamma
        broken = np.flatnonzero(~keep)
        if len(broken):
            Y, z = Y[broken[-1]+1:], z[broken[-1]+1:]    # only the unbroken contact run ending now
        if len(Y) < MIN_TICKS:
            return np.nan, np.inf, np.inf
        M, dof = self._matrix(len(Y))
        my = M@Y; info = float(Y@my)
        if info <= 0:
            return np.nan, np.inf, np.inf
        k = float(my@z)/info
        # M/noise^2 is the inverse covariance of d + e under the prior, so these are the GLS variance of k and the GLS residual test
        res = z-k*Y
        sigma = self.sigma_scale*p['noise']/math.sqrt(info)
        chi = self.chi_scale*float(res@(M@res))/p['noise']**2/dof
        self.dof = dof
        k = _clip_k(k)
        if np.isfinite(k) and sigma <= .25*k:
            self.k_ref = k
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); self.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi
