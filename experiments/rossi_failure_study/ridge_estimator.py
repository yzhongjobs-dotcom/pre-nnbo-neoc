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

    def __init__(self, nw=32, order=4, fc=14., memory=False, sigma_scale=1., chi_scale=1.):
        self.nw = int(nw); self.order = int(order); self.fc = float(fc); self.memory = bool(memory)
        self.sigma_scale = float(sigma_scale); self.chi_scale = float(chi_scale)

    def reset(self, policy):
        # The window grows from MIN_TICKS after contact or start, as the projection's does: no estimate must wait for a full window.
        self.rho = (2*math.sin(math.pi*self.fc*policy.Tc))**(-2*self.order)
        self.matrices = {}
        self.fit = _fit_record(degree=0, tw=self.nw*policy.Tc)
        self.private = None; self.count = 0; self.dof = max(self.nw-self.order-2, 1)
        if self.memory:
            f0 = float(policy.params.get('f_exc', 20.8))
            z = np.exp(-1j*2*math.pi*f0*policy.Tc); self.e = np.exp(-policy.Tc/POP_TAU)
            self.mag = abs(1+np.sum(POP_ALPHA*(1-z)/(1-self.e*z)))
            self.state = np.zeros(3); self.y_prev = None; self.y_ref = None; self.y_out = None

    def _matrix(self, m):
        if m not in self.matrices:
            D = np.diff(np.eye(m), self.order, axis=0)
            M = np.eye(m)-np.linalg.inv(np.eye(m)+self.rho*D.T@D)
            self.matrices[m] = (M, max(m-self.order-1, 1))     # residual dof of k and the prior's null space
        return self.matrices[m]

    def _log(self, o, n):
        if not self.memory:
            return o.log
        if self.private is None:
            self.private = np.zeros((4096, o.log.shape[1]))
        for j in range(self.count, n):
            row = o.log[j].copy(); y = float(row[2])
            if self.y_prev is None:
                self.y_prev = y; self.y_ref = y; self.y_out = y
            else:
                self.state = self.e*self.state+POP_ALPHA*(y-self.y_prev)
                self.y_prev = y
                self.y_out = self.y_ref+((y-self.y_ref)+float(np.sum(self.state)))/self.mag
            row[2] = self.y_out; self.private[j] = row
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
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); self.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi
