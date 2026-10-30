"""Causal constant-acceleration IMM for equivalent tissue displacement (SI).

Blom & Bar-Shalom (1988), doi:10.1109/9.1299: interaction/moment matching,
model-conditioned Kalman filtering, Bayesian model weights, output combination.
The model index changes the continuous white-jerk spectral density q [m^2/s^5].
This estimates equivalent contact drift, not necessarily physical surface motion.
"""
import math

import numpy as np


class IMMTracker:
    """Two CA models; ``single_model=True`` gives the matched KF ablation.

    ``reset(z0)`` consumes the first observation without advancing time.
    Subsequent ``update(z)`` calls advance exactly Tc. ``predict(L)`` is read-only.
    R is displacement measurement variance [m^2] and may be changed between ticks.
    Initial velocity/acceleration standard deviations are 0.03 m/s and 0.3 m/s^2.
    """

    def __init__(self, Tc=0.008, q_smooth=1e-5, q_maneuver=0.1,
                 p_stay=0.98, R=(10e-6)**2, single_model=False,
                 prob_floor=1e-12):
        vals = (Tc, q_smooth, q_maneuver, p_stay, R, prob_floor)
        if not all(np.isfinite(v) for v in vals):
            raise ValueError('IMM parameters must be finite')
        if Tc <= 0 or q_smooth < 0 or q_maneuver < 0:
            raise ValueError('Tc must be positive; q must be nonnegative')
        if not 0 <= p_stay <= 1 or not 0 < prob_floor < 0.5:
            raise ValueError('invalid p_stay or probability floor')
        self.Tc = float(Tc)
        self.prob_floor = float(prob_floor)
        self.R = R
        self.n_models = 1 if single_model else 2
        self.transition = (np.ones((1, 1)) if single_model else
                           np.array([[p_stay, 1-p_stay], [1-p_stay, p_stay]]))
        dt = self.Tc
        self.A = np.array([[1., dt, dt*dt/2], [0., 1., dt], [0., 0., 1.]])
        # Exact discretization of ds=v dt, dv=a dt, da=sqrt(q) dW.
        Q1 = np.array([[dt**5/20, dt**4/8, dt**3/6],
                       [dt**4/8, dt**3/3, dt**2/2],
                       [dt**3/6, dt**2/2, dt]])
        qs = [q_smooth] if single_model else [q_smooth, q_maneuver]
        self.Q = np.array([q*Q1 for q in qs])
        self.initialized = False
        self._mu = np.full(self.n_models, 1/self.n_models)
        self.x = np.zeros(3)
        self.P = np.zeros((3, 3))

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, value):
        if not np.isfinite(value) or value <= 0:
            raise ValueError('R must be a finite positive variance')
        self._R = float(value)

    @property
    def mu(self):
        """Posterior probabilities (copy), ordered smooth, maneuver."""
        return self._mu.copy()

    @staticmethod
    def _sym(P):
        return (P+P.T)*0.5

    def _combine(self):
        self.x = self._mu @ self.x_models
        delta = self.x_models-self.x
        self.P = self._sym(np.einsum('i,ijk->jk', self._mu, self.P_models)
                           + np.einsum('i,ij,ik->jk', self._mu, delta, delta))

    def reset(self, z0):
        if not np.isfinite(z0):
            raise ValueError('z0 must be finite')
        self._mu = np.full(self.n_models, 1/self.n_models)
        self.x_models = np.tile([float(z0), 0., 0.], (self.n_models, 1))
        P0 = np.diag([self.R, 0.03**2, 0.3**2])
        self.P_models = np.tile(P0, (self.n_models, 1, 1))
        self.initialized = True
        self._combine()
        return self.x.copy()

    def update(self, z):
        if not np.isfinite(z):
            raise ValueError('measurement must be finite')
        if not self.initialized:
            return self.reset(z)
        # 1. Interaction: c_j = sum_i p_ij mu_i, mu_i|j = p_ij mu_i/c_j.
        joint = self._mu[:, None]*self.transition
        prior = joint.sum(axis=0)
        mixing = joint/np.maximum(prior, np.finfo(float).tiny)
        mixed_x = mixing.T @ self.x_models
        mixed_P = np.empty_like(self.P_models)
        for j in range(self.n_models):
            d = self.x_models-mixed_x[j]
            mixed_P[j] = self._sym(np.einsum('i,ijk->jk', mixing[:, j], self.P_models)
                                  + np.einsum('i,ij,ik->jk', mixing[:, j], d, d))
        # 2. Model-conditioned KF, H=[1,0,0]. Joseph covariance update.
        log_likelihood = np.empty(self.n_models)
        for j in range(self.n_models):
            xp = self.A @ mixed_x[j]
            Pp = self._sym(self.A @ mixed_P[j] @ self.A.T+self.Q[j])
            innovation = float(z)-xp[0]
            S = max(float(Pp[0, 0]+self.R), np.finfo(float).tiny)
            K = Pp[:, 0]/S
            self.x_models[j] = xp+K*innovation
            J = np.eye(3)
            J[:, 0] -= K
            self.P_models[j] = self._sym(J @ Pp @ J.T+self.R*np.outer(K, K))
            log_likelihood[j] = -0.5*(math.log(2*math.pi)+math.log(S)+innovation**2/S)
        # 3. Bayes weights in log space; affine floor preserves sum=1 and mu>=floor.
        log_w = np.log(np.maximum(prior, np.finfo(float).tiny))+log_likelihood
        w = np.exp(log_w-log_w.max())
        w /= w.sum()
        self._mu = self.prob_floor+(1-self.n_models*self.prob_floor)*w
        # 4. Output moment matching includes the between-model covariance.
        self._combine()
        return self.x.copy()

    def predict(self, L):
        """CA extrapolation of each posterior, then mixture; returns (s, s_dot)."""
        if not self.initialized:
            raise RuntimeError('call reset or update before predict')
        if not np.isfinite(L) or L < 0:
            raise ValueError('prediction horizon must be finite and nonnegative')
        pos = self.x_models @ np.array([1., L, 0.5*L*L])
        vel = self.x_models @ np.array([0., 1., L])
        return float(self._mu @ pos), float(self._mu @ vel)
