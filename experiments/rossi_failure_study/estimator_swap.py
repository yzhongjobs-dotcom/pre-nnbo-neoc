"""One controller with a replaceable causal stiffness estimator.

Only ``_stiffness_est`` is specialized.  The inherited MotionAwareIMM/V1
motion tracker, admittance, feed-forward, excitation and recovery code remains
the execution path for every estimator.
"""
import math
import numpy as np

from v1_servo_id import ServoID
from estimators import ProjectionEstimator, BandPassEstimator, R2Estimator


class SwapPolicy(ServoID):
    """ServoID whose sole policy-level replacement is ``estimator``."""
    name = 'SwapPolicy'

    def __init__(self, estimator, **kwargs):
        kwargs.setdefault('kgood_rule', 'consistent')
        kwargs.setdefault('conservative', True)
        kwargs.setdefault('meas_bound', 1.0)
        self.estimator = estimator
        super().__init__(**kwargs)

    def reset(self, platform):
        super().reset(platform)
        self.accept_trace = []
        self.estimator.reset(self)

    def _stiffness_est(self, o, n):
        """MotionAware._stiffness_est with only ``fit`` supplied by estimator."""
        p = self.params
        k_est, sig_est, chi_est = self.estimator.estimate(o, n, self)
        fit = dict(self.estimator.fit)
        fit['k'] = k_est; fit['chi2'] = chi_est
        self.fit = fit
        if not np.isfinite(k_est) or fit.get('info', 0.) <= 0:
            return self.kbar, False
        k_hat = min(max(k_est, 200.), 3e5)
        sig = float(sig_est)
        self.k_hat = k_hat; self.rel = sig/k_hat
        precise = sig <= p['rel_max']*k_hat
        if precise:
            self.k_check = k_hat
        consistent = chi_est <= p['chi2_max']
        self.accept_trace.append((float(o.log[n-1, 0]), float(k_hat), sig,
                                  float(chi_est), bool(precise and consistent)))
        if p['kgood_rule'] == 'filter':
            if precise and consistent:
                self._filter_kgood(k_hat, sig/k_hat)
        elif precise and (consistent or p['kgood_rule'] == 'precise'):
            ok_meas = p['rel_meas'] is None or sig <= p['rel_meas']*k_hat
            if p['meas_bound'] is not None:
                dof = self.estimator.dof
                bound = (sig/k_hat)*(1.+math.sqrt(dof*max(0., chi_est-p['meas_bound'])))
                within = bound <= math.log(1.2)
                better = p['meas_better'] and np.isfinite(self.k_good) and bound < abs(math.log(k_hat/self.k_good))
                ok_meas = ok_meas and (within or better)
            if ok_meas:
                self.k_good = k_hat
                self.f_set = max(float(o.log[n-1, 5]), .3)
                self.kg_rel = sig/k_hat
                if p.get('oracle_meas') == 'stale':
                    self.k_good = max(float(o.log[n-1, 8]), 1.)
        elif precise and p['kgood_rule'] == 'raise' and np.isfinite(self.k_good) and k_hat > self.k_good:
            self.k_good = k_hat
        if not p['conservative']:
            return k_hat, True
        cand = min(max(k_hat+p['c']*sig, 200.), 3e5)
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
            return max(k_hat, self.kbar), False
        return self.kbar, False


def make_estimator(name, **params):
    if name in ('projection', 'Projection'):
        return ProjectionEstimator(**params)
    if name in ('bandpass', 'BandPass', 'butterworth'):
        return BandPassEstimator(**params)
    if name in ('r2', 'R2'):
        return R2Estimator(**params)
    raise ValueError(name)


def projection_policy(**kwargs):
    return SwapPolicy(ProjectionEstimator(nw=16, degree=3), **kwargs)
