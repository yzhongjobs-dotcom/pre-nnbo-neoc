"""Curvature candidates layered on the preserved effective lean estimator."""
import math
import time
import numpy as np
from joint_policy import JointPolicy
from contact_jet import fit


class CurvedContactPolicy(JointPolicy):
    def __init__(self, model='centered_curvature', require_resolved=False, **kw):
        super().__init__(rate=False, **kw)
        self.curve_model = model
        self.require_resolved = bool(require_resolved)

    def reset(self, platform):
        super().reset(platform)
        self.curvature_fits = []
        self.curvature_used = 0

    def _stiffness_est(self, o, n):
        feedback = super()._stiffness_est(o, n)
        if o.log[n-1, 0] < self.warm+.064+self.information_window*.008:
            return feedback
        tic = time.perf_counter()
        measured = o.log[:n, [0, 1, 2, 5, 11, 12, 13]]
        r = fit(measured, nw=self.information_window, tau=self.tau, fc=self.fc,
                noise=self.params['noise'], model=self.curve_model)
        if r is None:
            return feedback
        t_statistic = float(r['coefficients'][2]/np.sqrt(r['cov'][2, 2]))
        # Positive curvature branch motivated by convex elastic indentation.
        # This conditional test does not certify elasticity or a safe model.
        resolved = t_statistic > 1.96
        accepted = bool(200 <= r['k'] <= 3e5 and r['min_gain'] > 0 and
                        r['risk'] <= math.log(1.2) and r['chi2'] <= self.params['chi2_max'] and
                        (resolved or not self.require_resolved))
        self.curvature_fits.append(dict(time=float(o.log[n-1, 0]), k=r['k'], sigma=r['sigma'],
            curvature=float(r['coefficients'][2]/.001), curvature_t=t_statistic,
            accepted=accepted, resolved_positive=resolved, risk=r['risk'], chi2=r['chi2'],
            elapsed=time.perf_counter()-tic, lean_value_before=float(self.k_good)))
        if accepted:
            self.k_good = r['k']
            self.kg_rel = r['relative_sd']
            self.curvature_used += 1
        return feedback
