"""Limited integration test: SLS in-phase probe gain in the existing controller.

This tests the gain contribution only. It does NOT implement a frequency-aware
compliance inverse or pass a weakly identified relaxed modulus into control.
"""
import math
import numpy as np
from joint_policy import JointPolicy
from memory_contact import fit


class MemoryContactPolicy(JointPolicy):
    def __init__(self, memory_model='memoryc_040', **kw):
        super().__init__(rate=False, **kw)
        self.memory_model = memory_model

    def reset(self, platform):
        super().reset(platform)
        self.memory_fits = []
        self.memory_used = 0

    def _stiffness_est(self, o, n):
        feedback = super()._stiffness_est(o, n)
        if o.log[n-1, 0] < self.warm+.064+self.information_window*.008:
            return feedback
        m = o.log[:n, [0, 1, 2, 5, 11, 12, 13]]
        r = fit(m, nw=self.information_window, tau=self.tau, noise=self.params['noise'], model=self.memory_model)
        if r is None:
            return feedback
        accepted = bool(200 <= r['k'] <= 3e5 and r['physical'] and r['min_gain'] > 0 and
                        r['risk'] <= math.log(1.2) and r['chi2'] <= self.params['chi2_max'])
        self.memory_fits.append(dict(time=float(o.log[n-1, 0]), k=r['k'], sigma=r['sigma'], risk=r['risk'],
            accepted=accepted, coefficients=r['coefficients'].tolist(), chi2=r['chi2'],
            low_frequency_sd=r['low_frequency_sd'], low_frequency_gain=r['low_frequency_gain']))
        if accepted:
            self.k_good = r['k']
            self.kg_rel = r['sigma']/r['k']
            self.memory_used += 1
        return feedback
