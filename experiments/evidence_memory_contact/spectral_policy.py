"""Controller using the adaptive spectral window. Same interface as the lean candidate.

The original 16-tick projection keeps its own stream for fallback. A window result replaces the
measurement gain only when its posterior uncertainty, which now includes the modelled interference,
fits the feed-forward tolerance. Feedback can optionally use the same estimate under the original
conservative rule. Excitation, IMM, admittance and recovery are unchanged unless stated.
"""
import math
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'contact_information'))
sys.path.insert(0, str(HERE.parent/'method_redesign'))
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
from joint_policy import JointPolicy
from spectral_window import AdaptiveSpectralWindow


class SpectralContact(JointPolicy):
    def __init__(self, change_prior=.5, window=32, corners=(.5, 4., 8., 12., 16., 24.), fit_noise=.04, gate=math.log(1.2),
                 chi_max=None, feedback=False, amplitude_scale=1., warm=1., **kw):
        super().__init__(wave=None, new_fit=False, rate=False, amplitude_scale=amplitude_scale, warm=warm, **kw)
        self.window_fit = AdaptiveSpectralWindow(n=window, corners=corners, noise=fit_noise, change_prior=change_prior)
        self.window = window; self.gate = gate; self.chi_max = chi_max; self.use_feedback = feedback

    def reset(self, platform):
        super().reset(platform)
        self.spectral_trace = []; self.accepted = 0; self.trials = 0

    def _stiffness_est(self, o, n):
        self.k_good = self.baseline_gain
        feedback = super(JointPolicy, self)._stiffness_est(o, n)
        self.baseline_gain = self.k_good
        t = float(o.log[n-1, 0]); m = self.window
        if t < self.warm+.064+m*self.Tc or n < m+2 or np.any(o.log[n-m-2:n, 5] <= .5):
            return feedback
        r = self.window_fit.fit(o.log[n-m-1:n, 2], o.log[n-m-1:n, 5])
        if r is None:
            return feedback
        k, sig = r['k'], r['sigma']
        self.trials += 1
        rel = sig/max(k, 1.)
        ok = 200 <= k <= 3e5 and rel <= self.gate and (self.chi_max is None or r['chi2'] <= self.chi_max)
        self.spectral_trace.append((t, k, sig, r['corner'], r['chi2'], float(ok)))
        if ok:
            self.k_good = k; self.kg_rel = rel; self.accepted += 1
        if self.use_feedback and 200 <= k <= 3e5:
            p = self.params
            cand = float(np.clip(k+p['c']*sig, 200., 3e5))
            if rel <= p['rel_max']:
                return max(cand, (1-p['rho_d'])*self.kbar), True
            return self.kbar, False
        return feedback
