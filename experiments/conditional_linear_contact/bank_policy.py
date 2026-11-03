"""V2: single-prior conditional linear contact filter with a memory bank.

Diagnosis of V0: a narrow motion prior lets slow tracking motion read as a soft contact and
the gain collapses; per-tick likelihoods cannot identify the stiffness mode because the
probe response is far below the noise in one interval. So:
  - the contact filter always uses a broad motion prior (it must not be narrower than any
    motion the tissue can do); the IMM motion tracker stays separate and unchanged;
  - memory is chosen by a statistic with enough evidence: a slow and a fast filter run in
    parallel, and the slow one is reseeded when they disagree beyond their uncertainty.
Initialization comes from the first precise projection window, never from truth.
"""
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent/'method_redesign'))
import common
from common import K_MIN, K_MAX
from v1_servo_id import ServoID
from conditional_filter import interval_measurement

TC = .008
A = np.array([[1., TC, TC**2/2], [0., 1., TC], [0., 0., 1.]])
QCA = np.array([[TC**5/20, TC**4/8, TC**3/6], [TC**4/8, TC**3/3, TC**2/2], [TC**3/6, TC**2/2, TC]])


class ContactKF:
    """[k, D, r, a]: linear Kalman filter, measurement first then time update."""
    def __init__(self, q_motion, q_lnk, R):
        self.q, self.ql, self.R = q_motion, q_lnk, R
        self.x = None; self.P = None; self.nis = np.nan

    def start(self, k, D, rel=.5):
        self.x = np.array([k, D, 0., 0.])
        self.P = np.diag([(rel*k)**2, 1., 1e2, 1e4])

    def step(self, z, H, dy):
        P, x = self.P, self.x
        S = float(H@P@H+self.R)
        K = P@H/S; v = z-float(H@x)
        self.nis = v*v/S
        x = x+K*v
        I = np.eye(4)-np.outer(K, H)
        P = I@P@I.T+self.R*np.outer(K, K)
        F = np.eye(4); F[1, 0] = dy; F[1:, 1:] = A
        kref = float(np.clip(x[0], K_MIN, K_MAX))
        Q = np.zeros((4, 4)); Q[0, 0] = (kref*self.ql)**2*TC; Q[1:, 1:] = kref*kref*self.q*QCA
        self.x = F@x; self.P = F@P@F.T+Q; self.P = (self.P+self.P.T)*.5

    @property
    def k(self):
        return float(self.x[0])

    @property
    def sig(self):
        return math.sqrt(max(self.P[0, 0], 0.))


class BankContact(ServoID):
    name = 'conditional-bank-v2'

    def __init__(self, q_motion=None, q_slow=.1, q_fast=1., reseed=1.5, reseed_hold=8, gate_rel=math.log(1.2),
                 noise_R=.04, switch_time=1.0, feedback='filter', gain='filter', probe_rescale=True, **kw):
        super().__init__(**kw)
        self.bank_args = dict(q_motion=q_motion if q_motion is not None else kw.get('q_maneuver', 1.86),
                              q_slow=q_slow, q_fast=q_fast, reseed=reseed, reseed_hold=reseed_hold,
                              gate_rel=gate_rel, noise_R=noise_R, switch_time=switch_time,
                              feedback=feedback, gain=gain, probe_rescale=probe_rescale)

    def reset(self, platform):
        super().reset(platform)
        b = self.bank_args
        self.slow = ContactKF(b['q_motion'], b['q_slow'], b['noise_R']**2)
        self.fast = ContactKF(b['q_motion'], b['q_fast'], b['noise_R']**2)
        self.next_j = 1; self.ready = False; self.exceed = 0; self.lost = False
        self.reseeds = 0; self.published = 0; self.total = 0; self.bank_trace = []

    def _advance(self, o, n):
        b = self.bank_args
        F = float(o.log[n-1, 5])
        if F <= .5:
            self.lost = True; self.next_j = max(self.next_j, n-1)
            return
        if self.lost:
            self.lost = False; self.next_j = max(self.next_j, n-1)
            if self.ready:
                for f in (self.slow, self.fast):
                    k, sig = f.k, f.sig
                    f.start(k, float(o.log[n-1, 5]), rel=max(sig/max(k, 1.), .2))
        for j in range(self.next_j, n-1):
            item = interval_measurement(o.log, j, self.ns, self.eb, self.afc, self.tau)
            if item is None:
                continue
            H, z = item
            if not self.ready:
                continue
            dy = float(o.log[j+1, 2]-o.log[j, 2])
            self.slow.step(z, H, dy); self.fast.step(z, H, dy)
            for f in (self.slow, self.fast):
                f.x[0] = float(np.clip(f.x[0], K_MIN, K_MAX))
            # Disagreement of the slow memory with the recent evidence: a sustained
            # exceedance means real change, not one noisy interval.
            gap = abs(math.log(self.fast.k/self.slow.k))
            unc = math.hypot(self.fast.sig/self.fast.k, self.slow.sig/self.slow.k)
            self.exceed = self.exceed+1 if gap > b['reseed']*unc else 0
            if self.exceed >= b['reseed_hold']:
                self.slow.x, self.slow.P = self.fast.x.copy(), self.fast.P.copy()
                self.exceed = 0; self.reseeds += 1
        self.next_j = max(self.next_j, n-1)

    def step(self, o):
        previous = self.x_prev
        command = super().step(o)
        # The probe is scaled by 1/k_good. A stale high gain shrinks the actual probe and with it the
        # information needed to correct that gain. Scale by the fast filter when it shows the gain fell.
        if self.bank_args['probe_rescale'] and self.ready and self.mode == 0:
            dx = self.x_prev-previous
            k_good = self.k_good if np.isfinite(self.k_good) else self.k_seek
            gap = math.log(max(self.slow.k, 1.)/max(self.fast.k, 1.))
            unc = math.hypot(self.fast.sig/self.fast.k, self.slow.sig/self.slow.k)
            if gap > unc and self.fast.k < k_good:
                command += dx*(1./self.fast.k-1./k_good)
        return command

    def _stiffness_est(self, o, n):
        b = self.bank_args; p = self.params
        t = float(o.log[n-1, 0])
        self._advance(o, n)
        base = super()._stiffness_est(o, n)     # projection path: k_good, k_bar under the original rule
        if not self.ready:
            if t >= b['switch_time'] and np.isfinite(self.k_good):
                # First accepted projection value seeds the bank; no truth is read.
                D = float(o.log[n-1, 5])
                self.slow.start(self.k_good, D); self.fast.start(self.k_good, D)
                self.ready = True; self.next_j = n-1
            return base
        k, sig = self.slow.k, self.slow.sig
        rel = sig/max(k, 1.)
        self.total += 1
        if b['gain'] == 'filter' and rel <= b['gate_rel']:
            self.k_good = float(np.clip(k, K_MIN, K_MAX)); self.published += 1
            self.kg_rel = rel
        if b['feedback'] == 'filter':
            # Conservative feedback rule from the original method, on the filter output.
            cand = float(np.clip(k+p['c']*sig, K_MIN, K_MAX))
            if rel <= p['rel_max']:
                return max(cand, (1-p['rho_d'])*self.kbar), True
            return self.kbar, False
        return base
