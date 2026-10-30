"""Proposed controller (paper draft §三) and its reference/fixed variants.

MotionAware: every control period
 1. small known excitation, specified as a force amplitude and converted by 1/k_bar;
 2. one window fit z = k Y + B beta (cubic B, eliminated first) -> k_hat, sigma_k, drift g;
 3. conservative stiffness k_bar = max(k_hat + c sigma_k, (1-rho_d) k_bar_prev); windows that
    fail the validity test (chi2 too large, or relative precision sigma_k/k_hat above rel_max)
    may only raise k_bar;
 4. unit-stiffness admittance A1 computed in force units, command increments divided by k_bar;
 5. drift feed-forward: predicted drift increment across the command delay and servo lag,
    divided by k_bar; the constant part is left to the loop.
Feedback is divided by the conservative k_bar; the open-loop feed-forward and the excitation by the
unbiased k_hat. Contact loss: freeze k_bar, no feed-forward/excitation, approach at a limited speed,
resume normal control on re-contact (force >= 0.5 N) with the admittance state reset.

Oracle: same A1 and feed-forward structure, but k_bar = true local stiffness and the
feed-forward uses the true future surface height (privileged; reference only).
FixedEnd: same A1 with a constant k_bar (the method's own lower or upper bound), no
estimation, no feed-forward, no excitation.
"""
import json, math, os
import numpy as np
import common
from common import P, FD, K_MIN, K_MAX

A1 = json.load(open(os.path.join(common.HERE, 'A1_design.json')))


class _UnitAdmittance:
    """M1 w'' + D1 w' = e + ki I, exact ZOH per period; returns the force-domain increment."""
    def __init__(self, M1, D1, KPF1, ki, Tc):
        self.M1, self.D1, self.KPF1, self.ki, self.Tc = M1, D1, KPF1, ki, Tc
        self.eD = math.exp(-D1*Tc/M1)
        self.reset()

    def reset(self):
        self.wd = 0.; self.I = 0.; self.e_prev = None

    def increment(self, e):
        if self.e_prev is None:
            self.e_prev = e
        self.I = min(max(self.I+e*self.Tc, -20.), 20.)
        wd2 = self.wd*self.eD+(1-self.eD)*(e+self.ki*self.I)/self.D1
        dw = 0.5*(self.wd+wd2)*self.Tc
        self.wd = wd2
        du = dw+self.KPF1*(e-self.e_prev)
        self.e_prev = e
        return du


class MotionAware(common.DirectPolicy):
    name = 'Proposed'

    def __init__(self, M1=None, D1=None, KPF1=None, ki=None, a_force=.25, f_exc=20.8,
                 nw=16, degree=3, noise=.04, c=2., rho_d=.02, chi2_max=4., rel_max=.25,
                 feedforward=True, conservative=True, excitation=True, v_app=.006, adaptive_exc=False,
                 sector_lo=.5, tag=None):
        d = A1 if M1 is None else dict(M1=M1, D1=D1, KPF1=KPF1, ki=ki)
        self.params = dict(M1=d['M1'], D1=d['D1'], KPF1=d['KPF1'], ki=d['ki'] if ki is None else ki,
                           a_force=a_force, f_exc=f_exc, nw=nw, degree=degree, noise=noise, c=c,
                           rho_d=rho_d, chi2_max=chi2_max, rel_max=rel_max, feedforward=feedforward,
                           conservative=conservative, excitation=excitation, v_app=v_app,
                           adaptive_exc=adaptive_exc, sector_lo=sector_lo)
        if tag:
            self.name = tag

    # hooks overridden by the oracle
    def stiffness(self, o, n):
        p = self.params
        fit = common.window_fit(o.log, n, p['nw'], p['degree'], p['noise'], self.ns, self.eb,
                                self.afc, self.ndel)
        self.fit = fit
        if not np.isfinite(fit['k']) or fit['info'] <= 0:
            return self.kbar, False
        k_hat = min(max(fit['k'], K_MIN), K_MAX)
        sig = 1.0/math.sqrt(fit['info'])            # info = ||Y_perp||^2/noise^2
        self.k_hat = k_hat; self.rel = sig/k_hat
        if not p['conservative']:
            return k_hat, True
        cand = min(max(k_hat+p['c']*sig, K_MIN), K_MAX)
        precise = sig <= p['rel_max']*k_hat
        consistent = fit['chi2'] <= p['chi2_max']
        if precise:
            self.k_good = k_hat                      # last precise estimate (for FF/excitation)
        if self.first:
            if precise and consistent:
                self.first = False
                return cand, True
            return self.kbar, False
        if precise and consistent:                   # valid window
            return max(cand, (1-p['rho_d'])*self.kbar), True
        if precise:                                  # model mismatch (e.g. fast push/edge):
            return max(k_hat, self.kbar), False      # may raise to the estimate, never lower
        return self.kbar, False                      # too little information: hold

    def drift_pred(self, o, n):
        """predicted equivalent force drift at the time this command acts."""
        fit = self.fit
        if fit is None or not np.isfinite(fit['k']):
            return None
        h = (self.ndel+0.5)*self.Tc
        g, dg = common.drift_poly(fit, h)
        return g+self.tau*dg

    def reset(self, platform):
        super().reset(platform)
        p = self.params
        self.adm = _UnitAdmittance(p['M1'], p['D1'], p['KPF1'], p['ki'], self.Tc)
        self.kbar = K_MAX; self.k_hat = np.nan; self.k_good = np.nan; self.first = True; self.fit = None
        self.k_seek = 3e4                              # excitation scale before the first precise window
        self.exc_gain = 1.; self.rel = np.inf
        self.mode = 0; self.inband = 0
        self.q0 = None; self.q_prev = None; self.x_prev = 0.

    def step(self, o):
        p = self.params
        n, log = o.n, o.log
        F, t = log[n-1, 5], log[n-1, 0]
        sent = log[n-1, 11]
        e = FD-F
        # contact supervision
        if self.mode == 0 and F < .3:
            self.mode = 1; self.inband = 0
        if self.mode == 1:                                   # lost: approach
            if F >= .5:
                self.mode = 0; self.adm.reset(); self.q_prev = None; self.first_after = True
            else:
                self._trace(t, 1)
                return sent+p['v_app']*self.Tc
        k_bar, _ = self.stiffness(o, n)
        self.kbar = k_bar
        du = self.adm.increment(e)
        # drift feed-forward (normal mode only)
        dq = 0.
        if p['feedforward'] and self.mode == 0:
            q = self.drift_pred(o, n)
            if isinstance(q, tuple):                       # oracle: surface height [m]
                if self.q_prev is not None:
                    dq = self.kbar*(q[1]-self.q_prev)
                self.q_prev = q[1]
                q = None
            if q is not None:
                ye = log[n-1, 2]
                if self.q_prev is not None:
                    # each window's drift is referenced to its own end position y_e; remove the
                    # probe's own displacement so that only tissue-side drift is fed forward
                    k_ref = self.k_good if np.isfinite(self.k_good) else self.kbar
                    own = 0. if getattr(self, 'ff_absolute', False) or self.fit is None else k_ref*(ye-self.ye_prev)
                    dq = (q-self.q_prev)-own
                self.q_prev = q; self.ye_prev = ye
        else:
            self.q_prev = None
        # excitation (force units), only in normal contact
        if p['adaptive_exc']:
            # excite only as much as the sector requires: target relative precision
            # eps = (1-a)/(2c) (from c*sigma_k/k <= (1-a)/2); no new parameter
            eps = (1-p['sector_lo'])/(2*p['c'])
            want = min(max((self.rel/eps-0.5)/0.5, 0.), 1.) if np.isfinite(self.rel) else 1.
            self.exc_gain += 0.1*(want-self.exc_gain)
        x = p['a_force']*self.exc_gain*math.sin(2*math.pi*p['f_exc']*(t+self.Tc)) if (p['excitation'] and self.mode == 0) else 0.
        dx = x-self.x_prev; self.x_prev = x
        self._trace(t, self.mode)
        # feedback uses the conservative k_bar (stability); the open-loop feed-forward and the
        # excitation use the unbiased estimate (an error there only mis-scales compensation)
        if not np.isfinite(self.k_good):
            # information seeking: enlarge the excitation until a precise window is obtained
            self.k_seek = max(self.k_seek*0.85, K_MIN)
        k_un = self.k_good if np.isfinite(self.k_good) else self.k_seek
        return sent+du/self.kbar+(dx-dq)/k_un

    def _trace(self, t, mode):
        self.trace.append((t, self.kbar, self.k_hat, mode))


class Oracle(MotionAware):
    """privileged reference: true local stiffness for the loop and the true future surface
    height for the feed-forward (scaled by the true stiffness). Not deployable."""
    name = 'Oracle'

    def __init__(self, p_scene=None, **kw):
        kw.setdefault('excitation', False)
        super().__init__(**kw)
        self.p_scene = None if p_scene is None else np.asarray(p_scene, float)

    def stiffness(self, o, n):
        k = o.log[n-1, 8]
        self.fit = None
        k = (min(max(k, K_MIN), K_MAX) if k > 0 else self.kbar)
        self.k_good = k
        return k, True

    def drift_pred(self, o, n):
        x, t = o.log[n-1, 1], o.log[n-1, 0]
        h = (self.ndel+0.5)*self.Tc
        hh, ht, _ = P.surface(x, t+h, self.p_scene)
        return self.kbar*(hh+self.tau*ht), hh+self.tau*ht


class FixedEnd(MotionAware):
    """A1 at a constant stiffness bound; no estimation, feed-forward or excitation."""
    def __init__(self, k_fixed, **kw):
        kw.update(feedforward=False, excitation=False)
        super().__init__(**kw)
        self.k_fixed = k_fixed
        self.name = f'Fixed(k={k_fixed:g})'

    def stiffness(self, o, n):
        self.fit = None
        return self.k_fixed, True
