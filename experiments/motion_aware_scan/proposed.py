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
                 sector_lo=.5, kgood_rule='consistent', edge_rule='legacy', recovery=True, follow_tau=0., ff_band=None, kg_avg=2., rel_meas=None, kg_gamma=0., gain_noise=False, meas_bound=None, meas_better=False, recontact=False, enc_y=False, ff_jump=None, oracle_meas=False, oracle_fb=False, tag=None):
        d = A1 if M1 is None else dict(M1=M1, D1=D1, KPF1=KPF1, ki=ki)
        self.params = dict(M1=d['M1'], D1=d['D1'], KPF1=d['KPF1'], ki=d['ki'] if ki is None else ki,
                           a_force=a_force, f_exc=f_exc, nw=nw, degree=degree, noise=noise, c=c,
                           rho_d=rho_d, chi2_max=chi2_max, rel_max=rel_max, feedforward=feedforward,
                           conservative=conservative, excitation=excitation, v_app=v_app,
                           adaptive_exc=adaptive_exc, sector_lo=sector_lo, kgood_rule=kgood_rule,
                           edge_rule=edge_rule, recovery=recovery, follow_tau=follow_tau, ff_band=ff_band,
                           kg_avg=kg_avg, rel_meas=rel_meas, kg_gamma=kg_gamma, gain_noise=gain_noise, meas_bound=meas_bound, meas_better=meas_better, recontact=recontact, enc_y=enc_y, ff_jump=ff_jump,
                           oracle_meas=oracle_meas,
                           oracle_fb=oracle_fb)
        if tag:
            self.name = tag

    # hooks overridden by the oracle
    def stiffness(self, o, n):
        p = self.params
        if p.get('oracle_fb'):
            self._stiffness_est(o, n)                  # keep the estimator (and k_good) running
            return max(float(o.log[n-1, 8]), K_MIN), True   # diagnosis only: true stiffness for feedback
        return self._stiffness_est(o, n)

    def _stiffness_est(self, o, n):
        p = self.params
        fit = common.window_fit(o.log, n, p['nw'], p['degree'], p['noise'], self.ns, self.eb,
                                self.afc, self.ndel, p['enc_y'])
        self.fit = fit
        if not np.isfinite(fit['k']) or fit['info'] <= 0:
            return self.kbar, False
        k_hat = min(max(fit['k'], K_MIN), K_MAX)
        sig = 1.0/math.sqrt(fit['info'])            # info = ||Y_perp||^2/noise^2
        self.k_hat = k_hat; self.rel = sig/k_hat
        precise = sig <= p['rel_max']*k_hat
        if precise:
            self.k_check = k_hat                     # latest precise estimate, used only as a check
        consistent = fit['chi2'] <= p['chi2_max']
        # unbiased estimate for the open-loop channels (feed-forward/excitation); updated before the
        # feedback-stiffness choice so that every variant (incl. the ablation) shares this logic
        if p['kgood_rule'] == 'filter':
            if precise and consistent:
                self._filter_kgood(k_hat, sig/k_hat)
        elif precise and (consistent or p['kgood_rule'] == 'precise'):
            # the measurement channel has its own precision requirement (property 4: the feed-forward loop tolerates
            # only about -17% error in the measurement stiffness), stricter than the feedback gate
            ok_meas = p['rel_meas'] is None or sig <= p['rel_meas']*k_hat
            if p['meas_bound'] is not None:
                # property 1: bias from drift outside the subspace is |dk/k| <= ||r_perp||/(k||Y_perp||); the residual energy
                # beyond the noise floor (c0 * dof * noise^2) bounds ||r_perp||, so the window's worst-case relative error is
                # rel*(1 + sqrt(dof*max(0, chi2 - c0))); accept it for the measurement channel only if that stays inside
                # the feed-forward band (ln 1.2, property 4)
                dof = max(p['nw']-p['degree']-2, 1)
                bound = (sig/k_hat)*(1.+math.sqrt(dof*max(0., fit['chi2']-p['meas_bound'])))
                within = bound <= math.log(1.2)
                # also accept when even the worst case of the new window is closer than the current gain's discrepancy
                better = p['meas_better'] and np.isfinite(self.k_good) and bound < abs(math.log(k_hat/self.k_good))
                ok_meas = ok_meas and (within or better)
            if ok_meas:
                self.k_good = k_hat
                self.f_set = max(float(o.log[n-1, 5]), .3)    # force at which the measurement stiffness was identified
                self.kg_rel = sig/k_hat                        # relative uncertainty of the measurement gain
                if p.get('oracle_meas') == 'stale':
                    self.k_good = max(float(o.log[n-1, 8]), 1.)   # diagnosis only: exact value, estimator timing
        elif precise and p['kgood_rule'] == 'raise' and np.isfinite(self.k_good) and k_hat > self.k_good:
            # measurement stiffness: under-estimation destabilises the feed-forward loop (property 4), so a
            # precise window may raise it even when inconsistent (stiffness rising inside the window); it may
            # only be lowered on a consistent window -- the same asymmetry as the feedback stiffness
            self.k_good = k_hat
        if not p['conservative']:
            return k_hat, True
        cand = min(max(k_hat+p['c']*sig, K_MIN), K_MAX)
        if self.first:
            if precise and consistent:
                self.first = False
                return cand, True
            return self.kbar, False
        if precise and consistent:                   # valid window
            return max(cand, (1-p['rho_d'])*self.kbar), True
        if not consistent and p['edge_rule'] == 'raise':
            # model mismatch: local stiffness is changing inside the window; under-estimation is the
            # dangerous direction, so raise to the conservative bound (larger when less certain)
            return max(cand, self.kbar), False
        if precise:                                  # model mismatch (legacy rule): raise to estimate
            return max(k_hat, self.kbar), False
        return self.kbar, False                      # consistent but too little information: hold

    def _filter_kgood(self, k_hat, rel):
        """measurement stiffness to the precision the feed-forward band demands (property 4).
        Scalar Kalman filter on log k over consistent windows: per-tick estimates come from windows overlapping by
        N-1 samples, so each counts as 1/N of an independent window (variance rel^2 N); random-walk process noise
        sets an averaging memory of kg_avg windows; an innovation beyond 3 sigma is a stiffness jump (ridge entry)
        and resets the filter to the new window. The output is centred in the asymmetric stable band
        [0.75, 1.2] of kappa_g: k_good = exp(log k) / sqrt(0.75*1.2)."""
        N = self.params['nw']; z = math.log(k_hat); R = rel*rel*N
        if not np.isfinite(self._lk):
            self._lk, self._Pk = z, R
        else:
            Pp = self._Pk+R/(self.params['kg_avg']*N)**2
            S = Pp+R; inn = z-self._lk
            if inn*inn > 9.*S:
                self._lk, self._Pk = z, R
            else:
                K = Pp/S; self._lk += K*inn; self._Pk = (1-K)*Pp
        self.k_good = math.exp(self._lk)/math.sqrt(0.75*1.2)

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
        self.q0 = None; self.q_prev = None; self.x_prev = 0.; self.k_check = np.nan; self.ff_gain = 1.; self._lk = np.nan; self._Pk = np.nan; self.f_set = np.nan; self.kg_rel = np.nan

    def step(self, o):
        p = self.params
        n, log = o.n, o.log
        F, t = log[n-1, 5], log[n-1, 0]
        sent = log[n-1, 11]
        e = FD-F
        # contact supervision
        # recovery: True = full supervisor (approach mode); 'gate' = only switch off the motion
        # feed-forward and excitation and freeze k_bar while lost (k -> 0 leaves the FF loop's stable
        # range of k/k_ref, see ff_loop.py); False = none
        lost = p['recovery'] == 'gate' and F < .3
        if self.mode == 0 and F >= .5:
            self.y_contact = float(sent)                     # last command issued while in contact
        if self.mode == 0 and F < .3 and p['recovery'] is True:
            self.mode = 1; self.inband = 0; self.x_prev = 0.
            # surface velocity at loss (command direction), from the motion tracker if present
            self.v_follow = max(0., float(getattr(self, 'v_cmd_motion', 0.))); self.t_loss = t
        if self.mode == 1:                                   # lost: approach
            if F >= .5:
                self.mode = 0; self.adm.reset(); self.q_prev = None; self.x_prev = 0.
            else:
                self._trace(t, 1)
                v = p['v_app']
                yc = getattr(self, 'y_contact', None)
                if p['recontact'] and yc is not None and sent < yc:
                    # return at the speed limit to the last command that was in contact, then approach slowly
                    return min(sent+self.VMAX*self.Tc, yc)+v*self.Tc
                if p['follow_tau'] > 0:
                    # keep following the last estimated retraction, decaying, plus the approach speed
                    v += self.v_follow*math.exp(-(t-self.t_loss)/p['follow_tau'])
                return sent+min(v, self.VMAX)*self.Tc
        if not lost:
            k_bar, _ = self.stiffness(o, n)
            self.kbar = k_bar
        du = self.adm.increment(e)
        # drift feed-forward (normal mode only)
        dq = 0.
        if p['feedforward'] and self.mode == 0 and not lost:
            q = self.drift_pred(o, n)
            if isinstance(q, tuple):
                # displacement-domain drift prediction s [m] (oracle surface, or a motion tracker):
                # the feed-forward command increment is -Delta s
                if self.q_prev is not None:
                    k_un0 = self.k_good if np.isfinite(self.k_good) else self.kbar
                    dq = k_un0*(q[1]-self.q_prev)
                    if p['ff_jump'] is not None and np.isfinite(self.k_check) and np.isfinite(self.k_good):
                        # property 4: a discrepancy beyond x ff_jump between the latest precise estimate and the
                        # measurement gain cannot come from window noise (rel <= rel_max); it signals a stiffness jump
                        # (e.g. ridge entry), kappa_g is then certainly outside the stable band, so the feed-forward
                        # loop is opened until the measurement gain is updated by an accepted window
                        r = self.k_check/self.k_good
                        if r > p['ff_jump'] or r < 1./p['ff_jump']:
                            dq = 0.
                    if p['ff_band'] is not None and np.isfinite(self.k_check) and np.isfinite(self.k_good):
                        # property 4: the feed-forward loop gain scales with g (kappa_g - 1); keep g|kappa_g-1|
                        # inside the stable band [lo, hi] using the latest precise estimate as kappa_g check
                        lo, hi = p['ff_band']; r = self.k_check/self.k_good
                        g = (hi-1.)/(r-1.) if r > hi else ((1.-lo)/(1.-r) if r < lo else 1.)
                        self.ff_gain = g; dq *= g
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
        x = p['a_force']*self.exc_gain*math.sin(2*math.pi*p['f_exc']*(t+self.Tc)) if (p['excitation'] and self.mode == 0 and not lost) else 0.
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
