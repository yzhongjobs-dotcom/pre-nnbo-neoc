"""Contact-model candidates on one host controller, each adapted to the position interface as the source paper's physics implies.

Only the contact relation differs. Every dynamic form identifies its single parameter from the same online quantity, the phase phi of
the probe response, obtained by synchronous detection after the window prior M, so the force filter cancels. Rules shared by all:
measured force is never a regressor; only the speed-determined part of friction is used, so it is not collinear with the contact force.

  L    F = k y, the existing relation
  W    weak power law, five recursive branches fitted to (j w)^beta over .3-100 Hz, beta from phi through the discrete response at the
       probe frequency, so the response can carry any phase the probe reports up to the .5 rad limit
  S    Kelvin-Boltzmann, one branch at tau = 1/w where its phase range is largest, alpha from phi
  H    Hunt-Crossley, n = 1.26, phi attributed to the rate term, tau_h = tan(phi)/w, k(F) ~ F^.206 in the conversion
  Hz   T-RO normal law, k(F) ~ F^(1/3) in the conversion
  AS   AuSoScan: sliding pressure term m(v) of eqs. 14-15 with Table I k1, k2, friction mu F + eq. 16 resistance, F/T identifies mu
  TR   T-RO: Hz plus friction mu F [1 + (2 nu - 1) 3 a^2/(10 R^2)] + k_d v, F/T identifies mu and k_d
  LF   L plus the measured tangential force
  WF, WA, WT  W with the sliding part of LF, AS, TR
  HA, HT      H with the sliding part of AS, TR
Two hosts: the relaxed prior, and the tuned projection of the existing method (names 'P5:x').
Friction reaches the axial force through theta and the arm's cross compliance. That coupling c is unknown to every candidate and is
identified online by the same rule, only from sliding-state changes.
"""
import math
import time

import numpy as np
from scipy.optimize import nnls

import sweeps as S
import common
import physics_plant2 as PP2
from contact_compare import ModelPolicy, CTL, RIDGE, P5
from estimators import ProjectionEstimator
from ridge_estimator import RidgeProjection

OBS_COLS = (0, 1, 2, 5, 11, 12, 13)
FORMS = ('L', 'W', 'S', 'H', 'Hz', 'AS', 'TR', 'LF', 'WF', 'WA', 'WT', 'HA', 'HT')
NAMES = FORMS+tuple('P5:'+f for f in ('L', 'W', 'S', 'H', 'Hz', 'AS', 'TR', 'HA'))
OMEGA = 2.*math.pi*20.8
TC = .008
R, NU, K1, K2 = .03, .45, .001, .003
SLS_TAU = 1./OMEGA
POWER_TAU = np.array([.001, .004, .016, .064, .256])
MEMORY = ('W', 'S', 'WF', 'WA', 'WT')
HC = ('H', 'HA', 'HT')
SLIDING = ('AS', 'TR', 'LF', 'WF', 'WA', 'WT', 'HA', 'HT')
GAMMA = {'H': .206, 'HA': .206, 'HT': .206, 'Hz': 1/3, 'TR': 1/3, 'WT': 1/3}
PHASE_USERS = ('W', 'S', 'H', 'WF', 'WA', 'WT', 'HA', 'HT')


class Measured2:
    """Exposes t, encoder x and y, axial force, command, tangential force and moment only."""
    direct = True

    def __init__(self, inner):
        self.inner = inner

    def reset(self, platform):
        self.inner.reset(dict(PP2.PLATFORM))

    def __call__(self, o):
        a = np.full((o.n, PP2.LOGC), np.nan)
        a[:, OBS_COLS] = o.log[:o.n, OBS_COLS]
        return self.inner(PP2.OLD.Obs(a, o.n))


def branch(tau):
    e = np.exp(-TC/tau); b = tau/TC*(1-e); z = np.exp(-1j*OMEGA*TC)
    return e, b, b*(1-z)/(1-e*z)


def _power_table():
    """Branch weights of a power law of slope beta, flat below .3 Hz, and the phase each gives at the probe through the discrete branches."""
    f = np.geomspace(.3, 100., 60); w = 2*np.pi*f
    A = np.array([1j*w*t/(1+1j*w*t) for t in POWER_TAU]).T
    _, _, h = branch(POWER_TAU)
    betas = np.arange(0., .6001, .01); weights = []; phases = []
    for beta in betas:
        G = (f/.3)**beta*np.exp(.5j*beta*np.pi*(1-np.exp(-f/.3))); wt = 1/np.abs(G)
        a, _ = nnls(np.r_[(A.real.T*wt).T, (A.imag.T*wt).T], np.r_[(G.real-1)*wt, G.imag*wt])
        weights.append(a); phases.append(float(np.angle(1+np.sum(a*h))))
    return np.array(weights), np.array(phases)


POWER_WEIGHTS, POWER_PHASES = _power_table()


def memory_from_phase(kind, phi):
    """Branch weights whose response has phase phi at the probe frequency, under this form's assumption."""
    if kind == 'S':
        tau = np.array([SLS_TAU]); _, _, h = branch(tau)
        t = math.tan(phi); den = float(h[0].imag-t*h[0].real)
        return np.array([min(max(t/den, 0.), 20.)]) if den > 1e-9 else np.zeros(1), tau
    a = np.array([np.interp(phi, POWER_PHASES, POWER_WEIGHTS[:, q]) for q in range(len(POWER_TAU))])
    return a, POWER_TAU


class Preprocess:
    """y through the form's response, axial force less the sliding part; both causal, rows up to n-1."""
    def _setup(self, kind, policy):
        self.kind = kind; self.owner = policy; self.private = None; self.count = 0
        self.st = np.zeros(1 if kind == 'S' else len(POWER_TAU)); self.y_last = None; self.y_out = None; self.y0 = None

    def _log(self, o, n):
        kind = self.kind
        if kind not in MEMORY and kind not in SLIDING and kind not in HC:
            return o.log
        if self.private is None:
            self.private = np.zeros((4096, o.log.shape[1]))
        own = self.owner
        for j in range(self.count, n):
            row = o.log[j].copy(); y = float(row[2])
            if self.y_last is None:
                self.y_last = y; self.y0 = y; self.y_out = y
            else:
                dy = y-self.y_last; self.y_last = y
                if kind in MEMORY:
                    self.st = own.e_y*self.st+own.a_y*own.b_y*dy
                    self.y_out = self.y0+((y-self.y0)+float(np.sum(self.st)))/own.mag_y
                elif kind in HC:
                    self.y_out = self.y0+((y-self.y0)+own.tau_h*dy/TC)/own.mag_h
                else:
                    self.y_out = y
            row[2] = self.y_out
            if kind in SLIDING:
                row[5] = row[5]-own.cf[j]
            self.private[j] = row
        self.count = n
        return self.private[:n]


class ContactEstimator(Preprocess, RidgeProjection):
    def __init__(self, kind):
        RidgeProjection.__init__(self, memory=False, **RIDGE)
        self.kind = kind

    def reset(self, policy):
        RidgeProjection.reset(self, policy)
        self._setup(self.kind, policy)


class ContactProjection(Preprocess, ProjectionEstimator):
    """The existing tuned projection, reading the same preprocessed signals."""
    def __init__(self, kind):
        ProjectionEstimator.__init__(self, **P5)
        self.kind = kind; self.cache = {}

    def reset(self, policy):
        ProjectionEstimator.reset(self, policy)
        self._setup(self.kind, policy)

    def _matrix(self, m):
        if m not in self.cache:
            t = np.linspace(-1, 1, m); Q, _ = np.linalg.qr(np.vander(t, self.degree+1))
            self.cache[m] = (np.eye(m)-Q@Q.T, max(m-self.degree-2, 1))
        return self.cache[m]

    def estimate(self, o, n, policy):
        class View:
            pass
        v = View(); v.log = self._log(o, n); v.n = n
        return ProjectionEstimator.estimate(self, v, n, policy)


class CandidatePolicy(ModelPolicy):
    def __init__(self, kind, estimator, **kwargs):
        self.kind = kind
        super().__init__(estimator, **kwargs)

    def reset(self, platform):
        self.cf = np.zeros(4096); self.done = 0; self.c_t = 0.; self.fslow = None; self.f_last = None; self.mu_hat = .25; self.kd_hat = 0.
        self.c_hat = 0.; self.sxx = .5; self.sxy = 0.; self.P = np.diag([1., 100.]); self.fv_prev = None; self.cf_filter = 0.; self.m_prev = 1.
        self.phi = 0.; self.tau_h = 0.; self.phi_count = 0; self.phi_w = 0.; self.phi_s = 0.; self.v_prev = 0.
        super().reset(platform)
        self.decay = (1-self.afc)**self.ns
        self.set_memory(0.)
        if self.kind in MEMORY:
            self.zalpha, self.ztau = self.a_y, memory_from_phase(self.kind, 0.)[1]
            self.q = np.zeros(len(self.ztau)); self.e, self.b = self.e_y, self.b_y; self.r_c = 1.
        else:
            self.zalpha = None; self.q = np.zeros(0); self.r_c = 1.

    def set_memory(self, phi):
        self.phi = phi; self.tau_h = 0.; self.mag_h = 1.
        if self.kind in HC:
            diff = (1.-np.exp(-1j*OMEGA*TC))/TC        # the backward difference, not j omega
            t = math.tan(phi); self.tau_h = t/max(diff.imag-t*diff.real, 1.)
            self.mag_h = abs(1.+self.tau_h*diff)
        if self.kind in MEMORY:
            a, tau = memory_from_phase(self.kind, phi)
            e, b, h = branch(tau)
            self.a_y, self.e_y, self.b_y = a, e, b
            self.mag_y = abs(1+np.sum(a*h))
            self.zalpha = a; self.e = e; self.b = b; self.r_c = self.mag_y

    def update_phase(self, o, n):
        """Synchronous detection of the probe response through the window prior, on raw measurements."""
        est = self.estimator; p = self.params
        j0 = max(self.ndel+1, n-1-est.nw); m = n-2-j0+1
        if m < est.nw:
            return
        _, Y, z, keep = common._build(o.log, j0, m, 0, self.ns, self.eb, self.afc, self.ndel, o.log[n-1, 0], o.log[n-1, 2], est.nw*self.ns*1e-3, p['enc_y'])
        if not keep.all():
            return
        M, _ = est._matrix(m)
        carrier = np.exp(-1j*OMEGA*o.log[j0:j0+m, 0])
        yc = carrier@(M@Y); zc = carrier@(M@z)
        if abs(yc) < 1e-12:
            return
        self.phi_count += 1
        lam = math.exp(-TC/2.)
        self.phi_w = lam*self.phi_w+1.; self.phi_s = lam*self.phi_s+min(max(float(np.angle(zc/yc)), 0.), .5)
        self.set_memory(self.phi_s/self.phi_w)

    def sliding(self, o, n):
        """Force the sliding state predicts, cumulative and at the measured filter, for rows up to n-1. Row n-1 has no tangential force yet."""
        kind = self.kind
        while self.done < n:
            j = self.done; self.done += 1
            F = float(o.log[j, 5])
            if self.fslow is None:
                self.fslow = F
            self.fslow += (1-math.exp(-TC/1.))*(F-self.fslow)
            v = float((o.log[j, 1]-o.log[j-1, 1])/TC) if j > 0 else 0.
            s = math.tanh(v/.002); fx = float(o.log[j-1, 12]) if j > 1 else float('nan')
            k = self.k_good if np.isfinite(self.k_good) else None
            d0 = self.fslow/k if k else None
            Rp = getattr(self, 'R_probe', R)       # probe geometry is known to the Method
            a = math.sqrt(Rp*d0) if d0 and d0 > 0 else None
            if abs(v) > .002 and np.isfinite(fx) and self.fslow > .5:
                x = np.array([self.fslow*s, v if kind in ('TR', 'WT', 'HT') else 0.])
                P = self.P/math.exp(-TC/2.)
                g = P@x/(1+x@P@x); e = fx-(self.mu_hat*x[0]+self.kd_hat*x[1])
                self.mu_hat = float(np.clip(self.mu_hat+g[0]*e, 0., 1.)); self.kd_hat = float(np.clip(self.kd_hat+g[1]*e, 0., 50.))
                self.P = P-np.outer(g, x@P)
                self.P = (self.P+self.P.T)/2
                self.P[0, 0] = min(self.P[0, 0], 1.); self.P[1, 1] = min(self.P[1, 1], 100.)
                lim = .95*math.sqrt(self.P[0, 0]*self.P[1, 1]); self.P[0, 1] = self.P[1, 0] = min(max(self.P[0, 1], -lim), lim)
            pressure = 0.
            if kind in ('AS', 'WA', 'HA'):
                fv = self.mu_hat*self.fslow*s+4*self.fslow*(K1+K2)*v/(3*math.pi*Rp)
                if a:
                    pressure = self.fslow*(K2-K1)*abs(v)/(2*a)
            elif kind in ('TR', 'WT', 'HT'):
                fv = self.mu_hat*self.fslow*(1+(2*NU-1)*3*(a or 0.)**2/(10*Rp*Rp))*s+self.kd_hat*v
            else:
                # measured tangential force without the part that only follows the contact force
                fv = (fx-self.mu_hat*s*(F-self.fslow)) if np.isfinite(fx) else 0.
            if self.fv_prev is None:
                self.fv_prev = fv; self.p_prev = pressure; self.cf[j] = 0.; self.f_last = F; self.y_last_s = float(o.log[j, 2])
                continue
            df = fv-self.fv_prev; self.fv_prev = fv
            y = float(o.log[j, 2]); de = F-self.f_last-(k or 0.)*(y-self.y_last_s); self.f_last = F; self.y_last_s = y
            threshold = .1 if kind in ('LF', 'WF') else .02     # three sigma of the measured channel's tick difference; the model channels have no noise
            changing = abs(v-self.v_prev) > .0005; self.v_prev = v
            if changing and abs(df) > threshold:
                self.sxx = .95*self.sxx+df*df; self.sxy = .95*self.sxy+df*de
                self.c_hat = float(np.clip(self.sxy/self.sxx, -.5, .5))
            self.c_t += self.c_hat*df+(pressure-self.p_prev); self.p_prev = pressure
            if kind in ('LF', 'WF'):
                self.cf_filter = self.c_t          # the tangential channel already carries the force filter
            else:
                self.cf_filter = self.decay*self.cf_filter+(1-self.decay)*self.c_t
            self.cf[j] = self.cf_filter

    def _stiffness_est(self, o, n):
        if self.kind in SLIDING:
            self.sliding(o, n)
        out = super()._stiffness_est(o, n)
        if self.kind in PHASE_USERS and np.isfinite(getattr(self, 'rel', np.nan)) and self.rel <= self.params['rel_max']:
            self.update_phase(o, n)
        return out

    def _next_z(self):
        o, n = self._ctx; j = n-1
        F = float(o.log[j, 5]); fin = (F-self.decay*o.log[j-1, 5])/(1-self.decay)
        y0, y1 = o.log[j-1, 2], o.log[j, 2]
        y = float(y0+self._encoder_fraction*(y1-y0))
        dF = fin-self._fin_prev; dy = y-self._y_prev
        if self.kind in SLIDING:
            # same correction as the regression, taken through the same one-step filter inverse
            dF -= (self.cf[j]-self.decay*self.cf[j-1])/(1-self.decay)-(self.cf[j-1]-self.decay*self.cf[j-2])/(1-self.decay) if j >= 2 else 0.
        if min(F, o.log[j-1, 5]) <= .5:
            self.q = np.zeros_like(self.q); self.dc_prev = 0.
        kc = self.k_good
        g = GAMMA.get(self.kind)
        if g and np.isfinite(kc) and np.isfinite(getattr(self, 'f_set', np.nan)):
            kc = kc*(max(F, .3)/self.f_set)**g
        if not np.isfinite(kc) or kc <= 0:
            kc = self.kbar
        if self.kind in MEMORY:
            dFe = (dF-float(np.sum((self.e-1)*self.q)))/(1+float(np.sum(self.zalpha*self.b)))
            self.q = self.e*self.q+self.zalpha*self.b*dFe
            dc = dFe*self.r_c/kc
        else:
            dc = dF/kc
        if self.kind in HC:
            a = self.tau_h/TC
            dc = (dc*self.mag_h+a*self.dc_prev)/(1+a); self.dc_prev = dc
        self.Zd += dc-dy; self.Z = self.Zd
        return self.Zd


def build(name, overrides=None):
    p = dict(S.FINAL_C); p.update(CTL); p.update(overrides or {})
    if name.startswith('P5:'):
        kind = name[3:]
        pol = CandidatePolicy(kind, ContactProjection(kind), single_model=False, **p)
        pol.centre = False; pol.acquire = False
        return pol
    return CandidatePolicy(name, ContactEstimator(name), single_model=False, **p)


def run_one(name, sc, seed=0, sensor=True):
    pol = build(name)
    t0 = time.perf_counter()
    log = PP2.episode(sc, Measured2(pol), seed, sensor=sensor)
    return log, pol, time.perf_counter()-t0
