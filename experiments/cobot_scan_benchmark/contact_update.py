"""Development C1: profile contact drift out of stiffness identification.

Only measured force, encoder position and sent commands enter this policy.
The old C1 and its published task results remain unchanged. This implementation
uses the same robot predictor and objective, with a different estimator and a
rule-based, command-continuous contact recovery. No physical-source recovery or
closed-loop safety theorem is implied by the projection.
"""
import math
import time

import numpy as np
from numba import njit

import plant as P
from candidate1_model import Candidate1, _sync, _solve_ord, XLO, XHI, LDMIN, LDMAX


@njit(cache=True)
def profile_contact(log, n, nw, lam, ndel, eb, afc, ns, kprior,
                    relative_prior, noise, degree):
    """Weighted variable projection; no arrival prior on nuisance drift.

    z = k*Y + B*beta. Eliminate beta before fitting k. The Schur complement
    I_k = ||(I-P_B) sqrt(W)Y||^2/noise^2 is data information, not a prior-
    inflated confidence. A finite prior regularizes weak windows. Polynomial
    coefficients are force units; convert to equivalent velocity afterwards.
    """
    j0 = max(ndel + 1, n - 1 - nw)
    m = n - 2 - j0 + 1
    if m < degree + 3:
        return np.array([kprior, log[n-1, 5], 0., 0., 0., 0., 0.])
    te, ye = log[n-1, 0], log[n-1, 2]
    tw = nw * ns * .001
    a = 1. - afc
    X = np.zeros((m, degree + 1))
    Y = np.zeros(m)
    Z = np.zeros(m)
    valid = 0
    for jj in range(m):
        j = j0 + jj
        if min(log[j, 5], log[j+1, 5]) <= .5:
            continue
        valid += 1
        wt = math.sqrt(lam ** (n-2-j))
        yy = log[j, 2]
        tgt = log[j-ndel, 11]
        for mm in range(1, ns+1):
            yy = tgt + (yy-tgt)*eb
            w = afc*a**(ns-mm)*wt
            tau = (log[j, 0] + mm*.001-te)/tw
            Y[jj] += w*(yy-ye)
            for d in range(degree+1):
                X[jj, d] += w*tau**d
        Z[jj] = wt*(log[j+1, 5] - a**ns*log[j, 5])
    if valid < degree + 3:
        return np.array([kprior, log[n-1, 5], 0., 0., 0., 0., float(valid)])
    gram = X.T @ X
    for d in range(degree+1):
        gram[d, d] += 1e-12
    cy = np.linalg.solve(gram, X.T @ Y)
    cz = np.linalg.solve(gram, X.T @ Z)
    yr = Y-X@cy
    zr = Z-X@cz
    info = np.dot(yr, yr)/(noise*noise)
    prior_info = 1./(relative_prior*kprior)**2
    k = (np.dot(yr, zr)/(noise*noise) + prior_info*kprior)/(info+prior_info)
    k = min(max(k, 200.), 3e5)
    beta = cz-k*cy
    residual = Z-k*Y-X@beta
    rms = math.sqrt(np.dot(residual, residual)/max(valid-degree-2, 1))
    # Convert force polynomial to length drift. This is an equivalent contact
    # state; its value need not equal anatomical surface motion.
    vs = beta[1]/tw/k
    acc = 2*beta[2]/tw**2/k if degree >= 2 else 0.
    return np.array([k, beta[0], vs, acc, info, rms, float(valid)])


class ContactUpdate(Candidate1):
    def __init__(self, relative_prior=.3, degree=3, noise=.04,
                 approach=.006, recovery_damping=40., recovery_cp=.5,
                 projection=True, reference=True, motion_threshold=.004,
                 gate_release=.128, **kw):
        kw.pop('fb', None)
        for key in ('sup','rec_D','F_off','band_lo','band_hi','n_back',
                    'n_settle','u_rate','aw','vs_max','F_fb_hi'):
            kw.pop(key, None)
        kw['solver'] = 'ord'
        super().__init__(**kw)
        self.relative_prior = relative_prior
        self.degree, self.noise = degree, noise
        self.approach = approach
        self.recovery_damping, self.recovery_cp = recovery_damping, recovery_cp
        self.projection, self.reference = projection, reference
        self.motion_threshold,self.gate_release = motion_threshold,gate_release

    def reset(self, platform):
        super().reset(platform)
        self.mode = 0
        self.rc = 0
        self.inband = 0
        self.Ioff = 0.
        self.motion_before_loss = 0.
        self.details = []
        self.information = 0.
        self.project_until = -1.

    def _next_admittance(self, F, D, KP, ki, uff):
        ef = P.FD-F
        In = np.clip(self.rs[2]+ef*self.Tc, -20., 20.)
        vv = np.clip(self.rs[1]*math.exp(-D*self.Tc/self.M)
                     + (1-math.exp(-D*self.Tc/self.M))*(ef+ki*In)/D,
                     -self.VMAX, self.VMAX)
        return self.rs[0]+.5*(self.rs[1]+vv)*self.Tc

    def __call__(self, o):
        started = time.perf_counter()
        n, log = o.n, o.log
        if n == 0:
            D, KP = self.rule(self.kuse)
            self.x = np.array([math.log(D), KP/1e-4, 0., 0.])
            self.g = (self.M, D, 0., self.ki, KP, 0.)
            return self.g
        if self.rs is None:
            self.rs = np.array([log[0, 2], 0., 0., log[0, 2]])
        _sync(log, self.m, n, self.rs, *self.g[:4], self.Tc, self.VMAX)
        self.m = n
        F, y = log[n-1, 5], log[n-1, 2]
        if self.mode == 0 and F < .3:
            self.mode = 1
            self.motion_before_loss = np.clip(self.th[2], -.04, .04)
            self.Ioff = self.rs[2]
            self.inband = 0
        if self.mode:
            # Keep the contact-stiffness rule at re-contact. Fixed gains tuned
            # on a different tissue are not used. During loss close relative
            # to the last estimated drift, rather than removing all damping.
            D = np.clip(self.kuse/self.recovery_damping, 100., P.GMAX[1])
            KP = min(self.recovery_cp/self.kuse, P.GMAX[4])
            if F < .5:
                v = np.clip(-self.motion_before_loss+self.approach, -.04, .06)
                target = log[n-1, 11]+v*self.Tc
                self.inband = 0
            else:
                # A force-level correction advances the commanded reference,
                # so the first touch cannot introduce a proportional jump.
                v = np.clip((P.FD-F)/(self.kuse*.04), -.03, .03)
                target = log[n-1, 11]+v*self.Tc
                self.inband = self.inband+1 if 1. <= F <= 6. else 0
            target = np.clip(target, y-.003, y+.003)
            yn = self._next_admittance(F, D, KP, 0., 0.)
            uff = float(np.clip(target-yn-KP*(P.FD-F), -.02, .02))
            self.g = (self.M, D, 0., 0., KP, uff)
            self.th[:] = [math.log(self.kuse), F, 0., 0.]
            self.x = np.array([math.log(D), KP/1e-4, uff/1e-3, 0.])
            if self.inband >= 12:
                self.mode = 0
                self.rc = n-1
                self.Ioff = self.rs[2]
            self._record(o, started, 0., 0.)
            return self.g

        nw = self.Nw if self.rc == 0 else min(self.Nw, n-1-self.rc)
        legacy = None
        if not self.projection or self.projection == 'adaptive':
            from candidate1_model import _est
            legacy = self.th.copy()
            mu = legacy.copy()
            mu[2] += legacy[3]*self.Tc
            pw = np.array([1/self.sk**2,0.,1/self.sv**2,1/self.sa**2])
            info = np.zeros(2)
            valid = _est(log,n,nw,self.lam,self.ndel,self.eb,self.afc,self.ns,
                         legacy,mu,pw,self.nfree,self.est_lo,self.est_hi,
                         self.iters_est,self.sz,info)
            if valid == 0:
                legacy[:] = [mu[0],F,0.,0.]
            if abs(legacy[2]) > self.motion_threshold:
                self.project_until = o.t[-1]+self.gate_release
        project = bool(self.projection) and (self.projection != 'adaptive' or o.t[-1]<self.project_until)
        if project:
            fit = profile_contact(log, n, nw, self.lam, self.ndel,
                                  self.eb, self.afc, self.ns, self.kuse,
                                  self.relative_prior, self.noise, self.degree)
            k = max(fit[0], self.kuse*math.exp(-self.Tc/self.trel)) if self.trel else fit[0]
            self.kuse = k
            self.th[:] = [math.log(k), np.clip(fit[1], -5., 30.),
                          np.clip(fit[2], -.1, .1), np.clip(fit[3], -5., 5.)]
            self.information = fit[4]
        else:
            self.th[:] = legacy
            self.kuse = max(math.exp(self.th[0]), self.kuse*math.exp(-self.Tc/self.trel))

        # Conditional integral: remove the integral accumulated before/during
        # recovery without resetting the actual plant's hidden state.
        In = np.clip(self.rs[2]+(P.FD-F)*self.Tc, -20., 20.)
        ki = self.ki if self.Ioff == 0. else float(np.clip(self.ki*(In-self.Ioff)/In, 0., self.ki)) if abs(In)>1e-9 else 0.
        I0 = ki/self.ki*In-(P.FD-F)*self.Tc if self.ki else 0.
        nd = self.ndel
        tgt = log[n-1-nd,11] if n-1-nd >= 0 else log[0,2]
        q = np.array([log[j,11] if j>=0 else log[0,2] for j in range(n-nd,n)])
        st = np.array([self.rs[0],self.rs[1],I0,self.rs[3],y,F,tgt,
                       self.th[1],self.kuse,self.th[2],self.th[3],y,0.,0.,0.,0.])
        self.dith_a = np.clip(self.Fdith/self.kuse,self.amin,self.amax) if self.Fdith else 0.
        dith = np.array([self._dither(n+i) for i in range(self.H)])
        D, KP = self.rule(self.kuse)
        prior = np.array([math.log(D),KP/1e-4,0.,0.])
        prev = self.x.copy()
        prev[2] = np.clip(prev[2]+prev[3]*1e-2*min(self.Tc,self.tu)/1e-3,XLO[2],XHI[2])
        lo, hi = XLO.copy(), XHI.copy()
        lo[0], hi[0] = max(LDMIN,prev[0]-self.rate_slow),min(LDMAX,prev[0]+self.rate_fast)
        lo[1], hi[1] = max(0.,prev[1]-self.rateK_fast),min(10.,prev[1]+self.rateK_slow)
        if not self.reference:
            lo[2:] = 0.; hi[2:] = 0.; prev[2:] = 0.
        x = np.clip(prev, lo, hi)
        info = np.zeros(4)
        cost = _solve_ord(x,self.iters,self.H,self.per,st,q,dith,self.pr,
                          self.scen,prev,prior,self.rw,lo,hi,self.Jc,self.rprev,info)
        self.x = x
        self.g = (self.M,math.exp(x[0]),0.,ki,x[1]*1e-4,x[2]*1e-3+dith[0])
        self._record(o, started, cost, info[0])
        return self.g

    def _record(self, o, started, cost, rollouts):
        self.t_call.append(time.perf_counter()-started)
        self.details.append([o.t[-1],self.kuse,self.th[2],self.th[3],
                             self.information,self.mode,*self.g])
        self.nroll.append(rollouts)
        if o.n % 4 == 0:
            self.hist.append((o.t[-1],math.exp(self.th[0]),self.kuse,
                              self.th[2],self.g[1],self.g[4],self.g[5],cost))
