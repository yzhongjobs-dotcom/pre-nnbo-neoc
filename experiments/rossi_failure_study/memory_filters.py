"""Memory mechanisms over window estimates of ln k: one interface, used offline and inside the closed loop.

Each filter receives per update: ell (ln of the window estimate), R (its variance, already inflated for window overlap), ds (mm moved since the last update),
usable (precise and consistent window). It returns (m, P), the current ln k and its variance.
  window        no memory beyond the window (the baseline)
  kalman        q = q_t + q_s*|ds|, Huber innovation; q_s=0 gives a fixed random walk
  lepski        intersection of confidence intervals over dyadic histories (Lepski/Goldenshluger-Nemirovski)
  forgetting    exponential forgetting with lambda driven by the normalized innovation (Fortescue-Kershenbaum type)
  sparse        hierarchical sparse Bayes on changes: Delta ell ~ N(0, gamma_t) with gamma_t learned online by EM from a prior scale,
                in the spirit of Huang et al. CMAME 2020 (mechanism only; not a reproduction of that paper's dual state-parameter filter)
"""
import math

import numpy as np


class Window:
    name = 'window'

    def __init__(self, **_):
        self.m = None; self.P = None

    def update(self, ell, R, ds, usable):
        if usable:
            self.m, self.P = ell, R
        return self.m, self.P


class Kalman:
    name = 'kalman'

    def __init__(self, q_t=0., q_s=1e-3, dt=.048, huber=2., **_):
        self.q_t, self.q_s, self.dt, self.huber = q_t, q_s, dt, huber
        self.m = None; self.P = None

    def update(self, ell, R, ds, usable):
        if self.m is None:
            if usable:
                self.m, self.P = ell, R
            return self.m, self.P
        self.P += self.q_t*self.dt+self.q_s*ds
        if usable:
            S = self.P+R; nu = ell-self.m
            r = R
            if self.huber and abs(nu) > self.huber*math.sqrt(S):
                r = max((abs(nu)/self.huber)**2-self.P, R); S = self.P+r
            G = self.P/S
            self.m += G*nu; self.P *= 1.-G
        return self.m, self.P


class Lepski:
    name = 'lepski'

    def __init__(self, scales=(1, 2, 4, 8, 16), c=2., stride=4, **_):
        self.scales, self.c, self.stride = scales, c, stride
        self.h = []; self.m = None; self.P = None

    def update(self, ell, R, ds, usable):
        if not usable:
            return self.m, self.P
        self.h.append((ell, R))
        self.h = self.h[-self.scales[-1]*self.stride:]
        lo, hi = -math.inf, math.inf
        best = None
        for n in self.scales:
            part = self.h[-n*self.stride:] if n > 1 else self.h[-1:]
            if len(part) < (n*self.stride if n > 1 else 1):
                break
            w = np.array([1./p[1] for p in part]); x = np.array([p[0] for p in part])
            # overlap-inflated R already accounts for correlation of neighbouring windows
            mean = float(np.sum(w*x)/np.sum(w)); var = float(1./np.sum(w))
            half = self.c*math.sqrt(var)
            lo, hi = max(lo, mean-half), min(hi, mean+half)
            if lo > hi:
                break
            best = (mean, var)
        if best is not None:
            self.m, self.P = best
        return self.m, self.P


class Forgetting:
    name = 'forgetting'

    def __init__(self, lam_max=.995, lam_min=.8, gain=1., q_s=0., **_):
        self.lam_max, self.lam_min, self.gain, self.q_s = lam_max, lam_min, gain, q_s
        self.m = None; self.P = None

    def update(self, ell, R, ds, usable):
        self.pending = getattr(self, 'pending', 0.)+ds
        if not usable:
            return self.m, self.P
        if self.m is None:
            self.m, self.P = ell, R
            self.pending = 0.
            return self.m, self.P
        self.P += self.q_s*self.pending; self.pending = 0.
        S = self.P+R; nu = ell-self.m
        G = self.P/S; self.m += G*nu; self.P *= 1.-G
        # lambda falls when the normalized innovation exceeds its expected unit size
        e = max(nu*nu/S-1., 0.)
        lam = max(self.lam_min, self.lam_max-self.gain*(self.lam_max-self.lam_min)*e/(1.+e))
        self.P /= lam
        return self.m, self.P


class Sparse:
    name = 'sparse'

    def __init__(self, gamma0=1e-4, a=.5, b=1e-5, dt=.048, q_s=0., **_):
        self.gamma0, self.a, self.b, self.q_s = gamma0, a, b, q_s
        self.m = None; self.P = None; self.gamma = gamma0

    def update(self, ell, R, ds, usable):
        if self.m is None:
            if usable:
                self.m, self.P = ell, R
            return self.m, self.P
        if not usable:
            self.P += self.gamma0+self.q_s*ds
            return self.m, self.P
        nu = ell-self.m
        # EM for the change variance at this step, inverse-gamma hyperprior (a,b): a few passes of the local fixed point
        g = self.gamma
        for _ in range(3):
            Pp = self.P+g+self.q_s*ds; S = Pp+R; G = g/S
            post_mean, post_var = G*nu, g*(1.-G)
            g = (post_mean**2+post_var+2.*self.b)/(1.+2.*self.a+2.)
        Pp = self.P+g+self.q_s*ds; S = Pp+R; K = Pp/S
        self.m += K*nu; self.P = Pp*(1.-K); self.gamma = g
        return self.m, self.P


class HuangARD:
    """Parameter filter of Huang et al. CMAME 2020, eqs (10)-(17), scalar theta = ln k, on window estimates.

    The state filter is replaced by the D4 elimination inside the window (its covariance is in R). The reference-value layer (eqs 6, 18, 19) is off, as the paper
    does for load tracking (Sec. 4.2, variance set very large): contact stiffness varies both ways along the scan, so there is no sparse deviation from a reference.
    alpha is the change variance from the fixed point alpha = Z + (mu-mu^-)^2 (eq 14); rho is the noise precision a/b (eqs 15, 17), with a0=1, b0 small as in the paper.
    """
    name = 'huang_ard'

    def __init__(self, q_s=0., iters=8, learn_rho=True, **_):
        self.q_s, self.iters, self.learn_rho = q_s, iters, learn_rho
        self.m = None; self.P = None; self.a = 1.; self.b = 1.

    def update(self, ell, R, ds, usable):
        if self.m is None:
            if usable:
                self.m, self.P = ell, R
            return self.m, self.P
        rho = self.a/self.b if self.learn_rho else 1.
        P0 = self.P+self.q_s*ds
        if not usable:
            self.P = P0
            return self.m, self.P
        nu = ell-self.m; r = R/rho
        alpha = max(nu*nu-P0-r, 1e-12)
        for _ in range(self.iters):
            lam = P0+r
            Z = 1./(1./lam+1./alpha)
            mu_change = Z*nu/lam
            alpha = Z+mu_change**2
        Pm = P0+alpha; K = Pm/(Pm+r)
        self.m += K*nu; self.P = Pm*(1.-K)
        if self.learn_rho:
            self.a += .5; self.b += .5*nu*nu/(P0+alpha+R)
        return self.m, self.P


class ThresholdARD:
    """Huang's ARD with an explicit evidence threshold, derived from its closed form.

    With a flat hyperprior the fixed point of eq (14) is alpha* = max(nu^2-S, 0), so in a scalar low-SNR problem it fires at one sigma and leaves almost no memory.
    A spike-and-slab prior with change probability p moves the threshold to kappa^2 S, kappa^2 ~ 2 ln(1/p); here S uses the true window variance s2 = R/overlap, not the inflated R.
    The process model between changes is the physical prior q = q_t*dt + q_s*|ds|. The thresholds differ by direction because under-estimating is the unstable side:
    upward evidence moves the estimate sooner than downward evidence.
    """
    name = 'threshold_ard'

    def __init__(self, q_t=0., q_s=1e-4, dt=.048, kappa_up=2., kappa_down=3., overlap=4., **_):
        self.q_t, self.q_s, self.dt, self.ku, self.kd, self.ov = q_t, q_s, dt, kappa_up, kappa_down, overlap
        self.m = None; self.P = None

    def update(self, ell, R, ds, usable):
        if self.m is None:
            if usable:
                self.m, self.P = ell, R/self.ov
            return self.m, self.P
        self.P += self.q_t*self.dt+self.q_s*ds
        if not usable:
            return self.m, self.P
        s2 = R/self.ov
        nu = ell-self.m; S = self.P+s2
        kappa = self.ku if nu > 0 else self.kd
        alpha = max(nu*nu-kappa*kappa*S, 0.)
        Pm = self.P+alpha
        K = Pm/(Pm+R)      # between changes neighbouring windows overlap, so the update uses the inflated R
        if alpha > 0.:
            K = Pm/(Pm+s2)  # a detected change: this window alone carries the new level
        self.m += K*nu; self.P = Pm*(1.-K)
        return self.m, self.P


FILTERS = dict(threshold_ard=ThresholdARD, huang_ard=HuangARD, window=Window, kalman=Kalman, lepski=Lepski, forgetting=Forgetting, sparse=Sparse)


def attach(pol, kind, hop=6, overlap=None, raw_gate=True, a_scale=1., shape_fc=None, keep_sigma=False, **params):
    """Wrap a built policy's window estimator so that the Method uses (m, P) instead of the raw window, gates unchanged.

    The feed-forward jump check keeps reading the raw window. x is read from the observation log (col 1)."""
    est = pol.estimator; original = est.estimate
    flt = FILTERS[kind](**params)
    nw = getattr(est, 'nw', 24)
    overlap = nw/hop if overlap is None else overlap
    state = dict(tick=0, x=None, raw=None, out=None, ds=0.)
    p = pol.params
    if a_scale != 1.:
        p['a_force'] = p['a_force']*a_scale

    def estimate(o, n, policy, original=original):
        k, sig, chi = original(o, n, policy)
        if not np.isfinite(k) or k <= 0:
            return k, sig, chi
        k_ref = max(pol.k_good, .5*k) if getattr(pol, 'gate_on_prior', False) and np.isfinite(pol.k_good) else k
        if sig <= p['rel_max']*k_ref and k > 200.:
            state['raw'] = k     # the jump check keeps its own rule: the latest precise raw window
        x = float(o.log[n-1, 1])
        if state['x'] is not None:
            state['ds'] += abs(x-state['x'])*1e3
        state['x'] = x
        state['tick'] += 1
        if state['tick'] % hop == 0 or state['out'] is None:
            usable = chi <= p['chi2_max'] and (not raw_gate or sig <= p['rel_max']*k)
            m, P = flt.update(math.log(max(k, 200.)), (sig/k)**2*overlap, state['ds'], usable)
            state['ds'] = 0.
            state['out'] = (m, P)
        m, P = state['out']
        if m is None:
            return k, sig, chi
        kf = math.exp(m)
        if shape_fc is not None:
            # gain at the loop's upper crossover from the probe gain and online phase, first-order Bode, relative to the 14 deg design point the Method was selected at
            phi = float(getattr(pol, 'phi', math.radians(14.)))
            kf *= (shape_fc/20.8)**(2.*(phi-math.radians(14.))/math.pi)
        # keep_sigma: denoise the gain only. The window's relative sigma still sets the gates and the feedback margin k_bar = k + c*sigma that the loop was designed with;
        # letting the fused sigma into k_bar removes that margin (about 16%) and moves the loop towards the low-gain side.
        return kf, (kf*sig/k if keep_sigma else kf*math.sqrt(P)), chi

    est.estimate = estimate
    stiffness = pol._stiffness_est

    def stiffness_with_raw_check(o, n):
        out = stiffness(o, n)
        if state['raw'] is not None and np.isfinite(pol.k_check):
            pol.k_check = state['raw']
        return out

    pol._stiffness_est = stiffness_with_raw_check
    return flt


class HuangDual:
    """Both filters of Huang et al. run per control tick: theta = k with ARD on its change (eqs 10-17) and state x = D4-smooth d (eqs 20-26), cycled (Algorithm 1).

    Observation z = k*Y + x[0] + e (eq 2). x follows the integrated random walk that the D4 prior corresponds to, with state noise q learned (eq 24). Reference layer off.
    Scales: noise variance sigma2 comes from the Method's calibrated noise; the paper's learned precision multiplies it."""

    def __init__(self, sigma=.04, cycles=2, inner=4, q_s=0.):
        self.A = np.array([[1., 1., .5, 1./6.], [0., 1., 1., .5], [0., 0., 1., 1.], [0., 0., 0., 1.]])
        self.G = np.array([1./24., 1./6., .5, 1.])
        self.C = np.array([1., 0., 0., 0.])
        self.s2, self.cycles, self.inner, self.q_s = sigma**2, cycles, inner, q_s
        self.x = None; self.Px = None; self.k = None; self.Pk = None; self.q = 1e-4; self.a = 1.; self.b = 1.

    def start(self, k, sigk, z, Y):
        self.k, self.Pk = k, sigk**2
        self.x = np.array([z-k*Y, 0., 0., 0.]); self.Px = np.eye(4)*1.

    def step(self, z, Y, ds=0.):
        rho = self.a/self.b
        r = self.s2/rho
        xp = self.A@self.x; Px0 = self.A@self.Px@self.A.T
        Pk0 = self.Pk+self.q_s*ds
        k = self.k
        alpha = self.Pk*.01+1e-9; q = self.q
        cg = float(self.C@self.G)
        for _ in range(self.cycles):
            # state filter with the current k; q is the variance of the scalar driving noise, from the fixed point of eq (24)
            nu = z-k*Y-self.C@xp
            s0 = self.C@Px0@self.C+r+Y*Y*(Pk0+alpha)
            for _ in range(self.inner):
                S = s0+q*cg*cg
                w_mean = q*cg*nu/S; w_var = q*(1.-q*cg*cg/S)
                q = max(w_var+w_mean*w_mean, 1e-12)
            Pxm = Px0+q*np.outer(self.G, self.G)
            S = self.C@Pxm@self.C+r+Y*Y*(Pk0+alpha)
            K = Pxm@self.C/S
            x = xp+K*nu; Px = (np.eye(4)-np.outer(K, self.C))@Pxm
            # parameter filter with the current state prediction; alpha from eq (14)
            noise = r+self.C@Pxm@self.C
            nu_k = z-self.C@xp-self.k*Y
            # start the fixed point from the evidence of this sample, not from a small value it could not leave
            alpha = max((nu_k*nu_k-Y*Y*Pk0-noise)/max(Y*Y, 1e-30), Pk0*.01+1e-9)
            for _ in range(self.inner):
                lam = Y*Y*Pk0+noise
                Z = 1./(Y*Y/lam+1./alpha)
                mu = Z*Y*nu_k/lam
                alpha = max(Z+mu*mu, 1e-12)
            Pm = Pk0+alpha
            G = Pm*Y/(Y*Y*Pm+noise)
            k = self.k+G*nu_k; Pk = Pm*(1.-G*Y)
        s = Y*Y*(Pk0+alpha)+self.C@Pxm@self.C+self.s2
        self.a += .5; self.b += .5*(z-self.C@xp-self.k*Y)**2/s
        self.k, self.Pk, self.x, self.Px, self.q = float(k), float(Pk), x, Px, q
        return self.k, math.sqrt(max(self.Pk, 0.))


def attach_dual(pol, q_s=0., **params):
    """Per-tick HuangDual inside the Method. Uses the same preprocessed z_j, Y_j as the window (common._build); the window still supplies chi for the gates and starts the filter."""
    import common
    est = pol.estimator; original = est.estimate
    dual = HuangDual(sigma=pol.params['noise'], q_s=q_s, **params)
    state = dict(x=None, ye=None, raw=None, ds=0.)
    p = pol.params
    if a_scale != 1.:
        p['a_force'] = p['a_force']*a_scale

    def estimate(o, n, policy, original=original):
        k, sig, chi = original(o, n, policy)
        if not np.isfinite(k) or k <= 0:
            return k, sig, chi
        k_ref = max(pol.k_good, .5*k) if getattr(pol, 'gate_on_prior', False) and np.isfinite(pol.k_good) else k
        if sig <= p['rel_max']*k_ref and k > 200.:
            state['raw'] = k
        reg = est._log(o, n)
        j0 = n-2
        if state['ye'] is None:
            state['ye'] = float(reg[j0, 2])
        _, Y, z, keep = common._build(reg, j0, 1, 0, policy.ns, policy.eb, policy.afc, policy.ndel, reg[n-1, 0], state['ye'], est.nw*policy.Tc, p['enc_y'])
        x = float(o.log[n-1, 1])
        if state['x'] is not None:
            state['ds'] += abs(x-state['x'])*1e3
        state['x'] = x
        if dual.k is None:
            if sig <= p['rel_max']*k and chi <= p['chi2_max'] and keep[0]:
                dual.start(k, sig, float(z[0]), float(Y[0]))
            return k, sig, chi
        if keep[0]:
            kf, sf = dual.step(float(z[0]), float(Y[0]), state['ds']); state['ds'] = 0.
            kf = min(max(kf, 200.), 3e5)
            return kf, max(sf, 1e-3*kf), chi
        return dual.k, math.sqrt(max(dual.Pk, 0.)), chi

    est.estimate = estimate
    stiffness = pol._stiffness_est

    def stiffness_with_raw_check(o, n):
        out = stiffness(o, n)
        if state['raw'] is not None and np.isfinite(pol.k_check):
            pol.k_check = state['raw']
        return out

    pol._stiffness_est = stiffness_with_raw_check
    return dual
