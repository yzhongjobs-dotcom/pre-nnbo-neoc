"""Theory-specified memory filters for independent, non-overlapping D4 windows.

q is a variance RATE (ln(k)^2/s); Q=q*dt is a per-update variance.
No robustification, acceptance gate, overlap inflation, or parameter fitting.
"""
import math

import numpy as np

H = .192
S = .2
N = S*S*H
TAU_LONG = 2.
TAU_SHORT = .2
Q_RATE_LO = N/TAU_LONG**2
Q_RATE_HI = N/TAU_SHORT**2
V_REF = 20.  # mm/s
Q_DISTANCE = Q_RATE_LO/V_REF  # Q_lo/(H*v_ref), Q_lo=q_rate_lo*H


def probabilities(logp):
    p = np.exp(logp-np.max(logp))
    return p/p.sum()


def moments(w, m, p):
    mu = float(w@m)
    return mu, float(w@(p+(m-mu)**2))


class Kalman:
    def __init__(self, tau=TAU_LONG):
        self.q = N/tau**2
        self.m = self.P = None
        self.K = None

    def update(self, y, R, ds=0., dt=H):
        if self.m is None:
            self.m, self.P = y, R
        else:
            prior = self.P+self.q*dt
            self.K = prior/(prior+R)
            self.m += self.K*(y-self.m)
            self.P = (1.-self.K)*prior
        return self.m, self.P


class IMM:
    """Two genuine interacting modes, including conditional moment mixing."""
    def __init__(self):
        self.q = np.array([Q_RATE_LO, Q_RATE_HI])
        self.transition = np.array([[.98, .02], [.02, .98]])
        self.w = np.ones(2)/2
        self.means = self.variances = None
        self.m = self.P = None

    def update(self, y, R, ds=0., dt=H):
        if self.means is None:
            self.means = np.full(2, y)
            self.variances = np.full(2, R)
        else:
            joint = self.w[:, None]*self.transition
            predicted_w = joint.sum(axis=0)
            mixing = joint/predicted_w[None, :]
            mixed_m = mixing.T@self.means
            mixed_p = np.sum(mixing*(self.variances[:, None]+
                                      (self.means[:, None]-mixed_m)**2), axis=0)
            prior = mixed_p+self.q*dt
            innovation = y-mixed_m
            innovation_var = prior+R
            self.w = probabilities(np.log(predicted_w)-.5*(
                np.log(2*math.pi*innovation_var)+innovation**2/innovation_var))
            K = prior/innovation_var
            self.means = mixed_m+K*innovation
            self.variances = (1-K)*prior
        self.m, self.P = moments(self.w, self.means, self.variances)
        return self.m, self.P


class Lepski:
    def __init__(self):
        self.history = []
        self.scales = (1, 2, 4, 8, 16)
        self.selected = 1
        self.m = self.P = None

    def update(self, y, R, ds=0., dt=H):
        self.history.append((y, R))
        self.history = self.history[-16:]
        candidates = []
        for n in self.scales:
            if n > len(self.history):
                break
            v = 1./sum(1./r for _, r in self.history[-n:])
            m = v*sum(x/r for x, r in self.history[-n:])
            # Nested precision means have Cov(mean_m,mean_N)=v_N.
            if all(abs(m-mm) <= 2.5*math.sqrt(max(vv-v, 0.))
                   for _, mm, vv in candidates):
                self.selected, self.m, self.P = n, m, v
            # Compare to ALL shorter means, even if a shorter scale failed.
            candidates.append((n, m, v))
        return self.m, self.P


class EWA:
    def __init__(self):
        self.taus = (.2, .5, 1., 2., 4.)
        self.q = N/np.array(self.taus)**2
        self.w = np.ones(5)/5
        self.means = self.variances = None
        self.m = self.P = None

    def update(self, y, R, ds=0., dt=H):
        if self.means is None:
            self.means = np.full(5, y)
            self.variances = np.full(5, R)
        else:
            prior = self.variances+self.q*dt
            innovation = y-self.means
            innovation_var = prior+R
            # Fixed share before the next likelihood is equivalent to sharing
            # after the previous posterior. Experts themselves are not mixed.
            predicted_w = .98*self.w+.02/len(self.w)
            self.w = probabilities(np.log(predicted_w)-.5*(
                np.log(2*math.pi*innovation_var)+innovation**2/innovation_var))
            K = prior/innovation_var
            self.means += K*innovation
            self.variances = (1-K)*prior
        self.m, self.P = moments(self.w, self.means, self.variances)
        return self.m, self.P


class Forgetting:
    def __init__(self, distance=False):
        self.distance = distance
        self.m = self.P = None
        self.lam = .998
        self.Q = None

    def update(self, y, R, ds=0., dt=H):
        if self.m is None:
            self.m, self.P = y, R
            return self.m, self.P
        self.Q = Q_DISTANCE*abs(ds) if self.distance else Q_RATE_LO*dt
        self.P += self.Q
        innovation = y-self.m
        innovation_var = self.P+R  # uninflated measurement variance
        K = self.P/innovation_var
        self.m += K*innovation
        self.P *= 1-K
        e = max(innovation**2/innovation_var-1., 0.)
        self.lam = max(.8, .998-(.998-.8)*e/(1+e))
        # Preserve the old update-then-forget sequence; returned covariance
        # includes the inflation for the next update (not a pure posterior).
        self.P /= self.lam
        return self.m, self.P


class HuangARD:
    """Scalar eqs 14--17, D4 state replacement, reference layer disabled.

All covariances E and lambda share the precision beta. R is the base window
variance, P=E/beta the actual variance. First window bootstraps the mean;
E_0=R (unit covariance in whitened observation coordinates). For learning,
a_0=1, b_0=1e-10 as in the paper; fixed precision uses beta=1.

The exact positive-evidence solution avoids arbitrary iteration counts:
Gamma=R+E_prev, lambda=(beta_prev*nu^2-Gamma)_+.
Active branch: beta_new=beta_prev. Inactive branch: lambda=0 and
beta_new=(a_prev+.5)/(b_prev+.5*nu^2/Gamma). This solves eq14 and eq15
jointly, with a,b incremented ONCE per observation, not per inner iteration.
"""
    def __init__(self, learn=False):
        self.learn = learn
        self.a, self.b = 1., (1e-10 if learn else 1.)
        self.beta = self.a/self.b
        self.m = self.P = self.E = None
        self.active = False
        self.alpha = 0.
        self.last_residual = 0.

    def update(self, y, R, ds=0., dt=H):
        if self.m is None:
            self.m, self.E = y, R
            self.P = self.E/self.beta
            return self.m, self.P
        nu = y-self.m
        Gamma = R+self.E
        self.alpha = max(self.beta*nu**2-Gamma, 0.)  # scaled lambda
        self.active = self.alpha > 0.
        prior = self.E+self.alpha
        if self.learn:
            self.a += .5
            self.b += .5*nu**2/(R+prior)
            self.beta = self.a/self.b
        Z = self.alpha*Gamma/(self.alpha+Gamma)
        delta = self.alpha*nu/(self.alpha+Gamma)
        rhs = Z+self.beta*delta**2
        self.last_residual = abs(self.alpha-rhs)/max(1., self.alpha, rhs)
        K = prior/(prior+R)
        self.m += K*nu
        self.E = prior*R/(prior+R)
        self.P = self.E/self.beta
        return self.m, self.P


class BlockHuang:
    def __init__(self, learn=False):
        self.core = HuangARD(learn=learn)
        self.buffer = []
        self.m = self.P = None

    def update(self, y, R, ds=0., dt=H):
        self.buffer.append((y, R))
        if self.m is None:
            # Causal bootstrap hold until the first complete four-window block.
            self.m, self.P = y, R
        if len(self.buffer) == 4:
            mean = sum(x for x, _ in self.buffer)/4
            var = sum(r for _, r in self.buffer)/16
            self.m, self.P = self.core.update(mean, var, ds, 4*dt)
            self.buffer.clear()
        return self.m, self.P


class Window:
    def update(self, y, R, ds=0., dt=H):
        self.m, self.P = y, R
        return y, R


FACTORIES = dict(
    off_window24=Window,
    kal_long=lambda: Kalman(2.), kal_short=lambda: Kalman(.2),
    imm=IMM, lepski=Lepski, ewa=EWA,
    forget_time=Forgetting, forget_dist=lambda: Forgetting(distance=True),
    ard_fixed=HuangARD, ard_learn=lambda: HuangARD(learn=True),
    ard_block=lambda: BlockHuang(), ard_block_learn=lambda: BlockHuang(learn=True))


def attach(pol, kind):
    """Finite/info guard only; >=24 actual control ticks between used windows.

Raw D4 still executes at every host request (preserves its fit/side effects).
No raw estimate leaks into a held output after the causal bootstrap. ds is
integrated from the observed scan coordinate across ALL intervening ticks.
"""
    est = pol.estimator
    original = est.estimate
    flt = FACTORIES[kind]()
    state = dict(last=None, out=None, x=None, ds=0., updates=[],
                 invalid=0, max_hold_ticks=0)

    def estimate(o, n, policy):
        k, sig, chi = original(o, n, policy)
        x = float(o.log[n-1, 1])
        if state['x'] is not None:
            # Requests normally occur on each tick, but integrate any gap too.
            a = state['x'][0]
            state['ds'] += float(np.abs(np.diff(o.log[a:n, 1])).sum())*1e3
        state['x'] = (n-1, x)
        valid = (np.isfinite(k) and k > 0 and np.isfinite(sig) and sig > 0
                 and est.fit.get('info', 0.) > 0)
        due = state['last'] is None or n-state['last'] >= 24
        if valid and due:
            elapsed = H if state['last'] is None else (n-state['last'])*pol.Tc
            state['out'] = flt.update(math.log(max(k, 200.)), (sig/k)**2,
                                      state['ds'], elapsed)
            state['last'] = n
            state['ds'] = 0.
            state['updates'].append(n)
        elif not valid:
            state['invalid'] += 1
        if state['out'] is None:
            return k, sig, chi
        state['max_hold_ticks'] = max(state['max_hold_ticks'], n-state['last'])
        m, P = state['out']
        # Same numerical output range as the host; no change to the latent state.
        kf = math.exp(min(max(m, math.log(200.)), math.log(3e5)))
        return kf, kf*math.sqrt(max(P, 0.)), chi

    est.estimate = estimate
    flt.sampling = state
    return flt
