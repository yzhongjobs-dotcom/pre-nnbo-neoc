"""Window stiffness with a time-varying first-order term (functional-series expansion of a time-varying parameter).

The current estimator (common.window_fit) solves z = k Y + B beta by variable projection, i.e. it assumes k constant
in the 128 ms window (H1). Across a stiffness boundary, or when the force-dependent (Hertz-type) stiffness changes
during a push, k changes inside the window: the constant fit returns a window average, is inconsistent (chi2 > max)
and the measurement gain stays at an old value. Expanding k(tau) = k_m + c1 tau_c in the window (tau_c the filtered
row time, centred so that tau_c = 0 at mid-window and +1/2 at the window end) gives
    z = k_m Y + c1 (tau_c o Y) + B beta + eps,
still linear in all unknowns; B is eliminated by the same projection. The stiffness at the window end is
k_e = k_m + c1/2 with variance from the 2x2 projected information matrix. c1 = 0 reproduces the current estimator.
No contact model is used here; the expansion is generic in the cause of the change (path or force).
"""
import math
import numpy as np


def tv_fit(X, Y, Z, noise):
    """X: drift basis (rows x (p+1)) incl. constant column 0 and tau column 1, Y displacement, Z filter-inverse force.
    Returns dict(k_e, k_m, c1, cov (2x2 of k_m, c1), sd_e, chi2, k_const, sd_const, chi2_const)."""
    out = dict(k_e=np.nan, k_m=np.nan, c1=np.nan, sd_e=np.inf, chi2=np.inf, k_const=np.nan, sd_const=np.inf, chi2_const=np.inf)
    m, q = X.shape
    if m < q+4:
        return out
    tau = X[:, 1]/np.where(X[:, 0] != 0, X[:, 0], 1.)
    tc = tau+.5
    Q, _ = np.linalg.qr(X)
    P = lambda v: v-Q@(Q.T@v)
    y1, y2, zr = P(Y), P(tc*Y), P(Z)
    # constant model (current estimator)
    syy = float(y1@y1)
    if syy <= 0:
        return out
    k0 = float(y1@zr)/syy; r0 = zr-k0*y1
    out.update(k_const=k0, sd_const=noise/math.sqrt(syy), chi2_const=float(r0@r0)/noise**2/max(m-q-1, 1))
    A = np.column_stack([y1, y2]); G = A.T@A
    if np.linalg.cond(G) > 1e10:
        return out
    th = np.linalg.solve(G, A.T@zr); r = zr-A@th
    cov = np.linalg.inv(G)*noise**2
    e = np.array([1., .5])
    out.update(k_m=float(th[0]), c1=float(th[1]), k_e=float(e@th), sd_e=float(math.sqrt(e@cov@e)), cov=cov,
               chi2=float(r@r)/noise**2/max(m-q-2, 1))
    return out


def spline_fit(log, n, ns, eb, afc, ndel, noise=.04, nb=3, nw=16, degree=3, enc=True, prior=None):
    """Multi-block variable projection with a time-varying stiffness.

    Horizon = nb consecutive nw-tick blocks ending at n. Each block keeps its own degree-p drift (eliminated by
    projection, so motion invariance holds per block as in the current estimator); the stiffness is one continuous
    linear spline over the horizon with knots at the block boundaries (nb+1 knots, last knot = current stiffness).
    prior: optional (mean, information) on the knot values (arrival cost). nb=1 with a constant spline reproduces
    window_fit. Returns dict(k_e, sd_e, chi2, knots, cov) or None.
    """
    import common
    Xs, Ys, Zs, Ts, keepall = [], [], [], [], []
    for b in range(nb):
        nend = n-b*nw
        j0 = max(ndel+1, nend-1-nw); m = nend-2-j0+1
        if m < degree+3:
            return None
        X, Y, Z, keep = common._build(log, j0, m, degree, ns, eb, afc, ndel, log[nend-1, 0], log[nend-1, 2],
                                      nw*ns*1e-3, enc)
        if keep.sum() < degree+4:
            return None
        X, Y, Z = X[keep], Y[keep], Z[keep]
        tau = X[:, 1]/X[:, 0]-b                       # global horizon time in block units, 0 = now
        Xs.append(X); Ys.append(Y); Zs.append(Z); Ts.append(tau)
    m_tot = sum(len(z) for z in Zs); q = degree+1
    B = np.zeros((m_tot, nb*q)); Y = np.concatenate(Ys); Z = np.concatenate(Zs); tau = np.concatenate(Ts)
    r0 = 0
    for b in range(nb):
        B[r0:r0+len(Zs[b]), b*q:(b+1)*q] = Xs[b]; r0 += len(Zs[b])
    knots = -np.arange(nb, -1, -1, dtype=float)       # -nb ... 0
    Phi = np.clip(1-np.abs(tau[:, None]-knots[None, :]), 0, None)
    Qb, _ = np.linalg.qr(B)
    P = lambda V: V-Qb@(Qb.T@V)
    A = P(Phi*Y[:, None]); zr = P(Z)
    G = A.T@A/noise**2; g = A.T@zr/noise**2
    if prior is not None:
        G = G+prior[1]; g = g+prior[1]@prior[0]
    if np.linalg.cond(G) > 1e12:
        return None
    c = np.linalg.solve(G, g); cov = np.linalg.inv(G)
    r = zr-A@c; dof = max(m_tot-nb*q-(nb+1), 1)
    return dict(k_e=float(c[-1]), sd_e=float(np.sqrt(cov[-1, -1])), chi2=float(r@r)/noise**2/dof, knots=c, cov=cov)


class GainFilter:
    """Kalman filter on ln k fed by every projection window (methods 21.6).

    Measurement: ln k_hat of the current 16-tick window with variance (b)^2 N_eff, where b is the window's own risk
    score (relative sd inflated by the residual energy beyond the noise floor, property 1(d)), so a biased-looking
    window is down-weighted instead of rejected; overlapping windows count 1/N each (N = 16). State: level (and
    optionally rate). Process model: optional Hertz-family force term d ln k = gamma_F d ln F (gamma_F estimated as a
    state). With no process noise, no rate and no force term this is an average over accepted windows.
    """
    def __init__(self, q_level=1e-3, q_rate=0., gamma=False, q_gamma=1e-6, nw=16, dof=11, dt=.008):
        self.ql, self.qr, self.use_g, self.qg, self.N, self.dof, self.dt = q_level, q_rate, gamma, q_gamma, nw, dof, dt
        n = 1+(q_rate > 0)+bool(gamma); self.n = n; self.x = None; self.P = None; self.lF = None

    def step(self, k_hat, rel, chi2, F):
        lF = math.log(max(F, .3))
        if self.x is None:
            if not (np.isfinite(k_hat) and k_hat > 0): return None
            self.x = np.zeros(self.n); self.x[0] = math.log(k_hat); self.P = np.eye(self.n)*.25
            if self.use_g: self.x[-1] = 0.; self.P[-1, -1] = .25
            self.lF = lF; return self.x[0]
        # predict
        A = np.eye(self.n); u = np.zeros(self.n); Q = np.zeros((self.n, self.n)); Q[0, 0] = self.ql
        if self.qr > 0: A[0, 1] = self.dt; Q[1, 1] = self.qr
        if self.use_g: A[0, -1] = lF-self.lF; Q[-1, -1] = self.qg
        self.x = A@self.x; self.P = A@self.P@A.T+Q; self.lF = lF
        if np.isfinite(k_hat) and k_hat > 0 and np.isfinite(rel):
            b = rel*(1.+math.sqrt(self.dof*max(0., chi2-1.)))
            Rm = (b*b)*self.N
            H = np.zeros(self.n); H[0] = 1.
            S = H@self.P@H+Rm; K = self.P@H/S
            self.x = self.x+K*(math.log(k_hat)-self.x[0]); self.P = self.P-np.outer(K, H@self.P)
        return self.x[0]

    def sd(self):
        return math.sqrt(self.P[0, 0]) if self.P is not None else np.inf


class GainIMM:
    """Two-model IMM on x = [ln k, d ln k/dt] fed by every projection window (same measurement as GainFilter).
    Model 1 (still): level random walk, rate held at 0. Model 2 (changing): level + rate, rate random walk.
    Mode probabilities from the windows' own likelihoods; output = mixture level and its variance."""
    def __init__(self, q1=1e-3, q2l=1e-4, q2r=1., p_stay=.98, nw=16, dof=11, dt=.008):
        self.A = [np.array([[1., 0.], [0., 0.]]), np.array([[1., dt], [0., 1.]])]
        self.Q = [np.diag([q1, 1e-8]), np.diag([q2l, q2r])]
        self.T = np.array([[p_stay, 1-p_stay], [1-p_stay, p_stay]]); self.N, self.dof = nw, dof
        self.x = None

    def step(self, k_hat, rel, chi2, F=None):
        if self.x is None:
            if not (np.isfinite(k_hat) and k_hat > 0): return None
            self.x = [np.array([math.log(k_hat), 0.]) for _ in range(2)]; self.P = [np.diag([.25, 1.]) for _ in range(2)]
            self.mu = np.array([.5, .5]); return math.log(k_hat)
        c = self.mu@self.T; mix = (self.mu[:, None]*self.T)/np.maximum(c, 1e-300)
        xs, Ps = [], []
        for j in range(2):
            x0 = sum(mix[i, j]*self.x[i] for i in range(2))
            P0 = sum(mix[i, j]*(self.P[i]+np.outer(self.x[i]-x0, self.x[i]-x0)) for i in range(2))
            xs.append(self.A[j]@x0); Ps.append(self.A[j]@P0@self.A[j].T+self.Q[j])
        lik = np.ones(2)
        if np.isfinite(k_hat) and k_hat > 0 and np.isfinite(rel):
            b = rel*(1.+math.sqrt(self.dof*max(0., chi2-1.))); Rm = b*b*self.N; z = math.log(k_hat)
            for j in range(2):
                S = Ps[j][0, 0]+Rm; K = Ps[j][:, 0]/S; inn = z-xs[j][0]
                xs[j] = xs[j]+K*inn; Ps[j] = Ps[j]-np.outer(K, Ps[j][0, :])
                # likelihood of the overlapping window counts 1/N (same weighting as the update)
                lik[j] = math.exp(-.5*(inn*inn/S)/self.N)/math.sqrt(S)**(1./self.N)
        self.mu = c*lik; self.mu /= self.mu.sum(); self.mu = np.clip(self.mu, 1e-4, 1); self.mu /= self.mu.sum()
        self.x, self.P = xs, Ps
        return float(self.mu@np.array([x[0] for x in xs]))

    def sd(self):
        """sd of the mixture level (model-conditioned variances plus spread of the model means)"""
        if self.x is None: return np.inf
        m = self.mu@np.array([x[0] for x in self.x])
        return math.sqrt(sum(self.mu[j]*(self.P[j][0, 0]+(self.x[j][0]-m)**2) for j in range(2)))
