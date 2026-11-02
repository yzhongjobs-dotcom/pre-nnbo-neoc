"""Path-indexed modulus map by variable projection with an arrival cost (methods 21.6, item 2). PROTOTYPE.

Status (2026-11-02): reconstructed from the session record because it was deleted before its first commit. Not part
of any method. It was dropped after a short replay (2 scenes per family, p3_replay.py) and three quick fixes, so the
negative result limits this implementation, not finite-footprint contact modelling (handoff/03, handoff/04). It uses
footprint_kernels.py, whose geometry equals the plant's; E3 (step2_estimators.py) had already used mismatched geometry.

Window model. The current estimator regresses the filter-inverse force increment z on the probe displacement Y with
a smooth drift: z = k Y + B beta (variable projection, B a cubic time basis). If the local stiffness varies along the
scan inside the window, F_j = k_j (y_j - s_j) + g_j gives z = diag(k) Y + diag(k) B beta_1 + B beta_2: the drift
subspace becomes span[B, diag(kappa) B] with the within-window stiffness profile kappa_j = k_j / k_end. The profile
is predicted by the finite-footprint contact model from the modulus field E on a 1 mm grid along the known path:
k_j = K(x_j, d_j) E (footprint_kernels.py). The window residual
    r(theta) = (I - P_[B, kappa(theta) B]) (z - diag(k(theta)) Y),   theta = log E on the active cells,
is the separable nonlinear least-squares problem of Golub and Pereyra in its full form (the projection depends on
the nonlinear parameter); with one cell or kappa = 1 it is exactly the current projection estimator.
Across windows the map is carried in information form (arrival cost). Windows are evaluated every UPD ticks and
overlap, so each window's information is weighted UPD/N (as the existing measurement-gain filter does). Depth d_j
per row comes from the measured force and the current map (static elastic model, no relaxation factor). Up to
`iters` Gauss-Newton steps per update (Kaufman's approximation of the projection derivative), edge-preserving Huber
prior on neighbouring log-modulus differences (IRLS).
"""
import math
import numpy as np
import footprint_kernels as FK


class MapEstimator:
    def __init__(self, x0, x1, nw=16, degree=3, noise=.04, upd=4, h=1e-3, sd_inc=.5, sd_abs=3., forget=1.,
                 band=(-.045, .045), iters=3, kdeg=1, tol=.03, huber=.1):
        self.g = FK.Grid(x0, x1, h); self.nw, self.degree, self.noise, self.upd = nw, degree, noise, upd
        self.sd_inc, self.sd_abs, self.forget, self.band = sd_inc, sd_abs, forget, band
        self.iters, self.kdeg, self.tol, self.huber = iters, kdeg, tol, huber
        self.mu = None; self.Lam = None; self.dcache = {}; self.last = None

    # ---------------------------------------------------------------- model pieces
    def _init(self, k0, x, d):
        n = self.g.n
        K = FK.kernels(self.g, x, d)[2]
        self.mu = np.full(n, math.log(max(k0, 1.)/max(K.sum(), 1e-9)))
        D = np.zeros((n-1, n)); D[np.arange(n-1), np.arange(n-1)] = -1; D[np.arange(n-1), np.arange(1, n)] = 1
        self.D = D
        self.Lam = np.eye(n)/self.sd_abs**2              # data information (arrival cost) + weak absolute prior

    def _prior(self, th):
        """edge-preserving (Huber) smoothness prior on neighbouring log-modulus differences, IRLS weights at th"""
        dif = self.D@th
        wts = np.minimum(1., self.huber/np.maximum(np.abs(dif), 1e-12))
        H = self.D.T@(wts[:, None]*self.D)/self.sd_inc**2
        return H, float(np.sum(np.where(np.abs(dif) <= self.huber, dif**2, 2*self.huber*np.abs(dif)-self.huber**2)))/self.sd_inc**2

    def depth(self, x, F, E):
        """static elastic depth from the measured force and the map (a few secant steps on log d)"""
        d = 1e-3
        for _ in range(8):
            fe = FK.fe_m1_rows(self.g, x, d)[0]@E
            if fe <= 0: d *= 2; continue
            d2 = d*1.01; fe2 = FK.fe_m1_rows(self.g, x, d2)[0]@E
            s = math.log(fe2/fe)/math.log(1.01)
            d = float(np.clip(d*math.exp(math.log(max(F, .05)/fe)/max(s, .3)), 1e-6, 2e-2))
        return d

    def k_at(self, x, d, theta=None):
        E = np.exp(self.mu if theta is None else theta)
        return float(FK.kernels(self.g, x, d)[2]@E)

    def sd_lnk_at(self, x, d):
        K = FK.kernels(self.g, x, d)[2]; E = np.exp(self.mu); a = K*E/(K@E)
        C = np.linalg.inv(self.Lam+self._prior(self.mu)[0])   # marginal covariance of the map
        return float(np.sqrt(a@C@a))

    # ---------------------------------------------------------------- update
    def _basis(self, X, kap):
        """orthonormal basis of span[B, (kappa-1) B_01]: the stiffness-weighted drift columns are added only where the
        profile varies (rank-revealing; a flat profile gives exactly the current projection)"""
        Q1, _ = np.linalg.qr(X)
        C = (kap-1.)[:, None]*X[:, :self.kdeg+1]
        C = C-Q1@(Q1.T@C)
        U, sv, _ = np.linalg.svd(C, full_matrices=False)
        keep = sv > self.tol*np.linalg.norm(X[:, 0])
        return np.hstack([Q1, U[:, keep]]) if keep.any() else Q1

    def _resid(self, th, X, Y, Z, Krows):
        k = Krows@np.exp(th); Q = self._basis(X, k/k[-1])
        v = Z-k*Y; r = v-Q@(Q.T@v)
        return r, k, Q

    def update(self, log, n, X, Y, Z, rows, k_init=None):
        """X, Y, Z: window regression (rows already filtered); rows: log indices j of each row"""
        xs = .5*(log[rows, 1]+log[rows+1, 1]); Fs = log[rows+1, 5]
        if self.mu is None:
            if k_init is None or not np.isfinite(k_init): return None
            self._init(k_init, float(xs[-1]), max(float(Fs[-1])/k_init, 1e-5))
        w = self.upd/self.nw/self.noise**2
        th = self.mu.copy(); prior = self.mu.copy()
        for it in range(self.iters):
            E = np.exp(th)
            ds = np.array([self.depth(float(xs[i]), float(Fs[i]), E) for i in range(len(rows))])
            Krows = np.array([FK.kernels(self.g, float(xs[i]), float(ds[i]))[2] for i in range(len(rows))])
            act = np.where(Krows.sum(0) > 0)[0]
            if len(act) == 0: return None
            r, k, Q = self._resid(th, X, Y, Z, Krows)
            # beta of the stiffness-weighted drift is folded into the projection; derivative of the profile through
            # the projection is neglected (Kaufman); d r / d theta_c = -(I-P) diag(Y) K_c E_c
            G = Y[:, None]*(Krows[:, act]*E[act])
            J = -(G-Q@(Q.T@G))
            Lsub = self.Lam[np.ix_(act, act)]
            Hp, cp = self._prior(th); Hp = Hp[np.ix_(act, act)]
            H = Lsub+Hp+w*J.T@J
            gvec = w*J.T@r+Lsub@(th[act]-prior[act])+(self.D.T@(np.minimum(1., self.huber/np.maximum(np.abs(self.D@th), 1e-12))*(self.D@th)))[act]/self.sd_inc**2
            step = np.clip(np.linalg.solve(H, -gvec), -1., 1.)
            c0 = w*r@r+(th[act]-prior[act])@Lsub@(th[act]-prior[act])+cp
            for ls in (1., .5, .25):
                t2 = th.copy(); t2[act] += ls*step
                r2, _, _ = self._resid(t2, X, Y, Z, Krows)
                if w*r2@r2+(t2[act]-prior[act])@Lsub@(t2[act]-prior[act])+self._prior(t2)[1] < c0:
                    th = t2; break
            else:
                break
        self.mu = th
        self.Lam *= self.forget
        self.Lam[np.ix_(act, act)] += w*J.T@J
        r, k, Q = self._resid(th, X, Y, Z, Krows)
        dof = max(len(rows)-Q.shape[1]-1, 1)
        self.last = dict(chi2=float(r@r)/self.noise**2/dof, x=float(xs[-1]), d=float(ds[-1]))
        return self.last
