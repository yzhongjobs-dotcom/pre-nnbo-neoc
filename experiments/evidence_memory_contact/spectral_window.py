"""Finite-window frequency separation derived from a motion model.

    z = k_e Y + k_1 (t/T) Y + B beta + xi + eps,   Cov(eps) = R I,
    xi is triple-integrated white noise with spectral density q_F = R dt w_c^6.

xi has spectrum proportional to 1/w^6, so weights act as a third-order high-pass with corner w_c,
the only separation parameter. q_F is defined from the noise and the corner, so it needs no
stiffness estimate and the estimator stays linear. w_c -> 0 recovers projection on B (with the
slope term, the lean estimator). The fit is a precomputed quadratic form per window geometry,
as cheap as projection. The variance includes xi, so it reflects modelled interference.
"""
import math
import numpy as np

DT = .008
AF = 1-math.exp(-2*math.pi*30*.001)
W = AF*(1-AF)**(7-np.arange(8))
RHO = math.exp(-.008/.019)
FRAC = -np.expm1(-np.arange(8)*.001/.019)/(1-RHO)


def _kernel(n):
    """Covariance of triple-integrated unit white noise at the interval centres (filter gain W0)."""
    t = (np.arange(n)+.5)*DT
    s, u = np.meshgrid(t, t, indexing='ij')
    a = np.minimum(s, u); b = np.maximum(s, u)
    # Closed form of int_0^a (a-v)^2 (b-v)^2 / 4 dv.
    K = (a**5/20+(b-a)*a**4/8+(b-a)**2*a**3/12)
    return K*W.sum()**2


class SpectralWindow:
    def __init__(self, n=32, corner_hz=8., degree=3, slope=True, noise=.04):
        self.n, self.degree, self.slope, self.noise = n, degree, slope, noise
        self.corner = corner_hz
        R = noise**2
        wc = 2*math.pi*corner_hz
        q = R*DT*wc**6
        C = R*np.eye(n)+q*_kernel(n)
        m = np.arange(8)*.001
        tt = ((np.arange(n)[:, None]*DT+m)-n*DT)/(n*DT)
        self.tt = tt
        B = np.column_stack([(tt**p)@W for p in range(degree+1)])
        L = np.linalg.cholesky(C)
        Li = np.linalg.solve(L, np.eye(n))
        # Orthonormal whitened nuisance basis: profiling beta stays accurate at small corners.
        Q, _ = np.linalg.qr(Li@B)
        P = np.eye(n)-Q@Q.T
        self.M = Li.T@P@Li
        self.logdet = 2*float(np.log(np.diag(L)).sum())+2*float(np.log(np.abs(np.diag(np.linalg.qr(Li@B)[1]))).sum())
        self.B = B

    def regressors(self, y):
        """Filter-aligned displacement and time-weighted displacement, in the common domain."""
        n = self.n
        yy = y[:-1, None]+(y[1:]-y[:-1])[:, None]*FRAC-y[-1]
        return yy@W, (self.tt*yy)@W

    def fit(self, y, F):
        """y and F hold n+1 consecutive measured samples ending now."""
        z = F[1:]-(1-AF)**8*F[:-1]
        u0, u1 = self.regressors(y)
        M = self.M
        if self.slope:
            U = np.column_stack([u0, u1])
            A = U.T@M@U
            try:
                Ai = np.linalg.inv(A)
            except np.linalg.LinAlgError:
                return None
            th = Ai@(U.T@M@z)
            r = z-U@th
            quad = float(r@M@r)
            dof = self.n-self.degree-1-2
            return dict(k=float(th[0]), slope=float(th[1]), sigma=math.sqrt(max(Ai[0, 0], 0.)),
                        chi2=quad/dof, rel=None)
        a = float(u0@M@u0)
        if a <= 0:
            return None
        k = float(u0@M@z)/a
        r = z-k*u0
        return dict(k=k, slope=0., sigma=1/math.sqrt(a), chi2=float(r@M@r)/(self.n-self.degree-2), rel=None)


class AdaptiveSpectralWindow:
    """Corner chosen per window by restricted marginal likelihood.

    The evidence comes from the interference itself, whose energy is far above the sensor noise, so
    one window identifies it. Probe parameters are fixed effects, which keeps the choice from being
    driven by the gain. Nothing is carried between windows, so a wrong corner cannot persist.
    """
    def __init__(self, n=32, corners=(.5, 4., 8., 12., 16., 24.), degree=3, slope=True, noise=.04, mix=True, change_prior=None):
        self.bank = [SpectralWindow(n=n, corner_hz=c, degree=degree, slope=slope, noise=noise) for c in corners]
        self.corners = corners; self.n = n; self.slope = slope; self.mix = mix
        # Prior sd of the relative gain change across a window. None: free slope; 0: constant gain.
        self.change_prior = change_prior

    def fit(self, y, F):
        z = F[1:]-(1-AF)**8*F[:-1]
        u0, u1 = self.bank[0].regressors(y)
        U = np.column_stack([u0, u1]) if self.slope else u0[:, None]
        ll = []; sols = []
        for f in self.bank:
            A = U.T@f.M@U
            b = U.T@f.M@z
            if self.slope and self.change_prior is not None:
                k0 = max(float(b[0]/A[0, 0]), 200.)
                A = A.copy(); A[1, 1] += 1./(self.change_prior*k0)**2
            try:
                Ai = np.linalg.inv(A)
            except np.linalg.LinAlgError:
                return None
            th = Ai@b
            quad = float(z@f.M@z-b@th)
            # Profile the probe parameters: a log|A| term would reward corners that discard probe information.
            ll.append(-.5*(quad+f.logdet))
            sols.append((float(th[0]), float(Ai[0, 0]), quad))
        ll = np.asarray(ll); w = np.exp(ll-ll.max()); w /= w.sum()
        ks = np.array([s[0] for s in sols]); vs = np.array([s[1] for s in sols])
        if self.mix:
            k = float(w@ks); var = float(w@(vs+(ks-k)**2))
        else:
            i = int(np.argmax(w)); k, var = float(ks[i]), float(vs[i])
        i = int(np.argmax(w))
        dof = self.n-self.bank[0].degree-1-U.shape[1]
        return dict(k=k, sigma=math.sqrt(max(var, 0.)), chi2=sols[i][2]/dof, corner=float(w@np.asarray(self.corners)), rel=None)
