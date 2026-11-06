"""D⁴ with one structural break: the window may contain one abrupt change, at an unknown time, in the force not explained by the probe and in the stiffness.

Window model (same regression domain as the D⁴ host): z = k Y + Δk (Y∘H_τ) + a H_τ + b R_τ + d + ε, with d under the D⁴ prior (metric M).
H_τ is a step at τ, R_τ a ramp from τ. For each τ the generalised least squares solution is closed form; the break is kept only if it
lowers the GLS residual by more than the BIC price of its parameters (classical change-point regression, Hinkley 1970, Bai and Perron 1998,
inside the partially linear model of the D⁴ host). Without a break this is the D⁴ estimator. The stiffness reported is the one after the break.
"""
import math

import numpy as np

from contact_candidates import RIDGE, ContactEstimator, CandidatePolicy, Preprocess
import contact_candidates as C
import common
from estimators import _clip_k, _fit_record
from ridge_estimator import MIN_TICKS, RidgeProjection

MIN_AFTER_K = 6          # samples after the break needed to identify a stiffness change: about one probe cycle at 125 Hz


class BreakRidge(RidgeProjection):
    def __init__(self, penalty=1., **kw):
        super().__init__(**kw)
        self.penalty = float(penalty)

    def estimate(self, o, n, policy):
        p = policy.params; log = self._log(o, n); nw = self.nw
        self.fit = _fit_record(degree=0, tw=nw*policy.Tc)
        j0 = max(policy.ndel+1, n-1-nw); m = n-2-j0+1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        _, Y, z, keep = common._build(log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel, te, ye, nw*policy.ns*1e-3, p['enc_y'])
        broken = np.flatnonzero(~keep)
        if len(broken):
            Y, z = Y[broken[-1]+1:], z[broken[-1]+1:]
        m = len(Y)
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        M, dof = self._matrix(m)
        s2 = p['noise']**2
        # no break: the D⁴ estimator
        my = M@Y; info = float(Y@my)
        if info <= 0:
            return np.nan, np.inf, np.inf
        k0 = float(my@z)/info; r0 = z-k0*Y; q0 = float(r0@(M@r0))/s2
        best = (q0, k0, s2/info, dof, None)
        idx = np.arange(m)
        for tau, stiff in [(tau, stiff) for tau in range(2, m-1) for stiff in ((False, True) if m-tau >= MIN_AFTER_K else (False,))]:
            # a break in d alone, or in d and the stiffness: the stiffness change must earn its own price
            H = (idx >= tau).astype(float); R = np.maximum(idx-tau, 0).astype(float)
            cols = [Y, H, R] + ([Y*H] if stiff else [])
            G = np.array(cols).T
            A = G.T@M@G
            try:
                Ainv = np.linalg.inv(A)
            except np.linalg.LinAlgError:
                continue
            th = Ainv@(G.T@(M@z)); r = z-G@th; q = float(r@(M@r))/s2
            extra = G.shape[1]-1
            # BIC price: each extra parameter plus the choice of the break time
            if q+self.penalty*(extra+1)*math.log(m) < best[0]:
                k_after = th[0]+(th[3] if len(th) > 3 else 0.)
                g = np.zeros(len(th)); g[0] = 1.
                if len(th) > 3:
                    g[3] = 1.
                var = float(g@Ainv@g)*s2
                best = (q+self.penalty*(extra+1)*math.log(m), k_after, var, max(dof-extra, 1), (tau, q))
        crit, k, var, dof_b, brk = best
        q = brk[1] if brk else q0
        sigma = self.sigma_scale*math.sqrt(max(var, 1e-30))
        chi = self.chi_scale*q/dof_b
        self.dof = dof_b; self.last_break = None if brk is None else brk[0]
        k = _clip_k(k)
        if np.isfinite(k) and sigma <= .25*k:
            self.k_ref = k
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); self.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi


class ContactBreak(Preprocess, BreakRidge):
    def __init__(self, kind, penalty=1., **ridge):
        BreakRidge.__init__(self, memory=False, penalty=penalty, **dict(RIDGE, **ridge))
        self.kind = kind

    def reset(self, policy):
        BreakRidge.reset(self, policy)
        self._setup(self.kind, policy)


def build(sc, penalty=1., overrides=None, ridge=None):
    """The selected Method (H on the D⁴ host) with the break-aware window, through method3's interfaces."""
    import method3 as M
    import sweeps as S
    p = dict(S.FINAL_C); p.update(C.CTL); p.update(dict(M.TRACKER, **(overrides or {})))
    pol = CandidatePolicy('H', ContactBreak('H', penalty, **(ridge or {})), single_model=False, **p)
    return M.attach(pol, sc)
