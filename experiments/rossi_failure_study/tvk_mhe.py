"""Time-varying stiffness, joint banded IRLS with a rank-two correction.

Objective (standard Huber h):
 ||z-Ax-d||² + 2 rho s² sum h_c(Dp d/s)
                 + 2 lam_k noise² sum h_1(Dq x/delta_k),
 x=k/K, A=K[diag(Y-SY)+(SY)1'/m], s=noise/sqrt(rho).
 K is the cold quadratic constant-k estimate (a measured relative scale),
 fixed throughout the six IRLS steps. delta_k is a relative per-tick change
 (q=1) or change of per-tick slope (q=2). There is no positivity constraint;
 invalid endpoint estimates use the existing host's clipping/rejection rule.

IMPORTANT: q=2, lam_k -> infinity leaves an affine k, NOT a constant k.
 Moreover baseline chi divides a profiled GLS objective by m-p-1; the
 requested effective-df statistic divides by m-tr(H). These requirements
 cannot simultaneously reproduce the baseline. No hidden slope penalty or
 special infinity dispatch is used to obscure this discrepancy.
"""
import math
import numpy as np
from scipy.linalg import solveh_banded

import common
import contact_candidates as C
from d_robust import RobustRidge
from estimators import _clip_k, _fit_record
from ridge_estimator import MIN_TICKS


class TVKRidge(RobustRidge):
    def __init__(self, q=2, lam_k=.1, delta_k=.002, **kw):
        super().__init__(**kw)
        if q not in (1, 2) or not math.isfinite(lam_k) or not math.isfinite(delta_k) or lam_k <= 0 or delta_k <= 0:
            raise ValueError('q=1/2, positive finite lam_k and delta_k required')
        self.q = int(q); self.lam_k = float(lam_k); self.delta_k = float(delta_k)

    def reset(self, policy):
        super().reset(policy)
        self._joint_cache = {}
        self._probe_cache = {}

    def _joint_matrices(self, m):
        if m not in self._joint_cache:
            bw = 2*max(self.order, self.q)
            bases = []
            for order, parity in ((self.order, 1), (self.q, 0)):
                D = np.diff(np.eye(m), order, axis=0)
                basis = np.zeros((bw+1, 2*m, len(D)))
                for j in range(order+1):
                    basis[2*j, parity:2*(m-j):2] = (D[:, j:]*D[:, :m-j]).T
                bases.append((D, basis.reshape(-1, len(D))))
            self._joint_cache[m] = bw, bases
        return self._joint_cache[m]

    def solve_window(self, Y, z, noise, dense_audit=False):
        """Pure window solver, also exposed for reproducible synthetic checks.

        The full Hessian is T + U C U', where T is banded, U=[B'u,h],
        C=[[0,1],[1,u'u]]. h selects mean(x), B has rows [v_i,1].
        The inverse is T^-1 - T^-1 U (C^-1+U'T^-1 U)^-1 U'T^-1.
        Covariance uses the final solve's frozen weights, including the slope.
        chi uses z'(I-H)z (the same profiled numerator as RobustRidge),
        divided by noise²*(m-tr H). Both dof conventions are exposed.
        """
        Y = np.asarray(Y, float); z = np.asarray(z, float)
        m = len(Y); Y = Y-Y.mean()  # constant absorbed by unpenalized d nullspace
        D, dbasis = self._difference_bands(m)
        ab0 = (dbasis@np.ones(m-self.order)).reshape(self.order+1, m)
        ab0[0] += 1.
        sy, sz = solveh_banded(ab0, np.column_stack((Y,z)), lower=True,
                              check_finite=False).T
        my = Y-sy; info = float(Y@my)
        if info <= 1e-24:
            raise np.linalg.LinAlgError('no probe information')
        K = float(np.clip((my@z)/info, 200., 100000.))
        # The probe regressor ALWAYS uses the quadratic D4 smoother, even
        # when the disturbance prior is the p=3 comparison variant.
        if self.order == 4:
            probe_sy = sy
        else:
            if m not in self._probe_cache:
                dp = np.diff(np.eye(m),4,axis=0)
                rp = (2*math.sin(math.pi*self.fc*self.Tc))**-8
                full = rp*dp.T@dp+np.eye(m)
                band = np.zeros((5,m))
                for j in range(5):band[j,:m-j] = np.diag(full,-j)
                self._probe_cache[m] = band
            probe_sy = solveh_banded(self._probe_cache[m],Y,lower=True,check_finite=False)
        v = K*(Y-probe_sy); u = K*probe_sy
        bw, bases = self._joint_matrices(m)
        Dd, bd = bases[0]; Dk, bk = bases[1]
        h = np.zeros(2*m); h[::2] = 1./m
        bu = np.empty(2*m); bu[::2] = v*u; bu[1::2] = u
        U = np.column_stack((bu,h))
        cinv = np.array([[-float(u@u),1.],[1.,0.]])
        rhs = np.empty(2*m); rhs[::2] = v*z; rhs[1::2] = z
        rhs += h*float(u@z)
        obs = np.zeros((m,2*m)); ii = np.arange(m)
        obs[ii,2*ii] = v; obs[ii,2*ii+1] = 1.
        obs += np.outer(u,h)
        wd = np.ones(m-self.order); wk = np.ones(m-self.q)
        rk = self.lam_k*noise**2/self.delta_k**2
        for it in range(self.iters):
            ab = (self.rho*(bd@wd)+rk*(bk@wk)).reshape(bw+1,2*m)
            ab[0,::2] += v*v; ab[0,1::2] += 1.
            ab[1,::2] += v
            # Solve all final covariance/hat RHS with the very same factorization.
            rr = np.column_stack((rhs,U,obs.T)) if it==self.iters-1 else np.column_stack((rhs,U))
            solved = solveh_banded(ab, rr, lower=True, check_finite=False)
            tu = solved[:,1:3]
            small = cinv+U.T@tu
            invrr = solved-tu@np.linalg.solve(small, U.T@solved)
            theta = invrr[:,0]
            if it != self.iters-1:
                wd = np.minimum(1., self.c/np.maximum(np.abs(Dd@theta[1::2])/(noise/math.sqrt(self.rho)),1e-12))
                wk = np.minimum(1., 1./np.maximum(np.abs(Dk@theta[::2])/self.delta_k,1e-12))
        # Last k variance: solve one column, preserving full nuisance uncertainty.
        e = np.zeros(2*m); e[-2] = 1.
        te = solveh_banded(ab,e,lower=True,check_finite=False)
        ie = te-tu@np.linalg.solve(small,U.T@te)
        var = float(ie[-2])*K*K
        trace = float(np.sum(obs.T*invrr[:,3:]))
        if not np.isfinite(var) or var <= 0 or not np.isfinite(trace):
            raise np.linalg.LinAlgError('invalid endpoint covariance')
        dof = max(m-trace,1e-9)
        resid = z-obs@theta
        # Stable nonnegative profiled objective, avoiding z'z-rhs'theta cancellation.
        numerator = float(resid@resid+self.rho*np.sum(wd*(Dd@theta[1::2])**2)
                          +rk*np.sum(wk*(Dk@theta[::2])**2))
        result = dict(k=K*theta[::2], d=theta[1::2], sigma=self.sigma_scale*noise*math.sqrt(max(var,0.)),
                      chi=self.chi_scale*numerator/noise**2/dof,
                      legacy_chi=self.chi_scale*numerator/noise**2/max(m-self.order-1,1),
                      dof=dof, trace=trace, K=K, weights_d=wd, weights_k=wk)
        if dense_audit:
            full = obs.T@obs
            full[1::2,1::2] += self.rho*(Dd.T*wd)@Dd
            full[::2,::2] += rk*(Dk.T*wk)@Dk
            inv = np.linalg.inv(full); ref = inv@rhs
            result['dense_theta_error'] = float(np.max(np.abs(ref-theta)))
            result['dense_var_relative_error'] = float(abs(inv[-2,-2]-ie[-2])/abs(inv[-2,-2]))
            result['dense_trace_error'] = float(abs(np.trace(obs@inv@obs.T)-trace))
            result['normal_residual'] = float(np.linalg.norm(full@theta-rhs)/np.linalg.norm(rhs))
        return result

    def estimate(self, o, n, policy):
        p = policy.params; log = self._log(o,n); nw = self.nw
        self.fit = _fit_record(degree=0,tw=nw*policy.Tc)
        j0 = max(policy.ndel+1,n-1-nw); m = n-1-j0
        if m < MIN_TICKS: return np.nan,np.inf,np.inf
        _,Y,z,keep = common._build(log,j0,m,0,policy.ns,policy.eb,policy.afc,policy.ndel,
                                  log[n-1,0],log[n-1,2],nw*policy.ns*1e-3,p['enc_y'])
        broken = np.flatnonzero(~keep)
        if len(broken): Y,z = Y[broken[-1]+1:],z[broken[-1]+1:]
        if len(Y) < MIN_TICKS: return np.nan,np.inf,np.inf
        try: r = self.solve_window(Y,z,p['noise'])
        except np.linalg.LinAlgError: return np.nan,np.inf,np.inf
        self.last = r; self.dof = r['dof']; self.weights = r['weights_d']
        k = _clip_k(float(r['k'][-1])); sigma,chi = r['sigma'],r['chi']
        if np.isfinite(k) and sigma <= .25*k: self.k_ref = k
        self.fit = _fit_record(k,sigma,chi,degree=0,tw=nw*policy.Tc)
        self.fit['te'] = float(log[n-1,0])
        return k,sigma,chi


class ContactTVK(C.Preprocess,TVKRidge):
    def __init__(self, kind='H', **kw):
        TVKRidge.__init__(self,memory=False,**dict(C.RIDGE,**kw)); self.kind=kind

    def reset(self,policy):
        TVKRidge.reset(self,policy); self._setup(self.kind,policy)


def build(sc,cfg=None):
    import method3 as M
    import sweeps
    c = dict(cfg or {}); gate = c.pop('gate',.551)
    c.setdefault('nw',24); c.setdefault('order',4); c.setdefault('fc',14.)
    c.setdefault('sigma_scale',.796); c['chi_scale']=gate
    p = dict(sweeps.FINAL_C); p.update(C.CTL); p.update(M.TRACKER)
    return M.attach(C.CandidatePolicy('H',ContactTVK(**c),single_model=False,**p),sc)
