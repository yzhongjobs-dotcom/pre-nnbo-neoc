"""Literature mechanism adaptations, force/encoder only. See baseline_report.md.

k is in N/m, d is in the one-step force-filter equation domain (N).
The phase detector is a separate D4 observer and never supplies k.
"""
import math
from collections import deque
from functools import lru_cache
from types import SimpleNamespace

import numpy as np

import sweeps  # establishes the existing benchmark import paths
import common
import contact_candidates as C
import method3 as M
from estimators import BandPassEstimator, R2Estimator, _clip_k, _fit_record
from ridge_estimator import RidgeProjection


def kf_step(x, P, h, z, R, A, Q):
    """Linear Gaussian update; Joseph covariance, no clipping in the recursion."""
    xp = A @ x; pp = A @ P @ A.T + Q
    ph = pp @ h; s = float(h @ ph + R); e = float(z - h @ xp)
    gain = ph / s; b = np.eye(len(x)) - np.outer(gain, h)
    post = b @ pp @ b.T + R * np.outer(gain, gain)
    return xp + gain * e, (post + post.T)/2, e*e/s


def process_model(qk, qd, dt, d_order=1):
    """Shared exact discrete integrated-Wiener model in scaled k coordinates."""
    A=np.eye(d_order+1);Q=np.zeros_like(A);Q[0,0]=qk*qk*dt/1e6
    if d_order==1:Q[1,1]=qd*qd*dt
    elif d_order==2:
        A[1,2]=dt
        Q[1:,1:]=qd*qd*np.array([[dt**3/3,dt**2/2],[dt**2/2,dt]])
    else:raise ValueError(d_order)
    return A,Q


def mhe_solve(Y, z, R, mean=None, cov=None, qk=math.inf, qd=math.inf, dt=.008, d_order=1):
    """Exact unconstrained linear MHE; q=0 means a constant over the window.

    State is [k/1000,d]. Transition variances are qk^2*dt/1e6, qd^2*dt.
    Infinite transition precision is handled by exact variable elimination,
    not by ill-conditioned huge weights. Returns trajectory, endpoint variance,
    and weighted *measurement* residual sum (arrival/process costs excluded).
    """
    if d_order==2:
        return mhe_second_order(Y,z,R,mean,cov,qk,qd,dt)
    Y = np.asarray(Y)*1000.; z = np.asarray(z); n = len(z)
    if qk == 0 and qd == 0:
        H = np.c_[Y, np.ones(n)]
        J = H.T @ H/R; b = H.T @ z/R
        if mean is not None:
            prec = np.linalg.inv(cov); J += prec; b += prec @ mean
        P = np.linalg.inv(J); x = P @ b
        return np.tile(x, (n, 1)), float(P[0, 0]*1e6), float(np.sum((z-H@x)**2)/R)
    if min(qk, qd) <= 0 or not np.isfinite(qk+qd):
        raise ValueError('use both finite positive q or both zero')
    # Block tridiagonal Hessian, scalar half bandwidth 2; two RHS include
    # the selected inverse column for the endpoint k variance.
    J = np.zeros((3, 2*n)); rhs = np.zeros((2*n, 2))
    J[0, ::2] = Y*Y/R; J[0, 1::2] = 1/R
    J[1, :-1:2] = Y/R
    rhs[::2, 0] = Y*z/R; rhs[1::2, 0] = z/R
    for d, v in enumerate((qk*qk*dt/1e6, qd*qd*dt)):
        ix = np.arange(d, 2*n, 2); w = 1/v
        J[0, ix[:-1]] += w; J[0, ix[1:]] += w
        J[2, ix[:-1]] = -w
    if mean is not None:
        prec = np.linalg.inv(cov)
        J[0, :2] += np.diag(prec); J[1, 0] += prec[1, 0]
        rhs[:2, 0] += prec @ mean
    rhs[-2, 1] = 1.
    from scipy.linalg import cholesky_banded, cho_solve_banded
    sol = cho_solve_banded((cholesky_banded(J, lower=True, check_finite=False), True), rhs, check_finite=False)
    tr = sol[:, 0].reshape(n, 2)
    return tr, float(sol[-2, 1]*1e6), float(np.sum((z-Y*tr[:, 0]-tr[:, 1])**2)/R)


@lru_cache(maxsize=256)
def _mhe_process_precision(n,qk,qd,dt):
    A,Q=process_model(qk,qd,dt,2);W=np.linalg.inv(Q)
    J=np.zeros((6,3*n));left=A.T@W@A;off=-W@A
    for j in range(n-1):
        for a in range(3):
            for b in range(a+1):
                J[a-b,3*j+b]+=left[a,b]
                J[a-b,3*(j+1)+b]+=W[a,b]
            for b in range(3):J[3+a-b,3*j+b]+=off[a,b]
    return J


def mhe_second_order(Y,z,R,mean,cov,qk,qd,dt):
    """Block-banded exact quadratic MHE, state [k/1000,d,d_dot].

    The left arrival is already predictive. Process cost starts at sample 1.
    Uses exact non-diagonal Q; no artificial diagonal regularization.
    """
    from scipy.linalg import cholesky_banded,cho_solve_banded
    n=len(z);Y=np.asarray(Y)*1000.;z=np.asarray(z)
    if min(qk,qd)<=0:raise ValueError('second order requires finite positive q')
    J=_mhe_process_precision(n,qk,qd,dt).copy();rhs=np.zeros((3*n,2))
    # Lower-triangle band storage: J[i-j,j] = full Hessian[i,j].
    def block(b,offset):
        for i in range(3):
            for j in range(i+1):J[i-j,offset+j]+=b[i,j]
    if mean is not None:
        precision=np.linalg.inv(cov);block(precision,0);rhs[:3,0]=precision@mean
    J[0,::3]+=Y*Y/R;J[0,1::3]+=1/R;J[1,::3]+=Y/R
    rhs[::3,0]+=Y*z/R;rhs[1::3,0]+=z/R;rhs[-3,1]=1.
    sol=cho_solve_banded((cholesky_banded(J,lower=True,check_finite=False),True),rhs,check_finite=False)
    tr=sol[:,0].reshape(n,3)
    return tr,float(sol[-3,1]*1e6),float(np.sum((z-Y*tr[:,0]-tr[:,1])**2)/R)


class EquationEstimator:
    """Read each causal interval once, at a fixed measured reference position."""
    def reset(self, policy):
        self.yref = None; self.last_j = -1; self.fit = _fit_record(degree=0)
        self.dof = 22; self.nis = deque(maxlen=24); self.age = 0
        self._restart()

    def _restart(self):
        self.x = np.array([30., 0.] + ([0.] if getattr(self, 'd_order', 1) == 2 else []))
        self.P = np.diag([300.**2, 100.**2] + ([100.**2] if len(self.x) == 3 else []))
        self.nis.clear(); self.age = 0

    def pair(self, o, n, policy):
        if n < policy.ndel+4 or n-2 <= self.last_j:
            return None
        j = n-2; self.last_j = j
        if self.yref is None:
            self.yref = float(o.log[0, 2])
        _, y, z, keep = common._build(o.log, j, 1, 0, policy.ns, policy.eb, policy.afc,
                                      policy.ndel, o.log[n-1, 0], self.yref, policy.Tc,
                                      policy.params['enc_y'])
        if not keep[0]:
            self._restart(); return None
        return float(y[0]), float(z[0])

    def record(self, k, sig, chi, policy, n):
        self.fit = _fit_record(k, sig, chi, degree=0, tw=max(self.age, 1)*policy.Tc)
        return k, sig, chi


class JointKalman(EquationEstimator):
    """Haninger-style augmented random-walk states, adapted to measured force."""
    def __init__(self, q_k=300., q_d=.3, d_order=2, sigma_scale=.796, chi_scale=.551):
        self.q_k=q_k; self.q_d=q_d; self.d_order=d_order
        self.sigma_scale=sigma_scale; self.chi_scale=chi_scale

    def estimate(self, o, n, policy):
        pair = self.pair(o, n, policy)
        if pair is None:
            return self.record(np.nan, np.inf, np.inf, policy, n)
        y, z = pair; dt=policy.Tc; size=len(self.x)
        A,Q=process_model(self.q_k,self.q_d,dt,self.d_order)
        h=np.array([1000*y, 1.] + ([0.] if size == 3 else []))
        self.x, self.P, nis=kf_step(self.x,self.P,h,z,policy.params['noise']**2,A,Q)
        self.age+=1; self.nis.append(nis)
        return self.record(_clip_k(self.x[0]*1000), self.sigma_scale*math.sqrt(max(self.P[0,0],0))*1000,
                           self.chi_scale*float(np.mean(self.nis)), policy, n)


class MovingHorizon(EquationEstimator):
    """Zhang-style arrival + inverse-noise history; k,d may vary in the window.

    Arrival is the predictive marginal using only samples BEFORE the left edge.
    Reusing the last endpoint posterior at the left edge would double count data.
    No true surface, robot torques or dynamics are supplied.
    """
    def __init__(self, N=24, q_k=3000., q_d=.3, sigma_scale=.796, chi_scale=.551, d_order=1, chi_mode='smooth'):
        self.N=N; self.q_k=q_k; self.q_d=q_d; self.d_order=d_order; self.chi_mode=chi_mode
        self.sigma_scale=sigma_scale; self.chi_scale=chi_scale

    def _restart(self):
        super()._restart(); self.history=deque(maxlen=self.N)

    def estimate(self,o,n,policy):
        pair=self.pair(o,n,policy)
        if pair is None:
            return self.record(np.nan,np.inf,np.inf,policy,n)
        y,z=pair; dt=policy.Tc; R=policy.params['noise']**2
        A,Q=process_model(self.q_k,self.q_d,dt,self.d_order)
        # Prior belongs to THIS observation, before assimilating it.
        self.history.append((y,z,A@self.x,A@self.P@A.T+Q))
        self.x,self.P,nis=kf_step(self.x,self.P,np.array([1000*y,1.]+([0.] if self.d_order==2 else [])),z,R,A,Q)
        self.nis.append(nis)
        self.age+=1; self.dof=max(len(self.history)-len(self.x),1)
        if self.chi_mode=='nis':
            # Exact equivalent endpoint; match Kalman startup as well as its
            # 24-sample NIS. This diagnostic never enters development search.
            self.dof=22
            return self.record(_clip_k(self.x[0]*1000),self.sigma_scale*math.sqrt(max(self.P[0,0],0))*1000,
                               self.chi_scale*float(np.mean(self.nis)),policy,n)
        if len(self.history)<4:
            return self.record(np.nan,np.inf,np.inf,policy,n)
        Y=np.array([a[0] for a in self.history]); zz=np.array([a[1] for a in self.history])
        first=self.history[0]
        tr,var,sse=mhe_solve(Y,zz,R,first[2],first[3],self.q_k,self.q_d,dt,self.d_order)
        # Diagnostics only: expose the endpoint-equivalent Kalman NIS gate.
        chi=sse/self.dof
        return self.record(_clip_k(tr[-1,0]*1000),self.sigma_scale*math.sqrt(max(var,0)),
                           self.chi_scale*chi,policy,n)


def xue_update(theta, phi_previous, log_force, delta, gain=1.):
    """RAS2025 (15)-(16), gain=1 is exactly the printed law (no dt)."""
    e=float(log_force-theta@phi_previous)
    return theta+(gain*phi_previous*e/(1+phi_previous@phi_previous) if abs(e)>2*delta else 0.)


class XueHC:
    """Original log HC regressor, NOT the probe-frequency regression.

    Outputs tangent stiffness. Fixed sigma is a declared interface addition;
    optional log-force residual chi gating is also an addition, not Xue's claim.
    H-preprocessed displacement is used as requested, but has no known DC depth.
    """
    def __init__(self, delta=.05, gain=1., initial_k=3000., initial_n=1.26,
                 residual_gate=False, sigma_scale=.796, chi_scale=.551, contact_position=None):
        self.delta=delta;self.gain=gain;self.initial_k=initial_k;self.initial_n=initial_n
        self.residual_gate=residual_gate;self.sigma_scale=sigma_scale;self.chi_scale=chi_scale
        self.contact_position=contact_position

    def reset(self,policy):
        self.theta=np.array([math.log(self.initial_k),self.initial_n,0.])
        self.contact=self.contact_position;self.prev_phi=None;self.prev_y=None;self.last=-1
        self.res=deque(maxlen=24);self.dof=21;self.fit=_fit_record(degree=0)
        self.invalid_depth=0;self.updates=0

    def estimate(self,o,n,policy):
        row=o.log[n-1]; f=float(row[5]);y=float(row[2]); dt=policy.Tc
        v=0. if self.prev_y is None else (y-self.prev_y)/dt;self.prev_y=y
        if self.contact is None and f>.5:
            self.contact=y  # first detected contact, never from truth
        depth=y-self.contact if self.contact is not None else -1.
        if f<=.5 or depth<=1e-7:
            self.prev_phi=None;self.invalid_depth+=1
            return self._out(np.nan,np.inf,np.inf)
        lf=math.log(f)
        if self.prev_phi is not None:
            e=lf-float(self.theta@self.prev_phi)
            before=self.theta.copy();self.theta=xue_update(self.theta,self.prev_phi,lf,self.delta,self.gain)
            self.updates+=int(np.any(before!=self.theta))
            # Extra chi uses log-force measurement noise by the delta method.
            self.res.append((e*f/policy.params['noise'])**2)
        self.prev_phi=np.array([1.,math.log(depth),v])
        logkt=math.log(max(self.theta[1],1e-300))+self.theta[0]+(self.theta[1]-1)*math.log(depth)
        k=_clip_k(math.exp(min(max(logkt,-700.),700.))) if self.theta[1]>0 else np.nan
        chi=self.chi_scale*float(np.mean(self.res)) if self.residual_gate and self.res else 0.
        # A fixed positive sigma preserves the host interface; no covariance is claimed.
        return self._out(k,self.sigma_scale*1.,chi)

    def _out(self,k,sig,chi):
        self.fit=_fit_record(k,sig,chi,degree=0,tw=.192);return k,sig,chi


class ContactAdapter(C.Preprocess):
    """Common H preprocessing with an independent D4 phase-only matrix."""
    def __init__(self, inner):
        self.inner=inner;self.nw=24;self.phase=RidgeProjection(nw=24,fc=14.,order=4)

    def reset(self,policy):
        self._setup('H',policy);self.inner.reset(policy);self.phase.reset(policy)
        if isinstance(self.inner, (BandPassEstimator, R2Estimator)):
            self.inner.displacement=EncoderMatched(policy)

    def _matrix(self,m):
        return self.phase._matrix(m)

    @property
    def fit(self): return self.inner.fit

    @property
    def dof(self): return self.inner.dof

    def estimate(self,o,n,policy):
        view=SimpleNamespace(log=self._log(o,n),n=n)
        return self.inner.estimate(view,n,policy)


class EncoderMatched:
    """Same measured endpoints and intra-tick shape as common._build.

    Raw commands cannot interpolate an H-transformed displacement. Integrating
    its measured endpoints keeps Rossi's two BP channels in the same domain.
    """
    def __init__(self,p):
        self.p=p;self.last=0;self.filtered=None

    def update(self,o):
        p=self.p
        if self.filtered is None:self.filtered=float(o.log[0,2])
        for j in range(self.last,o.n-1):
            y0,y1=o.log[j:j+2,2];ebn=p.eb**p.ns
            for mm in range(1,p.ns+1):
                fr=(1-p.eb**(mm-1))/(1-ebn) if ebn<1 else (mm-1)/p.ns
                self.filtered+=p.afc*(y0+fr*(y1-y0)-self.filtered)
        self.last=o.n-1
        return self.filtered


def build(method, cfg, sc):
    cfg=dict(cfg);gate=cfg.pop('gate',.551)
    sigma_scale=cfg.pop('sigma_scale',.796)
    # All estimator comparisons share the selected host controller and tracker.
    if method in ('Huber-D4','D4','P5'):
        if method=='P5':
            C.P5=dict(dict(nw=16,degree=5,sigma_scale=sigma_scale,chi_scale=gate),**cfg)
            pol=C.build('P5:H',M.TRACKER)
            # Same host (P5's older build switches these off).
            pol.acquire=True;pol.centre=True
        else:
            C.RIDGE=dict(dict(nw=24,order=4,fc=14.,sigma_scale=sigma_scale,chi_scale=gate),**cfg)
            pol=C.build('H',M.TRACKER)
            if method=='Huber-D4':
                from d_robust import ContactRobust
                pol.estimator=ContactRobust('H',**C.RIDGE)
        return M.attach(pol,sc)
    kinds={'Rossi-ID':BandPassEstimator,'Rossi-R2':R2Estimator,'Kalman':JointKalman,
           'MHE':MovingHorizon,'Xue-const':XueHC,'Xue-chi':XueHC}
    if method.startswith('Xue'):
        cfg['residual_gate']=method=='Xue-chi'
    cls=kinds[method]
    if method=='Rossi-ID' and cfg.pop('filter','butterworth')=='r2':cls=R2Estimator
    inner=cls(sigma_scale=sigma_scale,chi_scale=gate,**cfg)
    p=dict(sweeps.FINAL_C);p.update(C.CTL);p.update(M.TRACKER)
    return M.attach(C.CandidatePolicy('H',ContactAdapter(inner),single_model=False,**p),sc)
