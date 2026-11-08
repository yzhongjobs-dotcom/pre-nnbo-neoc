"""D4 arrival-cost Kalman filter and Willsky--Jones sparse jump GLR.

The disturbance coordinates are BACKWARD DIFFERENCES, not continuous
derivatives: Fd=triu(ones(4,4)), Gd=ones(4), Var(w)=noise**2/rho.
This is exactly the finite-difference prior used by RidgeProjection.
Jump covariance is conditional on a selected time/type; it is NOT Bayesian
model averaging over undetected jumps. No robust measurement update is used.
"""
import math
from collections import deque
import numpy as np
from scipy.stats import chi2
from scipy.linalg import solveh_banded

import sweeps  # existing benchmark import-path initialization
import common
import contact_candidates as C
from ridge_estimator import RidgeProjection, MIN_TICKS
from estimators import _clip_k, _fit_record


def dynamics():
    f = np.eye(5); f[1:, 1:] = np.triu(np.ones((4, 4)))
    g = np.r_[0., np.ones(4)]
    return f, g


def batch_arrival(y, z, noise, rho, scale=3000.):
    """Exact diffuse initialization, flat k and cubic polynomial coefficients.

    Returns posterior at the final sample. This same routine permits an
    independent batch-vs-sequential equivalence check with q_k=0.
    """
    y = np.asarray(y); z = np.asarray(z); m = len(y)
    d = np.diff(np.eye(m), 4, axis=0)
    a = np.column_stack((scale*y, np.eye(m)))
    normal = a.T@a
    normal[1:, 1:] += rho*d.T@d
    cov = np.linalg.inv(normal)*noise**2
    theta = np.linalg.solve(normal, a.T@z)
    b = np.zeros((5, m+1)); b[0, 0] = 1.
    for q in range(4):
        b[1+q, 1:] = np.diff(np.eye(m), q, axis=0)[-1]
    return b@theta, b@cov@b.T


class JumpFilter:
    """Scalar innovations, vectorized bank of at most L candidate times.

    q_k is the per-tick standard deviation of k/scale (not variance).
    alpha bounds the per-tick union false-alarm probability under H0,
    conditional on the correct deterministic regressor/model and Gaussian noise.
    h is fixed from alpha/(number_of_families*L), never fitted to episodes.
    """
    def __init__(self, noise=.04, rho=None, q_k=.0005, L=96, alpha=1e-4,
                 impulse=True, detect=True, scale=3000., nw=24):
        if rho is None: rho = (2*math.sin(math.pi*14*.008))**-8
        if noise <= 0 or q_k < 0 or L < 2 or not 0 < alpha < 1:
            raise ValueError('invalid noise/q_k/L/alpha')
        self.noise = float(noise); self.rho = float(rho); self.q_k = float(q_k)
        self.L = int(L); self.alpha = float(alpha); self.impulse = bool(impulse)
        self.detect = bool(detect); self.scale = float(scale)
        self.F, g = dynamics(); self.Q = np.outer(g, g)*noise**2/rho
        self.Q[0, 0] = q_k**2
        tail = alpha/((3 if impulse else 2)*L)
        self.h1 = float(chi2.isf(tail, 1)); self.h2 = float(chi2.isf(tail, 2))
        self.x = None; self.P = None; self.start_y = []; self.start_z = []
        self.E = np.zeros((0, 5, 3)); self.J = np.zeros((0, 3, 3))
        self.b = np.zeros((0, 3)); self.theta = np.zeros(0, dtype=int)
        self.index = -1; self.events = []; self.nis = deque(maxlen=nw)
        self.last = {}; self.max_cov_asymmetry = 0.

    def clear_bank(self):
        self.E = np.zeros((0, 5, 3)); self.J = np.zeros((0, 3, 3))
        self.b = np.zeros((0, 3)); self.theta = np.zeros(0, dtype=int)

    def step(self, y, z):
        self.index += 1
        if self.x is None:
            self.start_y.append(y); self.start_z.append(z)
            if len(self.start_y) < MIN_TICKS: return np.nan, np.inf, np.inf
            try:
                self.x, self.P = batch_arrival(self.start_y, self.start_z,
                                               self.noise, self.rho, self.scale)
            except np.linalg.LinAlgError:
                self.start_y.pop(0); self.start_z.pop(0)
                return np.nan, np.inf, np.inf
            return self.output()
        h = np.array([self.scale*y, 1., 0., 0., 0.])
        xp = self.F@self.x; pp = self.F@self.P@self.F.T+self.Q
        ph = pp@h; v = float(h@ph+self.noise**2)
        gain = ph/v; innovation = float(z-h@xp)
        self.x = xp+gain*innovation
        # Joseph form remains positive semidefinite through long quiet runs.
        ik = np.eye(5)-np.outer(gain, h)
        self.P = ik@pp@ik.T+np.outer(gain, gain)*self.noise**2
        self.P = (self.P+self.P.T)*.5
        self.nis.append(innovation**2/v)
        self.last = dict(innovation=innovation, variance=v, detected=False)
        if self.detect:
            # Propagate true-minus-nominal state signature; append a new jump.
            ep = self.F@self.E[-(self.L-1):]
            new = np.zeros((1, 5, 3)); new[0, 0, 0] = 1.; new[0, 1, 1] = 1.
            ep = np.concatenate((ep, new))
            self.J = np.concatenate((self.J[-(self.L-1):], np.zeros((1, 3, 3))))
            self.b = np.concatenate((self.b[-(self.L-1):], np.zeros((1, 3))))
            self.theta = np.r_[self.theta[-(self.L-1):], self.index]
            signature = h@ep
            # Measurement impulse has no true persistent state jump. Its
            # subsequent negative signature undoes the nominal filter response.
            signature[-1, 2] += 1.
            self.E = ep-gain[None, :, None]*signature[:, None, :]
            self.J += signature[:, :, None]*signature[:, None, :]/v
            self.b += signature*(innovation/v)
            self._detect()
        return self.output()

    def _detect(self):
        j = self.J; b = self.b
        sd = b[:, 1]**2/np.maximum(j[:, 1, 1], 1e-300)
        si = b[:, 2]**2/np.maximum(j[:, 2, 2], 1e-300)
        # Schur complement tests k conditional on an unrestricted d step.
        schur = j[:, 0, 0]-j[:, 0, 1]**2/np.maximum(j[:, 1, 1], 1e-300)
        bk = b[:, 0]-j[:, 0, 1]*b[:, 1]/np.maximum(j[:, 1, 1], 1e-300)
        identifiable = schur > 1e-10*np.maximum(j[:, 0, 0], 1.)
        sk = np.where(identifiable, bk**2/np.maximum(schur, 1e-300), 0.)
        joint = sd+sk
        scores = np.column_stack((sd-self.h1,
                                  si-self.h1 if self.impulse else np.full(len(sd), -np.inf),
                                  np.where((sk>self.h1)&identifiable, joint-self.h2, -np.inf)))
        flat = int(np.argmax(scores)); row, mode = np.unravel_index(flat, scores.shape)
        self.last.update(max_stat=float(max(np.max(sd), np.max(si) if self.impulse else 0., np.max(joint))),
                         max_k_evidence=float(np.max(sk)))
        if scores[row, mode] <= 0: return
        cols = [1] if mode == 0 else ([2] if mode == 1 else [0, 1])
        jj = j[row][np.ix_(cols, cols)]; bb = b[row, cols]
        cov = np.linalg.inv(jj); amplitude = cov@bb
        err = self.E[row][:, cols]
        self.x += err@amplitude
        self.P += err@cov@err.T
        self.P = (self.P+self.P.T)*.5
        event = dict(index=self.index, theta=int(self.theta[row]),
                     kind=('d_step', 'd_impulse', 'k_d_step')[mode],
                     statistic=float((sd[row], si[row], joint[row])[mode]),
                     k_evidence=float(sk[row]), amplitude=amplitude.tolist(),
                     sigma_after=float(self.scale*math.sqrt(max(self.P[0, 0], 0.))))
        self.events.append(event); self.last.update(detected=True, event=event)
        self.clear_bank()

    def output(self):
        return (float(self.scale*self.x[0]),
                float(self.scale*math.sqrt(max(self.P[0, 0], 0.))),
                float(np.mean(self.nis)) if self.nis else 1.)


class JGLRRidge(RidgeProjection):
    def __init__(self, q_k=.0005, L=96, alpha=1e-4, impulse=True, **kw):
        super().__init__(**kw)
        self.q_k=q_k; self.L=L; self.alpha=alpha; self.impulse=impulse

    def reset(self, policy):
        super().reset(policy)
        self.core=None; self.last_row=None; self.probe_cache={}; self.events=[]
        self.times=[]; self.segment=0; self.regressor_trace=[]

    def _probe_weights(self, m):
        if m not in self.probe_cache:
            # Six-tick causal local mean, constrained to have an exact zero
            # at the known carrier. Unlike the causal D4 endpoint residual,
            # I-S has unit gain AND zero phase at the probe frequency.
            # D4's endpoint residual is only .392 at 20.8 Hz, phase +99 deg.
            ticks=np.arange(6)-5
            a=np.vstack((np.ones(6),np.cos(self.w0*self.Tc*ticks),
                          np.sin(self.w0*self.Tc*ticks)))
            low=a.T@np.linalg.solve(a@a.T,np.array([1.,0.,0.]))
            w=np.zeros(m);w[-6:]=-low;w[-1]+=1.
            w[-6:]-=w.sum()/6.;self.probe_cache[m]=w
        return self.probe_cache[m]

    def estimate(self, o, n, policy):
        p=policy.params; log=self._log(o,n); nw=self.nw
        self.fit=_fit_record(degree=0,tw=nw*policy.Tc)
        j0=max(policy.ndel+1,n-1-nw); m=n-1-j0
        if m<MIN_TICKS:return np.nan,np.inf,np.inf
        _,y,z,keep=common._build(log,j0,m,0,policy.ns,policy.eb,policy.afc,policy.ndel,
                                log[n-1,0],log[n-1,2],nw*policy.ns*1e-3,p['enc_y'])
        broken=np.flatnonzero(~keep)
        if len(broken):y,z=y[broken[-1]+1:],z[broken[-1]+1:]
        if len(y)<MIN_TICKS:
            self.core=None;self.last_row=None;return np.nan,np.inf,np.inf
        row=n-2
        if self.last_row is None or row!=self.last_row+1:
            self.core=JumpFilter(noise=p['noise'],rho=self.rho,q_k=self.q_k,
                                 L=self.L,alpha=self.alpha,impulse=self.impulse,nw=nw)
            self.times=[];self.segment+=1
        self.last_row=row
        yt=float(self._probe_weights(len(y))@(y-y[-1]))
        self.times.append(float(log[row,0]))
        old=len(self.core.events);k,sig,chi=self.core.step(yt,float(z[-1]))
        self.regressor_trace.append((float(log[row,0]),yt,float(z[-1])))
        if len(self.core.events)>old:
            ev=dict(self.core.events[-1],segment=self.segment,
                    time=self.times[-1],callback_time=float(log[n-1,0]),
                    theta_time=self.times[self.core.events[-1]['theta']])
            self.events.append(ev)
        # sigma is the filter posterior SD, without the legacy .796 shrinkage.
        chi*=self.chi_scale; k=_clip_k(k)
        if np.isfinite(k) and sig<=.25*k:self.k_ref=k
        self.fit=_fit_record(k,sig,chi,degree=0,tw=nw*policy.Tc)
        self.fit['te']=float(log[n-1,0]);return k,sig,chi


class ContactJGLR(C.Preprocess,JGLRRidge):
    def __init__(self,kind='H',**kw):
        JGLRRidge.__init__(self,memory=False,**dict(C.RIDGE,**kw));self.kind=kind

    def reset(self,policy):
        JGLRRidge.reset(self,policy);self._setup(self.kind,policy)


def build(sc,cfg=None):
    import method3 as M
    import sweeps
    c=dict(cfg or {});c.setdefault('nw',24);c.setdefault('fc',14.)
    c['sigma_scale']=1.;c['chi_scale']=.551
    p=dict(sweeps.FINAL_C);p.update(C.CTL);p.update(M.TRACKER)
    return M.attach(C.CandidatePolicy('H',ContactJGLR(**c),single_model=False,**p),sc)
