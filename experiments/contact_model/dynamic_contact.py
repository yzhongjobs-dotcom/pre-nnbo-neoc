"""Dynamic-consistent contact prototype. No plant, truth or scene imports.

Elastic pressure law remains FinitePatch/LinePatch. A one-relaxation branch
(passive in the frozen linear elastic case) is an explicit hypothesis, not
the known benchmark model/parameters:
    F = Fe + q; q_j = r q_(j-1) + alpha (Fe_j-Fe_(j-1)).
The inverse is stable for alpha>=0, 0<r<1. Hypotheses must be tested or
identified; this prototype does NOT claim to identify alpha/tau from one tone.
It does not prove global passivity with spatial changes/nonlinear elasticity.
"""
from dataclasses import dataclass
import math
import time
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter
from scipy.stats import chi2 as chi2_distribution
from finite_patch_model import FinitePatch, Geometry
from line_patch_model import LinePatch


@dataclass(frozen=True)
class Memory:
    alpha: float = .5
    tau: float = .2
    dt: float = .008

    @property
    def r(self):
        return math.exp(-self.dt/self.tau)

    def forward(self, elastic, initial_force=0., initial_elastic=0.):
        # DF2 state for F_j = b0 Fe_j + b1 Fe_prev + r F_prev.
        b0, b1 = 1+self.alpha, -(self.r+self.alpha)
        zi = np.array([b1*initial_elastic+self.r*initial_force])
        arr = np.asarray(elastic)
        return lfilter([b0,b1], [1.,-self.r], arr, zi=zi)[0]

    def inverse(self, forces, initial_force=None, initial_elastic=None):
        forces = np.asarray(forces, float)
        if initial_force is None: initial_force = float(forces[0])
        if initial_elastic is None: initial_elastic = initial_force
        b0 = 1/(1+self.alpha); b1 = -self.r/(1+self.alpha)
        pole = (self.r+self.alpha)/(1+self.alpha)
        zi = np.array([b1*initial_force+pole*initial_elastic])
        return lfilter([b0,b1], [1.,-pole], forces, zi=zi)[0]

    def operator(self, n):
        return lfilter([1+self.alpha,-(self.r+self.alpha)], [1.,-self.r], np.eye(n), axis=0)

    def excitation_gain(self, frequency=20.8):
        z = np.exp(-2j*np.pi*frequency*self.dt)
        return (1+self.alpha-(self.r+self.alpha)*z)/(1-self.r*z)


def project(B, y, tolerance=1e-7):
    # Normalise physical nuisance columns before rank selection. A slow
    # relaxation mode can be nearly cubic within a short window.
    norms = np.linalg.norm(B,axis=0)
    B = B[:,norms>1e-15]/norms[norms>1e-15]
    U,s,_ = np.linalg.svd(B,full_matrices=False)
    Q = U[:,s>tolerance*max(s[0],1.)]
    return y-Q@(Q.T@y), Q.shape[1]


class ContactFit:
    """A model-consistent force-domain residual, including arrival information.

    Elastic indentation is inverted from the causal elastic-force state. Force
    residual and surface nuisance are BOTH mapped through the same relaxation
    operator. The reported covariance/risk is nominal, not a safety guarantee.
    """
    def __init__(self, memory=Memory(), spatial=False, contact='foundation', moment=False,
                 arrival=True, noise=.04, max_nfev=20,data_confidence=True,blocks=4):
        self.memory, self.spatial, self.moment = memory, spatial, moment
        self.model = LinePatch(spatial) if contact=='line' else FinitePatch(Geometry(),spatial)
        self.arrival, self.noise, self.max_nfev = arrival, noise, max_nfev
        self.data_confidence=data_confidence
        self.previous = None
        self.blocks=int(blocks);self.history=[]
        self.retired_information=None;self.retired_vector=None

    def surface_profile(self,theta,obs,T,depth_init=None):
        """Profile smooth surface motion with the ORIGINAL forward force law.

        The inverse law is an initializer, never the final statistical objective.
        Forward residual is T(Fe_pred-Fe_observed); the carried memory initial
        state cancels because BOTH sides use the same causal inverse state.
        """
        x,y,Fe,B=(obs[k] for k in ('x','y','fe','B'))
        Fe=np.maximum(Fe,.02)
        d=self.model.inverse(x,Fe,theta) if depth_init is None else depth_init
        _,_,ke=self.model.response(x,d,theta)
        A=T@(ke[:,None]*B);u=T@(ke*(y-d))
        beta=np.linalg.lstsq(A,u,rcond=1e-8)[0]
        for _ in range(7):
            dep=np.maximum(y-B@beta,1e-9)
            pred,_,kp=self.model.response(x,dep,theta)
            r=T@(pred-Fe);J=-T@(kp[:,None]*B)
            change=np.linalg.lstsq(J,-r,rcond=1e-8)[0]
            old=float(r@r)
            for scale in (1.,.5,.25,.125):
                candidate=beta+scale*change
                fp=self.model.response(x,np.maximum(y-B@candidate,1e-9),theta,stiffness=False)[0]
                new=T@(fp-Fe)
                if float(new@new)<=old:
                    beta=candidate;break
            else:break
            if np.linalg.norm(scale*change)<1e-9:break
        dep=np.maximum(y-B@beta,1e-9)
        pred,moment,kp=self.model.response(x,dep,theta)
        _,rank=project(T@(kp[:,None]*B),np.zeros(len(y)))
        return T@(pred-Fe),pred,moment,beta,rank

    def fit(self, obs, fast_hint):
        tic=time.perf_counter()
        x,y,F,Fe,B=(obs[k] for k in ('x','y','f','fe','B'))
        Fe=np.maximum(Fe,.02); xe=obs['xend']; fe=float(Fe[-1])
        elastic_hint=np.clip(fast_hint/(1+self.memory.alpha),100.,3e5)
        lc=math.log(elastic_hint)
        if hasattr(self.model,'initial_log_scale'): lc=self.model.initial_log_scale(fe,elastic_hint)
        low=np.array([math.log(20.),-math.log(60),-6.,.1]) if self.spatial else np.array([math.log(20.)])
        high=np.array([math.log(1e11),math.log(60),6.,6.]) if self.spatial else np.array([math.log(1e11)])
        starts=[np.array([lc,math.log(3.),-.7,1.4]),np.array([lc,-math.log(3.),-.7,1.4])] if self.spatial else [np.array([lc])]
        prior=None; precision=None
        if self.history and obs['tend']-self.history[-1]['tend']>.3:
            self.history=[];self.retired_information=None;self.retired_vector=None
        self.history.append(obs)
        leaving=self.history.pop(0) if len(self.history)>self.blocks else None
        warm=None
        if self.arrival and self.previous is not None:
            prior,cov,oldx,oldt=self.previous; prior=prior.copy()
            if self.spatial: prior[2]+=(oldx-xe)/.01
            dt=max((obs['tend']-oldt)/.128,1.)
            cov=cov+np.diag([.12**2,.25**2,.04**2,.04**2] if self.spatial else [.12**2])*dt
            # Reopen a local patch after both old edges are outside support.
            if self.spatial and prior[2]+prior[3] < -3.:
                prior[1:]=(0.,.5,1.4); cov=np.diag([max(cov[0,0],.3**2),9.,16.,16.])
                self.retired_information=None;self.retired_vector=None
            ev,V=np.linalg.eigh(cov);precision=(V/np.sqrt(np.maximum(ev,1e-7)))@V.T
            starts=[prior,*starts]
            warm=prior.copy()
            if self.blocks>1:
                # Raw blocks are explicitly in the likelihood. Do not reuse
                # their posterior as another arrival cost (double counting).
                prior=None;precision=None
        T=self.memory.operator(len(F))
        def physical(z):
            return np.array([z[0],z[1],xe+.01*z[2],.01*z[3]]) if self.spatial else np.array([z[0]])
        def one_block(z,b):
            bx,bF,bFe=(b[k] for k in ('x','f','fe'));bFe=np.maximum(bFe,.02)
            bt=self.memory.operator(len(bF));th=physical(z);d=self.model.inverse(bx,bFe,th)
            reconstructed,_,_=self.model.response(bx,d,th)
            r,elastic_prediction,me,_,rank=self.surface_profile(th,b,bt,depth_init=d)
            # Inverse solvers have a declared finite indentation domain. A
            # saturated inverse is NOT a valid root; enforce the original force
            # equation so the optimiser cannot exploit a flat clipped depth.
            root_res=bt@(reconstructed-bFe)/self.noise
            rr=np.r_[r/self.noise,root_res]
            if self.moment:
                N=np.column_stack([np.ones(len(bF)),b['ft'],bF])
                predicted=(bF+r)/np.maximum(elastic_prediction,.02)*me
                rm,_=project(N,b['my']-predicted)
                rr=np.r_[rr,rm/.001]
            return rr,rank
        def components(z):
            parts=[one_block(z,b) for b in self.history]
            return np.concatenate([p[0] for p in parts]),sum(p[1] for p in parts)
        # Only a retired disjoint window enters the arrival information. The
        # information is maintained in GLOBAL cm coordinates for the boundary.
        npar=len(starts[0]);offset=np.zeros(npar)
        if self.spatial:offset[2]=xe/.01
        if self.retired_information is None:
            self.retired_information=np.zeros((npar,npar));self.retired_vector=np.zeros(npar)
        if self.blocks>1 and self.arrival:
            if np.any(self.retired_information):
                L=self.retired_information;mu=np.linalg.pinv(L,rcond=1e-10)@self.retired_vector
                dt=max((obs['tend']-(self.previous[3] if self.previous else obs['tend']))/.128,1.)
                Q=np.diag([.12**2,.25**2,.04**2,.04**2] if self.spatial else [.12**2])*dt
                propagated=L-L@np.linalg.solve(np.linalg.inv(Q)+L,L)
                self.retired_information=(propagated+propagated.T)/2
                self.retired_vector=self.retired_information@mu
            if leaving is not None:
                z0=warm if warm is not None else starts[0]
                r0,_=one_block(z0,leaving);J=np.empty((len(r0),npar))
                for j in range(npar):
                    h=1e-4;zp,zm=z0.copy(),z0.copy();zp[j]+=h;zm[j]-=h
                    J[:,j]=(one_block(zp,leaving)[0]-one_block(zm,leaving)[0])/(2*h)
                self.retired_information+=J.T@J
                self.retired_vector+=J.T@(J@(z0+offset)-r0)
        ev,V=np.linalg.eigh(self.retired_information)
        keep_arrival=ev>max(float(ev.max())*1e-10,1e-10)
        LA=np.sqrt(ev[keep_arrival])[:,None]*V[:,keep_arrival].T
        target_arrival=(V[:,keep_arrival].T@self.retired_vector)/np.sqrt(ev[keep_arrival])
        def data_arrival(z):return LA@(z+offset)-target_arrival
        def residual(z):
            rr,_=components(z)
            if prior is not None:return np.r_[rr,precision@(z-prior)]
            rr=np.r_[rr,data_arrival(z)]
            if self.spatial:rr=np.r_[rr,z[1]/3.,(z[3]-1.5)/5.]
            return rr
        best=None
        for st in starts:
            opt=least_squares(residual,np.clip(st,low+1e-6,high-1e-6),bounds=(low,high),max_nfev=self.max_nfev,
                              ftol=3e-4,xtol=3e-4,gtol=3e-4,x_scale='jac')
            val=float(opt.fun@opt.fun)
            if best is None or val<best[0]:best=(val,opt)
            # Once a warm start converges at the nominal data noise floor,
            # additional starts are not needed to force a lower training loss.
            if warm is not None and opt.success:
                rd,nr=components(opt.x);nu=max(sum(len(b['f']) for b in self.history)-nr-npar,1)
                if float(rd@rd)<=nu:break
        _,opt=best; th=physical(opt.x);rr,nuisance_rank=components(opt.x)
        H=opt.jac.T@opt.jac;cov=np.linalg.pinv(H,rcond=1e-10)
        # Root-consistency equations are constraints, not extra noisy samples.
        nobs=sum(len(b['f']) for b in self.history)*(2 if self.moment else 1)
        dof=max(nobs-nuisance_rank-len(opt.x),1)
        chi=float(rr@rr/dof)
        # Includes residual inflation, still not model-bias coverage.
        cov*=max(1.,chi)
        # Data information must be distinguished from regularisation. A prior
        # can make a Hessian invertible without observing the control gain.
        extra=data_arrival(opt.x) if prior is None else np.empty(0)
        JD=opt.jac[:len(rr)+len(extra),:]
        r_evidence=np.r_[rr,extra]
        U,sv,Vt=np.linalg.svd(JD,full_matrices=False)
        keep=sv>max(sv[0]*1e-8,1e-10)
        V=Vt.T
        data_cov=(V[:,keep]/sv[keep]**2)@V[:,keep].T*max(1.,chi)
        data_null=V[:,~keep]
        # First-order likelihood-profile displacement from the regularised
        # estimate: exact for the scalar linear limit, only local otherwise.
        shift=-(V[:,keep]/sv[keep])@(U[:,keep].T@r_evidence)
        if self.arrival:self.previous=opt.x.copy(),cov,xe,obs['tend']
        result=dict(theta=th,cov=cov,xe=xe,chi2=chi,dof=dof,success=bool(opt.success),
                    residual_rms=float(np.sqrt(np.mean(rr**2))),elapsed=time.perf_counter()-tic,
                    surface_beta=self.surface_profile(th,obs,T)[3],data_cov=data_cov,
                    data_null=data_null,data_shift=shift,data_rank=int(keep.sum()),
                    raw_blocks=len(self.history),arrival_rank=int(keep_arrival.sum()))
        prediction=self.predict(result,xe,fe)
        result.update(prediction)
        return result

    def predict(self, fit, x, elastic_force):
        th=fit['theta'];d=float(self.model.inverse(x,max(elastic_force,.02),th))
        ke=float(self.model.response(x,d,th)[2]);k0=ke*(1+self.memory.alpha)
        grad=np.empty(len(th))
        for i in range(len(th)):
            h=1e-4 if i<2 else 1e-6
            tp,tm=th.copy(),th.copy();tp[i]+=h;tm[i]-=h
            vals=[]
            for theta in (tp,tm):
                dd=self.model.inverse(x,max(elastic_force,.02),theta)
                vals.append(float(self.model.response(x,dd,theta)[2]))
            grad[i]=(math.log(max(vals[0],1.))-math.log(max(vals[1],1.)))/(2*h)
        if len(th)>1:grad[2:]*=.01  # covariance uses local cm coordinates
        C=fit['data_cov'] if self.data_confidence else fit['cov']
        rel=float(np.sqrt(max(grad@C@grad,0.)))
        unseen=bool(self.data_confidence and fit['data_null'].shape[1]>0 and
                    np.linalg.norm(fit['data_null'].T@grad)>1e-4*max(np.linalg.norm(grad),1e-12))
        bias=abs(float(grad@fit['data_shift'])) if self.data_confidence else 0.
        risk=rel*(1+math.sqrt(fit['dof']*max(fit['chi2']-1.,0.)))+bias
        if unseen:risk=float('inf')
        return dict(k=k0,elastic_k=ke,depth=d,relative_sd=rel,risk=risk,unobserved_gain=unseen,
                    regularization_log_shift=bias)

    def validate_contact(self, fit, obs, probability=.99):
        """Test a previously fitted contact law on the NEXT disjoint block.

        Contact parameters are frozen; only the unknown surface motion is
        profiled anew. Thus this is conditional contact validation, not a
        surface-motion forecast. The chi-square threshold is nominal: colored
        noise, model bias and weak excitation prevent a safety interpretation.
        """
        T=self.memory.operator(len(obs['f']))
        r,_,_,_,rank=self.surface_profile(fit['theta'],obs,T)
        dof=max(len(r)-rank,1)
        score=float(r@r/self.noise**2)
        threshold=float(chi2_distribution.ppf(probability,dof))
        reconstructed=self.model.response(obs['x'],self.model.inverse(obs['x'],obs['fe'],fit['theta']),
                                          fit['theta'],stiffness=False)[0]
        root_error=float(np.max(np.abs(reconstructed-obs['fe'])))
        return dict(valid=bool(np.isfinite(score) and score<=threshold and root_error<self.noise),
                    score=score,threshold=threshold,dof=dof,root_error=root_error,
                    force_rms=float(np.sqrt(np.mean(r*r))))

    def target_correction(self, fit, x, total_force, elastic_force, desired=4.):
        # A fast indentation change also changes the Maxwell branch, so the
        # target elastic force is Fe + (Fd-F)/(1+alpha), NOT Fd-q held frozen.
        target=max(elastic_force+(desired-total_force)/(1+self.memory.alpha),.02)
        d=self.model.inverse(x,max(elastic_force,.02),fit['theta'])
        dt=self.model.inverse(x,target,fit['theta'])
        return float(d-dt)


class Signals:
    """Streaming measurements only. Encoder endpoints and known force-filter spec."""
    def __init__(self,memory=Memory(),fc=30.,tau_servo=.019):
        self.memory,self.fc,self.tau_servo=memory,fc,tau_servo
        self.rows=[];self.previous=None;self.last_fe=None;self.last_force=None

    def append(self, row):
        # Seven allowed columns: t,x,y,filtered Fn,sent command,Ft,My.
        row=np.asarray(row,float)
        if row.shape!=(7,):raise ValueError('measurement-only seven-column interface')
        decay=math.exp(-2*math.pi*self.fc*self.memory.dt)
        if self.previous is None:
            f,ft,my=row[3],row[5],row[6];x,y=row[1],row[2]
        else:
            p=self.previous;denom=1-decay
            f=(row[3]-decay*p[3])/denom
            ft=(row[5]-decay*p[5])/denom;my=(row[6]-decay*p[6])/denom
            m=np.arange(8);af=1-math.exp(-2*math.pi*self.fc*.001)
            w=af*(1-af)**(7-m)
            fr=-np.expm1(-m*.001/self.tau_servo)/(-math.expm1(-self.memory.dt/self.tau_servo))
            q=float(w@fr/w.sum());x=p[1]+q*(row[1]-p[1]);y=p[2]+q*(row[2]-p[2])
        if self.last_fe is None or row[3]<.3:
            fe=f
        else:
            a=self.memory.alpha;r=self.memory.r
            fe=(f-r*self.last_force+(r+a)*self.last_fe)/(1+a)
        self.last_fe,self.last_force=float(fe),float(f)
        self.rows.append([row[0],x,y,f,fe,ft,my,row[3]])
        self.previous=row.copy()

    def window(self,n=16):
        if len(self.rows)<n:return None
        z=np.asarray(self.rows[-n:])
        if np.any(z[:,7]<=.5) or np.any(z[:,3]<=.1) or np.any(z[:,4]<=.02):return None
        tt=(z[:,0]-z[-1,0])/(n*self.memory.dt)
        return dict(tend=float(z[-1,0]),xend=float(z[-1,1]),fend=float(z[-1,3]),
                    x=z[:,1],y=z[:,2],f=z[:,3],fe=z[:,4],ft=z[:,5],my=z[:,6],
                    B=np.column_stack([tt**j for j in range(4)]))
