"""Explicit task adaptations of NNBO and the one-step MPVIC construction.

These retain source mechanisms but are NOT exact full-paper replications.
NNBO: Eq23/37/41 and Newton reference, with f=Fd+u for nonzero force cost.
Jin: uniform-acceleration affine gain input -> box QP, per Eq13-25;
fixed robot inertia, measured force held over one step, common rate/energy plant.
"""
import time
import numpy as np
from scipy.optimize import least_squares,minimize
from core import m,plant


def calibrate(p,zr,seed):
    """Synthetic pre-identification measurements at the initial material point.

    Prescribed pose/velocity probes, not a physically timed calibration robot.
    True p generates measurements only; controller receives the fitted HC.
    All methods may use this same declared initial calibration information.
    """
    rng=np.random.default_rng(seed);rows=[]
    for z in np.linspace(.7*zr,1.3*zr,8):
        for vz in (-.006,0.,.006):
            y=np.zeros(10);y[1]=z;y[3]=vz
            rows.append((z,vz,plant.contact(0.,y,p,129)[0][0]+rng.normal(0,.03)))
    data=np.array(rows);scale=max(zr,1e-4)
    # Dimensionless parameterization avoids ill-conditioning between k,b,n.
    def fun(v):return (np.exp(v[0])+np.exp(v[1])*data[:,1])*(data[:,0]/scale)**v[2]-data[:,2]
    fit=least_squares(fun,[np.log(3),np.log(10),1.2],bounds=([-4,-4,1.],[5,8,2.5]))
    n=fit.x[2];return np.array([np.exp(fit.x[0])/scale**n,np.exp(fit.x[1])/scale**n,n,0.]),{'samples':len(data),'force_rms':float(np.sqrt(np.mean(fun(fit.x)**2))),'mechanically_timed':False}


class NNBO:
    def __init__(self,hc,zr,cfg,omega=10.,zeta=1.,learning=5.):
        self.k,self.b,self.n,self.surface=hc;self.r=zr;self.cfg=cfg.copy();self.omega=omega;self.zeta=zeta;self.learning=learning
        damping=self.b*max(zr-self.surface,1e-5)**self.n
        stiffness=self.k*self.n*max(zr-self.surface,1e-5)**(self.n-1)
        # Stabilizing task-matched warm start from the local scalar ARE; a
        # declared adaptation, not the paper's [20,20,20] simulation setting.
        q=cfg[8]*(.3/.0015)**2;g=1/max(damping,1e-4);aa=stiffness/max(damping,1e-4)
        P=(-aa+np.sqrt(aa*aa+g*g*q))/(g*g)
        self.w=np.array([max(P,1e-6),0.,0.]);self.reference=zr;self.failures=0;self.timings=[]

    def decide(self,t,y,sensed,dt=.08):
        tic=time.perf_counter();delta=max(y[1]-self.surface,1e-6);e=y[1]-self.r;g=1/(self.b*delta**self.n)
        control=-.5*g*(2*self.w[0]*e+self.w[1]*self.r)
        edot=-self.k/self.b+g*(3.+control)
        nu=np.array([2*e*edot,self.r*edot,0.]);q=self.cfg[8]*(.3/.0015)**2
        H=self.w@nu+q*e*e+control*control
        self.w-=dt*self.learning*nu*H/(1+nu@nu)**2
        x=y[1];ok=False
        def eq(z):
            d=z-self.surface
            if d<=1e-7:return np.inf,np.nan
            gg=1/(self.b*d**self.n);grad=2*self.w[0]*(z-self.r)+self.w[1]*self.r
            return sensed[0]-3.+.5*gg*grad,gg*self.w[0]-.5*self.n*gg*grad/d
        for _ in range(30):
            value,der=eq(x)
            if abs(value)<1e-7:ok=True;break
            if not np.isfinite(der) or abs(der)<1e-12:break
            step=-value/der
            for j in range(24):
                trial=x+step*2.**(-j);tv,_=eq(trial)
                if np.isfinite(tv) and abs(tv)<abs(value):x=trial;break
            else:break
        if ok:self.reference=x
        else:self.failures+=1
        acceleration=self.omega**2*(self.reference-y[1])-2*self.zeta*self.omega*y[3]
        desired=sensed[0]+.5*acceleration+.5*y[3]
        # Position-reference interface expressed through the shared force plant.
        ff=desired-3-self.cfg[4]*(3-sensed[0])-y[9]+y[5]*(y[1]-self.r)+y[7]*y[3]
        self.timings.append(time.perf_counter()-tic)
        return np.r_[self.cfg[:4],ff]


class OneStepVIC:
    def __init__(self,cfg,horizon=.08,regularizer=1e-7):
        self.cfg=cfg.copy();self.T=horizon;self.reg=regularizer;self.timings=[]
    def decide(self,t,y,zr,sensed):
        tic=time.perf_counter();out=self.cfg[:4].copy()
        # Paper mechanism is instantiated on the normal interaction axis;
        # scanning uses the same already tuned tangential servo as NNBO.
        for idx,e,v,force,mass in ((1,y[1]-zr,y[3],3.-sensed[0],.5),):
            # Source uniform-acceleration discretization, gain variables K,D.
            T=self.T;B=np.array([[-.5*T*T*e/mass,-.5*T*T*v/mass],[-T*e/mass,-T*v/mass]])
            z=np.array([e+T*v+.5*T*T*force/mass,v+T*force/mass])
            Q=np.diag([self.cfg[8]/.0015**2,.04/.03**2]);R=self.reg*np.diag([1.,10.])
            H=B.T@Q@B+R;b=B.T@Q@z
            lo=np.array([m.GMIN[idx],m.GMIN[2+idx]]);hi=np.array([m.GMAX[idx],m.GMAX[2+idx]])
            # Exact active-set enumeration for a strictly convex 2-D box QP,
            # avoiding optimizer tolerance dependence on physical gain units.
            candidates=[np.clip(-np.linalg.solve(H,b),lo,hi)]
            for fixed in (0,1):
                free=1-fixed
                for bound in (lo[fixed],hi[fixed]):
                    gg=np.zeros(2);gg[fixed]=bound;gg[free]=np.clip(-(b[free]+H[free,fixed]*bound)/H[free,free],lo[free],hi[free]);candidates.append(gg)
            gains=min(candidates,key=lambda gg:float(gg@H@gg+2*b@gg))
            out[idx]=gains[0];out[idx+2]=gains[1]
        self.timings.append(time.perf_counter()-tic)
        return np.r_[out,0.]
