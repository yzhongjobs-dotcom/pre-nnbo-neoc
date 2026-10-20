"""Causal two-region Hertz moving-window estimation with surface/bias arrival.

Zero process-noise constant-velocity surface instance of the B0 model; not a
full per-sample stochastic-state MHE. Weakly identifiable components retain
their arrival priors. No true material or future observation is accepted.
"""
import time
import numpy as np
from scipy.optimize import least_squares
import physics as p

SCALE=np.array([1e5,1e5,.2,.2,3.,3.,.01,.001,.001,.02,.02,.00002])
LOW=np.array([.15,.15,0.,0.,0.,0.,-8.,-2.,-5.,-3.,-3.,-3.])
HIGH=np.array([4.,4.,2.5,2.5,4.,4.,8.,2.,5.,3.,3.,3.])
ARRIVAL=np.array([.07,.5,.12,.3,.8,.8,.35,.12,.25,.6,.6,.6])
NOISE=np.array([.01,.01,.00002])

class MHE:
    def __init__(self,depth,s0,speed,initial_wrench):
        a=float(p.halfwidth(depth)); e0=initial_wrench[0]*4*p.R/(np.pi*p.L*a*a)
        mu=np.clip((-initial_wrench[1]-2*speed)/max(initial_wrench[0],.05),.02,.4)
        physical=np.array([e0,e0,mu,mu,2.,2.,s0+2*a,0.,0.,0.,0.,0.])
        self.v=np.clip(physical/SCALE,LOW+1e-6,HIGH-1e-6)
        self.depth=depth;self.s0=s0;self.speed=speed;self.last=0.;self.log=[]

    def estimate(self,now,observations,use_moment=True):
        tic=time.perf_counter();obs=[a for a in observations if a[0]<=now-.02+1e-10 and a[0]>=now-.42]
        if len(obs)<6:return self.descriptor(now),{'seconds':0.,'samples':len(obs)}
        obs=obs[-36:];tt=np.array([a[0] for a in obs]);xx=np.array([a[1] for a in obs]);yy=np.array([a[2] for a in obs])
        end=tt[-1];arrival=self.v.copy();arrival[7]+=arrival[8]*(end-self.last);self.last=end
        def unpack(v):
            ph=v*SCALE
            return np.r_[ph[:7],ph[7]-ph[8]*end,ph[8],self.s0,self.depth,self.speed],ph[9:]
        channels=[0,1,2] if use_moment else [0,1]
        def residual(v):
            eta,beta=unpack(v);pred=p.wrench(tt,xx,eta)
            data=((pred+beta-yy)/NOISE)[:,channels].ravel()
            return np.r_[data,(v-arrival)/ARRIVAL]
        fit=least_squares(residual,np.clip(arrival,LOW+1e-8,HIGH-1e-8),bounds=(LOW,HIGH),max_nfev=18,ftol=1e-5,xtol=1e-5)
        self.v=fit.x; info={'seconds':time.perf_counter()-tic,'samples':len(obs),'rms_whitened':float(np.sqrt(np.mean(residual(fit.x)[:len(tt)*len(channels)]**2))),
                           'evaluations':fit.nfev,'success':bool(fit.success)}
        self.log.append(info);return self.descriptor(now),info

    def descriptor(self,now):
        ph=self.v*SCALE
        gamma=ph[7]+ph[8]*(now-self.last)
        return np.r_[ph[:7],gamma,ph[8],self.s0+self.speed*now,self.depth,self.speed]
