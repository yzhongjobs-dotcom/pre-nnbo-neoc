"""Causal slowly varying Hertz contact; a finite nodal value-critic carrier.

All algorithms solve the same frozen-current-descriptor finite-horizon family.
The physical plant continues to change between queries. No future descriptor is
used in controller prediction. Not a full R15 MHE/spatial-preview implementation.
"""
from pathlib import Path
import sys
import numpy as np

PREVIOUS=Path(__file__).resolve().parents[2]/'快速验证修复与闭环_2026-10-12'
sys.path.insert(0,str(PREVIOUS))
from bellman_carrier import Carrier,src,DREF,FD,T,MASS,DAMP,BR,RU,RATE


def force(d,eta):
    return eta*src.fast_force(0.,float(d),(1.,0.))


class UniformCarrier(Carrier):
    def __init__(self,nd=33,nv=17,steps=200):
        super().__init__(nd,nv,steps,focused=True)
        self.base_force=.02/(2*.03)*30000*np.pi*self.a*self.a/2

    def fn_nodes(self,t,eta):
        return self.base_force*eta[0]

    def model(self,eta):
        key=tuple(eta)
        if key in self.model_cache:return self.model_cache[key]
        d,v,k=self.x.T;f=self.base_force*eta[0]
        acc=(FD-f-k*(d-DREF)-(DAMP+BR)*v)/MASS
        geom=self.geometry(d+self.h*v,v+self.h*acc)
        stage=((d-DREF)/DREF)**2+(v/.008)**2+(f/FD-1)**2+.1*self.norm**2
        phi=((d-DREF)/DREF)**2+(f/FD-1)**2
        # Same stationary frozen model at every horizon slice, shared storage.
        result=([(geom,stage)]*self.N,phi)
        self.model_cache[key]=result
        return result

    def policy(self,V,t,state,eta):
        d,v,k=state[:3]
        f=force(d,eta[0])
        acc=(FD-f-k*(d-DREF)-(DAMP+BR)*v)/MASS
        geom=self.geometry(np.array([d+self.h*v]),np.array([v+self.h*acc]))
        # Every low-level call uses phase zero of the receding frozen horizon.
        return float(self.minimize(self.lines(V[1],geom),np.array([k]))[1][0])

    def greedy_system(self,V,eta):
        F=np.empty_like(V);P=[None]*self.N
        F[-1]=V[-1]-self.model(eta)[1]
        for i in range(self.N):
            b,P[i]=self.bellman(V[i+1],i,eta,True)
            F[i]=V[i]-b
        return F,P


def descriptor(t,duration,kind):
    phase=min(max(t/duration,0.),1.)
    if kind=='ramp':return 1.+.2*phase
    # Smooth rise then return, starts at the nominal initial descriptor.
    return 1.+.1*(1.-np.cos(2*np.pi*phase))


def advance(c,state,V,observed_eta,t,length,duration,kind,fixed=False,dt=.001):
    x=state.copy();rows=[];steps=round(length/dt);cost=0.
    for j in range(steps):
        now=t+j*dt;eta=descriptor(now,duration,kind)
        u=0. if fixed else c.policy(V,0.,x,(observed_eta,0.))
        f=force(x[0],eta);e=x[0]-DREF
        stage=(e/DREF)**2+(x[1]/.008)**2+(f/FD-1)**2+.1*((x[2]-150)/150)**2+RU*u*u
        rows.append([now,eta,observed_eta,*x,f,u,stage])
        def rhs(tt,y):
            ff=force(y[0],descriptor(tt,duration,kind));ee=y[0]-DREF
            return np.array([y[1],(FD-ff-y[2]*ee-(DAMP+BR)*y[1])/MASS,
                             u,.3*(DAMP+BR)*y[1]**2-.5*ee**2*u])
        aa=rhs(now,x);bb=rhs(now+dt/2,x+dt*aa/2)
        cc=rhs(now+dt/2,x+dt*bb/2);dd=rhs(now+dt,x+dt*cc)
        x+=dt*(aa+2*bb+2*cc+dd)/6;cost+=stage*dt
        if x[0]<=0:raise FloatingPointError('Lost physical contact')
    return x,np.asarray(rows),cost


def metrics(c,trace,cost,state,duration,kind):
    eta=descriptor(duration,duration,kind);f=force(state[0],eta)
    terminal=((state[0]-DREF)/DREF)**2+(f/FD-1)**2
    # Sampled verification is EVALUATOR work, never an acceptance step.
    ok=bool(np.all((trace[:,3]>=c.d[0])&(trace[:,3]<=c.d[-1])) and
            np.all(abs(trace[:,4])<=c.v[-1]) and
            np.all((trace[:,5]>=30-1e-8)&(trace[:,5]<=300+1e-8)) and
            np.all(trace[:,6]>=.001) and np.all(abs(trace[:,8])<=500+1e-8) and
            np.all((trace[:,7]>0)&(trace[:,7]<2*FD)))
    return dict(cost=float(cost+terminal),integrated_cost=float(cost),terminal_cost=float(terminal),
                force_rmse_N=float(np.sqrt(np.mean((trace[:,7]-FD)**2))),
                force_peak_N=float(np.max(trace[:,7])),minimum_energy=float(np.min(trace[:,6])),
                sampled_constraints_pass=ok,final_gain=float(state[2]))
