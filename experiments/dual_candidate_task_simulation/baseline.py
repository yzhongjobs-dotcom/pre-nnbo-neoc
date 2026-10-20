"""Jin TMECH2023 MS predictive VIC: declared planar scanning adaptation.

Source E41 eq13/14 (uniform-acceleration, frozen external force), eq21/28
(fixed inertia, free K/D), eq29 (tanh bounds), eq32 (soft state penalty),
eq36--42 and Algorithm1 (constant/variable stiffness tank conditioning).
Uses L-BFGS on the same transformed finite-horizon problem, rather than the
paper's iLQR numerical solver. No environment truth/preview is accepted.
It is a method-level reconstruction and task adaptation, NOT a reproduction
of the Panda experiments or an exact implementation of the author's solver.
"""
import time
import numpy as np
import casadi as ca
from scipy.optimize import minimize
import physics as p

class JinMS:
    def __init__(self,regularizer=1e-5,horizon_steps=10):
        self.reg=np.broadcast_to(np.asarray(regularizer,dtype=float),(2,)).copy();self.H=horizon_steps;self.T=.015;self.previous=None
        self.tank=np.ones(2)*.1;self.low=.01;self.high=1.;self.kc=p.GMIN[:2].copy();self.log=[]
        mu=ca.MX.sym('mu',4*self.H);s=ca.MX.sym('s',6);z=s;cost=0.;dt=self.T
        for i in range(self.H):
            gains=ca.DM(p.GMIN)+(ca.tanh(mu[i*4:i*4+4])+1)/2*ca.DM(p.GMAX-p.GMIN)
            chunks=[]
            for axis in range(2):
                e,v,f=z[axis*3:axis*3+3][0],z[axis*3+1],z[axis*3+2]
                k=gains[axis];d=gains[axis+2];m=p.MASS[axis]
                acc=(f-k*e-d*v)/m
                ne=e+dt*v+.5*dt*dt*acc;nv=v+dt*acc
                qe=1/(.002 if axis==0 else .001)**2;qv=.04/.015**2
                cost+=dt*(qe*ne**2+qv*nv**2+.1*f*f+self.reg[axis]*((k/m)**2+(d/m)**2+(1/m)**2))
                distance=ca.sqrt(ca.fmax(ca.fabs(ne)/.006-1,0)**2+ca.fmax(ca.fabs(nv)/.06-1,0)**2+1e-20)
                cost+=dt*50*(1-ca.exp(-distance/2))
                chunks.extend([ne,nv,f])
            z=ca.vertcat(*chunks)
        cost+=.1*((z[0]/.002)**2+(z[3]/.001)**2)
        self.fun=ca.Function('jin_value_gradient',[mu,s],[cost,ca.gradient(cost,mu)])

    def objective(self,mu,s):
        v,g=self.fun(mu,s);return float(v),np.asarray(g).ravel()

    def decide(self,x,sensed,speed,period=.02):
        tic=time.perf_counter()
        # Source external-force state; desired normal force is an explicit
        # fixed feedforward in this scanning task, hence f_n - F_d.
        net=np.array([sensed[1]-p.BR*(x[2]+speed),sensed[0]-p.FD-p.BR*x[3]])
        state=np.array([x[0],x[2],net[0],x[1],x[3],net[1]])
        if self.previous is None:
            ratio=np.clip(2*(x[4:8]-p.GMIN)/(p.GMAX-p.GMIN)-1,-.98,.98);start=np.tile(np.arctanh(ratio),self.H)
        else:start=np.r_[self.previous[4:],self.previous[-4:]]
        opt=minimize(lambda a:self.objective(a,state),start,jac=True,method='L-BFGS-B',options={'maxiter':50,'ftol':1e-9,'gtol':1e-6,'maxls':15})
        self.previous=opt.x.copy();g=p.GMIN+(np.tanh(opt.x[:4])+1)/2*(p.GMAX-p.GMIN)
        e=x[:2];v=x[2:4];kv=g[:2]-self.kc
        dissipation=g[2:]*v*v+kv*e*v
        sigma1=np.where(((dissipation>0)&(self.tank>self.high))|((dissipation<0)&(self.tank<self.low)),0.,1.)
        sigma2=np.where(dissipation>0,1.,sigma1)
        g[:2]=self.kc+sigma2*kv;self.tank+=period*sigma1*dissipation
        # Same physical gain/rate/tank interface, as fast as allowed; no extra
        # first-order actuator is imposed selectively on this comparator.
        action=np.clip((g-x[4:8])/(period*p.RATE),-1.,1.)
        self.log.append({'seconds':time.perf_counter()-tic,'success':bool(opt.success),'iterations':int(opt.nit),'source_tank_min':float(self.tank.min())})
        return action
