"""Reduced constrained moving-window parameter/state estimator, not B0 MHE.

Known robot pose/velocity history; local surface state and material parameters
jointly fit a delayed wrench window with a propagated arrival regularizer.
No true material, future samples, or true surface after initialization.
"""
import time
import numpy as np
from scipy.optimize import least_squares
from core import m, residuals


class Frontend:
    def __init__(self,zr,speed,moment=True):
        a=min(m.A,np.sqrt(2*m.R*max(zr,1e-6)))
        integral=m.L*(2*a*zr-a**3/(3*m.R))
        self.q=np.array([np.log(np.clip(3./max(integral,1e-10),3e5,3e7)),0,0,0,0,np.log(.008),.25])
        self.anchor=np.zeros(2);self.speed=speed;self.moment=moment
        self.times=[];self.states=[];self.forces=[];self.timings=[];self.fits=[]
        self.prior_sd=np.array([.15,.2,.15,.25,.5,.3,.06])

    def observe(self,t,y,force):
        self.times.append(t);self.states.append(y.copy());self.forces.append(force.copy())

    def update(self,i):
        j=i-2 # 20 ms actual sensing delay; no use of later force samples.
        if j<0 or i%5:return self.q,self.anchor
        tic=time.perf_counter();s=self.speed*self.times[j]+self.states[j][0]
        anchor=np.array([s,self.times[j]])
        prior=self.q.copy();dx=(s-self.anchor[0])/.01;dt=anchor[1]-self.anchor[1]
        prior[0]+=prior[1]*dx;prior[2]+=prior[3]*dx+prior[4]*dt
        low=np.array([12.,-2.,-20.,-3.,-10.,np.log(.001),.02])
        high=np.array([18.,2.,20.,3.,10.,np.log(.04),.8])
        prior=np.clip(prior,low+1e-8,high-1e-8)
        start=max(0,j-30);states=np.array(self.states[start:j+1]);times=np.array(self.times[start:j+1]);forces=np.array(self.forces[start:j+1])
        def fun(q):
            r=residuals(q,states,times,forces,anchor,self.speed,prior,self.prior_sd)
            if not self.moment:r[2:3*len(times):3]=0.
            return r
        fit=least_squares(fun,prior,bounds=(low,high),max_nfev=12,diff_step=1e-4,loss="soft_l1",f_scale=2.)
        self.q=fit.x;self.anchor=anchor
        self.timings.append(time.perf_counter()-tic);self.fits.append(float(np.sqrt(np.mean(fun(fit.x)[:3*len(times)]**2))))
        return self.q,self.anchor
