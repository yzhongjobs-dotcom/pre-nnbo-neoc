"""Independent constrained direct-shooting reference on the same physical OCP.

Not a replacement for NEOC/HJB: solely checks the value-critic policy quality.
Decision variables are actual stiffness knots with explicit rate constraints.
Forward analytic sensitivities differentiate the RK4 physical rollout.
"""
import math
import time
import numpy as np
from scipy.optimize import minimize,LinearConstraint,Bounds
from bellman_carrier import src,DREF,FD,T,MASS,DAMP,BR,RU

def force_slope(t,d,eta):
    a=src.fast_width(d)
    r=eta[1]*.003+.006-.008*t
    z=max(-1.,min(1.,r/a));theta=math.acos(z)
    fn=src.fast_force(t,d,eta)
    slope=(.02/(2*.03)*30000*.03*(math.pi+(eta[0]-1)*theta)/
           math.log((.02+math.sqrt(.02**2-a*a))/a))
    return fn,slope

class ShootingReference:
    def __init__(self,eta,segments=30,dt=.001):
        self.eta=eta;self.M=segments;self.dt=dt;self.times=np.linspace(0,T,segments+1)
        self.cached=None
    def gain(self,t,z):
        i=min(int(t/(T/self.M)),self.M-1);beta=min(max((t-self.times[i])/(T/self.M),0),1)
        k=np.r_[150.,z]
        grad=np.zeros(self.M);rate_grad=np.zeros(self.M)
        if i>0:grad[i-1]=1-beta;rate_grad[i-1]=-self.M/T
        grad[i]=beta;rate_grad[i]=self.M/T
        return (1-beta)*k[i]+beta*k[i+1],(k[i+1]-k[i])*self.M/T,grad,rate_grad
    def fun(self,z,save=False):
        if self.cached is not None and np.array_equal(z,self.cached[0]) and not save:
            return self.cached[1],self.cached[2]
        n=self.M;x=np.array([DREF,0.]);S=np.zeros((2,n));cost=0.;gradient=np.zeros(n);rows=[]
        h=self.dt;steps=round(T/h)
        def flow(t,x,S):
            k,u,dk,du=self.gain(t,z);f,fd=force_slope(t,x[0],self.eta)
            e=x[0]-DREF
            fx=np.array([[0.,1.],[-(fd+k)/MASS,-(DAMP+BR)/MASS]])
            ff=np.array([x[1],(FD-f-k*e-(DAMP+BR)*x[1])/MASS])
            ds=fx@S+np.outer([0.,-e/MASS],dk)
            stage=(e/DREF)**2+(x[1]/.008)**2+(f/FD-1)**2+.1*((k-150)/150)**2+RU*u*u
            lx=np.array([2*e/DREF**2+2*(f/FD-1)*fd/FD,2*x[1]/.008**2])
            grad=lx@S+.2*(k-150)/150**2*dk+2*RU*u*du
            return ff,ds,stage,grad,k,u,f
        for i in range(steps):
            t=i*h
            a,Sa,la,ga,k,u,f=flow(t,x,S)
            b,Sb,lb,gb,*_=flow(t+h/2,x+h*a/2,S+h*Sa/2)
            c,Sc,lc,gc,*_=flow(t+h/2,x+h*b/2,S+h*Sb/2)
            d,Sd,ld,gd,*_=flow(t+h,x+h*c,S+h*Sc)
            if save:rows.append([t,x[0],x[1],k,u,f,la])
            x+=h*(a+2*b+2*c+d)/6;S+=h*(Sa+2*Sb+2*Sc+Sd)/6
            cost+=h*(la+2*lb+2*lc+ld)/6;gradient+=h*(ga+2*gb+2*gc+gd)/6
            if x[0]<=0:raise FloatingPointError('Reference lost contact')
        f,fd=force_slope(T,x[0],self.eta);e=x[0]-DREF
        cost+=(e/DREF)**2+(f/FD-1)**2
        gradient+=np.array([2*e/DREF**2+2*(f/FD-1)*fd/FD,0.])@S
        self.cached=(z.copy(),float(cost),gradient.copy())
        return (cost,gradient,np.array(rows)) if save else (float(cost),gradient)
    def solve(self):
        D=np.eye(self.M)
        D[np.arange(1,self.M),np.arange(self.M-1)]=-1
        offset=np.zeros(self.M);offset[0]=150.
        rate=500*T/self.M
        limit=LinearConstraint(D,offset-rate,offset+rate)
        start=time.perf_counter()
        fit=minimize(lambda z:self.fun(z),np.full(self.M,150.),jac=True,method='SLSQP',
                     bounds=Bounds(np.full(self.M,30.),np.full(self.M,300.)),constraints=[limit],
                     options=dict(maxiter=150,ftol=1e-10,disp=False))
        val,grad,data=self.fun(fit.x,True)
        # Feasible local solution; no global-optimality claim for nonlinear shooting.
        return fit,data,dict(cost=float(val),success=bool(fit.success),message=fit.message,
                             iterations=int(fit.nit),solve_seconds=time.perf_counter()-start,
                             maximum_rate=float(np.max(abs(D@fit.x-offset))/(T/self.M)),
                             minimum_gain=float(min(fit.x)),maximum_gain=float(max(fit.x)))
