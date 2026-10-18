"""Two-DOF/four-gain stationary contact response, not the full B0 controller.

Known local mechanical coefficients, white or colored force disturbances.
Objective has displacement, contact-force, velocity and actuator-effort costs;
no gain preference, future disturbance realization or fitted image surrogate.
"""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov, solve_continuous_are
from scipy.optimize import minimize

MASS=np.array([1.,.25])
BE=np.diag([2.,1.])
SCALE=np.array([500.,500.,40.,40.])
LOW=np.array([5.,5.,2.,2.])/SCALE
HIGH=np.array([2000.,2000.,150.,150.])/SCALE
QX=np.eye(2)/.003**2
QV=np.eye(2)/.03**2
QF=np.eye(2)/.7**2
RU=np.eye(2)
WHITE=np.array([[.0009,.00027],[.00027,.0009]])
COLORED=np.array([[.04,.012],[.012,.04]])


def parameters(x,colored=False):
    """Features in [-1,1], fixed system family; values are illustrative."""
    x=np.asarray(x)
    kn=np.exp(np.log(300.)+(x[0]+1)*.5*np.log(900./300.))
    ratio=.20+(x[1]+1)*.15
    cross=.60*x[2]
    ke=np.array([[kn,cross*np.sqrt(ratio)*kn],
                 [cross*np.sqrt(ratio)*kn,ratio*kn]])
    decay=np.exp(np.log(4.)+(x[3]+1)*.5*np.log(40./4.)) if colored else None
    return ke,decay


class Problem:
    def __init__(self,x,colored=False,sensed=False):
        self.x=np.array(x);self.colored=colored
        self.sensed=sensed
        self.ke,self.decay=parameters(x,colored)
        self.base_n=6 if colored else 4
        self.n=self.base_n+(4 if sensed else 0)
        self.a0=np.zeros((self.n,self.n));self.a0[:2,2:4]=np.eye(2)
        self.a0[2:4,:2]=-self.ke/MASS[:,None]
        self.a0[2:4,2:4]=-BE/MASS[:,None]
        self.diffusion=np.zeros((self.n,self.n))
        if colored:
            self.a0[2:4,4:6]=np.diag(1/MASS)
            self.a0[4:6,4:6]=-self.decay*np.eye(2)
            self.diffusion[4:6,4:6]=2*self.decay*COLORED
        else:
            self.diffusion[2:4,2:4]=WHITE/MASS[:,None]/MASS[None,:]
        if sensed:
            # Finite-bandwidth observation error; never interpret continuous
            # white measurement noise as a finite instantaneous effort cost.
            rate=2*np.pi*30
            self.a0[self.base_n:,self.base_n:]=-rate*np.eye(4)
            noise_var=np.array([50e-6,50e-6,.002,.002])**2
            self.diffusion[self.base_n:,self.base_n:]=2*rate*np.diag(noise_var)
        self.contact_map=np.column_stack((self.ke,BE))
        self.basecost=np.zeros((self.n,self.n))
        self.basecost[:4,:4]=self.contact_map.T@QF@self.contact_map
        self.basecost[:2,:2]+=QX
        self.basecost[2:4,2:4]+=QV
        self.ais=[];self.lis=[]
        for i in range(4):
            li=np.zeros((2,self.n));li[i%2,i%2+(0 if i<2 else 2)]=SCALE[i]
            if sensed:li[i%2,self.base_n+i%2+(0 if i<2 else 2)]=SCALE[i]
            ai=np.zeros_like(self.a0);ai[2:4,:]=-li/MASS[:,None]
            self.ais.append(ai);self.lis.append(li)
        self.calls=0
        self.lyapunov_solves=0

    def matrices(self,g):
        l=sum(v*li for v,li in zip(g,self.lis))
        a=self.a0+sum(v*ai for v,ai in zip(g,self.ais))
        c=self.basecost+l.T@RU@l
        return a,c,l

    def evaluate(self,g,hessian=False):
        self.calls+=1
        self.lyapunov_solves+=10 if hessian else 2
        a,c,l=self.matrices(g)
        p=solve_continuous_lyapunov(a,-self.diffusion)
        v=solve_continuous_lyapunov(a.T,-c)
        cis=[li.T@RU@l+l.T@RU@li for li in self.lis]
        value=float(np.trace(c@p))
        grad=np.array([np.trace(ci@p)+2*np.trace(v@ai@p) for ai,ci in zip(self.ais,cis)])
        if not hessian:return value,grad
        h=np.zeros((4,4))
        for j,(aj,cj) in enumerate(zip(self.ais,cis)):
            pj=solve_continuous_lyapunov(a,-(aj@p+p@aj.T))
            vj=solve_continuous_lyapunov(a.T,-(aj.T@v+v@aj+cj))
            for i,(ai,ci) in enumerate(zip(self.ais,cis)):
                cij=self.lis[i].T@RU@self.lis[j]+self.lis[j].T@RU@self.lis[i]
                h[i,j]=np.trace(cij@p)+np.trace(ci@pj)+2*np.trace(vj@ai@p+v@ai@pj)
        return value,grad,(h+h.T)/2

    def metrics(self,g):
        a,c,l=self.matrices(g)
        p=solve_continuous_lyapunov(a,-self.diffusion)
        return dict(J=float(np.trace(c@p)),
            displacement_rms_m=np.sqrt(np.maximum(np.diag(p[:2,:2]),0)).tolist(),
            contact_force_rms_N=np.sqrt(np.maximum(np.diag(self.contact_map@p[:4,:4]@self.contact_map.T),0)).tolist(),
            actuator_rms_N=np.sqrt(np.maximum(np.diag(l@p@l.T),0)).tolist(),
            maximum_real_eigenvalue=float(np.max(np.linalg.eigvals(a).real)))

    def analytic(self):
        # Exact unbounded optimum in each independent white-force axis;
        # coupling/colored disturbances are deliberately ignored, not tuned out.
        ke=np.diag(self.ke);b=np.diag(BE)
        q=np.diag(QX)+np.diag(QF)*ke**2
        k=-ke+np.sqrt(ke**2+q/np.diag(RU))
        aa=(q+np.diag(RU)*k**2)/(ke+k)+(np.diag(QV)+np.diag(QF)*b*b)/MASS
        d=-b+np.sqrt(b*b+MASS*aa/np.diag(RU))
        return np.clip(np.r_[k,d]/SCALE,LOW,HIGH)

    def solve(self,starts=None):
        if starts is None:starts=[self.analytic()]
        out=[]
        for start in starts:
            sol=minimize(self.evaluate,np.clip(start,LOW,HIGH),jac=True,
                bounds=list(zip(LOW,HIGH)),method='L-BFGS-B',
                options=dict(ftol=1e-13,gtol=1e-9,maxiter=200,maxls=30))
            out.append(sol)
        return min(out,key=lambda s:s.fun)

    def care_reference(self):
        # Broader 2x4 feedback class, an exact white-noise lower bound. This
        # is NOT a diagonal-impedance implementation or colored-noise optimum.
        a=self.a0[:4,:4];b=np.zeros((4,2));b[2:,:]=np.diag(1/MASS)
        q=self.basecost[:4,:4]
        ric=solve_continuous_are(a,b,q,RU)
        l=np.linalg.solve(RU,b.T@ric)
        aa=self.a0.copy();aa[2:4,:4]-=l/MASS[:,None]
        ll=np.zeros((2,self.n));ll[:,:4]=l
        if self.sensed:
            ll[:,self.base_n:]=l
            aa[2:4,self.base_n:]-=l/MASS[:,None]
        cc=self.basecost+ll.T@RU@ll
        p=solve_continuous_lyapunov(aa,-self.diffusion)
        return float(np.trace(cc@p)),l


def projected_gradient(g,grad):
    return g-np.clip(g-grad,LOW,HIGH)


def newton(problem,start,frozen_hessian=None,maxiter=15,tol=1e-6):
    g=np.clip(np.array(start),LOW,HIGH);steps=0
    for steps in range(maxiter+1):
        if frozen_hessian is None:value,grad,h=problem.evaluate(g,True)
        else:value,grad=problem.evaluate(g);h=frozen_hessian
        residual=float(np.max(np.abs(projected_gradient(g,grad))))
        if residual<=tol or steps==maxiter:break
        active=((g<=LOW+1e-10)&(grad>0))|((g>=HIGH-1e-10)&(grad<0))
        free=~active;direction=np.zeros(4)
        hf=h[np.ix_(free,free)]
        if not free.any():break
        floor=max(0.,1e-7-np.linalg.eigvalsh(hf).min())
        direction[free]=np.linalg.solve(hf+floor*np.eye(free.sum()),-grad[free])
        accepted=False
        for fraction in .5**np.arange(20):
            trial=np.clip(g+fraction*direction,LOW,HIGH)
            if np.max(np.abs(trial-g))<1e-13:continue
            jt,_=problem.evaluate(trial)
            if jt<=value+1e-4*grad@(trial-g):g=trial;accepted=True;break
        if not accepted:break
    return g,dict(iterations=steps,projected_gradient=residual,converged=residual<=tol)
