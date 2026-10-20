"""Finite-horizon HJB critic and neighboring-root continuation.

Fixed degree-two basis/nodes; exact 4-D Hamiltonian minimisation eliminates
node actions. Newton matrix uses the envelope theorem (B0 appendix I.4).
This is a planar numerical HJB instance, not an exact HJB/safety certificate.
"""
import time
import numpy as np
from scipy.linalg import qr
import physics as p

def features(t,x):
    y=np.concatenate([np.asarray(t)[...,None]/p.HORIZON,(x-np.r_[np.zeros(4),p.GREF,.02])/p.XSCALE],axis=-1)
    vals=[np.ones(y.shape[:-1])]; grad=[np.zeros(y.shape)]
    for i in range(10):
        vals.append(y[...,i]); g=np.zeros_like(y);g[...,i]=1.;grad.append(g)
    for i in range(10):
        for j in range(i,10):
            vals.append(y[...,i]*y[...,j]);g=np.zeros_like(y);g[...,i]+=y[...,j];g[...,j]+=y[...,i];grad.append(g)
    phi=np.stack(vals,axis=-1); dp=np.stack(grad,axis=-2); factor=1-y[...,0]
    ps=phi*factor[...,None]
    pt=(-phi+factor[...,None]*dp[...,0])/p.HORIZON
    px=factor[...,None,None]*dp[...,1:]/p.XSCALE
    return ps,pt,px

class Critic:
    def __init__(self,eta,seed=15701):
        rng=np.random.default_rng(seed)
        t=rng.uniform(0,p.HORIZON*.98,1500)
        x=np.r_[np.zeros(4),p.GREF,.02]+rng.uniform(-1,1,(1500,9))*p.XSCALE*np.array([.6,.35,.2,.2,.3,.3,.3,.3,.3])
        ok=np.all(np.array([p.margins(tt,xx,eta)>0 for tt,xx in zip(t,x)]),axis=1)
        t=t[ok];x=x[ok];ps,_,_=features(t,x)
        _,_,idx=qr(ps.T,pivoting=True,mode='economic');idx=idx[:ps.shape[1]]
        self.t=t[idx];self.x=x[idx];self.ps,self.pt,self.px=features(self.t,self.x)
        self.w=np.zeros(len(idx));self.eta=eta.copy();self.C=None;self.log=[]

    def setup(self,eta):
        x=self.x;t=self.t
        self.base=p.terms(t,x,eta)
        self.pg=np.empty_like(x)
        for j in range(9):
            h=p.XSCALE[j]*1e-5; v=np.zeros(9);v[j]=h
            self.pg[:,j]=(p.terminal(x+v,eta)-p.terminal(x-v,eta))/(2*h)

    def residual(self,w):
        f0,ell,lo,hi,power,budget,_=self.base
        grad=self.pg+np.einsum('ikj,k->ij',self.px,w)
        lin=grad[:,4:8]*p.RATE;lin[:,:2]-=power*grad[:,8,None]
        a=p.qp(lin,lo,hi,power,budget)
        f=f0.copy();f[:,4:8]=a*p.RATE;f[:,8]-=np.sum(power*a[:,:2],axis=1)
        F=self.pt@w+ell+np.sum(p.RA*a*a,axis=1)+np.sum(grad*f,axis=1)
        J=self.pt+np.einsum('ikj,ij->ik',self.px,f)
        return F,J

    def solve(self,eta,start=None,maxiter=35):
        self.setup(eta);w=self.w.copy() if start is None else start.copy();history=[]
        for _ in range(maxiter):
            F,J=self.residual(w);norm=np.max(abs(F));history.append(float(norm))
            if norm<1e-6:break
            delta=np.linalg.lstsq(J,-F,rcond=1e-12)[0]
            best=w; value=norm
            for scale in (1.,.5,.25,.1,.03,.01):
                z=w+scale*delta;v=np.max(abs(self.residual(z)[0]))
                if v<value:best=z;value=v;break
            if value>=norm:break
            w=best
        F,J=self.residual(w)
        return w,{'residual':float(np.max(abs(F))),'condition':float(np.linalg.cond(J)),'iterations':len(history),'history':history}

    def anchor(self,eta):
        tic=time.perf_counter();w,info=self.solve(eta);self.w=w;self.eta=eta.copy()
        F,J=self.residual(w);C=np.empty((len(w),len(eta)))
        for j in range(len(eta)):
            h=max(abs(eta[j]),np.array([1e5,1e5,.1,.1,1.,1.,.01,.001,.001,.01,.001,.01])[j])*1e-5
            v=np.zeros(len(eta));v[j]=h
            self.setup(eta+v);fp=self.residual(w)[0];self.setup(eta-v);fm=self.residual(w)[0]
            C[:,j]=np.linalg.lstsq(J,-(fp-fm)/(2*h),rcond=1e-12)[0]
        self.setup(eta);self.C=C;info['seconds']=time.perf_counter()-tic
        return info

    def update(self,eta,method='neighbor'):
        tic=time.perf_counter()
        if method=='neighbor' and self.C is not None:start=self.w+self.C@(eta-self.eta)
        else:start=self.w
        w,info=self.solve(eta,start,maxiter=3 if method=='neighbor' else 35)
        self.w=w;self.eta=eta.copy();info['seconds']=time.perf_counter()-tic;self.log.append(info)
        return info

    def action(self,t,x,eta):
        _,_,px=features(np.array(t),x)
        pg=np.empty(9)
        for j in range(9):
            h=p.XSCALE[j]*1e-5;v=np.zeros(9);v[j]=h
            pg[j]=(p.terminal(x+v,eta)-p.terminal(x-v,eta))/(2*h)
        grad=pg+px.T@self.w
        _,_,lo,hi,power,budget,_=p.terms(t,x,eta)
        lin=grad[4:8]*p.RATE;lin[:2]-=power*grad[8]
        return p.qp(lin,lo,hi,power,budget)

if __name__=='__main__':
    eta=np.array([1e5,1.4e5,.15,.2,2.,3.,.025,0.,0.,0.,.001,.01])
    x,eta=p.initial(eta);c=Critic(eta);print(c.anchor(eta),flush=True)
    print('action',c.action(0,x,eta),flush=True)
