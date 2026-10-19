"""Trajectory-local quadratic HJB (iLQR) and neighboring feedback updates.

This is an explicit representation change from the failed global polynomial
collocation in critic.py. The backward pass returns a value gradient/Hessian
and feedback, not just an open-loop OCP action. No continuous HJB certificate.
"""
import time
import numpy as np
import casadi as ca
import physics as p

CENTRE=np.r_[np.zeros(4),p.GREF,.02]

def build(dt=.01):
    z=ca.SX.sym('z',9);u=ca.SX.sym('u',4);eta=ca.SX.sym('eta',12);t=ca.SX.sym('t')
    x=ca.DM(CENTRE)+z*ca.DM(p.XSCALE)
    def contact(tt,xx):
        delta=ca.fmax(1e-10,eta[10]-xx[1]-eta[7]-eta[8]*tt)
        a=ca.sqrt(p.R*delta)
        for _ in range(9):
            a=ca.fmin(p.LREF*.95,ca.fmax(1e-8,a));w=ca.sqrt(p.LREF**2-a*a)
            psi=(p.LREF**2-p.LREF*w+a*a*ca.log((p.LREF+w)/a))/(2*p.R)
            a=a-(psi-delta)/(a/p.R*ca.log((p.LREF+w)/a))
        a=ca.fmin(p.LREF*.95,ca.fmax(1e-8,a))
        rr=ca.fmin(a,ca.fmax(-a,eta[6]-eta[9]-eta[11]*tt-xx[0]));ww=ca.sqrt(ca.fmax(1e-30,a*a-rr*rr))
        i0=np.pi*a*a/2;iff=.5*(a*a*ca.acos(ca.fmin(1-1e-10,ca.fmax(-1+1e-10,rr/a)))-rr*ww)
        # Use the exact closed-coverage branches to avoid spurious derivatives.
        iff=ca.if_else(rr>=a,0,ca.if_else(rr<=-a,i0,iff))
        n0=p.L/(2*p.R)*eta[0]*(i0-iff);nf=p.L/(2*p.R)*eta[1]*iff
        fn=n0+nf;v=eta[11]+xx[2]
        ft=-ca.sign(v)*(eta[2]*n0+eta[3]*nf)-v*((1-iff/i0)*eta[4]+iff/i0*eta[5])
        return ca.vertcat(fn,ft)
    def rhs(tt,xx,aa):
        f=contact(tt,xx);g=xx[4:8];vel=xx[2:4];err=xx[:2]
        torque=-g[:2]*err-g[2:]*vel-ca.vertcat(0,p.FD)
        av=vel+ca.vertcat(eta[11],0);acc=(torque+ca.vertcat(f[1],f[0])-p.BR*av)/ca.DM(p.MASS)
        pd=ca.dot(g[2:],vel**2)+p.BR*ca.dot(av,av)
        # Forward projection is done by physics.qp at every physical step.
        dx=ca.vertcat(vel,acc,aa*ca.DM(p.RATE),p.CHI*pd-.5*ca.dot(err**2,aa[:2]*ca.DM(p.RATE[:2])))
        residual=ca.vertcat(err/ca.DM([.002,.001]),.2*vel/.015,(f[0]-p.FD)/.1,
                           np.sqrt(.08)*(g-ca.DM(p.GREF))/ca.DM(p.GSCALE),ca.sqrt(ca.DM(p.RA))*aa)
        return dx,residual
    dx,_=rhs(t,x,u);xm=x+.5*dt*dx;dx,r=rhs(t+.5*dt,xm,u)
    nextz=(x+dt*dx-ca.DM(CENTRE))/ca.DM(p.XSCALE)
    rx=ca.jacobian(r,z);ru=ca.jacobian(r,u)
    cost=dt*ca.dot(r,r)
    fn=ca.Function('local',[t,z,u,eta],[nextz,cost,ca.jacobian(nextz,z),ca.jacobian(nextz,u),
            ca.gradient(cost,z),ca.gradient(cost,u),2*dt*rx.T@rx,2*dt*ru.T@ru,2*dt*ru.T@rx])
    f=contact(p.HORIZON,x);rt=np.sqrt(.1)*ca.vertcat(x[0]/.002,x[1]/.001,(f[0]-p.FD)/.1);jx=ca.jacobian(rt,z)
    term=ca.Function('terminal',[z,eta],[ca.dot(rt,rt),2*jx.T@rt,2*jx.T@jx])
    return fn,term

class LocalHJB:
    def __init__(self,dt=.01):
        self.dt=dt;self.N=round(p.HORIZON/dt);self.fun,self.term=build(dt);self.log=[]

    def rollout(self,x0,eta,U,K=None,Zold=None,k=None,scale=1.):
        z=(x0-CENTRE)/p.XSCALE;Z=[z.copy()];V=[];total=0.;feasible=True
        for i in range(self.N):
            request=U[i].copy()
            if k is not None:request+=scale*k[i]+K[i]@(z-Zold[i])
            x=CENTRE+z*p.XSCALE
            _,_,lo,hi,power,budget,_=p.terms(i*self.dt,x,eta)
            a=p.qp(-2*p.RA*request,lo,hi,power,budget)
            vals=self.fun(i*self.dt,z,a,eta);z=np.asarray(vals[0]).ravel();total+=float(vals[1])
            V.append(a);Z.append(z.copy())
            if np.min(p.margins((i+1)*self.dt,CENTRE+z*p.XSCALE,eta)) < -1e-8:feasible=False
            if not np.isfinite(z).all() or np.max(abs(z))>100:return np.asarray(Z),np.asarray(V),1e20,False
        total+=float(self.term(z,eta)[0]);return np.asarray(Z),np.asarray(V),total,feasible

    def backward(self,Z,U,eta,reg=1e-4):
        vals=self.term(Z[-1],eta);vx=np.asarray(vals[1]).ravel();vxx=np.asarray(vals[2]);Vx=[vx];Vxx=[vxx];kk=[];KK=[]
        for i in range(self.N-1,-1,-1):
            data=self.fun(i*self.dt,Z[i],U[i],eta)
            _,_,A,B,lx,lu,lxx,luu,lux=[np.asarray(v) for v in data]
            lx=lx.ravel();lu=lu.ravel();qx=lx+A.T@vx;qu=lu+B.T@vx
            qxx=lxx+A.T@vxx@A;quu=luu+B.T@vxx@B+reg*np.eye(4);qux=lux+B.T@vxx@A
            feed=-np.linalg.solve(quu,qu);feedback=-np.linalg.solve(quu,qux)
            # A trust region stabilises nonlinear contact iterations; the
            # physical action is still the same exact gain/energy projection.
            feed=np.clip(feed,-.5,.5)
            vx=qx+feedback.T@quu@feed+feedback.T@qu+qux.T@feed
            vxx=qxx+feedback.T@quu@feedback+feedback.T@qux+qux.T@feedback
            vxx=.5*(vxx+vxx.T)
            Vx.append(vx);Vxx.append(vxx);kk.append(feed);KK.append(feedback)
        return np.array(kk[::-1]),np.array(KK[::-1]),np.array(Vx[::-1]),np.array(Vxx[::-1])

    def solve(self,x,eta,U=None,maxiter=20):
        tic=time.perf_counter()
        if U is None:U=np.zeros((self.N,4))
        Z,U,cost,feasible=self.rollout(x,eta,U);history=[cost]
        for _ in range(maxiter):
            k,K,_,_=self.backward(Z,U,eta);best=None
            for scale in (1.,.5,.25,.1,.03):
                trial=self.rollout(x,eta,U,K,Z,k,scale)
                if trial[2]<cost-1e-8 and (trial[3] or not feasible):best=trial;break
            if best is None:break
            Z,U,cost,feasible=best;history.append(cost)
        k,K,Vx,Vxx=self.backward(Z,U,eta)
        info={'cost':cost,'feasible':feasible,'iterations':len(history)-1,'seconds':time.perf_counter()-tic,
              'final_feedforward_max':float(np.max(abs(k))),'history':history}
        return {'Z':Z,'U':U,'K':K,'Vx':Vx,'Vxx':Vxx,'eta':eta.copy(),'x0':x.copy()},info

    def action(self,policy,phase,x):
        i=min(int(phase/self.dt),self.N-1);z=(x-CENTRE)/p.XSCALE
        return policy['U'][i]+policy['K'][i]@(z-policy['Z'][i])

class NeighborHJB:
    def __init__(self,solver):self.solver=solver;self.anchor_policy=None;self.sensitivity=None;self.cache=None
    def anchor(self,x,eta):
        tic=time.perf_counter();policy,info=self.solver.solve(x,eta);self.reanchor(policy)
        info['sensitivity_preparation_seconds']=time.perf_counter()-tic-info['seconds'];return policy,info
    def reanchor(self,policy):
        self.anchor_policy=policy
        self.cache=[tuple(np.asarray(v) for v in self.solver.fun(i*self.solver.dt,policy['Z'][i],policy['U'][i],policy['eta'])) for i in range(self.solver.N)]
    def predict(self,eta,x):
        """State/model directional neighboring-extremal Riccati sweep.

        Uses cached local quadratic data and new-minus-old model residuals;
        avoids 24 extra nonlinear solves to prepare a full FD sensitivity.
        It is an iLQR/Gauss--Newton neighboring approximation, not the exact
        implicit sensitivity or certificate of B0's fixed-node HJB root.
        """
        tic=time.perf_counter();pol=self.anchor_policy;N=self.solver.N
        newterm=self.solver.term(pol['Z'][-1],eta)
        dv=np.asarray(newterm[1]).ravel()-pol['Vx'][-1];dvectors=[dv.copy()];feeds=[];defects=[]
        for i in range(N-1,-1,-1):
            old=self.cache[i];new=tuple(np.asarray(v) for v in self.solver.fun(i*self.solver.dt,pol['Z'][i],pol['U'][i],eta))
            A=old[2];B=old[3];P=pol['Vxx'][i+1];lam=pol['Vx'][i+1]
            defect=(new[0]-old[0]).ravel();v=dv+P@defect
            qx=(new[4]-old[4]).ravel()+(new[2]-A).T@lam+A.T@v
            qu=(new[5]-old[5]).ravel()+(new[3]-B).T@lam+B.T@v
            h=old[7]+B.T@P@B+1e-4*np.eye(4)
            feed=-np.linalg.solve(h,qu);feeds.append(np.clip(feed,-.5,.5));defects.append(defect)
            dv=qx+pol['K'][i].T@qu;dvectors.append(dv.copy())
        feeds=feeds[::-1];defects=defects[::-1];dz=(x-pol['x0'])/p.XSCALE;Z=[pol['Z'][0]+dz];U=[]
        for i in range(N):
            du=feeds[i]+pol['K'][i]@dz;U.append(pol['U'][i]+du)
            dz=self.cache[i][2]@dz+self.cache[i][3]@du+defects[i];Z.append(pol['Z'][i+1]+dz)
        new={'Z':np.asarray(Z),'U':np.asarray(U),'K':pol['K'].copy(),'Vx':pol['Vx']+np.asarray(dvectors[::-1]),
             'Vxx':pol['Vxx'].copy(),'eta':eta.copy(),'x0':x.copy()}
        return new,{'prediction_seconds':time.perf_counter()-tic}

if __name__=='__main__':
    eta=np.array([1e5,1.4e5,.15,.2,2.,3.,.025,0.,0.,0.,.001,.01]);x,eta=p.initial(eta)
    solver=LocalHJB();pol,info=solver.solve(x,eta);print(info,flush=True)
