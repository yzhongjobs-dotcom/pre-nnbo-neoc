"""Candidate 1: execution-aware parametric OCP and local solution reuse.

Derivative-based finite-horizon realization replacing the unimplemented HJB
critic for this experiment; not a reproduction of NNBO or the complete B0.
"""
import sys,time
from pathlib import Path
import numpy as np
import casadi as ca
from scipy.optimize import minimize
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mhe_guided_impedance_pilot'))
from core import m,LOGLO,LOGHI,encode,decode


def context(t,y,zr,q,anchor,sensed,cfg):
    return np.r_[t,y,zr,q,anchor,sensed,cfg[:9]]


class Problem:
    def __init__(self,dt=None,horizon=.24,nq=33,rate_actions=False,gains_only=False,true_plant=False,closed_projection=True):
        if dt is None:dt=.002 if rate_actions else .004
        self.dt=dt;self.horizon=horizon
        self.rate_actions=rate_actions
        self.gains_only=gains_only
        envdim=14 if true_plant else 7
        y=ca.SX.sym('y',10);q=ca.SX.sym('q',envdim);an=ca.SX.sym('an',2);t=ca.SX.sym('t');speed=ca.SX.sym('speed')
        u=ca.DM(-m.A+(np.arange(nq)+.5)*2*m.A/nq);weight=m.L*2*m.A/nq
        dx=(speed*t+y[0]+u-an[0])/.01;k=ca.exp(ca.fmin(19,ca.fmax(11,q[0]+q[1]*dx)))
        hh=.001*(q[2]+q[3]*dx+q[4]*(t-an[1]));delta=y[1]-u*u/(2*m.R)-hh
        rate=y[3]+u*(speed+y[2])/m.R-.001*q[4]
        kx=k*q[1]/.01;hx=.1*q[3];tau=ca.exp(q[5]);mu=q[6];drag=3.5
        if true_plant:
            xx=speed*t+y[0]+u;shape=ca.tanh((xx-q[2])/q[3]);k=q[0]*ca.exp(q[1]*shape)
            hh=q[4]*(ca.sin(2*np.pi*xx/q[5]+q[6])-ca.sin(q[6]))+q[7]*(ca.sin(2*np.pi*q[8]*t+q[6])-ca.sin(q[6]))
            hd=q[7]*2*np.pi*q[8]*ca.cos(2*np.pi*q[8]*t+q[6]);delta=y[1]-u*u/(2*m.R)-hh
            rate=y[3]+u*(speed+y[2])/m.R-hd;kx=k*q[1]*(1-shape**2)/q[3];hx=q[4]*2*np.pi/q[5]*ca.cos(2*np.pi*xx/q[5]+q[6])
            tau=q[12];mu=q[10]*ca.exp(q[11]*ca.tanh((speed*t+y[0]-q[2])/q[3]));drag=q[13]
        elastic=k*delta;pressure=ca.if_else(delta>0,ca.fmax(0,elastic+k*tau*rate),0)
        fn=ca.sum1(pressure)*weight;moment=ca.dot(pressure,u)*weight
        ft=ca.sum1(ca.if_else(delta>0,.5*kx*delta**2-elastic*hx+(pressure-elastic)*u/m.R,0))*weight
        ft+=mu*fn*ca.tanh((speed+y[2])/.001)+drag*(speed+y[2])
        contact=ca.Function('contact',[t,y,q,an,speed],[ca.vertcat(fn,ft,-moment+.025*ft)])
        a=ca.SX.sym('a',5);zr=ca.SX.sym('zr');sensed=ca.SX.sym('sensed',3);cfg=ca.SX.sym('cfg',9);force=ca.SX.sym('f',3)
        clip=lambda z,lo,hi:ca.fmin(hi,ca.fmax(lo,z))
        target=ca.exp(ca.DM(LOGLO)+ca.DM(LOGHI-LOGLO)*(clip(a[:4],-1,1)+1)/2)
        if rate_actions:target=clip(y[4:8]+.035*ca.DM(m.GRATE)*clip(a[:4],-1,1),ca.DM(m.GMIN),ca.DM(m.GMAX))
        raw=(target-y[4:8])/.035;lo=ca.fmax(-ca.DM(m.GRATE),-m.ALPHA*(y[4:8]-ca.DM(m.GMIN)))
        hi=ca.fmin(ca.DM(m.GRATE),m.ALPHA*(ca.DM(m.GMAX)-y[4:8]));rate0=clip(raw,lo,hi)
        aa=ca.vertcat(.5*y[0]**2,.5*(y[1]-zr)**2,0,0)
        pd=y[6]*y[2]**2+y[7]*y[3]**2+.5*(y[2]**2+y[3]**2)
        budget=m.CHI*pd+m.ALPHA*(y[8]-m.E_MIN)
        if rate_actions and closed_projection:
            # The inverse target map makes raw rate already box-feasible.
            # Only two nonnegative stiffness coefficients enter tank power;
            # energy projection only reduces their rates, so only lower bounds
            # can activate. This exactly solves that 2-D Euclidean projection.
            multiplier=ca.fmax(0,ca.dot(aa,rate0)-budget)/ca.fmax(ca.dot(aa,aa),1e-30)
            vv=rate0-multiplier*aa
            at0=clip((budget-aa[1]*lo[1])/ca.fmax(aa[0],1e-30),lo[0],rate0[0])
            at1=clip((budget-aa[0]*lo[0])/ca.fmax(aa[1],1e-30),lo[1],rate0[1])
            v0=ca.if_else(vv[0]<lo[0],lo[0],ca.if_else(vv[1]<lo[1],at0,vv[0]))
            v1=ca.if_else(vv[0]<lo[0],at1,ca.if_else(vv[1]<lo[1],lo[1],vv[1]))
            limited=ca.vertcat(v0,v1,rate0[2:4])
            gainrate=ca.if_else(ca.dot(aa,rate0)>budget,limited,rate0)
        else:
            left=0.;right=1.
            for j in range(2):right=ca.fmax(right,ca.if_else(aa[j]>1e-16,ca.fmax(raw[j]-lo[j],0)/ca.fmax(aa[j],1e-16),0))
            for _ in range(28):
                mid=(left+right)/2;excess=ca.dot(aa,clip(raw-mid*aa,lo,hi))>budget
                left=ca.if_else(excess,mid,left);right=ca.if_else(excess,right,mid)
            gainrate=ca.if_else(ca.dot(aa,rate0)>budget,clip(raw-right*aa,lo,hi),rate0)
        edot=m.CHI*pd-ca.dot(aa,gainrate);edot=ca.if_else(ca.logic_and(y[8]>=m.E_MAX,edot>0),0,edot)
        ut=sensed[1]-y[4]*y[0]-y[6]*y[2]+.5*cfg[7]
        un=3+cfg[4]*(3-sensed[0])+y[9]-y[5]*(y[1]-zr)-y[7]*y[3]+(0 if gains_only else 3*clip(a[4],-1,1))
        ut=clip(ut,-12,12);un=clip(un,-12,12)
        # ki=0 in the shared experiment, explicitly not a general integral model.
        dy=ca.vertcat(y[2],y[3],(ut-force[1]-.5*(y[2]+cfg[7]))/.8,(un-force[0]-.5*y[3])/.5,gainrate,edot,0)
        stage=((force[0]-3)/.3)**2+(y[0]/.0015)**2+cfg[8]*((y[1]-zr)/.0015)**2+.04*(y[3]/.03)**2+.015*((ut/3)**2+((un-3)/3)**2)
        dynamics=ca.Function('dyn',[y,force,a,zr,sensed,cfg],[dy,stage])
        self.gain_rates=ca.Function('gain_rates',[y,a,zr],[gainrate])
        f=contact(t,y,q,an,speed);dy,_=dynamics(y,f,a,zr,sensed,cfg)
        ym=y+.5*dt*dy;fm=contact(t+.5*dt,ym,q,an,speed);dm,loss=dynamics(ym,fm,a,zr,sensed,cfg)
        step=ca.Function('step',[t,y,q,an,speed,a,zr,sensed,cfg],[y+dt*dm,dt*loss,f])
        anidx=12+envdim;sidx=anidx+2;cfgidx=sidx+3
        c=ca.MX.sym('c',cfgidx+9);act=ca.MX.sym('act',10);yy=c[1:11];ss=c[sidx:sidx+3]
        n=int(round(horizon/dt));lag=int(round(.02/dt));buf=[ss for _ in range(lag)];cost=0.
        for j in range(n):
            if rate_actions and j%int(round(.01/dt))==0:ss=c[cfgidx+6]*buf[j%lag]+(1-c[cfgidx+6])*ss
            feedback=ss if rate_actions else buf[j%lag]
            yy,loss,ff=step(c[0]+j*dt,yy,c[12:anidx],c[anidx:sidx],c[cfgidx+7],act[:5] if j<n//2 else act[5:],c[11],feedback,c[cfgidx:cfgidx+9])
            buf[j%lag]=ff;cost+=loss/horizon
        grad=ca.gradient(cost,act);hes=ca.jacobian(grad,act)
        self.value=ca.Function('value',[act,c],[cost])
        self.vg=ca.Function('vg',[act,c],[cost,grad])
        self.vgh=ca.Function('vgh',[act,c],[cost,grad,hes])
        self.cross=ca.Function('cross',[act,c],[ca.jacobian(grad,c)])
        self.contact=contact;self.step=step;self.dynamics=dynamics

    def fun(self,a,c):
        v,g=self.vg(a,c);return float(v),np.asarray(g).ravel()


class ReusedOptimizer:
    """Local active-free sensitivity prediction + bounded merit correction.

    A local numerical mechanism; no global optimality, safety certificate or
    automatic inheritance of the B0 proofs. All preparation is timed.
    """
    def __init__(self,problem,reuse=True,full=False):
        self.problem=problem;self.reuse=reuse;self.full=full
        self.a=None;self.c=None;self.S=None;self.calls=0;self.log=[]

    def decide(self,c,seed):
        tic=time.perf_counter();pr=self.problem;reason='correction'
        if self.a is None:
            init=seed.copy();reason='initial_anchor'
        elif self.reuse and self.S is not None:
            delta=self.S@(c-self.c)
            init=np.clip(self.a+np.clip(delta,-.25,.25),-1,1)
        else:init=np.r_[self.a[5:],self.a[5:]]
        # The previous held action is always an admissible merit competitor.
        hold=seed if self.a is None else np.r_[self.a[5:],self.a[5:]]
        if float(pr.value(hold,c))<float(pr.value(init,c)):init=hold
        opt=minimize(lambda a:pr.fun(a,c),init,jac=True,bounds=[(-1,1)]*10,method='L-BFGS-B',options={'maxiter':12 if self.full or self.a is None else 2,'maxls':6,'ftol':1e-6})
        a=opt.x.copy();v,g=pr.vg(a,c);g=np.asarray(g).ravel()
        active=((a<=-1+1e-5)&(g>0))|((a>=1-1e-5)&(g<0));free=np.flatnonzero(~active)
        pg=g.copy();pg[active]=0.;regularizer=0.
        if self.reuse and self.calls%4==0:
            _,_,H=pr.vgh(a,c);H=np.asarray(H)
            S=np.zeros((10,len(c)))
            if len(free):
                Hf=H[np.ix_(free,free)];regularizer=max(1e-4,-np.linalg.eigvalsh(Hf).min()+1e-3)
                cross=np.asarray(pr.cross(a,c))[free]
                S[free]=-np.linalg.solve(Hf+regularizer*np.eye(len(free)),cross)
            self.S=S
        self.a=a;self.c=c.copy();self.calls+=1
        self.log.append({'seconds':time.perf_counter()-tic,'predicted_cost':float(v),'projected_gradient':float(np.linalg.norm(pg,np.inf)),'iterations':int(opt.nit),'regularizer':regularizer,'reason':reason})
        if pr.rate_actions:
            from rate_core import decode_rate
            return decode_rate(a[:5])
        return decode(a[:5])


class FastNeighbor:
    """Residual-difference neighboring update with cached inverse curvature.

    Directional parameter response replaces a dense cross-Jacobian; secants
    refresh curvature after accepted corrections. Merit rejection and projected
    gradients handle imperfect anchors/active bounds. This is a quasi-Newton
    realization of the model route, not a novel general-purpose NLP solver.
    """
    def __init__(self,problem):
        self.problem=problem;self.a=None;self.g=None;self.inverse=None;self.log=[]
    def decide(self,c,seed):
        tic=time.perf_counter();pr=self.problem;calls=0
        if self.a is None:
            opt=minimize(lambda a:pr.fun(a,c),seed,jac=True,bounds=[(-1,1)]*10,method='L-BFGS-B',options={'maxiter':12,'maxls':8,'ftol':1e-6})
            a=opt.x;v,g,H=pr.vgh(a,c);H=np.asarray(H);vals,vecs=np.linalg.eigh(.5*(H+H.T));self.inverse=(vecs/np.maximum(vals,.005))@vecs.T
            g=np.asarray(g).ravel();v=float(v);calls=int(opt.nfev)+1;prediction=False
        else:
            a=self.a.copy();v,g=pr.fun(a,c);calls+=1
            change=g-self.g;direction=-self.inverse@change
            direction=np.clip(direction,-.2,.2);trial=np.clip(a+direction,-1,1)
            vv,gg=pr.fun(trial,c);calls+=1;prediction=vv<v
            if prediction:a,v,g=trial,vv,gg
            for _ in range(2):
                pg=g.copy();pg[((a<=-1+1e-6)&(g>0))|((a>=1-1e-6)&(g<0))]=0.
                if np.linalg.norm(pg,np.inf)<1e-4:break
                direction=-self.inverse@pg
                if direction@pg>=0:direction=-pg
                norm=max(1.,np.max(abs(direction))/.2);direction/=norm
                accepted=False
                for scale in (1.,.25,.0625):
                    trial=np.clip(a+scale*direction,-1,1);vv,gg=pr.fun(trial,c);calls+=1
                    if vv<v-1e-10:
                        s=trial-a;yd=gg-g;ys=yd@s
                        if ys>1e-9:
                            rho=1/ys;V=np.eye(10)-rho*np.outer(s,yd)
                            self.inverse=V@self.inverse@V.T+rho*np.outer(s,s)
                            vals,vecs=np.linalg.eigh(self.inverse);self.inverse=(vecs*np.clip(vals,.001,2000.))@vecs.T
                        a,v,g=trial,vv,gg;accepted=True;break
                if not accepted:break
        self.a=a.copy();self.g=g.copy()
        pg=g.copy();pg[((a<=-1+1e-6)&(g>0))|((a>=1-1e-6)&(g<0))]=0.
        self.log.append({'seconds':time.perf_counter()-tic,'predicted_cost':float(v),'projected_gradient':float(np.linalg.norm(pg,np.inf)),'gradient_calls':calls,'prediction_accepted':bool(prediction)})
        if pr.rate_actions:
            from rate_core import decode_rate
            return decode_rate(a[:5])
        return decode(a[:5])
