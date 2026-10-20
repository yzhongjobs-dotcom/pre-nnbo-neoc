"""Exploratory local finite-footprint model and executed-action rollouts.

Imports the previously corrected plant unchanged. No full B0/HJB claim.
q = [log k, log-k change/cm, h/mm, h change/mm/cm, hdot/mm/s, log tau, mu].
The local field is affine in h and log-affine in k; not the true sinusoid/tanh.
"""
import sys
from pathlib import Path
import numpy as np
from numba import njit

sys.path.append(str(Path(__file__).resolve().parents[1]/"within_scan_impedance_feedback"))
import model as m
import coupled_contact as plant

LO=np.array([40., 1., 2., 2.])
HI=m.GMAX.copy()
LOGLO=np.log(LO); LOGHI=np.log(HI)


@njit(cache=True)
def local_contact(t,y,q,anchor,speed,nq=33):
    s=speed*t+y[0];vx=speed+y[2];fn=0.;moment=0.;ft=0.;kn=0.
    tau=np.exp(q[5]);w=m.L*2*m.A/nq
    for j in range(nq):
        u=-m.A+(j+.5)*2*m.A/nq
        dx=(s+u-anchor[0])/.01
        k=np.exp(min(19.,max(11.,q[0]+q[1]*dx)))
        h=.001*(q[2]+q[3]*dx+q[4]*(t-anchor[1]))
        delta=y[1]-u*u/(2*m.R)-h
        if delta>0:
            rate=y[3]+u*vx/m.R-.001*q[4]
            elastic=k*delta;pressure=max(0.,elastic+k*tau*rate)
            fn+=pressure*w;moment+=pressure*u*w
            ft+=(.5*k*q[1]/.01*delta*delta-elastic*.1*q[3]+(pressure-elastic)*u/m.R)*w
            if pressure>0:kn+=k*w
    ft+=q[6]*fn*np.tanh(vx/.001)+3.5*vx
    return np.array([fn,ft,-moment+.025*ft]),kn


@njit(cache=True)
def decode(a):
    # Bounded targets; rate and energy limits are applied to actual gains later.
    out=np.zeros(5)
    out[:4]=np.exp(LOGLO+(LOGHI-LOGLO)*(np.clip(a[:4],-1.,1.)+1.)/2.)
    out[4]=3.*min(1.,max(-1.,a[4]))
    return out


def encode(target):
    return np.r_[2*(np.log(np.clip(target[:4],LO,HI))-LOGLO)/(LOGHI-LOGLO)-1, target[4]/3.]


@njit(cache=True)
def rollout(t,y,zr,q,anchor,speed,cfg,a,p,true_model=False,horizon=.24,dt=.002):
    """Two held targets, genuine gain dynamics and delayed simulated force feedback.

    Model or privileged plant is used for *training* costs only. During deployment
    future true p never reaches the actor or model optimizer.
    """
    yy=y.copy();n=int(round(horizon/dt));lag=max(1,int(round(.02/dt)))
    buf=np.zeros((lag,3));cost=0.
    f0=q[7:10] if len(q)>=10 else local_contact(t,yy,q,anchor,speed)[0]
    for j in range(lag):buf[j]=f0
    for j in range(n):
        tj=t+j*dt;target=decode(a[:5] if j<n//2 else a[5:])
        force=plant.contact(tj,yy,p,33)[0] if true_model else local_contact(tj,yy,q,anchor,speed)[0]
        sensed=buf[j%lag].copy();buf[j%lag]=force
        dy,_,_,_,_=m.dynamics(yy,force,target[:4],target[4],zr,sensed,cfg)
        mid=yy+.5*dt*dy
        fm=plant.contact(tj+.5*dt,mid,p,33)[0] if true_model else local_contact(tj+.5*dt,mid,q,anchor,speed)[0]
        dm,ut,un,_,_=m.dynamics(mid,fm,target[:4],target[4],zr,sensed,cfg)
        yy+=dt*dm
        cost+=dt*m.stage(mid,fm,ut,un,zr,cfg[8])
        if not np.isfinite(cost) or cost>1e8:return 1e8
    return cost/horizon


@njit(cache=True)
def cost_gradient(t,y,zr,q,anchor,speed,cfg,a,p,true_model=False):
    value=rollout(t,y,zr,q,anchor,speed,cfg,a,p,true_model)
    grad=np.zeros(10);eps=.003
    for j in range(10):
        aa=a.copy();bb=a.copy();aa[j]+=eps;bb[j]-=eps
        grad[j]=(rollout(t,y,zr,q,anchor,speed,cfg,aa,p,true_model)-rollout(t,y,zr,q,anchor,speed,cfg,bb,p,true_model))/(2*eps)
    return value,grad


@njit(cache=True)
def residuals(q,states,times,forces,anchor,speed,prior,prior_sd):
    n=len(times);r=np.zeros(n*3+7)
    scale=np.array([.05,.05,.00015])
    for i in range(n):
        f=local_contact(times[i],states[i],q,anchor,speed)[0]
        r[3*i:3*i+3]=(f-forces[i])/scale
    r[n*3:]=(q-prior)/prior_sd
    return r


def features(t,y,zr,q,anchor,speed,cfg):
    f,k=local_contact(t,y,q,anchor,speed)
    # Causal MHE coefficients, current robot states, executed gains, task weights.
    return np.r_[y[:4]/[.0015,.01,.02,.03],(y[1]-zr)/.0015,
                 (np.log(np.clip(y[4:8],LO,HI))-LOGLO)/(LOGHI-LOGLO),y[8]/.03,
                 (q-np.array([15.5,0,0,0,0,-4.8,.25]))/[1,1,10,1,2,1,.2],
                 (f-[3.,.8,0])/[1.,1.,.03],k/3000.,speed/.012,cfg[8]]


def analytic(t,y,zr,q,anchor,speed,cfg,omega=35.,zeta=.7):
    f,k=local_contact(t,y,q,anchor,speed);k=max(k,200.)
    h=.0001;yp=y.copy();yp[1]+=h
    k=max(200.,(local_contact(t,yp,q,anchor,speed)[0][0]-f[0])/h)
    # Temporal/spatial drift at held robot penetration, not tissue truth.
    drift=(local_contact(t+.005,y,q,anchor,speed)[0][0]-f[0])/.005
    ratio=cfg[8]*(.3/.0015)**2
    elastic=f[0]-k*np.exp(q[5])*y[3]
    zs=y[1]+(k*(3.-elastic)+ratio*(zr-y[1]))/(k*k+ratio)
    vs=-drift*k/(k*k+ratio)
    desired=f[0]+.5*(omega**2*(zs-y[1])+2*zeta*omega*(vs-y[3]))+.5*y[3]
    ff=desired-3.-cfg[4]*(3.-f[0])-y[9]+y[5]*(y[1]-zr)+y[7]*y[3]
    return np.r_[cfg[:4],np.clip(ff,-3,3)]
