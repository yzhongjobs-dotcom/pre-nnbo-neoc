"""B0-style gain-rate actions through the SAME bounded target actuator.

rate is implemented by inverse first-order target mapping, with target bounds
enforced BEFORE the existing rate/energy shield. No unbounded actuator request.
"""
import numpy as np
from numba import njit
from core import m,plant,local_contact


@njit(cache=True)
def target_for_rate(y,command):
    return np.minimum(m.GMAX,np.maximum(m.GMIN,y[4:8]+.035*command))


@njit(cache=True)
def advance_rate(t,y,p,command,ff,zr,sensed,cfg,duration,dt=.001,nq=65):
    yy=y.copy();metrics=np.zeros(8)
    for i in range(int(np.ceil(duration/dt-1e-10))):
        h=min(dt,duration-i*dt);f=plant.contact(t+i*dt,yy,p,nq)[0]
        d,_,_,_,_=m.dynamics(yy,f,target_for_rate(yy,command),ff,zr,sensed,cfg)
        mid=yy+.5*h*d;fm=plant.contact(t+i*dt+.5*h,mid,p,nq)[0]
        dm,ut,un,active,sat=m.dynamics(mid,fm,target_for_rate(mid,command),ff,zr,sensed,cfg)
        yy+=h*dm;metrics[0]+=h*m.stage(mid,fm,ut,un,zr,cfg[8])
        metrics[1]+=h*(fm[0]-3)**2;metrics[2]+=h*mid[0]**2;metrics[3]=max(metrics[3],fm[0])
        metrics[4]+=h*(fm[0]<.2);metrics[5]+=h*active;metrics[6]+=h*sat;metrics[7]=max(metrics[7],abs(mid[0]))
    return yy,metrics


@njit(cache=True)
def rollout_rate(t,y,zr,q,anchor,speed,cfg,a,p,true_model=True,horizon=.24,dt=.002):
    yy=y.copy();n=int(round(horizon/dt));lag=int(round(.02/dt));buf=np.zeros((lag,3));cost=0.
    sensed=q[7:10].copy();control_stride=int(round(.01/dt))
    for j in range(lag):buf[j]=q[7:10]
    for j in range(n):
        tj=t+j*dt;act=a[:5] if j<n//2 else a[5:];cmd=m.GRATE*np.clip(act[:4],-1,1);ff=3.*min(1.,max(-1.,act[4]))
        force=plant.contact(tj,yy,p,33)[0] if true_model else local_contact(tj,yy,q,anchor,speed)[0]
        if j%control_stride==0:sensed=cfg[6]*buf[j%lag]+(1-cfg[6])*sensed
        buf[j%lag]=force
        d,_,_,_,_=m.dynamics(yy,force,target_for_rate(yy,cmd),ff,zr,sensed,cfg)
        mid=yy+.5*dt*d;fm=plant.contact(tj+.5*dt,mid,p,33)[0] if true_model else local_contact(tj+.5*dt,mid,q,anchor,speed)[0]
        dm,ut,un,_,_=m.dynamics(mid,fm,target_for_rate(mid,cmd),ff,zr,sensed,cfg)
        yy+=dt*dm;cost+=dt*m.stage(mid,fm,ut,un,zr,cfg[8])/horizon
    return cost


@njit(cache=True)
def rate_gradient(t,y,zr,q,anchor,speed,cfg,a,p):
    value=rollout_rate(t,y,zr,q,anchor,speed,cfg,a,p);g=np.zeros(10);eps=.003
    for j in range(10):
        if j%5==4:continue
        aa=a.copy();bb=a.copy();aa[j]+=eps;bb[j]-=eps
        g[j]=(rollout_rate(t,y,zr,q,anchor,speed,cfg,aa,p)-rollout_rate(t,y,zr,q,anchor,speed,cfg,bb,p))/(2*eps)
    return value,g


def decode_rate(a):return np.r_[m.GRATE*np.clip(a[:4],-1,1),3*np.clip(a[4],-1,1)]
