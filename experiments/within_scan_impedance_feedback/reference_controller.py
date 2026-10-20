"""Simple equilibrium-reference compensation using the shared local response.

A principle-based analytic baseline, not an implementation of a named paper.
Same actuator, observation, gain dynamics and sensor/decision periods as model.py.
"""
import numpy as np
from numba import njit
import model as m


@njit(cache=True)
def decide(y,zref,cfg,resp,temporal):
    k=max(200.,resp[0]);visc=.008*k
    # Minimize the local static force/normal-pose quadratic in z, not a gain label.
    qratio=cfg[8]*(.3/.0015)**2
    elastic_force=resp[3]-visc*y[3]
    zstar=y[1]+(k*(m.FD-elastic_force)+qratio*(zref-y[1]))/(k*k+qratio)
    source=cfg[9]*(resp[1]*(cfg[7]+y[2])+resp[2])+(1-cfg[9])*temporal
    vstar=-source*k/(k*k+qratio)
    omega,zeta=cfg[10],cfg[11]
    desired=resp[3]+m.MASS[1]*(omega*omega*(zstar-y[1])+2*zeta*omega*(vstar-y[3]))+.5*y[3]
    # Express a reference/force correction through the same additive normal input.
    ff=desired-m.FD-cfg[4]*(m.FD-resp[3])-y[9]+y[5]*(y[1]-zref)+y[7]*y[3]
    out=np.zeros(5);out[:4]=cfg[:4];out[4]=min(3.,max(-3.,ff))
    return out


@njit(cache=True)
def episode(p,noise,cfg0,latency=0.,dt=m.DT,nq=33,geometry_prior=True):
    cfg=cfg0.copy();cfg[7]=p[9]
    y=m.initial(p,cfg[:4]);zr=y[1];count=len(noise)-1
    hist=np.zeros((count+1,8));trace=np.zeros((count,26));metrics=np.zeros(10);metrics[8]=1e3
    k=1000.
    if geometry_prior:
        k=1.5*m.FD/zr if zr<m.A*m.A/(2*m.R) else m.FD/(zr-m.A*m.A/(6*m.R))
        k=min(7000.,max(200.,k))
    old=np.array([k,0.,0.,3.,.6,zr,0.,0.]);sensed=np.array([3.,.6,0.])
    target=cfg[:4].copy();ff=0.;pending=target.copy();pendingff=0.;release=-1.
    for i in range(count):
        t=i*.01;oldy=y.copy();force,_=m.wrench(t,y,p,nq)
        hist[i,:4]=y[:4];hist[i,4:7]=force+noise[i];hist[i,7]=p[9]*t+y[0]
        resp,temp,_=m.response(hist,i,old);old=resp.copy()
        observation=hist[max(0,i-2),4:7].copy();observation[:2]=resp[3:5]
        sensed=cfg[6]*observation+(1-cfg[6])*sensed;resp[3:5]=sensed[:2]
        if i%4==0:
            proposal=decide(y,zr,cfg,resp,temp);pending=proposal[:4].copy();pendingff=proposal[4];release=t+latency
        if release>=0 and release<=t+1e-12:
            target=pending.copy();ff=pendingff;release=-1.
        if release>t and release<t+.01:
            first=release-t
            y,aa=m.advance(t,y,p,target,ff,zr,sensed,cfg,first,dt,nq)
            target=pending.copy();ff=pendingff;release=-1.
            y,bb=m.advance(t+first,y,p,target,ff,zr,sensed,cfg,.01-first,dt,nq)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:y,aa=m.advance(t,y,p,target,ff,zr,sensed,cfg,.01,dt,nq)
        metrics[:3]+=aa[:3];metrics[4:7]+=aa[4:7];metrics[3]=max(metrics[3],aa[3]);metrics[7]=max(metrics[7],aa[7])
        metrics[8]=min(metrics[8],y[8]-m.E_MIN);metrics[9]=max(metrics[9],max(np.max(m.GMIN-y[4:8]),np.max(y[4:8]-m.GMAX)))
        trace[i,0]=t;trace[i,1:11]=oldy;trace[i,11:14]=force;trace[i,14:17]=sensed
        trace[i,17:21]=target;trace[i,21]=ff;trace[i,22:25]=resp[:3];trace[i,25]=temp
    duration=count*.01;metrics[1]=np.sqrt(metrics[1]/duration);metrics[2]=np.sqrt(metrics[2]/duration);metrics[4:7]/=duration
    return metrics,trace
