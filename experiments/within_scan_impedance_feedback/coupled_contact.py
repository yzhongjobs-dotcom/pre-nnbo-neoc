"""Potential/Rayleigh-consistent translational contact sensitivity model.

Same finite footprint, spatial fields and actuator interface as model.py. Adds
material-point convective strain rate and conservative/dissipative tangential
coupling. Moment remains an idealized pressure-moment/lever observation channel;
not a full 6D probe/contact stress model or a calibrated tissue constitutive law.
"""
import numpy as np
from numba import njit
import model as m


@njit(cache=True)
def contact(t,y,p,nq=129):
    s=p[9]*t+y[0];vx=p[9]+y[2];vz=y[3]
    fn,moment,fel,fdamp,potential,diss,bodywork,kn=0.,0.,0.,0.,0.,0.,0.,0.
    for j in range(nq):
        u=-m.A+(j+.5)*2*m.A/nq
        x=s+u;k,h,hd=m.ground(x,t,p)
        penetration=y[1]-u*u/(2*m.R)-h
        if penetration>0:
            shape=np.tanh((x-p[2])/p[3])
            kx=k*p[1]*(1-shape*shape)/p[3]
            hx=p[4]*2*np.pi/p[5]*np.cos(2*np.pi*x/p[5]+p[6])
            rate=vz+u*vx/m.R-hd
            elastic=k*penetration
            pressure=max(0.,elastic+k*p[12]*rate)
            w=m.L*2*m.A/nq
            fn+=pressure*w;moment+=pressure*u*w
            # Exact derivative of the discrete moving-window elastic potential.
            # This includes finite-face elastic boundary flux, not just M/R.
            fel+=(.5*kx*penetration**2-elastic*hx)*w
            fdamp+=(pressure-elastic)*u/m.R*w
            potential+=.5*k*penetration**2*w
            diss+=(pressure-elastic)*rate*w
            bodywork+=pressure*hd*w
            if pressure>0:kn+=k*w
    mu=p[10]*np.exp(p[11]*np.tanh((s-p[2])/p[3]))
    friction=mu*fn*np.tanh(vx/.001)+p[13]*vx
    ft=fel+fdamp+friction
    diss+=friction*vx
    rawmy=-moment+.025*ft
    return np.array([fn,ft,rawmy]),potential,diss,bodywork,kn


@njit(cache=True)
def initial(p,g,nq=129):
    y=np.zeros(10);lo,hi=-.005,.04
    for _ in range(60):
        y[1]=(lo+hi)/2
        if contact(0.,y,p,nq)[0][0]<m.FD:lo=y[1]
        else:hi=y[1]
    y[1]=(lo+hi)/2;y[4:8]=g;y[8]=.01
    return y


@njit(cache=True)
def advance(t,y,p,target,ff,zr,sensed,cfg,duration,dt=.0005,nq=129):
    yy=y.copy();metrics=np.zeros(8)
    for i in range(int(np.ceil(duration/dt-1e-10))):
        h=min(dt,duration-i*dt)
        force=contact(t+i*dt,yy,p,nq)[0]
        d,ua,ub,active0,sat0=m.dynamics(yy,force,target,ff,zr,sensed,cfg)
        mid=yy+.5*h*d;fm=contact(t+i*dt+.5*h,mid,p,nq)[0]
        dm,ut,un,active,sat=m.dynamics(mid,fm,target,ff,zr,sensed,cfg)
        yy+=h*dm
        metrics[0]+=h*m.stage(mid,fm,ut,un,zr,cfg[8])
        metrics[1]+=h*(fm[0]-3.)**2;metrics[2]+=h*mid[0]**2
        metrics[3]=max(metrics[3],fm[0]);metrics[4]+=h*(fm[0]<.2)
        metrics[5]+=h*active;metrics[6]+=h*sat;metrics[7]=max(metrics[7],abs(mid[0]))
    return yy,metrics


@njit(cache=True)
def episode(p,noise,method,cfg0,policy=np.zeros(8),latency=0.,dt=.0005,nq=129):
    cfg=cfg0.copy();cfg[7]=p[9];y=initial(p,cfg[:4],nq);zr=y[1]
    n=len(noise)-1;hist=np.zeros((n+1,8));trace=np.zeros((n,26))
    kk=1.5*m.FD/zr if zr<m.A*m.A/(2*m.R) else m.FD/(zr-m.A*m.A/(6*m.R))
    old=np.array([min(7000.,max(200.,kk)),0.,0.,3.,.6,zr,0.,0.])
    sensed=np.array([3.,.6,0.]);target=cfg[:4].copy();ff=0.
    pending=target.copy();pendingff=0.;release=-1.;metrics=np.zeros(10);metrics[8]=1e3
    for i in range(n):
        t=i*.01;before=y.copy();force=contact(t,y,p,nq)[0]
        hist[i,:4]=y[:4];hist[i,4:7]=force+noise[i];hist[i,7]=p[9]*t+y[0]
        resp,temp,_=m.response(hist,i,old);old=resp.copy()
        ob=hist[max(0,i-2),4:7].copy()
        if method!=0:ob[:2]=resp[3:5]
        sensed=cfg[6]*ob+(1-cfg[6])*sensed;resp[3:5]=sensed[:2]
        if i%4==0:
            decision,_=m.decide(y,zr,cfg,resp,temp,method,policy)
            pending=decision[:4].copy();pendingff=decision[4];release=t+latency
        if release>=0 and release<=t+1e-12:target=pending.copy();ff=pendingff;release=-1.
        if release>t and release<t+.01:
            first=release-t;y,aa=advance(t,y,p,target,ff,zr,sensed,cfg,first,dt,nq)
            target=pending.copy();ff=pendingff;release=-1.
            y,bb=advance(t+first,y,p,target,ff,zr,sensed,cfg,.01-first,dt,nq)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:y,aa=advance(t,y,p,target,ff,zr,sensed,cfg,.01,dt,nq)
        metrics[:3]+=aa[:3];metrics[4:7]+=aa[4:7];metrics[3]=max(metrics[3],aa[3]);metrics[7]=max(metrics[7],aa[7])
        metrics[8]=min(metrics[8],y[8]-m.E_MIN);metrics[9]=max(metrics[9],max(np.max(m.GMIN-y[4:8]),np.max(y[4:8]-m.GMAX)))
        trace[i,0]=t;trace[i,1:11]=before;trace[i,11:14]=force;trace[i,14:17]=sensed
        trace[i,17:21]=target;trace[i,21]=ff;trace[i,22:25]=resp[:3];trace[i,25]=temp
    metrics[1]=np.sqrt(metrics[1]/(n*.01));metrics[2]=np.sqrt(metrics[2]/(n*.01));metrics[4:7]/=n*.01
    return metrics,trace
