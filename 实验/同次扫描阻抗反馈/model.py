"""First-pass finite-footprint scan, tangential/normal dynamics and executed gains.

Exploratory foundation model, not calibrated tissue, full B0 or a prior replication.
Positive normal coordinate points into tissue. End-effector rotation is fixed.
"""
import numpy as np
from numba import njit

FD = 3.
DT = .001
CONTROL_DT = .01
DECISION_DT = .04
SENSOR_DELAY = .02
A, L, R = .010, .020, .050
MASS = np.array([.8, .5])
GMIN = np.array([40., 0., 2., 2.])
GMAX = np.array([1600., 800., 150., 150.])
GRATE = np.array([1200., 1200., 120., 120.])
E_MIN, E_MAX, CHI, ALPHA = 1e-5, .03, .5, 30.
# x=[tangent position error,z,tangent velocity error,vz,Kt,Kn,Dt,Dn,E,force integral]
# p=[kappa,contrast,boundary,width,h_amp,wavelength,phase,breath_amp,freq,
#    speed,mu,mucontrast,viscous_time,tangent_drag]


def sample(n, seed, family=None):
    rng=np.random.default_rng(seed)
    out=[]
    for i in range(n):
        f=i%3 if family is None else family
        out.append([np.exp(rng.uniform(np.log(1.5e6),np.log(9e6))),
            rng.uniform(-.6,.6) if f else 0.,rng.uniform(.025,.041),
            rng.uniform(.003,.010),rng.uniform(.0002,.0008),rng.uniform(.045,.10),
            rng.uniform(-np.pi,np.pi),rng.uniform(.0002,.0008) if f==2 else 0.,
            rng.uniform(.15,.35),rng.choice(np.array([.008,.012,.016])),
            rng.uniform(.12,.4),rng.uniform(-.25,.25) if f else 0.,
            rng.uniform(.002,.015),rng.uniform(1.,6.)])
    return np.array(out)


@njit(cache=True)
def ground(s,t,p):
    k=p[0]*np.exp(p[1]*np.tanh((s-p[2])/p[3]))
    h=p[4]*(np.sin(2*np.pi*s/p[5]+p[6])-np.sin(p[6]))
    h+=p[7]*(np.sin(2*np.pi*p[8]*t+p[6])-np.sin(p[6]))
    hd=p[7]*2*np.pi*p[8]*np.cos(2*np.pi*p[8]*t+p[6])
    return k,h,hd


@njit(cache=True)
def wrench(t,y,p,nq=33):
    s=p[9]*t+y[0]
    fn,moment,kn=0.,0.,0.
    for j in range(nq):
        u=-A+(j+.5)*2*A/nq
        k,h,hd=ground(s+u,t,p)
        delta=y[1]-u*u/(2*R)-h
        if delta>0:
            elastic=k*delta
            normal=max(0.,elastic+k*p[12]*(y[3]-hd))
            w=L*2*A/nq
            fn+=normal*w
            moment+=normal*u*w
            kn+=k*w
    mu=p[10]*np.exp(p[11]*np.tanh((s-p[2])/p[3]))
    vx=p[9]+y[2]
    ft=mu*fn*np.tanh(vx/.001)+p[13]*vx
    # Moment referred to contact point; raw sensor My also contains lever*Ft.
    rawmy=-moment+.025*ft
    return np.array([fn,ft,rawmy]),kn


@njit(cache=True)
def initial(p,g):
    y=np.zeros(10)
    lo,hi=-.005,.04
    for _ in range(60):
        y[1]=(lo+hi)/2
        # Start at stationary tissue, but tangential scan is already established.
        value=wrench(0.,y,p)[0][0]
        if value<FD:lo=y[1]
        else:hi=y[1]
    y[1]=(lo+hi)/2
    y[4:8]=g
    y[8]=.01
    return y


@njit(cache=True)
def gain_rate(y,target,zref):
    raw=(target-y[4:8])/.035
    lo=np.maximum(-GRATE,-ALPHA*(y[4:8]-GMIN))
    hi=np.minimum(GRATE,ALPHA*(GMAX-y[4:8]))
    u=np.minimum(np.maximum(raw,lo),hi)
    a=np.array([.5*y[0]**2,.5*(y[1]-zref)**2,0.,0.])
    pd=y[6]*y[2]**2+y[7]*y[3]**2+.5*((y[2])**2+y[3]**2)
    budget=CHI*pd+ALPHA*(y[8]-E_MIN)
    active=np.dot(a,u)>budget
    if active:
        left,right=0.,1.
        for _ in range(48):
            if np.dot(a,np.minimum(np.maximum(raw-right*a,lo),hi))<=budget:break
            right*=2.
        for _ in range(28):
            mid=(left+right)/2
            v=np.minimum(np.maximum(raw-mid*a,lo),hi)
            if np.dot(a,v)>budget:left=mid
            else:right=mid
        u=np.minimum(np.maximum(raw-right*a,lo),hi)
    edot=CHI*pd-np.dot(a,u)
    if y[8]>=E_MAX and edot>0:edot=0.
    return u,edot,active


@njit(cache=True)
def dynamics(y,force,target,ff,zref,measured,cfg):
    rate,edot,active=gain_rate(y,target,zref)
    err=FD-measured[0]
    # cfg=[Kt,Kn,Dt,Dn,kp_force,ki_force,force_observer_smoothing]
    ut=measured[1]-y[4]*y[0]-y[6]*y[2]+.5*cfg[7]
    un=FD+cfg[4]*err+y[9]-y[5]*(y[1]-zref)-y[7]*y[3]+ff
    uta=min(12.,max(-12.,ut));una=min(12.,max(-12.,un))
    dy=np.zeros(10)
    dy[0]=y[2];dy[1]=y[3]
    dy[2]=(uta-force[1]-.5*(y[2]+cfg[7]))/MASS[0]
    dy[3]=(una-force[0]-.5*y[3])/MASS[1]
    dy[4:8]=rate;dy[8]=edot
    iy=cfg[5]*err
    if (y[9]>=5. and iy>0) or (y[9]<=-5. and iy<0):iy=0.
    if (un>12. and iy>0) or (un< -12. and iy<0):iy=0.
    dy[9]=iy
    return dy,uta,una,active,max(abs(ut),abs(un))>12.


@njit(cache=True)
def stage(y,force,ut,un):
    # Force, tangential path, normal velocity and executed command effort.
    return ((force[0]-FD)/.3)**2+(y[0]/.0015)**2+.04*(y[3]/.03)**2+.015*((ut/3.)**2+((un-FD)/3.)**2)


@njit(cache=True)
def advance(t,y,p,target,ff,zref,measured,cfg,duration,dt=DT,nq=33):
    yy=y.copy();metrics=np.zeros(8)
    count=int(np.ceil(duration/dt-1e-10))
    for i in range(count):
        h=min(dt,duration-i*dt)
        f,_=wrench(t+i*dt,yy,p,nq)
        d,ut,un,active,sat=dynamics(yy,f,target,ff,zref,measured,cfg)
        mid=yy+.5*h*d
        fm,_=wrench(t+i*dt+.5*h,mid,p,nq)
        dm,utm,unm,am,sm=dynamics(mid,fm,target,ff,zref,measured,cfg)
        yy+=h*dm
        metrics[0]+=h*stage(mid,fm,utm,unm)
        metrics[1]+=h*(fm[0]-FD)**2
        metrics[2]+=h*mid[0]**2
        metrics[3]=max(metrics[3],fm[0])
        metrics[4]+=h*(fm[0]<.2)
        metrics[5]+=h*am
        metrics[6]+=h*sm
        metrics[7]=max(metrics[7],abs(mid[0]))
    return yy,metrics


@njit(cache=True)
def response(history,i,old):
    # Aligned delayed force and robot history. Ordinary regularized local regression.
    # history columns: ex,z,ev,vz,Fn,Ft,My,x_absolute
    # response=[k,g_s,unmodeled_time_rate,Fn,Ft,z,x,vz]
    lag=2;j=max(0,i-lag)
    obs=history[j]
    mnormal=-obs[6]+.025*obs[5]
    gs=3*mnormal/(A*A)
    gs=min(1500.,max(-1500.,gs))
    gs=.65*old[1]+.35*gs
    # Differential affine regression over 0.35s; drift is a nuisance, not identified motion.
    k=old[0];drift=old[2]
    n=min(j,35)
    if n>=8:
        a,b,c,d,e=1e-7,0.,.02,1e-7*k,.02*drift
        for q in range(1,n+1):
            prev=history[j-q]
            dz=obs[1]-prev[1];elapsed=q*CONTROL_DT
            mm=-prev[6]+.025*prev[5]
            gg=.5*(gs+3*mm/(A*A))
            df=obs[4]-prev[4]-gg*(obs[7]-prev[7])
            a+=dz*dz;b+=dz*elapsed;c+=elapsed*elapsed
            d+=dz*df;e+=elapsed*df
        den=a*c-b*b
        nk=(d*c-b*e)/max(den,1e-15)
        nd=(a*e-b*d)/max(den,1e-15)
        k=.9*k+.1*min(7000.,max(200.,nk))
        drift=.8*drift+.2*min(8.,max(-8.,nd))
    # Separate purely temporal predictor, to compete rather than freeze force.
    n2=min(j,12)
    temporal=drift+gs*(history[i,7]-history[max(0,i-1),7])/CONTROL_DT
    if n2>=4:
        tm=n2/2;den=0.;num=0.
        for q in range(n2+1):
            tau=(q-tm)*CONTROL_DT
            prev=history[j-n2+q]
            residual=prev[4]-k*prev[1]
            num+=tau*residual;den+=tau*tau
        temporal=num/max(den,1e-12)
    elapsed=(i-j)*CONTROL_DT
    fn=obs[4]+k*(history[i,1]-obs[1])+gs*(history[i,7]-obs[7])+drift*elapsed
    ft=obs[5]*max(.2,min(3.,fn/max(obs[4],.2)))
    current=np.array([k,gs,drift,fn,ft,history[i,1],history[i,7],history[i,3]])
    return current,temporal,j


@njit(cache=True)
def library(y,zref,cfg,resp,method,temporal):
    # Fixed candidate layout across contexts; canonical gain/ref shapes, not oracle gains.
    # methods 0=fixed delayed,1=observer fixed,2=critical,3=analytic transport,
    # 4=temporal joint MPC,5=spatial gain MPC,6=spatial joint MPC,7=spatial ff only
    arr=np.zeros((25,5))
    arr[0,:4]=cfg[:4]
    k=resp[0]
    for index in range(1,25):
        j=index-1
        kfac=.5 if (j//12)==0 else 1.
        kn=0. if (j//6)%2==0 else 150.
        zeta=np.array([.35,.7,1.1])[(j//2)%3]
        fraction=0. if j%2==0 else 1.
        kt=cfg[0]*kfac
        dn=2*zeta*np.sqrt(MASS[1]*((1+cfg[4])*k+kn))
        dt=2*.8*np.sqrt(MASS[0]*kt)
        arr[index,:4]=np.array([kt,kn,dt,dn])
        source=temporal if method==4 else resp[1]*(cfg[7]+y[2])+resp[2]
        vel=-source/max(k,200.)
        arr[index,4]=fraction*dn*min(.02,max(-.02,vel))
        if method==5:arr[index,4]=0.
        if method==7:
            arr[index,:4]=cfg[:4]
            arr[index,4]=fraction*cfg[3]*min(.02,max(-.02,vel))
    arr[:,:4]=np.minimum(np.maximum(arr[:,:4],GMIN),GMAX)
    return arr


@njit(cache=True)
def predicted_force(tau,y,start,resp,temporal,mode):
    dx=y[0]-start[0] # both ex; add known moving reference displacement
    source=temporal*tau if mode==4 else resp[1]*(dx+tau*resp[8])+resp[2]*tau
    fn=resp[3]+resp[0]*(y[1]-start[1])+source+.008*resp[0]*(y[3]-start[3])
    fn=max(0.,fn)
    ft=resp[4]*fn/max(.2,resp[3])
    return np.array([fn,ft,0.])


@njit(cache=True)
def score_library(y,zref,cfg,resp,temporal,mode,arr):
    values=np.zeros(len(arr))
    # A deterministic local response model with the same PI and physical gain executor.
    # Not called a robust certificate; motion/parameter uncertainty not covered here.
    dt=.004;n=48
    ext=np.zeros(9);ext[:8]=resp;ext[8]=cfg[7]
    for c in range(len(arr)):
        yy=y.copy();tar=arr[c,:4];ff=arr[c,4]
        sensed=np.array([resp[3],resp[4],0.]);past=np.zeros((n+1,3))
        for step in range(n):
            tau=step*dt
            f=predicted_force(tau,yy,y,ext,temporal,mode)
            past[step]=f
            lagged=past[max(0,step-5)]
            # Match sensor filtering approximately at the predictor integration rate.
            alpha=1-(1-cfg[6])**(dt/CONTROL_DT)
            sensed=alpha*lagged+(1-alpha)*sensed
            d,ut,un,_,_=dynamics(yy,f,tar,ff,zref,sensed,cfg)
            mid=yy+.5*dt*d
            fm=predicted_force(tau+.5*dt,mid,y,ext,temporal,mode)
            dm,um,nm,_,_=dynamics(mid,fm,tar,ff,zref,sensed,cfg)
            yy+=dt*dm
            value=stage(mid,fm,um,nm)
            if not np.isfinite(value) or abs(yy[1])>.08:
                values[c]=1e10;break
            values[c]+=dt*value
    return values


@njit(cache=True)
def features(y,zref,cfg,resp,temporal):
    return np.array([y[0]/.002,(y[1]-zref)/.004,y[2]/.02,y[3]/.02,
        y[4]/800.,y[5]/500.,y[6]/60.,y[7]/60.,y[8]/.01,y[9]/3.,
        resp[0]/2000.,resp[1]/500.,resp[2]/3.,(resp[3]-FD)/.5,
        resp[4]/2.,temporal/3.])


@njit(cache=True)
def decide(y,zref,cfg,resp,temporal,method):
    if method<2:
        out=np.zeros(5);out[:4]=cfg[:4]
        return out,np.zeros(25)
    if method in (2,3):
        out=np.zeros(5);out[:4]=cfg[:4]
        out[3]=min(150.,max(2.,2*.8*np.sqrt(MASS[1]*((1+cfg[4])*resp[0]+cfg[1]))))
        if method==3:
            vel=-(resp[1]*(cfg[7]+y[2])+resp[2])/resp[0]
            out[4]=out[3]*min(.02,max(-.02,vel))
        return out,np.zeros(25)
    candidates=library(y,zref,cfg,resp,method,temporal)
    scores=score_library(y,zref,cfg,resp,temporal,method,candidates)
    return candidates[np.argmin(scores)],scores


@njit(cache=True)
def episode(p,noise,method,cfg0,latency=0.,dt=DT,nq=33,collect=False):
    cfg=cfg0.copy();cfg[7]=p[9]
    y=initial(p,cfg[:4]);zref=y[1]
    count=len(noise)-1
    hist=np.zeros((count+1,8))
    trace=np.zeros((count,26))
    ndecision=(count+3)//4
    xx=np.zeros((ndecision,16));qq=np.zeros((ndecision,25))
    metrics=np.zeros(10);metrics[8]=1e3;metrics[9]=0.
    response0=np.array([1000.,0.,0.,FD,.6,zref,0.,0.])
    measured=np.array([FD,.6,0.]);target=cfg[:4].copy();ff=0.
    pending=target.copy();pendingff=0.;release=-1.;decisions=0
    for i in range(count):
        t=i*CONTROL_DT
        oldy=y.copy()
        force,_=wrench(t,y,p,nq)
        hist[i,:4]=y[:4]
        hist[i,4:7]=force+noise[i]
        hist[i,7]=p[9]*t+y[0]
        resp,temp,j=response(hist,i,response0)
        response0=resp.copy()
        ob=hist[max(0,i-2),4:7]
        mm=ob.copy()
        if method!=0:mm[0:2]=resp[3:5]
        measured=cfg[6]*mm+(1-cfg[6])*measured
        resp[3:5]=measured[:2]
        if i%4==0:
            decision,scores=decide(y,zref,cfg,resp,temp,method)
            pending=decision[:4].copy();pendingff=decision[4]
            release=t+latency
            if collect:
                xx[decisions]=features(y,zref,cfg,resp,temp);qq[decisions]=scores
            decisions+=1
        if release>=0 and release<=t+1e-12:
            target=pending.copy();ff=pendingff;release=-1.
        if release>t and release<t+CONTROL_DT:
            first=release-t
            y,aa=advance(t,y,p,target,ff,zref,measured,cfg,first,dt,nq)
            target=pending.copy();ff=pendingff;release=-1.
            y,bb=advance(t+first,y,p,target,ff,zref,measured,cfg,CONTROL_DT-first,dt,nq)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:
            y,aa=advance(t,y,p,target,ff,zref,measured,cfg,CONTROL_DT,dt,nq)
        metrics[:3]+=aa[:3];metrics[4:7]+=aa[4:7]
        metrics[3]=max(metrics[3],aa[3]);metrics[7]=max(metrics[7],aa[7])
        metrics[8]=min(metrics[8],y[8]-E_MIN)
        metrics[9]=max(metrics[9],max(np.max(GMIN-y[4:8]),np.max(y[4:8]-GMAX)))
        trace[i,0]=t
        trace[i,1:11]=oldy
        trace[i,11:14]=force
        trace[i,14:17]=measured
        trace[i,17:21]=target
        trace[i,21]=ff
        trace[i,22:25]=resp[:3]
        trace[i,25]=temp
    duration=count*CONTROL_DT
    metrics[1]=np.sqrt(metrics[1]/duration)
    metrics[2]=np.sqrt(metrics[2]/duration)
    metrics[4:7]/=duration
    return metrics,trace,xx[:decisions],qq[:decisions]
