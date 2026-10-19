"""Policy-space rollout with causal observer prediction and local response models.

Transfers classical rollout/policy-space MPC; this implementation is exploratory,
not a new optimality theorem or a calibrated Bayesian contact posterior.
"""
from pathlib import Path
import sys
import numpy as np
from numba import njit

BASE = Path(__file__).resolve().parents[1] / '同次扫描阻抗反馈'
sys.path.insert(0, str(BASE))
import model as m
import coupled_contact as cm


@njit(cache=True)
def response_models(hist, i, resp):
    """Nine local hypotheses fitted ONLY to aligned, already-arrived wrenches.

    Returns k,c,spatial slope,time drift,Fn-now,Ft-now. Difference residuals
    share sensor noise, so normalized scores are weights, NOT posterior odds.
    """
    j=max(0,i-2);obs=hist[j];start=max(0,j-35)
    models=np.zeros((9,6));rss=np.zeros(9)
    gs=resp[1]
    for a in range(9):
        k=min(7000.,max(200.,resp[0]*np.array([.5,1.,2.])[a//3]))
        c=k*np.array([.002,.008,.015])[a%3]
        numerator=0.;denominator=1e-8
        for q in range(start,j):
            lag=(j-q)*.01
            residual=obs[4]-hist[q,4]-k*(obs[1]-hist[q,1])-c*(obs[3]-hist[q,3])-gs*(obs[7]-hist[q,7])
            numerator+=lag*residual;denominator+=lag*lag
        drift=min(8.,max(-8.,numerator/denominator)) if j>start else resp[2]
        error=0.
        for q in range(start,j):
            predicted=k*(obs[1]-hist[q,1])+c*(obs[3]-hist[q,3])+gs*(obs[7]-hist[q,7])+drift*(j-q)*.01
            error+=(obs[4]-hist[q,4]-predicted)**2
        rss[a]=error/max(1,j-start)
        fn=obs[4]+k*(hist[i,1]-obs[1])+c*(hist[i,3]-obs[3])+gs*(hist[i,7]-obs[7])+drift*(i-j)*.01
        models[a]=np.array([k,c,gs,drift,max(.0,fn),obs[5]*max(.2,min(3.,fn/max(obs[4],.2)))])
    # 0.03N declared sensor noise; no claim of independent difference residuals.
    weights=np.exp(-.5*(rss-rss.min())/(2*.03**2))
    weights/=weights.sum()
    return models,weights,rss


@njit(cache=True)
def force(tau,y,initial,par,speed):
    fn=max(0.,par[4]+par[0]*(y[1]-initial[1])+par[1]*(y[3]-initial[3])
           +par[2]*(speed*tau+y[0]-initial[0])+par[3]*tau)
    ft=par[5]*fn/max(.2,par[4])
    return np.array([fn,ft,-par[2]*m.A*m.A/3+.025*ft])


@njit(cache=True)
def action(y,zref,cfg,resp,temp,theta,scales):
    out,_=m.decide(y,zref,cfg,resp,temp,12,theta)
    out[1]=min(800.,max(0.,out[1]*scales[0]))
    out[3]=min(150.,max(2.,out[3]*scales[1]))
    out[4]*=scales[1]
    return out


@njit(cache=True)
def rollout(y,zref,cfg,hist,i,resp,sensed,temporal,par,theta,scales,horizon,consistent,release_delay,held):
    """Same future policy in every model; future actions use its own observations.

    Gain decisions 25Hz, observations 100Hz, plant integration 500Hz. Future
    sensor noise has zero conditional mean; uncertainty here is model variation.
    """
    yy=y.copy();sense=sensed.copy();old=resp.copy()
    offset=max(0,i-39);nold=i-offset+1;nobs=int(round(horizon/.01))
    history=np.zeros((nold+nobs+1,8));history[:nold]=hist[offset:i+1]
    # The simulator stores undelivered samples too. Replace them by model
    # predictions before the forecast can use them as future-arriving data.
    for q in range(nold):
        if offset+q>max(0,i-2):
            pose=y.copy();pose[:4]=history[q,:4]
            history[q,4:7]=force((offset+q-i)*.01,pose,y,par,cfg[7])
    target=np.zeros(5);total=0.;past=np.zeros((nobs+1,3))
    for step in range(nobs):
        tau=step*.01
        ff=force(tau,yy,y,par,cfg[7]);past[step]=ff
        if step>0:
            h=nold-1+step
            history[h,:4]=yy[:4];history[h,4:7]=ff;history[h,7]=hist[i,7]+cfg[7]*tau+yy[0]-y[0]
            rr,temp,_=m.response(history,h,old);old=rr.copy()
            ob=rr[3:5] if consistent else past[max(0,step-2),:2]
            sense[:2]=cfg[6]*ob+(1-cfg[6])*sense[:2]
            rr[3:5]=sense[:2]
        else:
            rr=resp.copy();temp=temporal
        if step%4==0:target=action(yy,zref,cfg,rr,temp,theta,scales)
        for sub in range(5):
            t=tau+sub*.002
            executing=held if t<release_delay else target
            f=force(t,yy,y,par,cfg[7]);d,_,_,_,_=m.dynamics(yy,f,executing[:4],executing[4],zref,sense,cfg)
            mid=yy+.001*d;fm=force(t+.001,mid,y,par,cfg[7])
            dm,ut,un,_,_=m.dynamics(mid,fm,executing[:4],executing[4],zref,sense,cfg)
            yy+=.002*dm
            if not np.isfinite(yy).all() or abs(yy[1])>.08:return 1e8
            total+=.002*m.stage(mid,fm,ut,un,zref,cfg[8])
    return total


@njit(cache=True)
def decide(y,zref,cfg,hist,i,resp,sensed,temporal,bank,scales,horizon,ensemble,consistent,release_delay,held):
    models,weights,rss=response_models(hist,i,resp)
    if not ensemble:
        # Nominal local response identical in meaning to the previous prototype.
        models[0]=np.array([resp[0],.008*resp[0],resp[1],resp[2],resp[3],resp[4]])
        weights[:]=0.;weights[0]=1.
    values=np.zeros(len(bank))
    for a in range(len(bank)):
        for b in range(len(models)):
            if weights[b]>1e-4:
                values[a]+=weights[b]*rollout(y,zref,cfg,hist,i,resp,sensed,temporal,models[b],bank[a],scales[a],horizon,consistent,release_delay,held)
    best=np.argmin(values)
    # Current branch uses deployment observations, not model-specific hidden state.
    return action(y,zref,cfg,resp,temporal,bank[best],scales[best]),best,values


@njit(cache=True)
def episode(p,noise,cfg0,bank,scales,horizon=.2,ensemble=True,consistent=True,
            latency=0.,policy_index=-1,dt=.0005,nq=129):
    cfg=cfg0.copy();cfg[7]=p[9];y=cm.initial(p,cfg[:4],nq);zr=y[1]
    n=len(noise)-1;hist=np.zeros((n+1,8));trace=np.zeros((n,27))
    kk=1.5*m.FD/zr if zr<m.A*m.A/(2*m.R) else m.FD/(zr-m.A*m.A/(6*m.R))
    old=np.array([min(7000.,max(200.,kk)),0.,0.,3.,.6,zr,0.,0.])
    sensed=np.array([3.,.6,0.]);target=cfg[:4].copy();ff=0.
    pending=target.copy();pendingff=0.;release=-1.;metrics=np.zeros(10);metrics[8]=1e3;selected=-1
    for i in range(n):
        t=i*.01;before=y.copy();f=cm.contact(t,y,p,nq)[0]
        hist[i,:4]=y[:4];hist[i,4:7]=f+noise[i];hist[i,7]=p[9]*t+y[0]
        resp,temp,_=m.response(hist,i,old);old=resp.copy()
        sensed[:2]=cfg[6]*resp[3:5]+(1-cfg[6])*sensed[:2];resp[3:5]=sensed[:2]
        if release>=0 and release<=t+1e-12:target=pending.copy();ff=pendingff;release=-1.
        if i%4==0 and release<0:
            if policy_index<0:
                held=np.empty(5);held[:4]=target;held[4]=ff
                out,selected,_=decide(y,zr,cfg,hist,i,resp,sensed,temp,bank,scales,horizon,ensemble,consistent,latency,held)
            else:
                selected=policy_index;out=action(y,zr,cfg,resp,temp,bank[selected],scales[selected])
            pending=out[:4].copy();pendingff=out[4];release=t+latency
        if release>=0 and release<=t+1e-12:target=pending.copy();ff=pendingff;release=-1.
        if release>t and release<t+.01:
            first=release-t;y,aa=cm.advance(t,y,p,target,ff,zr,sensed,cfg,first,dt,nq)
            target=pending.copy();ff=pendingff;release=-1.
            y,bb=cm.advance(t+first,y,p,target,ff,zr,sensed,cfg,.01-first,dt,nq)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:y,aa=cm.advance(t,y,p,target,ff,zr,sensed,cfg,.01,dt,nq)
        metrics[:3]+=aa[:3];metrics[4:7]+=aa[4:7];metrics[3]=max(metrics[3],aa[3]);metrics[7]=max(metrics[7],aa[7])
        metrics[8]=min(metrics[8],y[8]-m.E_MIN);metrics[9]=max(metrics[9],max(np.max(m.GMIN-y[4:8]),np.max(y[4:8]-m.GMAX)))
        trace[i,0]=t;trace[i,1:11]=before;trace[i,11:14]=f;trace[i,14:17]=sensed
        trace[i,17:21]=target;trace[i,21]=ff;trace[i,22:25]=resp[:3];trace[i,25]=temp;trace[i,26]=selected
    metrics[1]=np.sqrt(metrics[1]/(n*.01));metrics[2]=np.sqrt(metrics[2]/(n*.01));metrics[4:7]/=n*.01
    return metrics,trace
