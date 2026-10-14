"""Exploratory failure frontier, restoring R15 gain preference and constraints.

Shared actual gain dynamics, actuator clipping and continuous tank/rate protection.
Direct shooting is a diagnostic for the required decision horizon, not R15 HJB.
"""
from pathlib import Path
import sys,json,time,argparse
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'第一步_扫描优化必要性_2026-10-13'))
import scan_validation as old

BASE=old.CFG.copy();GP=np.array([250.,80.,8.,6.]);LOW=old.LOW;HIGH=old.HIGH;RATE=old.RATE;SPAN=old.SPAN
CASES=dict(loose=dict(coupling=0.,friction_scale=1.,reserve=.05,actuator=[.15,.3]),
           coupled=dict(coupling=.65,friction_scale=1.,reserve=.05,actuator=[.15,.3]),
           budget=dict(coupling=.65,friction_scale=1.,reserve=.000005,actuator=[.15,.3]),
           friction=dict(coupling=.65,friction_scale=1.6,reserve=.000005,actuator=[.12,.3]))
CASE=None;MI=None;ACT=None
MATERIAL=old.material

def material(x,kind):
    E,mu,kd=MATERIAL(x,kind)
    return E,mu*CASE['friction_scale'],kd

def rates(y,e,ev,target):
    g=y[:,4:8]
    lo=np.maximum(-RATE,-25*(g-LOW));hi=np.minimum(RATE,25*(HIGH-g))
    u=np.clip((target-g)/.04,lo,hi)
    p=.5*e**2
    diss=np.sum(g[:,2:]*ev**2,axis=1)+BASE['robot_damping']*np.sum(y[:,2:4]**2,axis=1)
    cap=.3*diss+10*(y[:,8]-BASE['energy_min'])
    excess=np.sum(p*u[:,:2],axis=1)-cap
    active=excess>0
    if np.any(active):
        ua=u[active,:2];pa=p[active];lower=lo[active,:2]
        lam=excess[active]/np.maximum(np.sum(pa**2*RATE[:2]**2,axis=1),1e-30)
        proposal=ua-lam[:,None]*pa*RATE[:2]**2
        # Exact two-variable projection: at most one lower bound is active
        # when the tank cap is nonnegative and zero gain rate is feasible.
        hit0=proposal[:,0]<lower[:,0]
        proposal[hit0,0]=lower[hit0,0]
        proposal[hit0,1]=(cap[active][hit0]-pa[hit0,0]*lower[hit0,0])/np.maximum(pa[hit0,1],1e-30)
        hit1=proposal[:,1]<lower[:,1]
        proposal[hit1,1]=lower[hit1,1]
        proposal[hit1,0]=(cap[active][hit1]-pa[hit1,1]*lower[hit1,1])/np.maximum(pa[hit1,0],1e-30)
        u[active,:2]=np.clip(proposal,lower,hi[active,:2])
    return u,diss

def dynamics(y,t,speed,target,env,kind,feedback=0.):
    fn,fr,e,ev,n,tan,delta,vt=old.geometry(y,t,speed,env,kind)
    u,diss=rates(y,e,ev,target);g=y[:,4:8]
    intended=(BASE['force_target']+feedback*(BASE['force_target']-fn))[:,None]*n-g[:,:2]*e-g[:,2:]*ev
    command=np.clip(intended,-ACT,ACT)
    fc=-fn[:,None]*n-fr[:,None]*tan
    acc=(command+fc-BASE['robot_damping']*y[:,2:4])@MI.T
    energy=.3*diss-.5*np.sum(e**2*u[:,:2],axis=1)
    rhs=np.column_stack((y[:,2:4],acc,u,energy))
    ell=((fn-BASE['force_target'])/.02)**2+(e[:,0]/.0001)**2+.2*(e[:,1]/.00025)**2
    ell+=.02*np.sum((ev/.008)**2,axis=1)+.01*np.sum((command/.2)**2,axis=1)
    ell+=.1*np.sum(((g-GP)/GP)**2,axis=1)+.0001*np.sum((u/RATE)**2,axis=1)
    penalty=1e4*np.maximum(-delta/.0001,0)**2+1e4*np.maximum((fn-.23)/.02,0)**2
    penalty+=1e4*np.sum((np.maximum(abs(intended)-ACT,0)/ACT)**2,axis=1)
    penalty+=1e4*np.maximum((.001-y[:,8])/.001,0)**2
    return rhs,ell,penalty,(fn,e,delta,command,u,vt,intended)

def configure(name):
    global CASE,MI,ACT
    CASE=CASES[name];rho=CASE['coupling'];M=np.array([[.2,rho*np.sqrt(.2*.3)],[rho*np.sqrt(.2*.3),.3]])
    MI=np.linalg.inv(M);ACT=np.array(CASE['actuator'])
    old.CFG.update(BASE);old.CFG['initial_energy']=.001+CASE['reserve']
    old.material=material;old.dynamics=dynamics

def steady_target(y,t,speed,env,feedback=0.):
    # Strong cheap quasi-static optimization. No dependence on future material.
    xr=speed*t;h,s=old.surface(xr);r=np.sqrt(1+s*s)
    n=np.array([-s,1.])/r;tan=np.array([1.,s])/r
    test=old.initial(speed,GP)[None,:];test[0,:2]=[xr,h+old.DNOM]
    f0=old.geometry(test,t,speed,env,'ramp')[0][0]
    ep=1e-7;test[0,1]+=ep
    kn=(old.geometry(test,t,speed,env,'ramp')[0][0]-f0)/ep
    E,mu,kd=env;vt=speed*r;direction=n+mu*tan
    def value(q,output=False):
        k=LOW[:2]+q*SPAN[:2]
        A=np.diag(k)+kn*np.outer(direction+feedback*n,n)
        e=np.linalg.solve(A,(BASE['force_target']-f0)*(1+feedback)*n-mu*f0*tan-(kd*vt+BASE['robot_damping']*vt)*tan)
        fn=f0+kn*n@e
        g=np.r_[k,GP[2:]]
        cost=((fn-.14)/.02)**2+(e[0]/.0001)**2+.2*(e[1]/.00025)**2+.1*np.sum(((g-GP)/GP)**2)
        command=(.14+feedback*(.14-fn))*n-k*e
        cost+=1e4*np.sum((np.maximum(abs(command)-ACT,0)/ACT)**2)
        return g if output else cost
    res=minimize(value,(np.clip(y[4:6],LOW[:2],HIGH[:2])-LOW[:2])/SPAN[:2],method='L-BFGS-B',bounds=[(0,1)]*2,
                 options=dict(maxiter=30,ftol=1e-9))
    return value(res.x,True)

def qp_target(y,t,speed,env,h=.04):
    # A one-step quadratic planner retaining gain preference and tank limits.
    fn,fr,e,ev,n,tan,delta,vt=old.geometry(y[None,:],t,speed,env,'ramp')
    e=e[0];ev=ev[0];n=n[0];tan=tan[0];g=y[4:8]
    base_cmd=.14*n-g[:2]*e-g[2:]*ev
    cmdJ=-np.column_stack((np.diag(e),np.diag(ev)))
    fc=-fn[0]*n-fr[0]*tan;acc0=MI@(base_cmd+fc-BASE['robot_damping']*y[2:4])
    aJ=MI@cmdJ*.5
    x0=e+h*ev+.5*h*h*acc0;X=.5*h*h*aJ
    v0=ev+h*acc0;V=h*aJ
    ep=1e-7;yp=y.copy();yp[1]+=ep
    kn=(old.geometry(yp[None,:],t,speed,env,'ramp')[0][0]-fn[0])/ep
    f0=fn[0]-.14+kn*n@(.5*h*h*acc0+h*y[2:4]);F=kn*n@X
    mats=[X/np.array([.0001,.00025/np.sqrt(.2)])[:,None],np.sqrt(.02)*V/.008,
          F[None,:]/.02,np.sqrt(.1)*np.diag(1/GP)]
    vals=[x0/np.array([.0001,.00025/np.sqrt(.2)]),np.sqrt(.02)*v0/.008,
          np.array([f0/.02]),np.sqrt(.1)*(g-GP)/GP]
    A=np.vstack(mats)*SPAN[None,:];b=np.concatenate(vals)
    H=2*A.T@A+1e-10*np.eye(4);grad=2*A.T@b
    lo=np.maximum(LOW-g,-h*RATE)/SPAN;hi=np.minimum(HIGH-g,h*RATE)/SPAN
    power=.5*np.r_[e**2,[0.,0.]]*SPAN
    dissip=np.sum(g[2:]*ev**2)+BASE['robot_damping']*np.sum(y[2:4]**2)
    cap=h*(.3*dissip+10*(y[8]-.001))
    res=minimize(lambda q:(.5*q@H@q+grad@q,H@q+grad),np.zeros(4),jac=True,
        method='SLSQP',bounds=list(zip(lo,hi)),constraints=[dict(type='ineq',fun=lambda q:cap-power@q,jac=lambda q:-power)],
        options=dict(maxiter=30,ftol=1e-8))
    return np.clip(g+res.x*SPAN,LOW,HIGH)

def optimize(y,t,speed,env,q):
    candidates=[q,np.tile((GP-LOW)/SPAN,2),np.tile((steady_target(y,t,speed,env)-LOW)/SPAN,2)]
    vals=old.shooting(candidates,y,t,speed,env,'ramp');start=np.asarray(candidates)[np.argmin(vals)]
    res=minimize(old.objective,start,args=(y,t,speed,env,'ramp'),jac=True,method='L-BFGS-B',bounds=[(0,1)]*8,
                 options=dict(maxiter=35,ftol=1e-7,gtol=3e-5,maxls=20))
    if res.fun>np.min(vals):return start,dict(success=False,message='kept feasible incumbent',iterations=int(res.nit))
    return res.x,dict(success=bool(res.success),message=str(res.message),iterations=int(res.nit))

def run(speed,mode,param=None,dt=.002,save=None):
    duration=BASE['scan_length']/speed;steps=round(duration/dt);period=round(.1/dt)
    y=old.initial(speed,GP)[None,:];target=GP[None,:];q=np.tile((GP-LOW)/SPAN,2)
    feedback=(param or {}).get('feedback',0.) if mode=='force' else 0.
    cost=0.;rows=[];seconds=0.;success=[]
    for j in range(steps):
        t=j*dt
        if j%period==0:
            env=material(y[0,0],'ramp');tic=time.perf_counter()
            if mode=='fixed':target=np.array(param['gains'])[None,:]
            elif mode in ['steady','force']:target=steady_target(y[0],t,speed,env,feedback)[None,:]
            elif mode=='qp':target=qp_target(y[0],t,speed,env)[None,:]
            elif mode=='mpc':
                q,info=optimize(y[0],t,speed,env,q);target=(LOW+q[:4]*SPAN)[None,:];success.append(info)
            seconds+=time.perf_counter()-tic
        _,ell,penalty,out=dynamics(y,t,speed,target,None,'ramp',feedback)
        fn,e,delta,command,u,vt,intended=out
        rows.append([t,*y[0],fn[0],*e[0],delta[0],*command[0],*u[0],vt[0],ell[0],penalty[0],*intended[0]])
        y,l,p=old.rk4(y,t,dt,speed,target,None,'ramp',feedback);cost+=l[0]
    a=np.array(rows)
    metrics=dict(cost=float(cost),force_rmse=float(np.sqrt(np.mean((a[:,10]-.14)**2))),
        track_rmse=float(np.sqrt(np.mean(a[:,11]**2))),normal_rmse=float(np.sqrt(np.mean(a[:,12]**2))),
        minimum_energy=float(a[:,9].min()),minimum_contact=float(a[:,13].min()),
        actuator_saturation_samples=int(np.sum(np.any(abs(a[:,23:25])>ACT+1e-9,axis=1))),
        maximum_rate_ratios=np.max(abs(a[:,16:20])/RATE,axis=0).tolist(),
        minimum_speed=float(a[:,20].min()),force_peak=float(a[:,10].max()),
        tank_violations=int(np.sum(a[:,9]<.001-1e-9)),gain_violations=int(np.sum((a[:,5:9]<LOW-1e-8)|(a[:,5:9]>HIGH+1e-8))),
        controller_seconds=seconds,solver_queries=len(success),solver_successes=sum(s['success'] for s in success))
    if save:
        np.savetxt(HERE/(save+'.csv'),a,delimiter=',',header='t,x,z,vx,vz,Kt,Kn,Dt,Dn,tank,Fn,ex,ez,delta,cmdx,cmdz,rateKt,rateKn,rateDt,rateDn,vt,ell,penalty,intendedx,intendedz',comments='')
        (HERE/(save+'_solve_status.json')).write_text(json.dumps(success,indent=2),encoding='utf8')
    return metrics,a

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cases',default='loose,coupled,budget,friction');args=parser.parse_args()
    contract=dict(cases=CASES,gain_preference=GP.tolist(),gain_penalty=.1,track_scale=.0001,
        comparison='all have current accurate material, geometry and physical gain/tank/actuator interfaces',
        status='exploratory failure search; not a formal benchmark or original HJB implementation')
    (HERE/'frozen_config.json').write_text(json.dumps(contract,indent=2),encoding='utf8')
    all_results={}
    for name in args.cases.split(','):
        configure(name);all_results[name]={}
        # Strong fixed candidates span the entire path, both speeds; explicit record.
        candidates=[GP]+[np.array([kt,kn,dt,dn]) for kt,kn,dt,dn in
            [(400,5,8,6),(800,5,8,6),(400,80,8,6),(800,80,8,6),(1200,80,8,6),(400,200,8,6),(800,200,15,15)]]
        tune=[]
        for g in candidates:
            r,_=run(.006,'fixed',dict(gains=g));tune.append(dict(gains=g.tolist(),metrics=r))
        best=min(tune,key=lambda a:a['metrics']['cost'])
        fixed=dict(gains=best['gains']);all_results[name]['fixed_search']=tune
        for speed in [.006,.003]:
            key='fast' if speed==.006 else 'slow';all_results[name][key]={}
            for mode in ['fixed','steady','qp','force','mpc']:
                p=fixed if mode=='fixed' else dict(feedback=3.) if mode=='force' else None
                r,a=run(speed,mode,p,save=name+'_'+key+'_'+mode)
                all_results[name][key][mode]=r
                (HERE/'results.json').write_text(json.dumps(all_results,indent=2),encoding='utf8')
                print(name,key,mode,json.dumps(r),flush=True)
    print('frontier search complete',flush=True)

if __name__=='__main__':main()
