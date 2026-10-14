"""Stage A: causal 2-D scanning, cylinder Hertz contact, actual gain dynamics.

Direct-shooting MPC is a finite-horizon local optimization reference, not an
HJB global optimum. Targets are decision variables; physical gains obey bounded
rates. Batched central differences provide an independently auditable gradient.
"""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
CFG=dict(mass=.2,robot_damping=1.,radius=.03,length=.02,lref=.02,
         force_target=.14,scan_length=.024,speeds=[.003,.006],dt=.002,
         query=.1,horizon=.4,shoot_dt=.008,blocks=2,
         gain_lower=[80.,5.,2.,2.],gain_upper=[1600.,800.,30.,30.],
         gain_rate=[4000.,2000.,100.,100.],initial_energy=.05,energy_min=.001,
         force_scale=.02,track_scale=.0005,normal_scale=.00025,
         velocity_scale=.008,effort_scale=.2,effort_weight=.01,
         rate_weight=.0001,force_ceiling=.23,track_tolerance=.0005,
         force_rmse_tolerance=.02,material_min=30000.,material_max=60000.,
         surface_amplitude=.0005,surface_wavelength=.048,
         optimizer_maxiter=35,optimizer_ftol=1e-7,optimizer_starts=2,
         development='ramp_slow',verification=['ramp_fast','cycle_slow','cycle_fast'],
         substantive_fraction=.10)
LOW=np.array(CFG['gain_lower']);HIGH=np.array(CFG['gain_upper'])
RATE=np.array(CFG['gain_rate']);SPAN=HIGH-LOW

# Tabulated monotone inverse of the exact cylindrical reference-plane relation.
AA=np.geomspace(1e-8,.012,30000)
ss=np.sqrt(CFG['lref']**2-AA**2)
DD=(CFG['lref']**2-CFG['lref']*ss+AA**2*np.log((CFG['lref']+ss)/AA))/(2*CFG['radius'])
def contact(delta,E):
    a=np.interp(np.maximum(delta,0.),DD,AA)
    return np.where(delta>0,CFG['length']*E*np.pi*a*a/(4*CFG['radius']),0.)
DNOM=float(np.interp(CFG['force_target'],contact(DD,30000.),DD))

def surface(x):
    w=2*np.pi/CFG['surface_wavelength'];a=CFG['surface_amplitude']
    return a*np.sin(w*x),a*w*np.cos(w*x)

def material(x,kind):
    p=np.clip(x/CFG['scan_length'],0.,1.)
    fraction=p if kind=='ramp' else .5*(1-np.cos(2*np.pi*p))
    return 30000.+30000.*fraction,.2+.3*fraction,2.+2.*fraction

def geometry(y,t,speed,env,kind):
    x,z,vx,vz=y[:,:4].T
    h,s=surface(x);r=np.sqrt(1+s*s)
    nx,nz=-s/r,1/r;tx,tz=1/r,s/r
    # Signed distance to the local tangent plane; small-slope cylinder projection.
    delta=(z-h)/r;vt=vx*tx+vz*tz
    E,mu,kd=material(x,kind) if env is None else env
    # Local osculating-plane Hertz: combine probe and known surface curvature.
    w=2*np.pi/CFG['surface_wavelength']
    curvature=-CFG['surface_amplitude']*w*w*np.sin(w*x)/(1+s*s)**1.5
    relative=1+CFG['radius']*curvature
    fn=relative*contact(delta/relative,E);fr=mu*fn*np.tanh(vt/.0002)+kd*vt
    xr=speed*t;hr,sr=surface(xr)
    e=np.column_stack((x-xr,z-hr-DNOM))
    ev=np.column_stack((vx-speed,vz-speed*sr))
    return fn,fr,e,ev,np.column_stack((nx,nz)),np.column_stack((tx,tz)),delta,vt

def dynamics(y,t,speed,target,env,kind,force_feedback=0.):
    fn,fr,e,ev,n,tan,delta,vt=geometry(y,t,speed,env,kind)
    gains=y[:,4:8];u=np.clip((target-gains)/.04,-RATE,RATE)
    command=(CFG['force_target']+force_feedback*(CFG['force_target']-fn))[:,None]*n
    command-=gains[:,:2]*e+gains[:,2:]*ev
    fcontact=-fn[:,None]*n-fr[:,None]*tan
    acc=(command+fcontact-CFG['robot_damping']*y[:,2:4])/CFG['mass']
    diss=np.sum(gains[:,2:]*ev**2,axis=1)+CFG['robot_damping']*np.sum(y[:,2:4]**2,axis=1)
    exchange=.5*np.sum(e**2*u[:,:2],axis=1)
    # Common tank protection; this is a physical constraint, not result acceptance.
    limited=(y[:,8]<=CFG['energy_min']+.0001)&(exchange>.3*diss)
    u[limited,:2]=np.minimum(u[limited,:2],0.)
    energy=.3*diss-.5*np.sum(e**2*u[:,:2],axis=1)
    rhs=np.column_stack((y[:,2:4],acc,u,energy))
    ell=((fn-CFG['force_target'])/CFG['force_scale'])**2
    ell+= (e[:,0]/CFG['track_scale'])**2+.2*(e[:,1]/CFG['normal_scale'])**2
    ell+= .02*np.sum((ev/CFG['velocity_scale'])**2,axis=1)
    ell+=CFG['effort_weight']*np.sum((command/CFG['effort_scale'])**2,axis=1)
    ell+=CFG['rate_weight']*np.sum((u/RATE)**2,axis=1)
    # Optimizer feasibility penalties; hard outcomes evaluated independently.
    penalty=1e4*np.maximum(-delta/.0001,0.)**2
    penalty+=1e4*np.maximum((fn-CFG['force_ceiling'])/.02,0.)**2
    penalty+=1e4*np.maximum((CFG['energy_min']-y[:,8])/.001,0.)**2
    return rhs,ell,penalty,(fn,e,delta,command,u,vt)

def rk4(y,t,h,speed,target,env,kind,feedback=0.):
    a,l,p,_=dynamics(y,t,speed,target,env,kind,feedback)
    b,lb,pb,_=dynamics(y+h*a/2,t+h/2,speed,target,env,kind,feedback)
    c,lc,pc,_=dynamics(y+h*b/2,t+h/2,speed,target,env,kind,feedback)
    d,ld,pd,_=dynamics(y+h*c,t+h,speed,target,env,kind,feedback)
    return y+h*(a+2*b+2*c+d)/6,h*(l+2*lb+2*lc+ld)/6,h*(p+2*pb+2*pc+pd)/6

def shooting(q,state,t,speed,env,kind,dt=None):
    q=np.atleast_2d(q);batch=len(q);y=np.repeat(state[None,:],batch,axis=0)
    targets=LOW+q.reshape(batch,CFG['blocks'],4)*SPAN
    h=CFG['shoot_dt'] if dt is None else dt;steps=round(CFG['horizon']/h)
    value=np.zeros(batch)
    for j in range(steps):
        block=min(CFG['blocks']-1,j*CFG['blocks']//steps)
        y,l,p=rk4(y,t+j*h,h,speed,targets[:,block],env,kind)
        value+=l+p
    # Short terminal tail, equal objective across all controls.
    _,l,p,_=dynamics(y,t+CFG['horizon'],speed,targets[:,-1],env,kind)
    return value+.05*(l+p)

def objective(q,state,t,speed,env,kind):
    eps=2e-5;eye=np.eye(len(q))*eps
    vals=shooting(np.vstack((q,q+eye,q-eye)),state,t,speed,env,kind)
    return float(vals[0]),(vals[1:1+len(q)]-vals[1+len(q):])/(2*eps)

def initial(speed,gains):
    h,s=surface(0.)
    return np.array([0.,h+DNOM,speed,speed*s,*gains,CFG['initial_energy']])

def run(kind,speed,mode,params,dt=None,save=None,maxiter=None):
    h=CFG['dt'] if dt is None else dt
    duration=CFG['scan_length']/speed;steps=round(duration/h);query=round(CFG['query']/h)
    gains=np.array(params['gains']);y=initial(speed,gains)[None,:]
    q=np.tile((gains-LOW)/SPAN,CFG['blocks']);target=gains[None,:]
    cost=0.;rows=[];solves=[];feedback=params.get('feedback',0.) if mode=='force' else 0.
    for j in range(steps):
        t=j*h
        if j%query==0:
            env=material(y[0,0],kind)
            if mode=='mpc':
                tic=time.perf_counter()
                # Rate saturation creates flat decision regions. Retain the fixed
                # incumbent and rank coarse alternatives before local refinement.
                candidates=[q,np.tile((gains-LOW)/SPAN,CFG['blocks']),
                    np.tile((y[0,4:8]-LOW)/SPAN,CFG['blocks'])]
                for kn in [5.,60.,400.]:
                    for dn in [2.,15.,30.]:
                        g=gains.copy();g[1]=kn;g[3]=dn
                        candidates.append(np.tile((g-LOW)/SPAN,CFG['blocks']))
                candidates=np.asarray(candidates)
                ranked=np.argsort(shooting(candidates,y[0].copy(),t,speed,env,kind))
                local=[]
                for start in ranked[:CFG['optimizer_starts']]:
                    local.append(minimize(objective,candidates[start],args=(y[0].copy(),t,speed,env,kind),jac=True,
                        method='L-BFGS-B',bounds=[(0.,1.)]*len(q),
                        options=dict(maxiter=maxiter or CFG['optimizer_maxiter'],ftol=CFG['optimizer_ftol'],gtol=2e-5,maxls=20)))
                res=min(local,key=lambda r:r.fun)
                q=res.x;target=(LOW+q[:4]*SPAN)[None,:]
                solves.append(dict(t=t,cost=float(res.fun),success=bool(res.success),
                    iterations=int(res.nit),message=str(res.message),seconds=time.perf_counter()-tic,
                    q=q.tolist(),state=y[0].tolist(),env=[float(a) for a in env],
                    local_objectives=[float(r.fun) for r in local]))
                # Causal warm start; no future material supplied to shooting.
            elif mode=='reactive':
                E,mu,kd=env
                # Material-aware impedance schedule with force-error response.
                fn=geometry(y,t,speed,None,kind)[0][0]
                normalized=(E-30000.)/30000.
                goal=gains.copy()
                goal[0]*=1+params['tangent_slope']*normalized
                goal[1]*=np.exp(-params['normal_slope']*normalized-params['force_gain']*(fn/CFG['force_target']-1))
                goal[2:]=params['damping_ratio']*2*np.sqrt(CFG['mass']*goal[:2])
                target=np.clip(goal,LOW,HIGH)[None,:]
        _,ell,penalty,out=dynamics(y,t,speed,target,None,kind,feedback)
        fn,e,delta,command,u,vt=out
        rows.append([t,*y[0],fn[0],e[0,0],e[0,1],delta[0],*command[0],*u[0],vt[0],ell[0],penalty[0]])
        y,l,p=rk4(y,t,h,speed,target,None,kind,feedback);cost+=float(l[0])
        if not np.all(np.isfinite(y)):raise FloatingPointError('nonfinite trajectory')
    trace=np.asarray(rows)
    # Columns: t, x,z,vx,vz,Kt,Kn,Dt,Dn,tank,Fn,ex,ez,delta,command[2],u[4],vt,ell,penalty.
    rmsef=float(np.sqrt(np.mean((trace[:,10]-CFG['force_target'])**2)))
    rmsex=float(np.sqrt(np.mean(trace[:,11]**2)))
    result=dict(cost=cost,force_rmse_N=rmsef,track_rmse_m=rmsex,
        normal_rmse_m=float(np.sqrt(np.mean(trace[:,12]**2))),
        force_peak_N=float(np.max(trace[:,10])),track_peak_m=float(np.max(abs(trace[:,11]))),
        min_contact_m=float(np.min(trace[:,13])),min_energy_J=float(np.min(trace[:,9])),
        contact_loss_samples=int(np.sum(trace[:,13]<=0)),force_violation_samples=int(np.sum(trace[:,10]>CFG['force_ceiling'])),
        gain_violation_samples=int(np.sum((trace[:,5:9]<LOW-1e-8)|(trace[:,5:9]>HIGH+1e-8))),
        rate_violation_samples=int(np.sum(abs(trace[:,16:20])>RATE+1e-8)),
        task_tolerances_pass=bool(rmsef<=CFG['force_rmse_tolerance'] and rmsex<=CFG['track_tolerance']),
        solve_seconds=float(sum(a['seconds'] for a in solves)),optimizer_successes=sum(a['success'] for a in solves),
        optimizer_queries=len(solves),final_state=y[0].tolist())
    if save:
        np.savetxt(HERE/(save+'.csv'),trace,delimiter=',',header='t,x,z,vx,vz,Kt,Kn,Dt,Dn,tank,Fn,ex,ez,delta,ux,uz,rateKt,rateKn,rateDt,rateDn,vt,ell,penalty',comments='')
        (HERE/(save+'_solves.json')).write_text(json.dumps(solves,indent=2),encoding='utf8')
    return result,trace,solves

def batch_cost(mode,parameters):
    """Evaluate independent tuning candidates together; no shared controller data."""
    speed=.003;h=CFG['dt'];steps=round(CFG['scan_length']/speed/h);query=round(CFG['query']/h)
    gains=np.array([p['gains'] for p in parameters]);y=np.array([initial(speed,g) for g in gains])
    target=gains.copy();total=np.zeros(len(parameters))
    feedback=np.array([p.get('feedback',0.) for p in parameters]) if mode=='force' else 0.
    for j in range(steps):
        t=j*h
        if mode=='reactive' and j%query==0:
            E,mu,kd=material(y[:,0],'ramp');normalized=(E-30000.)/30000.
            fn=geometry(y,t,speed,None,'ramp')[0]
            target=gains.copy()
            target[:,0]*=1+np.array([p['tangent_slope'] for p in parameters])*normalized
            target[:,1]*=np.exp(-np.array([p['normal_slope'] for p in parameters])*normalized-np.array([p['force_gain'] for p in parameters])*(fn/CFG['force_target']-1))
            target[:,2:]=np.array([p['damping_ratio'] for p in parameters])[:,None]*2*np.sqrt(CFG['mass']*target[:,:2])
            target=np.clip(target,LOW,HIGH)
        y,l,p=rk4(y,t,h,speed,target,None,'ramp',feedback);total+=l+p
    return total

def tune():
    records=[];best={}
    # Same development path only. Deterministic broad gain grid, then continuous refinement.
    for kt in [150.,600.,1400.]:
        for kn in [8.,60.,400.]:
            for ratio in [.5,1.,1.5]:
                g=np.clip([kt,kn,ratio*2*np.sqrt(.2*kt),ratio*2*np.sqrt(.2*kn)],LOW,HIGH)
                records.append(dict(mode='fixed',params=dict(gains=g.tolist())))
    costs=batch_cost('fixed',[a['params'] for a in records])
    for a,cost in zip(records,costs):a['tuning_cost']=float(cost)
    seed=min(records,key=lambda a:a['tuning_cost'])['params']['gains']
    def fixed_obj(q):
        ep=1e-4;eye=np.eye(4)*ep;points=np.vstack((q,q+eye,q-eye))
        costs=batch_cost('fixed',[dict(gains=(LOW+p*SPAN).tolist()) for p in points])
        return costs[0],(costs[1:5]-costs[5:])/(2*ep)
    res=minimize(fixed_obj,(np.array(seed)-LOW)/SPAN,jac=True,method='L-BFGS-B',bounds=[(0,1)]*4,options=dict(maxiter=20,ftol=1e-8))
    fixed=dict(gains=(LOW+res.x*SPAN).tolist())
    records.append(dict(mode='fixed',params=fixed,tuning_cost=float(res.fun)));best['fixed']=fixed
    # Material-aware schedule plus response to measured force. Shared starting gains.
    for ns in [0.,1.,3.]:
        for fg in [0.,2.,6.]:
            for ratio in [.7,1.2]:
                p=dict(gains=fixed['gains'],normal_slope=ns,force_gain=fg,tangent_slope=1.,damping_ratio=ratio)
                records.append(dict(mode='reactive',params=p))
    subset=[a for a in records if a['mode']=='reactive'];costs=batch_cost('reactive',[a['params'] for a in subset])
    for a,cost in zip(subset,costs):a['tuning_cost']=float(cost)
    best['reactive']=min(subset,key=lambda a:a['tuning_cost'])['params']
    for feedback in [0.,1.,3.,8.,20.]:
        p=dict(gains=fixed['gains'],feedback=feedback)
        records.append(dict(mode='force',params=p))
    subset=[a for a in records if a['mode']=='force'];costs=batch_cost('force',[a['params'] for a in subset])
    for a,cost in zip(subset,costs):a['tuning_cost']=float(cost)
    best['force']=min(subset,key=lambda a:a['tuning_cost'])['params']
    (HERE/'tuning.json').write_text(json.dumps(dict(records=records,selected=best),indent=2),encoding='utf8')
    print('tuning complete',best,flush=True)
    return best

def audit():
    # Force inverse and mechanical+spring+tank energy identity in stationary reference.
    inv=abs(float(contact(DNOM,30000.))-CFG['force_target'])
    state=initial(.003,[600.,60.,12.,6.]);env=material(state[0],'ramp')
    q=np.tile((state[4:8]-LOW)/SPAN,CFG['blocks']);q=np.clip(q+.1,.05,.95)
    value,g=objective(q,state,0.,.003,env,'ramp')
    rng=np.random.default_rng(15);d=rng.normal(size=len(q));d/=np.linalg.norm(d)
    diffs={}
    for ep in [1e-3,1e-4,1e-5]:
        fd=(shooting(q+ep*d,state,0.,.003,env,'ramp')[0]-shooting(q-ep*d,state,0.,.003,env,'ramp')[0])/(2*ep)
        diffs[str(ep)]=float(abs(fd-g@d)/max(1.,abs(fd)))
    y=state[None,:].copy();target=(LOW+q[:4]*SPAN)[None,:]
    rhs,_,_,out=dynamics(y,0.,.003,target,env,'ramp')
    fn,fr,e,ev,n,tan,delta,vt=geometry(y,0.,.003,env,'ramp');command=out[3];u=out[4]
    # E=kinetic + .5 e'Ke + tank. Moving reference/feedforward remain port powers.
    lhs=np.sum(CFG['mass']*y[:,2:4]*rhs[:,2:4])+np.sum(e*y[:,4:6]*ev)+.5*np.sum(e*e*u[:,:2])+rhs[0,8]
    v=y[:,2:4];refv=v-ev;fc=-fn[:,None]*n-fr[:,None]*tan
    diss=np.sum(y[:,6:8]*ev**2)+CFG['robot_damping']*np.sum(v**2)
    expected=np.sum(v*(CFG['force_target']*n+fc))-.7*diss-np.sum(refv*(y[:,4:6]*e+y[:,6:8]*ev))
    result=dict(force_inverse_abs_N=inv,directional_gradient_relative_errors=diffs,
                mechanical_tank_energy_identity_abs=float(abs(lhs-expected)))
    (HERE/'model_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print('audit',result,flush=True)
    assert inv<1e-6 and min(diffs.values())<1e-4 and abs(lhs-expected)<1e-10

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['all','tune','verify'],default='all');args=parser.parse_args()
    # Write the experiment contract before any comparisons run.
    (HERE/'frozen_config.json').write_text(json.dumps(dict(config=CFG,dnom=DNOM,
        information='current E,mu,kd; known reference geometry; frozen material in horizon',
        decisions='two target blocks, actual rate-limited gain/tank dynamics',
        criterion='>=10% cost and >=10% force or track improvement over strongest simple baseline; no extra violations; or additional task-tolerance pass',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf8')
    audit()
    if args.mode in ['all','tune']:best=tune()
    else:best=json.loads((HERE/'tuning.json').read_text(encoding='utf8'))['selected']
    if args.mode=='tune':return
    results={}
    for kind in ['ramp','cycle']:
        for speed in CFG['speeds']:
            name=kind+('_slow' if speed==CFG['speeds'][0] else '_fast');results[name]={}
            for mode in ['fixed','reactive','force','mpc']:
                p=best['fixed'] if mode=='mpc' else best[mode]
                r,trace,solves=run(kind,speed,mode,p,save=name+'_'+mode)
                results[name][mode]=r
                (HERE/'results.json').write_text(json.dumps(results,indent=2),encoding='utf8')
                print(name,mode,json.dumps(r),flush=True)
    print('all comparisons complete',flush=True)

if __name__=='__main__':main()
