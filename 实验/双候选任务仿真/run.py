"""Separate candidate experiments. No literature baseline is assumed here."""
import argparse,json,time
from pathlib import Path
import numpy as np
import physics as p
from estimator import MHE,NOISE
from local_hjb import LocalHJB,NeighborHJB

HERE=Path(__file__).resolve().parent

def cases(n,seed):
    rng=np.random.default_rng(seed);out=[]
    for i in range(n):
        kind=i%4;e0=rng.uniform(.8e5,1.2e5);ef=e0*rng.uniform(1.25,1.75)
        mu=rng.uniform(.10,.20);muf=mu*(rng.uniform(1.7,2.3) if kind in (1,3) else 1.)
        if kind==1:ef=e0
        speed=rng.uniform(.008,.012);boundary=rng.uniform(-.0008,.0008)
        eta=np.array([e0,ef,mu,muf,2.,3.,boundary,0.,0.,-.008,.001,speed]);x,eta=p.initial(eta)
        out.append({'id':i,'kind':kind,'eta':eta.tolist(),'x':x.tolist(),'surface_amplitude':rng.uniform(.00006,.00014) if kind in (2,3) else 0.,
                    'frequency':rng.uniform(.4,.8),'noise_multiplier':2. if kind==3 else 1.,'sensor_delay':.04 if kind==3 else .02})
    return out

def actual_eta(case,t,local=False):
    e=np.array(case['eta']);amp=case['surface_amplitude'];freq=case['frequency'];gamma=amp*np.sin(2*np.pi*freq*t)
    e[7]=gamma;e[8]=0.
    if local:e[8]=amp*2*np.pi*freq*np.cos(2*np.pi*freq*t);e[9]+=e[11]*t
    return e

def episode(case,method,seed,model=None,solver=None,collect=False,timed=False,no_material_preview=False,use_moment=True,legacy_old_region=False):
    rng=np.random.default_rng(seed);noise=rng.normal(size=(201,3))*NOISE*case['noise_multiplier']
    x=np.array(case['x']);base=np.array(case['eta']);sensed=p.wrench(0,x,actual_eta(case,0))+noise[0]
    mhe=MHE(base[10],base[9],base[11],sensed)
    traces=[];obs=[];pending=[];contexts=[];cost=0.;policy=None;neigh=None;timing=[];mhetime=[];reject=0;anchorinfo=None
    action=np.zeros(4);start=time.perf_counter();duration=1.6;dt=.002
    publication=None;action_publication=None;estimate_publication=None;policy_time=0.;estimate_time=0.;skipped=0;expired=0
    estimate=mhe.descriptor(0.)
    if method in ('candidate1','full_value','warm_value','held_value'):
        if solver is None:solver=LocalHJB()
        if method=='candidate1':
            neigh=NeighborHJB(solver);policy,anchorinfo=neigh.anchor(x,estimate)
        else:policy,anchorinfo=solver.solve(x,estimate)
    for j in range(round(duration/dt)):
        now=j*dt;true=actual_eta(case,now)
        if publication is not None and now>=publication[0]:
            _,policy,policy_time=publication;publication=None
        if action_publication is not None and now>=action_publication[0]:
            _,action=action_publication;action_publication=None
        if estimate_publication is not None and now>=estimate_publication[0]:
            _,estimate,estimate_time=estimate_publication;estimate_publication=None
        if j%5==0:pending.append((now,x.copy(),p.wrench(now,x,true)+noise[j//5]))
        while pending and pending[0][0]<=now-case['sensor_delay']+1e-10:obs.append(pending.pop(0))
        if obs:sensed=obs[-1][2]
        if j%50==0 and (publication is not None or estimate_publication is not None):skipped+=1
        if j%50==0 and publication is None and estimate_publication is None:
            new_estimate,mi=mhe.estimate(now,obs,use_moment=use_moment);mhetime.append(mi['seconds'])
            if timed and method in ('candidate2','no_preview'):estimate_publication=(now+mi['seconds'],new_estimate,now)
            else:estimate=new_estimate;estimate_time=now
            if no_material_preview:
                if legacy_old_region:
                    estimate=estimate.copy();estimate[[1,3,5]]=estimate[[0,2,4]]
                else:estimate=p.current_equivalent(x,estimate)
            if method in ('candidate1','full_value','warm_value') and j>0:
                tic=time.perf_counter()
                if method in ('full_value','warm_value'):
                    candidate,ci=solver.solve(x,estimate,policy['U'],maxiter=2 if method=='warm_value' else 12)
                    if ci['feasible']:
                        if not timed:policy=candidate;policy_time=now
                    else:reject+=1
                else:
                    candidate,ci=neigh.predict(estimate,x)
                    # Check the predicted feedback on the current state before
                    # publication; this is a sampled nominal check, not B0 cert.
                    Z,U,val,feasible=solver.rollout(x,estimate,candidate['U'],candidate['K'],candidate['Z'],np.zeros_like(candidate['U']))
                    if feasible:
                        # A neighboring prediction is an initial policy, not an
                        # exact optimal solution after a large descriptor jump.
                        candidate,ci=solver.solve(x,estimate,U,maxiter=2)
                        if ci['feasible']:
                            if not timed:policy=candidate;policy_time=now
                        else:reject+=1
                    else:
                        candidate,ci=solver.solve(x,estimate,policy['U'],maxiter=12);reject+=1
                        if ci['feasible']:
                            if not timed:policy=candidate;policy_time=now
                if method=='candidate1' and ci.get('feasible',False):neigh.reanchor(candidate)
                elapsed=time.perf_counter()-tic;timing.append(elapsed)
                if timed and ci.get('feasible',True):publication=(now+mi['seconds']+elapsed,candidate,now)
        if j%10==0:
            tic=time.perf_counter()
            current_estimate=estimate.copy();age=now-estimate_time
            current_estimate[7]+=current_estimate[8]*age;current_estimate[9]+=current_estimate[11]*age
            if method in ('candidate1','full_value','warm_value','held_value'):
                if now-policy_time>=p.HORIZON:action=np.zeros(4);expired+=1
                else:action=solver.action(policy,now-policy_time,x)
            elif method in ('candidate2','no_preview'):
                import torch
                from learning import features
                with torch.no_grad():new_action=model(torch.tensor(features(x,current_estimate,model.current_equivalent)[None],dtype=torch.float32)).numpy()[0]
                if timed:action_publication=(now+time.perf_counter()-tic,new_action)
                else:action=new_action
            elif method=='random':action=rng.uniform(-.2,.2,4)
            elif method=='jin_ms':
                new_action=model.decide(x,sensed,base[11])
                if timed:action_publication=(now+time.perf_counter()-tic,new_action)
                else:action=new_action
            else:action=np.zeros(4)
            if method in ('candidate2','no_preview','jin_ms'):timing.append(time.perf_counter()-tic)
            if collect:contexts.append((x.copy(),current_estimate.copy(),actual_eta(case,now,local=True)))
        dx,ell,executed=p.flow(now,x,true,action);margins=p.margins(now,x,true)
        traces.append(np.r_[now,x,p.wrench(now,x,true),executed,margins.min(),estimate])
        # Surface motion is regenerated at midpoint for the independent plant.
        mid=x+.5*dt*dx;dx,ell,executed=p.flow(now+.5*dt,mid,actual_eta(case,now+.5*dt),action)
        x=x+dt*dx;cost+=dt*ell
        if not np.isfinite(x).all() or np.max(abs(x[:4]))>1.:
            raise FloatingPointError((method,case['id'],now,x))
    a=np.asarray(traces);f=a[:,10];mins=a[:,17]
    endeta=actual_eta(case,duration)
    endeta[9]+=endeta[11]*(duration-p.HORIZON)
    terminal=p.terminal(x,endeta)
    metrics={'case':case['id'],'kind':case['kind'],'cost':float(cost+terminal),'running_cost':float(cost),'terminal_cost':float(terminal),
             'force_rmse_N':float(np.sqrt(np.mean((f-p.FD)**2))),'force_peak_N':float(f.max()),
             'tangent_rmse_mm':float(np.sqrt(np.mean(a[:,1]**2))*1000),'normal_rmse_mm':float(np.sqrt(np.mean(a[:,2]**2))*1000),
             'constraint_violating_samples':int(np.sum(mins< -1e-8)),'minimum_margin_mixed_units':float(mins.min()),
             'minimum_tank_J':float(a[:,9].min()),'decision_median_ms':float(np.median(timing)*1000) if timing else 0.,
             'decision_p95_ms':float(np.percentile(timing,95)*1000) if timing else 0.,'decision_total_seconds':float(sum(timing)),
             'mhe_total_seconds':float(sum(mhetime)) if method!='jin_ms' else 0.,'publication_rejections':reject,'initial_preparation':anchorinfo,
             'wall_seconds':time.perf_counter()-start,'computational_delay_applied':timed,
             'skipped_busy_updates':skipped,'expired_policy_queries':expired}
    if method=='jin_ms':
        metrics['source_tank_min_J']=min(z['source_tank_min'] for z in model.log)
        metrics['source_optimizer_failures']=sum(not z['success'] for z in model.log)
        metrics['source_optimizer_calls']=len(model.log)
    return metrics,a,contexts

def main(args):
    target=HERE/args.output
    if target.exists():raise FileExistsError(target)
    data=cases(args.n,args.seed);model=None;solver=None;setup_start=time.perf_counter()
    if args.softening:
        for case in data:
            if case['kind']!=1:case['eta'][1]=case['eta'][0]**2/case['eta'][1]
    if args.fast_motion:
        for case in data:
            if case['surface_amplitude']>0:case['surface_amplitude']*=2;case['frequency']*=2
    if args.method in ('candidate2','no_preview'):
        from learning import load
        model=load(args.model)
    if args.method in ('candidate1','full_value','warm_value','held_value'):solver=LocalHJB()
    batch_setup_seconds=time.perf_counter()-setup_start
    results=[];traces={};contexts=[]
    for case in data:
        construction=time.perf_counter()
        if args.method=='jin_ms':
            from baseline import JinMS
            model=JinMS(args.baseline_axis_regularizers if args.baseline_axis_regularizers is not None else args.baseline_regularizer)
        construction=time.perf_counter()-construction
        m,tr,ctx=episode(case,args.method,args.seed*100+case['id'],model,solver,args.collect,args.timed,args.no_preview,not args.no_moment,args.legacy_old_region_mask)
        m['controller_construction_seconds']=construction
        results.append(m);traces[str(case['id'])]=tr;contexts.extend(ctx)
        print(json.dumps({'method':args.method,**m}),flush=True)
    target.write_text(json.dumps({'method':args.method,'seed':args.seed,'n':args.n,'timed':args.timed,'batch_setup_seconds':batch_setup_seconds,'no_preview':args.no_preview,'no_moment':args.no_moment,
                                 'softening':args.softening,'fast_motion':args.fast_motion,'legacy_old_region_mask':args.legacy_old_region_mask,'cases':data,'results':results},indent=2),encoding='utf-8')
    np.savez_compressed(target.with_name(target.stem+'_traces.npz'),**traces)
    if contexts:
        xx,ee,tt=zip(*contexts);np.savez_compressed(target.with_name(target.stem+'_contexts.npz'),x=xx,estimate=ee,truth=tt)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--method',required=True);ap.add_argument('--n',type=int,default=4);ap.add_argument('--seed',type=int,default=16201)
    ap.add_argument('--model');ap.add_argument('--output',required=True);ap.add_argument('--collect',action='store_true')
    ap.add_argument('--timed',action='store_true');ap.add_argument('--no-preview',action='store_true');ap.add_argument('--no-moment',action='store_true')
    ap.add_argument('--softening',action='store_true');ap.add_argument('--fast-motion',action='store_true')
    ap.add_argument('--legacy-old-region-mask',action='store_true',help='reproduce rejected ablation only; not a valid no-preview method')
    ap.add_argument('--baseline-axis-regularizers',type=float,nargs=2)
    ap.add_argument('--baseline-regularizer',type=float,default=1e-5);main(ap.parse_args())
