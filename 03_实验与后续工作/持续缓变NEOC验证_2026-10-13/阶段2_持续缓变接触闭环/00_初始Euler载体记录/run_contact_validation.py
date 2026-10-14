from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
from continuous_contact import UniformCarrier,descriptor,advance,metrics,DREF,FD,T

HERE=Path(__file__).resolve().parent
METHODS=['local_NEOC','continuous_NEOC','cached_chord','one_Newton','full_solve','fixed_gain']
RADIUS=.025  # PREDECLARED diagnostic radius, NOT a proved deployment certificate
EPS=1e-5


def save_csv(name,rows,header):
    with (HERE/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)


def quality_audit(c):
    # Independent nonlinear shooting on the SAME uniform material problem.
    import independent_reference as ir
    from continuous_contact import force,src
    def slope(t,d,eta):
        a=src.fast_width(d)
        base=.02/(2*.03)*30000*.03*np.pi/np.log((.02+np.sqrt(.02**2-a*a))/a)
        return force(d,eta[0]),eta[0]*base
    ir.force_slope=slope  # local imported module; original source unchanged
    shoot=ir.ShootingReference((1.2,0.))
    z=np.full(shoot.M,150.);direction=np.sin(np.arange(shoot.M)+.3)
    _,g=shoot.fun(z);epsilon=1e-3
    fd=(shoot.fun(z+epsilon*direction)[0]-shoot.fun(z-epsilon*direction)[0])/(2*epsilon)
    ge=abs(fd-g@direction)/max(abs(fd),1e-10)
    assert ge<1e-4
    fit,shoot_trace,ref=shoot.solve();assert ref['success']
    save_csv('independent_uniform_shooting.csv',shoot_trace,
             ['time_s','d_m','v_m_s','gain','rate','force_N','stage'])
    audits=[]
    for cc in [c,UniformCarrier(49,49,200)]:
        V,_=cc.solve((1.2,0.));x=np.array([DREF,0.,150.,.02]);total=0.;chunks=[]
        # Frozen physical problem: descriptor exactly 1.2 throughout.
        # advance uses a zero-slope ramp with a supplied descriptor via duration
        # override below; implement the frozen rollout independently here.
        from continuous_contact import MASS,DAMP,BR,RU
        dt=.001
        for i in range(round(T/dt)):
            tt=i*dt;u=cc.policy(V,0.,x,(1.2,0.));f=force(x[0],1.2);e=x[0]-DREF
            stage=(e/DREF)**2+(x[1]/.008)**2+(f/FD-1)**2+.1*((x[2]-150)/150)**2+RU*u*u
            # Quality check uses finite-horizon phase, not receding phase.
            phase=min(int(tt/cc.h),cc.N-1)
            d,v,k=x[:3];acc=(FD-f-k*(d-DREF)-(DAMP+BR)*v)/MASS
            geom=cc.geometry(np.array([d+cc.h*v]),np.array([v+cc.h*acc]))
            fraction=min(max(tt/cc.h-phase,0.),1.)
            line=(1-fraction)*cc.lines(V[phase+1],geom)+fraction*cc.lines(V[min(phase+2,cc.N)],geom)
            u=float(cc.minimize(line,np.array([k]))[1][0])
            stage=(e/DREF)**2+(v/.008)**2+(f/FD-1)**2+.1*((k-150)/150)**2+RU*u*u
            def rhs(y):
                ff=force(y[0],1.2);ee=y[0]-DREF
                return np.array([y[1],(FD-ff-y[2]*ee-(DAMP+BR)*y[1])/MASS,
                                 u,.3*(DAMP+BR)*y[1]**2-.5*ee**2*u])
            aa=rhs(x);bb=rhs(x+dt*aa/2);dd=rhs(x+dt*bb/2);ee=rhs(x+dt*dd)
            x+=dt*(aa+2*bb+2*dd+ee)/6;total+=stage*dt
        value=total+((x[0]-DREF)/DREF)**2+(force(x[0],1.2)/FD-1)**2
        audits.append(dict(nd=len(cc.d),nv=len(cc.v),steps=cc.N,
                           cost=float(value),relative_gap=float(value/ref['cost']-1)))
        cc.model_cache.clear()
    return dict(shooting=ref,gradient_relative_error=float(ge),grids=audits)


def simulate(c,method,kind,duration,period):
    start=time.perf_counter();V,P=c.solve((1.,0.),True)
    initial_s=time.perf_counter()-start;prepare_s=0.;maintenance_s=0.;updates_s=0.;action_s=0.
    anchor=V.copy();anchor_eta=1.;last_eta=1.;Cw=None;dphi=None
    solves=1;derivatives=0;updates=0;residual_evals=0;checks=0
    if method=='local_NEOC':
        start=time.perf_counter();Cw,dphi=c.sensitivity(V,P,(1.,0.),eps=EPS)
        prepare_s+=time.perf_counter()-start;derivatives+=1
    elif method=='fixed_gain':
        # No optimization is needed by this comparator.
        initial_s=0.;solves=0
    x=np.array([DREF,0.,150.,.02]);cost=0.;traces=[];events=[]
    count=round(duration/period)
    for q in range(count):
        now=q*period;eta=float(descriptor(now,duration,kind));reset=False
        start=time.perf_counter()
        if q and method=='local_NEOC':
            if abs(eta-anchor_eta)>RADIUS+1e-12:
                V,P=c.solve((eta,0.),True);solves+=1
                anchor=V.copy();anchor_eta=eta
                Cw,dphi=c.sensitivity(V,P,(eta,0.),eps=EPS);derivatives+=1;reset=True
            else:
                phi=c.model((eta,0.))[1];oldphi=c.model((anchor_eta,0.))[1]
                V=anchor+(phi-oldphi)[None,:]+Cw*(eta-anchor_eta)
                V[-1]=phi
        elif q and method=='continuous_NEOC':
            # Refresh derivative at the actual previous approximate point.
            # These global passes COMPUTE the derivative, not acceptance checks.
            _,P=c.greedy_system(V,(last_eta,0.))
            Cw,dphi=c.sensitivity(V,P,(last_eta,0.),eps=EPS);derivatives+=1
            phi=c.model((eta,0.))[1];oldphi=c.model((last_eta,0.))[1]
            V=V+(phi-oldphi)[None,:]+Cw*(eta-last_eta);V[-1]=phi
        elif q and method in ('cached_chord','one_Newton'):
            phi=c.model((eta,0.))[1];oldphi=c.model((last_eta,0.))[1]
            base=V+(phi-oldphi)[None,:];base[-1]=phi
            if method=='one_Newton':
                F,newP=c.greedy_system(base,(eta,0.));residual_evals+=1
                V=base+c.linear_solve(F,newP);P=newP
            else:
                F=c.residual(base,(eta,0.));residual_evals+=1
                V=base+c.linear_solve(F,P)
                # Maintain the cache for the next query; counted, never free.
                _,P=c.greedy_system(V,(eta,0.))
        elif q and method=='full_solve':
            V,_=c.solve((eta,0.));solves+=1
        elapsed=time.perf_counter()-start
        if reset:maintenance_s+=elapsed
        else:updates_s+=elapsed
        events.append([now,eta,elapsed,reset,solves,derivatives])
        # Diagnostic ideal instant publication; computation deadlines separately
        # checked. Do not claim latency was simulated in the plant.
        start=time.perf_counter()
        x,chunk,increment=advance(c,x,V,eta,now,period,duration,kind,method=='fixed_gain')
        action_s+=time.perf_counter()-start;traces.append(chunk);cost+=increment
        last_eta=eta;updates+=1
        c.model_cache.clear()
        print(f'{kind} T={duration:g} h={period:g} {method} {q+1}/{count}',flush=True)
    trace=np.concatenate(traces);score=metrics(c,trace,cost,x,duration,kind)
    score.update(dict(method=method,path=kind,duration_s=duration,period_s=period,
        observed_max_descriptor_step=float(max(abs(events[i][1]-events[i-1][1]) for i in range(1,len(events)))),
        anchor_solves=solves,sensitivity_evaluations=derivatives,updates=updates,
        candidate_acceptance_checks=checks,residual_evaluations=residual_evals,
        initial_solve_seconds=initial_s,sensitivity_prepare_seconds=prepare_s,
        maintenance_seconds=maintenance_s,ordinary_update_seconds=updates_s,
        total_strategy_seconds=initial_s+prepare_s+maintenance_s+updates_s,
        rollout_with_action_seconds=action_s,
        maximum_query_seconds=max(row[2] for row in events),
        missed_query_deadlines=sum(row[2]>period for row in events)))
    name=f'{kind}_{duration:g}s_{period:g}s_{method}'
    save_csv('trace_'+name+'.csv',trace,['time_s','true_eta','observed_eta','d_m','v_m_s','gain','energy','force_N','rate','stage'])
    save_csv('events_'+name+'.csv',events,['time_s','eta','update_seconds','anchor_reset','anchor_solves','sensitivity_evaluations'])
    return score


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--skip-audit',action='store_true');args=parser.parse_args()
    c=UniformCarrier()
    cfg=dict(physical_model='Uniform cylinder Hertz contact, E(t)=30000 eta(t), slow real physical drift',
        descriptor_information='exact CURRENT eta sampled equally; no future eta, no MHE or preview',
        finite_horizon_s=T,root_type='fixed nodal finite Bellman family; piecewise smooth, not manuscript C2 collocation',
        radius=RADIUS,radius_status='predeclared diagnostic, NOT rigorously certified physical publish radius',
        grid=c.shape,time_steps=c.N,sensitivity_epsilon=EPS,
        latency='zero publication latency in quality rollouts; actual compute deadlines reported separately',
        full_solver='strong structured backward solve; warm start cannot reduce its fixed backward pass',
        comparisons=METHODS)
    result=dict(config=cfg,quality_audit=None,scenarios=[],status='RUNNING')
    if not args.skip_audit:
        print('Independent quality and mesh audit...',flush=True)
        result['quality_audit']=quality_audit(c)
        print(json.dumps(result['quality_audit']),flush=True)
    scenarios=[('ramp',3.,.25)] if args.pilot else [
        ('ramp',3.,.25),('ramp',6.,.25),('ramp',6.,.5),('cycle',6.,.25)]
    for kind,duration,period in scenarios:
        scores=[]
        for method in METHODS:
            scores.append(simulate(c,method,kind,duration,period))
        reference=next(x for x in scores if x['method']=='full_solve')
        for row in scores:
            row['cost_extra_fraction']=row['cost']/reference['cost']-1
            row['force_rmse_extra_fraction']=row['force_rmse_N']/reference['force_rmse_N']-1
            row['quality_gate_pass']=bool(row['cost_extra_fraction']<=.01 and
                row['force_rmse_extra_fraction']<=.05 and row['sampled_constraints_pass'])
        eligible=[x for x in scores if x['method'] in ('cached_chord','one_Newton','full_solve') and x['quality_gate_pass']]
        strongest=min(eligible,key=lambda x:x['total_strategy_seconds']) if eligible else None
        for row in scores:
            row['strongest_eligible_optimizer']=strongest['method'] if strongest else None
            row['total_strategy_time_reduction']=1-row['total_strategy_seconds']/strongest['total_strategy_seconds'] if strongest else None
            row['screen_pass']=bool(row['quality_gate_pass'] and strongest and row['total_strategy_time_reduction']>=.3)
        result['scenarios'].append(dict(path=kind,duration_s=duration,period_s=period,scores=scores))
        (HERE/('pilot_results.json' if args.pilot else 'results.json')).write_text(json.dumps(result,indent=2),encoding='utf-8')
        print('SCENARIO COMPLETE',json.dumps(scores),flush=True)
    result['status']='COMPLETED_DIAGNOSTIC_NOT_FULL_R15_CERTIFICATION'
    result['source_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in ['continuous_contact.py','run_contact_validation.py']}
    (HERE/('pilot_results.json' if args.pilot else 'results.json')).write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
