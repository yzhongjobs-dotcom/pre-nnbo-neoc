"""One-step rollout improvement of a cheap feedback policy, a strong rival.

This is an existing control principle, not proposed as a new R15 contribution.
It gets the same model/state and returns actual gain-rate actions. No optimal
family answers are queried by the controller. Candidate rollout is its own
planning computation and is timed. Evaluator optimum is separate.
"""
import json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize
import compiled_flow
from dynamic_tradeoff import Problem,dump,ROOT


@njit(cache=True)
def base_rate(y,pars,p,h,lo,hi):
    ak,ad,se,sv,bias,tau=pars
    gate=(y[0]/se)**2+(y[1]/sv)**2;gate=gate/(1+gate)
    sign=np.tanh(bias*y[0]*y[1]/(se*sv))
    factor=max(0.,min(gate*(1+sign)/2,1.))
    desired=np.array([p[13]+ak*factor,p[14]+ad*factor])
    z=np.empty(2)
    for i in range(2):
        desired[i]=max(lo[i],min(desired[i],hi[i]))
        zi=(desired[i]-y[2+i])/max(tau,h)/p[24+i]
        z[i]=max(-1.,min(zi,1.))
        z[i]=max((lo[i]-y[2+i])/(h*p[24+i]),min(z[i],(hi[i]-y[2+i])/(h*p[24+i])))
    reserve=max(y[4]-p[4],0.);needed=.5*y[0]**2*max(z[0]*p[24],0.)*h
    if needed>reserve and needed>0:z[0]*=reserve/needed
    return z


@njit(cache=True)
def future(first,y0,pars,th,p,n,sub,dt,h,lo,hi):
    y=y0.copy();cost=0.;minimum=1e99
    for j in range(n):
        z=first.copy() if j==0 else base_rate(y,pars,p,h,lo,hi)
        for i in range(2):
            z[i]=max((lo[i]-y[2+i])/(h*p[24+i]),min(z[i],(hi[i]-y[2+i])/(h*p[24+i])))
        for _ in range(sub):
            yn,ell,ya,yb,yc=compiled_flow.step(y,z,th,p,dt)
            for stage in (y,ya,yb,yc,yn):
                mar=compiled_flow.margins(stage,th,p)
                for v in mar:minimum=min(minimum,v)
            y=yn;cost+=ell
    cost+=p[23]*((y[0]/p[9])**2+p[19]*(th*y[0]/p[11])**2)
    return cost,minimum


@njit(cache=True)
def decide(y,pars,th,p,n,sub,dt,h,lo,hi):
    base=base_rate(y,pars,p,h,lo,hi)
    best=base.copy();cost,margin=future(base,y,pars,th,p,n,sub,dt,h,lo,hi)
    if margin < -1e-8:cost=1e99
    grid=np.array([-1.,-.5,0.,.5,1.])
    count=1
    for a in grid:
        for b in grid:
            z=np.array([a,b])
            for i in range(2):
                z[i]=max((lo[i]-y[2+i])/(h*p[24+i]),min(z[i],(hi[i]-y[2+i])/(h*p[24+i])))
            val,mar=future(z,y,pars,th,p,n,sub,dt,h,lo,hi);count+=1
            if mar>=-1e-8 and val<cost:cost=val;best=z.copy()
    return best,cost,count


def policy(p,pars,refine=False):
    y=p.y0.copy();zs=[];seconds=[];predicted=[]
    # JIT preparation measured separately; used for every method timing report.
    tic=time.perf_counter()
    decide(y,pars,p.theta,p.params,p.n,p.sub,p.dt,p.h,p.lo,p.hi)
    setup=time.perf_counter()-tic
    for j in range(p.n):
        tic=time.perf_counter()
        z,value,count=decide(y,pars,p.theta,p.params,p.n-j,p.sub,p.dt,p.h,p.lo,p.hi)
        if refine:
            # Refine only the current two actions; future actions remain the
            # fixed cheap policy. Keep the feasible grid incumbent on failure.
            cache={}
            def evaluate(q):
                key=tuple(q)
                if key not in cache:cache[key]=future(q,y,pars,p.theta,p.params,p.n-j,p.sub,p.dt,p.h,p.lo,p.hi)
                return cache[key]
            bounds=list(zip(np.maximum(-1.,(p.lo-y[2:4])/(p.h*p.rate)),
                            np.minimum(1.,(p.hi-y[2:4])/(p.h*p.rate))))
            res=minimize(lambda q:evaluate(q)[0],z,method='SLSQP',bounds=bounds,
                    constraints=[dict(type='ineq',fun=lambda q:evaluate(q)[1])],
                    options=dict(maxiter=35,ftol=1e-9))
            candidate,margin=evaluate(res.x)
            if res.success and margin>=-1e-8 and candidate<value:z,value=res.x,candidate
        seconds.append(time.perf_counter()-tic)
        if value>1e90:raise RuntimeError('No feasible rollout incumbent')
        zs.append(z);predicted.append(value)
        for _ in range(p.sub):y=compiled_flow.step(y,z,p.theta,p.params,p.dt)[0]
    return np.array(zs).ravel(),dict(setup_or_warmup_seconds=setup,
            decision_seconds=seconds,total_control_seconds=sum(seconds),
            worst_decision_seconds=max(seconds),candidate_count_per_decision=26)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--refine',action='store_true');args=parser.parse_args()
    audit=json.loads((ROOT/'pilot_compiled'/'audit.json').read_text())
    pars=np.array(audit['state_global_and_Powell']['parameters'])
    cases=[(100.,.0015,.008,None),(300.,.0015,.008,None),(600.,.002,-.012,None),
           (300.,.0008,.004,None),(750.,.0008,.004,None),
           (100.,.0015,.008,1.5e-5),(100.,.0015,-.012,1.5e-5)]
    out=ROOT/('rollout_rival_refined' if args.refine else 'rollout_rival');out.mkdir(exist_ok=True)
    dump(out/'config.json',dict(base_parameters=pars.tolist(),cases=cases,continuous_current_action_refinement=args.refine,
        target='Predicting one action plus a cheap feedback tail repairs coverage failure?',
        scope='Existing rollout principle; not NNBO, ENE or NEOC reproduction.'))
    results={}
    for j,(theta,e0,v0,reserve) in enumerate(cases):
        p=Problem(theta,initial_error=e0,initial_velocity=v0,
                    overrides={} if reserve is None else dict(energy_reserve=reserve))
        z,timing=policy(p,pars,args.refine);m=p.metrics(z)
        base=p.metrics(p.policy(pars,'state'));p.evaluate(z)
        roots=[]
        for start in [z,np.zeros(2*p.n)]:
            res,elapsed=p.solve(start,maxiter=220)
            roots.append(dict(success=bool(res.success),seconds=elapsed,metrics=p.metrics(res.x)))
        accepted=[r for r in roots if r['success'] and r['metrics']['minimum_margin']>=-1e-6]
        if not accepted:raise RuntimeError('No accepted independent optimization')
        best=min(accepted,key=lambda r:r['metrics']['cost'])
        results[str(j)]=dict(theta=theta,e0=e0,v0=v0,reserve=reserve,
              rollout=m,base=base,optimization=best,optimization_runs=roots,
              rollout_excess=m['cost']/best['metrics']['cost']-1,
              improvement_over_base=1-m['cost']/base['cost'],timing=timing)
        dump(out/'results.json',results)
        np.save(out/f'case{j}_rates.npy',z)
        np.savetxt(out/f'case{j}_trace.csv',p.rollout(z,True)[2],delimiter=',',header='t,e,v,K,D,E,uK,uD',comments='')
        print(j,'excess',results[str(j)]['rollout_excess'],'worst_seconds',timing['worst_decision_seconds'],flush=True)


if __name__=='__main__':main()
