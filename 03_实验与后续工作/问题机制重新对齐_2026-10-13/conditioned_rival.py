"""Material-conditioned cheap state feedback; offline preparation is charged.

Blend targets evaluated at the SAME actual state, never blend separate rollouts.
Anchors use only the common initial disturbance case, not held-out cases.
"""
import json,time
import numpy as np
import compiled_flow
from dynamic_tradeoff import Problem,dump,ROOT


def target(y,pars,p):
    ak,ad,se,sv,bias,tau=pars
    gate=(y[0]/se)**2+(y[1]/sv)**2;gate=gate/(1+gate)
    sign=np.tanh(bias*y[0]*y[1]/(se*sv))
    return p.pref+np.array([ak,ad])*np.clip(gate*(1+sign)/2,0,1),tau


def policy(p,thetas,parsets):
    # A continuous piecewise-linear material lookup, all anchors precomputed.
    i=np.clip(np.searchsorted(thetas,p.theta)-1,0,len(thetas)-2)
    w=np.clip((p.theta-thetas[i])/(thetas[i+1]-thetas[i]),0.,1.)
    y=p.y0.copy();rates=[];tic=time.perf_counter();decision_seconds=0.
    for j in range(p.n):
        start=time.perf_counter()
        a,taua=target(y,parsets[i],p);b,taub=target(y,parsets[i+1],p)
        desired=(1-w)*a+w*b;tau=(1-w)*taua+w*taub
        desired=np.clip(desired,p.lo,p.hi)
        z=np.clip((desired-y[2:4])/max(tau,p.h)/p.rate,-1,1)
        z=np.clip(z,(p.lo-y[2:4])/(p.h*p.rate),(p.hi-y[2:4])/(p.h*p.rate))
        reserve=max(y[4]-p.cfg['energy_min'],0.)
        need=.5*y[0]**2*max(z[0]*p.rate[0],0.)*p.h
        if need>reserve and need>0:z[0]*=reserve/need
        decision_seconds+=time.perf_counter()-start
        rates.append(z)
        for _ in range(p.sub):y=compiled_flow.step(y,z,p.theta,p.params,p.dt)[0]
    return np.array(rates).ravel(),decision_seconds


def main():
    pilot=json.loads((ROOT/'pilot_compiled'/'audit.json').read_text())
    family=json.loads((ROOT/'material_and_state_holdout'/'results.json').read_text())
    thetas=np.array([100.,300.,600.,900.])
    parsets=np.array([family['0']['retuned_state']['parameters'],
                     pilot['state_global_and_Powell']['parameters'],
                     family['1']['retuned_state']['parameters'],family['2']['retuned_state']['parameters']])
    cases=[(200.,.0015,.008),(450.,.0015,.008),(750.,.0015,.008),
           (300.,.0008,.004),(300.,.0025,-.008),(600.,.002,-.012),
           (200.,.002,-.012),(750.,.0008,.004)]
    folder=ROOT/'conditioned_rival';folder.mkdir(exist_ok=True)
    dump(folder/'config.json',dict(thetas=thetas.tolist(),parsets=parsets.tolist(),
       cases=cases,scope='Offline material-conditioned state controller; lookup is existing mechanism, not new contribution.',
       preparation='Four locally tuned anchor policies; tuning work is offline cost, not free or a zero-shot method.'))
    results={}
    for j,(th,e0,v0) in enumerate(cases):
        p=Problem(th,initial_error=e0,initial_velocity=v0);p.evaluate(np.zeros(2*p.n))
        z,seconds=policy(p,thetas,parsets);m=p.metrics(z)
        runs=[]
        for start in [z,np.zeros(2*p.n)]:
            res,elapsed=p.solve(start);mr=p.metrics(res.x)
            runs.append(dict(success=bool(res.success),seconds=elapsed,metrics=mr))
        accepted=[r for r in runs if r['success'] and r['metrics']['minimum_margin']>=-1e-6]
        if not accepted:raise RuntimeError('No accepted optimizer '+str(j))
        best=min(accepted,key=lambda r:r['metrics']['cost'])
        timings=[]
        for _ in range(5):
            zr,ts=policy(p,thetas,parsets);timings.append(ts)
            assert np.array_equal(zr,z)
        results[str(j)]=dict(theta=th,e0=e0,v0=v0,lookup_metrics=m,
               optimal=best,optimizer_runs=runs,
               lookup_excess=m['cost']/best['metrics']['cost']-1,
               controller_seconds_per_episode_median=float(np.median(timings)),
               controller_seconds_per_decision_median=float(np.median(timings)/p.n),
               evaluator_and_physical_simulation_excluded=True)
        dump(folder/'results.json',results)
        np.savetxt(folder/f'case{j}_trace.csv',p.rollout(z,True)[2],delimiter=',',
            header='t,e,v,K,D,E,uK,uD',comments='')
        print(j,th,'excess',results[str(j)]['lookup_excess'],flush=True)


if __name__=='__main__':main()
