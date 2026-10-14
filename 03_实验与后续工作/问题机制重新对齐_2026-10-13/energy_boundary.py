"""A second discriminator: actual tank activation, not just a small number."""
import json
import numpy as np
from dynamic_tradeoff import Problem,tune,dump,ROOT
from audit_and_family import repair_state,independent_integral


def main():
    folder=ROOT/'energy_boundary';folder.mkdir(exist_ok=True)
    cases=[dict(theta=100.,reserve=1.5e-5,e0=.0015,v0=.008),
           dict(theta=300.,reserve=1.5e-5,e0=.0025,v0=.012),
           dict(theta=300.,reserve=1.5e-6,e0=.0015,v0=.008),
           dict(theta=100.,reserve=1.5e-5,e0=.0015,v0=-.012)]
    all_results={}
    for j,c in enumerate(cases):
        p=Problem(c['theta'],initial_error=c['e0'],initial_velocity=c['v0'],
                   overrides=dict(energy_reserve=c['reserve']))
        p.evaluate(np.zeros(2*p.n))
        simple,z=tune(p,'state',maxiter=100)
        refined,zr=repair_state(p,simple['parameters'],global_search=True)
        if refined['metrics']['cost']<simple['metrics']['cost']:simple,z=refined,zr
        runs=[]
        for start in [np.zeros(2*p.n),z]:
            res,seconds=p.solve(start,maxiter=250);m=p.metrics(res.x)
            runs.append(dict(success=bool(res.success),message=str(res.message),
                              seconds=seconds,iterations=int(res.nit),metrics=m,rates=res.x.tolist()))
        okay=[r for r in runs if r['success'] and r['metrics']['minimum_margin']>=-1e-6]
        if not okay:
            all_results[str(j)]=dict(config=c,state=simple,runs=runs,status='NO_ACCEPTED_OPTIMUM')
            dump(folder/'results.json',all_results);continue
        best=min(okay,key=lambda r:r['metrics']['cost']);optz=np.array(best['rates'])
        independent,_=independent_integral(p,optz)
        # Difference is only interpreted when the simple policy is feasible.
        feasible=simple['metrics']['minimum_margin']>=-1e-6
        all_results[str(j)]=dict(config=c,state=simple,optimization=best,runs=runs,
            state_sampled_feasible=feasible,
            state_excess_if_feasible=simple['metrics']['cost']/best['metrics']['cost']-1 if feasible else None,
            independent_optimal_cost=independent,
            independent_relative_difference=abs(independent-best['metrics']['cost'])/independent,
            optimal_reserve_ratio=(best['metrics']['min_energy_J']-.001)/c['reserve'])
        for label,zout in [('optimization',optz),('state',z)]:
            tr=p.rollout(zout,True)[2]
            np.savetxt(folder/f'case{j}_{label}.csv',tr,delimiter=',',header='t,e,v,K,D,E,uK,uD',comments='')
        dump(folder/'results.json',all_results)
        print(j,'opt energy reserve ratio',all_results[str(j)]['optimal_reserve_ratio'],
                 'state feasible',feasible,'state excess',all_results[str(j)]['state_excess_if_feasible'],flush=True)


if __name__=='__main__':main()
