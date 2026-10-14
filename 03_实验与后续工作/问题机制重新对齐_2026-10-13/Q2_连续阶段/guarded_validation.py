"""Repair sampled CBF undershoot and verify actual continuous equations."""
import json,numpy as np
from contact_stages import Problem,ROOT,dump
from independent_audit import independent


def main():
    p=Problem();out=ROOT/'guarded_validation';out.mkdir(exist_ok=True)
    previous=json.loads((ROOT/'strong_rivals'/'results.json').read_text())
    runs=[]
    for seed in [np.array(previous['independent_shooting_from_rollout']['rates']),
                 np.array(previous['rollout']['rates'])]:
        res,seconds=p.solve(seed);m=p.metrics(res.x);ind=independent(p,res.x)
        states=ind.pop('segment_states')
        runs.append(dict(success=bool(res.success),message=str(res.message),iterations=int(res.nit),seconds=seconds,
                  metrics=m,independent_integral=ind,
                  independent_relative_cost_error=abs(ind['cost']-m['cost'])/m['cost'],rates=res.x.tolist()))
        np.save(out/f'opt{len(runs)-1}_rates.npy',res.x)
        np.savetxt(out/f'opt{len(runs)-1}_trace.csv',p.run(res.x,True)[2],delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
        dump(out/'results.json',dict(config=p.cfg,optimization_runs=runs))
        print('root',len(runs)-1,m,ind,flush=True)
    best=min([r for r in runs if r['success'] and r['independent_integral']['min_power_W']>=-1e-12],key=lambda r:r['metrics']['cost'])
    rows=[]
    for label,kind,pars in [('fixed','fixed',previous['best_fixed']['gains']),
        ('preview','preview',[previous['best_preview']['zeta'],previous['best_preview']['lead']])]:
        z,seconds=p.policy(p.targets(kind,pars));m=p.metrics(z);ind=independent(p,z);ind.pop('segment_states')
        rows.append(dict(label=label,parameters=pars,metrics=m,seconds=seconds,independent_integral=ind,
                        excess=m['cost']/best['metrics']['cost']-1))
    projected=json.loads((ROOT/'projected_rival'/'results.json').read_text())
    independent_projected=json.loads((ROOT/'audit.json').read_text())['projected_adaptive_integral']
    result=dict(config=p.cfg,optimization_runs=runs,best_optimization=best,rivals=rows,
            projected_rollout=projected['rollout'],projected_independent_integral=independent_projected,
            projected_excess=projected['rollout']['metrics']['cost']/best['metrics']['cost']-1,
            repaired_violation='Original held-rate sampled optimum had an independently observed CBF undershoot ~1.18e-8 W. Numerical margin added; old results retained.',
            scope='Numerical continuous-equation validation; sampled minimum, not a rigorous interval proof. No claim of causal realtime NEOC benefit yet.')
    dump(out/'results.json',result);print('projected excess',result['projected_excess'],flush=True)


if __name__=='__main__':main()
