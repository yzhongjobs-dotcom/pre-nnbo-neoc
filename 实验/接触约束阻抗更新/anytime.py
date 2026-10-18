"""Strong online comparator: execute the simple controller while refining.

Known-model predicted release state; no waiting with obsolete gains. A fixed
20ms release deadline and post-filter predicted task-cost comparison retain the
simple target on deadline miss. This is mature rollout/anytime control, not RTI
certification or a new theoretical contribution.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize
from scipy.interpolate import RBFInterpolator
import model as m
from screen import library,critical
from verify import rhs
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
DEADLINE=.020


@njit(cache=True)
def window(targets,state,param,duration):
    y=state.copy();cost=0.
    for phase in range(2):
        length=duration/2
        for i in range(int(np.ceil(length/m.DT-1e-12))):
            dt=min(m.DT,length-i*m.DT);target=targets[phase]
            a,ja,_,_,_,_=m.flow(y,target,param)
            b,jb,_,_,_,_=m.flow(y+dt*a/2,target,param)
            c,jc,_,_,_,_=m.flow(y+dt*b/2,target,param)
            d,jd,_,_,_,_=m.flow(y+dt*c,target,param)
            y+=dt*(a+2*b+2*c+d)/6
            cost+=dt*(ja+2*jb+2*jc+jd)/6
    return cost,y


@njit(cache=True)
def scores(plans,state,param,duration):
    out=np.empty(len(plans))
    for i in range(len(plans)):out[i]=window(plans[i],state,param,duration)[0]
    return out


def independent(plan,state,param,simple):
    value=np.r_[state,0.]
    switches=[0.,DEADLINE,DEADLINE+(.3-DEADLINE)/2,.3]
    for phase,(start,end) in enumerate(zip(switches[:-1],switches[1:])):
        target=simple if phase==0 else plan[phase-1]
        def f(t,y):
            dy,cost,*_=rhs(y[:9],target,param)
            return np.r_[dy,cost]
        sol=solve_ivp(f,(start,end),value,method='DOP853',rtol=2e-10,atol=2e-12,max_step=.001)
        assert sol.success
        value=sol.y[:,-1]
    return float(value[-1])


def main():
    raw=np.load(HERE/'ranking_paired.npz');base=json.loads((HERE/'results.json').read_text())
    fixed=np.array(base['fixed_target'])/m.SCALE;theta=np.array(base['critical_theta'])
    fit=np.load(HERE/'ranker_weights.npz')
    ranker=RBFInterpolator(fit['rbf_x'],fit['rbf_y'],kernel='thin_plate_spline',smoothing=float(fit['rbf_smoothing']))
    # Same 96 contexts: post-result strengthening, not new independent evidence.
    params,states=raw['params'],raw['states'];rows=[];arrays={};allplans={}
    for i,(p,s) in enumerate(zip(params,states)):
        simple=critical(p[None],theta)[0];start=np.tile(simple,(2,1))
        window(start,s,p,DEADLINE);scores(np.array([start]),s,p,.3-DEADLINE)
        for name in np.random.default_rng(7900+i).permutation(['library_with_immediate_simple','four_step_local_with_immediate_simple',
                     'RBF_warm_four_step_with_immediate_simple']):
            tic=time.perf_counter()
            # Includes future-state preparation and evaluation of the fallback.
            prefix,release=window(start,s,p,DEADLINE)
            fallback,_=window(start,release,p,.3-DEADLINE)
            best=start.copy();bestj=fallback;solmeta=None
            if name.startswith('library'):
                plans=library(release,p,fixed,theta)
                values=scores(plans,release,p,.3-DEADLINE)
                j=int(np.argmin(values))
                if values[j]<bestj:best=plans[j];bestj=values[j]
            else:
                guess=start
                if name.startswith('RBF'):
                    feature=(m.features(p[None],s[None])-fit['mean'])/fit['std']
                    guess=library(s,p,fixed,theta)[int(np.argmin(ranker(feature)))]
                    value,_=window(guess,release,p,.3-DEADLINE)
                    if value<bestj:best,bestj=guess.copy(),value
                class BudgetReached(Exception):pass
                last_evaluation=.0005;evaluations=0
                def objective(z):
                    nonlocal best,bestj
                    nonlocal last_evaluation,evaluations
                    # Reserve numerical-evaluation/return overhead; record any
                    # remaining misses rather than calling this a WCET guarantee.
                    reserve=max(.001,2*last_evaluation)
                    if time.perf_counter()-tic>=DEADLINE-reserve:raise BudgetReached
                    before=time.perf_counter()
                    target=m.normalized_targets(z)
                    j,_=window(target,release,p,.3-DEADLINE)
                    last_evaluation=time.perf_counter()-before;evaluations+=1
                    if j<bestj:best,bestj=target.copy(),j
                    return j
                try:
                    sol=minimize(objective,m.target_coordinates(best).ravel(),method='L-BFGS-B',
                        bounds=[(0,1)]*8,options=dict(maxiter=4,maxls=10,ftol=1e-10,eps=2e-5))
                    solmeta=dict(iterations=int(sol.nit),evaluations=int(sol.nfev),success=bool(sol.success),message=str(sol.message))
                except BudgetReached:
                    solmeta=dict(iterations=None,evaluations=evaluations,success=False,message='Planning budget reached; keep best evaluated feasible candidate')
            elapsed=time.perf_counter()-tic
            missed=elapsed>DEADLINE
            chosen=start if missed else best
            value=prefix+window(chosen,release,p,.3-DEADLINE)[0]
            allplans[(name,i)]=chosen
            rows.append(dict(case=i,method=name,J=float(value),
                improvement_over_immediate_simple_pct=float((1-value/raw['critical_J'][i])*100),
                worker_ms=elapsed*1000,missed_deadline=missed,
                predicted_prepublication_improvement_pct=float((1-(prefix+bestj)/raw['critical_J'][i])*100),
                solver=solmeta))
        if i%24==0:print('anytime',i,flush=True)
    summaries=[];checks=[]
    for name in sorted(set(r['method'] for r in rows)):
        chosen=[r for r in rows if r['method']==name]
        summaries.append(dict(method=name,mean_improvement_over_immediate_simple_pct=float(np.mean([r['improvement_over_immediate_simple_pct'] for r in chosen])),
            missed_deadlines=sum(r['missed_deadline'] for r in chosen),
            worker_median_ms=float(np.median([r['worker_ms'] for r in chosen])),
            worker_p95_ms=float(np.quantile([r['worker_ms'] for r in chosen],.95))))
        arrays[name+'_targets']=np.array([allplans[(name,i)] for i in range(len(states))])
        arrays[name+'_J']=np.array([r['J'] for r in chosen])
        for case in [0,1]:
            p,s=params[case],states[case];simple=critical(p[None],theta)[0]
            check=independent(allplans[(name,case)],s,p,simple)
            checks.append(dict(method=name,case=case,relative_cost_error=abs(chosen[case]['J']/check-1)))
    output=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
        protocol=dict(cases=96,source='Same seed7503 contexts, exploratory baseline strengthening after looking at earlier results',
            release_deadline_ms=20,execution='Immediate tuned critical-damping target, known-model predicted release state, candidate check, retain fallback on deadline miss',
            scope='One model update, exact model/state, no rolling MPC or unknown-parameter guarantee; local solver is not a fully tuned RTI baseline',
            RBF_cost='Actual RBF feature/target generation inside the measured worker interval'),
        methods=summaries,independent_checks=checks,cases=rows)
    assert max(r['relative_cost_error'] for r in checks)<.005
    assert min(r['improvement_over_immediate_simple_pct'] for r in rows)>-1e-8
    np.savez_compressed(HERE/'anytime_paired.npz',**arrays)
    (HERE/'anytime_results.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
