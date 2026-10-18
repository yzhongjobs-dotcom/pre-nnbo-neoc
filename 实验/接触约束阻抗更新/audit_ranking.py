"""Independent held-out cost/ranker checks and paired uncertainty summaries."""
from pathlib import Path
import hashlib,json
import numpy as np
import torch
from scipy.interpolate import RBFInterpolator
import model as m
from screen import library
from ranking import quadratic
from verify import independent

HERE=Path(__file__).resolve().parent


def bootstrap(values,seed):
    rng=np.random.default_rng(seed)
    boot=values[rng.integers(0,len(values),(10000,len(values)))].mean(1)
    return [float(v) for v in np.quantile(boot,[.025,.975])]


def main():
    torch.set_num_threads(1)
    raw=np.load(HERE/'ranking_paired.npz');fit=np.load(HERE/'ranker_weights.npz')
    base=json.loads((HERE/'results.json').read_text())
    fixed=np.array(base['fixed_target'])/m.SCALE;theta=np.array(base['critical_theta'])
    x=(m.features(raw['params'],raw['states'])-fit['mean'])/fit['std']
    rbf=RBFInterpolator(fit['rbf_x'],fit['rbf_y'],kernel='thin_plate_spline',smoothing=float(fit['rbf_smoothing']))
    choices={'quadratic':np.argmin(quadratic(x)@fit['quadratic'],axis=1),
             'rbf_rank':np.argmin(rbf(x),axis=1)}
    for seed in [11,22,33]:
        stored=torch.load(HERE/f'cost_actor_{seed}.pt',weights_only=False)
        net=torch.nn.Sequential(torch.nn.Linear(11,64),torch.nn.Tanh(),torch.nn.Linear(64,64),
                                torch.nn.Tanh(),torch.nn.Linear(64,30)).double()
        net.load_state_dict(stored['weights']);net.eval()
        with torch.no_grad():choices[f'cost_actor_{seed}']=net(torch.as_tensor(x,dtype=torch.float64)).argmax(1).numpy()
    reload_error=0.
    for name,indices in choices.items():
        for i,(s,p) in enumerate(zip(raw['states'],raw['params'])):
            target=library(s,p,fixed,theta)[indices[i]]
            reload_error=max(reload_error,float(np.max(abs(target-raw[name+'_targets'][i]))))
    assert reload_error<1e-12
    checks=[]
    for i in [0,1]:
        for name in ['critical','library','teacher','quadratic','rbf_rank','cost_actor_11','cost_actor_22','cost_actor_33']:
            j,_,emin=independent(raw[name+'_targets'][i],raw['states'][i],raw['params'][i])
            actual=float(raw[name+'_J'][i])
            checks.append(dict(case=i,method=name,relative_J_error=abs(actual/j-1),energy_margin_J=emin))
    summaries=[]
    for name in ['library','teacher','quadratic','rbf_rank','cost_actor_11','cost_actor_22','cost_actor_33']:
        improvement=(1-raw[name+'_delayed_J']/raw['critical_delayed_J'])*100
        vs_rbf=(1-raw[name+'_delayed_J']/raw['rbf_rank_delayed_J'])*100
        summaries.append(dict(method=name,mean_improvement_over_critical_pct=float(improvement.mean()),
            paired_bootstrap95_pct=bootstrap(improvement,7701),
            mean_improvement_over_RBF_pct=float(vs_rbf.mean()),
            paired_bootstrap95_over_RBF_pct=bootstrap(vs_rbf,7702)))
    # Algebraic rate-saturation counterexample, not a new sampled performance case.
    p=np.array([700.,.3]);s=np.r_[.001,.001,.01,-.01,np.array([300.,100.,20.,20.])/m.SCALE,.005]
    target1=s[4:8]+m.RATE*(m.HORIZON+.025+.01)
    target2=s[4:8]+m.RATE*(m.HORIZON+.025+.10)
    first=m.rollout(np.tile(target1,(2,1)),s,p)
    second=m.rollout(np.tile(target2,(2,1)),s,p)
    equivalent=dict(physical_target_difference=(m.SCALE*(target2-target1)).tolist(),
        J_difference=float(abs(first[0]-second[0])),maximum_final_state_difference=float(np.max(abs(first[2]-second[2]))),
        tank_active_duration_s=[float(first[1][5]),float(second[1][5])],
        interpretation='When the energy constraint is inactive and target-g0 >= RATE*(T+tau), distinct targets produce the same saturated actual gain rates throughout T. Gain-label error is not policy regret.')
    assert equivalent['J_difference']<1e-12
    output=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        ranker_reload_maximum_target_error=reload_error,independent_adaptive_ODE=checks,
        max_adaptive_relative_cost_error=max(v['relative_J_error'] for v in checks),
        paired_context_bootstrap=summaries,target_nonuniqueness_example=equivalent,
        scope='Context resampling only, not new material/subject or real-robot uncertainty; unmodeled plant discrepancy excluded')
    assert output['max_adaptive_relative_cost_error']<.005
    (HERE/'ranking_checks.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
