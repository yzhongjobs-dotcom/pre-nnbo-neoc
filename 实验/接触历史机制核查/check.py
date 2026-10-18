"""Illustrative passive-contact memory diagnostic, not a prior replication or controller benchmark."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parent
DT=0.0005
TAUS=np.array([0.05,0.5])
KMOD=np.array([500.,500.])
KINF=500.
DELTA=0.004
FD=3.
MASS=1.
HORIZON=0.3

def past(age):
    """Smooth monotone loading ending age seconds before now; same current displacement and force."""
    t=np.arange(-0.4,DT/2,DT)
    duration=.02
    z=np.clip((t+age+duration)/duration,0.,1.)
    ramp=(1-np.cos(np.pi*z))/2
    v=np.where((z>0)&(z<1),np.pi*np.sin(np.pi*z)/(2*duration),0.)
    # Unit-amplitude loading: solve z'_i = k_i v - z_i/tau_i.
    a=np.exp(-DT/TAUS)
    unit=np.zeros((len(t),2))
    for i in range(1,len(t)):
        # Exact first-order-hold input/displacement increment on each interval.
        vel=(ramp[i]-ramp[i-1])/DT
        unit[i]=a*unit[i-1]+KMOD*TAUS*(1-a)*vel
    amplitude=(FD-KINF*DELTA)/unit[-1].sum()
    delta=DELTA-amplitude+amplitude*ramp
    return t,delta,amplitude*unit,amplitude

def rhs(y,kd):
    delta,v,z1,z2,ka,da,integ=y.T
    force=KINF*delta+z1+z2
    # The diagnostic remains in compression; no extrapolation to contact loss.
    u=FD+ka*(DELTA-delta)-da*v
    acc=(u-force)/MASS
    ratek=np.clip((kd[:,0]-ka)/.03,-10000,10000)
    rated=np.clip((kd[:,1]-da)/.03,-2000,2000)
    cost=(force-FD)**2+.1*((delta-DELTA)/.001)**2+.01*(v/.01)**2
    return np.stack([v,acc,KMOD[0]*v-z1/TAUS[0],KMOD[1]*v-z2/TAUS[1],ratek,rated,cost],axis=1)

def rollout(z,kd,dt=DT):
    y=np.tile(np.r_[DELTA,0.,z,500.,40.,0.],(len(kd),1))
    min_force=1e10;min_delta=1e10
    for _ in range(round(HORIZON/dt)):
        a=rhs(y,kd);b=rhs(y+dt*a/2,kd);c=rhs(y+dt*b/2,kd);d=rhs(y+dt*c,kd)
        y+=dt*(a+2*b+2*c+d)/6
        min_force=min(min_force,float(np.min(KINF*y[:,0]+y[:,2]+y[:,3])))
        min_delta=min(min_delta,float(y[:,0].min()))
    return y[:,-1],y,min_force,min_delta

def observer(t,delta,measured,indices):
    a=np.exp(-DT/TAUS);forced=np.zeros((len(t),2))
    for i in range(1,len(t)):
        vel=(delta[i]-delta[i-1])/DT
        forced[i]=a*forced[i-1]+KMOD*TAUS*(1-a)*vel
    design=np.exp(-(t[indices]-t[0])[:,None]/TAUS)
    residual=measured-KINF*delta[indices]-forced[indices].sum(axis=1)
    initial=np.linalg.lstsq(design,residual,rcond=None)[0]
    current=np.exp(-(t[-1]-t[0])/TAUS)*initial+forced[-1]
    return current,float(np.linalg.cond(design))

def run():
    rng=np.random.default_rng(20261018)
    histories=[past(0.),past(.15)]
    future=np.linspace(0,.3,601)
    forces=np.array([KINF*DELTA+(h[2][-1]*np.exp(-future[:,None]/TAUS)).sum(axis=1) for h in histories])
    diff=np.abs(forces[0]-forces[1]);idx=int(np.argmax(diff))
    candidates=np.array([(k,d) for k in [30.,100.,300.,600.,1000.,1800.] for d in [5.,15.,30.,60.,100.,160.]])
    costs=[];verification=[];state_info=[];observer_results=[]
    for hi,(t,delta,z,amp) in enumerate(histories):
        js,ys,mf,md=rollout(z[-1],candidates)
        costs.append(js)
        idxs=np.arange(0,len(t),round(.01/DT))
        true_force=KINF*delta[idxs]+z[idxs].sum(axis=1)
        est,condition=observer(t,delta,true_force,idxs)
        state_info.append(dict(loading_end_age_s=[0.,.15][hi],ramp_s=.02,loading_amplitude_m=float(amp),initial_penetration_m=float(delta[0]),current_penetration_m=float(delta[-1]),current_force_N=float(KINF*delta[-1]+z[-1].sum()),current_internal_forces_N=z[-1].tolist(),min_future_force_across_candidates_N=mf,min_future_penetration_across_candidates_m=md,noiseless_observer_state_error_N=float(np.max(np.abs(est-z[-1]))),observer_design_condition=condition))
        rows=[]
        for noise in [.01,.05]:
            errs=[];regrets=[];choices=[]
            for _ in range(100):
                est,_=observer(t,delta,true_force+rng.normal(0,noise,len(idxs)),idxs)
                # Rank the shared feasible candidate set on estimated internal state.
                pred,_,_,_=rollout(est,candidates,dt=.001)
                pick=int(np.argmin(pred));errs.append(float(np.sqrt(np.mean((KINF*DELTA+(est*np.exp(-future[:,None]/TAUS)).sum(axis=1)-forces[hi])**2))))
                regrets.append(float((js[pick]-js.min())/js.min()));choices.append(pick)
            rows.append(dict(force_noise_std_N=noise,noise_draws=100,forecast_rmse_mean_N=float(np.mean(errs)),forecast_rmse_p95_N=float(np.quantile(errs,.95)),candidate_regret_mean_percent=100*float(np.mean(regrets)),candidate_regret_p95_percent=100*float(np.quantile(regrets,.95))))
        observer_results.append(rows)
        chosen=int(js.argmin());single=candidates[[chosen]];y0=np.r_[DELTA,0.,z[-1],500.,40.,0.]
        sol=solve_ivp(lambda tt,yy:rhs(yy[None,:],single)[0],(0,HORIZON),y0,method='DOP853',rtol=1e-10,atol=1e-12)
        half=rollout(z[-1],single,DT/2)[0][0]
        verification.append(dict(history=hi,candidate=single[0].tolist(),J_rk4=float(js[chosen]),J_dop853=float(sol.y[-1,-1]),dop853_relative_difference=float(abs(js[chosen]-sol.y[-1,-1])/sol.y[-1,-1]),halved_step_relative_difference=float(abs(half-js[chosen])/half)))
    costs=np.array(costs);best=costs.min(axis=1);paired=costs/best[:,None]-1
    common=int(np.argmin(paired.mean(axis=0)))
    result=dict(scope='Two constructed, reachable smooth loading histories in a known passive two-mode generalized Maxwell contact; diagnostic only. Not a robot/unknown-parameter estimator experiment, prior replication, learned-policy test, safety certificate, or top-journal contribution.',parameters=dict(dt_s=DT,taus_s=TAUS.tolist(),branch_stiffness_N_per_m=KMOD.tolist(),equilibrium_stiffness_N_per_m=KINF,mass_kg=MASS,current_penetration_m=DELTA,desired_force_N=FD,horizon_s=HORIZON),history_states=state_info,prescribed_hold=dict(max_future_force_separation_N=float(diff[idx]),at_s=float(future[idx]),force_separation_at_100ms_N=float(diff[200]),instantaneous_single_valued_predictor_minimax_error_lower_bound_N=float(diff[idx]/2)),candidate_decision=dict(candidates_KD=candidates.tolist(),cost_per_history=costs.tolist(),best_KD_per_history=candidates[np.argmin(costs,axis=1)].tolist(),best_cost_per_history=best.tolist(),best_common_KD=candidates[common].tolist(),best_common_candidate_mean_relative_regret_percent=float(paired[:,common].mean()*100),best_common_candidate_regret_percent=(paired[:,common]*100).tolist(),note='Optimal only among 36 fixed-target gain transitions with actual gain dynamics, in this local unconstrained diagnostic. No energy filter or complete direct-force baseline. Equal weights for two illustrative histories; no population interpretation.'),history_observer=observer_results,independent_integration=verification,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['history_states','prescribed_hold','history_observer','independent_integration']},ensure_ascii=False,indent=2))
    print(json.dumps({k:v for k,v in result['candidate_decision'].items() if k not in ['candidates_KD','cost_per_history']},ensure_ascii=False,indent=2))
if __name__=='__main__':run()
