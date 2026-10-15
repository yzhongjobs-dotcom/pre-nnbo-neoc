"""Repair the carrier, then perform actual NEOC and strongest-baseline probes.

This finite-horizon normal-impedance projection is a diagnostic of the idea,
not an exact discretization-independent implementation of every R15 theorem.
"""
from pathlib import Path
import csv
import gc
import hashlib
import json
import time
import numpy as np
from scipy.optimize import least_squares
from bellman_carrier import Carrier,rollout,src,DREF,FD,T,RATE
from independent_reference import ShootingReference

HERE=Path(__file__).resolve().parent
ETA=np.array([1.5,0.])

def save_csv(name,data,head):
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(data)

def array_hash(a):return hashlib.sha256(np.ascontiguousarray(a).view(np.uint8)).hexdigest()

def main():
    started=time.perf_counter()
    cfg=dict(date='2026-10-13',started_date='2026-10-12',
             scene='30->45 kPa, cylinder scan -6->6 mm at 8 mm/s',
             state=['indentation','normal_velocity','actual_normal_stiffness'],
             damping_fixed=8.,energy_initial=.02,energy_minimum=.001,
             model='same cylinder Hertz spatial kernel and inertial impedance equation',
             objective='finite 1.5 s, same normalized state/force/gain/rate costs as prior task',
             critic='fixed nodal value basis with exact descriptor-dependent terminal cost',
             anchor_descriptor=[1.,0.],test_descriptor=ETA.tolist(),
             sensitivity_epsilon=1e-5,sensitivity_type='one-sided soft->hard direction',
             trials='one physical soft->hard task; small-change check and solver diagnostics separate',
             constraint_check='sampled physical rollout, not full certified PublishOK',
             gates=dict(reference_cost_tolerance=.05,NEOC_cost_extra=.05,
                        force_rmse_extra_fraction=.02,latency_reduction=.3),
             novelty_scope='does not infer journal readiness or invalidate all parameter families')
    result=dict(status='RUNNING',config=cfg,source_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                for name in ['bellman_carrier.py','independent_reference.py','run_repaired_validation.py']},
                failures_preserved='previous stationary-critic and 5-state failed files unchanged',
                stages={},timing_scope='prototype same-host timings, no hard real-time claim')
    (HERE/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')

    print('Independent constrained shooting reference...',flush=True)
    shoot=ShootingReference(tuple(ETA))
    z=np.full(shoot.M,150.);direction=np.sin(np.arange(shoot.M)+.3)
    _,grad=shoot.fun(z);step=1e-3
    fd=(shoot.fun(z+step*direction)[0]-shoot.fun(z-step*direction)[0])/(2*step)
    gradient_error=abs(fd-grad@direction)/max(abs(fd),1e-10)
    assert gradient_error<1e-4
    fit,trace,reference=shoot.solve()
    assert reference['success'] and reference['maximum_rate']<=RATE+1e-8
    reference['gradient_direction_relative_error']=float(gradient_error)
    save_csv('independent_shooting_trace.csv',trace,['time_s','d_m','v_m_s','gain','rate','force_N','stage'])
    save_csv('independent_shooting_knots.csv',np.c_[shoot.times[1:],fit.x],['time_s','gain'])
    result['independent_reference']=reference

    # Preserved failed resolution, followed by physical-resolution repair.
    resolution=[]
    for nd,nv,steps,focused in [(17,17,100,False),(33,33,200,True),(49,49,200,True)]:
        c=Carrier(nd,nv,steps,focused=focused)
        begin=time.perf_counter();V,_=c.solve(tuple(ETA));solve_s=time.perf_counter()-begin
        trace,metric=rollout(c,V,tuple(ETA))
        entry=dict(nd=nd,nv=nv,steps=steps,focused=focused,shape=c.shape,
                   cost=metric['cost'],cost_gap_to_shooting=metric['cost']/reference['cost']-1,
                   constraints=metric['sampled_constraints_pass'],solve_seconds=solve_s)
        resolution.append(entry);print('Resolution check:',json.dumps(entry),flush=True)
        if focused and nd==49:
            full_reference=metric;full_reference_trace=trace
            chosen=c
        else:del c,V;gc.collect()
    result['resolution_checks']=resolution
    assert full_reference['sampled_constraints_pass']
    assert abs(full_reference['cost']/reference['cost']-1)<.05
    save_csv('full_root_oracle_trace.csv',full_reference_trace,
             ['time_s','d_m','v_m_s','gain','energy','force_N','rate','stage','accumulated_cost'])
    del V;gc.collect()
    c=chosen
    # Release the future model after quality checks so query preparation is measured.
    c.model_cache.clear();gc.collect()
    print('Preparing uniform anchor and directional root sensitivity...',flush=True)
    begin=time.perf_counter();anchor,P=c.solve((1.,0.),True);anchor_s=time.perf_counter()-begin
    anchor_res=float(np.max(abs(c.residual(anchor,(1.,0.)))))
    assert anchor_res<1e-10
    begin=time.perf_counter();Cw,dphi=c.sensitivity(anchor,P,(1.,0.),eps=1e-5);sensitivity_s=time.perf_counter()-begin
    pert=(1.+1e-5,0.)
    exact_pert,_=c.solve(pert)
    phi0=c.model((1.,0.))[1];phip=c.model(pert)[1]
    derivative=(exact_pert-anchor-(phip-phi0)[None,:])/1e-5
    derivative_error=float(np.linalg.norm(Cw-derivative)/max(np.linalg.norm(derivative),1e-12))
    print('Directional derivative check:',derivative_error,flush=True)
    result['stages']['shared_root']=dict(status='PASS_DISCRETE_ROOT_AND_PHYSICAL_REFERENCE',
        root_shape=anchor.shape,root_residual=anchor_res,anchor_seconds=anchor_s,
        sensitivity_seconds=sensitivity_s,directional_derivative_relative_error=derivative_error,
        anchor_sha256=array_hash(anchor),sensitivity_sha256=array_hash(Cw),
        continuous_HJB_optimality_claimed=False)
    del exact_pert,derivative
    c.model_cache.pop(pert,None);gc.collect()

    # A common causal wrench history. Truth is used only to generate measurements.
    print('Collecting common history and estimating the front material/boundary...',flush=True)
    history_trace,_=rollout(c,anchor,(1.,0.),true_eta=tuple(ETA))
    discrepancy=[]
    for row in history_trace:
        discrepancy.append(row[5]-src.fast_force(row[0],row[1],(1.,0.)))
    detect=int(np.flatnonzero(np.array(discrepancy)>.005*FD)[0])
    indices=detect+5*np.arange(12)
    observed=history_trace[indices];tt=observed[:,0];dd=observed[:,1]
    force,moment=src.contact(tt,dd,ETA,True)
    # Same measured force and moment history to both methods, no truth in fit.
    def residual(eta):
        ff,mm=src.contact(tt,dd,eta,True)
        return np.r_[(ff-force)/.003,(mm-moment)/2e-6]
    estimate=least_squares(residual,[1.1,.05],bounds=([.8,-.5],[2.,.5]),
                           ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100)
    assert estimate.success and np.linalg.norm(residual(estimate.x))<1e-6
    eta=tuple(estimate.x);event=float(tt[-1])
    estimate_record=dict(eta=list(eta),event_time=event,estimated_boundary_m=eta[1]*.003,
                         estimated_front_E_Pa=eta[0]*30000,nominal_lead_time_s=(eta[1]*.003+.006-.008*event)/.008,
                         fit_residual=float(np.linalg.norm(residual(estimate.x))),
                         samples=12,noise_added=False,geometry_surface_and_tangential_parameters_known=True,
                         full_MHE=False)
    result['common_estimate']=estimate_record
    save_csv('common_wrench_history.csv',np.c_[tt,dd,force,moment],
             ['time_s','indentation_m','normal_force_N','raw_pitch_moment_Nm'])
    print('Common estimate:',json.dumps(estimate_record),flush=True)

    begin=time.perf_counter();c.model(eta);model_s=time.perf_counter()-begin
    rows=[]
    methods=['NEOC','cached_chord','one_Newton','full_backward_solve']
    for method in methods:
        print('Testing:',method,flush=True)
        begin=time.perf_counter();candidate=c.candidate(method,eta,anchor,P,Cw,dphi)
        generation=time.perf_counter()-begin
        begin=time.perf_counter();F=c.residual(candidate,eta);res_s=time.perf_counter()-begin
        begin=time.perf_counter();trace,metric=rollout(c,candidate,eta,true_eta=tuple(ETA),event_time=event,oldV=anchor)
        rollout_s=time.perf_counter()-begin
        # Monotone interpolation: exact projected-grid root error <= sum of
        # per-time residual infinity norms. It is not a physical-ODE certificate.
        E=float(np.sum(np.max(abs(F),axis=1)))
        row=dict(method=method,metrics=metric,generation_seconds=generation,
                 shared_model_prepare_seconds=model_s,global_residual_check_seconds=res_s,
                 sampled_rollout_check_seconds=rollout_s,
                 with_standalone_checks_seconds=model_s+generation+res_s+rollout_s,
                 fused_exact_baseline_seconds=(model_s+generation+rollout_s if method=='full_backward_solve' else None),
                 maximum_discrete_residual=float(np.max(abs(F))),discrete_root_error_bound=E,
                 discrete_MDP_policy_gap_bound=2*E,
                 terminal_condition_error=float(np.max(abs(F[-1]))),
                 candidate_sha256=array_hash(candidate),NEOC_online_Newton_calls=0 if method=='NEOC' else None)
        rows.append(row)
        save_csv(f'trace_{method}.csv',trace,
                 ['time_s','d_m','v_m_s','gain','energy','force_N','rate','stage','accumulated_cost'])
        del F,candidate;gc.collect()
        print(json.dumps(row),flush=True)

    trace,fixed=rollout(c,anchor,(1.,0.),true_eta=tuple(ETA),fixed_gain=True)
    save_csv('trace_fixed_gain.csv',trace,
             ['time_s','d_m','v_m_s','gain','energy','force_N','rate','stage','accumulated_cost'])
    reference_row=next(r for r in rows if r['method']=='full_backward_solve')
    Jref=reference_row['metrics']['cost'];neoc=rows[0]
    for row in rows:
        row['relative_cost_extra_to_full']=row['metrics']['cost']/Jref-1
        row['force_rmse_extra_fraction']=(row['metrics']['force_rmse_N']-reference_row['metrics']['force_rmse_N'])/FD
    neoc_accuracy=(neoc['metrics']['sampled_constraints_pass'] and
                   neoc['relative_cost_extra_to_full']<=.05 and neoc['force_rmse_extra_fraction']<=.02)
    latency_reduction=1-neoc['with_standalone_checks_seconds']/reference_row['fused_exact_baseline_seconds']
    # Extra sensitivity preparation is charged, not hidden in the common anchor.
    neoc_total=sensitivity_s+neoc['with_standalone_checks_seconds']
    full_total=reference_row['fused_exact_baseline_seconds']
    benefit=(neoc_accuracy and latency_reduction>=.3 and neoc_total<full_total)
    result['stages']['parameter_update']=dict(results=rows,fixed_gain_metrics=fixed,
         NEOC_quality_pass=neoc_accuracy,latency_reduction_vs_fused_exact=latency_reduction,
         incremental_sensitivity_prepare_seconds=sensitivity_s,
         NEOC_single_event_total_seconds=neoc_total,full_single_event_total_seconds=full_total,
         fixed_cost_extra_to_full=fixed['cost']/Jref-1,
         fixed_cost_extra_to_independent_reference=fixed['cost']/reference['cost']-1,
         continue_to_preview_control=benefit)

    # A local check distinguishes the small-change premise from the 50% test.
    print('Checking the local 5% change premise...',flush=True)
    local_eta=(1.05,0.)
    local_neoc=c.candidate('NEOC',local_eta,anchor,P,Cw,dphi)
    _,lm=rollout(c,local_neoc,local_eta,true_eta=local_eta,event_time=event,oldV=anchor)
    del local_neoc;gc.collect()
    local_full=c.solve(local_eta)[0]
    _,lf=rollout(c,local_full,local_eta,true_eta=local_eta,event_time=event,oldV=anchor)
    _,lfix=rollout(c,anchor,(1.,0.),true_eta=local_eta,fixed_gain=True)
    result['local_change_check']=dict(descriptor=local_eta,NEOC=lm,full=lf,fixed=lfix,
                                     relative_NEOC_cost_extra=lm['cost']/lf['cost']-1)
    del local_full;gc.collect()
    result['stages']['preview_control']=dict(
        status='NOT_RUN_NO_ESTABLISHED_UPDATE_BENEFIT' if not benefit else 'ELIGIBLE_PENDING_CONTROL_COMPARISON',
        sensing_executed=True,spatial_preview_control_comparison_executed=False,
        note='Early wrench estimation is measured, not a claim of downstream preview-control benefit.')
    maxe=max(abs(c.d-DREF))
    result['energy_inactive_envelope_bound']=dict(
        lower_bound=.02-.5*maxe*maxe*RATE*T,required=.001,
        reason='Ignore nonnegative damping recovery; bound maximal stiffness injection while within grid envelope.')
    result['status']='CORE_BENEFIT_NOT_ESTABLISHED_AFTER_CARRIER_REPAIR' if not benefit else 'UPDATE_BENEFIT_FOUND'
    result['duration_seconds']=time.perf_counter()-started
    (HERE/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('Final status:',result['status'],flush=True)

if __name__=='__main__':main()
