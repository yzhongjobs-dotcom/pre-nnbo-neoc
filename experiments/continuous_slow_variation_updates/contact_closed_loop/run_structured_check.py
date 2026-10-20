from pathlib import Path
import hashlib
import json
import time
import numpy as np
import run_contact_validation as trial
from continuous_contact import UniformCarrier
from sensitivity_fast import fast_sensitivity

BASE=Path(__file__).resolve().parent
OUT=BASE/'01_reused_action_sensitivity'


def main():
    OUT.mkdir(exist_ok=True);c=UniformCarrier();checks=[]
    original=UniformCarrier.sensitivity
    for eta in [(1.,0.),(1.1,0.)]:
        V,P=c.solve(eta,True)
        begin=time.perf_counter();old,_=original(c,V,P,eta,eps=1e-5);slow_s=time.perf_counter()-begin
        begin=time.perf_counter();fast,dphi=fast_sensitivity(c,V,P,eta,eps=1e-5);fast_s=time.perf_counter()-begin
        pert=(eta[0]+1e-5,0.);exact,_=c.solve(pert)
        reference=(exact-V)/1e-5-dphi[None,:]
        diff=float(np.linalg.norm(fast-old)/max(np.linalg.norm(old),1e-12))
        gap=float(np.linalg.norm(fast-reference)/max(np.linalg.norm(reference),1e-12))
        checks.append(dict(eta=eta,fast_to_original_relative_error=diff,
            fast_to_complete_perturbed_root_relative_error=gap,
            original_derivative_seconds=slow_s,structured_derivative_seconds=fast_s))
        print(json.dumps(checks[-1]),flush=True)
        assert diff<.02 and gap<.02
        c.model_cache.clear()
    (OUT/'derivative_audit.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    UniformCarrier.sensitivity=fast_sensitivity
    trial.HERE=OUT
    result=dict(sensitivity='fixed-action envelope derivative; finite eta difference of transition',
        candidate_acceptance_checks=0,derivative_audit=checks,scenarios=[],status='RUNNING')
    for kind,duration,period in [('ramp',3.,.25),('ramp',6.,.25),('ramp',6.,.5),('cycle',6.,.25)]:
        scores=[trial.simulate(c,m,kind,duration,period) for m in trial.METHODS]
        ref=next(r for r in scores if r['method']=='full_solve')
        for row in scores:
            row['cost_extra_fraction']=row['cost']/ref['cost']-1
            row['force_rmse_extra_fraction']=row['force_rmse_N']/ref['force_rmse_N']-1
            row['quality_gate_pass']=bool(row['cost_extra_fraction']<=.01 and row['force_rmse_extra_fraction']<=.05 and row['sampled_constraints_pass'])
        eligible=[r for r in scores if r['method'] in ('cached_chord','one_Newton','full_solve') and r['quality_gate_pass']]
        strongest=min(eligible,key=lambda r:r['total_controller_seconds'])
        for row in scores:
            row['strongest_eligible_optimizer']=strongest['method']
            row['total_strategy_time_reduction']=1-row['total_controller_seconds']/strongest['total_controller_seconds']
            row['screen_pass']=bool(row['quality_gate_pass'] and row['total_strategy_time_reduction']>=.3)
        result['scenarios'].append(dict(path=kind,duration_s=duration,period_s=period,scores=scores))
        (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
        print('COMPLETE',json.dumps(scores),flush=True)
    result['status']='COMPLETED_DIAGNOSTIC_NOT_CERTIFIED'
    result['source_hashes']={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                           for name in ['continuous_contact.py','run_contact_validation.py','sensitivity_fast.py','run_structured_check.py']}
    (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
