"""Two scientific conditions; excitation/model ablations, no plant changes.

S2: local-material transitions; T3m: movement/material confounding. One reused
diagnosis episode each, no tuning/final cohort. Higher-amplitude comparison is
explicit: do not attribute a larger input budget to the estimator.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from concurrent.futures import ProcessPoolExecutor,as_completed
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0,str(HERE.parent/'method_redesign'))
import common,cohort_eval as CE
from v1_servo_id import ServoID
from local_contact import optimize_probe,sine_information
from joint_policy import JointPolicy


class MeasuredOnly:
    direct=True
    def __init__(self,inner):self.inner=inner
    def reset(self,p):self.inner.reset(dict(common.P.PLATFORMS['PC125']))
    def __call__(self,o):
        a=np.full((o.n,16),np.nan);a[:,[0,1,2,5,11,12,13]]=o.log[:o.n,[0,1,2,5,11,12,13]]
        return self.inner(common.P.Obs(a,o.n))


def episode(job):
    name,scene,wave,scale=job;fam,sid,p,T=scene
    args=dict(single_model=False,**CE.FINAL)
    if name=='baseline':pol=ServoID(**args)
    elif name=='amplitude_only':pol=JointPolicy(wave=None,new_fit=False,amplitude_scale=scale,**args)
    elif name=='long_constant':pol=JointPolicy(wave=None,new_fit=True,constant=True,**args)
    elif name=='varying_no_rate':pol=JointPolicy(wave=None,new_fit=True,rate=False,**args)
    elif name=='matched_gate_short':pol=JointPolicy(wave=None,new_fit=True,rate=False,replace_value='short',**args)
    elif name=='joint_lean':pol=JointPolicy(wave=wave,new_fit=True,rate=False,amplitude_scale=scale,**args)
    else:pol=JointPolicy(wave=wave if name in ('probe_only','joint') else None,
        new_fit=name in ('fit_only','joint'),amplitude_scale=scale if name in ('probe_only','joint') else 1.,**args)
    tic=time.perf_counter();log,_=common.run(MeasuredOnly(pol),p,T,sid)
    metrics=common.task_metrics(log,20.8)
    mask=log[:,0]>=1.5;t=log[mask,0];f=log[mask,3]
    frequencies=[20.8]
    if name in ('probe_only','joint','joint_lean') and wave is not None:
        spectrum=np.fft.rfft(wave)
        frequencies.extend([i/(len(wave)*.008) for i in range(1,len(spectrum)) if abs(spectrum[i])>1e-5*max(abs(spectrum))])
    B=np.column_stack([np.ones(len(t)),t-t.mean()]);tones=np.column_stack([fun(2*np.pi*freq*t) for freq in sorted(set(frequencies)) for fun in (np.sin,np.cos)])
    beta=np.linalg.lstsq(np.column_stack([B,tones]),f,rcond=None)[0]
    metrics['probe_coherent_true_force_rms']=float(np.sqrt(np.mean((tones@beta[2:])**2)))
    fits=getattr(pol,'contact_fits',[]);calls=np.asarray(pol.t_call)[126:]*1000
    # Scoring after the complete episode, unavailable to the controller.
    for f in fits:
        i=int(np.argmin(abs(log[:,0]-f['time'])));f['scoring_reference_k']=float(log[i,8])
        ahead=(log[:,0]>f['time'])&(log[:,0]<=f['time']+.128)
        if np.any(ahead):f['scoring_next_force_rmse']=float(np.sqrt(np.mean((log[ahead,3]-4.)**2)))
    return dict(model=name,family=fam,sid=sid,metrics=metrics,elapsed=time.perf_counter()-tic,
                fits=fits,used=getattr(pol,'used',0),
                call_ms_P50_P90_P99=np.percentile(calls,[50,90,99]).tolist(),
                probe_records=getattr(pol,'probe_records',[]))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--models',nargs='+',default=['baseline','amplitude_only','probe_only','fit_only','joint'])
    ap.add_argument('--families',nargs='+',default=['S2','T3m']);ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--probe-budget',choices=['command','response'],default='command')
    ap.add_argument('--include-original',action='store_true')
    ap.add_argument('--label',default='joint_pilot_v1');args=ap.parse_args()
    physical=args.probe_budget=='response'
    probe=optimize_probe(n=32,rate=True,physical_force_budget=physical);sigma=.04;target=math.log(1.2)
    probe['selected_kind']='periodic_fourier'
    if args.include_original:
        reference=sine_information(n=32,rate=True)[0]
        probe['original_sine_information']=reference;probe['periodic_information']=probe['information']
        if reference>=probe['information']:
            probe.update(wave=None,selected_kind='original_sine',information=reference)
    # Minimal scale for nominal endpoint precision; upper limit 2 is a declared
    # probe budget, NOT a guarantee on actual force excursions.
    scale=1. if physical else min(2.,max(1.,sigma/(target*np.sqrt(probe['information'])*CE.FINAL['a_force'])))
    scenes=[s for s in CE.cohort('diag',1) if s[0] in args.families]
    jobs=[(name,scene,probe['wave'],scale) for scene in scenes for name in args.models]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'local_contact.py',HERE/'joint_policy.py']}
    result=dict(protocol=dict(scientific_questions=['current gain during material transition','movement/material confounding'],
        families=args.families,scenarios_per_family=1,models=args.models,expected_rows=len(jobs),
        information_window=32,source_sha256=hashes,amplitude_scale=scale,probe=probe,
        probe_budget=args.probe_budget,
        original_probe_is_feasible_candidate=args.include_original,
        truth_for_scoring_only=True,tuning=False,plant_changed=False,final_used=False,
        computational_delay_simulated=False),rows=[])
    out=HERE/'results'/(args.label+'.json');out.parent.mkdir(exist_ok=True)
    if out.exists():raise ValueError('use a new label; preserve the original result')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(episode,j) for j in jobs]):
            r=future.result();result['rows'].append(r)
            out.write_text(json.dumps(result,indent=2),encoding='utf8')
            print(r['family'],r['model'],'RMSE',round(r['metrics']['rmse'],4),
                  'peak',round(r['metrics']['peak'],2),'used',r['used'],flush=True)
    print(out)


if __name__=='__main__':main()
