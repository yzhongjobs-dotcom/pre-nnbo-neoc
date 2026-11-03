"""Paired diagnosis-only integration of the corrected contact model.

No tuning or plant changes; all policies receive truth-masked observations.
Reports task RMSE, runaway, loss, model usage, timing and per-fit results.
"""
import argparse
import json
import math
import sys
import time
import hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0,str(HERE.parent/'method_redesign'))
import common,cohort_eval as CE
from v1_servo_id import ServoID
from dynamic_contact import Memory
from dynamic_contact_policy import ContactPolicy

CONFIGS={
 'baseline':None,
 'static_linear_gain':dict(spatial=False,memory=Memory(0.,.2),interface='gain'),
 'memory_linear_gain':dict(spatial=False,memory=Memory(.5,.2),interface='gain'),
 'memory_patch_gain':dict(spatial=True,memory=Memory(.5,.2),interface='gain'),
 'memory_patch_target':dict(spatial=True,memory=Memory(.5,.2),interface='target'),
 'memory_line_target':dict(spatial=True,memory=Memory(.5,.2),interface='target',contact='line'),
 'memory_patch_target_low':dict(spatial=True,memory=Memory(.2,.08),interface='target'),
 'memory_patch_target_high':dict(spatial=True,memory=Memory(.8,.6),interface='target'),
 'memory_patch_target_delay':dict(spatial=True,memory=Memory(.5,.2),interface='target',latency=.128)}
CONFIGS.update({
 'static_patch_target':dict(spatial=True,memory=Memory(0.,.2),interface='target'),
 'memory_patch_target_posterior':dict(spatial=True,memory=Memory(.5,.2),interface='target',data_confidence=False),
 'memory_linear_hold_posterior':dict(spatial=False,memory=Memory(.5,.2),interface='gain',fallback='hold',data_confidence=False)})
CONFIGS.update({
 'static_linear_checked':dict(spatial=False,memory=Memory(0.,.2),interface='gain',validate_next=True),
 'memory_linear_checked':dict(spatial=False,memory=Memory(.5,.2),interface='gain',validate_next=True),
 'memory_linear_checked_target':dict(spatial=False,memory=Memory(.5,.2),interface='target',validate_next=True),
 'memory_patch_checked_gain':dict(spatial=True,memory=Memory(.5,.2),interface='gain',validate_next=True),
 'memory_patch_checked_target':dict(spatial=True,memory=Memory(.5,.2),interface='target',validate_next=True),
 'static_patch_checked_target':dict(spatial=True,memory=Memory(0.,.2),interface='target',validate_next=True),
 'memory_line_checked_target':dict(spatial=True,memory=Memory(.5,.2),interface='target',contact='line',validate_next=True)})
for old,new in [('memory_linear_checked','linear_stable_gain'),
                ('memory_patch_checked_gain','patch_stable_gain'),
                ('memory_patch_checked_target','patch_stable_target'),
                ('static_patch_checked_target','static_patch_stable_target'),
                ('memory_line_checked_target','line_stable_target')]:
    CONFIGS[new]=dict(CONFIGS[old],stable_jacobian=True)
    CONFIGS[new+'_uncertainty']=dict(CONFIGS[new],propagate_contact=True)
    CONFIGS[new+'_level']=dict(CONFIGS[new],propagate_contact=True,uncertainty_domain='level')


class Masked:
    direct=True
    def __init__(self,inner):self.inner=inner
    def reset(self,p):
        # Controller receives its nominal calibration, never a changed plant's
        # actual servo/filter parameters. V1 may identify servo lag from data.
        self.inner.reset(dict(common.P.PLATFORMS['PC125']))
    def __call__(self,o):
        x=o.log[:o.n].copy();x[:,[3,4,6,7,8,14,15]]=np.nan
        return self.inner(common.P.Obs(x,o.n))


def episode(job):
    name,(fam,sid,p,T)=job
    pol=ServoID(single_model=False,**CE.FINAL) if name=='baseline' else ContactPolicy(single_model=False,**CE.FINAL,**CONFIGS[name])
    tic=time.perf_counter();log,_=common.run(Masked(pol),p,T,sid)
    metrics=common.task_metrics(log,20.8)
    records=getattr(pol,'fit_records',[])
    # All truth-dependent task/gain attribution occurs AFTER the episode.
    # These scoring fields are never returned through the observation wrapper.
    for record in records:
        use=(log[:,0]>record['time'])&(log[:,0]<=record['time']+.128)
        if np.any(use):
            record['scoring_next_window_force_rmse']=float(np.sqrt(np.mean((log[use,3]-common.FD)**2)))
            record['scoring_next_window_peak']=float(log[use,3].max())
    times=[r['elapsed'] for r in records if 'elapsed'in r]
    variances=np.asarray(getattr(pol,'observation_variances',[]))
    return dict(model=name,family=fam,sid=sid,metrics=metrics,
                elapsed=time.perf_counter()-tic,model_ticks=getattr(pol,'model_ticks',0),
                contact_ticks=getattr(pol,'contact_ticks',0),fits=records,
                contact_rejections=getattr(pol,'contact_rejections',[]),
                unobserved_corrections=getattr(pol,'unobserved_corrections',0),
                parameter_to_force_variance_P50_P90=(np.percentile(variances[:,0]/variances[:,1],[50,90]).tolist()
                                                       if len(variances) else None),
                call_ms_P50_P90_P99=(np.percentile(pol.t_call,[50,90,99])*1000).tolist(),
                fit_ms_P50_P90=(np.percentile(times,[50,90])*1000).tolist() if times else None)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--scenes',type=int,default=2)
    ap.add_argument('--workers',type=int,default=2);ap.add_argument('--families',nargs='+',default=list(CE.FAMS))
    ap.add_argument('--models',nargs='+',default=['baseline','static_linear_gain','memory_linear_gain','memory_patch_gain','memory_patch_target'])
    ap.add_argument('--label',default='dynamic_contact_closed_loop_v1');args=ap.parse_args()
    cohort=[s for s in CE.cohort('diag',args.scenes) if s[0]in args.families]
    jobs=[(name,s) for s in cohort for name in args.models]
    rows=[]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),HERE/'dynamic_contact.py',HERE/'dynamic_contact_policy.py',
             HERE/'finite_patch_model.py',HERE/'line_patch_model.py']}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(episode,jobs,chunksize=1):
            rows.append(row);m=row['metrics']
            print(row['family'],row['sid'],row['model'],'RMSE',round(m['rmse'],4),'peak',round(m['peak'],2),
                  'model ticks',row['model_ticks'],flush=True)
            out=dict(protocol=dict(cohort='diag',scenes_per_family=args.scenes,models=args.models,
                                   truth_masked=True,tuning=False,plant_changed=False,final_used=False,
                                   wall_time_included_as_plant_delay=False,
                                   expected_rows=len(jobs),source_sha256=hashes,
                                   configs={n:({k:(vars(v) if isinstance(v,Memory) else v) for k,v in CONFIGS[n].items()}
                                               if CONFIGS[n] else None) for n in args.models},
                                   controller_calibration='nominal PC125'),rows=rows)
            (HERE/'results'/(args.label+'.json')).write_text(json.dumps(out,indent=2,default=float),encoding='utf-8')
    for fam in args.families:
        for name in args.models:
            r=[v for v in rows if v['family']==fam and v['model']==name]
            vals=[v['metrics']['rmse'] for v in r]
            print(fam,name,'RMSE median/mean',np.round([np.median(vals),np.mean(vals)],4),
                  'runaways',sum(v['metrics']['peak']>50 or v['metrics']['diverged'] for v in r),
                  'long loss',sum(v['metrics']['loss_max']>.5 for v in r))


if __name__=='__main__':main()
