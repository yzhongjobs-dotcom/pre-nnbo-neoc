"""Reproducible paired closed-loop validation of the event-aware estimator.

No changes to the simulator or published controller. Only common.window_fit
is replaced temporarily during each policy's stiffness-estimation call.
The scored truth is never visible inside the controller (measured-only wrapper).
"""
import argparse, collections, json, math, os, sys, time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
for p in (HERE,ROOT/'motion_aware_scan',ROOT/'method_redesign',ROOT/'contact_information'):
    sys.path.insert(0,str(p))
import common
from cohort_eval import FINAL, FAMS, cohort
from v1_servo_id import ServoID
from joint_policy import JointPolicy
from event_projection import event_fit

class MeasuredOnly:
    direct=True
    def __init__(self,inner):self.inner=inner
    def reset(self,platform):self.inner.reset(platform)
    def __call__(self,obs):
        a=np.full((obs.n,obs.log.shape[1]),np.nan)
        cols=[0,1,2,5,11,12,13]
        a[:,cols]=obs.log[:obs.n,cols]
        return self.inner(common.P.Obs(a,obs.n))

class EventMixin:
    def reset(self,platform):
        super().reset(platform)
        self.events=collections.Counter()
        self.event_samples=[]
    def _stiffness_est(self,o,n):
        original=common.window_fit
        p=self.params
        base=original(o.log,n,p['nw'],p['degree'],p['noise'],self.ns,
                      self.eb,self.afc,self.ndel,p['enc_y'])
        result=event_fit(o.log,n,p['nw'],p['degree'],p['noise'],self.ns,
                         self.eb,self.afc,self.ndel,p['enc_y'],_base=base)
        et=result['event_type']
        self.events[et]+=1
        if et!='base' and len(self.event_samples)<50:
            self.event_samples.append((float(o.log[n-1,0]),et,float(result['k']),float(base['k'])))
        def substituted(*args,**kwargs):
            return result
        common.window_fit=substituted
        try:
            return super()._stiffness_est(o,n)
        finally:
            common.window_fit=original
class EventV1(EventMixin,ServoID):pass
class EventLean(EventMixin,JointPolicy):pass

def build(config):
    if config=="v1":return ServoID(single_model=False,**FINAL)
    if config=="event_v1":return EventV1(single_model=False,**FINAL)
    if config=="lean":return JointPolicy(wave=None,new_fit=True,rate=False,amplitude_scale=1.,**FINAL)
    if config=="event_lean":return EventLean(wave=None,new_fit=True,rate=False,amplitude_scale=1.,**FINAL)
    raise ValueError(config)

def run_one(config,scene,lag=1.):
    fam,sid,p,T=scene
    pol=build(config)
    plat=dict(common.P.PLATFORMS['PC125'])
    plat['tau_a']*=lag
    start=time.perf_counter()
    log,_=common.run(MeasuredOnly(pol),p,T,sid,platform=plat)
    if log is None:
        return dict(config=config,family=fam,sid=sid,rmse=float('nan'),peak=float('nan'),
                    loss=float('nan'),error='plant returned None',elapsed=time.perf_counter()-start)
    m=common.task_metrics(log,20.8)
    ret=dict(config=config,family=fam,sid=sid,rmse=float(m['rmse']),peak=float(m['peak']),
             loss=float(m['loss_max']),diverged=bool(m['diverged']),
             elapsed=time.perf_counter()-start)
    if hasattr(pol,'events'):
        ret['events']=dict(pol.events)
        ret['event_samples']=pol.event_samples
    # truth only used post-hoc; no access by policy.
    return ret

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--scenes',type=int,default=3)
    pa.add_argument('--lag',type=float,default=1.)
    pa.add_argument('--families',nargs='+',default=list(FAMS))
    pa.add_argument('--configs',nargs='+',default=['v1','event_v1','lean','event_lean'])
    pa.add_argument('--label',default='diagnostic')
    args=pa.parse_args()
    if args.scenes<1 or args.scenes>50:raise ValueError("scenes out of range")
    scenes=[s for s in cohort('diag',args.scenes) if s[0] in args.families]
    rows=[]
    for scene in scenes:
        for config in args.configs:
            start=time.perf_counter()
            try:
                row=run_one(config,scene,args.lag)
            except Exception as exc:
                row=dict(config=config,family=scene[0],sid=scene[1],error=repr(exc))
            rows.append(row)
            print(f'{scene[0]} {scene[1]} {config}: '
                  f'RMSE={row.get("rmse")} peak={row.get("peak")} '
                  f'elapsed={time.perf_counter()-start:.2f}s '
                  f'error={row.get("error","")}',flush=True)
    report={'protocol':vars(args),'summary':{},'rows':rows}
    for fam in args.families:
        for config in args.configs:
            r=[x for x in rows if x.get('family')==fam and x.get('config')==config]
            good=[x for x in r if np.isfinite(x.get('rmse',float('nan')))]
            if not good:
                report['summary'][fam+'|'+config]={'error':'no successful episodes'};continue
            vals=np.array([x['rmse'] for x in good])
            report['summary'][fam+'|'+config]={
                'median':float(np.median(vals)), 'mean':float(np.mean(vals)),
                'failed':len(r)-len(good)+sum(x['peak']>50 for x in good),
                'loss_gt_half':sum(x['loss']>.5 for x in good),
                'events':dict(sum((collections.Counter(x.get('events',{})) for x in r),collections.Counter()))
            }
    path=HERE/f'{args.label}.json'
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf8')
    print('SUMMARY_JSON_START',flush=True)
    print(json.dumps(report['summary'],ensure_ascii=False,allow_nan=False),flush=True)
    print('SUMMARY_JSON_END',flush=True)

if __name__=='__main__':main()
