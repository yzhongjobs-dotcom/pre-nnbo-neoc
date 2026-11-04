"""Post-hoc, same-input estimator verification on original simulated trajectories.
Truth is removed from estimator data and used only to score estimates.
"""
import argparse, json, math, os, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
for p in (HERE,ROOT/'motion_aware_scan',ROOT/'method_redesign'):
    sys.path.insert(0,str(p))
import common
from cohort_eval import FINAL,FAMS,cohort
from v1_servo_id import ServoID
from event_projection import event_fit
from validate_closed_loop import MeasuredOnly

def score(a):
    if len(a)==0:return dict(n=0)
    a=np.asarray(a,float)
    return dict(n=len(a),median=float(np.median(a)),p90=float(np.percentile(a,90)),
                low20=float(np.mean(a>math.log(1.2)))) # a signed log(k_true/k_hat)

def compute(scene,lag,stride):
    fam,sid,p,T=scene
    pol=MeasuredOnly(ServoID(single_model=False,**FINAL))
    plat=dict(common.P.PLATFORMS['PC125']);plat['tau_a']*=lag
    log,_=common.run(pol,p,T,sid,platform=plat)
    measured=log.copy()
    # Cut every privileged field; observable columns are the same as closed-loop policy.
    cols=(0,1,2,5,11,12,13)
    measured[:]=np.nan;measured[:,cols]=log[:,cols]
    ns=8;afc=1-math.exp(-2*math.pi*plat['fc']*.001)
    eb=math.exp(-.001/plat['tau_a'])
    out=[]
    for n in range(max(32,int(1./plat['Tc'])),len(log)-1,stride):
        # Encode the same causal history as would be available at time n.
        l=measured[:n]
        b=common.window_fit(l,n,16,3,.04,ns,eb,afc,plat['ndel'],True)
        e=event_fit(l,n,16,3,.04,ns,eb,afc,plat['ndel'],True,_base=b)
        truth=float(log[n-1,8])
        if truth<=0 or not np.isfinite(b['k']) or b['k']<=0 or not np.isfinite(e['k']) or e['k']<=0:
            continue
        j=max(0,n-9)
        previous=max(float(log[j,8]),1.)
        changing=abs(math.log(truth/previous))>.15
        # Actual true tissue surface velocity (scoring label, not a method input).
        v=float(common.P.surface(float(log[n-1,1]),float(log[n-1,0]),p)[1])
        moving=abs(v)>.03
        label=("ridge" if changing else "motion" if moving else "stable")
        out.append(dict(time=float(log[n-1,0]), label=label,
                        old=float(math.log(truth/b['k'])),new=float(math.log(truth/e['k'])),
                        event=e['event_type'],base_chi2=float(b['chi2'])))
    return dict(family=fam,sid=sid,rows=out)

def agg(data):
    ret={}
    for fam in FAMS:
        for label in ("stable","motion","ridge"):
            samples=[x for ep in data if ep['family']==fam for x in ep['rows'] if x['label']==label]
            if not samples:continue
            k=f'{fam}|{label}'
            ret[k]=dict(
                n=len(samples),
                old=dict(median=float(np.median([abs(x['old']) for x in samples])),
                         p90=float(np.percentile([abs(x['old']) for x in samples],90)),
                         under20=float(np.mean([x['old']>math.log(1.2) for x in samples]))),
                new=dict(median=float(np.median([abs(x['new']) for x in samples])),
                         p90=float(np.percentile([abs(x['new']) for x in samples],90)),
                         under20=float(np.mean([x['new']>math.log(1.2) for x in samples]))),
                event_fraction=float(np.mean([x['event']!='base' for x in samples])))
    return ret

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--scenes',type=int,default=2)
    pa.add_argument('--lag',type=float,default=1.)
    pa.add_argument('--stride',type=int,default=4)
    pa.add_argument('--label',default='same_input')
    args=pa.parse_args()
    rows=[]
    for scene in cohort('diag',args.scenes):
        ep=compute(scene,args.lag,args.stride)
        rows.append(ep)
        print(scene[0],scene[1],len(ep['rows']),flush=True)
    summary=agg(rows)
    (HERE/f'{args.label}.json').write_text(json.dumps(dict(summary=summary,episodes=rows),indent=1),encoding='utf8')
    print('ESTIMATE_SUMMARY',json.dumps(summary),flush=True)

if __name__=='__main__':main()
