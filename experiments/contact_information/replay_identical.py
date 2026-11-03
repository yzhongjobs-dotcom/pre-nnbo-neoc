"""Same-input replay of endpoint gain and held-out conditional force prediction.

Reuses baseline measurement files. No new simulations, no tuning. Logged
instantaneous reference K is scoring only and is NOT an optimal control gain.
One-step prediction conditions on the next available encoder; next force is
hidden from prediction. Windows are trajectory diagnostics, not independent n.
"""
import hashlib
import json
import math
import argparse
import sys
from pathlib import Path
import numpy as np
from local_contact import fit,features,sine_information,optimize_probe

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'motion_aware_scan'))
import common,cohort_eval as CE


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--label',default='identical_replay_v1');args=ap.parse_args()
    rows=[];inputs={}
    for fam in ('S2','T3m','T2m-stiff'):
        path=next((HERE.parent/'contact_model'/'results').glob(f'finite_patch_input_{fam}_*.npz'))
        inputs[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        a=np.load(path);m=a['measured'];log=np.full((len(m),16),np.nan)
        log[:,[0,1,2,5,11,12,13]]=m
        for end in range(192,len(m)-1,16):
            tau=float(a['baseline_tau'][end-1]);ref=float(a['scoring_truth_k'][end-1]);p=CE.FINAL
            for name in ('original16','constant32','varying32','varying_rate32'):
                if name=='original16':
                    r=common.window_fit(log,end,16,3,.04,8,math.exp(-.001/tau),1-math.exp(-2*np.pi*30*.001),1,True)
                    if r['info']<=0 or not np.isfinite(r['k']):continue
                    r.update(sigma=1/math.sqrt(r['info']),dof=max(r['valid']-5,1),coefficients=np.array([r['k']]),
                             start_gain=r['k'],tw=.128,ye=m[end-1,2])
                    r['risk']=r['sigma']/max(r['k'],1)*(1+math.sqrt(r['dof']*max(r['chi2']-1,0)))
                else:
                    r=fit(m[:end],nw=32,tau=tau,rate=name=='varying_rate32',constant=name=='constant32',return_aux=True)
                if r is None or r['k']<=0:continue
                accepted=bool(200<=r['k']<=3e5 and r['start_gain']>0 and r['risk']<=math.log(1.2) and r['chi2']<=p['chi2_max'])
                # Construct the held-out interval in exactly the fitted domain.
                B,U=features(np.array([m[end-1,0]]),np.array([m[end-1,2]]),np.array([m[end,2]]),
                             r['te'],r['ye'],r['tw'],tau,30.)
                if name in ('original16','constant32'):B=B[:,:4];U=U[:,:1]
                elif name=='varying32':U=U[:,:2]
                zp=float((B@r['beta']+U@r['coefficients'])[0])
                fp=math.exp(-2*np.pi*30*.008)*m[end-1,3]+zp
                rows.append(dict(family=fam,end=end,time=float(m[end-1,0]),model=name,k=r['k'],
                    sigma=r['sigma'],risk=r['risk'],accepted=accepted,reference_for_scoring=ref,
                    reference_log_error=abs(math.log(r['k']/ref)),conditional_force_error=fp-m[end,3]))
    summaries=[]
    for fam in ('S2','T3m','T2m-stiff'):
        for model in ('original16','constant32','varying32','varying_rate32'):
            rr=[r for r in rows if r['family']==fam and r['model']==model];accepted=[r for r in rr if r['accepted']]
            summary=dict(family=fam,model=model,windows=len(rr),accepted=len(accepted),
                reference_log_P50_P90=np.percentile([r['reference_log_error'] for r in rr],[50,90]).tolist(),
                accepted_reference_log_P50_P90=(np.percentile([r['reference_log_error'] for r in accepted],[50,90]).tolist() if accepted else None),
                conditional_force_RMSE=float(np.sqrt(np.mean([r['conditional_force_error']**2 for r in rr]))))
            summaries.append(summary);print(fam,model,summary,flush=True)
    probe=optimize_probe(n=32,physical_force_budget=True)
    information=dict(original_sine32=sine_information()[0],selected_periodic32=probe['information'])
    out=dict(protocol=dict(input_hashes=inputs,all_models_same_measurements=True,truth_for_scoring_only=True,
                           window_stride=16,windows_independent=False,tuning=False,final_used=False),
             predicted_information=information,summary=summaries,rows=rows)
    target=HERE/'results'/(args.label+'.json')
    if target.exists():raise ValueError('preserve existing replay')
    target.write_text(json.dumps(out,indent=2),encoding='utf8')
    print('information',information)


if __name__=='__main__':main()
