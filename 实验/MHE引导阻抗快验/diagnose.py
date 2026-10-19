"""Does a low one-step wrench error imply a useful decision gradient?"""
import json
import numpy as np
from scipy.optimize import minimize
from core import cost_gradient,rollout,encode
from run import ROOT


def main():
    meta=json.loads((ROOT/'first_pilot.json').read_text());data=np.load(ROOT/'first_pilot_traces.npz');cfg=np.array(meta['cfg']);rows=[]
    for case,p0 in enumerate(meta['parameters']):
        p=np.array(p0);cfg[7]=p[9];tr=data[f'{case}_fixed'];zr=tr[0,2]
        for i in range(40,len(tr)-24,40):
            r=tr[i];t=r[0];y=r[1:11];q=r[17:24];j=i-i%5-2
            anchor=np.array([p[9]*tr[j,0]+tr[j,1],tr[j,0]])
            # Reconstruct exactly the causal smoothed noisy sensing stream.
            rng=np.random.default_rng(meta['config']['seed']*100+case);noise=rng.normal(size=(len(tr),3))*[.03,.03,.0001];sensed=np.array([3.,.8,0.])
            for h in range(i+1):sensed=cfg[6]*(tr[max(0,h-2),11:14]+noise[max(0,h-2)])+(1-cfg[6])*sensed
            qc=np.r_[q,sensed];a=np.tile(encode(np.r_[cfg[:4],0.]),2)
            v,g=cost_gradient(t,y,zr,qc,anchor,p[9],cfg,a,p,False)
            vt,gt=cost_gradient(t,y,zr,qc,anchor,p[9],cfg,a,p,True)
            cs=float(g@gt/max(1e-12,np.linalg.norm(g)*np.linalg.norm(gt)))
            fn=lambda aa:cost_gradient(t,y,zr,qc,anchor,p[9],cfg,aa,p,False)
            opt=minimize(fn,a,jac=True,bounds=[(-1,1)]*10,method='L-BFGS-B',options={'maxiter':15})
            after=rollout(t,y,zr,qc,anchor,p[9],cfg,opt.x,p,True)
            rows.append({'case':case,'t':t,'gradient_cosine':cs,'predicted_before':v,'predicted_after':float(opt.fun),'true_before':vt,'true_after':after})
    out={'rows':rows,'negative_gradient_fraction':float(np.mean([r['gradient_cosine']<0 for r in rows])),'median_gradient_cosine':float(np.median([r['gradient_cosine'] for r in rows])),'true_cost_increase_fraction':float(np.mean([r['true_after']>r['true_before'] for r in rows]))}
    (ROOT/'pilot_diagnosis.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='rows'})


if __name__=='__main__':main()
