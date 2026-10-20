"""Task preference study uses NNBO Q/R=5e5,1e6 as declared numerical references."""
import json
import numpy as np
from experiment import ROOT,BASE,episode
from core import m


def main():
    cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);params=m.sample(6,9101);records=[]
    for weight in (12.5,25.):
        cfg[8]=weight
        for method,values in [('nnbo_adapted',[1.,2.,4.,8.,16.]),('mpvic_os_adapted',[1e-8,1e-7,1e-6,1e-5,1e-4])]:
            for v in values:
                costs=[]
                for i,p in enumerate(params):
                    r,_=episode(p,910100+i,method,cfg,omega=v if method=='nnbo_adapted' else 10.,qp_reg=v if method=='mpvic_os_adapted' else 1e-7)
                    costs.append(r['metrics'][0])
                row={'weight':weight,'method':method,'parameter':v,'costs':costs,'mean':float(np.mean(costs))};records.append(row);print(row,flush=True)
    best={str(w):{m:min([r for r in records if r['method']==m and r['weight']==w],key=lambda r:r['mean'])['parameter'] for m in ('nnbo_adapted','mpvic_os_adapted')} for w in (12.5,25.)}
    (ROOT/'taskweight_baseline_tuning.json').write_text(json.dumps({'seed':9101,'n':6,'weights':[12.5,25.],'records':records,'best':best},indent=2))


if __name__=='__main__':main()
