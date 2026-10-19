"""Train-only tuning, retaining native cheap-feedback update frequency."""
import json
import numpy as np
from experiment import ROOT,BASE,episode
from core import m


def main():
    cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);params=m.sample(6,9101)
    records=[]
    for method,values in [('nnbo_adapted',[.5,1.,2.,4.,8.]),('mpvic_os_adapted',[1e-7,1e-6,1e-5,1e-4,1e-3])]:
        for v in values:
            costs=[];fr=[]
            for i,p in enumerate(params):
                r,_=episode(p,910100+i,method,cfg,omega=v if method=='nnbo_adapted' else 10.,qp_reg=v if method=='mpvic_os_adapted' else 1e-7)
                costs.append(r['metrics'][0]);fr.append(r['metrics'][1])
            row={'method':method,'parameter':v,'costs':costs,'force_rms':fr,'mean':float(np.mean(costs))};records.append(row);print(row,flush=True)
    best={m:min([r for r in records if r['method']==m],key=lambda r:r['mean'])['parameter'] for m in ('nnbo_adapted','mpvic_os_adapted')}
    (ROOT/'domain_baseline_tuning_v3.json').write_text(json.dumps({'seed':9101,'n':6,'native_stride':1,'normal_axis_vic':True,'records':records,'best':best},indent=2))


if __name__=='__main__':main()
