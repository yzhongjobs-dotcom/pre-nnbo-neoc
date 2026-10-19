"""Independent axis regularizers are allowed by Jin's Q/R formulation."""
import json
from pathlib import Path
import numpy as np
from run import cases,episode
from baseline import JinMS

def scores(tr):
    tangent=1.6*np.mean((tr[:,1]/.002)**2+.04*(tr[:,3]/.015)**2)
    normal=1.6*np.mean((tr[:,2]/.001)**2+.04*(tr[:,4]/.015)**2+((tr[:,10]-.6)/.1)**2)
    return [float(tangent),float(normal)]

if __name__=='__main__':
    target=Path(__file__).with_name('baseline_axis_tuning.json')
    if target.exists():raise FileExistsError(target)
    cs=cases(8,16101);out=[]
    for reg in [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4]:
        rows=[]
        for c in cs:
            m,tr,_=episode(c,'jin_ms',1610100+c['id'],model=JinMS(reg));rows.append({'case':c['id'],'scores':scores(tr),'constraints':m['constraint_violating_samples']})
        out.append({'regularizer':reg,'cases':rows,'means':np.mean([r['scores'] for r in rows],axis=0).tolist()});print(out[-1],flush=True)
    chosen=[min(out,key=lambda z:z['means'][i])['regularizer'] for i in range(2)];verified=[]
    for c in cs:
        m,tr,_=episode(c,'jin_ms',1610100+c['id'],model=JinMS(chosen));verified.append({'case':c['id'],'scores':scores(tr),'constraints':m['constraint_violating_samples']})
    target.write_text(json.dumps({'selection':'eight training scenes only; choose each axis weight from the same fixed logarithmic grid; verify combined controller',
         'results':out,'selected':chosen,'combined_training_validation':verified},indent=2))
