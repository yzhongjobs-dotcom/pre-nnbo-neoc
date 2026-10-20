"""Select the Jin-MS task regularizer on eight training scenes only."""
import json
from pathlib import Path
import numpy as np
from run import cases,episode
from baseline import JinMS

if __name__=='__main__':
    target=Path(__file__).with_name('baseline_tuning.json')
    if target.exists():raise FileExistsError(target)
    cs=cases(8,16101);out=[]
    for reg in [1e-6,1e-5,1e-4]:
        rows=[]
        for c in cs:
            metrics,tr,_=episode(c,'jin_ms',1610100+c['id'],model=JinMS(reg))
            score=1.6*np.mean((tr[:,1]/.002)**2+(tr[:,2]/.001)**2+.04*np.sum((tr[:,3:5]/.015)**2,axis=1)+((tr[:,10]-.6)/.1)**2)
            rows.append({'case':c['id'],'physical_running_score':float(score),'force_rmse_N':metrics['force_rmse_N'],'constraints':metrics['constraint_violating_samples']})
        out.append({'regularizer':reg,'cases':rows,'mean_physical_running_score':float(np.mean([r['physical_running_score'] for r in rows]))});print(out[-1],flush=True)
    target.write_text(json.dumps({'selection':'training scenes only; physical running score, shared scales fixed in physics.py',
        'results':out,'selected':min(out,key=lambda z:z['mean_physical_running_score'])['regularizer']},indent=2))
