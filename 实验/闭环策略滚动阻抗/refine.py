"""Refine the new control trajectories without changing their policies."""
from pathlib import Path
import hashlib,json
import numpy as np
import controller as c
from run import noises
HERE=Path(__file__).resolve().parent
fit=json.loads((HERE/'pilot.json').read_text());protocol=json.loads((HERE/'protocol.json').read_text())
data=np.load(HERE/'paired.npz');ps=data['test_params'];ns=noises(ps,protocol['noise_seed'])
cfg=np.array(fit['cfg']);bank=np.array(fit['bank']);scales=np.array(fit['scales']);rows=[]
for case in [0,1,2]:
    for name,index,ensemble in [('constant_gain',0,False),('learned_joint',2,False),('ensemble_closed',-1,True)]:
        latency=protocol['timing']['ensemble_closed']['p95_ms']/1000 if index<0 else .000003
        value,_=c.episode(ps[case],ns[case],cfg,bank,scales,.16,ensemble,True,latency,
                          policy_index=index,dt=.00025,nq=257)
        old=float(data[name][case,0]);row={'case':case,'method':name,'original_J':old,
                'fine_J':float(value[0]),'relative_difference':float(value[0]/old-1)}
        rows.append(row);print(row,flush=True)
out={'rows':rows,'scope':'Three fixed test cases; physical plant refinement only, no new independent cases or global error bound.',
     'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['controller.py','refine.py','protocol.json','paired.npz']}}
(HERE/'refinement.json').write_text(json.dumps(out,indent=2)+'\n')
