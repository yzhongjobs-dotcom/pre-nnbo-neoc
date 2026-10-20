from pathlib import Path
import hashlib,json
import numpy as np
import model as m
import reference_controller as rc
from run import noises
HERE=Path(__file__).resolve().parent
p=m.sample(18,8101);ns=noises(p,8102)
base=np.array(json.loads((HERE/'final_policies.json').read_text())['cfg'])
choices=[]
settings=[(omega,zeta,1.) for omega in [4.,8.,12.,20.,32.,50.] for zeta in [.7,1.,1.4]]
settings+= [(omega,zeta,blend) for omega in [20.,32.,50.,80.] for zeta in [.4,.7,1.] for blend in [0.,.5,1.]]
for omega,zeta,blend in settings:
   cfg=np.r_[base,omega,zeta];cfg[9]=blend
   value=np.mean([rc.episode(pp,nn,cfg)[0][0] for pp,nn in zip(p,ns)])
   choices.append((float(value),omega,zeta,blend))
   print('reference',omega,zeta,blend,value,flush=True)
best=min(choices);cfg=np.r_[base,best[1:3]];cfg[9]=best[3]
v=m.sample(18,8201);vn=noises(v,8202)
values=np.array([rc.episode(pp,nn,cfg)[0] for pp,nn in zip(v,vn)])
out={'cfg':cfg.tolist(),'train_J':best[0],'validation_J':float(values[:,0].mean()),
     'train_seed':8101,'validation_seed':8201,'choices':choices,
     'scope':'Static-equilibrium quadratic reference feedback; full published controller not claimed; cfg[9] is spatial/temporal blend for this baseline',
     'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','reference_controller.py','tune_reference.py']}}
(HERE/'reference_tuning.json').write_text(json.dumps(out,indent=2)+'\n')
print('reference selected',out,flush=True)
