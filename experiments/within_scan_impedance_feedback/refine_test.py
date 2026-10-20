"""Same independent test cases at finer numerical resolution, not new test data."""
from pathlib import Path
import json,hashlib
import numpy as np
import model as m
from run import noises
HERE=Path(__file__).resolve().parent
protocol=json.loads((HERE/'test_protocol.json').read_text())
final=json.loads((HERE/'final_policies.json').read_text())
data=np.load(HERE/'paired.npz');ps=data['test_params'];ns=noises(ps,protocol['noise_seed'])
base=json.loads((HERE/'joint_validation.json').read_text())
specs=[('analytic_inverse',12,np.array(final['cfg']),np.array(final['policies']['analytic_inverse'])),
       ('learned_feedback',12,np.array(final['cfg']),np.array(final['policies']['learned_feedback'])),
       ('spatial_joint',6,np.array(base['cfgs']['spatial_joint']),np.zeros(8))]
arrays={};rows={}
for name,method,cfg,policy in specs:
 values=[]
 for p,n in zip(ps,ns):
  metric,tr,_,_=m.episode(p,n,method,cfg,policy=policy,geometry_prior=True,latency=protocol['delay_replay_s'][name],dt=.00025,nq=129)
  values.append(np.r_[metric,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
 values=np.array(values);arrays[name]=values
 coarse=data['informed_prior_'+name]
 rows[name]={'mean_J':float(values[:,0].mean()),'force_rms_N':float(values[:,1].mean()),
    'normal_pose_rms_mm':float(values[:,10].mean()*1000),
    'relative_mean_J_change':float(values[:,0].mean()/coarse[:,0].mean()-1),
    'max_instance_relative_J_change':float(np.max(np.abs(values[:,0]/coarse[:,0]-1)))}
 print('fine',name,rows[name],flush=True)
rng=np.random.default_rng(8492);ids=rng.integers(0,48,size=(4000,48))
a=arrays['analytic_inverse'];b=arrays['learned_feedback']
draw=1-b[ids,0].mean(1)/a[ids,0].mean(1)
summary={'rows':rows,'learned_vs_analytic':{'relative_mean_J_benefit':float(1-b[:,0].mean()/a[:,0].mean()),
 'paired_95pct':np.quantile(draw,[.025,.975]).tolist()},
 'scope':'All same 48 cases, 129 vs 33 spatial nodes and 0.25ms vs 1ms integration; same causal observations and delay replay, not new independent data',
 'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','refine_test.py','test_protocol.json']}}
(HERE/'refinement.json').write_text(json.dumps(summary,indent=2)+'\n')
np.savez_compressed(HERE/'refinement.npz',**arrays)
print('refined comparison',summary['learned_vs_analytic'],flush=True)
