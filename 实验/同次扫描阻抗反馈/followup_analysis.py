"""Report the stronger simple comparison and the follow-up reversal explicitly."""
from pathlib import Path
import hashlib,json
import numpy as np
import model as m
from run import noises
HERE=Path(__file__).resolve().parent
puregain=json.loads((HERE/'pure_gain_fit.json').read_text());cfg=np.array(puregain['cfg'])
new=np.load(HERE/'pure_gain_paired.npz')

def comparison(a,b):
 rng=np.random.default_rng(8591);ids=rng.integers(0,len(a),size=(5000,len(a)))
 r=1-a[ids,0].mean(1)/b[ids,0].mean(1)
 return {'relative_mean_J_benefit':float(1-a[:,0].mean()/b[:,0].mean()),
   'paired_95pct':np.quantile(r,[.025,.975]).tolist(),
   'worst_instance_absolute_regression':float(np.max(a[:,0]-b[:,0])),
   'better_instances':int(np.sum(a[:,0]<b[:,0])),'n':len(a)}

out={'fresh24_gain_learning_vs_continuous_constant':comparison(new['learned_gain_only'],new['continuous_gain_only']),
 'fresh24_joint_learning_vs_continuous_constant':comparison(new['joint_learned_unchanged'],new['continuous_gain_only']),
 'scope':'Fresh24 is independent of pure-gain fitting; old48 below is a transparent retrospective stronger-baseline check, not another independent experiment.'}
base=np.load(HERE/'paired.npz');ps=base['test_params'];ns=noises(ps,8402)
latency=json.loads((HERE/'results.json').read_text())['timing']['learned_feedback']['p95_ms']/1000
values=[]
for p,n in zip(ps,ns):
 val,tr,_,_=m.episode(p,n,12,cfg,policy=np.array(puregain['policies']['continuous_gain_only']),geometry_prior=True,latency=latency)
 values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
values=np.array(values)
out['retrospective48_joint_learning_vs_new_constant']=comparison(base['informed_prior_learned_feedback'],values)
out['retrospective48_new_constant_mean_J']=float(values[:,0].mean())
# Diagnostic cost components use force RMS and sampled pose RMS; not a new utility.
duration=np.round(.064/new['test_params'][:,9]/.01)*.01
for key in ['continuous_gain_only','learned_gain_only','joint_learned_unchanged']:
 a=new[key]
 out[key+'_component_diagnostic']={'force_integral_component_mean':float(np.mean(duration*(a[:,1]/.3)**2)),
   'sampled_normal_pose_component_mean':float(np.mean(duration*(a[:,10]/.0015)**2))}
out['source_sha256']={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','followup_analysis.py','pure_gain_fit.json']}
(HERE/'followup_analysis.json').write_text(json.dumps(out,indent=2)+'\n')
np.savez_compressed(HERE/'retrospective_simple.npz',constant_gain_only=values)
print(json.dumps(out,indent=2))
