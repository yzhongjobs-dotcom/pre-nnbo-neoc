"""Power-gradient verification and held-policy transfer to corrected mechanics."""
from pathlib import Path
import hashlib,json
import numpy as np
import model as m
import coupled_contact as cm
from run import noises
HERE=Path(__file__).resolve().parent

def main():
 rng=np.random.default_rng(8610);ps=m.sample(48,8611);checks=[]
 for i,p in enumerate(ps):
  y=cm.initial(p,np.array([1000.,80.,45.25,18.]));t=rng.uniform(0,.064/p[9])
  y[0]=rng.uniform(-.001,.001);y[1]+=rng.uniform(-.0005,.0005);y[2:4]=rng.uniform(-.003,.003,2)
  force,U,diss,pbody,_=cm.contact(t,y,p)
  eps=2e-7;gradient=[]
  for j in [0,1]:
   a=y.copy();b=y.copy();a[j]+=eps;b[j]-=eps
   gradient.append((cm.contact(t,a,p)[1]-cm.contact(t,b,p)[1])/(2*eps))
  a=y.copy();b=y.copy();a[0]-=p[9]*eps;b[0]+=p[9]*eps
  ut=(cm.contact(t+eps,a,p)[1]-cm.contact(t-eps,b,p)[1])/(2*eps)
  Udot=gradient[0]*(p[9]+y[2])+gradient[1]*y[3]+ut
  residual=force[1]*(p[9]+y[2])+force[0]*y[3]-Udot-diss-pbody
  checks.append([residual,diss])
 check={'cases':48,'seed':8611,'max_abs_power_residual_W':float(np.max(np.abs(np.array(checks)[:,0]))),
   'min_dissipation_W':float(np.min(np.array(checks)[:,1])),
   'scope':'Discrete reduced environment potential/velocity-power identity, not real tissue, 6D wrench or closed-loop safety proof'}
 assert check['max_abs_power_residual_W']<1e-5 and check['min_dissipation_W']>=-1e-10
 print('power check',check,flush=True)
 final=json.loads((HERE/'final_policies.json').read_text());pure=json.loads((HERE/'pure_gain_fit.json').read_text())
 cfg=np.array(final['cfg']);timing=json.loads((HERE/'results.json').read_text())['timing']
 specs=[('fixed_PI',0,np.zeros(8)),('analytic_inverse',12,np.array(final['policies']['analytic_inverse'])),
   ('continuous_gain_only',12,np.array(pure['policies']['continuous_gain_only'])),
   ('learned_feedback',12,np.array(final['policies']['learned_feedback'])),
   ('learned_gain_only',12,np.array(pure['policies']['learned_gain_only'])),
   ('spatial_gain',5,np.zeros(8)),('spatial_joint',6,np.zeros(8))]
 source={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','coupled_contact.py','coupled_check.py','final_policies.json','pure_gain_fit.json']}
 protocol={'stage':'Physical-model sensitivity after primary/follow-up results; policies remain frozen',
  'test_seed':8601,'noise_seed':8602,'n':24,'trained_on':'Original local-damping model only; no retraining for this plant',
  'changes':'Material-point convective rate plus gradient of discrete elastic potential and dissipative tangential force',
  'not_claimed':'Full 6D probe mechanics or real-material validation','source_sha256':source}
 (HERE/'coupled_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 ps=m.sample(24,8601);ns=noises(ps,8602);arrays={'test_params':ps};summary={}
 for name,method,theta in specs:
  latency=timing[name if name in timing else 'learned_feedback']['p95_ms']/1000
  values=[]
  for case,(p,n) in enumerate(zip(ps,ns)):
   metric,trace=cm.episode(p,n,method,cfg,theta,latency)
   values.append(np.r_[metric,np.sqrt(np.mean((trace[:,2]-trace[0,2])**2))])
   if name in ['continuous_gain_only','learned_gain_only','learned_feedback'] and case==2:arrays[name+'_trace2']=trace
  values=np.array(values);arrays[name]=values
  summary[name]={'mean_J':float(values[:,0].mean()),'force_rms_N':float(values[:,1].mean()),
    'normal_pose_rms_mm':float(values[:,10].mean()*1000),'tangent_rms_mm':float(values[:,2].mean()*1000),
    'max_force_N':float(values[:,3].max()),'loss_fraction':float(values[:,4].mean()),'saturation_fraction':float(values[:,6].mean())}
  print('corrected plant',name,summary[name],flush=True)
 for name in summary:summary[name]['relative_benefit_vs_continuous_fixed']=float(1-summary[name]['mean_J']/summary['continuous_gain_only']['mean_J'])
 out={'summary':summary,'power_checks':check,'source_sha256':source,
   'scope':'24 fresh cases of a physically coupled reduced plant with frozen policies; no new method superiority or full-model validation claim'}
 (HERE/'coupled_results.json').write_text(json.dumps(out,indent=2)+'\n')
 np.savez_compressed(HERE/'coupled_paired.npz',**arrays)
 print('coupled transfer complete',flush=True)
if __name__=='__main__':main()
