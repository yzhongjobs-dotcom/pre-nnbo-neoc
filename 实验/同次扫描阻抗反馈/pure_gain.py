"""Follow-up after the joint-action test: isolate gain learning, no normal FF action.

Uses a NEW test seed after fits/protocol are frozen. Previous 48 cases are not
relabelled as independent verification of this post-test method refinement.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize
import model as m
from run import noises
from learn_policy import pad
HERE=Path(__file__).resolve().parent


@njit(cache=True)
def task_loss(theta,ps,ns,counts,cfg):
    value=0.
    for i in range(len(ps)):
        value+=m.episode(ps[i],ns[i,:counts[i]],12,cfg,policy=theta,geometry_prior=True)[0][0]/len(ps)
    return value


def summarize(a,b):
    rng=np.random.default_rng(8590);ids=rng.integers(0,len(a),size=(4000,len(a)))
    draws=1-a[ids,0].mean(1)/b[ids,0].mean(1)
    return {'mean_J':float(a[:,0].mean()),'force_rms_N':float(a[:,1].mean()),
      'normal_pose_rms_mm':float(a[:,10].mean()*1000),'worst_peak_force_N':float(a[:,3].max()),
      'relative_mean_J_benefit':float(1-a[:,0].mean()/b[:,0].mean()),
      'paired_95pct':np.quantile(draws,[.025,.975]).tolist()}


def main():
    start=time.perf_counter();final=json.loads((HERE/'final_policies.json').read_text());cfg=np.array(final['cfg'])
    ps=m.sample(18,8101);ns,counts=pad(noises(ps,8102));vp=m.sample(18,8201);vn,vc=pad(noises(vp,8202))
    templates={};optim={}
    for name,theta,ids in [
       ('continuous_gain_only',np.array(final['policies']['continuous_constant']),[0,3]),
       ('analytic_gain_only',np.array(final['policies']['analytic_inverse']),[3])]:
        theta=theta.copy();theta[5]=0.;theta[6]=1.
        calls=[0]
        def fun(x):
            tt=theta.copy();tt[ids]=x;calls[0]+=1
            return task_loss(tt,ps,ns,counts,cfg)
        bounds=[(np.log(1),np.log(800)) if j==0 else (np.log(2),np.log(150)) for j in ids]
        fit=minimize(fun,theta[ids],method='Powell',bounds=bounds,options={'maxiter':4,'maxfev':100,'xtol':.02,'ftol':1e-5})
        theta[ids]=fit.x;templates[name]=theta;optim[name]={'evaluations':calls[0],'train_J':float(fit.fun)}
        print('pure baseline',name,fit.fun,task_loss(theta,vp,vn,vc,cfg),flush=True)
    mean=np.array(final['policies']['learned_feedback']);mean[5]=0.;mean[6]=1.
    rng=np.random.default_rng(8307);sd=np.array([.2,.25,.2,.2,.2,0.,0.,.2])
    lo=np.array([np.log(1),-2,-1,np.log(2),-.75,0.,1.,-1.]);hi=np.array([np.log(800),1,1,np.log(150),.75,0.,1.,1.])
    best=mean.copy();bestval=task_loss(best,vp,vn,vc,cfg);history=[]
    for iteration in range(6):
        pop=np.clip(mean+rng.normal(size=(24,8))*sd,lo,hi)
        pop[0]=mean;pop[1]=templates['continuous_gain_only'];pop[2]=templates['analytic_gain_only']
        losses=np.array([task_loss(tt,ps,ns,counts,cfg) for tt in pop]);order=np.argsort(losses)
        elite=pop[order[:5]];mean=.3*mean+.7*elite.mean(0);sd=np.maximum(.025,.3*sd+.7*elite.std(0));sd[[5,6]]=0
        val=task_loss(pop[order[0]],vp,vn,vc,cfg)
        if val<bestval:bestval=val;best=pop[order[0]].copy()
        history.append({'iteration':iteration,'train_J':float(losses.min()),'validation_J':float(val)})
        print('pure CEM',iteration,losses.min(),val,flush=True)
    templates['learned_gain_only']=best
    output={'stage':'Post-primary-test method refinement with fresh follow-up cases',
      'cfg':cfg.tolist(),'policies':{k:v.tolist() for k,v in templates.items()},
      'train_seed':8101,'validation_seed':8201,'actor_evaluations':144,'simple_optimizer':optim,
      'history':history,'validation_J':{k:float(task_loss(v,vp,vn,vc,cfg)) for k,v in templates.items()},
      'training_seconds':time.perf_counter()-start,'normal_feedforward_action':0.,
      'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','pure_gain.py']}}
    (HERE/'pure_gain_fit.json').write_text(json.dumps(output,indent=2)+'\n')
    protocol={'fit_sha256':hashlib.sha256((HERE/'pure_gain_fit.json').read_bytes()).hexdigest(),
      'test_seed':8501,'noise_seed':8502,'n':24,'task':'Same joint force/pose loss; pure-gain methods have no normal feedforward',
      'comparison':'Fresh cases after all pure-gain policies are fitted; original joint policy also replayed unchanged',
      'scope':'Exploratory follow-up, not preregistration before primary findings or physical validation'}
    (HERE/'pure_gain_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    # First use of the new seed occurs only after the above freeze.
    tp=m.sample(24,8501);tn=noises(tp,8502)
    timings=json.loads((HERE/'results.json').read_text())['timing']
    latency=timings['learned_feedback']['p95_ms']/1000
    policies=dict(templates);policies['joint_learned_unchanged']=np.array(final['policies']['learned_feedback'])
    policies['joint_analytic_unchanged']=np.array(final['policies']['analytic_inverse'])
    arrays={'test_params':tp};table={}
    for name,theta in policies.items():
        values=[]
        for p,n in zip(tp,tn):
            val,tr,_,_=m.episode(p,n,12,cfg,policy=theta,geometry_prior=True,latency=latency)
            values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
        arrays[name]=np.array(values)
    for name,a in arrays.items():
        if name=='test_params':continue
        table[name]=summarize(a,arrays['analytic_gain_only']);print('new pure test',name,table[name],flush=True)
    result={'summary_vs_analytic_gain_only':table,
      'learned_gain_vs_joint_analytic':summarize(arrays['learned_gain_only'],arrays['joint_analytic_unchanged']),
      'scope':'Fresh 24 synthetic scans, same physical model; single fitted policy; no added sensing or hidden future',
      'protocol_sha256':hashlib.sha256((HERE/'pure_gain_protocol.json').read_bytes()).hexdigest(),
      'model_sha256':output['source_sha256']['model.py']}
    (HERE/'pure_gain_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(HERE/'pure_gain_paired.npz',**arrays)
if __name__=='__main__':main()
