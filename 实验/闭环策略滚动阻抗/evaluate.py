"""Validation then a frozen fresh-case test; do not select from the test."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import controller as c
import model as m
import coupled_contact as cm
from run import noises
from pilot import context,measured

HERE=Path(__file__).resolve().parent


def summary(a,b):
    rng=np.random.default_rng(9190);ids=rng.integers(0,len(a),size=(6000,len(a)))
    benefit=1-a[ids,0].mean(1)/b[ids,0].mean(1)
    return {'mean_J':float(a[:,0].mean()),'force_rms_N':float(a[:,1].mean()),
       'normal_pose_rms_mm':float(a[:,10].mean()*1000),'max_force_N':float(a[:,3].max()),
       'loss_fraction':float(a[:,4].mean()),'saturation_fraction':float(a[:,6].mean()),
       'energy_active_fraction':float(a[:,5].mean()),
       'relative_mean_J_benefit':float(1-a[:,0].mean()/b[:,0].mean()),
       'paired95pct':np.quantile(benefit,[.025,.975]).tolist()}


def batch(ps,ns,cfg,bank,scales,index,horizon,ensemble,consistent,latency):
    values=[];traces={}
    for case,(p,n) in enumerate(zip(ps,ns)):
        val,tr=c.episode(p,n,cfg,bank,scales,horizon,ensemble,consistent,latency,policy_index=index)
        values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
        if case<3:traces[str(case)]=tr
    return np.array(values),traces


def main():
    pilot=json.loads((HERE/'pilot.json').read_text());cfg=np.array(pilot['cfg']);bank=np.array(pilot['bank']);scales=np.array(pilot['scales'])
    vp=m.sample(18,8201);vn=noises(vp,8202)
    variants=[('nominal_closed',.16,False,True),('ensemble_closed',.16,True,True)]
    timings={};all_contexts=[]
    for case,where in [(0,30),(1,150),(2,300)]:all_contexts.append(context(vp[case],vn[case],cfg,bank[0],where))
    for name,h,e,ob in variants:
        samples=[]
        for ctx in all_contexts:samples.extend(measured(ctx,bank,scales,h,e,ob)['samples_ms'])
        timings[name]={'median_ms':float(np.median(samples)),'p95_ms':float(np.quantile(samples,.95)),'n':len(samples)}
        print('timing',name,timings[name],flush=True)
    validation={};va={}
    for index,name in enumerate(['constant_gain','analytic_joint','learned_joint','learned_gain']):
        a,_=batch(vp,vn,cfg,bank,scales,index,.16,False,True,.000003)
        validation[name]=float(a[:,0].mean());va[name]=a
        print('validation',name,validation[name],flush=True)
    for name,h,e,ob in variants:
        a,_=batch(vp,vn,cfg,bank,scales,-1,h,e,ob,timings[name]['p95_ms']/1000)
        validation[name]=float(a[:,0].mean());va[name]=a
        print('validation',name,validation[name],flush=True)
    baseline=min(['constant_gain','analytic_joint','learned_joint','learned_gain'],key=validation.get)
    source={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['controller.py','pilot.py','pilot.json','evaluate.py']}
    protocol={'stage':'Written after validation and before fresh test generation','test_seed':9101,'noise_seed':9102,'n':24,
              'selected_simple_or_learned_baseline':baseline,'validation_J':validation,'timing':timings,
              'variants':variants,'source_sha256':source,
              'scope':'Fixed corrected reduced model and old joint objective. Not a physically calibrated task or complete prior reproduction. All methods have causal deployment information.'}
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');np.savez_compressed(HERE/'validation.npz',**va)
    ps=m.sample(24,9101);ns=noises(ps,9102);arrays={'test_params':ps};traces={}
    for index,name in enumerate(['constant_gain','analytic_joint','learned_joint','learned_gain']):
        arrays[name],tr=batch(ps,ns,cfg,bank,scales,index,.16,False,True,.000003)
        for case,v in tr.items():traces[name+'_'+case]=v
    for name,h,e,ob in variants:
        arrays[name],tr=batch(ps,ns,cfg,bank,scales,-1,h,e,ob,timings[name]['p95_ms']/1000)
        for case,v in tr.items():traces[name+'_'+case]=v
    # Previous target-gain MPC on the same corrected plant and tuned common cfg.
    latency=json.loads((c.BASE/'results.json').read_text())['timing']['spatial_joint']['p95_ms']/1000
    old=[]
    for p,n in zip(ps,ns):
        val,tr=cm.episode(p,n,6,cfg,latency=latency)
        old.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
    arrays['previous_target_MPC']=np.array(old)
    rows={name:summary(a,arrays[baseline]) for name,a in arrays.items() if name!='test_params'}
    out={'summary_vs_validation_selected_baseline':rows,'baseline':baseline,
         'ensemble_vs_nominal':summary(arrays['ensemble_closed'],arrays['nominal_closed']),
         'ensemble_vs_analytic':summary(arrays['ensemble_closed'],arrays['analytic_joint']),
         'ensemble_vs_constant':summary(arrays['ensemble_closed'],arrays['constant_gain']),
         'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
         'source_sha256':source,
         'scope':'One trained policy family and one independent 24-case synthetic batch. Intervals omit training uncertainty; prototype timing excludes complete I/O/OS.'}
    np.savez_compressed(HERE/'paired.npz',**arrays);np.savez_compressed(HERE/'traces.npz',**traces)
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
