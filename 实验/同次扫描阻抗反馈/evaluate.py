"""Freeze protocol before generating independent cases; no test-based retuning."""
from pathlib import Path
import hashlib,json,time
import numpy as np
import model as m
import reference_controller as rc
from run import noises
HERE=Path(__file__).resolve().parent


def comparisons(final):
    config=json.loads((HERE/'joint_validation.json').read_text())
    specs=[]
    for name,method in [('fixed_PI',0),('observer_fixed',1),('critical',2),('temporal_MPC',4),
        ('spatial_gain',5),('spatial_joint',6),('spatial_ff_only',7),('transport_fixed',8),
        ('task_matched',9),('task_matched_transport',10)]:
        oldname='temporal_joint' if name=='temporal_MPC' else name
        cfg=np.array(config['cfgs'][oldname])
        specs.append((name,method,cfg,np.zeros(8)))
    for name,theta in final['policies'].items():specs.append((name,12,np.array(final['cfg']),np.array(theta)))
    reference=json.loads((HERE/'reference_tuning.json').read_text())
    specs.append(('equilibrium_reference',13,np.array(reference['cfg']),np.zeros(8)))
    actor=np.array(final['policies']['learned_feedback'])
    no_state=actor.copy();no_state[[2,7]]=0.
    fixed_gain=actor.copy();fixed_gain[[1,2,4,7]]=0.
    specs.append(('learned_no_state',12,np.array(final['cfg']),no_state))
    specs.append(('learned_fixed_gains',12,np.array(final['cfg']),fixed_gain))
    return specs


def timing(specs):
    # Recreate representative validation states using the final common information.
    ps=m.sample(18,8201);ns=noises(ps,8202);contexts=[]
    for case in [0,1,2]:
        cfg=specs[0][2].copy();cfg[7]=ps[case,9]
        _,tr,_,_=m.episode(ps[case],ns[case],6,cfg,geometry_prior=True)
        hist=np.zeros((len(tr)+1,8));hist[:-1,:4]=tr[:,1:5]
        hist[:-1,4:7]=tr[:,11:14]+ns[case][:-1]
        hist[:-1,7]=tr[:,0]*ps[case,9]+tr[:,1]
        for i in np.linspace(20,len(tr)-1,24).astype(int):
            y=tr[i,1:11].copy();old=np.r_[tr[i-1,22:25],np.zeros(5)]
            contexts.append((y,tr[0,2],cfg.copy(),hist,i,old,tr[i,14:17].copy()))
    timings={}
    for name,method,cfg0,policy in specs:
        records=[]
        for y,zr,cc,hist,i,old,sensed in contexts[:3]:
            cfg=cfg0.copy();cfg[7]=cc[7]
            resp,temp,_=m.response(hist,i,old);resp[3:5]=sensed[:2]
            if method==13:rc.decide(y,zr,cfg,resp,temp)
            else:m.decide(y,zr,cfg,resp,temp,method,policy)
        for _ in range(3):
          for y,zr,cc,hist,i,old,sensed in contexts:
            cfg=cfg0.copy();cfg[7]=cc[7]
            tic=time.perf_counter_ns()
            if method==0:resp=np.zeros(8);temp=0.
            else:
                resp,temp,_=m.response(hist,i,old);resp[3:5]=sensed[:2]
            if method==13:rc.decide(y,zr,cfg,resp,temp)
            else:m.decide(y,zr,cfg,resp,temp,method,policy)
            records.append((time.perf_counter_ns()-tic)*1e-9)
        timings[name]=np.array(records)
        print('timing',name,'median/p95 ms',np.percentile(records,[50,95])*1000,flush=True)
    return timings


def summarize(metrics,baseline):
    rng=np.random.default_rng(8490);n=len(metrics)
    ids=rng.integers(0,n,size=(4000,n))
    benefit=baseline[:,0]-metrics[:,0]
    ratios=1-metrics[ids,0].mean(1)/baseline[ids,0].mean(1)
    return {'mean_J':float(metrics[:,0].mean()),'median_J':float(np.median(metrics[:,0])),
      'force_rms_N':float(metrics[:,1].mean()),'tangent_rms_mm':float(metrics[:,2].mean()*1000),
      'normal_pose_rms_mm':float(metrics[:,10].mean()*1000),'worst_peak_force_N':float(metrics[:,3].max()),'contact_loss_fraction':float(metrics[:,4].mean()),
      'active_energy_fraction':float(metrics[:,5].mean()),'actuator_saturation_fraction':float(metrics[:,6].mean()),
      'min_energy_margin':float(metrics[:,8].min()),'max_gain_violation':float(metrics[:,9].max()),
      'benefit_vs_fixed_ratio_of_mean_J':float(1-metrics[:,0].mean()/baseline[:,0].mean()),
      'paired_benefit_vs_fixed_95pct':np.quantile(ratios,[.025,.975]).tolist(),
      'mean_absolute_J_benefit_vs_fixed':float(benefit.mean())}


def main():
    final=json.loads((HERE/'final_policies.json').read_text());specs=comparisons(final)
    times=timing(specs);delays={name:float(np.percentile(v,95)) for name,v in times.items()}
    assert max(delays.values())<m.DECISION_DT,'Representatively measured delay exceeds decision interval; handle overruns before test.'
    source={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','evaluate.py','finalize_policy.py','learn_policy.py','run.py','reference_controller.py','reference_tuning.json','final_policies.json','force_validation.json']}
    protocol={'stage':'Frozen before test cases are generated','test_seed':8401,'noise_seed':8402,'n':48,
      'families':['homogeneous material/static surface','heterogeneous/static','heterogeneous/moving'],
      'primary_task':'Joint force and normal/tangent pose objective, exact scales in model.stage',
      'primary_initialization':'Geometry/initial force prior, shared by every method; no material truth',
      'sensitivities':['constant 1000 N/m prior with unchanged policies','force-only task with unchanged train-selected configs'],
      'methods':[name for name,_,_,_ in specs],
      'main_pairs':[['spatial_joint','analytic_inverse'],['learned_feedback','analytic_inverse'],['learned_feedback','continuous_constant'],['learned_feedback','equilibrium_reference']],
      'timing':'Warmed validation states, required response update plus decision; replay fixed method-specific p95, not WCET or full sensor/OS latency',
      'sensor_delay_s':m.SENSOR_DELAY,'gain_decision_period_s':m.DECISION_DT,'delay_replay_s':delays,
      'tuning':'Training/validation only; no further choices after these test cases are generated; bootstrap conditions on this trained policy',
      'source_sha256':source,'limitations':['Uncalibrated tissue foundation','Ideal compensated two-axis force interface',
          'Initial contact already established','No complete published-prior replication','No real robot or imaging evidence']}
    (HERE/'test_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    # The first access to independent parameter/noise seeds occurs only here.
    ps=m.sample(48,8401);ns=noises(ps,8402)
    arrays={'test_params':ps};summaries={};retained={}
    for geometry in [True,False]:
        label='informed_prior' if geometry else 'constant_prior'
        allvalues=[]
        for name,method,cfg,policy in specs:
            values=[]
            for i,(p,n) in enumerate(zip(ps,ns)):
                if method==13:val,tr=rc.episode(p,n,cfg,latency=delays[name],geometry_prior=geometry)
                else:val,tr,_,_=m.episode(p,n,method,cfg,latency=delays[name],policy=policy,geometry_prior=geometry)
                values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
                if geometry and i in [0,2] and name in ['fixed_PI','analytic_inverse','spatial_joint','learned_feedback','equilibrium_reference']:
                    retained[f'{name}_case{i}']=tr
            values=np.array(values);arrays[f'{label}_{name}']=values;allvalues.append(values)
            print(label,name,float(values[:,0].mean()),float(values[:,1].mean()),flush=True)
        summaries[label]={name:summarize(values,allvalues[0]) for (name,_,_,_),values in zip(specs,allvalues)}
        # Explicit pairwise comparison of scientific interest, not only vs the fixed rule.
        lookup={s[0]:v for s,v in zip(specs,allvalues)}
        comparisons_summary={}
        for first,second in protocol['main_pairs']:
            value=summarize(lookup[first],lookup[second])
            comparisons_summary[first+'_vs_'+second]={'relative_mean_J_benefit':value['benefit_vs_fixed_ratio_of_mean_J'],
                'paired_95pct':value['paired_benefit_vs_fixed_95pct']}
        summaries[label]['main_pair_comparisons']=comparisons_summary
    # Secondary force-only task remains visible, not substituted for the joint task.
    config=json.loads((HERE/'force_validation.json').read_text());force_summary={};base=None
    for name,method in [('fixed_PI',0),('transport_fixed',8),('spatial_gain',5),('spatial_joint',6)]:
        cfg=np.array(config['cfgs'][name]);values=[]
        for p,n in zip(ps,ns):
            val,tr,_,_=m.episode(p,n,method,cfg,geometry_prior=True,latency=delays[name])
            values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
        values=np.array(values);arrays['force_'+name]=values
        if base is None:base=values
        force_summary[name]=summarize(values,base)
        print('force only',name,force_summary[name]['mean_J'],force_summary[name]['force_rms_N'],flush=True)
    output={'scope':'Independent synthetic first-pass scans, no journal-level conclusion',
      'joint':summaries,'force_only':force_summary,'timing':{name:{'median_ms':float(np.median(t)*1000),
         'p95_ms':delays[name]*1000,'max_ms':float(t.max()*1000),'queries':len(t)} for name,t in times.items()},
      'protocol_sha256':hashlib.sha256((HERE/'test_protocol.json').read_bytes()).hexdigest(),'source_sha256':source}
    for f,h in source.items():assert hashlib.sha256((HERE/f).read_bytes()).hexdigest()==h,'Source changed during frozen test'
    (HERE/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    np.savez_compressed(HERE/'paired.npz',**arrays,**retained)
    np.savez_compressed(HERE/'timing_samples.npz',**times)
    print('independent tests complete',flush=True)
if __name__=='__main__':main()
