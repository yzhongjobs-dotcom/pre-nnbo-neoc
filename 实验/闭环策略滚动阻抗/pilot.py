"""Train-only pilot for policy rollout; independent evaluation is a later stage."""
from pathlib import Path
import hashlib,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize
import controller as c
import model as m
import coupled_contact as cm
from learn_policy import pad
from run import noises

HERE=Path(__file__).resolve().parent


@njit(cache=True)
def loss(theta,ps,ns,counts,cfg):
    bank=np.empty((1,8));bank[0]=theta;scales=np.ones((1,2));value=0.
    for i in range(len(ps)):
        value+=c.episode(ps[i],ns[i,:counts[i]],cfg,bank,scales,policy_index=0,dt=.001,nq=33)[0][0]/len(ps)
    return value


def context(p,n,cfg,theta,where=30):
    metric,trace=cm.episode(p,n,12,cfg,theta)
    i=min(where,len(trace)-1);hist=np.zeros((i+1,8));hist[:,:4]=trace[:i+1,1:5]
    hist[:,4:7]=trace[:i+1,11:14]+n[:i+1];hist[:,7]=trace[:i+1,0]*p[9]+hist[:,0]
    zr=trace[0,2];kk=1.5*3/zr if zr<m.A*m.A/(2*m.R) else 3/(zr-m.A*m.A/(6*m.R))
    old=np.array([min(7000.,max(200.,kk)),0.,0.,3.,.6,zr,0.,0.]);sense=np.array([3.,.6,0.])
    for j in range(i+1):
        resp,temp,_=m.response(hist,j,old);old=resp.copy()
        sense[:2]=cfg[6]*resp[3:5]+(1-cfg[6])*sense[:2];resp[3:5]=sense[:2]
    held=np.r_[trace[max(0,i-1),17:21],trace[max(0,i-1),21]]
    cc=cfg.copy();cc[7]=p[9]
    return trace[i,1:11],zr,cc,hist,i,resp,sense,temp,held


def measured(ctx,bank,scales,horizon,ensemble,consistent):
    y,zr,cfg,hist,i,resp,sense,temp,held=ctx
    # JIT compilation is excluded; observer and required decision included.
    c.decide(y,zr,cfg,hist,i,resp,sense,temp,bank,scales,horizon,ensemble,consistent,0.,held)
    samples=[]
    for _ in range(12):
        start=time.perf_counter_ns()
        rr,tt,_=m.response(hist,i,resp)
        rr[3:5]=sense[:2]
        c.decide(y,zr,cfg,hist,i,rr,sense,tt,bank,scales,horizon,ensemble,consistent,0.,held)
        samples.append((time.perf_counter_ns()-start)/1e6)
    return {'median_ms':float(np.median(samples)),'p95_ms':float(np.quantile(samples,.95)),'samples_ms':samples}


def main():
    start=time.perf_counter();old=json.loads((c.BASE/'final_policies.json').read_text())
    pure=json.loads((c.BASE/'pure_gain_fit.json').read_text());cfg=np.array(old['cfg'])
    ps=m.sample(18,8101);ns,counts=pad(noises(ps,8102))
    policies={'constant_gain':np.array(pure['policies']['continuous_gain_only']),
              'analytic_joint':np.array(old['policies']['analytic_inverse']),
              'learned_joint':np.array(old['policies']['learned_feedback']),
              'learned_gain':np.array(pure['policies']['learned_gain_only'])}
    fits={}
    for name,ids in [('constant_gain',[0,3]),('analytic_joint',[0,3,5])]:
        initial=policies[name].copy();calls=[0]
        def fun(x):
            theta=initial.copy();theta[ids]=x;calls[0]+=1
            return loss(theta,ps,ns,counts,cfg)
        bounds=[(0.,np.log(800)) if j==0 else (np.log(2),np.log(150)) if j==3 else (-1.,2.) for j in ids]
        fitted=minimize(fun,initial[ids],method='Powell',bounds=bounds,
                        options={'maxfev':35,'maxiter':3,'xtol':.02,'ftol':1e-5})
        candidate=initial.copy();candidate[ids]=fitted.x
        # Retain a train-selected old fit if finite-budget optimization regresses.
        baseline=loss(initial,ps,ns,counts,cfg)
        if fitted.fun<baseline:policies[name]=candidate
        fits[name]={'calls':calls[0],'old_train_J':float(baseline),'new_train_J':float(fitted.fun),
                    'optimizer_success':bool(fitted.success),'message':str(fitted.message)}
        print('fit',name,fits[name],flush=True)
    names=['constant_gain','analytic_joint','learned_joint','learned_gain']
    bank=np.array([policies[k] for k in names]);scales=np.ones((4,2))
    # Small feedback-policy family, rather than a larger target-gain grid.
    for scale in [[.5,1.],[2.,1.],[1.,.5],[1.,2.]]:
        bank=np.vstack([bank,policies['learned_joint']]);scales=np.vstack([scales,scale]);names.append('learned_scaled_'+str(scale))
    pp=ps[:6];pn=noises(ps,8102)[:6]
    ctx=context(pp[0],pn[0],cfg,policies['constant_gain'],30)
    variants=[('nominal_raw',.16,False,False),('nominal_closed',.16,False,True),
              ('ensemble_closed',.16,True,True),('ensemble_long',.32,True,True)]
    timings={};table={};arrays={'pilot_params':pp}
    for name,horizon,ensemble,consistent in variants:
        ti=measured(ctx,bank,scales,horizon,ensemble,consistent);timings[name]=ti
        latency=ti['p95_ms']/1000;values=[]
        print('timing',name,ti['median_ms'],ti['p95_ms'],flush=True)
        for p,n in zip(pp,pn):
            val,tr=c.episode(p,n,cfg,bank,scales,horizon,ensemble,consistent,latency)
            values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
        a=np.array(values);arrays[name]=a;table[name]={'mean_J':float(a[:,0].mean()),'force_rms_N':float(a[:,1].mean()),
             'normal_pose_rms_mm':float(a[:,10].mean()*1000)}
        print('pilot',name,table[name],flush=True)
    for index,name in enumerate(['constant_gain','analytic_joint','learned_joint','learned_gain']):
        values=[]
        for p,n in zip(pp,pn):
            val,tr=c.episode(p,n,cfg,bank,scales,policy_index=index,latency=.000003)
            values.append(np.r_[val,np.sqrt(np.mean((tr[:,2]-tr[0,2])**2))])
        a=np.array(values);arrays[name]=a;table[name]={'mean_J':float(a[:,0].mean()),'force_rms_N':float(a[:,1].mean()),
             'normal_pose_rms_mm':float(a[:,10].mean()*1000)}
        print('pilot baseline',name,table[name],flush=True)
    out={'stage':'Train-only mechanism pilot; NOT independent evidence or new superiority claim',
         'train_seed':8101,'noise_seed':8102,'n_train':18,'pilot_cases':'First six of training set, two per family',
         'cfg':cfg.tolist(),'policies':{k:v.tolist() for k,v in policies.items()},'bank':bank.tolist(),
         'scales':scales.tolist(),'bank_names':names,'fits':fits,'variants':variants,'timing':timings,
         'table':table,'seconds':time.perf_counter()-start,
         'scope':'Same corrected reduced plant, same task; local response weights are heuristic, not calibrated posterior; no full B0/PPO replication',
         'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['controller.py','pilot.py']},
         'dependencies_sha256':{f:hashlib.sha256((c.BASE/f).read_bytes()).hexdigest() for f in ['model.py','coupled_contact.py','final_policies.json','pure_gain_fit.json','run.py']}}
    (HERE/'pilot.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(HERE/'pilot.npz',**arrays)
    print('pilot complete',flush=True)


if __name__=='__main__':main()
