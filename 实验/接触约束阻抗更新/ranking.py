"""Task-cost learning instead of imitating non-unique gain targets.

Mature decision-focused finite-action learning, compared to ordinary quadratic
regression and RBF on exactly the same rollout labels. All choices fixed on
validation data, followed by a fresh 96-context test (seed 7503).
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
from scipy.interpolate import RBFInterpolator
import torch
from torch import nn
import model as m
from screen import library,critical,teacher

HERE=Path(__file__).resolve().parent


def quadratic(x):
    # Ordinary degree-two regression, a deliberately simple learned alternative.
    return np.column_stack((np.ones(len(x)),x,
        np.array([x[:,i]*x[:,j] for i in range(x.shape[1]) for j in range(i,x.shape[1])]).T))


def main():
    torch.set_num_threads(1)
    old=json.loads((HERE/'results.json').read_text());data=np.load(HERE/'paired.npz')
    fixed=np.array(old['fixed_target'])/m.SCALE;theta=np.array(old['critical_theta'])
    xm=m.features(data['train_params'],data['train_states']).mean(0)
    xs=m.features(data['train_params'],data['train_states']).std(0)
    x=(m.features(data['train_params'],data['train_states'])-xm)/xs
    xv=(m.features(data['validation_params'],data['validation_states'])-xm)/xs
    # Shift/scale per-context costs; preserves the candidate ordering.
    trainj=data['train_library_values'];valj=data['validation_library_values']
    y=trainj/trainj.min(1)[:,None]-1
    yv=valj/valj.min(1)[:,None]-1
    ymean=y.mean(0);yscale=np.maximum(y.std(0),1e-4)
    def regret(pred):return float(np.mean(yv[np.arange(len(yv)),np.argmin(pred,axis=1)]))
    fits={};metadata={}
    t0=time.perf_counter();candidates=[]
    xx=quadratic(x);vv=quadratic(xv)
    for lam in [.001,.01,.1,1.,10.]:
        weights=np.linalg.solve(xx.T@xx+lam*np.eye(xx.shape[1]),xx.T@y)
        candidates.append((regret(vv@weights),lam,weights))
    score,lam,weights=min(candidates,key=lambda v:v[0]);fits['quadratic']=weights
    metadata['quadratic']=dict(regularization=lam,validation_relative_regret=score,fit_seconds=time.perf_counter()-t0)
    t0=time.perf_counter();candidates=[]
    for smoothing in [.001,.01,.1,1.,10.]:
        rbf=RBFInterpolator(x,y,kernel='thin_plate_spline',smoothing=smoothing)
        candidates.append((regret(rbf(xv)),smoothing,rbf))
    score,smoothing,rbf=min(candidates,key=lambda v:v[0]);fits['rbf_rank']=rbf
    metadata['rbf_rank']=dict(smoothing=smoothing,validation_relative_regret=score,fit_seconds=time.perf_counter()-t0)
    nets=[]
    for seed in [11,22,33]:
        torch.manual_seed(seed);tic=time.perf_counter()
        net=nn.Sequential(nn.Linear(x.shape[1],64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,y.shape[1])).double()
        opt=torch.optim.Adam(net.parameters(),lr=.001,weight_decay=1e-4)
        tx=torch.as_tensor(x,dtype=torch.float64);ty=torch.as_tensor(y,dtype=torch.float64)
        tv=torch.as_tensor(xv,dtype=torch.float64)
        best=float('inf');bestweights=None
        for step in range(2400):
            ids=torch.randint(len(x),(96,));pred=net(tx[ids])
            # Expected physical regret under a soft candidate policy. No gain MSE.
            loss=(torch.softmax(pred,dim=1)*ty[ids]).sum(1).mean()
            opt.zero_grad();loss.backward();opt.step()
            if step%30==0:
                with torch.no_grad():which=net(tv).argmax(1).numpy()
                score=float(np.mean(yv[np.arange(len(yv)),which]))
                if score<best:best=score;bestweights={k:v.detach().clone() for k,v in net.state_dict().items()}
        net.load_state_dict(bestweights);net.eval();nets.append(net)
        metadata[f'cost_actor_{seed}']=dict(validation_relative_regret=best,fit_seconds=time.perf_counter()-tic)
        torch.save(dict(weights=bestweights,mean=xm,std=xs),HERE/f'cost_actor_{seed}.pt')
    np.savez_compressed(HERE/'ranker_weights.npz',quadratic=weights,mean=xm,std=xs,rbf_x=x,rbf_y=y,
                        rbf_smoothing=np.array(smoothing))
    # Fresh data generated only after all fits/hyperparameter choices above.
    params,states=m.sample(96,7503)
    features=(m.features(params,states)-xm)/xs
    protocol=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
        train_contexts=192,validation_contexts=48,new_test_contexts=96,new_test_seed=7503,
        labels='Same 30 executed-library physical costs for all rankers; no exact optimizer labels needed',
        scope='A post-screen method iteration; fresh test contexts, not independent physical-task or full-system validation',
        timing='Random-order, warmed single-query timing incl feature construction, candidate construction and selection',
        decision_rules='Hyperparameters/checkpoints selected only by validation task regret, not action error')
    (HERE/'ranking_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    names=['fixed','critical','library','teacher','quadratic','rbf_rank','cost_actor_11','cost_actor_22','cost_actor_33']
    arrays={'params':params,'states':states};plans={name:[] for name in names};times={name:[] for name in names};metas=[]
    for i,(p,s) in enumerate(zip(params,states)):
        candidates=library(s,p,fixed,theta)
        # Compilation/cache warmup, including the candidate scorer.
        m.candidate_cost(candidates,s,p)
        with torch.no_grad():
            for net in nets:net(torch.as_tensor(features[i:i+1],dtype=torch.float64))
        for name in np.random.default_rng(7600+i).permutation(names):
            tic=time.perf_counter()
            if name=='fixed':plan=np.tile(fixed,(2,1))
            elif name=='critical':plan=np.tile(critical(p[None],theta)[0],(2,1))
            else:
                xx=(m.features(p[None],s[None])-xm)/xs
                candidates=library(s,p,fixed,theta)
                if name in ['library','teacher']:
                    costs=m.candidate_cost(candidates,s,p)
                    plan=candidates[int(np.argmin(costs))]
                    if name=='teacher':plan,meta=teacher(s,p,plan);metas.append(dict(index=i,**meta))
                elif name=='quadratic':plan=candidates[int(np.argmin(quadratic(xx)@weights))]
                elif name=='rbf_rank':plan=candidates[int(np.argmin(rbf(xx)))]
                else:
                    net=nets[[11,22,33].index(int(name.split('_')[-1]))]
                    with torch.no_grad():choice=int(net(torch.as_tensor(xx,dtype=torch.float64)).argmax())
                    plan=candidates[choice]
            times[name].append((time.perf_counter()-tic)*1000);plans[name].append(plan)
        if i%24==0:print('fresh',i,flush=True)
    rows=[]
    for name in names:
        js=[];delayed=[];diag=[]
        for plan,s,p,ms in zip(plans[name],states,params,times[name]):
            j,d,_=m.rollout(plan,s,p);js.append(j);diag.append(d)
            delayed.append(m.rollout(plan,s,p,delay=min(ms/1000,.3))[0])
        arrays[name+'_targets']=np.array(plans[name]);arrays[name+'_ms']=np.array(times[name])
        arrays[name+'_J']=np.array(js);arrays[name+'_delayed_J']=np.array(delayed);arrays[name+'_diagnostics']=np.array(diag)
    for name in names:
        j=arrays[name+'_J'];jd=arrays[name+'_delayed_J'];ms=arrays[name+'_ms'];d=arrays[name+'_diagnostics']
        rows.append(dict(method=name,mean_J=float(j.mean()),
            mean_improvement_over_critical_pct=float(np.mean(1-j/arrays['critical_J'])*100),
            measured_release_improvement_over_critical_pct=float(np.mean(1-jd/arrays['critical_delayed_J'])*100),
            gap_to_local_teacher_pct=float(np.mean(j/arrays['teacher_J']-1)*100),
            inference_median_ms=float(np.median(ms)),inference_p95_ms=float(np.quantile(ms,.95)),
            worst_relative_loss_to_critical_pct=float(np.max(j/arrays['critical_J']-1)*100),
            tank_active_cases=int(np.sum(d[:,5]>0)),maximum_force_N=float(np.max(d[:,1]))))
    out=dict(protocol=protocol,fit=metadata,methods=rows,teachers=metas)
    np.savez_compressed(HERE/'ranking_paired.npz',**arrays)
    (HERE/'ranking_results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':main()
