"""Value and compute screen for optimal impedance, with frozen test splits.

The learned map is an explicitly model-conditioned supervised variant of route
2, not PPO or an observation-only controller. Exact current model information is
shared by all methods. This tests subproblems of B0, not certified deployment.
"""
from pathlib import Path
import copy
import hashlib
import json
import time
import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.optimize import minimize
import torch
from torch import nn
import model as m

HERE=Path(__file__).resolve().parent
FAMILIES=[('white',False,False),('colored',True,False),('colored_observation_noise',True,True)]
TRAIN,VALID,TEST=320,64,96
EPOCHS=2500


def fit_actor(x,y,xv,yv,seed):
    torch.manual_seed(seed)
    torch.set_num_threads(1)
    net=nn.Sequential(nn.Linear(x.shape[1],64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,4)).double()
    mean=np.log(y).mean(axis=0);sd=np.log(y).std(axis=0);sd=np.maximum(sd,.05)
    tx=torch.tensor(x);ty=torch.tensor((np.log(y)-mean)/sd)
    tv=torch.tensor(xv);vy=torch.tensor((np.log(yv)-mean)/sd)
    optimizer=torch.optim.Adam(net.parameters(),lr=.002)
    schedule=torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,EPOCHS,eta_min=2e-5)
    bestloss=np.inf;best=None;bestepoch=0;history=[]
    gen=torch.Generator().manual_seed(seed+1000)
    for epoch in range(EPOCHS):
        idx=torch.randint(len(tx),(128,),generator=gen)
        optimizer.zero_grad();loss=torch.mean((net(tx[idx])-ty[idx])**2)
        loss.backward();optimizer.step();schedule.step()
        if epoch%50==0 or epoch==EPOCHS-1:
            with torch.no_grad():vl=float(torch.mean((net(tv)-vy)**2))
            history.append([epoch,float(loss),vl])
            if vl<bestloss:bestloss=vl;best=copy.deepcopy(net.state_dict());bestepoch=epoch
    net.load_state_dict(best);net.eval()
    torch.save(dict(state_dict=best,mean=mean.tolist(),sd=sd.tolist(),seed=seed,
                    input_dim=x.shape[1],selected_epoch=bestepoch),HERE/f'actor_{fit_actor.family}_{seed}.pt')
    def predict(xx):
        with torch.no_grad():out=net(torch.tensor(np.asarray(xx))).numpy()
        return np.clip(np.exp(out*sd+mean),m.LOW,m.HIGH)
    return predict,dict(seed=seed,selected_epoch=bestepoch,validation_log_gain_MSE=bestloss,history=history)


def main():
    torch.set_num_threads(1)
    protocol=dict(question='Rapid update of useful near-optimal impedance after environment changes',
        model='2-DOF preloaded linearized contact, 4 diagonal impedance gains, stationary H2 objective',
        train=TRAIN,validation=VALID,test=TEST,seed_base=2801,families=[f[0] for f in FAMILIES],
        gain_order=['K_normal','K_roll_normalized','D_normal','D_roll_normalized'],
        gain_low=(m.LOW*m.SCALE).tolist(),gain_high=(m.HIGH*m.SCALE).tolist(),
        masses=m.MASS.tolist(),position_cost_scale_m=.003,velocity_cost_scale_m_s=.03,
        contact_force_cost_scale_N=.7,actuation_cost_scale_N=1.,
        observation_noise_std=[50e-6,50e-6,.002,.002],observation_noise_bandwidth_Hz=30.,
        actor_seeds=[11,22,33],actor_updates=EPOCHS,teacher='L-BFGS-B; two starts on test; local numerical reference, not global proof',
        information='Exact current environment coefficients and disturbance family shared; no state-estimation claim',
        pilot='Three hand-selected points were used to check derivatives before random splits; no held-out test used to change bounds or cost.',
        timing_revision='Method evaluation order randomized after initial first-actor cache effect; model, splits and training unchanged.',
        exclusions=['finite horizon and actual gain-rate/energy dynamics','contact establishment/loss and hard-force safety',
                    'full B0 MHE/HJB/NEOC','ultrasound image quality','real robot timing'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest())
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    allresults={};arrays={};totalstart=time.perf_counter()
    for family_index,(label,colored,sensed) in enumerate(FAMILIES):
        rng=np.random.default_rng(2801+family_index);dim=4 if colored else 3
        xx=rng.uniform(-1,1,(TRAIN+VALID+TEST,dim));gains=[];statuses=[]
        teacher_start=time.perf_counter()
        for i,x in enumerate(xx):
            p=m.Problem(x,colored,sensed)
            starts=[p.analytic()]
            if i>=TRAIN+VALID:starts.append(np.ones(4))
            sol=p.solve(starts);gains.append(sol.x)
            statuses.append(float(np.max(np.abs(m.projected_gradient(sol.x,sol.jac)))))
        teacher_time=time.perf_counter()-teacher_start;gg=np.array(gains)
        xtr,xv,xt=xx[:TRAIN],xx[TRAIN:TRAIN+VALID],xx[-TEST:]
        ytr,yv,yt=gg[:TRAIN],gg[TRAIN:TRAIN+VALID],gg[-TEST:]
        problems=[m.Problem(x,colored,sensed) for x in xtr]
        def aggregate(g):
            vals=[p.evaluate(g) for p in problems]
            return float(np.mean([v[0] for v in vals])),np.mean([v[1] for v in vals],axis=0)
        t=time.perf_counter()
        fixedsol=minimize(aggregate,np.mean(ytr,axis=0),jac=True,
            bounds=list(zip(m.LOW,m.HIGH)),method='L-BFGS-B',options=dict(gtol=1e-8,ftol=1e-13,maxiter=150))
        fixed_time=time.perf_counter()-t
        t=time.perf_counter();rbf=RBFInterpolator(xtr,np.log(ytr),kernel='thin_plate_spline',degree=1,smoothing=1e-8)
        rbf_time=time.perf_counter()-t
        actors=[];training=[]
        fit_actor.family=label
        for seed in [11,22,33]:
            t=time.perf_counter();pred,meta=fit_actor(xtr,ytr,xv,yv,seed)
            meta['training_s']=time.perf_counter()-t;training.append(meta);actors.append(pred)
        methods=['fixed_tuned','analytic_schedule','teacher_RBF','actor_11','actor_22','actor_33','RBF_plus_Newton','actor_11_plus_Newton']
        values=np.zeros((TEST,len(methods)));times=np.zeros_like(values);residuals=np.zeros_like(values)
        poses=np.zeros((TEST,len(methods),4));ref=[];care=[];detailed=[]
        for i,(x,gref) in enumerate(zip(xt,yt)):
            p=m.Problem(x,colored,sensed);ref.append(p.evaluate(gref)[0])
            tc=time.perf_counter();jc,_=p.care_reference();care.append([jc,time.perf_counter()-tc])
            def predict(method):
                if method==0:return fixedsol.x
                if method==1:return p.analytic()
                if method==2:return np.clip(np.exp(rbf(x[None,:])[0]),m.LOW,m.HIGH)
                if method<6:return actors[method-3](x)
                initial=predict(2 if method==6 else 3)
                return m.newton(p,initial,maxiter=25)[0]
            for method in rng.permutation(len(methods)):
                # One-query inference, no batch-time division or setup hiding.
                t=time.perf_counter();g=predict(method);times[i,method]=time.perf_counter()-t
                value,grad=p.evaluate(g);values[i,method]=value;poses[i,method]=g
                residuals[i,method]=np.max(np.abs(m.projected_gradient(g,grad)))
            if i<3:detailed.append(dict(test_case=i,reference=p.metrics(gref),fixed=p.metrics(fixedsol.x),
                analytic=p.metrics(p.analytic()),actor=p.metrics(poses[i,3])))
        ref=np.array(ref);care=np.array(care)
        summary={}
        for i,name in enumerate(methods):
            gap=values[:,i]/ref-1
            summary[name]=dict(mean_J=float(values[:,i].mean()),
                mean_relative_gap=float(gap.mean()),p95_relative_gap=float(np.quantile(gap,.95)),max_relative_gap=float(gap.max()),
                mean_paired_improvement_over_fixed=float(np.mean(1-values[:,i]/values[:,0])),
                median_inference_ms=float(1000*np.median(times[:,i])),p95_inference_ms=float(1000*np.quantile(times[:,i],.95)),
                max_projected_gradient=float(residuals[:,i].max()))
        allresults[label]=dict(methods=summary,teacher_generation_s=teacher_time,
            teacher_max_projected_gradient=max(statuses),fixed_optimization_s=fixed_time,
            RBF_setup_s=rbf_time,actor_training=training,fixed_gains=(fixedsol.x*m.SCALE).tolist(),
            CARE=dict(mean_J=float(care[:,0].mean()),median_ms=float(1000*np.median(care[:,1])),
                mean_diagonal_reference_excess=float(np.mean(ref/care[:,0]-1)),
                status='exact broader-feedback lower bound' if not colored and not sensed else 'white-model full-state feedback, NOT a lower bound or optimum here'),
            selected_physical_metrics=detailed)
        arrays[label+'_features']=xx;arrays[label+'_teacher_gains']=gg;arrays[label+'_test_J']=values
        arrays[label+'_test_reference_J']=ref;arrays[label+'_test_times_s']=times;arrays[label+'_test_gains']=poses
        print(label,json.dumps(summary),flush=True)
    np.savez_compressed(HERE/'paired.npz',**arrays)
    (HERE/'results.json').write_text(json.dumps(dict(protocol=protocol,results=allresults,
        elapsed_s=time.perf_counter()-totalstart),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
