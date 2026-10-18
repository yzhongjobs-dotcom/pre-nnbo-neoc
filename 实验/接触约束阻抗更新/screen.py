"""Exploratory known-model contact screen. Freeze splits before fitting.

Strong alternatives include tuned fixed targets, tuned environment-conditioned
critical damping and a small rollout library. Teachers are local numerical
solutions; none is a global optimum or a full B0 controller.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc
from scipy.interpolate import RBFInterpolator
import torch
from torch import nn
import model as m

HERE = Path(__file__).resolve().parent


def critical(params, theta):
    kt = np.exp(theta[0]);kn = np.exp(theta[1])*params[:, 0]
    dt = 2*np.exp(theta[2])*np.sqrt(m.MASS[0]*kt)*np.ones(len(params))
    dn = 2*np.exp(theta[3])*np.sqrt(m.MASS[1]*(kn+params[:, 0]))
    return np.clip(np.column_stack((kt*np.ones(len(params)),kn,dt,dn))/m.SCALE,m.LOW,m.HIGH)


def fit_baselines(params, states):
    tic = time.perf_counter()
    def objective(z):
        g = m.normalized_targets(z)
        targets = np.tile(g[None], (len(states),2,1))
        return float(m.batch_cost(targets,states,params).mean())
    starts = np.vstack((m.target_coordinates(np.array([300.,100.,25.,15.])/m.SCALE),
                        qmc.Sobol(4,scramble=True,seed=7301).random_base2(5)))
    values = np.array([objective(z) for z in starts])
    sols = [minimize(objective,starts[j],method='L-BFGS-B',bounds=[(0,1)]*4,
                     options=dict(maxiter=50,ftol=1e-10,eps=2e-5)) for j in np.argsort(values)[:2]]
    fixed = m.normalized_targets(min(sols,key=lambda s:s.fun).x)[0]
    fixed_seconds = time.perf_counter()-tic
    def co(theta):
        g = critical(params,theta)
        return float(m.batch_cost(np.repeat(g[:,None,:],2,axis=1),states,params).mean())
    bounds = [(np.log(5.),np.log(1600.)),(np.log(.005),np.log(2.)),
              (np.log(.2),np.log(3.)),(np.log(.2),np.log(3.))]
    tic=time.perf_counter()
    starts = [np.log([300.,.1,1.,1.]), np.log([100.,.01,.5,.5])]
    sols = [minimize(co,z,method='L-BFGS-B',bounds=bounds,
                     options=dict(maxiter=60,ftol=1e-10,eps=2e-5)) for z in starts]
    best = min(sols,key=lambda s:s.fun)
    return fixed,best.x,dict(fixed_seconds=fixed_seconds,critical_seconds=time.perf_counter()-tic,
        fixed_train_J=objective(m.target_coordinates(fixed)),critical_train_J=float(best.fun))


def library(state,param,fixed,theta):
    # Physically interpretable constant targets, no access to teacher solutions.
    candidates=[state[4:8],fixed,critical(param[None],theta)[0]]
    for kt in [40.,300.,1200.]:
        for fraction in [.005,.05,.5]:
            for zeta in [.5,1.,2.]:
                kn = fraction*param[0]
                g=np.array([kt,kn,2*zeta*np.sqrt(m.MASS[0]*kt),
                            2*zeta*np.sqrt(m.MASS[1]*(kn+param[0]))])/m.SCALE
                candidates.append(np.clip(g,m.LOW,m.HIGH))
    return np.repeat(np.asarray(candidates)[:,None,:],2,axis=1)


def teacher(state,param,start):
    def cost(z):
        return float(m.rollout(m.normalized_targets(z),state,param)[0])
    tic=time.perf_counter()
    z0=m.target_coordinates(start).ravel()
    sol=minimize(cost,z0,method='L-BFGS-B',bounds=[(0,1)]*8,
                 options=dict(maxiter=65,ftol=1e-10,eps=2e-5,maxls=25))
    targets = m.normalized_targets(sol.x) if sol.fun <= cost(z0) else start
    return targets,dict(seconds=time.perf_counter()-tic,success=bool(sol.success),
                        iterations=int(sol.nit),evaluations=int(sol.nfev),message=str(sol.message))


def fit_actor(x,y,xv,yv,seed):
    torch.set_num_threads(1);torch.manual_seed(seed)
    mean=x.mean(0);std=np.maximum(x.std(0),1e-8)
    xs=torch.as_tensor((x-mean)/std,dtype=torch.float64)
    ys=torch.as_tensor(y,dtype=torch.float64)
    vs=torch.as_tensor((xv-mean)/std,dtype=torch.float64)
    vy=torch.as_tensor(yv,dtype=torch.float64)
    net=nn.Sequential(nn.Linear(x.shape[1],64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,8),nn.Sigmoid()).double()
    opt=torch.optim.Adam(net.parameters(),lr=.002)
    tic=time.perf_counter();best=float('inf');weights=None
    for step in range(1800):
        ids=torch.randint(len(xs),(min(96,len(xs)),))
        loss=((net(xs[ids])-ys[ids])**2).mean()
        opt.zero_grad();loss.backward();opt.step()
        if step%30==0:
            with torch.no_grad():v=float(((net(vs)-vy)**2).mean())
            if v<best:best=v;weights={k:v.detach().clone() for k,v in net.state_dict().items()}
    net.load_state_dict(weights);net.eval()
    torch.save(dict(weights=weights,mean=mean,std=std),HERE/f'actor_{seed}.pt')
    return net,mean,std,dict(seed=seed,training_seconds=time.perf_counter()-tic,validation_action_MSE=best)


def main():
    protocol=dict(status='Exploratory new carrier; splits fixed before fitting, not a top-journal validation',
        splits={'train':[192,7201],'validation':[48,7202],'test':[96,7203]},
        task='Tangential/normal motion after a local environment-model revision; 3/2 unilateral normal spring and smooth sliding friction',
        information='All methods receive the same exact current state, actual gains, tank and updated local model; no estimation advantage',
        horizon_s=m.HORIZON,dt_s=m.DT,force_target_N=m.FD,sliding_speed_m_s=m.SPEED,
        slope_range_N_m=[300,1500],old_to_new_slope_ratio=[.75,1.25],friction_range=[.1,.6],
        force_diagnostic_range_N=[.2,6.],actuator_clip_N=12.,
        constraints='Common continuous-feedback gain/rate/tank filter and actuator saturation; force thresholds are diagnostics, NOT enforced hard safety',
        methods=['hold actual gains','training-tuned fixed target','training-tuned critical damping schedule',
                 '30-constant-target rollout library','two-block local numerical teacher','same-label RBF','3 supervised state-conditioned actors','actor plus library rollout selection'],
        timings='Single-query warm runtime incl target generation; teacher cost includes optimization, library all predictions; independent rollout evaluation outside inference',
        exclusions=['cylindrical Hertz footprint and wrench','MHE and unknown local parameters','full robot','ultrasound images','stochastic disturbances','hard force guarantee','global optimality'],
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'model.py',Path(__file__)]})
    (HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    splits={name:m.sample(n,seed) for name,(n,seed) in protocol['splits'].items()}
    m.rollout(np.tile(splits['train'][1][0,4:8],(2,1)),splits['train'][1][0],splits['train'][0][0])
    params,states=splits['train'];fixed,theta,tuning=fit_baselines(params,states)
    print('baseline fits',fixed*m.SCALE,theta,tuning,flush=True)
    arrays={};teacher_meta=[];inference={};targets={}
    for split,(params,states) in splits.items():
        labs=[];libplans=[];libtime=[];teachtime=[];allplans=[];allvalues=[]
        for i,(p,s) in enumerate(zip(params,states)):
            tic=time.perf_counter();plans=library(s,p,fixed,theta)
            values=m.candidate_cost(plans,s,p);best=plans[int(np.argmin(values))]
            libtime.append((time.perf_counter()-tic)*1000);libplans.append(best)
            allplans.append(plans);allvalues.append(values)
            label,meta=teacher(s,p,best);labs.append(label);teachtime.append(meta['seconds']*1000+libtime[-1])
            teacher_meta.append(dict(split=split,index=i,**meta))
            if i%32==0:print(split,i,meta['seconds'],flush=True)
        arrays[split+'_params']=params;arrays[split+'_states']=states;arrays[split+'_teacher']=np.array(labs)
        arrays[split+'_library_targets']=np.array(allplans);arrays[split+'_library_values']=np.array(allvalues)
        if split=='test':
            targets['library']=np.array(libplans);targets['teacher']=np.array(labs)
            inference['library']=libtime;inference['teacher']=teachtime
    np.savez_compressed(HERE/'paired.npz',**arrays)
    x=m.features(arrays['train_params'],arrays['train_states'])
    xv=m.features(arrays['validation_params'],arrays['validation_states'])
    xt=m.features(arrays['test_params'],arrays['test_states'])
    y=m.target_coordinates(arrays['train_teacher']).reshape(-1,8)
    yv=m.target_coordinates(arrays['validation_teacher']).reshape(-1,8)
    xm=x.mean(0);xs=x.std(0)
    tic=time.perf_counter();rbf=RBFInterpolator((x-xm)/xs,y,kernel='thin_plate_spline',smoothing=.01)
    rbf_setup=time.perf_counter()-tic
    trainmeta=[];nets=[]
    for seed in [11,22,33]:
        net,mean,std,meta=fit_actor(x,y,xv,yv,seed);trainmeta.append(meta);nets.append((net,mean,std))
    params,states=splits['test']
    for key in ['hold','fixed','critical','rbf','actor_11','actor_22','actor_33','actor_library']:
        targets[key]=[];inference[key]=[]
    for i,(p,s) in enumerate(zip(params,states)):
        # Warm all actor paths; inference ordering randomized to reduce order bias.
        with torch.no_grad():
            for net,mean,std in nets:net(torch.as_tensor((xt[i:i+1]-mean)/std,dtype=torch.float64))
        for name in np.random.default_rng(8000+i).permutation(['hold','fixed','critical','rbf','actor_11','actor_22','actor_33','actor_library']):
            tic=time.perf_counter()
            if name=='hold':g=np.tile(s[4:8],(2,1))
            elif name=='fixed':g=np.tile(fixed,(2,1))
            elif name=='critical':g=np.tile(critical(p[None],theta)[0],(2,1))
            elif name=='rbf':g=m.normalized_targets(np.clip(rbf((xt[i:i+1]-xm)/xs)[0],0,1))
            else:
                net,mean,std=nets[0 if name=='actor_library' else [11,22,33].index(int(name.split('_')[1]))]
                with torch.no_grad():out=net(torch.as_tensor((xt[i:i+1]-mean)/std,dtype=torch.float64)).numpy()[0]
                g=m.normalized_targets(out)
                if name=='actor_library':
                    plans=np.concatenate((library(s,p,fixed,theta),g[None]))
                    g=plans[int(np.argmin(m.candidate_cost(plans,s,p)))]
            elapsed=(time.perf_counter()-tic)*1000
            targets[name].append(g);inference[name].append(elapsed)
    rows=[]
    for name,plans in targets.items():
        plans=np.asarray(plans);cost=[];diag=[];delayed=[]
        for g,s,p,ms in zip(plans,states,params,inference[name]):
            j,d,_=m.rollout(g,s,p);cost.append(j);diag.append(d)
            delayed.append(m.rollout(g,s,p,delay=min(ms/1000,m.HORIZON))[0])
        arrays['test_'+name+'_targets']=plans;arrays['test_'+name+'_J']=np.array(cost)
        arrays['test_'+name+'_delayed_J']=np.array(delayed)
        arrays['test_'+name+'_diagnostics']=np.array(diag);arrays['test_'+name+'_ms']=np.array(inference[name])
    ref=arrays['test_teacher_J']
    for name in targets:
        j=arrays['test_'+name+'_J'];jd=arrays['test_'+name+'_delayed_J'];d=arrays['test_'+name+'_diagnostics'];ms=np.array(inference[name])
        rows.append(dict(method=name,mean_J=float(j.mean()),mean_gap_to_local_teacher_pct=float(np.mean(j/ref-1)*100),
            mean_improvement_over_fixed_pct=float(np.mean(1-j/arrays['test_fixed_J'])*100),
            measured_release_mean_J=float(jd.mean()),
            measured_release_improvement_over_fixed_pct=float(np.mean(1-jd/arrays['test_fixed_delayed_J'])*100),
            inference_median_ms=float(np.median(ms)),inference_p95_ms=float(np.quantile(ms,.95)),
            force_below_point2_cases=int(np.sum(d[:,0]<.2)),force_above6_cases=int(np.sum(d[:,1]>6.)),
            maximum_force_N=float(d[:,1].max()),minimum_energy_margin_J=float(d[:,4].min()),
            tank_filter_active_cases=int(np.sum(d[:,5]>0)),max_gain_violation=float(d[:,7].max()),
            raw_actuator_saturation_cases=int(np.sum(d[:,6]>12.))))
    np.savez_compressed(HERE/'paired.npz',**arrays)
    result=dict(protocol=protocol,tuning=tuning,fixed_target=(fixed*m.SCALE).tolist(),critical_theta=theta.tolist(),
        teachers=teacher_meta,learning=trainmeta,rbf_setup_seconds=rbf_setup,methods=rows)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':main()
