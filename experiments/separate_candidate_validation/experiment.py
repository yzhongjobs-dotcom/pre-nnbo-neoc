"""Separate candidate results in a common causal scanning carrier.

Named paper baselines are marked task adaptations, not full replications.
This runner does not overwrite the previous engineering-precheck data.
"""
import argparse,json,time,sys
import casadi as ca
from pathlib import Path
import numpy as np
import torch
from optimization import Problem,ReusedOptimizer,FastNeighbor,context
from core import m,plant,features,encode,decode,LOGLO,LOGHI,LO,HI,local_contact
from frontend import Frontend
from run import Actor,save_contexts
from baselines import NNBO,OneStepVIC,calibrate
from rate_core import advance_rate,decode_rate
from compatible import CompatibleMPC

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'mhe_guided_impedance_pilot'
torch.set_num_threads(1)


def raw_features(y,zr,speed,cfg,sensed):
    x=np.zeros(23);x[:4]=y[:4]/[.0015,.01,.02,.03];x[4]=(y[1]-zr)/.0015
    x[5:9]=(np.log(np.clip(y[4:8],LO,HI))-LOGLO)/(LOGHI-LOGLO);x[9]=y[8]/.03
    x[17:20]=(sensed-[3.,.8,0])/[1.,1.,.03];x[21]=speed/.012;x[22]=cfg[8]
    return x


def episode(p,seed,method,cfg,problem=None,actor=None,steps=600,timed=False,dt=.001,omega=10.,stride=None,qp_reg=1e-7):
    cfg=cfg.copy();cfg[7]=p[9];y=plant.initial(p,cfg[:4],65);zr=y[1]
    front=Frontend(zr,p[9]);noise=np.random.default_rng(seed).normal(size=(steps,3))*[.03,.03,.0001]
    action_rng=np.random.default_rng(seed+1212)
    hc,cal=calibrate(p,zr,seed+777);controller=None
    if method=='nnbo_adapted':controller=NNBO(hc,zr,cfg,omega)
    if stride is None:stride=1 if method in ('nnbo_adapted','mpvic_os_adapted','candidate2','static_vic') else 8
    if method=='mpvic_os_adapted':controller=OneStepVIC(cfg,horizon=.01*stride,regularizer=qp_reg)
    if method in ('candidate1','ordinary_warm','full_mpc'):controller=ReusedOptimizer(problem,reuse=method=='candidate1',full=method=='full_mpc')
    if method=='candidate1_fast':controller=FastNeighbor(problem)
    if method=='rate_neighbor':controller=FastNeighbor(problem)
    if method=='rate_mpc':controller=ReusedOptimizer(problem,reuse=False,full=True)
    if method=='rate_compatible':controller=CompatibleMPC(problem)
    if method=='rate_bestfit':controller=CompatibleMPC(problem,best_only=True)
    rate_mode=method.startswith('rate_') or (method=='candidate2' and getattr(actor,'rate_actions',False))
    needs_front=method not in ('nnbo_adapted','mpvic_os_adapted','fixed_aux','random_gains') and not (method=='candidate2' and getattr(actor,'no_front',False))
    target=cfg[:4].copy();ff=0.;sensed=np.array([3.,.8,0.]);pending=None;release=0.
    if rate_mode:target=np.zeros(4)
    advance=advance_rate if rate_mode else plant.advance
    metrics=np.zeros(8);trace=[];decision_ms=[];release_delay=[];state_min=np.inf;gain_violation=0.;skipped=0
    for i in range(steps):
        t=i*.01;force=plant.contact(t,y,p,65)[0];front.observe(t,y,force+noise[i]);tic=time.perf_counter()
        q,anchor=front.update(i);mhe_time=time.perf_counter()-tic
        sensed=cfg[6]*front.forces[max(0,i-2)]+(1-cfg[6])*sensed
        if pending is not None and t>=release-1e-12:target=pending[:4].copy();ff=pending[4];pending=None
        if i%stride==0:
            if pending is not None:skipped+=1
            else:
                tic=time.perf_counter()
                if method=='fixed_aux':proposal=np.r_[cfg[:4],0.]
                elif method=='random_gains':
                    gg=np.exp(action_rng.uniform(np.log([100.,1.,5.,2.]),np.log([1600.,800.,150.,100.])))
                    if action_rng.random()<.5:gg[1]=action_rng.uniform(1.,800.)
                    proposal=np.r_[gg,0.]
                elif method=='nnbo_adapted':proposal=controller.decide(t,y,sensed,dt=.01*stride)
                elif method=='mpvic_os_adapted':proposal=controller.decide(t,y,zr,sensed)
                elif method=='static_vic':
                    # Static linear-contact minimizer: K/(1+kp)=q_z/(q_F*k).
                    # Principle-based variable-gain comparator, not a paper
                    # replication. Uses the same causal contact-slope estimate.
                    slope=max(200.,local_contact(t,y,q,anchor,p[9])[1])
                    proposal=np.r_[cfg[:4],0.];proposal[1]=np.clip((1+cfg[4])*cfg[8]*(.3/.0015)**2/slope,1.,800.)
                elif method in ('rate_compatible','rate_bestfit'):proposal=controller.decide(t,y,zr,sensed,cfg,front,i)
                elif method=='candidate2':
                    x=raw_features(y,zr,p[9],cfg,sensed) if getattr(actor,'no_front',False) else features(t,y,zr,q,anchor,p[9],cfg)
                    with torch.no_grad():a=actor(torch.tensor(x,dtype=torch.float32)[None]).numpy()[0]
                    proposal=decode_rate(a[:5]) if rate_mode else decode(a[:5])
                else:
                    c=context(t,y,zr,q,anchor,sensed,cfg);initial=np.tile(encode(np.r_[cfg[:4],0.]),2)
                    if rate_mode:initial=np.zeros(10)
                    proposal=controller.decide(c,initial)
                elapsed=time.perf_counter()-tic;decision_ms.append(1000*elapsed)
                delay=elapsed+(mhe_time if needs_front else 0.) if timed else 0.
                pending=proposal;release=t+delay;release_delay.append(delay)
                if delay==0.:target=pending[:4].copy();ff=pending[4];pending=None
        old=y.copy()
        if pending is not None and t<release<t+.01:
            d=release-t;y,aa=advance(t,y,p,target,ff,zr,sensed,cfg,d,dt,65)
            target=pending[:4].copy();ff=pending[4];pending=None
            y,bb=advance(t+d,y,p,target,ff,zr,sensed,cfg,.01-d,dt,65)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:y,aa=advance(t,y,p,target,ff,zr,sensed,cfg,.01,dt,65)
        metrics[:3]+=aa[:3];metrics[3]=max(metrics[3],aa[3]);metrics[4:7]+=aa[4:7];metrics[7]=max(metrics[7],aa[7])
        state_min=min(state_min,y[8]-m.E_MIN);gain_violation=max(gain_violation,float(np.max(m.GMIN-y[4:8])),float(np.max(y[4:8]-m.GMAX)))
        trace.append(np.r_[t,old,force,q,anchor,target,ff,sensed])
    duration=steps*.01;metrics[0]/=duration;metrics[1:3]=np.sqrt(metrics[1:3]/duration);metrics[4:7]/=duration
    extra={}
    if method=='nnbo_adapted':extra={'reference_failures':controller.failures,'final_critic':controller.w.tolist(),'hc':hc.tolist()}
    if isinstance(controller,(ReusedOptimizer,FastNeighbor,CompatibleMPC)):extra['optimizer']=controller.log
    out={'metrics':metrics.tolist(),'decision_ms':decision_ms,'mhe_ms':(np.array(front.timings)*1000).tolist(),'uses_frontend':needs_front,'publication_delay_s':release_delay,'energy_margin_min':state_min,'gain_violation':gain_violation,'skipped_updates':skipped,'calibration':{**cal,'used':method=='nnbo_adapted'},**extra}
    return out,np.array(trace)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--methods',default='nnbo_adapted,mpvic_os_adapted,candidate1,candidate2');ap.add_argument('--n',type=int,default=6);ap.add_argument('--seed',type=int,default=9201);ap.add_argument('--name',default='validation_v1');ap.add_argument('--actor',default='candidate2_v1.pt');ap.add_argument('--steps',type=int,default=600);ap.add_argument('--timed',action='store_true');ap.add_argument('--omega',type=float,default=10.);ap.add_argument('--dt',type=float,default=.001);ap.add_argument('--stride',type=int);ap.add_argument('--qp-reg',type=float,default=1e-7);ap.add_argument('--weight',type=float,default=1.);ap.add_argument('--compiled');ap.add_argument('--collect',action='store_true');args=ap.parse_args()
    if (ROOT/f'{args.name}.json').exists():raise FileExistsError('Choose a new run name; results are immutable.')
    cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);cfg[8]=args.weight;actor=None
    if 'candidate2' in args.methods:
        ck=torch.load(ROOT/args.actor,weights_only=True)
        if ck.get('kind','').startswith('guided'):
            from learning import GuidedActor
            actor=GuidedActor(ck['input'],ck['gains_only'],ck.get('task_log',False),ck.get('no_front',False))
        else:actor=Actor(ck['input'])
        actor.load_state_dict(ck['state']);actor.eval();actor.rate_actions=ck.get('rate_actions',False)
    if 'rate_' in args.methods and any(k in args.methods for k in ('candidate1','ordinary_warm','full_mpc')):raise ValueError('Run rate and target optimizer experiments separately.')
    problem=Problem(rate_actions='rate_' in args.methods,gains_only='rate_' in args.methods) if 'rate_' in args.methods or any(k in args.methods for k in ('candidate1','ordinary_warm','full_mpc')) else None
    if problem is not None and args.compiled:
        manifest=json.loads((Path(args.compiled).parent/'manifest.json').read_text())
        assert manifest['rate_actions']==problem.rate_actions and not manifest['true_plant']
        if problem.rate_actions:assert manifest['prediction_version']=='sampled_force_v2'
        problem.value=ca.external('value',args.compiled);problem.vg=ca.external('vg',args.compiled)
    params=m.sample(args.n,args.seed);results=[];traces={};contexts=[]
    for i,p in enumerate(params):
        for method in args.methods.split(','):
            result,tr=episode(p,args.seed*100+i,method,cfg,problem,actor,args.steps,args.timed,args.dt,args.omega,args.stride,args.qp_reg)
            results.append({'case':i,'method':method,**result});traces[f'{i}_{method}']=tr
            if args.collect:
                cc=cfg.copy();cc[7]=p[9];zr=tr[0,2]
                for row in tr[::8]:
                    t=row[0];y=row[1:11];q=row[14:21];an=row[21:23];sensed=row[28:31]
                    contexts.append((t,y.copy(),zr,q.copy(),an.copy(),p.copy(),sensed.copy(),features(t,y,zr,q,an,p[9],cc)))
            print(i,method,'cost/forceRMS/tangentRMS/peak',np.round(result['metrics'][:4],6),'decision_ms',round(np.median(result['decision_ms']),3),flush=True)
    out={'config':vars(args),'cfg':cfg.tolist(),'parameters':params.tolist(),'results':results}
    (ROOT/f'{args.name}.json').write_text(json.dumps(out,indent=2));np.savez_compressed(ROOT/f'{args.name}_traces.npz',**traces)
    if contexts:save_contexts(contexts,ROOT/f'{args.name}_contexts.npz')


if __name__=='__main__':main()
