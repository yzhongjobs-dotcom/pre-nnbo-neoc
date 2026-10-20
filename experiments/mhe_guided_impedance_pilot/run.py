"""Causal closed-loop experiments, model optimization and task-guided NN.

Run --help. Seeds/splits explicit. Optimizer is direct shooting, not HJB/NEOC.
The NN is a small feedforward actor trained with numerical task-cost gradients;
not PPO, not a reproduction of YOPO, and not a global optimality certificate.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import torch
from torch import nn
from core import m,plant,local_contact,features,analytic,encode,decode,cost_gradient,rollout
from frontend import Frontend

ROOT=Path(__file__).resolve().parent
torch.set_num_threads(1)
CFG=np.array([1000.,80.,45.,18.,.8,0.,.5,.012,1.])


def episode(p,seed,method="fixed",cfg=None,actor=None,collect=False,moment=True,steps=600,dt=.001,latency=0.):
    cfg=(CFG if cfg is None else cfg).copy();cfg[7]=p[9]
    y=plant.initial(p,cfg[:4],65);zr=y[1];front=Frontend(zr,p[9],moment)
    rng=np.random.default_rng(seed);noise=rng.normal(size=(steps,3))*[.03,.03,.0001]
    target=cfg[:4].copy();ff=0.;sensed=np.array([3.,.8,0.]);metrics=np.zeros(8)
    contexts=[];trace=[];decision_times=[];warm=None;pending=None;release=0.
    for i in range(steps):
        t=i*.01;f=plant.contact(t,y,p,65)[0];front.observe(t,y,f+noise[i])
        q,anchor=front.update(i)
        sensed=cfg[6]*front.forces[max(0,i-2)]+(1-cfg[6])*sensed
        if pending is not None and t>=release-1e-10:
            target=pending[:4].copy();ff=pending[4];pending=None
        if i%8==0:
            x=features(t,y,zr,q,anchor,p[9],cfg)
            if collect:contexts.append((t,y.copy(),zr,q.copy(),anchor.copy(),p.copy(),sensed.copy(),x))
            tic=time.perf_counter()
            base=np.tile(encode(analytic(t,y,zr,q,anchor,p[9],cfg)),2)
            if method=="fixed":proposal=np.r_[cfg[:4],0.]
            elif method=="analytic":proposal=decode(base[:5])
            elif method.startswith("nn"):
                with torch.no_grad():a=actor(torch.tensor(x,dtype=torch.float32)[None]).numpy()[0]
                if method=="nn_gains":a[[4,9]]=0.
                proposal=decode(a[:5])
            elif method in ("mpc","mpc_gains"):
                qc=np.r_[q,sensed]
                def fun(a):return cost_gradient(t,y,zr,qc,anchor,p[9],cfg,a,p,False)
                init=base if warm is None else np.r_[warm[5:],warm[5:]]
                # Compare shifted warm start with a competent analytic seed.
                if fun(base)[0]<fun(init)[0]:init=base
                bounds=[(-1.,1.)]*10
                if method=="mpc_gains":
                    init[[4,9]]=0.;bounds[4]=bounds[9]=(0.,0.)
                opt=minimize(fun,init,jac=True,bounds=bounds,method="L-BFGS-B",options={"maxiter":8,"ftol":1e-5,"maxls":8})
                warm=opt.x.copy();proposal=decode(warm[:5])
            else:raise ValueError(method)
            elapsed=time.perf_counter()-tic;decision_times.append(elapsed)
            # Same declared publication delay for a sensitivity comparison. A
            # pending result is never overwritten before publication.
            if pending is None:
                pending=proposal;release=t+latency
                if latency==0.:target=pending[:4].copy();ff=pending[4];pending=None
        old=y.copy()
        # Split the integration at a publication time, instead of silently
        # rounding it to the 10 ms servo boundary.
        if pending is not None and t<release<t+.01:
            d=release-t;y,aa=plant.advance(t,y,p,target,ff,zr,sensed,cfg,d,dt,65)
            target=pending[:4].copy();ff=pending[4];pending=None
            y,bb=plant.advance(t+d,y,p,target,ff,zr,sensed,cfg,.01-d,dt,65)
            aa[:3]+=bb[:3];aa[3]=max(aa[3],bb[3]);aa[4:7]+=bb[4:7];aa[7]=max(aa[7],bb[7])
        else:y,aa=plant.advance(t,y,p,target,ff,zr,sensed,cfg,.01,dt,65)
        metrics[:3]+=aa[:3];metrics[3]=max(metrics[3],aa[3]);metrics[4:7]+=aa[4:7];metrics[7]=max(metrics[7],aa[7])
        predicted=local_contact(t,old,q,anchor,p[9])[0]
        trace.append(np.r_[t,old,f,predicted,q,target,ff,anchor,sensed])
    duration=steps*.01;metrics[0]/=duration;metrics[1:3]=np.sqrt(metrics[1:3]/duration);metrics[4:7]/=duration
    return {"metrics":metrics.tolist(),"mhe_ms":(1000*np.array(front.timings)).tolist(),"decision_ms":(1000*np.array(decision_times)).tolist(),"fit_rms":front.fits},np.array(trace),contexts


class Actor(nn.Module):
    def __init__(self,n):
        super().__init__();self.net=nn.Sequential(nn.Linear(n,64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,10),nn.Tanh())
    def forward(self,x):return self.net(x)


def save_contexts(contexts,path):
    np.savez_compressed(path,t=np.array([c[0] for c in contexts]),y=np.array([c[1] for c in contexts]),zr=np.array([c[2] for c in contexts]),q=np.array([c[3] for c in contexts]),anchor=np.array([c[4] for c in contexts]),p=np.array([c[5] for c in contexts]),sensed=np.array([c[6] for c in contexts]),x=np.array([c[7] for c in contexts]))


def train(data,path,iterations=1000,seed=1041,true_model=True):
    torch.manual_seed(seed);rng=np.random.default_rng(seed);d=np.load(data);x=torch.tensor(d["x"],dtype=torch.float32)
    actor=Actor(x.shape[1]);opt=torch.optim.Adam(actor.parameters(),lr=.002);logs=[]
    # Imitate a feasible analytic seed briefly, then train against task losses.
    seeds=[]
    for i in range(len(x)):
        cfg=CFG.copy();cfg[7]=d["p"][i,9]
        seeds.append(np.tile(encode(analytic(d["t"][i],d["y"][i],d["zr"][i],d["q"][i],d["anchor"][i],cfg[7],cfg)),2))
    labels=torch.tensor(np.array(seeds),dtype=torch.float32)
    for _ in range(250):
        ids=rng.integers(len(x),size=64);opt.zero_grad();loss=((actor(x[ids])-labels[ids])**2).mean();loss.backward();opt.step()
    tic=time.perf_counter()
    for it in range(iterations):
        ids=rng.integers(len(x),size=32);out=actor(x[ids]);a=out.detach().numpy().astype(np.float64);g=[];vals=[]
        for j,i in enumerate(ids):
            cfg=CFG.copy();cfg[7]=d["p"][i,9];qc=np.r_[d["q"][i],d["sensed"][i]]
            value,grad=cost_gradient(d["t"][i],d["y"][i],d["zr"][i],qc,d["anchor"][i],cfg[7],cfg,a[j],d["p"][i],true_model)
            # Relative cost gradient keeps high-cost examples from dominating.
            g.append(grad/max(.1,value));vals.append(value)
        opt.zero_grad();out.backward(torch.tensor(np.array(g),dtype=torch.float32)/len(ids));nn.utils.clip_grad_norm_(actor.parameters(),5.);opt.step()
        if it%100==0:
            row={"iteration":it,"batch_cost":float(np.mean(vals)),"seconds":time.perf_counter()-tic};logs.append(row);print(row,flush=True)
    torch.save({"state":actor.state_dict(),"input":x.shape[1],"seed":seed,"iterations":iterations,"true_training_plant":true_model},path)
    Path(str(path)+".json").write_text(json.dumps(logs,indent=2))


def main():
    ap=argparse.ArgumentParser();ap.add_argument("mode",choices=["pilot","collect","train","eval"]);ap.add_argument("--n",type=int,default=6);ap.add_argument("--seed",type=int,default=9101);ap.add_argument("--steps",type=int,default=600);ap.add_argument("--iterations",type=int,default=1000);ap.add_argument("--name",default="pilot");ap.add_argument("--methods",default="fixed,analytic,mpc");ap.add_argument("--actor",default="actor.pt");ap.add_argument("--data",default="contexts.npz");ap.add_argument("--cfg");ap.add_argument("--latency",type=float,default=0.);args=ap.parse_args()
    cfg=np.array(json.loads(Path(args.cfg).read_text())["cfg"]) if args.cfg else CFG.copy()
    # Training config is explicitly kept consistent; do not silently train with
    # the default configuration then report a different low-level controller.
    CFG[:]=cfg
    if args.mode=="train":return train(ROOT/args.data,ROOT/args.actor,args.iterations)
    actor=None
    if "nn" in args.methods:
        ck=torch.load(ROOT/args.actor,weights_only=True);actor=Actor(ck["input"]);actor.load_state_dict(ck["state"]);actor.eval()
    results=[];contexts=[];traces={};params=m.sample(args.n,args.seed)
    for i,p in enumerate(params):
        for method in args.methods.split(","):
            result,trace,ctx=episode(p,args.seed*100+i,method,cfg,actor,args.mode=="collect",steps=args.steps,latency=args.latency)
            results.append({"case":i,"method":method,**result});traces[f"{i}_{method}"]=trace;contexts+=ctx
            print(i,method,np.round(result["metrics"][:4],6),"MHEms",round(np.median(result["mhe_ms"]),2),"decisionms",round(np.median(result["decision_ms"]),2),flush=True)
    out={"config":vars(args),"cfg":cfg.tolist(),"parameters":params.tolist(),"results":results}
    (ROOT/f"{args.name}.json").write_text(json.dumps(out,indent=2));np.savez_compressed(ROOT/f"{args.name}_traces.npz",**traces)
    if contexts:save_contexts(contexts,ROOT/args.data)


if __name__=="__main__":main()
