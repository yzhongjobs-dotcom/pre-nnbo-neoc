"""Candidate 2 independent training: task gradients through actual gain rates.

v2 isolates gains instead of allowing feedforward to mask gain quality. A
constant feasible controller initializes training; no material truth as input.
"""
import argparse,json,time,sys
from pathlib import Path
import numpy as np
import torch
import casadi as ca
from torch import nn
from numba import njit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'mhe_guided_impedance_pilot'))
from core import rollout,encode
from rate_core import rate_gradient

ROOT=Path(__file__).resolve().parent;torch.set_num_threads(1)


class GuidedActor(nn.Module):
    def __init__(self,n,gains_only=True,task_log=True,no_front=False):
        super().__init__();self.net=nn.Sequential(nn.Linear(n,64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,10),nn.Tanh())
        self.register_buffer('mask',torch.tensor([1,1,1,1,0 if gains_only else 1]*2,dtype=torch.float32))
        self.task_log=task_log
        self.no_front=no_front
    def forward(self,x):
        if self.no_front:
            x=x.clone();x[...,10:17]=0.;x[...,20]=0.
        if self.task_log:
            x=x.clone();x[...,-1]=torch.log(x[...,-1].clamp(min=1e-6))/np.log(25.)
        return self.net(x)*self.mask


@njit(cache=True)
def gradient(t,y,zr,q,anchor,speed,cfg,a,p,gains_only=True):
    value=rollout(t,y,zr,q,anchor,speed,cfg,a,p,True);g=np.zeros(10);eps=.003
    for j in range(10):
        if gains_only and j%5==4:continue
        aa=a.copy();bb=a.copy();aa[j]+=eps;bb[j]-=eps
        g[j]=(rollout(t,y,zr,q,anchor,speed,cfg,aa,p,True)-rollout(t,y,zr,q,anchor,speed,cfg,bb,p,True))/(2*eps)
    return value,g


def train(data,output,iterations=700,gains_only=True,seed=9141,no_front=False,rate_actions=False,resume=None,compiled=None):
    d=np.load(data);cfg=np.array(json.loads((ROOT.parent/'mhe_guided_impedance_pilot'/'tuned_fixed.json').read_text())['cfg'])
    torch.manual_seed(seed);rng=np.random.default_rng(seed);x=torch.tensor(d['x'],dtype=torch.float32);actor=GuidedActor(x.shape[1],gains_only,no_front=no_front)
    if resume:
        previous=torch.load(resume,weights_only=True);actor.load_state_dict(previous['state'])
    if no_front:x[:,17:20]=torch.tensor((d['sensed']-[3.,.8,0])/[1.,1.,.03],dtype=torch.float32)
    opt=torch.optim.Adam(actor.parameters(),lr=.0015)
    for _ in range(0 if resume else 200):
        ids=rng.integers(len(x),size=64);xx=x[ids].clone();weights=rng.choice([1.,12.5,25.],size=64);xx[:,-1]=torch.tensor(weights,dtype=torch.float32)
        labels=[]
        for j,w in enumerate(weights):
            target=cfg[:4].copy();target[1]=min(800.,target[1]*w)
            if rate_actions:
                rate=np.clip((target-d['y'][ids[j],4:8])/.24/np.array([1200.,1200.,120.,120.]),-.8,.8);labels.append(np.tile(np.r_[rate,0.],2))
            else:labels.append(np.tile(encode(np.r_[target,0.]),2))
        opt.zero_grad();loss=((actor(xx)-torch.tensor(np.array(labels),dtype=torch.float32))**2).mean();loss.backward();opt.step()
    log=[];tic=time.perf_counter()
    if compiled:
        manifest=json.loads((Path(compiled).parent/'manifest.json').read_text());assert manifest['true_plant'] and manifest['rate_actions'] and manifest['prediction_version']=='sampled_force_v2'
    cost_ad=ca.external('vg',str(compiled)) if compiled else None
    for it in range(iterations):
        ids=rng.integers(len(x),size=32);weights=rng.choice([1.,12.5,25.],size=32);xx=x[ids].clone();xx[:,-1]=torch.tensor(weights,dtype=torch.float32)
        out=actor(xx);actions=out.detach().numpy().astype(np.float64);grads=[];values=[]
        for j,i in enumerate(ids):
            cc=cfg.copy();cc[7]=d['p'][i,9];cc[8]=weights[j];qc=np.r_[d['q'][i],d['sensed'][i]]
            if cost_ad:
                context=np.r_[d['t'][i],d['y'][i],d['zr'][i],d['p'][i],d['anchor'][i],d['sensed'][i],cc]
                vv,gg=cost_ad(actions[j],context);v=float(vv);g=np.asarray(gg).ravel()
            elif rate_actions:v,g=rate_gradient(d['t'][i],d['y'][i],d['zr'][i],qc,d['anchor'][i],cc[7],cc,actions[j],d['p'][i])
            else:v,g=gradient(d['t'][i],d['y'][i],d['zr'][i],qc,d['anchor'][i],cc[7],cc,actions[j],d['p'][i],gains_only)
            grads.append(g/max(.1,v));values.append(v)
        opt.zero_grad();out.backward(torch.tensor(np.array(grads),dtype=torch.float32)/len(ids));nn.utils.clip_grad_norm_(actor.parameters(),5.);opt.step()
        if it%100==0:
            row={'iteration':it,'batch_cost':float(np.mean(values)),'seconds':time.perf_counter()-tic};log.append(row);print(row,flush=True)
    torch.save({'kind':'guided_gains' if gains_only else 'guided_joint','input':x.shape[1],'state':actor.state_dict(),'gains_only':gains_only,'task_log':True,'no_front':no_front,'rate_actions':rate_actions,'task_weights':[1.,12.5,25.],'iterations':iterations,'seed':seed,'cfg':cfg.tolist(),'data':str(data),'resume':str(resume),'gradient':'compiled AD' if compiled else 'finite difference'},output)
    Path(str(output)+'.json').write_text(json.dumps({'cfg':cfg.tolist(),'training_contexts':len(x),'seconds':time.perf_counter()-tic,'log':log},indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',default='training72_contexts.npz');ap.add_argument('--output',default='candidate2_v2.pt');ap.add_argument('--iterations',type=int,default=700);ap.add_argument('--joint',action='store_true');ap.add_argument('--no-front',action='store_true');ap.add_argument('--rates',action='store_true');ap.add_argument('--resume');ap.add_argument('--compiled');args=ap.parse_args()
    if (ROOT/args.output).exists():raise FileExistsError('Do not overwrite a checkpoint.')
    if args.compiled and not args.rates:raise ValueError('The compiled training kernel is for pure rate actions.')
    train(ROOT/args.data,ROOT/args.output,args.iterations,not args.joint,no_front=args.no_front,rate_actions=args.rates,resume=ROOT/args.resume if args.resume else None,compiled=args.compiled)
