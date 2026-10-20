"""Candidate 2: causal estimate-conditioned task-gradient gain-rate learning.

The simulation field is used ONLY to calculate the offline loss. The actor
receives measured motion/execution states and estimated, timestamped fields.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from torch import nn
import physics as p

torch.set_num_threads(1)
HERE=Path(__file__).resolve().parent
ESCALE=np.array([1e5,1e5,.2,.2,3.,3.,.01,.001,.001,.01,.001,.01])
CENTRE=np.r_[np.zeros(4),p.GREF,.02]

def features(x,eta,no_preview=False):
    ee=p.current_equivalent(x,eta) if no_preview else eta.copy()
    ee[...,6]-=ee[...,9];ee[...,9]=0.
    return np.concatenate([(x-CENTRE)/p.XSCALE,ee/ESCALE],axis=-1)

class Actor(nn.Module):
    def __init__(self,no_preview=False):
        super().__init__();self.no_preview=no_preview;self.current_equivalent=no_preview
        self.net=nn.Sequential(nn.Linear(21,64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,4),nn.Tanh())
        nn.init.zeros_(self.net[-2].weight);nn.init.zeros_(self.net[-2].bias)
    def forward(self,f):
        if self.no_preview:
            f=f.clone();f[:,10]=f[:,9];f[:,12]=f[:,11];f[:,14]=f[:,13];f[:,15]=0.
        return self.net(f)

def contact(t,x,eta):
    delta=eta[:,10]-x[:,1]-eta[:,7]-eta[:,8]*t;dd=torch.clamp(delta,min=1e-10,max=float(p.psi(p.LREF*.95)))
    a=torch.sqrt(p.R*dd)
    for _ in range(7):
        a=torch.clamp(a,1e-8,p.LREF*.95);w=torch.sqrt(p.LREF**2-a*a)
        psi=(p.LREF**2-p.LREF*w+a*a*torch.log((p.LREF+w)/a))/(2*p.R)
        a=a-(psi-dd)/(a/p.R*torch.log((p.LREF+w)/a))
    a=torch.clamp(a,1e-8,p.LREF*.95)
    rr=eta[:,6]-eta[:,9]-eta[:,11]*t-x[:,0];ratio=torch.clamp(rr/a,-1+1e-7,1-1e-7)
    rc=ratio*a;ww=torch.sqrt(torch.clamp(a*a-rc*rc,min=1e-25));i0=np.pi*a*a/2
    iff=.5*(a*a*torch.acos(ratio)-rc*ww);iff=torch.where(rr>=a,0.,torch.where(rr<=-a,i0,iff))
    n0=p.L/(2*p.R)*eta[:,0]*(i0-iff);nf=p.L/(2*p.R)*eta[:,1]*iff;fn=n0+nf
    v=eta[:,11]+x[:,2];ft=-torch.sign(v)*(eta[:,2]*n0+eta[:,3]*nf)-v*((1-iff/i0)*eta[:,4]+iff/i0*eta[:,5])
    return torch.stack([torch.where(delta>0,fn,0.),torch.where(delta>0,ft,0.)],dim=-1)

def flow(t,x,eta,request):
    tensor=lambda a:torch.as_tensor(a,dtype=x.dtype,device=x.device)
    g=x[:,4:8];err=x[:,:2];vel=x[:,2:4];rate=tensor(p.RATE)
    lo=torch.maximum(tensor(-np.ones(4)),-p.ALPHA*(g-tensor(p.GMIN))/rate)
    hi=torch.minimum(tensor(np.ones(4)),p.ALPHA*(tensor(p.GMAX)-g)/rate)
    aa=torch.maximum(lo,torch.minimum(hi,request));power=.5*err**2*rate[:2]
    av=vel+torch.stack([eta[:,11],torch.zeros_like(t)],dim=-1)
    pd=torch.sum(g[:,2:]*vel**2,dim=-1)+p.BR*torch.sum(av**2,dim=-1)
    budget=p.CHI*pd+p.ALPHA*(x[:,8]-p.EMIN);A=power[:,0];B=power[:,1]
    active=torch.sum(power*aa[:,:2],dim=-1)>budget
    # Inactive branches must remain numerically well conditioned too: torch
    # evaluates both branches, and 0*inf would otherwise poison gradients.
    A=torch.where(active,A,torch.ones_like(A));B=torch.where(active,B,torch.ones_like(B))
    asafe=torch.clamp(A,min=1e-12);bsafe=torch.clamp(B,min=1e-12)
    low=torch.maximum(lo[:,0],(budget-B*hi[:,1])/asafe);high=torch.minimum(hi[:,0],(budget-B*lo[:,1])/asafe)
    opt=(request[:,0]*B*B+A*(budget-B*request[:,1]))/(A*A+B*B+1e-30)
    r0=torch.maximum(low,torch.minimum(high,opt));r1=(budget-A*r0)/bsafe
    r0=torch.where(B<1e-18,torch.minimum(aa[:,0],budget/asafe),r0);r1=torch.where(B<1e-18,aa[:,1],r1)
    r0=torch.where(A<1e-18,aa[:,0],r0);r1=torch.where(A<1e-18,torch.minimum(aa[:,1],budget/bsafe),r1)
    executed=torch.stack([torch.where(active,r0,aa[:,0]),torch.where(active,r1,aa[:,1]),aa[:,2],aa[:,3]],dim=-1)
    f=contact(t,x,eta);torque=-g[:,:2]*err-g[:,2:]*vel-tensor([0.,p.FD]);acc=(torque+f[:,[1,0]]-p.BR*av)/tensor(p.MASS)
    edot=p.CHI*pd-torch.sum(power*executed[:,:2],dim=-1)
    dx=torch.cat([vel,acc,executed*rate,edot[:,None]],dim=-1)
    loss=(err[:,0]/.002)**2+(err[:,1]/.001)**2+.04*torch.sum((vel/.015)**2,dim=-1)
    loss+=((f[:,0]-p.FD)/.1)**2+.08*torch.sum(((g-tensor(p.GREF))/tensor(p.GSCALE))**2,dim=-1)
    loss+=torch.sum(tensor(p.RA)*executed**2,dim=-1)
    return dx,loss

def task_loss(x,eta,action,dt=.01):
    loss=torch.zeros(len(x),dtype=x.dtype);t=torch.zeros(len(x),dtype=x.dtype)
    for _ in range(round(p.HORIZON/dt)):
        dx,_=flow(t,x,eta,action);mid=x+.5*dt*dx;dx,ell=flow(t+.5*dt,mid,eta,action)
        x=x+dt*dx;t=t+dt;loss=loss+dt*ell
    fn=contact(t,x,eta)[:,0]
    return loss+.1*((x[:,0]/.002)**2+(x[:,1]/.001)**2+((fn-p.FD)/.1)**2)

def train(args):
    torch.manual_seed(args.seed);rng=np.random.default_rng(args.seed)
    data=np.load(HERE/args.data);ft=torch.tensor(features(data['x'],data['estimate'],args.no_preview),dtype=torch.float32)
    x=torch.tensor(data['x'],dtype=torch.float32);truth=torch.tensor(data['truth'],dtype=torch.float32)
    model=Actor(args.no_preview);optim=torch.optim.Adam(model.parameters(),lr=.0007);history=[];tic=time.perf_counter()
    for j in range(args.iterations):
        ix=rng.integers(0,len(x),48);a=model(ft[ix]);loss=task_loss(x[ix],truth[ix],a).mean()
        optim.zero_grad();loss.backward()
        if not torch.isfinite(loss) or any(not torch.isfinite(v.grad).all() for v in model.parameters() if v.grad is not None):
            raise FloatingPointError(f'nonfinite training loss/gradient at step {j}')
        torch.nn.utils.clip_grad_norm_(model.parameters(),10.);optim.step()
        if j%100==0:
            row={'step':j,'loss':float(loss.detach()),'seconds':time.perf_counter()-tic};history.append(row);print(json.dumps(row),flush=True)
    path=HERE/args.output
    if path.exists():raise FileExistsError(path)
    torch.save({'state_dict':model.state_dict(),'no_preview':args.no_preview,'seed':args.seed,'iterations':args.iterations,
                'data':args.data,'feature_version':'current_contact_equivalent_v2','seconds':time.perf_counter()-tic,'history':history},path)

def load(name):
    ckpt=torch.load(HERE/name,map_location='cpu',weights_only=True);model=Actor(ckpt['no_preview']);model.load_state_dict(ckpt['state_dict'])
    model.current_equivalent=model.no_preview and ckpt.get('feature_version')=='current_contact_equivalent_v2'
    model.eval();return model

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',default='training.npz');ap.add_argument('--output',required=True)
    ap.add_argument('--seed',type=int,default=16111);ap.add_argument('--iterations',type=int,default=600);ap.add_argument('--no-preview',action='store_true');train(ap.parse_args())
