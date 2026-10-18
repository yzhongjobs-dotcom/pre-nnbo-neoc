"""Direct task-loss structured policy search, versus continuous-tuned simple policy.

Cross-entropy search is a mature optimizer, not a novelty claim. No gain-label MSE,
no privileged state in actor, no neural network necessity claim.
"""
from pathlib import Path
import json,time,hashlib
import numpy as np
from numba import njit
from scipy.optimize import minimize
import model as m
from run import noises
HERE=Path(__file__).resolve().parent


def pad(ns):
    out=np.zeros((len(ns),max(map(len,ns)),3));counts=[]
    for i,n in enumerate(ns):out[i,:len(n)]=n;counts.append(len(n))
    return out,np.array(counts)


@njit(cache=True)
def loss(theta,ps,ns,counts,cfg):
    total=0.
    for i in range(len(ps)):
        metric=m.episode(ps[i],ns[i,:counts[i]],12,cfg,policy=theta)[0]
        total+=metric[0]/len(ps)
    return total


def main():
    config=json.loads((HERE/'joint_validation.json').read_text())
    cfg=np.array(config['cfgs']['fixed_PI'])
    ps=m.sample(18,8101);ns,counts=pad(noises(ps,8102))
    vp=m.sample(18,8201);vn,vcounts=pad(noises(vp,8202))
    base=np.array([np.log(max(1.,cfg[1])),0.,0.,np.log(cfg[3]),0.,.75,1.,0.])
    loss(base,ps[:1],ns[:1],counts[:1],cfg)
    start=time.perf_counter();calls=0
    def constant_obj(x):
        nonlocal calls
        theta=base.copy();theta[[0,3,5,6]]=x
        value=loss(theta,ps,ns,counts,cfg);calls+=1
        return value
    fit=minimize(constant_obj,base[[0,3,5,6]],method='Powell',
       bounds=[(np.log(1.),np.log(800.)),(np.log(2.),np.log(150.)),(0.,1.5),(0.,1.)],
       options={'maxiter':5,'maxfev':200,'xtol':.02,'ftol':1e-5})
    constant=base.copy();constant[[0,3,5,6]]=fit.x
    consttrain=loss(constant,ps,ns,counts,cfg);constval=loss(constant,vp,vn,vcounts,cfg)
    print('continuous simple',consttrain,constval,constant.tolist(),flush=True)
    rng=np.random.default_rng(8301)
    low=np.array([np.log(1.),-2.,-1.,np.log(2.),-.75,0.,0.,-1.])
    high=np.array([np.log(800.),1.,1.,np.log(150.),.75,1.5,1.,1.])
    mean=constant.copy();sd=np.array([.35,.6,.25,.25,.3,.25,.25,.25])
    best=constant.copy();bestv=constval;history=[]
    for it in range(10):
        population=np.clip(mean+rng.normal(size=(28,8))*sd,low,high)
        population[0]=constant;population[1]=mean
        scores=np.array([loss(t,ps,ns,counts,cfg) for t in population])
        order=np.argsort(scores);elite=population[order[:6]]
        mean=.3*mean+.7*elite.mean(0)
        sd=np.maximum(.035,.3*sd+.7*elite.std(0))
        vloss=loss(population[order[0]],vp,vn,vcounts,cfg)
        if vloss<bestv:bestv=vloss;best=population[order[0]].copy()
        history.append({'iteration':it,'train_best':float(scores.min()),'validation':float(vloss)})
        print('CEM',it,scores.min(),vloss,'bestval',bestv,flush=True)
    out={'method':'Direct full-scan task-loss cross-entropy search of a bounded structured feedback policy',
      'not_claimed':'PPO, full neural actor or optimality certificate','train_seed':8101,'validation_seed':8201,
      'train_scans':18,'validation_scans':18,'constant_training_evaluations':calls,'actor_training_evaluations':280,
      'constant_theta':constant.tolist(),'actor_theta':best.tolist(),'cfg':cfg.tolist(),
      'constant_train_J':float(consttrain),'constant_validation_J':float(constval),'actor_validation_J':float(bestv),
      'validation_used':'one train-best candidate each iteration; initial continuous simple policy remains eligible',
      'history':history,'training_seconds':time.perf_counter()-start,
      'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','learn_policy.py']}}
    (HERE/'learned_policy.json').write_text(json.dumps(out,indent=2)+'\n')
    print('saved policies',flush=True)
if __name__=='__main__':main()
