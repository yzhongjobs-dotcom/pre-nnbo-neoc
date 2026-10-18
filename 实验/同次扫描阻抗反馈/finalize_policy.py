"""Final train/validation iteration with the stronger, shared geometry prior."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from numba import njit
from scipy.optimize import minimize
import model as m
from run import noises
from learn_policy import pad
HERE=Path(__file__).resolve().parent


@njit(cache=True)
def objective(theta,ps,ns,counts,cfg):
    total=0.
    for i in range(len(ps)):
        total+=m.episode(ps[i],ns[i,:counts[i]],12,cfg,policy=theta,geometry_prior=True)[0][0]/len(ps)
    return total


def main():
    old=json.loads((HERE/'learned_policy.json').read_text());cfg=np.array(old['cfg'])
    ps=m.sample(18,8101);ns,counts=pad(noises(ps,8102));vp=m.sample(18,8201);vn,vc=pad(noises(vp,8202))
    start=time.perf_counter();policies={};evaluations={}
    base=np.array(old['constant_theta']);objective(base,ps[:1],ns[:1],counts[:1],cfg)
    for name,theta,ids,bounds in [
      ('continuous_constant',base,[0,3,5,6],[(np.log(1.),np.log(800.)),(np.log(2.),np.log(150.)),(0,1.5),(0,1)]),
      ('analytic_inverse',np.array([np.log((1+cfg[4])*(.3/.0015)**2/1000),-1.,0.,base[3],0.,base[5],base[6],0.]),
       [3,5,6],[(np.log(2.),np.log(150.)),(0,1.5),(0,1)])]:
        calls=[0]
        def fun(x):
            tt=theta.copy();tt[ids]=x;calls[0]+=1
            return objective(tt,ps,ns,counts,cfg)
        fit=minimize(fun,theta[ids],method='Powell',bounds=bounds,
           options={'maxiter':4,'maxfev':130,'xtol':.02,'ftol':1e-5})
        tt=theta.copy();tt[ids]=fit.x;policies[name]=tt;evaluations[name]=calls[0]
        print(name,objective(tt,ps,ns,counts,cfg),objective(tt,vp,vn,vc,cfg),tt.tolist(),flush=True)
    mean=np.array(old['actor_theta']);sd=np.array([.15,.2,.15,.15,.15,.15,.15,.15])
    rng=np.random.default_rng(8305)
    low=np.array([np.log(1.),-2.,-1.,np.log(2.),-.75,0.,0.,-1.]);high=np.array([np.log(800.),1.,1.,np.log(150.),.75,1.5,1.,1.])
    best=mean.copy();bestv=objective(best,vp,vn,vc,cfg);history=[]
    for it in range(5):
        pop=np.clip(mean+rng.normal(size=(28,8))*sd,low,high)
        pop[0]=mean;pop[1]=policies['continuous_constant'];pop[2]=policies['analytic_inverse']
        scores=np.array([objective(tt,ps,ns,counts,cfg) for tt in pop]);order=np.argsort(scores)
        elite=pop[order[:6]];mean=.3*mean+.7*elite.mean(0);sd=np.maximum(.025,.3*sd+.7*elite.std(0))
        vv=objective(pop[order[0]],vp,vn,vc,cfg)
        if vv<bestv:bestv=vv;best=pop[order[0]].copy()
        history.append({'iteration':it,'train_best_J':float(scores.min()),'validation_J':float(vv)})
        print('strong-prior CEM',it,scores.min(),vv,flush=True)
    policies['learned_feedback']=best
    output={'stage':'Final exploratory train/validation, independent test not yet generated',
      'geometry_prior':'Initial 3N/pose and calibrated initial surface frame, not material truth; shared across methods',
      'train_seed':8101,'validation_seed':8201,'train_noise_seed':8102,'validation_noise_seed':8202,
      'train_n':18,'validation_n':18,'cfg':cfg.tolist(),'policies':{k:v.tolist() for k,v in policies.items()},
      'validation_J':{k:float(objective(v,vp,vn,vc,cfg)) for k,v in policies.items()},
      'training_J':{k:float(objective(v,ps,ns,counts,cfg)) for k,v in policies.items()},
      'additional_optimization_evaluations':evaluations,'actor_additional_evaluations':140,
      'history':history,'seconds':time.perf_counter()-start,
      'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','finalize_policy.py','learn_policy.py']}}
    (HERE/'final_policies.json').write_text(json.dumps(output,indent=2)+'\n')
    print('final validation',output['validation_J'],flush=True)
if __name__=='__main__':main()
