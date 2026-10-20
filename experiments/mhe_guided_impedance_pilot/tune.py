"""Tune the strong fixed-gain/direct-force baseline before learning, on train only."""
import json
import time
import numpy as np
from scipy.optimize import differential_evolution
from numba import njit
from core import m,plant
from run import ROOT,CFG,episode


@njit(cache=True)
def fixed_cost(p,noise,cfg0,steps=600):
    cfg=cfg0.copy();cfg[7]=p[9];y=plant.initial(p,cfg[:4],65);zr=y[1]
    history=np.zeros((steps,3));sensed=np.array([3.,.8,0.]);cost=0.
    for i in range(steps):
        t=.01*i;history[i]=plant.contact(t,y,p,65)[0]+noise[i]
        sensed=cfg[6]*history[max(0,i-2)]+(1-cfg[6])*sensed
        y,metrics=plant.advance(t,y,p,cfg[:4],0.,zr,sensed,cfg,.01,.001,65)
        cost+=metrics[0]
    return cost/(steps*.01)


def main():
    params=m.sample(12,9101);noises=np.array([np.random.default_rng(910100+i).normal(size=(600,3))*[.03,.03,.0001] for i in range(12)])
    initial=np.array([fixed_cost(p,noises[i],CFG) for i,p in enumerate(params)])
    count=[0];tic=time.perf_counter()
    def fun(v):
        cfg=CFG.copy();cfg[:4]=np.exp(v[:4]);cfg[4]=v[4];cfg[6]=v[5]
        vals=np.array([fixed_cost(p,noises[i],cfg) for i,p in enumerate(params)])
        count[0]+=1
        return float(np.mean(vals/np.maximum(.03,initial)))
    bounds=[(np.log(100),np.log(1600)),(np.log(1),np.log(800)),(np.log(5),np.log(150)),(np.log(3),np.log(150)),(0,2.5),(.1,.95)]
    result=differential_evolution(fun,bounds,popsize=5,maxiter=12,seed=9191,polish=False,workers=1)
    cfg=CFG.copy();cfg[:4]=np.exp(result.x[:4]);cfg[4]=result.x[4];cfg[6]=result.x[5]
    final=np.array([fixed_cost(p,noises[i],cfg) for i,p in enumerate(params)])
    out={"cfg":cfg.tolist(),"train_seed":9101,"n":12,"initial_cost":initial.tolist(),"final_cost":final.tolist(),"objective":float(result.fun),"evaluations":count[0],"seconds":time.perf_counter()-tic}
    (ROOT/"tuned_fixed.json").write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
    # Independent Python closed loop including MHE should match the tuned score.
    rr,_,_=episode(params[0],910100,"fixed",cfg)
    assert abs(rr["metrics"][0]-final[0])<1e-10,(rr["metrics"][0],final[0])


if __name__=="__main__":main()
