"""Compare early convergence against the unchanged fixed-iteration roots.

No benchmark imports or object samples. Synthetic positive contact parameters
test algebra/numerics only, including domain-saturated cases, not task benefit.
"""
import json
import time
from pathlib import Path
import numpy as np
from finite_patch_model import FinitePatch
from line_patch_model import LinePatch


def old_inverse(model,x,f,theta):
    f=np.maximum(np.asarray(f,float),.005)
    if isinstance(model,LinePatch):
        a=np.clip(np.sqrt(f/np.exp(theta[0])),1e-8,.08)
        for _ in range(16):
            ff,_,da=model._at_a(x,a,theta);h=(ff-f)/np.maximum(da,1e-6)
            a=np.clip(a-np.clip(h,-.6*a,.6*a),1e-9,.08)
        return model.depth_of_a(a)
    d=np.clip((f/np.exp(theta[0]))**(1/model.g.exponent),1e-7,.06)
    for _ in range(16):
        ff,_,k=model.response(x,d,theta);step=(ff-f)/np.maximum(k,1.)
        d=np.clip(d-np.clip(step,-.6*d,.6*d),1e-8,.08)
    return d


def main():
    rng=np.random.default_rng(730112);rows=[]
    for cls in (FinitePatch,LinePatch):
        for spatial in (False,True):
            model=cls(spatial=spatial);cases=[]
            for i in range(50):
                x=np.linspace(-.02,.02,16);f=rng.uniform(.5,12.,16)
                scale=rng.uniform(8.,14.) if cls is FinitePatch else rng.uniform(10.,15.)
                theta=np.array([scale,rng.uniform(-2.,2.),rng.uniform(-.02,.01),rng.uniform(.003,.025)]) if spatial else np.array([scale])
                cases.append((x,f,theta))
            times=[];errors=[]
            for fn in (lambda m,x,f,t:old_inverse(m,x,f,t),lambda m,x,f,t:m.inverse(x,f,t)):
                tic=time.perf_counter();values=[fn(model,*c) for c in cases];times.append(time.perf_counter()-tic)
                errors.append(values)
            delta=max(float(np.max(np.abs(a-b))) for a,b in zip(*errors))
            assert delta<2e-13,(cls,spatial,delta)
            rows.append(dict(model=cls.__name__,spatial=spatial,cases=len(cases),
                             max_depth_difference=delta,fixed_seconds=times[0],early_seconds=times[1]))
    result=dict(protocol=dict(seed=730112,truth_or_scene_inputs=False),rows=rows)
    target=Path(__file__).resolve().parent/'results'/'contact_root_audit.json'
    target.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
