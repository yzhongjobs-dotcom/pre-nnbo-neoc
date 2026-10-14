"""Check projection/Bregman Hamiltonian gap and prediction-order diagnostics."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
import multinode as v

P=Path(__file__).resolve().parent
def main():
    xs=np.linspace(0,2,2001);rows=[]
    for anchor,sign in [(1.,1.),(1.1,-1.)]:
        wa,ca=v.exact_derivatives(anchor)
        for distance in [.01,.02,.04,.08]:
            theta=anchor+sign*distance;truth,_=v.exact_derivatives(theta)
            pred=wa+ca*sign*distance
            p=CubicSpline(v.NODES,np.r_[0,truth],bc_type='natural')(xs)
            ph=CubicSpline(v.NODES,np.r_[0,pred],bc_type='natural')(xs)
            u=np.clip(-p/2,-v.UB,v.UB);uh=np.clip(-ph/2,-v.UB,v.UB)
            gap=uh*uh+p*uh-u*u-p*u
            bound=(ph-p)**2/4
            assert np.min(gap)>-1e-12 and np.max(gap-bound)<1e-12
            rows.append(dict(anchor=anchor,theta=theta,distance=distance,
                node_gradient_error=float(np.max(abs(pred-truth))),
                maximum_hamiltonian_regret=float(np.max(gap)),
                maximum_bregman_bound=float(np.max(bound)),
                regret_over_distance_fourth=float(np.max(gap)/distance**4)))
    rng=np.random.default_rng(25);p=rng.normal(size=10000);ph=rng.normal(size=10000)
    u=np.clip(-p/2,-v.UB,v.UB);uh=np.clip(-ph/2,-v.UB,v.UB)
    gap=uh**2+p*uh-u**2-p*u;bound=(ph-p)**2/4
    assert np.max(gap-bound)<1e-12 and np.min(gap)>-1e-12
    result=dict(rows=rows,random_projection_cases=10000,
        maximum_gap_minus_bound=float(np.max(gap-bound)),scope='local Hamiltonian gap, not total time-varying optimality')
    (P/'regret_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
