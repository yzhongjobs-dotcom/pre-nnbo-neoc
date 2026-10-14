"""Envelope derivative: reuse old minimizing actions for F_eta evaluation.

No new greedy minimization is needed to evaluate the partial derivative at a
regular minimizing branch. We audit against the full residual difference and
perturbed complete root; interpolation branch changes remain a stated limit.
"""
import numpy as np


def fast_sensitivity(self,V,P,eta,eps=1e-5):
    pert=np.array(eta,dtype=float);pert[0]+=eps
    old_rows,old_phi=self.model(eta);new_rows,new_phi=self.model(pert)
    dphi=(new_phi-old_phi)/eps
    F=np.empty_like(V);F[-1]=-dphi
    for i in range(self.N):
        oldgeom,oldstage=old_rows[i];newgeom,newstage=new_rows[i]
        _,j,beta=P[i]
        v0=self.apply_policy(V[i+1],(oldgeom,j,beta))
        v1=self.apply_policy(V[i+1],(newgeom,j,beta))
        F[i]=-(self.h*(newstage-oldstage)+v1-v0)/eps
    Cv=self.linear_solve(F,P)
    return Cv-dphi[None,:],dphi
