"""Reference controllers that need no estimator: fixed gains and known-environment oracles.

TruthK: continuous truth-stiffness schedule (oracle B3b). Uses the true effective contact
stiffness under the footprint, delayed by `lat` s from its onset (not from probe-centre
crossing): D = max(Dmin, k/wI), KPF = cP/k, fixed ki, M = 0.3 kg, K = 0.
"""
import numpy as np

class TruthK:
    oracle=True
    def __init__(self,lat,wI,cP,ki,Dmin=20.,M=0.3):
        self.lat=lat; self.wI=wI; self.cP=cP; self.ki=ki; self.Dmin=Dmin; self.M=M
    def reset(self,pl): self.Tc=pl['Tc']; self.k=300.
    def __call__(self,o):
        n=o.n
        if n>0:
            i=max(n-1-int(round(self.lat/self.Tc)),0)
            k=o.k_true[i]
            if k>=300.: self.k=k        # hold last contact stiffness during contact loss
        k=self.k
        return (self.M,max(self.Dmin,k/self.wI),0.,self.ki,self.cP/k)
