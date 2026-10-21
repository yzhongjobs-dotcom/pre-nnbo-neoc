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

class TruthKPreview(TruthK):
    """Oracle with spatial preview: the stiffness the footprint will meet `preview` s ahead
    at the current indentation (uses the scene p; oracle only). preview <= 0 falls back to
    the delayed TruthK."""
    def __init__(self,p,preview,wI,cP,ki,Dmin=20.,M=0.3):
        super().__init__(0.,wI,cP,ki,Dmin,M); self.p=np.asarray(p,float); self.preview=preview
    def __call__(self,o):
        import plant as P
        n=o.n
        if n>0:
            x=o.x[n-1]+self.p[5]*self.preview; d=o.log[n-1,6]
            k=P.k_eff(x,d,self.p) if d>0 else 0.
            if k>=300.: self.k=max(k,o.k_true[n-1]) if self.preview>0 else k
        k=self.k
        return (self.M,max(self.Dmin,k/self.wI),0.,self.ki,self.cP/k)

def d_star(x,p):
    """static indentation giving FD at footprint centre x (elastic part)."""
    import plant as P
    lo,hi=1e-7,0.03
    for _ in range(40):
        mid=0.5*(lo+hi)
        if P.footprint(x,mid,p)[0]>P.FD: hi=mid
        else: lo=mid
    return 0.5*(lo+hi)

class TruthKFF(TruthK):
    """TruthK gains plus a stiffness-map feedforward on the command: u_ff = d*(x + v*preview) - d*(x0),
    the change of the indentation needed for FD (oracle; knows the stiffness map, not the surface motion)."""
    def __init__(self,p,preview,lat,wI,cP,ki,Dmin=20.,M=0.3,ff_gain=1.0):
        super().__init__(lat,wI,cP,ki,Dmin,M); self.p=np.asarray(p,float); self.preview=preview; self.fg=ff_gain
        import plant as P; self.d0=d_star(P.X0,self.p)
    def __call__(self,o):
        g=list(super().__call__(o))
        x=o.x[o.n-1] if o.n>0 else self.p[5]*0
        if o.n==0:
            import plant as P; x=P.X0
        u=self.fg*(d_star(x+self.p[5]*self.preview,self.p)-self.d0)
        return tuple(g)+(u,)
