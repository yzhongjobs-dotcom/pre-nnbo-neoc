"""Event-aware structured projection; candidate validation, not an adopted controller.

All design columns come from logged encoder and filtered-force measurements.
No truth, patient-motion state, future measurement, or plant internals are used.
"""
import math
from functools import lru_cache
import numpy as np
from scipy.stats import chi2
import common

def _compute(log,n,nw,degree,ns,eb,afc,ndel,enc):
    j0=max(ndel+1,n-1-nw)
    m=n-1-j0
    if m<degree+8:
        return None
    te=float(log[n-1,0]);ye=float(log[n-1,2]);tw=nw*ns*.001
    B,Y,z,keep=common._build(log,j0,m,degree,ns,eb,afc,ndel,te,ye,tw,enc)
    if keep.sum()<degree+8:
        return None
    idx=np.arange(j0,j0+m)[keep]
    B,Y,z=B[keep],Y[keep],z[keep]
    mm=np.arange(ns)
    dt=.001
    w=afc*(1-afc)**(ns-1-mm)
    ebn=eb**ns
    frac=(1-eb**mm)/(1-ebn) if ebn<1 else mm/ns
    if enc:
        y=log[idx,2,None]+(log[idx+1,2]-log[idx,2])[:,None]*frac[None,:]-ye
    else:
        # Consistent with common._build on the non-encoder branch.
        y=log[idx,2,None]*eb**mm[None,:]+log[idx-ndel,11,None]*(1-eb**mm[None,:])-ye
    tm=log[idx,0,None]+mm[None,:]*dt
    return B,Y,z,tm,y,w,tw,te

def _estimate(A,z,noise,kind,te,tau,tw,degree,desired):
    # Stable QR/SVD instead of normal equations. Exact data rank check.
    try:
        U,s,Vt=np.linalg.svd(A,full_matrices=False)
        if s[-1]<=max(1e-13, s[0]*1e-9):
            return None
        th=Vt.T@((U.T@z)/s)
        r=z-A@th
        e=float(r@r)
        dof=max(len(z)-len(th),1)
        c=np.zeros(len(th)); c[0]=1.
        if kind in ("k","both","slope"):
            c[1]=1. if kind!="slope" else 0.5
        inv=Vt.T@np.diag(1/s**2)@Vt
        sd=float(noise*math.sqrt(max(c@inv@c,0.)))
        k=float(c@th)
        if not np.isfinite(k) or not np.isfinite(sd):
            return None
        # [k0, delta-k, B, motion columns] or [k0, B, motion columns].
        bi=2 if kind in ("k","both","slope") else 1
        beta=th[bi:bi+degree+1].copy()
        if kind in ("motion","both"):
            a,b=th[-2:]
            beta[0]+=a+b*(te-tau)/tw
            if len(beta)>1:beta[1]+=b
        if kind=="k":
            # A last nuisance step absorbs the force-level shift caused by
            # a changed reference surface when k changes.
            beta[0]+=th[-1]
        return dict(k=k,sd=sd,rss=e,chi2=e/(noise*noise*dof),
                    beta=beta,kind=kind, rank=len(th), dof=dof)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return None

def event_fit(log,n,nw,degree,noise,ns,eb,afc,ndel,enc=False,
              alpha=.01,mode="combined",_base=None):
    base=common.window_fit(log,n,nw,degree,noise,ns,eb,afc,ndel,enc) if _base is None else _base
    # Preserve the exact original value/uncertainty/gating on normal windows.
    def fallback():
        return dict(base, event_type="base", event_tau=None)
    if not np.isfinite(base["k"]) or base["info"]<=0 or degree!=3:
        return fallback()
    design=_compute(log,n,nw,degree,ns,eb,afc,ndel,enc)
    if design is None:return fallback()
    B,Y,z,tm,y,w,tw,te=design
    if len(z)<12:return fallback()
    A0=np.column_stack([Y,B])
    # Compare all candidates to exactly the same original fit.
    rss0=float(np.dot(z-A0@np.concatenate(([base["k"]],base["beta"])),
                      z-A0@np.concatenate(([base["k"]],base["beta"]))))
    # Include finite-width transition and a continuous local-slope explanation.
    candidates=[]
    ntests=max(1,(len(z)-7)*2*3)
    # The candidate grid uses only timestamps already measured by the current call.
    grid=np.unique(tm[:,0])[4:-3]
    for tau in grid:
        for width in (0.,.024):
            S=(tm>=tau).astype(float) if width==0 else np.clip((tm-tau)/width,0.,1.)
            HS=S@w
            HYS=(S*y)@w
            R=(np.maximum(tm-tau,0.)/tw)@w
            if mode in ("combined","movement"):
                candidates.append(("motion",np.column_stack([Y,B,HS,R]),tau,2))
            if mode in ("combined","stiffness"):
                candidates.append(("k",np.column_stack([Y,HYS,B,HS]),tau,2))
            if mode=="combined":
                candidates.append(("both",np.column_stack([Y,HYS,B,HS,R]),tau,3))
    if mode in ("combined","stiffness"):
        trend=((tm-te)/tw*y)@w
        candidates.append(("slope",np.column_stack([Y,trend,B]),te-tw,1))
    if not candidates:return fallback()
    best=None
    for kind,A,tau,nadded in candidates:
        result=_estimate(A,z,noise,kind,te,tau,tw,degree,0)
        if result is None:continue
        gain=max(0.,(rss0-result["rss"])/(noise*noise))
        # Conservative family-wise screening, candidate-position scanning included.
        thresh=float(chi2.isf(alpha/len(candidates),nadded))
        if not np.isfinite(thresh) or gain<thresh:continue
        # Complexity penalty and rank check, beyond best raw residual.
        score=result["rss"]/(noise*noise)+result["rank"]*math.log(len(z))
        if best is None or score<best[0]:
            best=(score,result,tau)
    if best is None:return fallback()
    result=best[1]
    k,sd=result["k"],result["sd"]
    # Protection is part of the scientific test, not a post-hoc acceptance relax.
    if k<common.K_MIN or k>common.K_MAX or k<=0:return fallback()
    if sd/max(k,1.)>min(.20,alpha*12+0.07):return fallback()
    if result["chi2"]>min(4.,max(float(base["chi2"]),2.7)):return fallback()
    # Standard conditional covariance ignores the selection, so guard by
    # inflating the risk using extra model complexity, not shrinking it.
    risk_sd=sd*math.sqrt(max(1.,result["chi2"]))
    return dict(base,k=k,beta=result["beta"],
                info=1/max(risk_sd**2,1e-100),chi2=result["chi2"],
                event_type=result["kind"],event_tau=float(best[2]))
