"""Compiled scalar evaluation of the same closed Hertz equations.

Avoid Python/NumPy scalar dispatch dominating nominal feasibility checks.
Vector observation fitting still uses the independent NumPy implementation.
"""
import math
import numpy as np
from numba import njit
from physics import R,L,LREF,FD,BR,CHI,EMIN,GMIN,GMAX,GREF,GSCALE,RATE,ALPHA,MASS

@njit(cache=True)
def wrench(t,x,e):
    delta=e[10]-x[1]-e[7]-e[8]*t
    if delta<=0:return np.zeros(3)
    # Same bound as physics.halfwidth, evaluated exactly once.
    top=LREF*.95;w=math.sqrt(LREF*LREF-top*top)
    dmax=(LREF*LREF-LREF*w+top*top*math.log((LREF+w)/top))/(2*R)
    dd=min(max(delta,1e-10),dmax);a=min(max(math.sqrt(R*dd),1e-8),LREF*.9)
    for _ in range(9):
        w=math.sqrt(LREF*LREF-a*a);psi=(LREF*LREF-LREF*w+a*a*math.log((LREF+w)/a))/(2*R)
        deriv=a/R*math.log((LREF+w)/a);a=min(max(a-(psi-dd)/max(deriv,1e-12),1e-8),LREF*.95)
    rr=min(max(e[6]-e[9]-e[11]*t-x[0],-a),a);w=math.sqrt(max(0,a*a-rr*rr))
    i0=math.pi*a*a/2;i2=math.pi*a**4/8
    iff=.5*(a*a*math.acos(rr/a)-rr*w);i1=w**3/3;i2f=a**4*math.acos(rr/a)/8+rr*w*(a*a-2*rr*rr)/8
    c=L/(2*R);n0=c*e[0]*(i0-iff);nf=c*e[1]*iff;fn=n0+nf;v=e[11]+x[2]
    sign=1. if v>0 else (-1. if v<0 else 0.)
    ft=-sign*(e[2]*n0+e[3]*nf)-v*((1-iff/i0)*e[4]+iff/i0*e[5])
    my=-c*(e[1]-e[0])*i1-sign*c/(2*R)*(e[2]*e[0]*(i2-i2f)+e[3]*e[1]*i2f)
    my-=v/(2*R*i0)*(e[4]*(i2-i2f)+e[5]*i2f)
    return np.array([fn,ft,my])

@njit(cache=True)
def terms(t,x,e):
    f=wrench(t,x,e);g=x[4:8];err=x[:2];vel=x[2:4]
    tau=-g[:2]*err-g[2:]*vel;tau[1]-=FD
    av=vel.copy();av[0]+=e[11];acc=(tau+np.array([f[1],f[0]])-BR*av)/MASS
    pd=np.sum(g[2:]*vel**2)+BR*np.sum(av**2)
    power=.5*err**2*RATE[:2];lo=np.maximum(-np.ones(4),-ALPHA*(g-GMIN)/RATE)
    hi=np.minimum(np.ones(4),ALPHA*(GMAX-g)/RATE);budget=CHI*pd+ALPHA*(x[8]-EMIN)
    f0=np.zeros(9);f0[:2]=vel;f0[2:4]=acc;f0[8]=CHI*pd
    ell=(err[0]/.002)**2+(err[1]/.001)**2+.04*np.sum((vel/.015)**2)
    ell+=((f[0]-FD)/.1)**2+.08*np.sum(((g-GREF)/GSCALE)**2)
    return f0,ell,lo,hi,power,budget,tau
