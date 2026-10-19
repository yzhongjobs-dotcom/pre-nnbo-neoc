"""Planar torque-interface instance of B0 equations 5, 11--15, 28--34.

Synthetic engineering values, not a calibrated robot or clinical tissue model.
State: tangential/vertical reference error, their velocities, four gains, tank.
Descriptor: E0, Ef, mu0, muf, kd0, kdf, boundary, gamma, gamma_dot,
            current centre s, fixed reference depth, reference speed.
Vertical error is positive upwards; indentation = reference_depth - e_z.
"""
import numpy as np

R=.03; L=.02; LREF=.03; FD=.6; HORIZON=.3
MASS=np.array([.8,.5]); BR=.5; CHI=.3; EMIN=.001
GMIN=np.array([50.,30.,4.,4.]); GMAX=np.array([1500.,1000.,100.,100.])
GREF=np.array([500.,250.,25.,20.]); GSCALE=np.array([500.,300.,40.,30.])
RATE=np.array([3000.,2000.,300.,300.]); ALPHA=20.
XSCALE=np.r_[.002,.001,.015,.015,GSCALE,.02]
RA=np.ones(4)*.04

def psi(a):
    aa=np.clip(a,1e-12,LREF*(1-1e-9)); w=np.sqrt(LREF**2-aa**2)
    return (LREF**2-LREF*w+aa**2*np.log((LREF+w)/aa))/(2*R)

def halfwidth(delta):
    dd=np.clip(delta,1e-10,float(psi(LREF*.95)))
    a=np.clip(np.sqrt(R*dd),1e-8,LREF*.9)
    for _ in range(9):
        w=np.sqrt(LREF**2-a*a); deriv=a/R*np.log((LREF+w)/a)
        a=np.clip(a-(psi(a)-dd)/np.maximum(deriv,1e-12),1e-8,LREF*.95)
    return a

def wrench(t,x,eta):
    x=np.asarray(x); eta=np.asarray(eta)
    e0,ef,mu0,muf,kd0,kdf,b,gamma,gd,s0,depth,speed=np.moveaxis(eta,-1,0)
    d=depth-x[...,1]-gamma-gd*t; a=halfwidth(d)
    s=s0+speed*t+x[...,0]; rr=np.clip(b-s,-a,a); w=np.sqrt(np.maximum(0,a*a-rr*rr))
    i0=np.pi*a*a/2; i2=np.pi*a**4/8
    iff=.5*(a*a*np.arccos(rr/a)-rr*w); i1=w**3/3
    i2f=a**4*np.arccos(rr/a)/8+rr*w*(a*a-2*rr*rr)/8
    c=L/(2*R); n0=c*e0*(i0-iff); nf=c*ef*iff; fn=n0+nf
    v=speed+x[...,2]; sign=np.sign(v); frac=iff/i0
    ft=-sign*(mu0*n0+muf*nf)-v*((1-frac)*kd0+frac*kdf)
    my=-c*(ef-e0)*i1-sign*c/(2*R)*(mu0*e0*(i2-i2f)+muf*ef*i2f)
    my-=v/(2*R*i0)*(kd0*(i2-i2f)+kdf*i2f)
    out=np.stack([fn,ft,my],axis=-1)
    return np.where(np.asarray(d)[...,None]>0,out,0.)

def terms(t,x,eta):
    x=np.asarray(x); f=wrench(t,x,eta); g=x[...,4:8]
    err=x[...,:2]; vel=x[...,2:4]
    torque=-g[...,:2]*err-g[...,2:]*vel; torque[...,1]-=FD
    absolute_v=vel+np.stack([np.broadcast_to(np.asarray(eta)[...,11],vel.shape[:-1]),np.zeros(vel.shape[:-1])],axis=-1)
    accel=(torque+f[...,[1,0]]-BR*absolute_v)/MASS
    pd=np.sum(g[...,2:]*vel**2,axis=-1)+BR*np.sum(absolute_v**2,axis=-1)
    power=.5*err**2*RATE[:2]
    lo=np.maximum(-1.,-ALPHA*(g-GMIN)/RATE); hi=np.minimum(1.,ALPHA*(GMAX-g)/RATE)
    budget=CHI*pd+ALPHA*(x[...,8]-EMIN)
    f0=np.concatenate([vel,accel,np.zeros_like(g),CHI*pd[...,None]],axis=-1)
    ell=(err[...,0]/.002)**2+(err[...,1]/.001)**2+.04*np.sum((vel/.015)**2,axis=-1)
    ell+=((f[...,0]-FD)/.1)**2+.08*np.sum(((g-GREF)/GSCALE)**2,axis=-1)
    return f0,ell,lo,hi,power,budget,torque

def qp(linear,lo,hi,power,budget):
    """Exact diagonal Hamiltonian QP, including the energy-gradient channel."""
    raw=-np.asarray(linear)/(2*RA); out=np.clip(raw,lo,hi)
    A=power[...,0]; B=power[...,1]; active=np.sum(power*out[...,:2],axis=-1)>budget
    aa=np.maximum(A,1e-25); bb=np.maximum(B,1e-25)
    # Minimise weighted distance on the energy-boundary line segment.
    low=np.maximum(lo[...,0],(budget-B*hi[...,1])/aa)
    high=np.minimum(hi[...,0],(budget-B*lo[...,1])/aa)
    opt=(RA[0]*raw[...,0]*B*B+RA[1]*A*(budget-B*raw[...,1]))/(RA[0]*B*B+RA[1]*A*A+1e-50)
    r0=np.clip(opt,low,high); r1=(budget-A*r0)/bb
    # Axis-aligned limiting cases have no division by a vanishing coefficient.
    r0=np.where(B<1e-18,np.minimum(out[...,0],budget/aa),r0)
    r1=np.where(B<1e-18,out[...,1],r1)
    r0=np.where(A<1e-18,out[...,0],r0)
    r1=np.where(A<1e-18,np.minimum(out[...,1],budget/bb),r1)
    out[...,0]=np.where(active,r0,out[...,0]); out[...,1]=np.where(active,r1,out[...,1])
    return out

def flow(t,x,eta,a):
    f0,ell,lo,hi,power,budget,_=terms(t,x,eta)
    # Project all requests by the same metric; inputs are normalised gain rates.
    aa=qp(-2*RA*np.asarray(a),lo,hi,power,budget)
    dx=f0.copy(); dx[...,4:8]=aa*RATE
    dx[...,8]-=np.sum(power*aa[...,:2],axis=-1)
    return dx,ell+np.sum(RA*aa*aa,axis=-1),aa

def step(t,x,eta,a,dt=.001):
    f,_,_=flow(t,x,eta,a); mid=x+.5*dt*f
    f,c,aa=flow(t+.5*dt,mid,eta,a)
    return x+dt*f,c,aa

def terminal(x,eta):
    f=wrench(HORIZON,x,eta)
    return .1*((x[...,0]/.002)**2+(x[...,1]/.001)**2+((f[...,0]-FD)/.1)**2)

def margins(t,x,eta):
    f=wrench(t,x,eta); *_,tau=terms(t,x,eta)
    delta=eta[10]-x[1]-eta[7]-eta[8]*t
    return np.r_[f[0]-.05,1.5-f[0],1.-abs(f[1]),.006-abs(x[:2]),
                 .06-abs(x[2:4]),8.-abs(tau),delta-.00005,.004-delta,
                 eta[11]+x[2]-.0005,.04-eta[11]-x[2],x[4:8]-GMIN,GMAX-x[4:8],x[8]-EMIN]

def initial(eta):
    from scipy.optimize import brentq
    q=eta.copy(); q[10]=brentq(lambda d:float(wrench(0,np.r_[0.,0.,0.,0.,GREF,.02],np.r_[q[:10],d,q[11]])[0])-FD,.00001,.004)
    return np.r_[0.,0.,0.,0.,GREF,.02],q

def current_equivalent(x,eta):
    """Uniform model matching CURRENT contact; remove only spatial preview.

    Do not keep the old region's material forever after the boundary passes.
    Uses the causal estimate, never the true field or future observations.
    """
    e=np.asarray(eta).copy();xx=np.asarray(x)
    a=halfwidth(e[...,10]-xx[...,1]-e[...,7]);r=np.clip(e[...,6]-e[...,9]-xx[...,0],-a,a)
    frac=(a*a*np.arccos(r/a)-r*np.sqrt(np.maximum(0,a*a-r*r)))/(np.pi*a*a)
    eq=e[...,0]*(1-frac)+e[...,1]*frac
    mu=(e[...,2]*e[...,0]*(1-frac)+e[...,3]*e[...,1]*frac)/np.maximum(eq,1.)
    kd=e[...,4]*(1-frac)+e[...,5]*frac
    e[...,0]=eq;e[...,1]=eq;e[...,2]=mu;e[...,3]=mu;e[...,4]=kd;e[...,5]=kd
    return e

# Keep independent vector equations available for fitting and verification.
numpy_wrench=wrench
numpy_terms=terms
def wrench(t,x,eta):
    if np.ndim(t)==0 and np.ndim(x)==1 and np.ndim(eta)==1:
        from fast_physics import wrench as scalar
        return scalar(float(t),np.asarray(x,dtype=float),np.asarray(eta,dtype=float))
    return numpy_wrench(t,x,eta)
def terms(t,x,eta):
    if np.ndim(t)==0 and np.ndim(x)==1 and np.ndim(eta)==1:
        from fast_physics import terms as scalar
        return scalar(float(t),np.asarray(x,dtype=float),np.asarray(eta,dtype=float))
    return numpy_terms(t,x,eta)
