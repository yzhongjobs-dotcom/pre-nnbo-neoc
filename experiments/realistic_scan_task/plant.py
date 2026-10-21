"""Planar scanning plant for the redesigned task (realistic_scan_task).

Plant (truth) is deliberately richer than any controller-side model:
- planar 3R arm (x, z, pitch), joint-torque interface, gravity-model and
  payload errors, uncompensated joint friction, 2-step (2 ms) command delay
  and first-order torque bandwidth;
- B0 cylindrical finite-footprint contact (half-width from indentation via
  the B0 psi relation) with spatially varying modulus E(x) and friction mu(x),
  standard-linear-solid relaxation, Hunt-Crossley damping, Stribeck friction;
- respiratory surface motion;
- wrist F/T (Fx, Fz, My) with noise and bias random walk, encoder noise.
Synthetic engineering values; not a calibrated robot or tissue.
Units: SI. z is up; the probe presses downwards.
"""
import math
import numpy as np
from numba import njit

# ---- arm -------------------------------------------------------------------
LINK=np.array([0.45,0.40,0.12]); MASS=np.array([4.0,3.0,1.2]); GRAV=9.81
INER=MASS*LINK**2/12.0
# ---- probe / contact -------------------------------------------------------
R=0.03; L=0.04; LREF=0.03; FD=4.0
# gain bounds and rate limits for the four adapted gains Kx, Kz, Dx, Dz
GMIN=np.array([100.,10.,5.,5.]); GMAX=np.array([3000.,2000.,300.,300.])
GRATE=np.array([20000.,10000.,2000.,2000.])   # per second
KROT=40.0; DROT=1.5
DT=2e-4; CTRL_EVERY=5   # 5 kHz plant, 1 kHz impedance loop
# sample-index layout of the parameter vector p (see scenarios.py)
# 0 E1,1 E2,2 E3,3 xb1,4 xb2,5 wb,6 mu1,7 mu2,8 mu3,9 A_resp,10 f_resp,11 phase,
# 12 speed,13 payload,14 fc_scale,15 noise_scale,16 delay_steps(ctrl periods),17 unused
NP=18

@njit(cache=True)
def sig(u):
    if u>30: return 1.0
    if u<-30: return 0.0
    return 1.0/(1.0+math.exp(-u))

@njit(cache=True)
def field(x,p):
    """Piecewise-smooth modulus and friction at world position x (three regions)."""
    s1=sig((x-p[3])/p[5]); s2=sig((x-p[4])/p[5])
    E=p[0]*(1-s1)+p[1]*(s1-s2)+p[2]*s2
    mu=p[6]*(1-s1)+p[7]*(s1-s2)+p[8]*s2
    return E,mu

@njit(cache=True)
def psi(a):
    aa=min(max(a,1e-9),LREF*(1-1e-9)); w=math.sqrt(LREF*LREF-aa*aa)
    return (LREF*LREF-LREF*w+aa*aa*math.log((LREF+w)/aa))/(2*R)

@njit(cache=True)
def halfwidth(d):
    top=LREF*0.95
    dd=min(max(d,1e-10),psi(top)); a=min(max(math.sqrt(R*dd),1e-8),LREF*0.9)
    for _ in range(12):
        w=math.sqrt(LREF*LREF-a*a); der=a/R*math.log((LREF+w)/a)
        a=min(max(a-(psi(a)-dd)/max(der,1e-12),1e-8),top)
    return a

NQ=12
@njit(cache=True)
def footprint(xc,d,p):
    """Elastic normal force, first moment and mean friction coefficient."""
    if d<=0: return 0.,0.,0.,0.
    a=halfwidth(d); fn=0.; m1=0.; fm=0.
    for i in range(1,NQ+1):
        th=i*math.pi/(NQ+1); xi=a*math.cos(th); w=math.pi/(NQ+1)*math.sin(th)**2*a*a
        E,mu=field(xc+xi,p); pr=L/(2*R)*E*w
        fn+=pr; m1+=xi*pr; fm+=mu*pr
    return fn,m1,fm/max(fn,1e-12),a

# ---- arm kinematics/dynamics ---------------------------------------------
@njit(cache=True)
def fk(q):
    c1=q[0];c12=q[0]+q[1];c123=c12+q[2]
    x=LINK[0]*math.cos(c1)+LINK[1]*math.cos(c12)+LINK[2]*math.cos(c123)
    z=LINK[0]*math.sin(c1)+LINK[1]*math.sin(c12)+LINK[2]*math.sin(c123)
    return x,z,c123

@njit(cache=True)
def jac(q):
    a=np.array([q[0],q[0]+q[1],q[0]+q[1]+q[2]]); J=np.zeros((3,3))
    for j in range(3):
        for k in range(j,3):
            J[0,j]-=LINK[k]*math.sin(a[k]); J[1,j]+=LINK[k]*math.cos(a[k])
        J[2,j]=1.0
    return J

@njit(cache=True)
def mass_matrix(q,mscale,payload):
    a=np.array([q[0],q[0]+q[1],q[0]+q[1]+q[2]]); M=np.zeros((3,3))
    for i in range(3):
        Jv=np.zeros((2,3))
        for j in range(i+1):
            for k in range(j,i+1):
                r=LINK[k]*(0.5 if k==i else 1.0)
                Jv[0,j]-=r*math.sin(a[k]); Jv[1,j]+=r*math.cos(a[k])
        mi=MASS[i]*mscale
        for j in range(3):
            for k in range(3):
                M[j,k]+=mi*(Jv[0,j]*Jv[0,k]+Jv[1,j]*Jv[1,k])
                if j<=i and k<=i: M[j,k]+=INER[i]*mscale
    if payload>0:
        Je=jac(q)
        for j in range(3):
            for k in range(3): M[j,k]+=payload*(Je[0,j]*Je[0,k]+Je[1,j]*Je[1,k])
    return M

@njit(cache=True)
def gravity(q,mscale,payload):
    a=np.array([q[0],q[0]+q[1],q[0]+q[1]+q[2]]); g=np.zeros(3)
    for i in range(3):
        for j in range(i+1):
            dz=0.
            for k in range(j,i+1):
                r=LINK[k]*(0.5 if k==i else 1.0); dz+=r*math.cos(a[k])
            g[j]+=MASS[i]*mscale*GRAV*dz
    Je=jac(q)
    for j in range(3): g[j]+=payload*GRAV*Je[1,j]
    return g

@njit(cache=True)
def coriolis(q,qd,payload):
    h=1e-6; c=np.zeros(3); M0=mass_matrix(q,1.0,payload)
    dM=np.zeros((3,3,3))
    for k in range(3):
        qp=q.copy(); qp[k]+=h; dM[k]=(mass_matrix(qp,1.0,payload)-M0)/h
    for i in range(3):
        s=0.
        for k in range(3):
            for j in range(3):
                s+=dM[k][i,j]*qd[k]*qd[j]-0.5*dM[i][j,k]*qd[j]*qd[k]
        c[i]=s
    return c

@njit(cache=True)
def ik(x,z,th,qguess):
    q=qguess.copy()
    for _ in range(50):
        X,Z,T=fk(q); e=np.array([x-X,z-Z,th-T])
        if np.abs(e).max()<1e-12: break
        q+=np.linalg.solve(jac(q),e)
    return q

# ---- simulation state ------------------------------------------------------
# s: 0-2 q, 3-5 qd, 6-8 applied torque, 9 F1 (SLS), 10 Fe_prev, 11-13 F/T bias,
#    14 time, 15 force-PI integral, 16-18 vel filter, 19 xr0, 20 zr, 21 thr,
#    22-24 delayed torque buffer slot0, 25-27 slot1, 28 tangential reference speed
NS=40

@njit(cache=True)
def truth_wrench(s,p):
    q=s[0:3]; qd=s[3:6]; x,z,th=fk(q); J=jac(q); v=J@qd
    t=s[14]; surf=p[9]*math.sin(2*math.pi*p[10]*t+p[11]); surfd=p[9]*2*math.pi*p[10]*math.cos(2*math.pi*p[10]*t+p[11])
    d=surf-(z-R); dd=surfd-v[1]
    fe,m1,mu,a=footprint(x,d,p)
    return fe,m1,mu,a,d,dd,v,x,z,th

@njit(cache=True)
def step_chunk(s,p,gains,mode,ctrl,nsteps,rng_state,log,logi):
    """Advance the plant nsteps of DT. gains: actual (Kx,Kz,Dx,Dz) held this chunk.
    mode 0: impedance with fixed force feedforward FD.
    mode 1: impedance + explicit force PI on Fz (ctrl[0]=kp, ctrl[1]=ki).
    log rows: t, x, z, th, vx, vz, Fn_true, Ft_true, My_true, Fx_meas, Fz_meas, My_meas,
              delta, a, E_under, xr."""
    np.random.seed(rng_state)
    alpha=0.4; tau=0.3; lam=0.3; ns=p[15]
    fc=np.array([0.4,0.3,0.12])*p[14]; fv=np.array([0.6,0.5,0.2])
    for it in range(nsteps):
        q=s[0:3].copy(); qd=s[3:6].copy(); t=s[14]
        fe,m1,mu,a,d,dd,v,x,z,th=truth_wrench(s,p)
        # standard linear solid relaxation on the elastic normal force
        s[9]+=alpha*(fe-s[10])-DT*s[9]/tau; s[10]=fe
        fn=fe+s[9]
        if fn>0: fn*=max(0.,1.0+lam*dd)
        else: fn=0.
        if d<=0: fn=0.; s[9]=0.
        scale=fn/fe if fe>1e-12 else 0.
        vs=v[0]; strib=mu*(1.0+0.5*math.exp(-(vs/0.004)**2))
        ft=-strib*fn*math.tanh(vs/0.0005)-3.0*vs*(1.0 if fn>0 else 0.)
        # horizontal component of pressure on the curved probe face (B0 first moment / R):
        # an asymmetric footprint (e.g. stiffer tissue ahead) resists or assists motion
        ft+=-scale*m1/R
        my=scale*m1+R*ft
        Fext=np.array([ft,fn,my])
        # ---- 1 kHz controller with measurement noise and delay ----
        if int(round(t/DT))%CTRL_EVERY==0:
            qm=q+1e-5*ns*np.random.randn(3)
            Xm,Zm,Tm=fk(qm); Jm=jac(qm)
            vm=Jm@qd+np.array([2e-4,2e-4,2e-3])*ns*np.random.randn(3)
            for k in range(3): s[16+k]+=0.5*(vm[k]-s[16+k])
            s[11:14]+=np.array([2e-4,2e-4,1e-5])*ns*np.random.randn(3)
            Fm=Fext+s[11:14]+np.array([0.05,0.08,0.003])*ns*np.random.randn(3)
            xr=s[19]+s[28]*t
            Kx,Kz,Dx,Dz=gains[0],gains[1],gains[2],gains[3]
            ff=FD
            if mode==1:
                ef=FD-Fm[1]; s[15]+=ef*DT*CTRL_EVERY; ff=FD+ctrl[0]*ef+ctrl[1]*s[15]
            Fcmd=np.array([Kx*(xr-Xm)+Dx*(s[28]-s[16]),-ff+Kz*(s[20]-Zm)-Dz*s[17],
                           KROT*(s[21]-Tm)-DROT*s[18]])
            tau_cmd=Jm.T@Fcmd+gravity(qm,1.0,0.0)   # nominal model: no payload, nominal masses
            # two-period transport delay
            s[22:25]=s[25:28]; s[25:28]=tau_cmd
            log[logi[0],0]=t; log[logi[0],1]=Xm; log[logi[0],2]=Zm; log[logi[0],3]=Tm
            log[logi[0],4]=s[16]; log[logi[0],5]=s[17]; log[logi[0],6]=fn; log[logi[0],7]=ft
            log[logi[0],8]=my; log[logi[0],9]=Fm[0]; log[logi[0],10]=Fm[1]; log[logi[0],11]=Fm[2]
            log[logi[0],12]=d; log[logi[0],13]=a; log[logi[0],14]=field(x,p)[0]; log[logi[0],15]=xr
            logi[0]+=1
        tau_cmd_applied=s[22:25] if p[16]>=2 else s[25:28]
        s[6:9]+=DT*2*math.pi*100.0*(tau_cmd_applied-s[6:9])
        # ---- true arm dynamics (payload, mass error, joint friction) ----
        Mt=mass_matrix(q,1.01,p[13]); Jt=jac(q)
        tau_f=-(fc*np.tanh(qd/0.01)+fv*qd)
        rhs=s[6:9]+Jt.T@Fext-coriolis(q,qd,p[13])-gravity(q,1.01,p[13])+tau_f
        qdd=np.linalg.solve(Mt,rhs)
        s[3:6]=qd+DT*qdd; s[0:3]=q+DT*s[3:6]; s[14]=t+DT
    return 0

def initial_state(p,x0=0.42,gains=None):
    """Static contact equilibrium at FD with the probe pointing down."""
    s=np.zeros(NS); th=-math.pi/2
    from scipy.optimize import brentq
    surf0=p[9]*math.sin(p[11])
    dstar=brentq(lambda d: footprint(x0,d,p)[0]-FD,1e-6,0.02)
    z0=surf0-dstar+R
    q=ik(x0,z0,th,np.array([0.9,-1.6,-0.9-math.pi/2+0.7]))
    s[0:3]=q; s[10]=FD; s[19]=x0; s[20]=z0; s[21]=th; s[28]=p[12]
    s[6:9]=jac(q).T@np.array([0.,-FD,0.])+gravity(q,1.0,0.0)
    s[22:25]=s[6:9]; s[25:28]=s[6:9]
    return s

def rate_limit(g,target,period):
    lo=np.maximum(GMIN,g-GRATE*period); hi=np.minimum(GMAX,g+GRATE*period)
    return np.clip(target,lo,hi)

def episode(p,policy,T=7.0,mode=0,ctrl=(0.,0.),seed=0,x0=0.42,update=0.01):
    """policy(t, obs, gains)->target gains; called every `update` s with the
    latest logged measurements (obs = log rows so far). Returns log and metrics."""
    s=initial_state(p,x0); n_ctrl=int(round(T/DT/CTRL_EVERY))+2
    log=np.zeros((n_ctrl,16)); logi=np.zeros(1,dtype=np.int64)
    chunk=int(round(update/DT)); g=np.clip(policy(0.0,log[:0],None),GMIN,GMAX); gl=[]
    ctrl=np.asarray(ctrl,dtype=float); k=0
    while s[14]<T-1e-9:
        step_chunk(s,p,g,mode,ctrl,chunk,seed*100003+k,log,logi); k+=1
        g=rate_limit(g,np.asarray(policy(s[14],log[:logi[0]],g)),update); gl.append(g.copy())
        if not np.isfinite(s).all() or abs(s[3:6]).max()>50:
            return log[:logi[0]],None
    return log[:logi[0]],np.asarray(gl)
