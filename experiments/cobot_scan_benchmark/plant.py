"""Position-controlled cobot scanning plant (UR5/xArm class) for the task protocol.

Robot (kinematic, as in the Wang TASE2024 simulation model, F09:338-372): the
joint position servo is stiff, so contact forces do not deflect the arm; the
actual Cartesian position follows the command after `ndel` control periods of
transport delay through a first-order lag tau_a (UR5e TPI 0.1 per 2 ms ~ 19 ms,
F09:353-362) or, as a bracket, mostly dead time. Lateral x follows the scan path
(position controlled); the normal axis is admittance + force controlled.

Interaction controller (normal axis, penetration coordinate y = -z, every Tc):
  ef = FD - F_filt,  I = int ef dt (anti-windup +-20 N s)
  M yc'' + D yc' + K (yc - y_ref) = ef + ki I          (exact ZOH), |yc'| <= VMAX
  y_cmd = yc + KPF ef                                    (direct proportional force path)
Gains g = (M, D, K, ki, KPF) are set by a policy at its update period.

Contact: B0 cylindrical finite footprint (R 30 mm, L 40 mm, B0 psi relation)
with spatially varying modulus, standard-linear-solid relaxation, Hunt-Crossley
damping, Stribeck friction; bone-like ridges with partial aperture coverage and
optional bump; respiration and patient push/retract surface events.
F/T sampled at 1 kHz with noise and bias random walk, first-order low-pass.
Synthetic engineering values, not a calibrated robot or tissue. SI units.
"""
import math
import numpy as np
from numba import njit

R=0.03; L=0.04; LREF=0.03; FD=4.0; X0=0.42
GMIN=np.array([0.25,20.,0.,0.,0.]); GMAX=np.array([3.0,2500.,1000.,30.,1e-3])
# parameter vector p
# 0 E_soft, 1 mu, 2 A_resp, 3 f_resp, 4 phase_resp, 5 speed,
# 6 n_ridges, 7 ridge pitch, 8 ridge half-width, 9 first ridge centre (x - X0),
# 10 E_ridge, 11 ridge edge width, 12 ridge bump height, 13 aperture coverage (0..1),
# 14 noise scale, 15 n_events, 16+4k (k<3): event onset, amplitude (+ toward probe), ramp, hold
NP=28

@njit(cache=True)
def sig(u):
    if u>30: return 1.0
    if u<-30: return 0.0
    return 1.0/(1.0+math.exp(-u))

@njit(cache=True)
def ridges(x,p):
    """bone membership w (0..1), bump height and its slope at x."""
    w=0.; h=0.; hs=0.
    for k in range(int(p[6])):
        dx=x-(X0+p[9]+k*p[7]); hw=p[8]; ew=max(p[11],1e-5)
        w+=sig((dx+hw)/ew)*(1.0-sig((dx-hw)/ew))
        if abs(dx)<hw and p[12]!=0:
            ph=math.pi*dx/(2*hw)
            h+=p[12]*math.cos(ph)**2
            hs+=-p[12]*math.pi/(2*hw)*2*math.cos(ph)*math.sin(ph)
    return min(w,1.0),h,hs

@njit(cache=True)
def modulus(x,p):
    w,h,hs=ridges(x,p); c=p[13]*w
    return p[0]*(1-c)+p[10]*c

@njit(cache=True)
def surface(x,t,p):
    """surface height (+ toward the probe), its time derivative and spatial slope."""
    ph=2*math.pi*p[3]*t+p[4]
    h=p[2]*math.sin(ph); ht=p[2]*2*math.pi*p[3]*math.cos(ph)
    w,hr,hs=ridges(x,p); h+=hr; hx=hs
    for k in range(int(p[15])):
        t0=p[16+4*k]; amp=p[17+4*k]; tr=max(p[18+4*k],1e-3); hold=p[19+4*k]
        a=min(max((t-t0)/tr,0.),1.); b=min(max((t-t0-tr-hold)/tr,0.),1.)
        h+=amp*(a*a*(3-2*a)-b*b*(3-2*b))
        if 0<a<1: ht+=amp*6*a*(1-a)/tr
        if 0<b<1: ht-=amp*6*b*(1-b)/tr
    return h,ht,hx

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
    """elastic normal force, first moment, mean friction coefficient, half-width."""
    if d<=0: return 0.,0.,0.,0.
    a=halfwidth(d); fn=0.; m1=0.
    for i in range(1,NQ+1):
        th=i*math.pi/(NQ+1); xi=a*math.cos(th); w=math.pi/(NQ+1)*math.sin(th)**2*a*a
        pr=L/(2*R)*modulus(xc+xi,p)*w
        fn+=pr; m1+=xi*pr
    return fn,m1,p[1],a

# state s: 0 t, 1 x_act, 2 y_act, 3 yc, 4 vyc, 5 I, 6 F_filt, 7 bias, 8 F1 (SLS), 9 fe_prev,
#          10 d_prev, 11 y_ref, 12 tick counter, 13.. delay buffer (NB slots)
NB=8; NS=13+NB
# log columns per control tick
LOGC=12   # t, x, y, Fn_true, Ft_true, F_filt, d, E_under, k_eff, D, KPF, y_cmd

@njit(cache=True)
def k_eff(x,d,p):
    hd=2e-6
    if d<=hd: return 0.
    return (footprint(x,d+hd,p)[0]-footprint(x,d-hd,p)[0])/(2*hd)*1.4   # SLS instantaneous (1+alpha)

@njit(cache=True)
def advance(s,p,g,nticks,NC,ndel,tau_a,fc,VMAX,DT,noise,log,nlog):
    """advance nticks control periods (NC plant steps each) with gains g; returns new nlog
    (negative on divergence). noise: (n,2) standard normals per 1 kHz F/T sample."""
    alpha=0.4; tauS=0.3; lam=0.3; ns=p[14]; v=p[5]
    M=g[0]; D=g[1]; K=g[2]; ki=g[3]; KPF=g[4]
    Tc=NC*DT; NFT=max(int(round(1e-3/DT)),1); afc=1.0-math.exp(-2*math.pi*fc*1e-3)
    fn=0.; ft=0.
    for tick in range(nticks):
        for j in range(NC):
            t=s[0]; it=int(round(t/DT))
            xcmd=X0+v*max(t-ndel*Tc,0.)
            if tau_a>0: s[1]+=DT/tau_a*(xcmd-s[1])
            else: s[1]=xcmd
            h,ht,hx=surface(s[1],t,p)
            d=h+s[2]                       # y = -z_bottom
            dd=(d-s[10])/DT; s[10]=d
            fe,m1,mu,a=footprint(s[1],d,p)
            s[8]+=alpha*(fe-s[9])-DT*s[8]/tauS; s[9]=fe
            fn=fe+s[8]
            if fn>0: fn*=max(0.,1.0+lam*dd)
            else: fn=0.
            if d<=0:
                fn=0.; s[8]=0.
            strib=mu*(1.0+0.5*math.exp(-(v/0.004)**2))
            ft=-strib*fn*math.tanh(v/0.0005)-3.0*v*(1.0 if fn>0 else 0.)
            ft+=-(fn/fe if fe>1e-12 else 0.)*m1/R-fn*hx
            if it%NFT==0:
                kk=min(it//NFT,noise.shape[0]-1)
                s[7]+=2e-4*ns*noise[kk,0]
                s[6]+=afc*(fn+s[7]+0.08*ns*noise[kk,1]-s[6])
            if j==0:
                ef=FD-s[6]
                s[5]=min(max(s[5]+ef*Tc,-20.),20.)
                u=ef+ki*s[5]-K*(s[3]-s[11])
                e=math.exp(-D*Tc/M)
                vn=s[4]*e+(1-e)*u/D
                vn=min(max(vn,-VMAX),VMAX)
                s[3]+=0.5*(s[4]+vn)*Tc; s[4]=vn
                for q in range(ndel): s[13+q]=s[14+q]
                s[13+ndel]=s[3]+KPF*ef
                if nlog<log.shape[0]:
                    log[nlog,0]=t; log[nlog,1]=s[1]; log[nlog,2]=s[2]; log[nlog,3]=fn; log[nlog,4]=ft
                    log[nlog,5]=s[6]; log[nlog,6]=d; log[nlog,7]=modulus(s[1],p); log[nlog,8]=k_eff(s[1],d,p)
                    log[nlog,9]=D; log[nlog,10]=KPF; log[nlog,11]=s[13+ndel]
                    nlog+=1
            yt=s[13]
            if tau_a>0: s[2]+=DT/tau_a*(yt-s[2])
            else: s[2]=yt
            s[0]=t+DT
            if not math.isfinite(s[2]) or abs(s[2])>0.5: return -nlog
    return nlog

PLATFORMS={
 'PC125':dict(Tc=0.008,ndel=1,tau_a=0.019,fc=30.,VMAX=0.25),   # UR5/xArm 125 Hz (F07:353, F08:168-172), UR5e lag
 'PC125dead':dict(Tc=0.008,ndel=3,tau_a=0.003,fc=30.,VMAX=0.25), # latency bracket: mostly dead time
 'PC500':dict(Tc=0.002,ndel=1,tau_a=0.019,fc=30.,VMAX=0.25),    # UR5e RTDE 500 Hz (F09:500)
}

def init_state(p,dt):
    lo,hi=1e-7,0.03
    for _ in range(60):
        mid=0.5*(lo+hi)
        if footprint(X0,mid,p)[0]>FD: hi=mid
        else: lo=mid
    h0=surface(X0,0.,p)[0]; y=0.5*(lo+hi)-h0
    s=np.zeros(NS); s[1]=X0; s[2]=y; s[3]=y; s[6]=FD; s[9]=FD; s[10]=h0+y; s[11]=y; s[13:]=y
    return s

class Obs:
    """what a causal policy may read: time, measured (filtered) force, commanded and
    encoder positions, and its own gains; truth columns are for oracles only."""
    def __init__(self,log,n): self.log=log; self.n=n
    @property
    def t(self): return self.log[:self.n,0]
    @property
    def force(self): return self.log[:self.n,5]
    @property
    def y(self): return self.log[:self.n,2]
    @property
    def y_cmd(self): return self.log[:self.n,11]
    @property
    def x(self): return self.log[:self.n,1]
    # truth (oracle only)
    @property
    def k_true(self): return self.log[:self.n,8]

def episode(p,policy,T=10.0,platform='PC125',seed=0,update=None,dt=5e-4):
    """policy: either a fixed gain vector (array-like, 5) or an object with
    .reset(platform_dict) and .__call__(obs)->gains, called every `update` s
    (default: every control period). Returns log (control ticks x LOGC) or None."""
    pl=PLATFORMS[platform] if isinstance(platform,str) else platform
    NC=int(round(pl['Tc']/dt)); nt=int(round(T/pl['Tc']))
    log=np.zeros((nt+2,LOGC)); s=init_state(np.asarray(p,float),dt); p=np.asarray(p,float)
    noise=np.random.default_rng(seed).standard_normal((int(T*1000)+10,2))
    args=(NC,int(pl['ndel']),float(pl['tau_a']),float(pl['fc']),float(pl['VMAX']),dt,noise)
    if not callable(policy):
        g=np.clip(np.asarray(policy,float),GMIN,GMAX)
        n=advance(s,p,g,nt,*args,log,0)
        return None if n<0 else log[:n]
    policy.reset(pl)
    per=max(int(round((update or pl['Tc'])/pl['Tc'])),1); n=0
    while n<nt:
        g=np.clip(np.asarray(policy(Obs(log,n)),float),GMIN,GMAX)
        n=advance(s,p,g,min(per,nt-n),*args,log,n)
        if n<0: return None
    return log[:n]
