"""Task families of the protocol (README) and the shared metrics.

T1 ladder   : soft base with 1-2 bone-like ridges (ratio r), breathing; timeliness instrument.
T2 push     : patient push/retract on homogeneous soft or stiff tissue; non-inferiority /
              false-adaptation test.
T3 spine    : spine sagittal scan (3-5 spinous-process ridges, 10-30x), breathing, 1-2 patient
              motion events at random times; headline scanning task.
N1 homog    : homogeneous tissue 10 kPa-2 MPa with breathing; negative control.
Ranges are synthetic engineering values; ridge geometry and patient-motion speed are not
sourced from anatomy or patient data and are treated as dose variables.
"""
import numpy as np
import plant as P

def _base(rng,E,A,speed,T):
    p=np.zeros(P.NP)
    p[0]=E; p[1]=rng.uniform(0.05,0.4); p[2]=A; p[3]=rng.uniform(0.2,0.33); p[4]=rng.uniform(0,2*np.pi)
    p[5]=speed; p[13]=1.0; p[14]=1.0
    return p

def _ridges(p,rng,n,ratio,pitch,hw,first,bump,cover):
    p[6]=n; p[7]=pitch; p[8]=hw; p[9]=first; p[10]=p[0]*ratio; p[11]=rng.uniform(0.0005,0.002)
    p[12]=bump; p[13]=cover

def _events(p,rng,k,T,amp,ramp,on_range=None):
    lo,hi=on_range if on_range else (1.5,T-2.0)
    ts=np.sort(rng.uniform(lo,hi,k))
    for i in range(k):
        if i>0: ts[i]=max(ts[i],ts[i-1]+1.0+2*ramp)
        a=amp*(1 if rng.random()<0.5 else -1)
        p[16+4*i:20+4*i]=[ts[i],a,ramp,rng.uniform(0.5,1.5)]
    p[15]=k

def make(task,rng,dose=None):
    """dose: dict overriding the swept variables (A, speed, ratio, ramp, amp)."""
    d=dict(dose or {})
    if task=='T1':
        T=8.0; E=rng.uniform(1e4,4e4); speed=d.get('speed',rng.choice([0.010,0.020,0.030]))
        p=_base(rng,E,d.get('A',rng.uniform(0.002,0.005)),speed,T)
        t_edge=rng.uniform(1.5,3.0); hw=rng.uniform(0.004,0.0125)
        n=int(rng.integers(1,3))
        _ridges(p,rng,n,d.get('ratio',rng.uniform(10,30)),hw*2+rng.uniform(0.01,0.02),hw,speed*t_edge+hw,0.0,1.0)
    elif task=='T2':
        T=8.0; stiff=d.get('stiff',rng.random()<0.5)
        E=rng.uniform(5e5,1.5e6) if stiff else rng.uniform(1e4,4e4)
        p=_base(rng,E,d.get('A',rng.uniform(0.002,0.005)),0.015,T)
        _events(p,rng,int(rng.integers(1,4)),T,d.get('amp',rng.uniform(0.003,0.008)),d.get('ramp',rng.uniform(0.15,0.5)))
    elif task=='T3':
        T=10.0; E=rng.uniform(1e4,4e4); speed=d.get('speed',rng.uniform(0.010,0.020))
        p=_base(rng,E,d.get('A',rng.uniform(0.002,0.005)),speed,T)
        hw=rng.uniform(0.004,0.010); pitch=rng.uniform(0.025,0.035)
        first=rng.uniform(0.005,0.02)+hw
        n=int(np.ceil((speed*T+0.03-first)/pitch))+1
        _ridges(p,rng,min(n,5),d.get('ratio',rng.uniform(10,30)),pitch,hw,first,
                d.get('bump',0.0),d.get('cover',1.0))
        _events(p,rng,int(rng.integers(1,3)),T,d.get('amp',rng.uniform(0.003,0.008)),d.get('ramp',rng.uniform(0.15,0.5)))
    elif task=='N1':
        T=6.0; E=d.get('E',np.exp(rng.uniform(np.log(1e4),np.log(2e6))))
        p=_base(rng,E,d.get('A',rng.uniform(0.002,0.005)),0.015,T)
    else: raise ValueError(task)
    return p,T

def cohort(task,n,seed,dose=None):
    rng=np.random.default_rng(seed); return [make(task,rng,dose) for _ in range(n)]

# splits by seed (README): train 1000+, dev 2000+, test 3000+
SPLIT={'train':1000,'dev':2000,'test':3000}

TOL_F=0.5; FMIN=0.5; FPEAK=P.FD+3.0
def metrics(log,warm=1.0,Tc=0.008):
    """two-sided physical metrics (true contact force) and the tuning composite J."""
    if log is None or log.shape[0]<10:
        return dict(J=1e3,rmse=np.nan,p95=np.nan,maxerr=np.nan,inband=0.,peak=np.nan,t_over=np.nan,
                    loss=1.0,loss_max=np.nan,osc=np.nan,shear=np.nan,diverged=1)
    m=log[:,0]>=warm; f=log[m,3]; e=f-P.FD
    loss=f<FMIN
    # longest contiguous loss
    lm=0; c=0
    for b in loss:
        c=c+1 if b else 0; lm=max(lm,c)
    # oscillation: RMS of force above 10 Hz (difference from a 10 Hz moving average)
    w=max(int(round(0.1/Tc)),1); k=np.ones(w)/w
    osc=float(np.sqrt(np.mean((f-np.convolve(f,k,mode='same'))[w:-w]**2))) if len(f)>3*w else np.nan
    J=float(np.mean((e/TOL_F)**2)+20*np.mean(loss)+np.mean(np.maximum(0,f-FPEAK)**2))
    return dict(J=J,rmse=float(np.sqrt(np.mean(e**2))),p95=float(np.percentile(abs(e),95)),maxerr=float(np.max(abs(e))),
                inband=float(np.mean(abs(e)<=TOL_F)),peak=float(f.max()),t_over=float(np.sum(f>FPEAK)*Tc),
                loss=float(np.mean(loss)),loss_max=float(lm*Tc),osc=osc,shear=float(np.max(abs(log[m,4]))),diverged=0)
