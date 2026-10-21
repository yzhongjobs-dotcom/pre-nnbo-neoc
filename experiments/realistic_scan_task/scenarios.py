"""Representative tasks for the redesigned scan simulation.

T1 boundary: scan across a soft<->stiff boundary (ratio 2-6) with friction change.
T2 inclusion: scan over a stiff inclusion (ratio 3-8, width 8-25 mm).
T3 homogeneous: scan homogeneous tissue of varying stiffness under respiration
   and robot-side model error (payload/friction), no spatial change.
All with respiratory surface motion. Effective contact moduli 20-300 kPa give
linearised contact stiffness of roughly 1-8 N/mm at 4 N (synthetic, to be
calibrated against platform/phantom data).
"""
import numpy as np
import plant as P

def make(task,rng):
    p=np.zeros(P.NP)
    E1=rng.uniform(2e4,8e4); mu1=rng.uniform(0.08,0.25)
    speed=rng.uniform(0.010,0.020); x0=0.42
    tb=rng.uniform(2.8,3.6)          # time at which probe centre reaches the first boundary
    xb=x0+speed*tb
    if task=='T1':
        ratio=rng.uniform(2,6); E2=E1*ratio if rng.random()<0.7 else E1/ratio*3
        p[:9]=[E1,E2,E2,xb,10.,rng.uniform(0.001,0.006),mu1,mu1*rng.uniform(0.6,2.0),0.]
        p[8]=p[7]
    elif task=='T2':
        ratio=rng.uniform(3,8); width=rng.uniform(0.008,0.025)
        p[:9]=[E1,E1*ratio,E1,xb,xb+width,rng.uniform(0.001,0.004),mu1,mu1*rng.uniform(0.8,1.5),mu1]
    elif task=='T3':
        E1=rng.uniform(2e4,3e5)
        p[:9]=[E1,E1,E1,10.,11.,0.003,mu1,mu1,mu1]
    else: raise ValueError(task)
    p[9]=rng.uniform(0.002,0.006); p[10]=rng.uniform(0.2,0.33); p[11]=rng.uniform(0,2*np.pi)
    p[12]=speed; p[13]=rng.uniform(0.0,0.1); p[14]=rng.uniform(0.5,1.5); p[15]=1.0; p[16]=2
    return p

def cohort(task,n,seed):
    rng=np.random.default_rng(seed); return [make(task,rng) for _ in range(n)]

TOL_F=0.5; TOL_X=0.001; FMIN=0.5; FMAX=P.FD+3.0; FTMAX=2.5
def score(log,warm=1.0):
    """Task-tolerance-normalised score (lower is better) and components."""
    if log is None: return {'score':1e3}
    m=log[:,0]>=warm; f=log[m,6]; ex=log[m,1]-log[m,15]
    ef=(f-P.FD)/TOL_F; exn=ex/TOL_X
    loss=np.mean(f<FMIN); over=np.mean(np.maximum(0,f-FMAX)**2)
    shear=np.mean(np.maximum(0,abs(log[m,7])-FTMAX)**2)
    sc=np.mean(ef**2)+np.mean(exn**2)+20*loss+over+shear
    return {'score':float(sc),'force_rmse':float(np.sqrt(np.mean((f-P.FD)**2))),'force_max_err':float(np.max(abs(f-P.FD))),
            'lateral_rmse_mm':float(np.sqrt(np.mean(ex**2))*1e3),'contact_loss_frac':float(loss),'overforce':float(over),'shear_excess':float(shear),'ft_max':float(np.max(abs(log[m,7])))}
