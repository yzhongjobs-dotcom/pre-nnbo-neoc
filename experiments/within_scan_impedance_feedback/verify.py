"""Independent pressure quadrature and continuous closed-loop integration checks."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq
import model as m
HERE=Path(__file__).resolve().parent


def pressure_quad(t,y,p):
    x=p[9]*t+y[0]
    def scalar(u):
        # Written from physical definitions, without calling model.ground/wrench.
        k=p[0]*np.exp(p[1]*np.tanh((x+u-p[2])/p[3]))
        h=p[4]*(np.sin(2*np.pi*(x+u)/p[5]+p[6])-np.sin(p[6]))
        h+=p[7]*(np.sin(2*np.pi*p[8]*t+p[6])-np.sin(p[6]))
        hd=p[7]*2*np.pi*p[8]*np.cos(2*np.pi*p[8]*t+p[6])
        delta=y[1]-u*u/(2*m.R)-h
        return m.L*max(0.,k*(delta+p[12]*(y[3]-hd))) if delta>0 else 0.
    fn=quad(scalar,-m.A,m.A,epsabs=1e-8,limit=200,points=np.linspace(-m.A,m.A,21))[0]
    moment=quad(lambda u:u*scalar(u),-m.A,m.A,epsabs=1e-10,limit=200,points=np.linspace(-m.A,m.A,21))[0]
    mu=p[10]*np.exp(p[11]*np.tanh((x-p[2])/p[3]))
    ft=mu*fn*np.tanh((p[9]+y[2])/.001)+p[13]*(p[9]+y[2])
    return np.array([fn,ft,-moment+.025*ft])


def pressure_discrete(t,y,p,nq):
    u=-m.A+(np.arange(nq)+.5)*2*m.A/nq
    x=p[9]*t+y[0]
    k=p[0]*np.exp(p[1]*np.tanh((x+u-p[2])/p[3]))
    h=p[4]*(np.sin(2*np.pi*(x+u)/p[5]+p[6])-np.sin(p[6]))
    h+=p[7]*(np.sin(2*np.pi*p[8]*t+p[6])-np.sin(p[6]))
    hd=p[7]*2*np.pi*p[8]*np.cos(2*np.pi*p[8]*t+p[6])
    penetration=y[1]-u*u/(2*m.R)-h
    load=np.where(penetration>0,np.maximum(0,k*(penetration+p[12]*(y[3]-hd))),0.)
    fn=m.L*2*m.A/nq*load.sum();moment=m.L*2*m.A/nq*(u*load).sum()
    mu=p[10]*np.exp(p[11]*np.tanh((x-p[2])/p[3]))
    ft=mu*fn*np.tanh((p[9]+y[2])/.001)+p[13]*(p[9]+y[2])
    return np.array([fn,ft,-moment+.025*ft])


def gain_independent(y,target,zref):
    a=np.array([.5*y[0]**2,.5*(y[1]-zref)**2,0.,0.])
    lo=np.maximum(-m.GRATE,-m.ALPHA*(y[4:8]-m.GMIN))
    hi=np.minimum(m.GRATE,m.ALPHA*(m.GMAX-y[4:8]))
    raw=(target-y[4:8])/.035
    power=y[6]*y[2]**2+y[7]*y[3]**2+.5*(y[2]**2+y[3]**2)
    bound=m.CHI*power+m.ALPHA*(y[8]-m.E_MIN)
    u=np.clip(raw,lo,hi)
    if a@u>bound:
        f=lambda lam:a@np.clip(raw-lam*a,lo,hi)-bound
        right=1.
        while f(right)>0:right*=10.
        u=np.clip(raw-brentq(f,0.,right,xtol=1e-10)*a,lo,hi)
    edot=m.CHI*power-a@u
    if y[8]>=m.E_MAX and edot>0:edot=0.
    return u,edot


def independent_advance(t,y,p,target,ff,zref,measured,cfg,duration=.01,nq=33):
    def rhs(tt,yy):
        force=pressure_discrete(tt,yy,p,nq)
        rate,edot=gain_independent(yy,target,zref)
        err=m.FD-measured[0]
        ut=measured[1]-yy[4]*yy[0]-yy[6]*yy[2]+.5*cfg[7]
        un=m.FD+cfg[4]*err+yy[9]-yy[5]*(yy[1]-zref)-yy[7]*yy[3]+ff
        ut=np.clip(ut,-12,12);un=np.clip(un,-12,12)
        iy=cfg[5]*err
        if (yy[9]>=5 and iy>0) or (yy[9]<=-5 and iy<0):iy=0.
        out=np.r_[yy[2],yy[3],(ut-force[1]-.5*(yy[2]+p[9]))/.8,
                  (un-force[0]-.5*yy[3])/.5,rate,edot,iy]
        return out
    return solve_ivp(rhs,(t,t+duration),y,method='DOP853',rtol=2e-10,atol=1e-12).y[:,-1]


def main():
    rng=np.random.default_rng(8351);ps=m.sample(30,8352)
    cfg=np.array(json.loads((HERE/'validation.json').read_text())['selected_cfg'])
    formula=[];quaderrors=[];integrations=[];rateerrors=[];active=0
    for i,p in enumerate(ps):
        cfg[7]=p[9];y=m.initial(p,cfg[:4]);zr=y[1]
        t=rng.uniform(0,.064/p[9]);y[1]+=rng.uniform(-.0006,.0006);y[0]=rng.uniform(-.0005,.0005)
        y[2:4]=rng.uniform(-.004,.004,2)
        target=np.clip(y[4:8]*rng.uniform(.5,1.5,4)+[0,50,0,0],m.GMIN,m.GMAX)
        if i<10:y[8]=m.E_MIN+1e-6;y[0]=.01
        ff=rng.uniform(-.1,.1);measured=pressure_discrete(t,y,p,33)
        fast=m.wrench(t,y,p)[0];explicit=pressure_discrete(t,y,p,33)
        formula.append(np.max(np.abs(fast-explicit)))
        exact=pressure_quad(t,y,p)
        quaderrors.append((np.abs(fast-exact)).tolist())
        r,e=m.gain_rate(y,target,zr)[:2];ri,ei=gain_independent(y,target,zr)
        rateerrors.append(max(np.max(np.abs(r-ri)),abs(e-ei)))
        active+=bool(m.gain_rate(y,target,zr)[2])
        a=m.advance(t,y,p,target,ff,zr,measured,cfg,.01)[0]
        b=independent_advance(t,y,p,target,ff,zr,measured,cfg)
        integrations.append((np.abs(a-b)).tolist())
    out={'seed':8351,'cases':30,'formula_wrench_max_absolute_error':max(formula),
         'quadrature33_vs_adaptive_max_abs_Fn_Ft_My':np.max(quaderrors,axis=0).tolist(),
         'independent_projection_max_abs_error':max(rateerrors),'active_projection_cases':active,
         'DOP853_one_control_interval_max_state_abs_error':np.max(integrations,axis=0).tolist(),
         'scope':'Pressure quadrature, gain executor and independently coded 10ms dynamics; not whole-controller or real-tissue validation',
         'model_sha256':hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest()}
    (HERE/'numerical_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
