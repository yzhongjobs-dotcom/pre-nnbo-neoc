"""Local contact information and constrained probe design. No plant imports.

Force-domain local response, after the specified force filter:
 z = k_end Y + k_slope (t/T)Y + d_scaled T*velocity + B4 beta + noise.
The quartic drift follows from a cubic motion times a linear gain. The rate
coefficient is a local response nuisance, NOT an identified SLS viscosity.
Covariance/risk are conditional on this local model, not safety certificates.
"""
import math
import numpy as np


def residualize(A,v):
    U,s,_=np.linalg.svd(A,full_matrices=False)
    Q=U[:,s>max(float(s[0])*1e-10,1e-14)]
    return v-Q@(Q.T@v)


def features(t,y0,y1,te,ye,tw,tau=.019,fc=30.):
    """Eight 1-ms samples; endpoint encoder interpolation, common noise domain."""
    m=np.arange(8)*.001;af=1-math.exp(-2*math.pi*fc*.001)
    w=af*(1-af)**(7-np.arange(8));rho=math.exp(-.008/tau)
    frac=-np.expm1(-m/tau)/(1-rho);speed=np.exp(-m/tau)/(tau*(1-rho))
    yy=y0[:,None]+(y1-y0)[:,None]*frac-ye
    tt=(t[:,None]+m-te)/tw
    U=np.column_stack([yy@w,(tt*yy)@w,tw*((y1-y0)[:,None]*speed)@w])
    B=np.column_stack([(tt**p)@w for p in range(5)])
    return B,U


def design(measured,nw=16,tau=.019,fc=30.):
    """Only seven measured channels t,x,y,filtered F,command,Ft,My."""
    if measured.ndim!=2 or measured.shape[1]!=7:raise ValueError('seven measured channels required')
    n=len(measured);j=np.arange(max(1,n-1-nw),n-1)
    if len(j)<13:return None
    good=(measured[j,3]>.5)&(measured[j+1,3]>.5)
    j=j[good]
    if len(j)<13:return None
    t=measured[j,0];te=float(measured[-1,0]);tw=nw*.008
    B,U=features(t,measured[j,2],measured[j+1,2],te,measured[-1,2],tw,tau,fc)
    decay=math.exp(-2*math.pi*fc*.008)
    z=measured[j+1,3]-decay*measured[j,3]
    return B,U,z


def fit(measured,nw=16,tau=.019,fc=30.,noise=.04,rate=True,constant=False):
    data=design(measured,nw,tau,fc)
    if data is None:return None
    B,U,z=data
    if constant:B=B[:,:4];U=U[:,:1]
    else:U=U if rate else U[:,:2]
    ur=residualize(B,U);zr=residualize(B,z)
    norms=np.linalg.norm(ur,axis=0)
    if np.any(norms<1e-12):return None
    A=ur/norms;u,s,vt=np.linalg.svd(A,full_matrices=False)
    if s[-1]<s[0]*1e-8:return None
    coef=np.linalg.lstsq(A,zr,rcond=1e-8)[0]/norms
    inv=(vt.T/s**2)@vt;cov=noise**2*inv/norms[:,None]/norms[None,:]
    r=zr-ur@coef;dof=max(len(z)-B.shape[1]-U.shape[1],1)
    chi=float(r@r/(dof*noise**2));cov*=max(1.,chi)
    sd=float(np.sqrt(cov[0,0]));k=float(coef[0])
    risk=sd/max(k,1.)*(1+math.sqrt(dof*max(chi-1.,0.)))
    return dict(k=k,sigma=sd,relative_sd=sd/max(k,1.),risk=risk,chi2=chi,dof=dof,
                coefficients=coef,cov=cov,valid=len(z),start_gain=float(coef[0] if constant else coef[0]-coef[1]),
                residual_rms=float(np.sqrt(np.mean(r*r))),condition=float(s[0]/s[-1]))


def periodic_encoder(wave,tau=.019,delay=1):
    rho=math.exp(-.008/tau);n=len(wave)
    # Solve periodic fixed point analytically, then propagate one period.
    u=np.roll(np.asarray(wave,float),delay)
    first=(1-rho)*sum(rho**(n-1-i)*u[i] for i in range(n))/(1-rho**n)
    y=np.empty(n+1);y[0]=first
    for i in range(n):y[i+1]=rho*y[i]+(1-rho)*u[i]
    return y[:-1]


def predicted_force_rms(wave,tau=.019):
    """Nominal unit-gain, UNFILTERED physical-force ripple at 1 kHz.

    Known servo only. Unknown contact/closed-loop effects remain to be measured;
    this is neither a plant query nor a bound on actual force excursions.
    """
    y=periodic_encoder(wave,tau);n=len(y)
    frac=-np.expm1(-np.arange(8)*.001/tau)/(1-math.exp(-.008/tau))
    sub=y[:,None]+(np.roll(y,-1)-y)[:,None]*frac
    return float(np.sqrt(np.mean((sub-sub.mean())**2)))


def sine_reference_rms(frequency=20.8,tau=.019):
    # Exact steady-state first-order servo with ZOH, sampled every 1 ms.
    w=2*np.pi*frequency*.008;rho=math.exp(-.008/tau)
    G=(1-rho)*np.exp(-1j*w)/(np.exp(1j*w)-rho)
    frac=-np.expm1(-np.arange(8)*.001/tau)/(1-rho)
    A=G*(1+frac*(np.exp(1j*w)-1))
    return float(np.sqrt(np.mean(abs(A)**2)/2))


def conditional_information(wave,tau=.019,fc=30.,rate=True):
    """Worst sliding-phase information for END gain, profiling slope/rate/B4."""
    n=len(wave);y=periodic_encoder(wave,tau);info=[]
    t=np.arange(n)*.008;te=n*.008
    for phase in range(n):
        idx=(np.arange(n)+phase)%n;y0=y[idx];y1=y[(idx+1)%n]
        B,U=features(t,y0,y1,te,y[(phase+n)%n],te,tau,fc)
        other=np.column_stack([B,U[:,1:]]) if rate else np.column_stack([B,U[:,1]])
        r=residualize(other,U[:,0]);info.append(float(r@r))
    return min(info),info


def optimize_probe(count=192,seed=73621,rate=True,n=16,physical_force_budget=False):
    """C-optimal (conditional END gain) Fourier design, no task data or truth.

    Unit amplitude; peak<=1, RMS<=sine RMS, periodic slew<=original sine's
    largest 8-ms increment. Feasible scaling transfers the Fourier-constraint
    mechanism of Huang et al. TRO2025; the nuisance-profile criterion is ours.
    """
    rng=np.random.default_rng(seed);phase=np.arange(n)*2*np.pi/n
    harmonics=np.arange(1,round(n*32*.008)+1)
    H=np.sin(harmonics[:,None]*phase)
    max_slew=2*math.sin(math.pi*20.8*.008)
    best=None;reference=sine_reference_rms()
    for i in range(count):
        coef=np.eye(len(harmonics))[i] if i<len(harmonics) else rng.normal(size=len(harmonics))
        w=coef@H
        scale=max(np.max(np.abs(w)),np.sqrt(np.mean(w*w))/.7071067811865476,
                  np.max(np.abs(w-np.roll(w,1)))/max_slew)
        if physical_force_budget:scale=max(scale,predicted_force_rms(w)/reference)
        w=w/scale
        val,phases=conditional_information(w,rate=rate)
        if best is None or val>best['information']:
            best=dict(wave=w.tolist(),information=val,phase_information=phases,
                      fourier=coef.tolist(),scale=float(scale),seed=seed,candidates=count,
                      rms=float(np.sqrt(np.mean(w*w))),peak=float(np.max(np.abs(w))),
                      max_slew=float(np.max(np.abs(w-np.roll(w,1)))),
                      physical_force_budget=physical_force_budget,predicted_force_rms=predicted_force_rms(w),
                      reference_force_rms=reference)
    return best
