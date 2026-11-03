"""Shared moving-surface state for elastic AND memory contact contributions.

Same SLS as memory_contact; replaces free force-polynomial coefficients with
one local constant-acceleration (or constant-jerk) surface trajectory. Its
contribution is evaluated through BOTH contact branches, with unknown initial
memory. Variable projection eliminates motion coefficients, not physics.
"""
import math
import numpy as np
from scipy.optimize import least_squares
from local_contact import residualize
from memory_contact import features, fit as unrestricted_fit


def basis_parts(t,y0,y1,te,ye,tw,tau,tau_m,origin,order,q0=0.):
    """Precompute bilinear factors once, rather than convolving in each solve."""
    B,U,qend=features(t,y0,y1,te,ye,tw,tau,tau_m,origin,q0=q0)
    m=np.arange(8)*.001
    w=-math.expm1(-2*math.pi*30*.001)*np.exp(-2*math.pi*30*.001*(7-np.arange(8)))
    v=(t[:,None]+m-origin)/tw
    local_time=(t[:,None]+m-te)/tw
    elapsed=t[:,None]+m-origin
    # Memory of a surface polynomial v^p, initially zero; independent exp
    # column profiles unknown initial relative memory, including encoder state.
    memory=[np.zeros_like(v)]
    for p in range(1,order+1):
        # ∫ exp(-(t-u)/tau) d(u/tw)^p, polynomial convolution recursion.
        memory.append(p*tau_m/tw*(v**(p-1)-memory[p-1]))
        # Previous recurrence includes the initial-value term only for p=1.
        if p==1:
            memory[p]-=tau_m/tw*np.exp(-elapsed/tau_m)
    derivatives=np.zeros((3,len(t),order+2))
    for p in range(order+1):
        derivatives[0,:,p]=-(v**p)@w
        derivatives[1,:,p]=-(local_time*v**p)@w
        derivatives[2,:,p]=-memory[p]@w
    base=np.zeros((len(t),order+2));base[:,-1]=np.exp(-elapsed/tau_m)@w
    return base,derivatives,U,B,qend


def basis(theta,t,y0,y1,te,ye,tw,tau,tau_m,origin,order):
    base,derivatives,U,B,qend=basis_parts(t,y0,y1,te,ye,tw,tau,tau_m,origin,order)
    D=base+np.einsum('i,ijk->jk',theta,derivatives)
    return D,U@np.asarray(theta),U,B,qend


def fit(measured,nw=32,tau=.019,noise=.04,tau_m=.04,order=2,max_nfev=30):
    if measured.ndim!=2 or measured.shape[1]!=7:
        raise ValueError('seven measured channels required')
    n=len(measured);j=np.arange(max(1,n-1-nw),n-1)
    if len(j)<16 or np.any(measured[j,3]<=.5) or np.any(measured[j+1,3]<=.5):
        return None
    te,ye,tw,origin=float(measured[-1,0]),float(measured[-1,2]),nw*.008,float(measured[j[0],0])
    t,y0,y1=measured[j,0],measured[j,2],measured[j+1,2]
    z=measured[j+1,3]-math.exp(-2*math.pi*30*.008)*measured[j,3]
    initial=unrestricted_fit(measured,model='memoryc_040',nw=nw,tau=tau,noise=noise)
    if initial is None:
        return None
    # All relaxation hypotheses here are explicit; default initialization uses
    # the existing 40-ms fit, never a plant parameter or scoring stiffness.
    p0=np.asarray(initial['coefficients'])
    endpoints=np.clip([p0[0]-p0[1],p0[0],p0[2]], [200,200,0], [3e5,3e5,3e5])
    transform=np.array([[0.,1.,0.],[-1.,1.,0.],[0.,0.,1.]])
    base,derivatives,U,B,q=basis_parts(t,y0,y1,te,ye,tw,tau,tau_m,origin,order)
    cached={}
    def components(parameters):
        if 'parameters' in cached and np.array_equal(parameters,cached['parameters']):
            return cached['parts']
        theta=transform@parameters*1e4
        D=base+np.einsum('i,ijk->jk',theta,derivatives);mean=U@theta
        norms=np.linalg.norm(D,axis=0)
        Dplus=np.linalg.pinv(D/norms,rcond=1e-10)/norms[:,None]
        beta=Dplus@(z-mean);e=mean+D@beta-z
        raw=U+np.einsum('ijk,k->ji',derivatives,beta)
        profile_jac=raw-D@(Dplus@raw)-Dplus.T@np.einsum('ijk,j->ki',derivatives,e)
        parts=(theta,D,mean,beta,e,raw,profile_jac)
        cached.update(parameters=parameters.copy(),parts=parts)
        return parts
    def residual(parameters):
        return components(parameters)[4]/noise
    def jacobian(parameters):
        return components(parameters)[6]@transform*(1e4/noise)
    solution=least_squares(residual,endpoints/1e4,bounds=([.02,.02,0.],[30.,30.,30.]),
                           jac=jacobian,max_nfev=max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-6)
    theta,D,mean,beta,e,raw,profile_jac=components(solution.x)
    # Full local derivative includes the change of the SHARED surface response
    # with contact parameters. Profiling only U would overstate information.
    J=residualize(D,raw)
    norm=np.linalg.norm(J,axis=0)
    if min(norm)<1e-13:
        return None
    _,s,vt=np.linalg.svd(J/norm,full_matrices=False)
    if s[-1]<s[0]*1e-8:
        return None
    e=mean+D@beta-z;dof=max(len(z)-D.shape[1]-3,1)
    chi=float(e@e/(dof*noise**2));variance=noise**2*max(1.,chi)
    covariance=variance*((vt.T/s**2)@vt)/norm[:,None]/norm[None,:]
    omega=2*math.pi*20.8;u=(omega*tau_m)**2
    functional=np.array([1.,0.,u/(1+u)])
    gain=float(functional@theta);sigma=math.sqrt(max(float(functional@covariance@functional),0.))
    # Convert the fitted shared surface back into the original force-domain
    # nuisance span for forward prediction with the existing exact convolution.
    force_beta=np.linalg.lstsq(B,D@beta,rcond=1e-12)[0]
    physical=bool(theta[0]>0 and theta[0]-theta[1]>0 and theta[2]>=0)
    full=np.column_stack([raw,D]);fn=np.linalg.norm(full,axis=0)
    _,fs,fv=np.linalg.svd(full/fn,full_matrices=False)
    full_cov=variance*((fv.T/fs**2)@fv)/fn[:,None]/fn[None,:]
    return dict(k=gain,sigma=sigma,coefficients=theta,cov=covariance,beta=force_beta,
        te=te,ye=ye,tw=tw,origin=origin,qend=q,tau_m=tau_m,variance=variance,
        chi2=chi,dof=dof,physical=physical,min_gain=float(min(theta[0],theta[0]-theta[1])),
        risk=sigma/max(gain,1)*(1+math.sqrt(dof*max(chi-1,0))),
        surface_coefficients=beta,order=order,optimizer_success=bool(solution.success),
        nfev=int(solution.nfev),residual_rms=float(np.sqrt(np.mean(e**2))),
        surface_to_force_span_error=float(np.linalg.norm(B@force_beta-D@beta)),
        full_joint_cov=full_cov,fit_n=n,noise=noise)


def forecast_reference(r,measured,tau):
    """Causal post-fit innovation checks, conditional on observed encoders."""
    j=np.arange(r['fit_n']-1,len(measured)-1)
    if not len(j):
        return dict(score=float('inf'),count=0)
    base,derivatives,U,B,q=basis_parts(measured[j,0],measured[j,2],measured[j+1,2],
        r['te'],r['ye'],r['tw'],tau,r['tau_m'],r['origin'],r['order'],q0=r['qend'])
    D=base+np.einsum('i,ijk->jk',r['coefficients'],derivatives)
    beta=np.asarray(r['surface_coefficients'])
    pred=U@r['coefficients']+D@beta
    raw=U+np.einsum('ijk,k->ji',derivatives,beta)
    J=np.column_stack([raw,D])
    variance=r['noise']**2+np.einsum('ij,jk,ik->i',J,r['full_joint_cov'],J)
    observed=measured[j+1,3]-math.exp(-2*math.pi*30*.008)*measured[j,3]
    return dict(score=float(np.mean(((pred-observed)**2/variance)[-4:])),count=min(len(j),4))
