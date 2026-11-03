"""Independent physics and differential checks, separate from task evidence."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from dynamic_compliance import inverse_step,step_compliance,reconstruct_increment
from structured_memory import fit,basis_parts

HERE=Path(__file__).resolve().parent


def main():
    checks=[]
    ke,km,tm=12000.,6000.,.04
    for horizon in (.008,.032,.256):
        expected=1/ke+(1/(ke+km)-1/ke)*math.exp(-horizon*ke/((ke+km)*tm))
        actual=step_compliance([ke,0,km],tm,horizon)
        np.testing.assert_allclose(actual,expected,rtol=2e-5)
    # Exact FOH formula vs independent ODE integration, with nonzero state.
    E0,E1,c0=.13,-.07,.000012
    exact=solve_ivp(lambda t,c:[(E1-E0)/(.008*(ke+km))+
        (E0+(E1-E0)*t/.008)/((ke+km)*tm)-ke/((ke+km)*tm)*c[0]],
        [0,.008],[c0],rtol=1e-12,atol=1e-16).y[0,-1]
    np.testing.assert_allclose(inverse_step(c0,E0,E1,[ke,0,km],tm),exact,atol=1e-14)
    # Known displacement drives an independent forward SLS ODE. The inverse
    # receives force samples, not its true memory state. Matched initial c is
    # supplied only for this algebra test, never to the field controller.
    t=np.arange(501)*.008
    c=lambda z:.0001*np.sin(2*np.pi*2*z)+.00002*np.sin(2*np.pi*20.8*z)
    cd=lambda z:.0001*2*np.pi*2*np.cos(2*np.pi*2*z)+.00002*2*np.pi*20.8*np.cos(2*np.pi*20.8*z)
    state=solve_ivp(lambda z,r:[cd(z)-r[0]/tm],[0,t[-1]],[.00003],t_eval=t,rtol=1e-10,atol=1e-14)
    derivative=1000.;k=ke+derivative*t
    force_error=k*c(t)+km*state.y[0]
    inverse=[c(0.)]
    for j in range(1,len(t)):
        inverse.append(inverse_step(inverse[-1],force_error[j-1],force_error[j],
                                    [k[j],derivative*.256,km],tm))
    dynamic_rmse=float(np.sqrt(np.mean((inverse-c(t))**2)))
    scalar_rmse=float(np.sqrt(np.mean((force_error/(k+km)-c(t))**2)))
    assert dynamic_rmse<scalar_rmse/5
    checks.append(dict(independent_forward_time_varying_SLS=True,
                       matched_initial_condition_only=True,dynamic_rmse_m=dynamic_rmse,scalar_rmse_m=scalar_rmse))
    # Independently generate moving-surface physical force at 1 ms and apply
    # the known sensor filter. The fit sees only the seven measurement columns.
    t=np.arange(76)*.008;n=64;te=t[n-1];tw=.256;ts=.019
    y=.002+.00004*np.sin(2*np.pi*20.8*t)+.00002*np.sin(2*np.pi*7.3*t)+.0001*t
    measured=np.c_[t,.02*t,y,np.full(len(t),4.),np.zeros((len(t),3))]
    af=-math.expm1(-2*math.pi*30*.001);q=.000041
    surface=lambda z:.00175+.00007*z+.00003*z*z
    for j in range(len(t)-1):
        speed=(y[j+1]-y[j])/(ts*(-math.expm1(-.008/ts)))
        sol=solve_ivp(lambda h,r:[speed*np.exp(-h/ts)-(.00007+.00006*(t[j]+h))-r[0]/tm],
            [0,.008],[q],t_eval=np.r_[np.arange(8)*.001,.008],rtol=1e-11,atol=1e-14)
        force=measured[j,3]
        for i,h in enumerate(np.arange(8)*.001):
            yy=y[j]+speed*ts*(-math.expm1(-h/ts));k=ke+1000*(t[j]+h-te)/tw
            raw=k*(yy-surface(t[j]+h))+km*sol.y[0,i]
            force+=af*(raw-force)
        measured[j+1,3]=force;q=sol.y[0,-1]
    r=fit(measured[:n],tau=ts)
    np.testing.assert_allclose(r['coefficients'],[ke,1000,km],rtol=2e-5,atol=.05)
    assert r['surface_to_force_span_error']<1e-8
    # Profiled Jacobian at an intentionally nonzero residual: compare with
    # numerical differentiation, including motion-basis parameter dependence.
    j=np.arange(n-33,n-1);theta=np.array([12500.,800.,5500.])
    base,dD,U,B,q=basis_parts(t[j],y[j],y[j+1],te,y[n-1],tw,ts,tm,t[j[0]],2)
    z=measured[j+1,3]-math.exp(-2*math.pi*30*.008)*measured[j,3]+.02*np.cos(3.3*t[j])
    def profile(p):
        D=base+np.einsum('i,ijk->jk',p,dD);norm=np.linalg.norm(D,axis=0)
        pinv=np.linalg.pinv(D/norm,rcond=1e-10)/norm[:,None]
        beta=pinv@(z-U@p);e=U@p+D@beta-z
        raw=U+np.einsum('ijk,k->ji',dD,beta)
        jac=raw-D@(pinv@raw)-pinv.T@np.einsum('ijk,j->ki',dD,e)
        return e,jac
    error,jac=profile(theta);numerical=[]
    for k in range(3):
        plus=theta.copy();minus=theta.copy();plus[k]+=.01;minus[k]-=.01
        numerical.append((profile(plus)[0]-profile(minus)[0])/.02)
    numerical=np.column_stack(numerical)
    relative=float(np.linalg.norm(numerical-jac)/np.linalg.norm(numerical))
    assert relative<1e-5
    checks.append(dict(physical_shared_motion_recovery=r['coefficients'].tolist(),
                       independent_sensor_filter=True,profile_jacobian_relative_error=relative))
    uncertainty=reconstruct_increment(r,measured[:n],uncertainty=True)
    rng=np.random.default_rng(39177);noise_samples=[]
    for _ in range(1000):
        noisy=measured[:n].copy();noisy[:,3]+=rng.normal(0,.04,n)
        noise_samples.append(reconstruct_increment(r,noisy)['current'])
    observed=float(np.std(noise_samples,ddof=1))
    np.testing.assert_allclose(observed,uncertainty['sd_force'],rtol=.1)
    checks.append(dict(fixed_parameter_sensor_noise_influence=True,
        reference_sd_m=uncertainty['sd_force'],monte_carlo_sd_m=observed,
        joint_parameter_noise_correlation_not_assumed_independent=True,
        unobserved_initial_state_bias_not_covered=True))
    output=HERE/'results/shared_response_verification.json'
    output.write_text(json.dumps(dict(checks=checks,frozen_inverse_closed_form=True,
        exact_FOH_vs_ODE=True,not_task_performance_evidence=True),indent=2),encoding='utf8')
    print(output.read_text())


if __name__=='__main__':main()
