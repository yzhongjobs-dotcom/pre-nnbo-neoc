"""Independent derivative/variance checks for the used contact observation.

Synthetic contact algebra, never task evidence. Covariance excludes uncertainty
in the assumed memory parameters/geometry. Verifies functional null-space checks
and the scalar analytic limit, not global statistical coverage.
"""
import json
import math
from pathlib import Path
import numpy as np
from dynamic_contact import ContactFit,Memory


def main():
    linear=ContactFit(memory=Memory(.5,.2))
    var=.01**2;k=12000.
    fit=dict(theta=np.array([math.log(k)]),data_cov=np.array([[var]]),
             data_null=np.zeros((1,0)),data_shift=np.zeros(1))
    cur=np.array([1.,.01,.002,4.2,4.1]);prev=np.array([.992,.009,.002,3.9,4.])
    result=linear.correction_increment(fit,cur,prev);expected=(4.2-3.9)/(1.5*k)
    np.testing.assert_allclose(result['increment'],expected,rtol=1e-10)
    np.testing.assert_allclose(result['gradient'],[-expected],rtol=1e-8)
    np.testing.assert_allclose(result['variance'],expected**2*var,rtol=1e-8)
    null=dict(fit,data_null=np.ones((1,1)))
    assert linear.correction_increment(null,cur,prev)['unobserved']
    same=prev.copy();same[3]=cur[3]
    assert not linear.correction_increment(null,cur,same)['unobserved']
    rng=np.random.default_rng(73601)
    sampled=expected*np.exp(-rng.normal(0.,np.sqrt(var),100000))
    empirical=float(np.var(sampled));assert abs(empirical/result['variance']-1)<.02
    # Shared contact parameters across a sequence create correlated errors.
    # Covariance of the cumulative observation uses the SUM of Jacobians,
    # not a per-increment variance and not a sum of independent variances.
    forces=np.linspace(4.,4.8,9);increments=np.diff(forces)/(1.5*k)
    accumulated=float(forces[-1]-forces[0])/(1.5*k)
    sequence_sample=accumulated*np.exp(-rng.normal(0.,np.sqrt(var),100000))
    level_variance=accumulated**2*var
    assert abs(float(np.var(sequence_sample))/level_variance-1)<.02
    wrong_independent_variance=float(np.sum(increments**2)*var)
    assert abs(level_variance/wrong_independent_variance-8)<1e-10
    rows=[]
    for contact in ('foundation','line'):
        model=ContactFit(spatial=True,contact=contact)
        theta=np.array([math.log(18000. if contact=='foundation' else 5e5),math.log(3.),.003,.012])
        fit=dict(theta=theta,data_cov=np.eye(4)*.01,data_null=np.zeros((4,0)),data_shift=np.zeros(4))
        inc=model.correction_increment(fit,cur,prev);g=[]
        for i in range(4):
            h=1e-4 if i<2 else 1e-6;tp,tm=theta.copy(),theta.copy();tp[i]+=h;tm[i]-=h
            values=[]
            for th in (tp,tm):
                d=dict(fit,theta=th)
                values.append(model.target_correction(d,cur[1],cur[3],cur[4])-
                              model.target_correction(d,prev[1],prev[3],prev[4]))
            g.append((values[0]-values[1])/(2*h)*(1. if i<2 else .01))
        err=float(np.linalg.norm(inc['gradient']-g)/np.linalg.norm(g));assert err<1e-4,(contact,err)
        rows.append(dict(contact=contact,implicit_gradient=inc['gradient'].tolist(),
                         finite_difference_gradient=g,relative_difference=err))
    out=dict(scalar=dict(increment=expected,propagated_variance=result['variance'],
                         empirical_variance=empirical,null_space_and_zero_increment_checks=True),
             nonlinear_gradient_checks=rows)
    out['shared_parameter_sequence']=dict(level_variance=level_variance,
        empirical_variance=float(np.var(sequence_sample)),
        independent_increment_sum=wrong_independent_variance,
        missed_correlation_factor=level_variance/wrong_independent_variance)
    path=Path(__file__).resolve().parent/'results'/'contact_functional_audit.json'
    path.write_text(json.dumps(out,indent=2),encoding='utf8');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
