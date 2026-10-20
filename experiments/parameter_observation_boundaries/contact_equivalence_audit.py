"""Audit R15 Appendix D.2 using the complete Fn/Ft/My model.

This verifies a model-internal observational equivalence and its prediction
consequence. It is neither a new identifiability theorem nor a controller test.
No force limits, desired gain trajectories, or favorable controller costs used.
"""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent
R, L, A, E0, EF, V0 = .03, .02, .003, 30000., 45000., .003
CH = L / (2 * R)
MODEL_A = np.array([.20, .25, 2., 3.])  # mu0, muf, kd0, kdf
DMU = np.array([.020, -.025])

def integrals(a, r):
    r = np.clip(r, -a, a)
    g = np.sqrt(max(a*a-r*r, 0.))
    i0 = np.pi*a*a/2
    iff = .5*(a*a*np.arccos(r/a)-r*g)
    i1 = g**3/3
    i2 = np.pi*a**4/8
    i2f = a**4/8*np.arccos(r/a)+r*g/8*(a*a-2*r*r)
    return i0, iff, i1, i2, i2f

def observation(p, s, v, a=A):
    i0, iff, i1, i2, i2f = integrals(a, -s)
    mu0, muf, kd0, kdf = p
    sigma = np.sign(v)
    n0, nf = CH*E0*(i0-iff), CH*EF*iff
    alpha = iff/i0
    fn = n0+nf
    ft = -sigma*(mu0*n0+muf*nf)-v*((1-alpha)*kd0+alpha*kdf)
    my = -CH*(EF-E0)*i1 - sigma*CH/(2*R)*(mu0*E0*(i2-i2f)+muf*EF*i2f)
    my -= v/(2*R*i0)*(kd0*(i2-i2f)+kdf*i2f)
    return np.array([fn, ft, my])

def quadrature_observation(p, s, v, a=A):
    # Integrate each material separately; do not use the closed-form moments.
    nodes, weights = leggauss(100)
    values = np.zeros(3)
    boundary = np.clip(-s, -a, a)
    for lo, hi, e, mu, kd in [(-a,boundary,E0,p[0],p[2]),
                              (boundary,a,EF,p[1],p[3])]:
        if hi <= lo:
            continue
        # x=a*sin(phi) removes endpoint square-root singularities.
        plo, phi = np.arcsin(lo/a), np.arcsin(hi/a)
        angle = (phi-plo)/2*nodes+(phi+plo)/2
        x = a*np.sin(angle)
        jac = a*np.cos(angle)*(phi-plo)/2
        shape = a*np.cos(angle)
        qn = CH*e*shape
        qt = -np.sign(v)*mu*qn - v*kd*shape/(np.pi*a*a/2)
        values += np.array([np.sum(weights*jac*qn),
                            np.sum(weights*jac*qt),
                            np.sum(weights*jac*(x*x/(2*R)*qt-x*qn))])
    return values

def regressor(s, v, a=A):
    base = observation(np.zeros(4), s, v, a)
    return np.column_stack([observation(np.eye(4)[j],s,v,a)-base
                            for j in range(4)])

def main():
    i0 = np.pi*A*A/2
    shift = np.r_[DMU, -np.sign(V0)*CH*i0*np.array([E0,EF])/V0*DMU]
    model_b = MODEL_A + shift
    assert np.all(model_b > 0)
    positions = np.linspace(-1.2*A,1.2*A,121)
    traces = np.array([[s,*observation(MODEL_A,s,V0),
                        *observation(model_b,s,V0)] for s in positions])
    differences = np.abs(traces[:,1:4]-traces[:,4:7]).max(axis=0)
    history = np.vstack([regressor(s,V0) for s in positions])
    # Fixed scaling used only to report numerical rank with comparable columns.
    row_scale = np.tile([1.,1.,1/A], len(positions))
    column_scale = np.array([1.,1.,CH*E0*i0/V0,CH*EF*i0/V0])
    scaled = row_scale[:,None]*history*column_scale[None,:]
    singular = np.linalg.svd(scaled, compute_uv=False)
    rank = int(np.sum(singular > singular[0]*1e-10))
    assert rank == 2
    future = []
    for factor in [.5,1.,2.,3.]:
        ya, yb = observation(MODEL_A,0.,V0*factor), observation(model_b,0.,V0*factor)
        future.append(dict(velocity_factor=factor, observation_a=ya.tolist(),
                           observation_b=yb.tolist(), difference=(yb-ya).tolist()))
    footprint = []
    for factor in [.8,1.,1.2]:
        ya, yb = observation(MODEL_A,0.,V0,A*factor), observation(model_b,0.,V0,A*factor)
        footprint.append(dict(halfwidth_factor=factor,
                              observation_difference=(yb-ya).tolist()))
    excited = np.vstack([regressor(s,V0*f) for f in [1.,2.] for s in positions])
    sr = np.tile([1.,1.,1/A],2*len(positions))
    se = np.linalg.svd(sr[:,None]*excited*column_scale[None,:],compute_uv=False)
    excited_rank = int(np.sum(se > se[0]*1e-10))
    assert excited_rank == 4
    qerror = np.max([np.abs(observation(p,s,V0)-quadrature_observation(p,s,V0))
                    for p in [MODEL_A,model_b] for s in positions],axis=0)
    assert differences.max() < 1e-12 and qerror.max() < 1e-12
    assert abs(future[-1]['difference'][1]) > 1e-3
    out=ROOT/'contact_equivalence_audit'
    out.mkdir(exist_ok=True)
    np.savetxt(out/'identical_history.csv',traces,delimiter=',',
               header='s,Fn_A,Ft_A,My_A,Fn_B,Ft_B,My_B',comments='')
    result=dict(model=dict(R=R,L=L,a=A,E0=E0,Ef=EF,v0=V0,
                           known_elasticity_geometry_boundary=True,
                           model_a=MODEL_A.tolist(),model_b=model_b.tolist()),
        original_source='R15 appendix D.2 and main equations (12)-(15). Existing degeneracy; this audit adds complete-wrench replication and counterfactual predictions.',
        history=dict(samples=len(positions), maximum_difference_Fn_Ft_My=differences.tolist(),
                     regressor_rank=rank,scaled_singular_values=singular.tolist()),
        independent_quadrature_max_error=qerror.tolist(),
        changed_speed_prediction=future,changed_footprint_prediction=footprint,
        two_speed_identification=dict(rank=excited_rank,scaled_singular_values=se.tolist()),
        scope='An exact structural limitation of this reduced contact model under constant footprint and velocity. No controller failure, physiological parameter calibration, online timing, noise robustness, or novelty claimed.')
    (out/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
