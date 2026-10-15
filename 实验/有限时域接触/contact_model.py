"""Cylindrical contact helpers extracted without changing their numerical formulas.
Source: 03_实验与后续工作/最小闭环验证_2026-10-12/run_experiment.py at b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c.
The failed collocation solver and its development artifacts are intentionally omitted.
"""
import math
import numpy as np

CFG = dict(R=.03, L=.02, lref=.02, E0=30000., Ef=45000., boundary=0.,
           a0=.003, speed=.008, s0=-.006, T=1.5, dt=.001,
           mass=.2, Br=.2, k0=150., damping0=8., energy0=.02,
           energy_min=.001, chi=.3, kmin=30., kmax=300.,
           dmin=2., dmax=20., krate=500., drate=50., alpha=10.,
           action_weight=.01, gain_weight=.1, seed=20261012,
           node_candidates=512, newton_tolerance=1e-9,
           newton_max_iterations=40, sensitivity_step=1e-5,
           candidate_bellman_residual_allow=.05,
           history_samples=12, history_interval=.005,
           timing_repeats=9,
           thresholds=dict(force_rmse_fraction=.1, cost_relative_extra=.05,
                           force_rmse_extra_fraction=.02, latency_reduction=.3))

def psi(a):
    h = np.sqrt(CFG['lref']**2-a*a)
    return (CFG['lref']**2-CFG['lref']*h+a*a*np.log((CFG['lref']+h)/a))/(2*CFG['R'])

DREF = float(psi(CFG['a0']))

FD = np.pi*CFG['L']*CFG['E0']*CFG['a0']**2/(4*CFG['R'])

def halfwidth(delta):
    """Vector bracketed Newton for the monotone cylindrical indentation law."""
    delta = np.asarray(delta)
    if np.any(delta <= 0) or np.any(delta >= psi(CFG['lref']*.999)):
        raise ValueError('Indentation outside the calibrated contact branch')
    lo = np.full_like(delta, 1e-12, dtype=float)
    hi = np.full_like(delta, CFG['lref']*.999, dtype=float)
    a = np.minimum(CFG['a0']*np.sqrt(delta/DREF), hi*.95)
    for _ in range(15):
        h = np.sqrt(CFG['lref']**2-a*a)
        log = np.log((CFG['lref']+h)/a)
        val = (CFG['lref']**2-CFG['lref']*h+a*a*log)/(2*CFG['R'])
        err = val-delta
        if np.max(np.abs(err)) < 2e-15:
            break
        lo = np.where(err < 0, a, lo)
        hi = np.where(err >= 0, a, hi)
        nxt = a-err/(a/CFG['R']*log)
        update = np.where((nxt > lo)&(nxt < hi), nxt, .5*(lo+hi))
        a = np.where(np.abs(err)<2e-15, a, update)
    if np.max(np.abs(psi(a)-delta)) > 1e-12:
        raise ValueError('Geometry inversion did not converge')
    return a

def contact(t, delta, eta, moments=False):
    """eta=(Ef/E0,b/a0); raw My includes known tangential traction."""
    a = halfwidth(delta)
    s = CFG['s0']+CFG['speed']*np.asarray(t)
    r = eta[1]*CFG['a0']-s
    z = np.clip(r/a, -1., 1.)
    theta = np.arccos(z)
    h = np.sqrt(np.maximum(a*a-r*r, 0.))
    area = .5*(a*a*theta-np.clip(r,-a,a)*h)
    efront = eta[0]*CFG['E0']
    const = CFG['L']/(2*CFG['R'])
    fn = const*(CFG['E0']*np.pi*a*a/2+(efront-CFG['E0'])*area)
    slope_a = const*a*(CFG['E0']*np.pi+(efront-CFG['E0'])*theta)
    ah = np.sqrt(CFG['lref']**2-a*a)
    da_dd = CFG['R']/(a*np.log((CFG['lref']+ah)/a))
    slope_d = slope_a*da_dd
    if not moments:
        return fn, slope_d
    m1 = const*(efront-CFG['E0'])*h**3/3
    i2 = np.pi*a**4/8
    i2front = a**4/8*(theta+z*(1-2*z*z)*np.sqrt(np.maximum(1-z*z,0.)))
    tangential_moment = -.2/(2*CFG['R'])*const*(CFG['E0']*i2+(efront-CFG['E0'])*i2front)
    tangential_moment -= CFG['speed']*1.2/(2*CFG['R'])*a*a/4
    return fn, -m1+tangential_moment

def fast_width(d):
    if not 0<d<float(psi(CFG['lref']*.999)):raise ValueError('Lost contact / invalid indentation')
    a=CFG['a0']*math.sqrt(d/DREF)
    for _ in range(7):
        h=math.sqrt(CFG['lref']**2-a*a); lg=math.log((CFG['lref']+h)/a)
        err=(CFG['lref']**2-CFG['lref']*h+a*a*lg)/(2*CFG['R'])-d
        if abs(err)<2e-15:return a
        a-=err/(a/CFG['R']*lg)
    return a

def fast_force(t,d,eta):
    a=fast_width(d); r=eta[1]*CFG['a0']-CFG['s0']-CFG['speed']*t
    z=max(-1.,min(1.,r/a)); h=math.sqrt(max(a*a-r*r,0.))
    area=.5*(a*a*math.acos(z)-max(-a,min(a,r))*h)
    return CFG['L']/(2*CFG['R'])*(CFG['E0']*math.pi*a*a/2+(eta[0]-1)*CFG['E0']*area)
