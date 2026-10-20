"""Pointwise moving-contact kinematics audit, NOT a new performance experiment.

At a fixed material coordinate X, a parabolic probe indentation has derivative
vz + (X-s)*vx/R - h_t. The benchmark kept only vz-h_t in its damping term.
This checks a missing physical term, not a full finite-face contact solution.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent
A,L,R,kappa,tau,vx=.01,.02,.05,6e6,.01,.012
FD=3.
z=FD/(2*A*L*kappa)+A*A/(6*R)
def base(u):return kappa*(z-u*u/(2*R))
def moving(u):return kappa*(z-u*u/(2*R)+tau*u*vx/R)
fn0=L*quad(base,-A,A,epsabs=1e-12)[0]
fn1=L*quad(moving,-A,A,epsabs=1e-12)[0]
m0=L*quad(lambda u:u*base(u),-A,A,epsabs=1e-14)[0]
m1=L*quad(lambda u:u*moving(u),-A,A,epsabs=1e-14)[0]
expected=L*kappa*tau*vx/R*(2*A**3/3)
assert min(moving(-A),moving(A))>0
assert abs(m1-expected)<1e-13
result={'scope':'Illustration of an omitted moving-footprint term, not real tissue or full controller falsification',
 'geometry':{'halfwidth_m':A,'length_m':L,'radius_m':R},'kappa_N_m3':kappa,'viscosity_time_s':tau,'speed_m_s':vx,'normal_position_m':z,
 'Fn_local_damping_only_N':fn0,'Fn_with_convective_rate_N':fn1,
 'normal_moment_local_only_Nm':m0,'normal_moment_with_convective_rate_Nm':m1,
 'apparent_gradient_if_mapped_by_3M_A2_N_m':3*m1/A**2,
 'apparent_spatial_force_change_over_predictor_horizon_N':3*m1/A**2*vx*.192,
 'qualifications':['Uniform medium and constant motion: resultant force has no spatial change although moment is nonzero',
  'Existing drift observer may absorb part of the apparent gradient; not a proved whole-controller prediction error',
  'Finite-face edge transitions, nonlocal elasticity and full power balance are not established by this calculation',
  'Velocity-dependent pressure is already represented in relevant prior such as AuSoScan; no novelty claim'],
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'contact_kinematics_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
