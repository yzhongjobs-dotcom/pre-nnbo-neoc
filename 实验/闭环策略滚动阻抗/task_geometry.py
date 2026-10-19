"""Static task geometry, not a bound on the nonlinear finite-horizon problem."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent


def main():
    qf=1/.3**2;qz=1/.0015**2;kp=.5
    rows=[]
    for k in [600.,1000.,1500.,3000.,6000.]:
        delta=.001
        def cost(z):return qf*(k*(z-delta))**2+qz*z*z
        opt=minimize_scalar(cost,bounds=(-.01,.01),method='bounded',options={'xatol':1e-14})
        ratio=qz/(qf*k*k);kn=(1+kp)*qz/(qf*k)
        zstar=delta/(1+ratio);predicted=qz*delta**2/(1+ratio)
        assert np.isclose(opt.fun,predicted,rtol=1e-10)
        constant=58.3;alpha=constant/((1+kp)*k)
        static_excess=(alpha-ratio)**2/(ratio*(1+alpha)**2)
        rows.append({'contact_slope_N_m':k,'analytic_Kn_N_m':kn,'zstar_m':zstar,
                     'fractional_reduction_vs_perfect_force_tracking':ratio/(1+ratio),
                     'relative_static_excess_of_Kn58_3':static_excess})
    out={'rows':rows,'scope':'Frozen linear spring, quasi-static equilibrium, same world-frame reference; ignores motion, noise, delay, work, gain dynamics and energy. NOT an upper bound on dynamic control benefit.',
         'task_scale_ratio_N_m':float(np.sqrt(qz/qf)),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'task_geometry.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
