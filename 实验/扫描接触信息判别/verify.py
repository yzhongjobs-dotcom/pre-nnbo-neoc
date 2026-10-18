"""Independent force integrals, finite differences, twin worlds and replay ODE."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
import screen as s

HERE = Path(__file__).resolve().parent


def independent_wrench(x, location, time, p):
    def density(u):
        xx = location+u
        k = p[0]*np.exp(p[1]*np.tanh((xx-p[2])/p[3]))
        h = (p[9]+p[4]*np.sin(2*np.pi*xx/p[5]+p[6])
             +p[7]*np.sin(2*np.pi*p[8]*time+p[6])+p[10]*(xx/.08)**2)
        return s.L*k*max(x[0]+x[1]*u/s.A-u*u/.1-h, 0.)
    return np.array([quad(density, -s.A, s.A, epsabs=1e-10, limit=200)[0],
                     quad(lambda u: u/s.A*density(u), -s.A, s.A, epsabs=1e-10, limit=200)[0]])


def main():
    data = np.load(HERE/'paired.npz')
    report = json.loads((HERE/'results.json').read_text())
    errors, jerrors, independent_loop = [], [], []
    rng = np.random.default_rng(992)
    for family in ('uniform_static','heterogeneous_static','heterogeneous_moving'):
        ps = data[family+'_parameters']
        for i in range(4):
            trace = data[family+'_trace_fixed_PI']
            for row in trace[rng.integers(0, len(trace), 8)]:
                x = np.array([row[2], s.A*row[3]])
                tm, location = row[0]+s.DT, row[1]+s.SPEED*s.DT
                exact = independent_wrench(x, location, tm, ps[i])
                coarse, jac = s.wrench(x, location, tm, ps[i], 61)
                errors.append(float(np.max(np.abs(coarse-exact))))
                eps = 1e-7
                fd = np.column_stack([(independent_wrench(x+eps*np.eye(2)[k], location, tm, ps[i])
                        -independent_wrench(x-eps*np.eye(2)[k], location, tm, ps[i]))/(2*eps)
                        for k in range(2)])
                jerrors.append(float(np.linalg.norm(fd-jac)/max(1.,np.linalg.norm(fd))))
        # Separately implemented fixed PI closed loop; independent quadrature
        # supplies its own feedback. No replay of original commands/states.
        fi = ('uniform_static','heterogeneous_static','heterogeneous_moving').index(family)
        noise = s.noise_for(8000+fi,128)[0]
        state = np.r_[data[family+'_initial_poses'][0],np.zeros(2)]
        integral = np.zeros(2)
        observations = []
        cfg = report['selected']['fixed_PI']
        positions = []
        loss = 0.
        for q in range(1000):
            tm = q*s.DT
            measured = independent_wrench(state[:2],s.SPEED*tm,tm,ps[0])+noise[q]
            observations.append(measured)
            err = np.array([s.FD,0.])-observations[max(0,q-2)]
            integral = np.clip(integral+s.DT*err,-2.,2.)
            cmd = (cfg['bandwidth']*err+cfg['integral']*cfg['bandwidth']**2*integral)/[1500.,500.]
            cmd = np.clip(cmd,[-.020,-.0075],[.020,.0075])
            sol = solve_ivp(lambda t,y: np.r_[y[2:],(cmd-y[2:])/[.020,.030]],
                [tm,tm+s.DT],state,method='DOP853',rtol=1e-11,atol=1e-13,dense_output=True)
            for offset in (.0025,.0075):
                tt = tm+offset
                ff = independent_wrench(sol.sol(tt)[:2],s.SPEED*tt,tt,ps[0])
                loss += .005*((ff[0]-s.FD)**2+ff[1]**2)
            state = sol.y[:,-1]
            positions.append(state[:2])
        original = data[family+'_trace_fixed_PI']
        original_x = np.column_stack((original[:,2],s.A*original[:,3]))
        original_loss = float(s.score(data[family+'_metrics'][0,0]))
        independent_loop.append(dict(family=family,max_position_difference=float(np.max(np.abs(np.array(positions)-original_x))),
            original_loss=original_loss,independent_loss=loss,relative_loss_difference=loss/original_loss-1))
    # Two flat uniform surfaces in the actual kappa/height range, full footprint.
    a, L, F, radius = s.A, s.L, s.FD, .05
    area, cbar = 2*a*L, a*a/(6*radius)
    k1, k2 = 3e6, 6e6
    h1, h2 = -F/(4*k1*area), F/(4*k1*area)
    z0 = h1+F/(k1*area)+cbar
    def analytic(k, h, z):
        return area*k*(z-h-cbar)
    twins = dict(kappas=[k1,k2], surface_heights=[h1,h2], identical_pose=z0,
        forces=[analytic(k1,h1,z0),analytic(k2,h2,z0)], moments=[0.,0.],
        min_indentation=[z0-h1-a*a/(2*radius),z0-h2-a*a/(2*radius)],
        local_normal_slopes=[area*k1,area*k2],
        forces_after_common_0_1mm_step=[analytic(k1,h1,z0+.0001),analytic(k2,h2,z0+.0001)])
    # Same controller, observation cadence and noise; refine only mechanics/metrics.
    refinement = {}
    for family_index, family in enumerate(('uniform_static','heterogeneous_static','heterogeneous_moving')):
        ps, xs = data[family+'_parameters'][:8], data[family+'_initial_poses'][:8]
        nn = s.noise_for(8000+family_index,128)[:8]
        rows = {}
        for m,name in enumerate(s.NAMES):
            cfg = report['selected'][name]
            fine = s.batch(ps,xs,nn,m,cfg['bandwidth'],cfg['integral'],nq=241,substeps=20)
            coarse = data[family+'_metrics'][m,:8]
            rows[name] = dict(mean_loss_coarse=float(s.score(coarse).mean()),
                mean_loss_fine=float(s.score(fine).mean()),
                relative_change_mean=float((s.score(fine).mean()-s.score(coarse).mean())/s.score(coarse).mean()),
                max_absolute_peak_difference=float(np.max(np.abs(fine[:,2]-coarse[:,2]))))
        refinement[family] = rows
    results = dict(source_sha256=hashlib.sha256((HERE/'screen.py').read_bytes()).hexdigest(),
        verify_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        independent_quad_cases=len(errors), max_normalized_wrench_abs_error=max(errors),
        max_J_relative_difference_vs_independent_fd=max(jerrors),
        independent_fixed_PI_closed_loop=independent_loop, twins=twins, refinement=refinement,
        note="Independent PI loop uses DOP853 actuators, adaptive force quadrature and separately coded feedback. Refined other loops keep same sensing/controller cadence. Models are elastic/synthetic.")
    (HERE/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__ == '__main__':
    main()
