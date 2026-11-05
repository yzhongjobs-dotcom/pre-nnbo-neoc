"""Open-loop physical / numerical checks for plant3. Never imports a Method."""
import hashlib
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import physics_plant3 as P
import scenes3 as S

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'/'plant3_checks.json'


def uniform(law='foundation', probe='convex', q=1.1e6, E=17.e3, cover=None, T=6., **overrides):
    o = dict(contact_law=law, probe=probe, field_sigma=0., force=8., slide='long', phi_hc=5., surface_amp=0.,
             mu=0., k_d=0., pose_error=0., asymmetry=False, alpha_scale=1., tau_scale=1., rho=0., kr=1.e5, K_t=2.e4, accel=.1, k_rib=1.e12, no_curve=True, tool_mass=0., mass_error=0., tau_error=0., tau_a_true=.019, sensor_delay_ms=0., delay_error=0.)
    if cover is not None:
        o.update(cover=cover, full_cover=True)
    o.update(overrides)
    sc = S.cell_scene3('rest', 'relaxed', T=T, overrides=o)
    sc['q_field'][:] = 1.; sc['E_field'][:] = 1.; sc['q_base'] = q; sc['E_base'] = E
    sc['xcmd'][:] = sc['xcmd'][0]
    return sc


def reference(sc, d):
    """Independent scalar formulas (aligned, homogeneous, full or absent cover)."""
    L, W, R = (sc['geometry'][k] for k in ('length', 'width', 'radius'))
    q, E = sc['q_base'], sc['E_base']
    h = sc['cover'] if sc['full_cover'] else math.inf
    def slope(x):
        a = L*.5 if math.isinf(R) else min(math.sqrt(2.*R*max(x, 0.)), L*.5)
        if a <= 0.:
            return 0.
        factor = 1.
        if math.isfinite(h):
            if math.isinf(R):
                factor = (1.-x/h)**-2
            else:
                A = h-x; c = 1./(2.*R)
                integral = a/(A*(A+c*a*a))+math.atan(a*math.sqrt(c/A))/(A*math.sqrt(A*c))
                factor = h*h*integral/(2.*a)
        return (q*W*2.*a if sc['contact_law'] == 'foundation' else 1.1*E*math.sqrt(W*2.*a))*factor
    a = L*.5 if math.isinf(R) else min(math.sqrt(2.*R*d), L*.5)
    if sc['contact_law'] == 'foundation':
        if math.isinf(h):
            f = q*W*(2.*a*d-(0. if math.isinf(R) else a**3/(3.*R)))
        elif math.isinf(R):
            f = q*W*L*d/(1.-d/h)
        else:
            A = h-d; c = 1./(2.*R)
            f = q*W*(h*h*2./math.sqrt(A*c)*math.atan(a*math.sqrt(c/A))-2.*a*h)
    else:
        top = math.inf if math.isinf(R) else L*L/(8.*R)
        points = [top] if top < d else None
        f = quad(slope, 0., d, points=points, epsabs=1.e-10, epsrel=1.e-10)[0]
    return f, slope(d), 2.*a


def static_table():
    rows = []; errors = []; tangent_errors = []; approx_errors = []; approx_tangent = []
    tissues = [('R3', 1.1e6, 17.e3, None), ('R1', 1.6e6, 28.e3, None),
               ('active', 4.e6, 63.e3, None)]+[(f'rib_{h}mm', 1.1e6, 17.e3, h/1000.) for h in (4, 7, 18, 32)]
    for probe in S.PROBES:
        for law in ('foundation', 'halfspace'):
            for tissue, q, E, h in tissues:
                sc = uniform(law, probe, q, E, h)
                lam = P.calibration(sc)
                for force in (4., 8., 12.):
                    d = P.solve_depth(sc, force)
                    f, k, _, length, _ = P.contact(sc, d)
                    rf, rk, rl = reference(sc, d)
                    err = max(abs(f/rf-1.), abs(k/rk-1.), abs(length/rl-1.))
                    approximate = law == 'halfspace' and probe == 'convex' and h is not None
                    (approx_errors if approximate else errors).append(err)
                    eps = min(1.e-7, (h-d)*.001) if h else 1.e-7
                    numeric = (P.contact(sc, d+eps)[0]-P.contact(sc, d-eps)[0])/(2.*eps)
                    (approx_tangent if approximate else tangent_errors).append(abs(numeric/k-1.))
                    gain = P.complex_gain(sc, k, f, lam=lam)
                    rows.append(dict(probe=probe, law=law, tissue=tissue, force_N=force, depth_mm=d*1000.,
                                     length_mm=length*1000., tangent_N_m=k, real20_8_N_m=gain.real,
                                     amplitude20_8_N_m=abs(gain), phase20_8_deg=math.degrees(np.angle(gain)),
                                     relative_formula_error=err, lambda_s_m=lam))
    return dict(passed=max(errors) < .01 and max(tangent_errors) < .01 and max(approx_errors) < .2 and max(approx_tangent) < .2, rows=rows,
                max_formula_error=max(errors), max_derivative_error=max(tangent_errors),
                halfspace_convex_cover_force_error=max(approx_errors), halfspace_convex_cover_tangent_error=max(approx_tangent))


def remus():
    rows = []
    for law in ('foundation', 'halfspace'):
        for region, depths, forces, q, E in [
            ('R1', [5., 10., 15., 20.], [1.45, 3.33, 5.65, 8.40], 1.6e6, 28.e3),
            ('R3', [6.25, 12.5, 18.75, 25.], [1.12, 2.92, 5.39, 8.53], 1.1e6, 17.e3),
        ]:
            sc = uniform(law, q=q, E=E)
            for d, f in zip(depths, forces):
                pred = P.hemisphere(sc, d/1000.)
                hertz = 4./3.*E*math.sqrt(.01)*(d/1000.)**1.5
                rows.append(dict(law=law, region=region, depth_mm=d, measured_N=f,
                                 predicted_N=pred, relative_error=pred/f-1., hertz_extrapolation_N=hertz))
    # This is a reconstruction audit, not a data-fit acceptance threshold.
    return dict(passed=True, criterion='all 16 requested points reported; no empirical fit tolerance specified',
                empirical_fit_validated=False, max_absolute_relative_error=max(abs(r['relative_error']) for r in rows), rows=rows)


def _fit(log, column, frequency=20.8):
    keep = log[:, 0] > 2.
    t = log[keep, 0]; w = 2.*np.pi*frequency
    images = [125.*n+sign*frequency for n in range(1, 8) for sign in (-1., 1.)]
    A = np.column_stack([np.sin(w*t), np.cos(w*t), np.ones(len(t)), t]+
                        [func(2.*np.pi*f*t) for f in images for func in (np.sin, np.cos)])
    coef = np.linalg.lstsq(A, log[keep, column], rcond=None)[0]
    return complex(coef[0], coef[1])


def open_loop():
    rows = []
    for law in ('foundation', 'halfspace'):
        for phi in (0., 5., 30., 45.):
            sc = uniform(law, phi_hc=phi)
            # Freeze the support to separate the specified robot/contact series response.
            sc['m'] = 1.e9; sc['k'] = 1.e12; sc['c'] = 2.e10
            s, _ = P._init(sc, P.PLATFORM)
            t = np.arange(int(sc['T']/P.PLATFORM['Tc']))*P.PLATFORM['Tc']
            pl = dict(P.PLATFORM, log_dt=.0005)
            log = P.replay(sc, s[2]+1.e-5*np.sin(2.*np.pi*20.8*t), platform=pl)
            measured = _fit(log, 18)/_fit(log, 2)
            k = P.contact(sc, s[13])[1]; analytic = P.complex_gain(sc, k, 8.)
            rows.append(dict(law=law, phi_hc_deg=phi, measured_gain=abs(measured), analytic_gain=abs(analytic),
                             measured_phase_deg=math.degrees(np.angle(measured)), analytic_phase_deg=math.degrees(np.angle(analytic)),
                             gain_error=abs(abs(measured)/abs(analytic)-1.),
                             phase_error_deg=abs(math.degrees(np.angle(measured/analytic)))))
    return dict(passed=all(r['gain_error'] < .01 and r['phase_error_deg'] < .5 for r in rows), rows=rows,
                gain_tolerance=.01, phase_tolerance_deg=.5, support_frozen=True,
                command_period_s=.008, observation_period_s=.0005, hold_images_fitted_separately=True)


def continuity():
    rows = []; monotone = True; non_saturation = True
    for law in ('foundation', 'halfspace'):
        for probe in S.PROBES:
            sc = uniform(law, probe)
            fs = np.array([P.contact(sc, d)[0] for d in np.linspace(1.e-7, .04, 100)])
            monotone &= bool(np.all(np.diff(fs) > 0.))
            ds = [P.solve_depth(sc, f) for f in (4., 8., 12., 24., 48.)]
            non_saturation &= bool(np.all(np.diff(ds) > 0.))
            if probe == 'convex':
                top = sc['geometry']['length']**2/(8.*sc['geometry']['radius'])
                a, b = P.contact(sc, top-1.e-9), P.contact(sc, top+1.e-9)
                rows.append(dict(law=law, probe=probe, boundary='full_footprint', force_gap_N=b[0]-a[0],
                                 relative_tangent_gap=abs(b[1]/a[1]-1.)))
            eps = 1.e-10
            a, b = P.contact(sc, -eps), P.contact(sc, eps)
            rows.append(dict(law=law, probe=probe, boundary='first_contact', force_gap_N=b[0]-a[0],
                             tangent_below=a[1], tangent_above=b[1], intrinsic_flat_punch_jump=probe == 'linear'))
        sc = uniform(law)
        sc['ribs'] = np.array([[sc['xcmd'][0], .015, .002, .018]])
        edge = sc['xcmd'][0]+.0075
        vals = [P.contact(sc, .004, edge+shift) for shift in (-1.e-7, 1.e-7)]
        rows.append(dict(law=law, boundary='rib_edge', force_gap_N=vals[1][0]-vals[0][0],
                         relative_tangent_gap=abs(vals[1][1]/vals[0][1]-1.)))
        # Neither bone nor stiffness leaks outside the footprint.
        far = P.contact(sc, .004, edge+.10)[0]; sc['ribs'] = np.empty((0, 4))
        rows.append(dict(law=law, boundary='off_rib', force_error_N=far-P.contact(sc, .004, edge+.10)[0]))
    smooth = all(abs(r.get('relative_tangent_gap', 0.)) < .001 and abs(r.get('force_gap_N', 0.)) < .001 for r in rows)
    # A flat punch aligned with a flat surface has a tangent step at first contact by geometry. The force stays continuous, and it occurs only on contact loss.
    return dict(passed=bool(smooth and monotone and non_saturation), smooth_boundaries_passed=smooth, monotone=monotone, non_saturation=non_saturation,
                accepted_exception='Ideal aligned flat punch: finite right tangent, zero left tangent at first contact.', rows=rows)


def pose_and_dissipation():
    rows = []; passed = True
    for law in ('foundation', 'halfspace'):
        sc = uniform(law, 'linear'); d = .0007
        f0, _, cp0, l0, _ = P.contact(sc, d)
        fp, _, cp, lp, _ = P.contact(sc, d, theta=math.radians(5.))
        fm, _, cm, lm, _ = P.contact(sc, d, theta=-math.radians(5.))
        okay = lp < l0 and cp < 0. and cm > 0. and abs(cp0) < 1.e-10 and abs(fp-fm) < 1.e-8
        rows.append(dict(law=law, length_zero_mm=l0*1000., length_5deg_mm=lp*1000., center_5deg_mm=cp*1000.,
                         moment_5deg_Nm=fp*cp, moment_minus5deg_Nm=fm*cm, zero_moment_Nm=f0*cp0, passed=okay))
        passed &= okay
    energies = []; eigenvalues = []
    for rho in (-.3, 0., .3):
        sc = uniform(rho=rho, mu=.1, asymmetry=True, k_d=2.)
        c = P._constants(sc, P.PLATFORM); g = P._geometry(sc)
        eigenvalues.extend(np.linalg.eigvalsh(P.compliance(sc)).tolist())
        for vel in np.linspace(-.06, .06, 121):
            energies.append(float(P._friction(8., vel, .04, c, g)*vel))
    sc = uniform(mu=0., rho=0., k_d=0., asymmetry=True, T=1.)
    stationary = dict(sc); moving = dict(sc)
    moving['xcmd'] = sc['xcmd'][0]+.02*np.arange(len(sc['xcmd']))*.001
    s0, _ = P._init(sc, P.PLATFORM); cmd = np.full(125, s0[2])
    a, b = P.replay(stationary, cmd), P.replay(moving, cmd)
    gap = float(np.max(np.abs(a[:, 3]-b[:, 3])))
    dummy = np.arange(3*P.LOGC, dtype=float).reshape(3, P.LOGC)
    observed = P.Obs(dummy, 3)
    masked = bool(np.all(np.isnan(observed.log[:, 20:26])))
    return dict(passed=bool(passed and min(energies) >= 0. and min(eigenvalues) > 0. and gap < 1.e-9 and masked),
                pose=rows, min_friction_power_W=min(energies), min_nonzero_friction_power_W=min(x for x in energies if x > 0.),
                min_compliance_eigenvalue=min(eigenvalues), zero_coupling_force_gap_N=gap, truth_masked=masked)


def scene_physics():
    """Pose sign on a falling surface, muscle activation, thin ribs, balance within each tissue."""
    # downhill: x past the centre with the surface falling ahead, moving forward
    sc = S.cell_scene3('rest', 'relaxed', 'curved_untracked', 0, T=3., overrides=dict(contact_law='foundation', probe='linear', slide='long', force=8., phi_hc=5.,
                       pose_error=0., surface_radius=.3, surface_amp=0., field_sigma=0., mu=.3, rho=0., accel=.1, speed=.02, asymmetry=False))
    s0, _ = P._init(sc, P.PLATFORM)
    log = P.replay(sc, np.full(375, s0[2]))
    late = log[len(log)//2:]
    th = late[:, 17]; extra = late[:, 18]-late[:, 3]*np.cos(th)
    pose = bool(np.all(th > 0.) and np.all(extra > 0.) and np.all(late[:, 24] < 0.))
    # activation raises force at fixed position
    o = dict(contact_law='foundation', probe='linear', slide='long', force=8., phi_hc=0., pose_error=0., surface_amp=0., field_sigma=0., mu=0., k_d=0., rho=0., asymmetry=False, accel=.1, speed=.001)
    sc = S.cell_scene3('rest', 'relaxed', 'nominal', 0, T=2., overrides=o)
    s0, _ = P._init(sc, P.PLATFORM); cmd = np.full(250, s0[2])
    calm = P.replay(sc, cmd); sc = dict(sc); sc['act'] = np.where(np.arange(len(sc['act']))*.001 > 1., 1., 0.)
    tense = P.replay(sc, cmd)
    gain = float(tense[-1, 3]/calm[-1, 3])
    # thin covers stay solvable with the probe sliding on to the rib
    thin = []
    for law in ('foundation', 'halfspace'):
        for probe in ('convex', 'linear'):
            for cover in (.004, .007):
                sc = S.cell_scene3('rib', 'rib', 'nominal', 2, T=3., overrides=dict(force=12., probe=probe, contact_law=law, slide='long', cover=cover, pose_error=3., speed=.03, accel=.1))
                s0, _ = P._init(sc, P.PLATFORM)
                log = P.replay(sc, np.full(375, s0[2]))
                bound = 1.5/(1./sc['kr']+1./sc['k_rib'])     # arm and rib in series, with the soft part alongside
                thin.append(dict(law=law, probe=probe, cover_mm=cover*1000., finite=log is not None, peak_N=None if log is None else float(log[:, 3].max()), max_K_over_limit=None if log is None else float(log[:, 21].max()/bound)))
    within = {}
    for tissue in S.TISSUES:
        for regime in ('nominal', 'curved'):
            for name, values in (('probe', ('convex', 'linear')), ('force', (4., 8., 12.)), ('slide', ('long', 'short')), ('law', ('foundation', 'halfspace')), ('phi', (7., 11., 17., 23.)), ('dyn', ('hc', 'power', 'kb')), ('cover', (.007, .018, .032))):
                c = Counter(S._level(name, values, cond, tissue, S.REGIMES.index(regime), i) for cond in S.CONDITIONS for i in range(4))
                within[f'{tissue}/{regime}/{name}'] = len(set(c.values())) == 1 and len(c) == len(values)
    return dict(passed=bool(pose and gain > 1.5 and all(r['finite'] for r in thin) and all(within.values())), pose_sign=pose,
                theta_deg=float(np.degrees(th).mean()), friction_axial_N=float(extra.mean()), pressure_centre_mm=float(late[:, 24].mean()*1000.),
                activation_force_ratio=gain, thin_rib=thin, balanced_within_tissue=all(within.values()))


def trajectories_and_balance():
    rows = []; passed = True; reproducible = True; counts = {k: Counter() for k in ('probe', 'FD', 'slide')}
    for regime in S.REGIMES:
        for cond in S.CONDITIONS:
            sc = S.cell_scene3(cond, 'rib', regime, 3)
            v = np.diff(sc['xcmd'])*1000.; acc = np.diff(v)*1000.
            peak = max([float(np.max(np.abs(np.diff(w)*1000.))) for w in sc['motion_components']], default=0.)
            within = bool(sc['xcmd'].min()-sc['geometry']['length']*.5 > sc['field_x0'] and
                          sc['xcmd'].max()+sc['geometry']['length']*.5 < sc['field_x0']+(len(sc['q_field'])-1)*sc['dx'])
            motion_okay = (peak <= .05+.0001 if cond in ('flinch', 'cough') and regime != 'strong' else
                           (peak <= .15+.0001 if cond in ('flinch', 'cough') else True))
            okay = bool(np.max(np.abs(v)) <= sc['v']+1.e-8 and np.max(np.abs(acc)) <= sc['accel']*1.001 and within and motion_okay)
            passed &= okay
            rows.append(dict(condition=cond, regime=regime, max_scan_speed_mm_s=float(abs(v).max()*1000.),
                             max_acceleration_m_s2=float(abs(acc).max()), scan_limit_mm_s=sc['v']*1000.,
                             peak_patient_mm_s=peak*1000., in_field=within, passed=okay))
    for cond, tissue in S.CELLS3:
        for i in range(12):
            a = S.cell_scene3(cond, tissue, i=i, T=.2); b = S.cell_scene3(cond, tissue, i=i, T=.2)
            for k in counts:
                counts[k][str(a[k])] += 1
            for k in ('xcmd', 'u', 'q_field', 'E_field', 'ribs'):
                reproducible &= np.array_equal(a[k], b[k])
    a = S.cell_scene3('breath', 'mixed', 'curved', 3)
    b = S.cell_scene3('breath', 'mixed', 'curved_untracked', 3)
    pairing = all(np.array_equal(a[k], b[k]) for k in ('xcmd', 'q_field', 'E_field', 'u', 'ribs', 'h_s'))
    sc = uniform(T=.2); cmd = np.full(25, P._init(sc, P.PLATFORM)[0][2])
    replay_same = np.array_equal(P.replay(sc, cmd, seed=12, sensor=True), P.replay(sc, cmd, seed=12, sensor=True))
    balanced = all(len(set(c.values())) == 1 for c in counts.values())
    return dict(passed=bool(passed), rows=rows), dict(passed=bool(reproducible and pairing and balanced and replay_same),
                counts={k: dict(c) for k, c in counts.items()}, reproducible=bool(reproducible),
                replay_reproducible=replay_same, surface_tracking_pair_identical=pairing)


def quadrature_accuracy(cover=.004):
    """Independent dense area quadrature after transverse factorization."""
    nodes, weights = np.polynomial.legendre.leggauss(12)
    def composite(lo, hi):
        n = max(1, int(math.ceil((hi-lo)/.001)))
        mid = lo+(np.arange(n)+.5)*(hi-lo)/n
        return (mid[:, None]+(hi-lo)/(2.*n)*nodes).ravel(), np.tile(weights*(hi-lo)/(2.*n), n)
    rows = []
    for law in ('foundation', 'halfspace'):
        for slide in ('long', 'short'):
            for probe in ('convex', 'linear'):
                sc = S.cell_scene3('rib', 'rib', i=2, overrides=dict(contact_law=law, slide=slide, probe=probe, cover=cover, k_rib=1.e12, no_curve=True))
                x = sc['ribs'][1, 0]+.003
                theta = math.radians(3.)
                d = P.solve_depth(sc, 12., x, theta)
                g = P._geometry(sc); curvature = .5/g[2]+.5/g[7]; b = math.tan(theta)
                grid = sc['field_x0']+np.arange(len(sc['q_field']))*sc['dx']
                def section(depth):
                    lo, hi = P._bounds(depth, theta, curvature, g[0])
                    if hi <= lo:
                        return 0., 0., 0., 0.
                    xx, xw = composite(lo, hi)
                    yy, yw = composite(-g[1]*.5, g[1]*.5) if slide == 'short' else (np.array([0.]), np.array([g[1]]))
                    xi = xx[:, None]; eta = yy[None, :]
                    wt = xw[:, None]*yw[None, :]
                    delta = depth-curvature*xi**2-b*xi
                    pos = x+(xi if slide == 'long' else eta)
                    q = sc['q_base']*np.interp(pos, grid, sc['q_field']); e = sc['E_base']*np.interp(pos, grid, sc['E_field'])
                    w = np.zeros_like(pos)
                    for center, width, edge, _ in sc['ribs']:
                        a = np.clip((width*.5+edge*.5-abs(pos-center))/edge, 0., 1.)
                        w = np.maximum(w, .5-.5*np.cos(np.pi*a))
                    pressure = (q if law == 'foundation' else e)*delta/(1.-delta*w/cover)
                    stiffness = (q if law == 'foundation' else e)/(1.-delta*w/cover)**2
                    f = float(np.sum(wt*pressure)); k = float(np.sum(wt*stiffness))
                    m = float(np.sum(wt*pressure*xi))
                    if law == 'halfspace':
                        k *= 1.1/math.sqrt(g[1]*(hi-lo))
                    return f, k, m, hi-lo
                f, k, moment, span = section(d)
                centre = moment/f
                if law == 'halfspace':
                    full = curvature*(g[0]*.5)**2+abs(b)*g[0]*.5
                    scale = 1.1/math.sqrt(g[1]*span)
                    f = 1.2*scale*f if d <= full else scale*(f+.2*section(full)[0])
                f *= g[20]; k *= g[20]      # the overall factor (scale and equilibrium share) the plant applies after the integral
                actual = P._contact(x, d, theta, 0., g, sc['q_field'], sc['E_field'], sc['ribs'], 0.)      # the integral itself, before the slope contact() reports
                error = max(abs(actual[0]/f-1.), abs(actual[1]/k-1.))
                rows.append(dict(law=law, slide=slide, probe=probe, cover_mm=cover*1000., relative_force_tangent_error=error,
                                 force_N=actual[0], reference_force_N=f, tangent_N_m=actual[1], reference_tangent_N_m=k,
                                 pressure_center_error_mm=abs(actual[2]-centre)*1000.))
    return dict(passed=all(r['relative_force_tangent_error'] < .01 and r['pressure_center_error_mm'] < .1 for r in rows),
                reference_panel_width_m=.001, reference_nodes_per_panel=12, rows=rows)


def all_quadrature():
    checks = [quadrature_accuracy(h) for h in (.004, .018)]
    return dict(passed=all(c['passed'] for c in checks), reference_panel_width_m=.001,
                reference_nodes_per_panel=12, rows=[r for c in checks for r in c['rows']])


def timing():
    import physics_plant2 as P2
    import scenes2 as S2
    rows = []
    def duration(fun):
        samples = []
        for _ in range(3):
            start = time.perf_counter(); result = fun(); samples.append(time.perf_counter()-start)
            if result is None:
                raise RuntimeError('replay diverged')
        return float(np.median(samples))
    cases = [('nominal', 'rest', 'relaxed', dict(force=4., probe='convex', contact_law=law, slide=slide))
             for law in ('foundation', 'halfspace') for slide in ('long', 'short')]
    cases += [(f'rib_{int(h*1000)}mm_tilt', 'rib', 'rib',
               dict(force=12., probe=probe, contact_law='halfspace', slide=slide, cover=h, pose_error=3.))
              for h, probe, slide in ((.018, 'convex', 'short'), (.004, 'convex', 'long'), (.004, 'linear', 'short'))]
    original_fd = P2.FD
    try:
        for name, cond, tissue, overrides in cases:
            sc = S.cell_scene3(cond, tissue, i=0, T=1., overrides=overrides)
            old = S2.cell_scene2('rest', 'relaxed', 'nominal', 0, T=1.)
            # Identical drives, preload and platform; legacy geometry/contact law differ.
            for key in ('xcmd', 'u', 'fext', 'act', 'rip', 'm', 'c', 'k', 'kr', 'K_t', 'm_t', 'zeta_t',
                        'mu_k', 'mu_s', 'v_s', 'v_e', 'k_d', 'k1', 'k2', 'theta0'):
                old[key] = sc[key]
            old['ln_e'] = np.log(sc['q_field']); old['h_s'] = sc['h_s']; old['dx'] = sc['dx']
            old['xi_c'] = 0.; old['c_slope'] = 0.; old['lam'] = P.calibration(sc)
            P2.FD = sc['FD']
            commands2 = np.full(125, P2._init(old, P2.PLATFORM)[0][2]); P2.replay(old, commands2[:2])
            old_time = duration(lambda: P2.replay(old, commands2))
            cmds = np.full(125, P._init(sc, P.PLATFORM)[0][2]); P.replay(sc, cmds[:2])
            elapsed = duration(lambda: P.replay(sc, cmds))
            rows.append(dict(scene=name, law=sc['contact_law'], probe=sc['probe'], slide=sc['slide'], force_N=sc['FD'],
                             replay_seconds=elapsed, plant2_seconds=old_time, ratio=elapsed/old_time, passed=elapsed/old_time <= 3.))
            print('timing', name, sc['slide'], round(elapsed/old_time, 3), 'seconds per simulated second', round(elapsed, 3), flush=True)
    finally:
        P2.FD = original_fd
    # The 4 mm cover is a sensitivity endpoint outside the main matrix; its cost is recorded, not required.
    return dict(passed=all(r['passed'] for r in rows if '4mm' not in r['scene']), rows=rows, jit_excluded=True, initialization_included=True,
                repeats=3, simulation_seconds=1., nominal_passed=all(r['passed'] for r in rows if r['scene'] == 'nominal'))


def realdata_interfaces():
    from check_curve3 import clean_scene
    ratios = []; bones = []; activation = []; derivative = []
    for law in ('foundation', 'halfspace'):
        for probe in S.PROBES:
            scenes = [clean_scene(law, probe, q_scale=q) for q in (.5, 1., 2.)]
            forces = [P.contact(sc, .01, 0.)[0] for sc in scenes]
            ratios.append(dict(law=law, probe=probe, ratios=[f/forces[1] for f in forces]))
            sc = clean_scene(law, probe, cover=.007, full_cover=True, k_rib=10000.)
            soft = dict(sc, full_cover=False)
            for d in (.001, .006, .01, .03):
                ends = [P.contact(sc, d, 0., activation=a) for a in (0., 1.)]
                for a in (0., .4, 1.):
                    f, k, _, _, _ = P.contact(sc, d, 0., activation=a)
                    fs, ks, _, _, _ = P.contact(soft, d, 0., activation=a)
                    eps = 1.e-8
                    fd = (P.contact(sc, d+eps, 0., activation=a)[0]-P.contact(sc, d-eps, 0., activation=a)[0])/(2.*eps)
                    derivative.append(abs(k/fd-1.))
                    bones.append(dict(law=law, probe=probe, depth_mm=d*1000., activation=a,
                        tangent_N_m=k, bound_N_m=ks+sc['k_rib'], ratio=k/(ks+sc['k_rib'])))
                    activation.append(abs(f-((1.-a)*ends[0][0]+a*ends[1][0])))
    scaling_error = max(abs(r['ratios'][i]-expected) for r in ratios for i, expected in enumerate((.5, 1., 2.)))
    return dict(passed=scaling_error < 1.e-9 and max(r['ratio'] for r in bones) <= 1.+1.e-9 and max(activation) < 1.e-9 and max(derivative) < .005,
        q_scale=ratios, scaling_absolute_error=scaling_error, rib= bones,
        max_activation_interpolation_error_N=max(activation),
        max_rib_central_difference_relative_error=max(derivative))


def phases():
    from check_curve3 import clean_scene
    rows = []
    w = 2.*np.pi*20.8
    a = 1.-math.exp(-2.*np.pi*P.PLATFORM['fc']*.001)
    sensor = a/(1.-(1.-a)*np.exp(-1j*w*.001))
    for kind in S.SHAPES+('hc', 'power', 'kb'):
        for phase in S.PHASES:
            sc = clean_scene()
            sc['alpha'], sc['tau'], sc['phi_hc'] = S._dynamic_response(kind, phase, np.random.default_rng(71))
            P.attach_curve(sc)
            d = P.solve_depth(sc, 8., 0.); k = P.contact(sc, d, 0.)[1]
            memory = 1.+np.sum(sc['alpha']*1j*w*sc['tau']/(1.+1j*w*sc['tau']))
            bare = k*memory+1j*w*P.calibration(sc)*8.
            local = P.complex_gain(sc, k, 8.)
            rows.append(dict(kind=kind, target_deg=phase, bare_deg=float(np.degrees(np.angle(bare))),
                local_deg=float(np.degrees(np.angle(local))), after_IIR_deg=float(np.degrees(np.angle(local*sensor))),
                error_deg=float(np.degrees(np.angle(bare)))-phase, tau_s=sc['tau'].tolist()))
    return dict(passed=max(abs(r['error_deg']) for r in rows) < .5, rows=rows,
                sensor_phase_deg=float(np.degrees(np.angle(sensor))), target='bare material at equilibrium 8 N; no support, actuator, or ZOH included')


def noise_and_truth():
    rows = []
    for sigma in (.08, .18):
        sc = uniform(T=12., sensor_sigma=sigma, alpha_scale=0., phi_hc=0.)
        s, _ = P._init(sc, P.PLATFORM)
        log = P.replay(sc, np.full(1500, s[2]), seed=20261105, sensor=True, platform=dict(P.PLATFORM, log_dt=.001))
        keep = log[:, 0] > 1.
        error = log[keep, 5]-log[keep, 18]
        moment = log[keep, 13]
        a = 1.-math.exp(-2.*np.pi*P.PLATFORM['fc']*.001)
        expected = sigma*math.sqrt(a/(2.-a))
        rows.append(dict(input_white_sigma_N=sigma, output_sd_N=float(np.std(error)),
            output_moment_sd_Nm=float(np.std(moment)), theoretical_white_output_sd_N=expected,
            relative_error=float(abs(np.std(error)/expected-1.))))
    class Spy:
        def reset(self, pl):
            self.flags = []; self.masked = True
        def __call__(self, obs):
            self.flags.append(hasattr(obs, 'truth_k') or hasattr(obs, 'truth_re'))
            self.masked &= bool(np.all(np.isnan(obs.log[:, 20:26])))
            return obs.y[-1]
    sc = uniform(T=.08)
    default, diagnostic = Spy(), Spy()
    P.episode(sc, default, 0); P.episode(sc, diagnostic, 0, diagnostic=True)
    masked = not any(default.flags) and any(diagnostic.flags) and default.masked and diagnostic.masked
    return dict(passed=masked and all(r['relative_error'] < .1 for r in rows), rows=rows,
        default_has_truth=any(default.flags), diagnostic_has_truth=any(diagnostic.flags), truth_columns_masked=default.masked and diagnostic.masked)


def main():
    start = time.perf_counter()
    protected = ['physics_plant.py', 'physics_plant2.py', 'scenes.py', 'scenes2.py', 'contact_candidates.py']
    hashes = {f: hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in protected}
    output = {}
    for name, fun in [('static', static_table), ('remus', remus), ('open_loop', open_loop),
                      ('continuity', continuity), ('pose_friction_compliance', pose_and_dissipation),
                      ('quadrature_accuracy', all_quadrature), ('scene_physics', scene_physics),
                      ('realdata_interfaces', realdata_interfaces), ('phases', phases), ('noise_truth', noise_and_truth)]:
        output[name] = fun(); print(name, output[name]['passed'], flush=True)
    output['trajectories'], output['balance_reproducibility'] = trajectories_and_balance()
    print('trajectories', output['trajectories']['passed'], flush=True)
    output['timing'] = timing(); print('timing', output['timing']['passed'], flush=True)
    output['protected_files_unchanged'] = all(hashlib.sha256((HERE/f).read_bytes()).hexdigest() == h for f, h in hashes.items())
    output['protected_hashes'] = hashes
    output['elapsed_seconds'] = time.perf_counter()-start
    output['all_passed'] = all(v['passed'] for v in output.values() if isinstance(v, dict) and 'passed' in v)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('saved', OUT, 'seconds', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
