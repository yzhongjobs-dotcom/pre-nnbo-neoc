"""Static scene proposal audit; no plant imports, simulation or model selection.

Units inside contact functions: mm, N, N/mm^3. This independent foundation
approximation transfers indenter force/volume to finite probe geometry; it is
not a Hertz solution or an identified human constitutive law.
"""
import cmath
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT = HERE / 'results' / 'scene_settings_check.json'
FREQUENCY = 20.8
FORCES = [4., 8., 12.]
PROBES = {
    'convex_idealized': dict(length_mm=69., width_mm=17., radius_mm=60.),
    'linear_L4_20t_RS': dict(length_mm=48.43, width_mm=11.56, radius_mm=None),
}
TISSUES = {
    'relaxed_lateral': dict(q_N_mm3=.0011, cover_mm=None),
    'relaxed_anterior': dict(q_N_mm3=.0016, cover_mm=None),
    'active_equivalent': dict(q_N_mm3=.004, cover_mm=None),
    'rib_cover_18mm': dict(q_N_mm3=.0011, cover_mm=18.),
    'rib_cover_7mm': dict(q_N_mm3=.0011, cover_mm=7.),
}
ALPHA = [.35, .35, .35]
TAU = [.008, .08, .8]


def hemisphere_volume(d, r=10.):
    if d <= r:
        return math.pi * (r*d*d-d**3/3.)
    return 2.*math.pi*r**3/3. + math.pi*r*r*(d-r)


def contact(d, probe, tissue):
    if d <= 0.:
        return 0., 0., 0.
    length, width, radius = (probe[k] for k in ('length_mm', 'width_mm', 'radius_mm'))
    q, h = tissue['q_N_mm3'], tissue['cover_mm']
    if radius is None:
        area = length*width
        force, stiffness, contact_width = q*area*d, q*area, length
    else:
        # Parabolic face along the long axis; full short-axis contact assumed.
        a = min(math.sqrt(2.*radius*d), length/2.)
        force = q*width*(2.*a*d-a**3/(3.*radius))
        stiffness, contact_width = q*width*2.*a, 2.*a
    if h is not None:
        # Positive, monotone finite-cover proxy; diverges at the bone boundary.
        factor = 1./(1.-d/h)
        stiffness = stiffness*factor + force/h*factor**2
        force *= factor
    return force, stiffness, contact_width


def solve(force, probe, tissue):
    lo, hi = 0., tissue['cover_mm'] or 1000.
    hi *= 1.-1e-10
    for _ in range(100):
        mid = (lo+hi)/2.
        if contact(mid, probe, tissue)[0] < force:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2.


def main():
    protected = ['physics_plant2.py', 'scenes.py', 'scenes2.py', 'contact_candidates.py']
    hashes = {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in protected}
    # Published Table 5 points; no values extracted from plotted pixels.
    points = [('R1', 10., 3.33), ('R1', 15., 5.65),
              ('R3', 12.5, 2.92), ('R3', 18.75, 5.39)]
    calibration = []
    for region, d, force in points:
        q = .0016 if region == 'R1' else .0011
        predicted = q*hemisphere_volume(d)
        calibration.append(dict(region=region, depth_mm=d, measured_force_N=force,
                                inferred_q_N_mm3=force/hemisphere_volume(d),
                                proposed_q_N_mm3=q, predicted_force_N=predicted,
                                relative_error=predicted/force-1.))
    holdout = []
    for region, d, force in [('R1', 5., 1.45), ('R1', 20., 8.40),
                              ('R3', 6.25, 1.12), ('R3', 25., 8.53)]:
        q = .0016 if region == 'R1' else .0011
        predicted = q*hemisphere_volume(d)
        holdout.append(dict(region=region, depth_mm=d, measured_force_N=force,
                            predicted_force_N=predicted, relative_error=predicted/force-1.))
    omega = 2.*math.pi*FREQUENCY
    memory = 1.+sum(a*(1j*omega*t)/(1.+1j*omega*t) for a, t in zip(ALPHA, TAU))
    rows = []
    checks = dict(force_balance=True, footprint_bounded=True, positive_stiffness=True,
                  monotone_force=True, monotone_q_stiffness_at_fixed_force=True,
                  no_force_saturation=True, continuous_at_full_footprint=True,
                  derivative_consistent=True, bone_not_crossed=True,
                  remus_calibration_within_20_percent=True)
    for probe_name, probe in PROBES.items():
        for tissue_name, tissue in TISSUES.items():
            ref_d = solve(8., probe, tissue)
            ref_k = contact(ref_d, probe, tissue)[1]*1000.
            # Define HC strength by its separate phase at 8 N, not scan speed.
            lam = math.tan(math.radians(5.))*ref_k/(omega*8.)
            for force in FORCES:
                d = solve(force, probe, tissue)
                actual, k_mm, width = contact(d, probe, tissue)
                k = k_mm*1000.
                gain = k*memory+1j*omega*lam*force
                robot_gain = gain/(1.+gain/1e5)
                rows.append(dict(probe=probe_name, tissue=tissue_name, force_N=force,
                                 depth_mm=d, contact_long_mm=width,
                                 contact_short_mm=probe['width_mm'], stiffness_N_mm=k_mm,
                                 gain20_8_N_mm=abs(gain)/1000., phase20_8_deg=math.degrees(cmath.phase(gain)),
                                 robot_series_gain20_8_N_mm=abs(robot_gain)/1000.,
                                 robot_series_phase20_8_deg=math.degrees(cmath.phase(robot_gain)),
                                 lambda_s_m=lam, full_footprint=abs(width-probe['length_mm'])<1e-8))
                checks['force_balance'] &= abs(actual-force)<1e-8
                checks['footprint_bounded'] &= width <= probe['length_mm']+1e-8
                checks['positive_stiffness'] &= k_mm>0.
                checks['bone_not_crossed'] &= tissue['cover_mm'] is None or d<tissue['cover_mm']
                eps = 1e-5
                numeric = (contact(d+eps, probe, tissue)[0]-contact(d-eps, probe, tissue)[0])/(2.*eps)
                checks['derivative_consistent'] &= abs(numeric/k_mm-1.)<1e-5
            max_d = min(40., .95*(tissue['cover_mm'] or 1000.))
            fs = [contact(max_d*i/100., probe, tissue)[0] for i in range(101)]
            checks['monotone_force'] &= all(b>a for a, b in zip(fs, fs[1:]))
            # The inverse must remain finite above all proposed load levels.
            d24, d48 = solve(24., probe, tissue), solve(48., probe, tissue)
            checks['no_force_saturation'] &= d48>d24 and abs(contact(d48, probe, tissue)[0]-48.)<1e-7
        if probe['radius_mm'] is not None:
            edge = (probe['length_mm']/2.)**2/(2.*probe['radius_mm'])
            t = TISSUES['relaxed_lateral']
            below, above = contact(edge-1e-7, probe, t), contact(edge+1e-7, probe, t)
            checks['continuous_at_full_footprint'] &= abs(above[0]-below[0])<1e-6 and abs(above[1]-below[1])<1e-6
        for force in FORCES:
            ks = []
            for q in [.001, .0011, .0016, .004, .005]:
                t = dict(q_N_mm3=q, cover_mm=None)
                ks.append(contact(solve(force, probe, t), probe, t)[1])
            checks['monotone_q_stiffness_at_fixed_force'] &= all(b>a for a, b in zip(ks, ks[1:]))
    checks['remus_calibration_within_20_percent'] = all(abs(p['relative_error'])<.2 for p in calibration)
    # Analytic scene-drive bounds: these are prescribed support motion, not
    # the actual body response after mass/damping/contact force are included.
    n = 3
    cosine_factor = 2.*n*(1.-1./(2.*n))**((2.*n-1.)/2.)/math.sqrt(2.*n)
    drives = dict(old_breath_peak_speed_upper_mm_s=cosine_factor*math.pi*8.*1.2/(3.*.85),
                  proposed_breath_sinusoid_peak_speed_mm_s=math.pi*8./4.2,
                  old_slow_minjerk_speed_mm_s=[1.875*2./1., 1.875*15./.5],
                  old_fast_minjerk_speed_mm_s=[1.875*2./.4, 1.875*12./.15],
                  old_flinch_initial_speed_envelope_mm_s=4.*5./.03,
                  old_cough_single_pulse_initial_speed_envelope_mm_s=4.*3./.02,
                  note='Conservative twitch bounds use rise<=decay: normalized peak >= 1/4; not measured speeds.')
    result = dict(model='Finite foundation proposal, not existing plant or human validation',
                  frequency_Hz=FREQUENCY, forces_N=FORCES, probes=PROBES, tissues=TISSUES,
                  remus_source='https://doi.org/10.3389/fbioe.2024.1384062 Table 5',
                  calibration=calibration, same_study_unused_points=holdout,
                  unused_point_max_abs_relative_error=max(abs(p['relative_error']) for p in holdout),
                  alpha=ALPHA, tau_s=TAU,
                  hc_reference_phase_deg=5., hc_reference_force_N=8.,
                  dynamic_definition='dF/dd = k * (1 + sum(alpha*s*tau/(1+s*tau))) + s*lambda*F; bare contact, positive force lead',
                  memory_gain_ratio=abs(memory), memory_phase_deg=math.degrees(cmath.phase(memory)),
                  rows=rows, checks=checks, drive_bounds=drives, protected_file_sha256=hashes)
    assert all(checks.values()), checks
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(rows=len(rows), checks=checks, calibration=calibration,
                          depth_range_mm=[min(r['depth_mm'] for r in rows), max(r['depth_mm'] for r in rows)],
                          memory_gain_ratio=abs(memory), memory_phase_deg=result['memory_phase_deg'],
                          drive_bounds=drives), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
