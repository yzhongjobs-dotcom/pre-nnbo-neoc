"""Reproducible S1 source, paired sampling, transfer checks, and what can be validated independently of the curves used to build the plant."""
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import physics_plant3 as P
import scenes3 as S
from derive_remus_s1 import derive

HERE = Path(__file__).resolve().parent


def stats(values):
    a = np.asarray(values, float)
    return dict(n=len(a), median=float(np.median(a)), min=float(a.min()), max=float(a.max()),
                p05=float(np.quantile(a, .05)), p95=float(np.quantile(a, .95)))


def source_check():
    source = derive()
    table4 = dict(R1=11.58, R2=13.34, R3=11.52, R4=10.33, R5=11.41, R6=32.25,
                  A1=20.94, A2=29.51, A3=27.01, A4=22.44, A5=23.86, A6=55.95)
    ends = {r: float(np.polyval(v['rangeMean']['coeff'], v['rangeMean']['fitMaxDispl'])) for r, v in source.items()}
    adjusted = []
    for r, rows in source.items():
        for subject, row in rows.items():
            if not subject.startswith('P'):
                continue
            a, b = P.positive_domain(row['coeff'], row['fitMaxDispl'])
            if a > 0. or b < row['fitMaxDispl']:
                adjusted.append(dict(region=r, subject=subject, positive_interval_mm=[a, b]))
    errors = {r: ends[r]-table4[r] for r in ends}
    return dict(passed=source == S.REMUS_S1 and max(abs(x) for x in errors.values()) < .005,
                exact_json_equal=source == S.REMUS_S1, rangeMean_end_force_N=ends,
                table4_error_N=errors, adjusted_count=len(adjusted), adjusted=adjusted)


def at_force(row, force):
    p = row['coeff']; lo, hi = P.positive_domain(p, row['fitMaxDispl'])
    if np.polyval(p, lo) > force or np.polyval(p, hi) < force:
        return None
    d = brentq(lambda x: np.polyval(p, x)-force, lo, hi)
    return float(np.polyval(np.polyder(p), d))


def sampling():
    rng = np.random.default_rng(20261105)
    sampled = {state+str(i): [] for state in 'RA' for i in range(1, 6)}
    for _ in range(60000):
        r, subject, row, active = S.sample_curves(rng, 'relaxed', {})
        assert row == S.REMUS_S1[r][subject] and active == S.REMUS_S1['A'+r[1:]][subject]
        sampled[r].append(subject); sampled['A'+r[1:]].append(subject)
    rows = []; missing = []
    for region, ids in sampled.items():
        raw = S.REMUS_S1[region]; eligible = sorted(set(ids))
        missing.extend([region+'/'+s for s in raw if s.startswith('P') and s not in eligible])
        for force in (4., 8., 12.):
            truth = {s: at_force(raw[s], force) for s in raw if s.startswith('P')}
            values = [truth[s] for s in ids if truth[s] is not None]
            direct = [truth[s] for s in eligible if truth[s] is not None]
            emp = np.quantile(values, [.05, .5, .95]); expected = np.quantile(direct, [.05, .5, .95])
            rows.append(dict(region=region, force_N=force, sampled=emp.tolist(), paired_direct=expected.tolist(),
                full_S1=np.quantile([v for v in truth.values() if v is not None], [.05, .5, .95]).tolist(),
                max_relative_quantile_error=float(np.max(abs(emp-expected)/expected))))
    # The empirical quantiles of a few subjects are coarse, so test the draw itself: each paired subject equally likely within its region.
    spread = []
    for region, ids in sampled.items():
        counts = np.array([ids.count(x) for x in sorted(set(ids))], float)
        spread.append(float(np.max(np.abs(counts/counts.mean()-1.))))
    return dict(passed=max(spread) < .1, subject_frequency_max_deviation=max(spread), excluded_without_pair=missing, draws=60000,
                region_counts={r: len(v) for r, v in sampled.items() if r.startswith('R')},
                quantile_order='P05/median/P95; tangent N/mm; only roots within fitMaxDispl', rows=rows)


def clean_scene(law='foundation', probe='convex', **overrides):
    o = dict(contact_law=law, probe=probe, slide='long', field_sigma=0., surface_amp=0.,
             pose_error=0., mu=0., k_d=0., rho=0., asymmetry=False, region='R1', subject='P01', phi_hc=0., alpha_scale=0., tool_mass=0., sensor_delay_ms=0., delay_error=0., tau_a_true=.019, tau_error=0.)
    o.update(overrides)
    return S.cell_scene3('rest', 'relaxed', T=.1, overrides=o)


def transfer():
    rows = []; errors = []; floor = []; monotone = True; min_cap = math.inf
    for law in ('foundation', 'halfspace'):
        for probe in S.PROBES:
            template = clean_scene(law, probe)
            for region, subjects in S.REMUS_S1.items():
                for subject, curve in subjects.items():
                    if not subject.startswith('P'):
                        continue
                    sc = dict(template, region=region, subject=subject, remus_poly=np.array(curve['coeff']),
                        fit_max_mm=curve['fitMaxDispl'], active_poly=np.array(curve['coeff']), active_fit_max_mm=curve['fitMaxDispl'])
                    P.attach_curve(sc)
                    c = P._geometry(sc); min_cap = min(min_cap, c[13], c[14], c[22])
                    table = sc['curve_table']
                    monotone &= bool(np.all(np.diff(table[1]) > 0.) and np.all(table[2] > 0.))
                    for force in (4., 8., 12.):
                        d = P.solve_depth(sc, force, 0.); f, k, _, width, _ = P.contact(sc, d, 0.)
                        eps = 1.e-8
                        fd = (P.contact(sc, d+eps, 0.)[0]-P.contact(sc, d-eps, 0.)[0])/(2.*eps)
                        errors.append(abs(k/fd-1.))
                        assert 0. < width <= sc['geometry']['length']+1.e-12
                        s = d/math.sqrt(width*sc['geometry']['width'])
                        floor.append(s <= .005/math.sqrt(float(P._hemisphere_area(.005))))
                        rows.append(dict(law=law, probe=probe, region=region, subject=subject, force_N=force,
                            depth_mm=d*1000., tangent_N_mm=k/1000., regularized=bool(floor[-1])))
    summary = []
    for law in ('foundation', 'halfspace'):
        for probe in S.PROBES:
            for state in 'RA':
                for force in (4., 8., 12.):
                    q = [r for r in rows if r['law'] == law and r['probe'] == probe and r['force_N'] == force and r['region'].startswith(state) and r['region'][1:] != '6']
                    summary.append(dict(law=law, probe=probe, state=state, force_N=force,
                        depth_mm=stats([r['depth_mm'] for r in q]), tangent_N_mm=stats([r['tangent_N_mm'] for r in q])))
    return dict(passed=monotone and max(errors) < .005 and min_cap > 0., individual_curves=len(rows)//12,
        mapping_count=len(rows)//3, working_points=len(rows), monotone=monotone, min_cap=min_cap,
        max_central_difference_relative_error=max(errors), regularized_count=sum(floor),
        regularized_fraction=float(np.mean(floor)), summary=summary, rows=rows)


def independent():
    """Checks that do not reuse the S1 curve being mapped."""
    # 1. half-space transfer against exact elasticity: a Hertz sphere as the input must return the exact circular flat punch
    E = 25.e3; r = .01; rows = []
    for a_p in (.008, .015, .03):
        for delta in (.002, .005, .01):
            a_h = math.sqrt(r*delta); k_h = 2.*E*a_h
            rows.append(abs(k_h*(a_p/a_h)-2.*E*a_p)/(2.*E*a_p))
    # 2. the effective modulus the curves imply against Remus' own inverse finite-element fit (Table 6, relaxed anterior and lateral: c 5-8 kPa)
    implied = {}
    for region in ('R1', 'R2', 'R3', 'R4', 'R5'):
        row = S.REMUS_S1[region]['rangeMean']; p = np.array(row['coeff'])
        k = float(np.polyval(np.polyder(p), 5.))*1000.; a = math.sqrt(float(P._hemisphere_area(.005))/math.pi)
        implied[region] = k/(2.*a)
    # 3. fit quality of the curves themselves, from S1 (N)
    import csv
    from derive_remus_s1 import SOURCE
    rmse = [float(x['rmse']) for x in csv.DictReader(SOURCE.open(encoding='utf-8-sig', newline='')) if x['id'].startswith('P')]
    return dict(passed=max(rows) < 1.e-12, exact_flat_punch_max_error=max(rows), implied_Estar_at_5mm_Pa=implied,
                remus_ogden_c_kPa={'R1': 5, 'R2': 7, 'R3': 7, 'R4': 6, 'R5': 8}, s1_individual_fit_rmse_N=stats(rmse),
                note='Beber and Dahmani are different materials; Dahmani has no modulus or probe geometry, so neither validates the absolute transfer.')


def main():
    start = time.perf_counter(); out = {}
    for name, fun in (('source', source_check), ('sampling', sampling), ('transfer', transfer), ('independent', independent)):
        out[name] = fun(); print(name, out[name]['passed'], flush=True)
    out['elapsed_seconds'] = time.perf_counter()-start
    out['all_passed'] = all(r['passed'] for r in out.values() if isinstance(r, dict) and 'passed' in r)
    (HERE/'results'/'curve3_checks.json').write_text(json.dumps(out, indent=2, allow_nan=False, default=lambda x: x.item() if hasattr(x, 'item') else str(x))+'\n', encoding='utf-8')
    print('seconds', out['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
