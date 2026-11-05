"""Probe inertia and calibration mismatch: does the measurement model hold, does the compensation remove the term, and what does the Method lose to what it cannot know?"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
W = 2.*math.pi*20.8
QUIET = dict(contact_law='halfspace', probe='linear', slide='long', force=8., field_sigma=0., surface_amp=0., pose_error=0., mu=0., k_d=0., rho=0.,
             asymmetry=False, region='R3', subject='P04', phi_hc=0., alpha_scale=0., accel=.1, speed=.0001)
CASES = [(c, t, 'nominal', i) for c, t in (('breath', 'relaxed'), ('rest', 'relaxed'), ('breath', 'active'), ('shift_fast', 'relaxed'), ('breath', 'mixed'), ('rib', 'relaxed')) for i in range(4)]
CHANGES = {
    'ideal system, no probe mass': dict(tool_mass=0., mass_error=0., tau_error=0., tau_a_true=.019, sensor_delay_ms=0., delay_error=0.),
    'all system factors': dict(),
    'probe mass only, exact calibration': dict(mass_error=0., tau_error=0., sensor_delay_ms=0., delay_error=0.),
    'mass error only 5%': dict(mass_error=.05, tau_error=0., sensor_delay_ms=0., delay_error=0.),
    'mass error only -5%': dict(mass_error=-.05, tau_error=0., sensor_delay_ms=0., delay_error=0.),
    'servo lag error only 5%': dict(mass_error=0., tau_error=.05, sensor_delay_ms=0., delay_error=0.),
    'servo lag error only -5%': dict(mass_error=0., tau_error=-.05, sensor_delay_ms=0., delay_error=0.),
    'servo lag error only 10% (stress)': dict(mass_error=0., tau_error=.1, sensor_delay_ms=0., delay_error=0.),
    'force offset 2 ms, calibrated exactly': dict(mass_error=0., tau_error=0., sensor_delay_ms=2., delay_error=0.),
    'force offset 1 ms, calibration +.25 ms': dict(mass_error=0., tau_error=0., sensor_delay_ms=1., delay_error=.25),
    'force offset 1 ms, calibration -.25 ms': dict(mass_error=0., tau_error=0., sensor_delay_ms=1., delay_error=-.25),
    'force offset 1 ms, not calibrated': dict(mass_error=0., tau_error=0., sensor_delay_ms=1., delay_error=-1.),
    'force offset 2 ms, no probe mass, calibrated': dict(tool_mass=0., mass_error=0., tau_error=0., sensor_delay_ms=2., delay_error=0.),
    'no inertia compensation': dict(uncompensated=True, mass_error=0., tau_error=0., sensor_delay_ms=0., delay_error=0.),
}


def amplitude(t, x):
    a = np.c_[np.sin(W*t), np.cos(W*t), np.ones_like(t), t]
    c = np.linalg.lstsq(a, x, rcond=None)[0]
    return complex(c[0], c[1])


def measurement():
    """Open loop, sine command, no noise: reading minus contact force against m*a, then what the compensator leaves."""
    rows = []
    for mass in (.15, .3, .45):
        sc = S3.cell_scene3('rest', 'relaxed', 'nominal', 0, T=4., base=491000, overrides=dict(QUIET, tool_mass=mass, mass_error=0., tau_error=0., tau_a_true=.019, sensor_delay_ms=0., delay_error=0.))
        s0, _ = P3._init(sc, P3.PLATFORM)
        tick = np.arange(int(sc['T']/.008))*.008
        cmd = s0[2]+2.e-4*np.sin(W*tick)
        log = P3.replay(sc, cmd)
        ref = dict(sc, tool_mass=0.)
        base = P3.replay(ref, cmd)
        m = log[:, 0] > 1.
        extra = amplitude(log[m, 0], log[m, 5]-base[m, 5])
        y = amplitude(log[m, 0], log[m, 2])
        # the inertial force through the 30 Hz filter
        a = 1.-math.exp(-2.*math.pi*P3.PLATFORM['fc']*.001); filt = a/(1.-(1.-a)*np.exp(-1j*W*.001))*np.sinc(20.8*.008)     # filter and averaging over one period
        # the compensator run on the same log
        class Stub:
            step = staticmethod(lambda o: 0.); reset = staticmethod(lambda pl: None)
            ndel = P3.PLATFORM['ndel']; ns = 8; afc = a; tau = .019; eb = math.exp(-1.e-3/.019)
        comp = M.ForceSensing(Stub, mass, 0.); comp.reset(dict(P3.PLATFORM, tau_a=.019))
        class O:
            pass
        o = O(); o.n = len(log); o.log = log.copy()
        comp(o)
        left = amplitude(log[m, 0], o.log[m, 5]-base[m, 5])
        rows.append(dict(mass_kg=mass, y_amp_um=abs(y)*1e6, inertial_N=abs(extra), expected_mw2y_N=mass*W*W*abs(y), ratio=abs(extra)/(mass*W*W*abs(y)*abs(filt)),
                         after_compensation_N=abs(left), remaining_fraction=abs(left)/abs(extra),
                         apparent_stiffness_N_mm=mass*W*W/1e3, contact_change_um=float(np.max(np.abs(log[m, 6]-base[m, 6])))*1e6))
    return rows


def job(a):
    case, change = a
    c, t, g, i = case; spec = dict(CHANGES[change])
    unc = spec.pop('uncompensated', False)
    sc = S3.cell_scene3(c, t, g, i, base=492000, overrides=dict(force=8., dynamics='power', phase=11., **spec))
    pol = M.build('L', sc)
    if unc:
        pol.step = pol.step.step; pol.reset = type(pol).reset.__get__(pol)
    log = P3.episode(sc, pol, 950+i)
    m = M.metric(log, pol, sc); m.update(case=case, change=change)
    return m


def main():
    out = dict(measurement=measurement())
    print('open loop: probe mass, y amplitude um, inertial reading N, against m w2 y through the filter, left after compensation N, fraction, apparent stiffness N/mm, change in contact depth um')
    for r in out['measurement']:
        print(' ', r['mass_kg'], round(r['y_amp_um'], 1), round(r['inertial_N'], 3), round(r['ratio'], 3), round(r['after_compensation_N'], 4), round(r['remaining_fraction'], 3), round(r['apparent_stiffness_N_mm'], 2), round(r['contact_change_um'], 2))
    jobs = [(c, n) for c in CASES for n in CHANGES]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(job, jobs, chunksize=1))
    out['rows'] = rows
    base = {tuple(r['case']): r for r in rows if r['change'] == 'ideal system, no probe mass'}
    print('\nchange | axial RMSE median/mean | peak>15 | peak>36.5 | numeric | accept | gain error | L against ideal')
    for n in CHANGES:
        q = [r for r in rows if r['change'] == n]; ok = [r for r in q if not r['failed']]
        d = np.mean([np.log1p(r['J'])-np.log1p(base[tuple(r['case'])]['J']) for r in q])
        print(f"{n:38s} {np.median([r['rmse'] for r in ok]):.3f}/{np.mean([r['rmse'] for r in ok]):.3f}  {sum(r['peak']>15 for r in ok):2d} {sum(r['peak']>36.5 for r in ok):2d} {sum(r['failed'] for r in q):2d}  {np.mean([r['accept'] for r in ok]):.2f}  {np.nanmedian([abs(r['est_abs']) for r in ok]):.3f}  {d:+.3f}")
    (HERE/'results'/'system3_checks.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
