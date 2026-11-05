"""Read-only, paired diagnosis of D4 stiffness-window bias (stdout only).

Run: python -X utf8 -B d_bias_source.py [--workers 4]
Eight prescribed conditions, six scene seeds each; sensor seeds follow d_compare.
No result files, policy changes, or source patches. Extra diagnostics replay the
actual two-pole encoder motion, without rerunning or changing the closed loop.
"""
import os
import sys

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
# Keep JIT cache writes outside the project, including in spawned workers.
import tempfile
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(tempfile.gettempdir(), 'd_bias_source_numba'))

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import time

import numpy as np

import d_compare as D
import method3 as M
import physics_plant3 as P3
import scenes3 as S3
import common

HERE = Path(__file__).resolve().parent
BASE = dict(slide='long', force=8., dynamics='neural_white', phase=14.)
ELASTIC = dict(phi_hc=0., alpha_scale=0.)
CASES = (
    ('1 default H', {}, 'H'),
    ('2 elastic H', ELASTIC, 'H'),
    ('3 phase28 H', dict(phase=28.), 'H'),
    ('4 viscous H', dict(dynamics='viscous'), 'H'),
    ('5 rigid arm H', dict(kr=1.e9), 'H'),
    ('6 exact system H', dict(tool_mass=0., sensor_delay_ms=0., mass_error=0., tau_error=0., delay_error=0.), 'H'),
    ('7 noiseless H', dict(noise_scale=0.), 'H'),
    ('8 elastic L', ELASTIC, 'L'),
)


def quantiles(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if not len(a):
        return dict(n=0, median=None, q25=None, q75=None, iqr=None)
    q = np.percentile(a, [25, 50, 75])
    return dict(n=len(a), median=float(q[1]), q25=float(q[0]), q75=float(q[2]), iqr=float(q[2]-q[0]))


def encoder_replay(log, sc):
    """Identical 0.5 ms actuator recurrence; 1 ms force-sample positions per tick."""
    dt = .0005
    ns = int(round(P3.PLATFORM['Tc']/.001))
    t2 = P3.ACTUATOR_T2
    t1 = sc['tau_a_true']-t2
    ndel = P3.PLATFORM['ndel']
    y = float(log[0, 2]); velocity = 0.
    positions = np.empty((len(log), ns))
    inertia = np.empty(len(log)); acceleration_history = np.zeros(4)
    inertia_filter = 0.; afc = 1.-math.exp(-2.*math.pi*P3.PLATFORM['fc']*.001)
    lag = int(sc['sensor_delay_ms']); frac = sc['sensor_delay_ms']-lag
    errors = []
    for j in range(len(log)):
        errors.append(y-log[j, 2])
        inertia[j] = inertia_filter
        command = log[j-ndel, 11] if j >= ndel else log[0, 2]
        for k in range(2*ns):
            acceleration = (command-y-(t1+t2)*velocity)/(t1*t2)
            if k % 2 == 0:
                positions[j, k//2] = y
                acceleration_history[1:] = acceleration_history[:-1]
                acceleration_history[0] = acceleration
                sensed = sc['tool_mass']*((1.-frac)*acceleration_history[lag]+frac*acceleration_history[min(lag+1, 3)])
                inertia_filter += afc*(sensed-inertia_filter)
            velocity += dt*acceleration
            y += dt*velocity
    return positions, inertia, float(np.max(np.abs(errors)))


def noise_replay(sc, seed, length):
    """Independent linear sensor-noise component, with the plant's exact RNG order."""
    noise = P3._noise(np.random.default_rng(seed), int(sc['T']*1000.)+10, True)
    afc = 1.-math.exp(-2.*math.pi*P3.PLATFORM['fc']*.001)
    drift = 0.; filtered = 0.; ns = int(round(P3.PLATFORM['Tc']/.001))
    out = np.empty(length); index = 0
    for j in range(length):
        out[j] = filtered
        for _ in range(ns):
            drift += 2.e-4*sc['noise_scale']*noise[index, 0]
            filtered += afc*(drift+sc['sensor_sigma']*sc['noise_scale']*noise[index, 1]-filtered)
            index += 1
    return out


def harmonic_amplitude(t, x):
    # Remove slow scene drift; this is only a line-response diagnostic.
    tt = (t-t.mean())/max(np.ptp(t), 1.e-9)
    a = np.c_[np.cos(D.C.OMEGA*t), np.sin(D.C.OMEGA*t), np.ones(len(t)), tt, tt**2, tt**3]
    b = np.linalg.lstsq(a, x, rcond=None)[0]
    return complex(b[0], -b[1])


def one_job(case_index, i):
    name, changes, kind = CASES[case_index]
    sc = S3.cell_scene3('rest', 'relaxed', 'nominal', i, base=530000, overrides=dict(BASE, **changes))
    if case_index in (1, 7):
        assert not np.any(sc['alpha']) and sc['phi_hc'] == 0.
        assert 'memory_start' not in sc and 'equilibrium' not in sc
        assert P3.calibration(sc) == 0.
    pol = D.make('D4', (0, 0, 0), sc)
    if kind == 'L':
        # D.make sets the selected D4 structure and gate; only the contact form differs.
        pol = M.build('L', sc, D.TRACK[0])
    assert pol.kind == kind and pol.params['enc_y']
    history = {}
    estimate = pol.estimator.estimate

    def record_estimate(o, n, policy):
        history[round(float(o.log[n-1, 0]), 8)] = (policy.eb, policy.phi, policy.ndel)
        return estimate(o, n, policy)

    pol.estimator.estimate = record_estimate
    start = time.perf_counter()
    log = P3.episode(sc, pol, 900+i)
    metric = M.metric(log, pol, sc)
    if log is None:
        return dict(case=case_index, i=i, failed=True, seconds=time.perf_counter()-start)
    assert np.array_equal(log[:, 8], log[:, 20])
    true_y, true_inertia, replay_error = encoder_replay(log, sc)
    assert replay_error < 1.e-9, replay_error
    compensation = getattr(pol.step, 'comp', np.zeros(len(log)))[:len(log)]
    residual_inertia = true_inertia-compensation
    sensor_noise = noise_replay(sc, 900+i, len(log))
    # Preserve the exact H preprocessing history (its phase varies online).
    # For L, recreate ForceSensing's two input-channel operations.
    reg_log = getattr(pol.estimator, 'private', None)
    if reg_log is None:
        reg_log = log.copy(); reg_log[:, 5] -= compensation
        if hasattr(pol.step, 'delay') and pol.step.delay > 0.:
            shift = pol.step.delay*.001/pol.Tc
            reg_log[1:, 2] = log[1:, 2]-shift*np.diff(log[:, 2])
    ns = pol.ns; a = 1.-pol.afc
    weights = pol.afc*a**np.arange(ns-1, -1, -1)
    actual_Y = true_y@weights
    rows = []; eligible = 0; accepted = 0; eligible_bias = []; regression_error = 0.
    for tk, kh, sig, chi, ok in pol.accept_trace:
        j = int(np.clip(np.searchsorted(log[:, 0], tk), 0, len(log)-1))
        if tk <= 1. or log[j, 21] <= 200.:
            continue
        eligible += 1
        eligible_bias.append(math.log(kh/log[j, 21]))
        if not ok:
            continue
        accepted += 1
        ratios = [math.log(kh/log[j, c]) if log[j, c] > 0. else float('nan') for c in (21, 20, 9)]
        eb, phi, ndel = history[round(float(tk), 8)]
        # Exact source-window indexing of RidgeProjection. Diagnose the raw
        # encoder interpolation separately from H's phase preprocessing.
        n = j+1; j0 = max(ndel+1, n-1-pol.estimator.nw); m = n-1-j0
        _, y, _, keep = common._build(log, j0, m, 0, ns, eb, pol.afc, ndel,
                                      tk, log[j, 2], pol.estimator.nw*pol.Tc, True)
        matrix, _ = pol.estimator._matrix(m)
        exact = actual_Y[j0:j0+m]-weights.sum()*log[j, 2]
        slope = float(y@matrix@exact/(y@matrix@y)) if keep.all() else float('nan')
        support = complex(sc['k']-sc['m']*D.C.OMEGA**2, sc['c']*D.C.OMEGA)
        kt = log[j, 21]*np.exp(1j*log[j, 22])
        support_shift = math.log(abs(1./(1.+kt/support)))
        # From the plant's depth constraint: d=cos(theta)*(y+body+surface
        # -(Fax+m_tool*y'')/kr-Cnt*Ft). At a narrow line the numerator gains
        # 1+m_tool*w^2/kr, even after sensor inertial-force compensation.
        inertia_shift = math.log1p(sc['tool_mass']*D.C.OMEGA**2/sc['kr'])
        _, yreg, zreg, _ = common._build(reg_log, j0, m, 0, ns, eb, pol.afc, ndel,
                                       tk, reg_log[j, 2], pol.estimator.nw*pol.Tc, True)
        my = matrix@yreg; denominator = float(yreg@my)
        reconstructed_k = float(my@zreg/denominator)
        regression_error = max(regression_error, abs(reconstructed_k-kh))
        assert abs(reconstructed_k-kh) < 1.e-5*max(kh, 1.), (tk, reconstructed_k, kh)
        decay = a**ns
        def component_k(x):
            z = x[j0+1:j0+m+1]-decay*x[j0:j0+m]
            return float(my@z/denominator)
        k_residual = component_k(residual_inertia)
        k_noise = component_k(sensor_noise)
        cleaned_k = kh-k_residual-k_noise
        sensor_log = math.log(kh/(kh-k_residual)) if kh-k_residual > 0. else float('nan')
        noise_log = math.log(kh/(kh-k_noise)) if kh-k_noise > 0. else float('nan')
        clean_bias = math.log(cleaned_k/log[j, 21]) if cleaned_k > 0. else float('nan')
        rows.append(ratios+[math.log(slope) if slope > 0. else float('nan'), phi,
                            math.log(log[j, 21]/log[j, 20]), math.log(log[j, 21]/log[j, 9]), support_shift,
                            inertia_shift, ratios[0]-math.log(slope)-inertia_shift-support_shift if slope > 0. else float('nan'),
                            sensor_log, noise_log, clean_bias,
                            clean_bias-math.log(slope)-inertia_shift-support_shift if slope > 0. else float('nan'), float(tk)])
    warm = log[:, 0] > 1.
    t = log[warm, 0]
    force_line = harmonic_amplitude(t, log[warm, 18])
    y_line = harmonic_amplitude(t, log[warm, 2])
    body_line = harmonic_amplitude(t, log[warm, 4])
    return dict(case=case_index, i=i, failed=False, eligible=eligible, accepted=accepted,
                all_accept=metric['accept'], rmse=metric['rmse'], rows=rows,
                seconds=time.perf_counter()-start, replay_error=replay_error,
                line_gain=float(abs(force_line/y_line)), truth21=float(np.median(log[warm, 21])),
                body_over_y=float(abs(body_line/y_line)), phi_end=float(pol.phi),
                tau_true=sc['tau_a_true'], tau_hat=sc['tau_a_hat'], tau_identified=pol.tau,
                mass=sc['tool_mass'], mass_hat=sc['tool_mass_hat'],
                delay=sc['sensor_delay_ms'], delay_hat=sc['sensor_delay_hat'],
                eligible_bias=quantiles(eligible_bias), kr=sc['kr'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--cases', type=int, nargs='+', default=list(range(1, 9)))
    parser.add_argument('--replicas', type=int, default=6)
    args = parser.parse_args()
    chosen = json.loads((HERE/'results'/'d_compare.json').read_text(encoding='utf-8'))['chosen']['D4']
    assert chosen == dict(structure=D.STRUCT['D4'][0], chi=D.CHI[0], tracker=D.TRACK[0])
    print('CONFIG', json.dumps(chosen), flush=True)
    print('Window pooling: accepted, t>1 s, truth21>200; acceptance denominator: all eligible accept_trace entries.', flush=True)
    print('RMSE: median of per-episode M.metric axial RMSE (t>=1). IQR=q75-q25.', flush=True)
    start = time.perf_counter(); jobs = [(c-1, i) for c in args.cases for i in range(args.replicas)]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(one_job, c, i): (c, i) for c, i in jobs}
        for future in as_completed(futures):
            r = future.result(); results.append(r)
            print('RUN', json.dumps({k: v for k, v in r.items() if k != 'rows'}, allow_nan=True), flush=True)
    print('\ncase | accepted/eligible | ln(k/21) median [q25,q75] IQR | ln(k/20) | ln(k/9) | axial RMSE median [q25,q75]', flush=True)
    for c in sorted({c for c, _ in jobs}):
        rr = sorted([r for r in results if r['case'] == c], key=lambda r: r['i'])
        good = [r for r in rr if not r['failed']]
        a = np.asarray([row for r in good for row in r['rows']])
        def fmt(v):
            q = quantiles(v)
            return f"{q['median']:+.5f} [{q['q25']:+.5f},{q['q75']:+.5f}] IQR {q['iqr']:.5f}" if q['n'] else 'NA'
        na = sum(r['accepted'] for r in good); ne = sum(r['eligible'] for r in good)
        print(f"{CASES[c][0]} | {na}/{ne} ({na/max(ne,1):.3%}) | "+' | '.join(fmt(a[:, j]) for j in range(3))+' | '+fmt([r['rmse'] for r in good]), flush=True)
        diagnostics = {key: quantiles(a[:, j]) for j, key in enumerate(
            ('bias21', 'bias20', 'bias9', 'raw_encoder_interpolation_log_slope', 'phi_at_estimate',
             'log21_over20', 'log21_over9', 'continuous_support_log_shift', 'arm_inertia_log_shift',
             'bias_minus_interpolation_inertia_support', 'sensor_compensation_residual_log_contribution',
             'sensor_noise_log_contribution', 'bias_after_subtracting_sensor_residual_and_noise',
             'clean_bias_minus_interpolation_inertia_support'))}
        diagnostics['seed_bias21_medians'] = [quantiles([v[0] for v in r['rows']])['median'] for r in good]
        diagnostics['rmse_each'] = [r['rmse'] for r in good]
        diagnostics['full_trace_accept_mean'] = float(np.mean([r['all_accept'] for r in good]))
        print('DIAGNOSTICS', CASES[c][0], json.dumps(diagnostics, allow_nan=True), flush=True)
        if c != 0 and any(r['case'] == 0 for r in results):
            reference = {(r['i'], round(v[-1], 8)): v[0] for r in results if r['case'] == 0 for v in r.get('rows', [])}
            paired = [v[0]-reference[(r['i'], round(v[-1], 8))] for r in good for v in r['rows']
                      if (r['i'], round(v[-1], 8)) in reference]
            print('COMMON_ACCEPTED_BIAS_DIFFERENCE_VS_DEFAULT', CASES[c][0], json.dumps(quantiles(paired)), flush=True)
    print(f'TOTAL_SECONDS {time.perf_counter()-start:.1f}', flush=True)


if __name__ == '__main__':
    main()
