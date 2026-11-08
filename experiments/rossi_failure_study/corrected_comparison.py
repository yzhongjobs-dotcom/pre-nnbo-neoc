"""Correct cached preprocessing, verify legacy replay, then a fixed paired grid.

python -X utf8 experiments/rossi_failure_study/corrected_comparison.py --workers 8
Simulation only; no parameter search, commits, or pushes.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace

import overlay_compare as O
import numpy as np
import rib_transition_eval as R
import memory_segment_eval as SE
import contact_candidates as C
from ridge_estimator import RidgeProjection

HERE = Path(__file__).resolve().parent
OUT = HERE/'results/corrected_comparison.json'
BASE, SENSOR_BASE = 874500, 910000
SETTINGS = (('neural_white', 28.), ('viscous', 14.))
MEMORY = ('kal_long', 'kal_short', 'imm', 'lepski', 'ewa', 'forget_time', 'ard_fixed')
METHODS = ('Huber-D4', 'D4', 'Rossi-ID', 'Kalman', 'MHE', 'off_none') + MEMORY + O.B.CTRL_METHODS + ('oracle',)
SCENES = ('A', 'A_slow_thin', 'B', 'A+B', 'cough', 'combo')


def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [clean(x) for x in v]
    if isinstance(v, (float, np.floating)):
        return float(v) if np.isfinite(v) else None
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def save(data):
    OUT.write_text(json.dumps(clean(data), ensure_ascii=False, allow_nan=False,
                              separators=(',', ':')), encoding='utf-8')


def cases():
    selected = []
    overlay, original = O.cells(), O.B.cells()
    for label, cond, indices, variant in (
        ('A', 'breath_rib', (0, 2, 6), 'nominal'),
        ('A_slow_thin', 'breath_rib', (1, 4, 5), 'slow_thin'),
        ('B', 'breath_fast', (0, 1, 2), 'nominal'),
        ('A+B', 'breath_rib_fast', (0, 2, 6), 'nominal'),
        ('cough', 'cough', (0, 3, 7), 'nominal'),
        ('combo', 'combo', (0, 3, 7), 'nominal')):
        source = original if cond in ('cough', 'combo') else overlay
        for dy, ph in SETTINGS:
            for i in indices:
                cell, = [c for c in source if c[0] == cond and c[2:] == (dy, ph, i)]
                selected.append(dict(index=len(selected), group=label, cell=cell,
                                     variant=variant, seed=SENSOR_BASE+len(selected)))
    assert len(selected) == 36 and len({c['seed'] for c in selected}) == 36
    return selected


def fixture(est, owner, log):
    est.reset(owner)
    obs = SimpleNamespace(log=log)
    for n in range(1, len(log)+1):
        est._log(obs, n)
    return est.private[:len(log), 2].copy()


def numerical_checks():
    y = .018+np.cumsum(np.r_[0., np.sin(np.arange(1, 1000)*.71)*6e-6])
    log = np.zeros((len(y), 26)); log[:, 2] = y; log[:, 5] = 8.
    results = []
    # HC, every memory form (including the sliding variants), and both Ridge
    # preprocessing stages. Other signal columns do not affect this fixture.
    for kind in C.HC+C.MEMORY:
        for scale in (0., .7, 2.):
            own = SimpleNamespace(tau_h=scale*.008, mag_h=1+scale,
                a_y=np.full(1 if kind == 'S' else len(C.POWER_TAU), scale),
                e_y=np.full(1 if kind == 'S' else len(C.POWER_TAU), .6),
                b_y=np.full(1 if kind == 'S' else len(C.POWER_TAU), .8),
                mag_y=1+scale, cf=np.zeros(len(y)))
            old, new = C.Preprocess(), C.Preprocess()
            old.h_increment = old.memory_increment = False
            old._setup(kind, own); new._setup(kind, own)
            obs = SimpleNamespace(log=log)
            for n in range(1, len(y)+1):
                old._log(obs, n); new._log(obs, n)
            a, b = old.private[:len(y), 2], new.private[:len(y), 2]
            err = float(max(abs(np.diff(a)-np.diff(b))))
            assert err < 3e-17, (kind, err)
            results.append(dict(path='Contact:'+kind, scale=scale,
                                max_increment_difference_m=err))
    owner = SimpleNamespace(Tc=.008, params={})
    for memory, lam in ((True, 0.), (False, .05), (True, .05)):
        old = RidgeProjection(memory=memory, lam=lam, increment=False)
        new = RidgeProjection(memory=memory, lam=lam)
        old.reset(owner); new.reset(owner)
        old.k_ref = new.k_ref = 30000.
        obs = SimpleNamespace(log=log)
        for n in range(1, len(y)+1):
            old._log(obs, n); new._log(obs, n)
        a, b = old.private[:len(y), 2], new.private[:len(y), 2]
        err = float(max(abs(np.diff(a)-np.diff(b))))
        assert err < 3e-17, (memory, lam, err)
        results.append(dict(path='Ridge', memory=memory, lam=lam,
                            max_increment_difference_m=err))
    # A step after a finite motion, with zero current/previous velocity and
    # zero branch state: isolates coefficient artifacts from physical relaxation.
    rest = np.zeros((1003, 26)); rest[:, 2] = np.r_[np.linspace(0., .018, 1000), [.018]*3]
    rest[:, 5] = 8.; obs = SimpleNamespace(log=rest)
    steps = []
    for kind in C.HC+C.MEMORY:
        own = SimpleNamespace(tau_h=0., mag_h=1., a_y=np.zeros(1 if kind == 'S' else 5),
            e_y=np.zeros(1 if kind == 'S' else 5), b_y=np.ones(1 if kind == 'S' else 5),
            mag_y=1., cf=np.zeros(len(rest)))
        pair = [C.Preprocess(), C.Preprocess()]
        pair[0].h_increment = pair[0].memory_increment = False
        for est in pair:
            est._setup(kind, own); est._log(obs, 1002)
        before = [est.y_out for est in pair]
        cached = pair[1].private[:1002].copy()
        own.tau_h = .01; own.mag_h = own.mag_y = 1.001; own.a_y[:] = 2.
        for est in pair:
            est._log(obs, 1003)
        delta = [est.y_out-a for est, a in zip(pair, before)]
        assert delta[1] == 0. and np.array_equal(cached, pair[1].private[:1002])
        steps.append(dict(path='Contact:'+kind, old_m=delta[0], corrected_m=delta[1], cached_rows_unchanged=True))
    for memory, lam in ((True, 0.), (False, .05), (True, .05)):
        pair = [RidgeProjection(memory=memory, lam=lam, increment=False, alpha=np.zeros(3)),
                RidgeProjection(memory=memory, lam=lam, alpha=np.zeros(3))]
        for est in pair:
            est.reset(owner)
            est.k_ref = 30000.; est._log(obs, 1002)
        before = [est.private[1001, 2] for est in pair]
        cached = pair[1].private[:1002].copy()
        for est in pair:
            if memory:
                est.mag *= 1.001
            est.k_ref = 15000.; est._log(obs, 1003)
        delta = [est.private[1002, 2]-a for est, a in zip(pair, before)]
        assert delta[1] == 0. and np.array_equal(cached, pair[1].private[:1002])
        steps.append(dict(path='Ridge', memory=memory, lam=lam, old_m=delta[0],
                          corrected_m=delta[1], cached_rows_unchanged=True))
    return dict(constant_coefficients=results, rest_coefficient_steps=steps,
        max_increment_difference_m=max(r['max_increment_difference_m'] for r in results),
        memory_rest='Zero branch state and zero current/previous dy; remaining branch relaxation at rest is physical, not an exactly-zero response.')


def legacy_check():
    O.setup()
    cell = ('breath_rib', 'rib', 'neural_white', 28., 2)
    stored = next(r for r in O.ORIGINAL_LOAD('overlay_test.json')['rows']
                  if r['method'] == 'Huber-D4' and O.B.key(r) == cell)
    sc = O.B.scene(cell, 674500, 0.)
    pol = O.B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    pol.estimator.h_increment = pol.estimator.memory_increment = False
    seed = O.sensor_seed(cell)
    log = SE.P3.episode(sc, pol, seed)
    metric = SE.M.metric(log, pol, sc)
    with np.load(HERE/'results'/stored['log_file']) as z:
        assert np.array_equal(log[:, z['cols']], z['log']), 'Legacy trajectory differs'
        count, columns = z['log'].shape
    for field in ('J', 'rmse', 'peak'):
        assert metric[field] == stored[field], (field, metric[field], stored[field])
    return dict(passed=True, cell=cell, base=674500, seed=seed,
        stored_file=stored['log_file'], bitwise_equal_rows=count, columns=columns,
        J=metric['J'], rmse=metric['rmse'], peak=metric['peak'],
        J_delta=0., rmse_delta=0., peak_delta=0.)


def job(task):
    method, case = task
    O.setup()
    cell = tuple(case['cell'])
    sc = SE.MC.scene_for(cell, BASE, case['variant'])
    assert sc['tool_mass'] > 0. and sc['tool_mass_hat'] > 0.
    host = 'Huber-D4' if method in ('off_none', 'oracle')+MEMORY else method
    pol = O.B.build(host, O.chosen()[host], sc)
    identifier = host not in O.B.CTRL_METHODS
    if identifier:
        assert getattr(pol.estimator, 'h_increment', True)
    if method in ('off_none', 'oracle')+MEMORY:
        SE.gates_off(pol)
    sampling = R.make_filter(pol, method) if method in MEMORY else None
    oracle_truth = []
    if method == 'oracle':
        original_estimate = pol.estimator.estimate

        def oracle(o, n, policy):
            original_estimate(o, n, policy)
            k = float(o.truth_k)
            pol.estimator.fit = dict(pol.estimator.fit, info=1.)
            oracle_truth.append((n-1, k))
            return 1.2*k, 0., 0.

        pol.estimator.estimate = oracle
    channels = ('k_good', 'k_bar') if identifier else ('k_good',) if method in ('Rossi-control', 'DO-MPAC') else ()
    records = []
    step = pol.step

    def capture(o):
        command = step(o)
        if identifier:
            kg, kb = pol.k_good, pol.kbar
        else:
            kg = getattr(pol, 'k_hat', getattr(pol, 'k_h0', np.nan)); kb = np.nan
        records.append((float(o.log[o.n-1, 0]), kg, kb))
        return command

    pol.step = capture
    tic = time.perf_counter()
    log = R.oracle_run(sc, pol, case['seed']) if method == 'oracle' else SE.P3.episode(sc, pol, case['seed'])
    if log is None:
        raise RuntimeError(('Plant failed', method, case))
    rec = np.asarray(records)
    assert len(rec) == len(log) and np.allclose(rec[:, 0], log[:, 0], rtol=0., atol=1e-9)
    metric = SE.M.metric(log, pol, sc)
    # R.score accepts two estimate channels. Absent controller channels use an
    # offline placeholder solely to avoid all-NaN reductions, then are removed
    # from every metric. These values are never passed to any controller.
    scoring = rec.copy()
    for c, name in enumerate(('k_good', 'k_bar')):
        if name not in channels:
            scoring[:, c+1] = np.maximum(log[:, 21], 1.)
    row = dict(method=method, config=method, **case, base=BASE,
        L=float(np.log1p(metric['J'])), force_rmse=metric['rmse'], peak=metric['peak'],
        catastrophic=bool(metric['peak'] > 36.5), seconds=time.perf_counter()-tic,
        tool_mass=sc['tool_mass'], tool_mass_hat=sc['tool_mass_hat'],
        estimate_channels=channels,
        estimate_meaning='fixed nominal k_h0, not an identified gain' if method == 'DO-MPAC' else
                         'no stiffness estimate' if method == 'Xing-DAHIC' else
                         'own stiffness estimate' if method == 'Rossi-control' else 'Method k_good/k_bar')
    row.update(R.score(log, scoring, sc))
    row['plateau'] = {k: v for k, v in row['plateau'].items() if k in channels or k == 'force_band_rms'}
    for event in row['transitions']:
        event['channels'] = {k: v for k, v in event['channels'].items() if k in channels}
        a = event['start']; b = int(np.searchsorted(log[:, 0], event['t1']+.5+1e-9))
        event['force']['peak_true_N'] = float(max(log[a:b, 18]))
    if sampling is not None:
        gaps = np.diff(sampling['updates'])
        assert not len(gaps) or gaps.min() >= 6
        row['sampling'] = dict(updates=len(sampling['updates']), invalid=sampling['invalid'],
            min_gap_ticks=int(gaps.min()) if len(gaps) else None,
            max_gap_ticks=int(gaps.max()) if len(gaps) else None,
            gaps_over_six=int((gaps > 6).sum()))
    if method == 'oracle':
        difference = max((abs(k-log[n, 21]) for n, k in oracle_truth), default=0.)
        assert difference < 1e-8
        row['oracle_max_truth_delta'] = difference
    if identifier:
        assert pol.estimator.h_increment
        row['corrected_preprocessing'] = True
    else:
        row['corrected_preprocessing'] = 'Own controller equations; no H/MEMORY/Ridge preprocessing'
    return clean(row)


def summary(rows):
    groups = []
    for scene in SCENES:
        for method in METHODS:
            rr = [r for r in rows if r['group'] == scene and r['method'] == method]
            groups.append(dict(group=scene, method=method, n=len(rr),
                L_mean=float(np.mean([r['L'] for r in rr])),
                L_median=float(np.median([r['L'] for r in rr])),
                force_rmse_median=float(np.median([r['force_rmse'] for r in rr])),
                max_peak=max(r['peak'] for r in rr), above_15=sum(r['peak'] > 15. for r in rr)))
    transitions, plateau = [], []
    for method in METHODS:
        rr = [r for r in rows if r['method'] == method]
        for ch in ('k_good', 'k_bar'):
            pp = [r['plateau'][ch] for r in rr if ch in r['plateau']]
            plateau.append(dict(method=method, channel=ch,
                episode_count=sum(p['n'] > 0 for p in pp), ticks=sum(p['n'] for p in pp),
                bias=R.dist([p['median'] for p in pp]),
                robust_sd=R.dist([p['robust_sd'] for p in pp])))
        rr = [r for r in rr if r['group'] in ('A', 'A_slow_thin', 'A+B')]
        for direction in ('rising', 'falling'):
            detected = [e for r in rr for e in r['transitions']
                        if e['size_class'] == 'large' and e['direction'] == direction]
            events = [e for r in rr if not r['catastrophic'] for e in r['transitions']
                      if e['size_class'] == 'large' and e['direction'] == direction]
            item = dict(method=method, direction=direction, n=len(events),
                detected=len(detected), too_few=len(events) < 5,
                excluded_catastrophic=len(detected)-len(events),
                excluded_prior_contact_loss=sum(1 for r in rr for e in r['excluded_prior_contact_loss']
                    if e['size_class'] == 'large' and e['direction'] == direction),
                force={k: R.dist([e['force'][k] for e in events])
                       for k in ('peak_abs', 'peak_true_N', 'positive', 'negative')}, channels={})
            for ch in ('k_good', 'k_bar'):
                cc = [e['channels'][ch] for e in events if ch in e['channels']]
                item['channels'][ch] = {k: R.dist([c[k] for c in cc]) for k in
                    ('lag', 'negative', 'under_magnitude', 'raw_negative', 'raw_under_magnitude')}
                item['channels'][ch].update(
                    settling=R.dist([c['settling']['time'] for c in cc]),
                    settling_censored=sum(c['settling']['censored'] for c in cc),
                    settling_short_horizon=sum(c['settling']['observed_horizon'] < 2.-1e-8 for c in cc),
                    lag_at_upper_bound=sum(c['lag_at_upper_bound'] for c in cc),
                    bias_fallback=sum(c['bias_fallback'] for c in cc))
            transitions.append(item)
    return groups, transitions, plateau


def code_audit():
    return [
        dict(path='contact_candidates.Preprocess._log/HC', status='fixed, default h_increment=True',
             live='Huber-D4/D4 plus Rossi-ID/Kalman/MHE via ContactAdapter; all gates-off variants and oracle; H/HA/HT and P5 share fix'),
        dict(path='contact_candidates.Preprocess._log/MEMORY', status='fixed, default memory_increment=True',
             live='W/S/WF/WA/WT; none of this comparison uses these contact forms; memory filters in the comparison are log-gain filters, not contact MEMORY'),
        dict(path='ridge_estimator.RidgeProjection._log/memory', status='fixed, default increment=True',
             live='Its own optional branch; comparison Huber/D4 override _log with Preprocess; other comparison identifiers use ContactAdapter'),
        dict(path='ridge_estimator.RidgeProjection._log/lam', status='fixed, default increment=True',
             live='Optional force/k_ref-dependent HC normalization, inactive in the requested configurations'),
        dict(path='ridge_estimator.RidgeProjection.estimate/gamma', status='unchanged',
             reason='Explicit window-end static-stiffening reference; window regressors recalculated, no transformed level cached; gamma=0 in requested methods'),
        dict(path='contact_candidates.CandidatePolicy._next_z; contact_compare.ModelPolicy._next_z; sliding', status='unchanged',
             reason='Current inverse-response or force/sliding increments integrated once; no episode-origin renormalization'),
        dict(path='id_baselines.EncoderMatched; estimators._StreamingBandPass; id_baselines.EquationEstimator', status='unchanged',
             reason='Causal filtering/recursion, fixed references/filter coefficients; receives corrected H rows from ContactAdapter'),
        dict(path='method3.ForceSensing', status='unchanged',
             reason='Recalculates observation view using fixed calibrated delay and mass; acceleration compensation is causal; not online level renormalization'),
        dict(path='ctrl_baselines.RossiController/DOMPACController/XingController', status='bypasses H/MEMORY/Ridge',
             reason='Complete independent controllers, their own equations; Rossi fixed-coefficient BP/EWRLS, DO fixed nominal gain, Xing no gain estimate'),
        dict(path='low_estimate_forensics.H_tests; h_fix_eval.numerical_checks/metadata', status='historical diagnostics unchanged',
             reason='Offline old-form conditional reconstruction and explicit old/new ablation metadata; not comparison execution paths. h_fix_eval old default assumptions are superseded; use this script for current checks')]


def hashes():
    names = ('corrected_comparison.py', 'contact_candidates.py', 'ridge_estimator.py',
        'id_baselines.py', 'd_robust.py', 'estimators.py', 'ctrl_baselines.py',
        'method3.py', 'contact_compare.py', 'estimator_swap.py',
        'rib_transition_eval.py', 'memory_segment_eval.py', 'memory_filters2.py',
        'memory_closed_loop.py', 'overlay_compare.py', 'baseline_compare.py',
        'scenes3.py', 'physics_plant3.py', 'results/baseline_test_r2.json')
    import comparators
    import common
    paths = [(n, HERE/n) for n in names]
    paths += [(str(Path(m.__file__).relative_to(HERE.parent.parent)), Path(m.__file__))
              for m in (comparators, common)]
    return {n: hashlib.sha256(p.read_bytes()).hexdigest() for n, p in paths}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--checks-only', action='store_true')
    args = ap.parse_args()
    tic = time.perf_counter()
    source_hashes = hashes()
    if args.resume:
        data = json.loads(OUT.read_text(encoding='utf-8'))
        assert data['metadata']['source_sha256'] == source_hashes, 'Source changed; refuse resume'
    else:
        data = dict(metadata=dict(base=BASE, sensor_seed_base=SENSOR_BASE,
            cases=cases(), methods=METHODS, chosen=O.chosen(), workers=args.workers,
            simulation_only=True, probe_mass=True, tuning=False, source_sha256=source_hashes,
            sensor_index='Sequential 0..35 in listed scene, setting, i order, identical across methods',
            force_metrics='method3.metric: t>=1s true axial force; L=ln(1+J), J capped at1000, force RMSE against8N; peak=max true F, >15 counted per episode',
            memory='rib_transition_eval.make_filter unchanged: >=6 ticks, R=4*(sigma/k)^2; IMM/EWA likelihood power1/4; forget_time includes Q_RATE_LO time floor; no forget_dist variant',
            oracle='1.2*instantaneous true |K(20.8)| as window output through gates-off host, sigma=chi=0; unchanged side effects; truth verified at each call',
            transitions='rib_transition_eval.score unchanged; pooled A nominal, A slow_thin, A+B, large abs(delta log k)>=.3; rising/falling separately. Exact truth for scoring, .1s boxcar for detection, rate>.2/s, gaps<.1s merged.',
            transition_exclusions='Original rules: prior true-contact loss excludes event; episode force peak>36.5 excludes all its transitions. All exclusions and catastrophic episodes retained in JSON.',
            transition_metrics='lag grid0:.008:1.496 plus1.5s; local .5s pretransition median log bias removed (episode fallback recorded); settling abs(local error)<=.1 for .2s, cap2s; interval[t0,t1+.5]. Under-estimate = -min(log(hat/k)-local bias), signed; max across events is largest. Raw error also retained. Force peak_abs=max|F-8|; true F max also retained.',
            plateau='SE plateau masks; ln(hat/k) bias=per-episode median, sd=(Q75-Q25)/1.349; reported medians over episodes with plateau, including catastrophes. DO-MPAC is nominal fixed k_h0 only; Xing has no stiffness estimate.',
            too_few='n<5 marked; counts refer to method-specific closed-loop trajectories, not matched events'),
            numerical_checks=numerical_checks(), code_path_audit=code_audit(), episodes=[])
        data['legacy_check'] = legacy_check()
        save(data)
        print('checks passed', data['numerical_checks']['max_increment_difference_m'], data['legacy_check'], flush=True)
    if args.checks_only:
        return
    present = {(r['method'], r['index']) for r in data['episodes']}
    tasks = [(m, c) for m in METHODS for c in cases() if (m, c['index']) not in present]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        pending = [ex.submit(job, task) for task in tasks]
        for i, future in enumerate(as_completed(pending), 1):
            data['episodes'].append(future.result())
            if i % 12 == 0 or i == len(tasks):
                save(data)
                print('episodes', len(data['episodes']), '/', len(METHODS)*36,
                      'seconds', round(time.perf_counter()-tic, 1), flush=True)
    data['episodes'].sort(key=lambda r: (METHODS.index(r['method']), r['index']))
    assert len(data['episodes']) == 36*len(METHODS) == 612
    assert all(sum(r['method'] == m for r in data['episodes']) == 36 for m in METHODS)
    assert hashes() == source_hashes, 'Code changed during run'
    data['scene_summary'], data['transition_summary'], data['plateau_summary'] = summary(data['episodes'])
    data['seconds_this_invocation'] = time.perf_counter()-tic
    data['validation'] = dict(passed=True, episodes=612, methods=17, episodes_per_method=36,
        episodes_per_scene_method=6, corrected_identifier_episodes=504,
        independent_controller_episodes=108, memory_spacing_episodes=252,
        catastrophes=sum(r['catastrophic'] for r in data['episodes']),
        numerical_checks_passed=True, legacy_episode_bitwise_equal=True)
    save(data)
    print('saved', OUT, flush=True)


if __name__ == '__main__':
    main()
