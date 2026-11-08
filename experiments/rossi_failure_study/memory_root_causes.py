"""Simulation-only root-cause tests, frozen theory parameters, development seeds.

python -X utf8 experiments/rossi_failure_study/memory_root_causes.py --workers 8
Only results/memory_root_causes.json is written. --resume accepts identical code
and reuses completed configurations (including the checked gates-on baseline).
"""
import argparse
import hashlib
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Must precede numpy, to inherit the host's single-thread BLAS environment.
import overlay_compare as O
import numpy as np
from scipy.ndimage import uniform_filter1d

import memory_filters2 as F
import memory_segment_eval as SE

HERE = Path(__file__).resolve().parent
OUT = HERE/'results/memory_root_causes.json'
CONFIGS = ('on_none', 'off_none', *F.FACTORIES)
REPS = 80


def code_hashes():
    return {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ('memory_filters2.py', 'memory_root_causes.py', 'memory_segment_eval.py',
             'memory_filters.py')}


def analytic_kalman(tau):
    a = F.H/tau
    K = (math.sqrt(a**4+4*a*a)-a*a)/2
    sd = F.S*math.sqrt(K/(2-K))
    # All 24 causal held ticks, including the update tick, equally weighted.
    offsets = np.arange(24)*.008
    time_lag = F.H*(1-K)/K+float(offsets.mean())
    return dict(K=K, stationary_sd=sd, held_time_lag=time_lag,
                held_ramps={str(r): dict(lag=r*time_lag,
                    rmse=math.sqrt(sd**2+r*r*(time_lag**2+float(offsets.var()))))
                            for r in (.29, .73)},
                continuous=dict(sd=math.sqrt(F.N/(2*tau)),
                    ramp29_lag=.29*tau, ramp73_lag=.73*tau,
                    ramp29_rmse=math.sqrt(F.N/(2*tau)+(.29*tau)**2),
                    ramp73_rmse=math.sqrt(F.N/(2*tau)+(.73*tau)**2)))


def feed(name, truth, noises):
    flt = F.FACTORIES[name]()
    out = np.empty(len(truth))
    active = []
    selected = []
    beta = []
    residual = 0.
    for j, (x, z) in enumerate(zip(truth, noises)):
        out[j], _ = flt.update(float(x+z), F.S**2, F.H*F.V_REF, F.H)
        core = getattr(flt, 'core', flt)
        complete = not hasattr(flt, 'buffer') or not flt.buffer
        if j >= 300 and complete:
            if hasattr(core, 'active'):
                active.append(core.active)
                beta.append(core.beta)
                residual = max(residual, core.last_residual)
            if hasattr(flt, 'selected'):
                selected.append(flt.selected)
    return out, dict(active=active, beta=beta, selected=selected,
                     max_eq14_relative_residual=residual)


def synthetic_job(name):
    start = time.perf_counter()
    stationary = []
    activations, betas, selected = [], [], []
    max_residual = 0.
    ramp_sums = {rho: [0., 0., 0] for rho in (.29, .73)}
    # Phase-averaged step response on the actual 8ms held-output grid.
    npost, npre = 300, 300
    ticks = np.arange(npost*24)
    step_curves = {d: np.zeros(len(ticks)) for d in (.3, 1.1)}
    offsets = np.arange(24)*.008
    for rep in range(REPS):
        rng = np.random.default_rng(920000+rep)
        z = rng.normal(0., F.S, 2000)
        m, diag = feed(name, np.zeros(2000), z)
        stationary.append(m[300:])
        activations.extend(diag['active'])
        betas.extend(diag['beta'])
        selected.extend(diag['selected'])
        max_residual = max(max_residual, diag['max_eq14_relative_residual'])
        # Common independent window noises for both ramps and all methods.
        rz = rng.normal(0., F.S, 800)
        for rho in ramp_sums:
            tr = rho*F.H*np.maximum(np.arange(800)-300, 0)
            m, _ = feed(name, tr, rz)
            err = tr[600:, None]+rho*offsets[None, :]-m[600:, None]
            accum = ramp_sums[rho]
            accum[0] += float(err.sum())
            accum[1] += float((err**2).sum())
            accum[2] += err.size
        phase = rep % 24  # fixed balanced-to-within-one sampling phases
        for d in step_curves:
            # First post-change measurement occurs phase ticks after t=0.
            tr = np.r_[np.zeros(npre), np.full(npost, d)]
            sz = rng.normal(0., F.S, npre+npost)
            m, _ = feed(name, tr, sz)
            indices = npre+np.floor_divide(ticks-phase, 24)
            step_curves[d] += m[indices]/REPS
    st = np.concatenate(stationary)
    sd = float(st.std())
    result = dict(config=name, reps=REPS, stationary_sd=sd,
                  stationary_mean=float(st.mean()), equivalent_tau=F.N/(2*sd*sd),
                  ramps={str(r): dict(lag=v[0]/v[2], rmse=math.sqrt(v[1]/v[2]))
                         for r, v in ramp_sums.items()}, steps={})
    for d, curve in step_curves.items():
        hits = np.flatnonzero(curve >= d/2)
        result['steps'][str(d)] = dict(delay50=float(hits[0]*.008) if len(hits) else None,
                                       horizon=npost*F.H)
    if activations:
        result['ard'] = dict(activation_fraction=float(np.mean(activations)),
            beta_median=float(np.median(betas)), max_eq14_relative_residual=max_residual)
    if selected:
        result['scale_mean'] = float(np.mean(selected))
        result['longest_scale_fraction'] = float(np.mean(np.array(selected) == 16))
    result['seconds'] = time.perf_counter()-start
    return result


def property_checks():
    rng = np.random.default_rng(930000)
    variances = np.geomspace(.01, .1, 16)
    values = rng.normal(size=(100000, 16))*np.sqrt(variances)
    mean = lambda n: (values[:, -n:]/variances[-n:]).sum(1)/(1/variances[-n:]).sum()
    v = lambda n: 1./(1/variances[-n:]).sum()
    pairwise = []
    for m, n in ((1, 2), (2, 8), (4, 16)):
        measured = float(np.var(mean(m)-mean(n)))
        prediction = v(m)-v(n)
        pairwise.append(dict(m=m, n=n, measured=measured, theory=prediction,
                             ratio=measured/prediction))
    time_filter, dist_filter = F.Forgetting(), F.Forgetting(distance=True)
    maxdiff = 0.
    for y in rng.normal(size=300):
        a = time_filter.update(float(y), .04, F.H*20, F.H)
        b = dist_filter.update(float(y), .04, F.H*20, F.H)
        maxdiff = max(maxdiff, abs(a[0]-b[0]), abs(a[1]-b[1]))
    # One-sample, calibrated, near-zero P test, distinct from the adaptive run.
    z = rng.normal(size=1000000)
    alpha = np.maximum(z*z-1, 0)
    delta = alpha/(alpha+1)*z
    block = rng.normal(0, .2, (100000, 4)).mean(axis=1)
    # Mechanism checks, independent of the closed-loop scores.
    mechanisms = {}
    for name in ('imm', 'ewa'):
        mechanisms[name] = {}
        for rho in (0., .29, .73):
            model = F.FACTORIES[name]()
            weights = []
            nz = np.random.default_rng(940000).normal(0., .2, 2000)
            for j, noise in enumerate(nz):
                model.update(rho*F.H*max(j-300, 0)+noise, .04)
                if j >= 1800:
                    weights.append(model.w.copy())
            mechanisms[name][str(rho)] = np.mean(weights, axis=0).tolist()
    # Diagnostic only: do not change the benchmark's paper initialization.
    prior_sensitivity = []
    for b0 in (1e-10, 1.):
        model = F.HuangARD(learn=True)
        model.b = b0
        model.beta = model.a/model.b
        output = []
        for noise in np.random.default_rng(920000).normal(0., .2, 2000):
            output.append(model.update(float(noise), .04)[0])
        prior_sensitivity.append(dict(a0=1., b0=b0,
            stationary_sd=float(np.std(output[300:])), beta_final=model.beta))
    # Isolate the floor against the original update-then-forget rule, with
    # IDENTICAL independent measurements and uninflated R in both variants.
    import memory_filters as old
    floor_ablation = {}
    for floor in (False, True):
        curves = {d: np.zeros(150) for d in (.3, 1.1)}
        for rep in range(80):
            nz = np.random.default_rng(950000+rep).normal(0., .2, 2150)
            for d in curves:
                model = F.Forgetting() if floor else old.Forgetting(lam_max=.998)
                for j, noise in enumerate(nz):
                    obs = float(noise+(d if j >= 2000 else 0.))
                    m, _ = (model.update(obs, .04) if floor else model.update(obs, .04, 0., True))
                    if j >= 2000:
                        curves[d][j-2000] += m/80
        floor_ablation[str(floor)] = {}
        for d, curve in curves.items():
            hits = np.flatnonzero(curve >= d/2)
            floor_ablation[str(floor)][str(d)] = float(hits[0]*F.H) if len(hits) else None
    return dict(nested_difference_variance=pairwise,
                mechanism_mean_weights=mechanisms,
                ARD_initial_precision_diagnostic=prior_sensitivity,
                forgetting_floor_ablation_update_aligned_delay50=floor_ablation,
                distance_time_equality_at_20mm_s_max_abs=maxdiff,
                distance_Q_per_update_by_speed={str(speed): F.Q_DISTANCE*F.H*speed
                                                for speed in (5., 10., 20., 40.)},
                floor_time_Q_per_update=F.Q_RATE_LO*F.H,
                calibrated_single_sample_ARD=dict(activation=float(np.mean(alpha > 0)),
                    update_variance_over_R=float(np.var(delta)),
                    expected_activation=.317310507862914, expected_variance=.333261882350746),
                four_window_mean=dict(measured_variance=float(np.var(block)), theory=.04/4),
                exact_kalman={str(t): analytic_kalman(t) for t in (2., .2)},
                section2_continuous_optimum=[dict(rho=r, tau=(F.N/(4*r*r))**(1/3),
                    rmse=math.sqrt(F.N/(2*(F.N/(4*r*r))**(1/3))+
                                   (r*(F.N/(4*r*r))**(1/3))**2)) for r in (.05, .29, .73)])


def signed_rate(log):
    """Same smoothing/contact runs as SE.segment_masks; sign for lag only."""
    t, k = log[:, 0], log[:, 21]
    contact = (k > 0.) & (log[:, 3] > 0.) & np.isfinite(k)
    rate = np.full(len(t), np.nan)
    starts = np.flatnonzero(contact & ~np.r_[False, contact[:-1]])
    ends = np.flatnonzero(contact & ~np.r_[contact[1:], False])+1
    width = max(1, int(round(.2/float(np.median(np.diff(t))))))
    for a, b in zip(starts, ends):
        if b-a >= 2:
            smooth = uniform_filter1d(k[a:b], size=width, mode='nearest')
            rate[a:b] = np.gradient(np.log(smooth), t[a:b])
    return rate


def job(task):
    config, cell, variant, q = task
    if config == 'on_none':
        # The original job and original baseline check, without any wrapper.
        return SE.job((config, cell, variant, q, 1))
    O.setup()
    sc = SE.MC.scene_for(cell, SE.MC.DEV_BASE, variant)
    pol = O.B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    SE.gates_off(pol)
    state = dict(raw=np.nan)
    estimate = pol.estimator.estimate

    def raw_capture(o, n, policy):
        out = estimate(o, n, policy)
        state['raw'] = out[0]
        return out

    pol.estimator.estimate = raw_capture
    flt = F.attach(pol, config) if config != 'off_none' else None
    records = []
    step = pol.step

    def capture(o):
        state['raw'] = np.nan
        command = step(o)
        records.append((float(o.log[o.n-1, 0]), pol.k_good, pol.kbar, state['raw']))
        return command

    pol.step = capture
    log = SE.P3.episode(sc, pol, SE.MC.SEED+q)
    if log is None:
        raise RuntimeError(f'Plant failed: {task}')
    rec = np.asarray(records)
    assert len(rec) == len(log)
    assert np.allclose(rec[:, 0], log[:, 0], rtol=0., atol=1e-9)
    metric = SE.M.metric(log, pol, sc)
    eligible, masks = SE.segment_masks(log, sc)
    with np.errstate(divide='ignore', invalid='ignore'):
        errors = np.log(rec[:, 1:4]/log[:, 21, None])
    peak_j = int(np.argmax(np.where(log[:, 0] >= 1., log[:, 18], -np.inf)))
    peak_segment = next((s for s in SE.SEGMENTS if masks[s][peak_j]), 'unscored')
    row = dict(config=config, scene=cell[0], variant=variant,
               dynamics=cell[2], phase=cell[3], i=cell[4], seed=SE.MC.SEED+q,
               L=float(np.log1p(metric['J'])), peak=float(metric['peak']),
               peak_segment=peak_segment, eligible_ticks=int(eligible.sum()), segments={})
    # Intentionally identical scoring to SE.job, including all force samples.
    for s in SE.SEGMENTS:
        mask = masks[s]
        force_error = log[mask, 18]-8.
        row['segments'][s] = dict(
            ticks=int(mask.sum()), share=float(mask.sum()/eligible.sum()),
            e_g=SE.distribution(errors[mask, 0], fractions=True),
            e_bar=SE.distribution(errors[mask, 1]), e_raw=SE.distribution(errors[mask, 2]),
            force_rmse=float(np.sqrt(np.mean(force_error**2))) if len(force_error) else None,
            force_p99=float(np.percentile(abs(force_error), 99)) if len(force_error) else None,
            peak_over_15=int(row['peak'] > 15. and peak_segment == s),
            _force_abs=abs(force_error).tolist())
    rate = signed_rate(log)
    row['transition_lag'] = {
        direction: SE.distribution(errors[masks['transition'] & condition, 0])
        for direction, condition in (('rising', rate > 0), ('falling', rate < 0))}
    if flt is not None:
        sampling = flt.sampling
        gaps = np.diff(sampling['updates'])
        assert not len(gaps) or min(gaps) >= 24
        row['sampling'] = dict(updates=len(sampling['updates']),
            min_gap_ticks=int(min(gaps)) if len(gaps) else None,
            max_gap_ticks=int(max(gaps)) if len(gaps) else None,
            max_hold_ticks=sampling['max_hold_ticks'], invalid=sampling['invalid'])
    row['observed_scan_path_mm_s'] = float(np.abs(np.diff(log[:, 1])).sum()*1000/
                                          (log[-1, 0]-log[0, 0]))
    return row


def aggregate_config(rows, config):
    # Reuse the original aggregation without changing SE.CONFIGS or any host.
    rr = [dict(r, config='off_none') for r in rows]
    out = SE.aggregate(rr)
    for group in out:
        group['config'] = config
        selected = [r for r in rows if (r['scene'], r['variant']) ==
                    (group['scene'], group['variant'])]
        group['peak_over_15'] = sum(r['peak'] > 15. for r in selected)
        if config != 'on_none':
            group['transition_lag'] = {d: SE.mean_present(
                [r['transition_lag'][d]['median'] for r in selected])
                for d in ('rising', 'falling')}
    return out


def save(saved):
    OUT.write_text(json.dumps(saved, ensure_ascii=False, indent=1, allow_nan=False),
                   encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    start = time.perf_counter()
    hashes = code_hashes()
    episodes = [(c, 'nominal', q) for q, c in enumerate(O.cells())
                if c[0] in ('breath_rib', 'breath_fast')]
    episodes += [(c, 'slow_thin', q) for q, c in enumerate(O.cells()) if c[0] == 'breath_rib']
    assert len(episodes) == 84
    if args.resume:
        saved = json.loads(OUT.read_text(encoding='utf-8'))
        assert saved['metadata']['code_sha256'] == hashes, 'Code changed; resume rejected.'
    else:
        saved = dict(metadata=dict(code_sha256=hashes, dev_base=SE.MC.DEV_BASE,
            sensor_seed_base=SE.MC.SEED, fresh_seeds=False, workers=args.workers,
            episode_counts=dict(A_nominal=36, A_slow_thin=36, B=12),
            hop_ticks=24, dt=F.H, R='(sigma/k)^2, no overlap inflation',
            q_rate_lo=F.Q_RATE_LO, q_rate_hi=F.Q_RATE_HI, q_distance=F.Q_DISTANCE,
            q_reference_noise_sd=F.S, q_fixed_across_speeds_and_observed_R=True,
            additional_sample_hold_latency='0..23 ticks (0..0.184s); conservative bound <0.192s; gaps from unavailable windows reported separately',
            block='arithmetic mean of 4 consecutive nonoverlapping windows, sum(R)/16; block outputs every .768s',
            host='SE.gates_off unchanged; on_none uses original SE.job and SE.check_baseline',
            off_none='original every-tick no-memory baseline; off_window24 is additional matched-rate control',
            scoring='SE.segment_masks/distribution/aggregate unchanged; force uses every tick without trimming',
            lag='mean of episode median e_g during transition, rising/falling by sign of same smoothed k derivative',
            synthetic='80 common-noise replicates, independent .2 sd measurements every .192s; holds scored on 8ms grid; stationary 2000 updates discard 300; ramps last 200 of 800 updates; step mean-response 50% delay with phase cycling mod24',
            ARD='E covariance scaled by beta, actual P=E/beta; fixed beta=1; learn a0=1,b0=1e-10; exact joint evidence maximum; D4 state replacement/reference layer off'),
            properties=property_checks(), synthetic=[], episodes=[], summary=[],
            pooled_summary=[], completed=[])
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        present = {r['config'] for r in saved['synthetic']}
        futures = {ex.submit(synthetic_job, name): name for name in F.FACTORIES if name not in present}
        for f in as_completed(futures):
            result = f.result()
            saved['synthetic'].append(result)
            print('synthetic', result['config'], 'sd', round(result['stationary_sd'], 5), flush=True)
            save(saved)
        # Freeze the properties before any host comparisons.
        syn = {r['config']: r for r in saved['synthetic']}
        saved['properties']['block_memory_ratio'] = {
            'fixed_precision': syn['ard_block']['equivalent_tau']/syn['ard_fixed']['equivalent_tau'],
            'learned_precision': syn['ard_block_learn']['equivalent_tau']/syn['ard_learn']['equivalent_tau'],
            'predicted': 4.}
        for name, tau in (('kal_long', 2.), ('kal_short', .2)):
            theory = saved['properties']['exact_kalman'][str(tau)]
            assert abs(syn[name]['stationary_sd']/theory['stationary_sd']-1.) < .04
            for rho in ('.29', '.73'):
                r = str(float(rho))
                assert abs(syn[name]['ramps'][r]['lag']-theory['held_ramps'][r]['lag']) < .01
        assert saved['properties']['distance_time_equality_at_20mm_s_max_abs'] < 1e-12
        for config in CONFIGS:
            if config in saved['completed']:
                continue
            if config != 'on_none':
                assert saved['baseline_check']['passed']
            futures = [ex.submit(job, (config, c, v, q)) for c, v, q in episodes]
            rows = []
            for f in as_completed(futures):
                rows.append(f.result())
                if len(rows) % 12 == 0:
                    print(config, len(rows), '/84', flush=True)
            if config == 'on_none':
                saved['baseline_check'] = SE.check_baseline(rows)
                print('baseline_check', {k: v for k, v in saved['baseline_check'].items()
                                          if k != 'checks'}, flush=True)
                if not saved['baseline_check']['passed']:
                    save(saved)
                    raise SystemExit('Gates-on baseline failed; no candidates run.')
            saved['summary'].extend(aggregate_config(rows, config))
            # Compact overall report: same per-episode masks and scoring,
            # all 84 episodes combined, episode-equal distributions as before.
            pooled = aggregate_config([dict(r, scene='breath_rib', variant='nominal')
                                       for r in rows], config)[0]
            pooled['scene'], pooled['variant'] = 'all_84', 'dev_only'
            saved['pooled_summary'].append(pooled)
            for row in rows:
                for part in row['segments'].values():
                    del part['_force_abs']
            saved['episodes'].extend(sorted(rows, key=lambda r: (r['variant'], r['scene'], r['seed'])))
            saved['completed'].append(config)
            saved['seconds'] = saved.get('seconds', 0.)+time.perf_counter()-start
            start = time.perf_counter()
            save(saved)
        old = json.loads((HERE/'results/memory_segment_eval.json').read_text(encoding='utf-8'))
        key = lambda r: (r['scene'], r['variant'], r['dynamics'], r['phase'], r['i'])
        refs = {key(r): r for r in old['episodes'] if r['config'] == 'off_none'}
        rr = [r for r in saved['episodes'] if r['config'] == 'off_none']
        saved['off_none_reproduction'] = dict(episodes=len(rr),
            max_abs_L_delta=max(abs(r['L']-refs[key(r)]['L']) for r in rr),
            max_abs_peak_delta=max(abs(r['peak']-refs[key(r)]['peak']) for r in rr))
        save(saved)
    print('saved', OUT, 'episodes', len(saved['episodes']), flush=True)


if __name__ == '__main__':
    main()
