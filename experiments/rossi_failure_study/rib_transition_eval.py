"""Transition-level simulation evaluation; only this file and OUT are written.

python -X utf8 experiments/rossi_failure_study/rib_transition_eval.py --workers 8
The untouched on_none must pass SE.check_baseline before any other method runs.
No tuning, fresh seeds, selection, git operations, or truth access by estimators
other than the explicitly labelled oracle. --resume checks source hashes.
"""
import os
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import inspect
import json
import math
import textwrap
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import overlay_compare as O  # sets BLAS thread limits before numpy
import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import butter, sosfiltfilt
import memory_segment_eval as SE
import memory_filters2 as F
import memory_filters as MF  # original overlap convention, not its old filters

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results/rib_transition_eval.json'
CONFIGS = ('on_none', 'off_none', 'oracle', 'kal_long', 'kal_short', 'imm',
           'lepski', 'ewa', 'forget_time', 'forget_dist', 'ard_fixed',
           'D4', 'Rossi-ID', 'Kalman', 'MHE')
MEMORY = CONFIGS[3:11]
HOP = 6


def spans(mask):
    return list(zip(np.flatnonzero(mask & ~np.r_[False, mask[:-1]]),
                    np.flatnonzero(mask & ~np.r_[mask[1:], False])+1))


def dist(values):
    a = np.asarray([x for x in values if x is not None], float)
    a = a[np.isfinite(a)]
    if not len(a):
        return dict(n=0, median=None, p10=None, p90=None, min=None, max=None)
    return dict(n=len(a), median=float(np.median(a)),
                p10=float(np.percentile(a, 10)), p90=float(np.percentile(a, 90)),
                min=float(a.min()), max=float(a.max()))


def transitions(log):
    """Exact 0.1 s centered boxcar of log truth, within contact runs only."""
    t, k = log[:, 0], log[:, 21]
    dt = float(np.median(np.diff(t)))
    contact = (k > 0) & np.isfinite(k) & (log[:, 3] > 0)
    ell = np.full(len(t), np.nan)
    radius = int(math.ceil(.05/dt))
    centres = np.arange(-radius, radius+1)*dt
    kernel = np.maximum(0., np.minimum(centres+dt/2, .05)-
                        np.maximum(centres-dt/2, -.05))/.1
    assert abs(kernel.sum()-1) < 1e-10
    rate = ell.copy()
    for a, b in spans(contact):
        ell[a:b] = convolve1d(np.log(k[a:b]), kernel, mode='nearest')
        if b-a >= 2:
            rate[a:b] = np.gradient(ell[a:b], t[a:b])
    runs = spans(contact & (t >= 1.) & (abs(rate) > .2))
    merged = []
    for a, b in runs:
        # Count the duration of the intervening sub-threshold ticks.
        if merged and (a-merged[-1][1])*dt < .1-1e-10 and contact[merged[-1][1]:a].all():
            merged[-1][1] = b
        else:
            merged.append([a, b])
    losses = np.flatnonzero(~contact & np.r_[False, np.maximum.accumulate(contact)[:-1]])
    first_loss = int(losses[0]) if len(losses) else len(t)
    result, excluded = [], []
    for a, b in merged:
        delta = float(ell[b-1]-ell[a])
        if abs(delta) < .15:
            continue
        row = dict(start=int(a), end=int(b-1), t0=float(t[a]), t1=float(t[b-1]),
                   duration=float(t[b-1]-t[a]), delta=delta, size=abs(delta),
                   direction='rising' if delta > 0 else 'falling',
                   size_class='large' if abs(delta) >= .3 else 'small')
        (excluded if first_loss < a else result).append(row)
    return ell, contact, result, excluded


def make_filter(pol, kind):
    """F2 theory constants; actual elapsed time and observed path in mm.

    Same nw/hop overlap convention as MF.attach. IMM/EWA gains use R=4s²;
    likelihoods use prior+s², then power 1/4. This keeps gain correction and
    correlated evidence tempering separate rather than inflating it twice.
    Raw windows still execute at every request, preserving estimator side effects.
    """
    flt = F.FACTORIES[kind]()
    if kind in ('imm', 'ewa'):
        original_update = F.IMM.update if kind == 'imm' else F.EWA.update
        flt.update = SE.bind_source(flt, original_update, [
            ('innovation_var = prior+R', 'innovation_var = prior+R/4.'),
            ('-.5*(', '-.5*.25*('),
            ('K = prior/innovation_var', 'K = prior/(prior+R)')])
    est, state = pol.estimator, dict(last=None, out=None, index=None, ds=0.,
                                    updates=[], invalid=0)
    original = est.estimate

    def estimate(o, n, policy):
        k, sig, chi = original(o, n, policy)
        if state['index'] is not None:
            state['ds'] += float(abs(np.diff(o.log[state['index']:n, 1])).sum())*1000
        state['index'] = n-1
        valid = (np.isfinite(k) and k > 0 and np.isfinite(sig) and sig > 0
                 and est.fit.get('info', 0.) > 0)
        due = state['last'] is None or n-state['last'] >= HOP
        if valid and due:
            dt = HOP*pol.Tc if state['last'] is None else (n-state['last'])*pol.Tc
            state['out'] = flt.update(math.log(max(k, 200.)), 4*(sig/k)**2,
                                      state['ds'], dt)
            state['last'] = n
            state['ds'] = 0.
            state['updates'].append(n)
        elif not valid:
            state['invalid'] += 1
        if state['out'] is None:
            return k, sig, chi
        m, p = state['out']
        kf = math.exp(np.clip(m, math.log(200.), math.log(300000.)))
        return kf, kf*math.sqrt(max(p, 0.)), chi

    est.estimate = estimate
    return state


def oracle_run(sc, pol, seed):
    """Private plant runner supplies instantaneous truth to oracle only.

    A shadow tick on copied state gives the exact truth at the current logged
    tick, avoiding the normal diagnostic's one-tick-old truth. The actuator has
    at least one tick of delay, so the current command cannot affect this value.
    The real tick, sensor stream and host run unchanged. Exact agreement with
    the real tick is asserted after the episode; shadow work is not timed as a
    feasible estimator. Nothing is exposed to any non-oracle policy.
    """
    source = textwrap.dedent(inspect.getsource(SE.P3._run))
    old = """            if diagnostic and n > 0:
                obs.truth_k = float(log[n-1, 21]); obs.truth_re = float(log[n-1, 20]); obs.truth_k1 = float(log[n-1, 9]); obs.truth_body = float(log[n-1, 4])"""
    new = """            if diagnostic:
                shadow = np.zeros((int(round(pl['Tc']/log_dt))+2, LOGC))
                _tick(s.copy(), buf.copy(), float(buf[-1]), c, geom,
                      sc['q_field'], sc['E_field'], sc['ribs'], sc['h_s'],
                      sc['xcmd'], sc['u'], sc['fext'], sc['act'], sc['rip'],
                      alpha, tau, noise, kk, shadow, 0)
                obs.truth_k = float(shadow[0, 21])"""
    assert source.count(old) == 1
    namespace = {}
    exec(compile(source.replace(old, new), '<oracle-only runner>', 'exec'),
         SE.P3._run.__globals__, namespace)
    return namespace['_run'](sc, seed, True, None, policy=pol, diagnostic=True)


def hold_time(error, t, end, tolerance):
    """First start of a continuous >=0.2 s in-band run, up to 2 s.

    Durations include a tick's support; a start at t1 can have zero recovery.
    Censored values are capped at 2, including right censoring at episode end;
    the separate observed horizon makes that censoring visible.
    """
    dt = float(np.median(np.diff(t)))
    need = int(math.ceil(.2/dt-1e-8))
    start = int(np.searchsorted(t, end-1e-9))
    good = np.isfinite(error) & (abs(error) <= tolerance)
    for a, b in spans(good[start:]):
        a, b = a+start, b+start
        if b-a >= need and t[a]-end <= 2.+1e-9:
            return dict(time=max(0., float(t[a]-end)), censored=False,
                        observed_horizon=min(2., float(t[-1]+dt-end)))
    return dict(time=2., censored=True,
                observed_horizon=max(0., min(2., float(t[-1]+dt-end))))


def force_metrics(error, band, t, a, b, recovery_end, plateau_band):
    z, bz = error[a:b], band[a:b]
    rms = float(np.sqrt(np.mean(bz*bz)))
    return dict(peak_abs=float(max(abs(z))), positive=float(max(z)),
                negative=float(min(z)), band_rms=rms,
                band_ratio=rms/plateau_band if plateau_band > 0 else None,
                recovery=hold_time(error, t, recovery_end, .5))


def score(log, rec, sc):
    t = log[:, 0]
    _, contact, events, excluded = transitions(log)
    dt = float(np.median(np.diff(t)))
    _, masks = SE.segment_masks(log, sc)
    plateau = masks['plateau']
    truth_raw = np.full(len(t), np.nan)
    truth_raw[contact] = np.log(log[contact, 21])
    # Smooth only for locating transitions. Scoring uses ell=ln|K| itself:
    # otherwise even an exact oracle receives spurious errors near curvature.
    ell = truth_raw
    with np.errstate(divide='ignore', invalid='ignore'):
        hat = np.log(rec[:, 1:3])
    error = log[:, 18]-8.
    band = sosfiltfilt(butter(4, [15., 35.], btype='bandpass',
                             fs=1/dt, output='sos'), error)
    plateau_band = float(np.sqrt(np.mean(band[plateau]**2))) if plateau.any() else 0.
    episode_bias = [float(np.nanmedian((hat[:, c]-ell)[contact & (t >= 1.)]))
                    for c in range(2)]
    plateau_result = {name: SE.distribution((hat[:, c]-truth_raw)[plateau])
                      for c, name in enumerate(('k_good', 'k_bar'))}
    plateau_result['force_band_rms'] = plateau_band
    scored = []
    plateau_runs = spans(plateau)
    for event in events:
        a, z = event['start'], event['end']
        b = int(np.searchsorted(t, event['t1']+.5+1e-9))
        window_end = event['t1']+.5
        row = dict(event, interval_ticks=b-a,
                   interval_truncated=bool(t[-1]+dt < window_end-1e-9),
                   contact_loss_in_interval=bool(not contact[a:b].all()),
                   recovery_window_next_transition=any(
                       e['t0'] > event['t1'] and e['t0'] < event['t1']+2.2
                       for e in events))
        row['channels'] = {}
        pre = (t >= event['t0']-.5-1e-9) & (t < event['t0']-1e-9)
        idx = np.arange(a, b)
        tau_grid = np.r_[np.arange(0., 1.5, dt), 1.5]
        # Do not interpolate across a contact boundary or endpoint extension.
        shifted = np.array([np.interp(t[idx]-tau, t, ell, left=np.nan, right=np.nan)
                            for tau in tau_grid])
        for c, name in enumerate(('k_good', 'k_bar')):
            pre_error = (hat[:, c]-ell)[pre]
            fallback = (len(pre_error) < int(round(.5/dt))-1 or
                        not np.isfinite(pre_error).all())
            bias = episode_bias[c] if fallback else float(np.median(pre_error))
            local = hat[:, c]-ell-bias
            lag_errors = np.mean(abs(hat[idx, c][None, :]-shifted-bias), axis=1)
            finite = np.flatnonzero(np.isfinite(lag_errors))
            j = int(finite[np.argmin(lag_errors[finite])]) if len(finite) else None
            raw = (hat[:, c]-truth_raw)[a:b]
            row['channels'][name] = dict(
                bias=bias, bias_fallback=fallback,
                lag=float(tau_grid[j]) if j is not None else None,
                lag_at_upper_bound=bool(j is not None and j == len(tau_grid)-1),
                lag_mae=float(lag_errors[j]) if j is not None else None,
                settling=hold_time(local, t, event['t1'], .1),
                negative=float(np.nanmin(local[a:b])),
                positive=float(np.nanmax(local[a:b])),
                under_magnitude=float(-np.nanmin(local[a:b])),
                raw_negative=float(np.nanmin(raw)),
                raw_under_magnitude=float(-np.nanmin(raw)))
        row['force'] = force_metrics(error, band, t, a, b, event['t1'], plateau_band)
        # Choose the nearest interval of IDENTICAL tick length, entirely inside
        # one contiguous plateau run; never pool disconnected plateau samples.
        length = b-a
        choices = []
        for pa, pb in plateau_runs:
            if pb-pa >= length:
                start = int(np.clip(a, pa, pb-length))
                choices.append(start)
        if choices:
            pa = min(choices, key=lambda x: (abs(x-a), x))
            pb = pa+length
            virtual_end = float(t[pa]+event['t1']-event['t0'])
            matched = force_metrics(error, band, t, pa, pb, virtual_end, plateau_band)
            row['matched_plateau'] = dict(t0=float(t[pa]), t1=float(t[pb-1]),
                interval_ticks=length, force=matched,
                excess={k: row['force'][k]-matched[k] if row['force'][k] is not None and matched[k] is not None else None
                        for k in ('peak_abs', 'positive', 'negative', 'band_rms', 'band_ratio')},
                recovery_excess=row['force']['recovery']['time']-matched['recovery']['time'])
        else:
            row['matched_plateau'] = None
        scored.append(row)
    return dict(transitions=scored, excluded_prior_contact_loss=excluded,
                plateau=plateau_result, plateau_ticks=int(plateau.sum()))


def job(task):
    config, cell, variant, q = task
    O.setup()
    sc = SE.MC.scene_for(cell, SE.MC.DEV_BASE, variant)
    method = config if config in CONFIGS[11:] else 'Huber-D4'
    pol = O.B.build(method, O.chosen()[method], sc)
    if config in ('off_none', 'oracle')+MEMORY:
        SE.gates_off(pol)
    sampling = make_filter(pol, config) if config in MEMORY else None
    oracle_truth = []
    if config == 'oracle':
        original = pol.estimator.estimate

        def oracle_estimate(o, n, policy):
            # Preserve the window estimator's side effects, replacing its output.
            original(o, n, policy)
            k = float(o.truth_k)
            pol.estimator.fit = dict(pol.estimator.fit, info=1.)
            oracle_truth.append((n-1, k))
            return 1.2*k, 0., 0.

        pol.estimator.estimate = oracle_estimate
    records = []
    original_step = pol.step

    def capture(o):
        command = original_step(o)
        records.append((float(o.log[o.n-1, 0]), pol.k_good, pol.kbar))
        return command

    pol.step = capture
    seed = SE.MC.SEED+q
    log = oracle_run(sc, pol, seed) if config == 'oracle' else SE.P3.episode(sc, pol, seed)
    if log is None:
        raise RuntimeError(f'Plant failed: {task}')
    rec = np.asarray(records)
    assert len(log) == len(rec) and np.allclose(log[:, 0], rec[:, 0], atol=1e-9, rtol=0)
    metric = SE.M.metric(log, pol, sc)
    row = dict(config=config, cell=list(cell), scene=cell[0], variant=variant,
               dynamics=cell[2], phase=cell[3], i=cell[4], seed=seed,
               L=float(np.log1p(metric['J'])), peak=float(metric['peak']),
               catastrophic=bool(metric['peak'] > 36.5))
    # Keep detection and plateau records, but all transition metrics from a
    # catastrophic episode are marked excluded by aggregation.
    row.update(score(log, rec, sc))
    if sampling is not None:
        gaps = np.diff(sampling['updates'])
        assert not len(gaps) or gaps.min() >= HOP
        row['sampling'] = dict(updates=len(sampling['updates']), invalid=sampling['invalid'],
            min_gap_ticks=int(gaps.min()) if len(gaps) else None,
            max_gap_ticks=int(gaps.max()) if len(gaps) else None,
            gaps_over_six=int((gaps > HOP).sum()))
    if config == 'oracle':
        differences = [abs(k-log[n, 21]) for n, k in oracle_truth]
        assert max(differences, default=0.) < 1e-8
        row['oracle_max_truth_delta'] = max(differences, default=0.)
    return row


def summarize_events(rows, size_class, direction):
    events = [e for r in rows if not r['catastrophic'] for e in r['transitions']
              if e['size_class'] == size_class and e['direction'] == direction]
    out = dict(n=len(events), episodes=sum(any(e['size_class'] == size_class and
               e['direction'] == direction for e in r['transitions'])
               for r in rows if not r['catastrophic']), channels={})
    for channel in ('k_good', 'k_bar'):
        cc = [e['channels'][channel] for e in events]
        out['channels'][channel] = {k: dist([c[k] for c in cc]) for k in
            ('bias', 'lag', 'lag_mae', 'negative', 'positive', 'under_magnitude',
             'raw_negative', 'raw_under_magnitude')}
        out['channels'][channel].update(
            settling=dist([c['settling']['time'] for c in cc]),
            settling_censored=sum(c['settling']['censored'] for c in cc)/len(cc) if cc else None,
            settling_short_horizon=sum(c['settling']['observed_horizon'] < 2.-1e-8 for c in cc),
            bias_fallback=sum(c['bias_fallback'] for c in cc),
            lag_at_upper_bound=sum(c['lag_at_upper_bound'] for c in cc))
    out['force'] = {k: dist([e['force'][k] for e in events]) for k in
                    ('peak_abs', 'positive', 'negative', 'band_rms', 'band_ratio')}
    out['force'].update(recovery=dist([e['force']['recovery']['time'] for e in events]),
        recovery_censored=sum(e['force']['recovery']['censored'] for e in events)/len(events) if events else None,
        recovery_short_horizon=sum(e['force']['recovery']['observed_horizon'] < 2.-1e-8 for e in events))
    matched = [e for e in events if e['matched_plateau'] is not None]
    out['matched_plateau'] = dict(n=len(matched), missing=len(events)-len(matched),
        force={k: dist([e['matched_plateau']['force'][k] for e in matched]) for k in
               ('peak_abs', 'positive', 'negative', 'band_rms', 'band_ratio')},
        recovery=dist([e['matched_plateau']['force']['recovery']['time'] for e in matched]),
        excess={k: dist([e['matched_plateau']['excess'][k] for e in matched]) for k in
                ('peak_abs', 'positive', 'negative', 'band_rms', 'band_ratio')},
        recovery_excess=dist([e['matched_plateau']['recovery_excess'] for e in matched]))
    out.update(interval_truncated=sum(e['interval_truncated'] for e in events),
               contact_loss_in_interval=sum(e['contact_loss_in_interval'] for e in events),
               recovery_window_next_transition=sum(e['recovery_window_next_transition'] for e in events))
    return out


def summarize(rows):
    groups = []
    counts = []
    for variant in ('nominal', 'slow_thin'):
        for config in CONFIGS:
            rr = [r for r in rows if r['variant'] == variant and r['config'] == config]
            if not rr:
                continue
            group = dict(variant=variant, config=config, episodes=len(rr),
                catastrophes=[dict(cell=r['cell'], seed=r['seed'], peak=r['peak'])
                              for r in rr if r['catastrophic']],
                plateau={ch: {k: dist([r['plateau'][ch][k] for r in rr])
                               for k in ('median', 'robust_sd')}
                         for ch in ('k_good', 'k_bar')}, transitions={})
            group['plateau_noncatastrophic'] = {ch: {k: dist([
                r['plateau'][ch][k] for r in rr if not r['catastrophic']])
                for k in ('median', 'robust_sd')} for ch in ('k_good', 'k_bar')}
            for size_class in ('large', 'small'):
                for direction in ('rising', 'falling'):
                    key = size_class+'_'+direction
                    group['transitions'][key] = summarize_events(rr, size_class, direction)
                    all_events = [e for r in rr for e in r['transitions']
                                  if (e['size_class'], e['direction']) == (size_class, direction)]
                    counts.append(dict(variant=variant, config=config, category=key,
                        detected=len(all_events), retained=group['transitions'][key]['n'],
                        excluded_catastrophic=len(all_events)-group['transitions'][key]['n'],
                        excluded_prior_contact_loss=sum(1 for r in rr for e in r['excluded_prior_contact_loss']
                            if (e['size_class'], e['direction']) == (size_class, direction)),
                        duration=dist([e['duration'] for e in all_events]),
                        size=dist([e['size'] for e in all_events])))
            groups.append(group)
    return groups, counts


def hashes():
    paths = [Path(__file__), Path(O.__file__), Path(O.B.__file__), Path(SE.__file__),
             Path(SE.MC.__file__), Path(SE.M.__file__), Path(SE.P3.__file__),
             Path(SE.ES.__file__), Path(SE.C.__file__), Path(SE.CC.__file__),
             Path(F.__file__), Path(MF.__file__), HERE/'id_baselines.py',
             HERE/'d_robust.py', HERE/'results/baseline_test_r2.json',
             HERE/'results/memory_closed_loop_true_base.json']
    return {str(p.relative_to(HERE)) if p.is_relative_to(HERE) else str(p):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def save(saved):
    OUT.write_text(json.dumps(saved, ensure_ascii=False, allow_nan=False,
                              separators=(',', ':')), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--baseline-only', action='store_true')
    args = parser.parse_args()
    start = time.perf_counter()
    code_hashes = hashes()
    cells = [(q, c) for q, c in enumerate(O.cells()) if c[0] == 'breath_rib']
    assert len(cells) == 36 and [q for q, _ in cells] == list(range(36))
    episodes = [(c, v, q) for v in ('nominal', 'slow_thin') for q, c in cells]
    if args.resume:
        saved = json.loads(OUT.read_text(encoding='utf-8'))
        assert saved['metadata']['code_sha256'] == code_hashes, 'Changed code; resume rejected'
    else:
        saved = dict(metadata=dict(code_sha256=code_hashes, configs=list(CONFIGS),
            base=SE.MC.DEV_BASE, sensor_seed_base=SE.MC.SEED, fresh_seeds=False,
            probe_mass=True, workers=args.workers, cells_per_variant=36,
            host='SE.gates_off for off_none/oracle/memory; on_none and B.build comparison estimators unchanged',
            update='Memory: every 6 actual control ticks when valid, no block or 24-tick hold; original estimators called each host request',
            baseline_rate='on_none/off_none and adopted D4/Rossi-ID/Kalman/MHE remain untouched; oracle at every host request',
            gain_variance='R=4*(sigma/k)^2 (nw/hop=24/6, MF.attach convention)',
            likelihood='IMM/EWA: innovation variance prior+(sigma/k)^2, likelihood power 1/4; gains use prior+4*(sigma/k)^2',
            parameters='F2 factories and constants unchanged; actual dt and path; no tuning',
            q_rate_lo=F.Q_RATE_LO, q_rate_hi=F.Q_RATE_HI, q_distance=F.Q_DISTANCE,
            oracle='Window output 1.2*instantaneous true |K(20.8)|, sigma=0, chi=0, info=1; through unchanged gates-off host including CENTRE and k_bar decline bound; not a feasible estimator',
            centering_factor=SE.CC.CENTRE,
            truth='0.1 s centered boxcar of ln|K(20.8)| within contact runs for detection ONLY; exact unsmoothed ell for offset, lag, local error, raw log ratio and plateau',
            transitions='t>=1, contact, abs(dell/dt)>.2/s, merge subthreshold gaps <.1s, abs(endpoint delta)>=.15; large>=.3',
            prior_contact_loss='After first true contact, any true-contact loss before transition start excludes that transition; initial no-contact acquisition is not a loss',
            catastrophe='Episode peak true axial force after t>=1 >36.5N: retained and listed, excluded from all transition aggregates and paired comparisons; plateau all-episode and noncatastrophic versions both stored',
            local_offset='Median of own error on full .5s before start; otherwise eligible episode median; counted separately per channel',
            lag='Exhaustive grid 0:.008:1.496 plus 1.5s; linear interpolation of exact log truth; invalid lags omitted; ties choose smallest',
            settling='First start at/after t1 with .2s continuous in-band; <=.1 log unit, cap2s; end-of-episode censor horizon recorded',
            force='True axial force log[:,18]-8; fourth-order Butterworth 15-35Hz zero-phase full-episode offline filter; interval band RMS / pooled episode plateau band RMS',
            force_recovery='First start at/after t1 with abs(F-8)<=.5N for .2s, cap2s; next transitions do not reset horizon; overlap counted',
            matched_plateau='Nearest contiguous SE plateau interval of identical tick length in same episode (ties earlier); no stitched ticks or cross-episode substitution; virtual t1 at same offset, missing matches counted',
            distributions='Each transition has equal weight; med/P90 incl capped censored times. Signed negative extremes also P10; under_magnitude=-min(error), so P90 is adverse tail. Plateau distributions give each episode equal weight.'),
            episodes=[], completed=[])
    def run_config(ex, config):
        present = {(r['variant'], r['seed']) for r in saved['episodes'] if r['config'] == config}
        tasks = [(config, c, v, q) for c, v, q in episodes if (v, SE.MC.SEED+q) not in present]
        futures = [ex.submit(job, task) for task in tasks]
        for i, future in enumerate(as_completed(futures), 1):
            saved['episodes'].append(future.result())
            if i % 12 == 0 or i == len(tasks):
                print(config, i, '/', len(tasks), flush=True)
                save(saved)
        saved['completed'].append(config)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        if 'on_none' not in saved['completed']:
            run_config(ex, 'on_none')
        baseline = [r for r in saved['episodes'] if r['config'] == 'on_none']
        assert len(baseline) == 72
        saved['baseline_check'] = SE.check_baseline(baseline)
        print('baseline_check', {k: v for k, v in saved['baseline_check'].items() if k != 'checks'}, flush=True)
        save(saved)
        assert saved['baseline_check']['passed'], 'STOP: on_none failed baseline reproduction'
        if args.baseline_only:
            return
        for config in CONFIGS[1:]:
            if config not in saved['completed']:
                run_config(ex, config)
    saved['episodes'].sort(key=lambda r: (CONFIGS.index(r['config']), r['variant'], r['seed']))
    saved['summary'], saved['transition_counts'] = summarize(saved['episodes'])
    saved['seconds_this_invocation'] = time.perf_counter()-start
    save(saved)
    print('saved', OUT, 'episodes', len(saved['episodes']), flush=True)


if __name__ == '__main__':
    main()
