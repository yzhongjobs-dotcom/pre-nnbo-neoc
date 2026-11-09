"""Verify frozen scenes/configuration and one old row for each saved method."""
import hashlib
import json
import argparse
from concurrent.futures import ProcessPoolExecutor

import identification_completion as IC
import numpy as np


def digest(value):
    return hashlib.sha256(json.dumps(IC.CC.clean(value), sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def synthetic_transition_tests():
    def log_for(t, ell):
        log = np.zeros((len(t), IC.P3.LOGC))
        log[:, 0], log[:, 3], log[:, 21], log[:, 26] = t, 1., 1., 1000.*np.exp(ell)
        return log

    t = np.arange(0., 4.008, .008)
    ell = np.interp(t, [0., 1.2, 1.8, 2.4, 3., 3.6, 4.], [0., 0., 1., 0., 1., 0., 0.])
    _, _, events, excluded = IC.R.transitions(log_for(t, ell))
    assert not excluded and [e['direction'] for e in events] == ['rising', 'falling']*2
    assert all((e['delta'] > 0) == (j % 2 == 0) for j, e in enumerate(events))

    t = np.arange(0., 4.001, .001)
    ell = np.interp(t, [0., 1.2, 2., 2.11, 2.9, 4.], [0., 0., .8, .8, 1.59, 1.59])
    smoothed, _, ramp, excluded = IC.R.transitions(log_for(t, ell))
    runs = IC.R.spans((t >= 1.) & (np.gradient(smoothed, t) > .2))
    assert len(runs) == 2
    gap = float(t[runs[1][0]]-t[runs[0][1]])
    assert .045 <= gap <= .055, gap
    assert not excluded and len(ramp) == 1 and ramp[0]['direction'] == 'rising'

    t = np.arange(0., 3., .008)
    for limit in (.9, 1., 1.1, 1.199):
        result = IC.R.hold_time(np.zeros(len(t)), t, 1., IC.LOG12, limit=limit)
        assert result['censored'] and result['censored_by'] == 'next_change'
        assert abs(result['observed_horizon']-max(0., limit-1.)) < 1e-12
    result = IC.R.hold_time(np.zeros(len(t)), t, 1., IC.LOG12, limit=1.2)
    assert not result['censored'] and result['time'] == 0.
    result = IC.R.hold_time(np.ones(len(t)), t, 1., IC.LOG12, limit=np.inf)
    assert result['censored_by'] == 'horizon' and result['observed_horizon'] == 2.
    result = IC.R.hold_time(np.zeros(len(t)), t, 2.9, IC.LOG12, limit=np.inf)
    assert result['censored_by'] == 'episode_end'
    return dict(triangle_flanks=4, signs=['rising', 'falling']*2,
                ramp_events=1, ramp_subthreshold_gap_s=gap, bounded_hold_cases=7)


def independent_recovery(error, t, end, next_t0):
    """Brute-force tick windows, independent of hold_time/spans."""
    dt = float(np.median(np.diff(t)))
    need = int(np.ceil(.2/dt-1e-8))
    limit = np.inf if next_t0 is None else next_t0
    boundary = min(limit, float(t[-1]+dt))
    horizon = max(0., min(2., boundary-end))
    reason = 'horizon' if horizon >= 2.-1e-9 else 'next_change' if limit <= t[-1]+dt else 'episode_end'
    if horizon > 0:
        candidates = np.flatnonzero((t >= end-1e-9) & (t <= end+2.+1e-9))
        for j in candidates:
            if j+need > len(t) or t[j+need-1]+dt > boundary+1e-9:
                continue
            values = error[j:j+need]
            if np.isfinite(values).all() and (abs(values) <= IC.LOG12).all():
                return dict(time=max(0., float(t[j]-end)), censored=False,
                            observed_horizon=horizon, censored_by=None)
    return dict(time=horizon, censored=True, observed_horizon=horizon, censored_by=reason)


def saved_results(saved):
    closed = json.loads((IC.HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))
    keyed = {(r['method'], r['index']): r for r in closed['rows']}
    assert len(keyed) == len(closed['rows'])
    for name in IC.TABLE1_METHODS:
        indices = range(36, 42) if name == IC.CC.MARGIN else range(42)
        assert all((name, i) in keyed for i in indices), name
    data = json.loads((IC.HERE/'results/identification_replay.json').read_text(encoding='utf-8'))
    assert len(data['cases']) == 84
    errors, source_deltas, metric_deltas, under_deltas = [], [], [], []
    real_checks = []
    for case in data['cases']:
        meta, log, ticks = IC.load_case(IC.TICKS/f"{case['source']}_{case['index']:02d}.npz")
        assert meta['methods'] == list(IC.SOURCES[case['source']])
        errors.append(meta['command_replay_error'])
        prior = keyed[case['source'], case['index']]
        source_deltas += [abs(meta['source_L']-prior['L']), abs(meta['source_rmse']-prior['force_rmse']), abs(meta['source_peak']-prior['peak'])]
        contact = IC.SE.contact_mask(log)
        ell = np.full(len(log), np.nan)
        ell[contact] = np.log(log[contact, 26])
        for name, result in case['methods'].items():
            for channel, key in (('raw', 'raw'), ('k_good', 'good')):
                values = ticks[name+'__'+key].astype(float)
                if channel == 'raw':
                    values = np.where(ticks[name+'__valid'], values, np.nan)
                with np.errstate(divide='ignore', invalid='ignore'):
                    error = np.log(values)-ell
                for event in result['transitions']:
                    in_change = (log[:, 0] >= event['t0']-1e-9) & (log[:, 0] <= event['t1']+1e-9) & np.isfinite(error)
                    direct = float(np.mean(abs(error[in_change]))) if in_change.any() else None
                    reported = event['channels'][channel]['change_mae']
                    assert (direct is None) == (reported is None)
                    if direct is not None:
                        metric_deltas.append(abs(direct-reported))
                    under = event['channels'][channel]['raw_under_magnitude']
                    assert under is None or under >= 0.
                    next_t0 = np.inf if event['next_t0'] is None else event['next_t0']
                    window = (log[:, 0] >= event['t0']) & (log[:, 0] < min(event['t1']+.5, next_t0)) & np.isfinite(error)
                    direct_under = max(0., -float(error[window].min())) if window.any() else None
                    assert (direct_under is None) == (under is None)
                    if under is not None:
                        under_deltas.append(abs(direct_under-under))
                if case['source'] == 'Huber-D4' and case['index'] in (36, 37, 39) and name in IC.MAIN[:6]+('Rossi-ID',):
                    event = next(e for e in result['transitions'] if e['size_class'] == 'large')
                    direct = independent_recovery(error, log[:, 0], event['t1'], event['next_t0'])
                    reported = event['channels'][channel]['recovery_raw']
                    assert reported['censored'] == direct['censored']
                    assert reported['censored_by'] == direct['censored_by']
                    assert abs(reported['time']-direct['time']) < 1e-12
                    assert abs(reported['observed_horizon']-direct['observed_horizon']) < 1e-12
                    real_checks.append(dict(index=case['index'], t0=event['t0'], method=name,
                        channel=channel, recovery=direct))
    assert max(errors) < 1e-10
    assert max(source_deltas) == 0.
    assert max(metric_deltas, default=0.) < 1e-12
    assert max(under_deltas, default=0.) < 1e-12
    for row in data['tracking']:
        for ch in IC.CHANNELS:
            p = row[ch]
            assert p['recovered']+p['unrecovered']+p['truncated']+p['truncated_next_change'] == row['n']
    assert len(real_checks) == 3*7*2
    saved['transition_tests'] = dict(synthetic=synthetic_transition_tests(),
        real_A_large_events=3, real_channel_method_checks=len(real_checks), real_checks=real_checks)
    saved['replay_checks'] = dict(cases=84, command_error_bound=1e-10,
        command_max_abs_error=max(errors), source_metric_max_abs_delta=max(source_deltas),
        change_mae_checks=len(metric_deltas), change_mae_max_abs_delta=max(metric_deltas, default=0.),
        raw_under_checks=len(under_deltas), raw_under_max_abs_delta=max(under_deltas, default=0.),
        recovery_counts_partition_events=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved-only', action='store_true', help='Validate completed data without rerunning any plant')
    parser.add_argument('--synthetic-only', action='store_true', help='Run detector/window tests without reading or simulating episodes')
    args = parser.parse_args()
    if args.synthetic_only:
        print('PASS synthetic', synthetic_transition_tests(), flush=True)
        return
    IC.O.setup()
    path = IC.HERE/'results/arms_reproduction_checks.json'
    saved = json.loads(path.read_text(encoding='utf-8'))
    cases = IC.CC.cases()
    assert IC.CC.clean(cases[:36]) == saved['old_cases']
    assert hashlib.sha256(json.dumps(IC.O.chosen(), sort_keys=True).encode()).hexdigest() == saved['chosen_sha256']
    for row in saved['scene_hashes']:
        c = cases[row['index']]
        assert digest(IC.SE.MC.scene_for(tuple(c['cell']), IC.CC.BASE, c['variant'])) == row['sha256']
    saved['all_36_scenes_unchanged'] = True
    saved['group_hash_checks'] = [next(r for r in saved['scene_hashes'] if r['group'] == g) for g in IC.CC.SCENES[:-1]]
    shared = ('xcmd', 'q_field', 'E_field', 'h_s', 'u', 'act', 'fext', 'rip', 'k_rib')
    for c in cases[36:]:
        original = IC.SE.MC.scene_for(tuple(c['cell']), IC.CC.BASE, 'nominal')
        large = IC.SE.MC.scene_for(tuple(c['cell']), IC.CC.BASE, 'large_rib')
        for key in shared:
            assert digest(original[key]) == digest(large[key]), key
        assert digest(original['ribs'][:, :2]) == digest(large['ribs'][:, :2])
    saved['A_large_shared_field_motion_rib_centres_widths'] = True
    old = json.loads((IC.HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']
    for r in old:
        key = r['method']+'__'+str(r['index'])
        if key in saved['row_hashes']:
            assert hashlib.sha256(json.dumps(r, sort_keys=True).encode()).hexdigest() == saved['row_hashes'][key]
    saved['stored_rows_preserved'] = len(saved['row_hashes'])
    if args.saved_only:
        saved_results(saved)
        path.write_text(json.dumps(saved, ensure_ascii=False, indent=2), encoding='utf-8')
        print('PASS: preserved scenes/rows/configuration and completed replay metrics', saved['replay_checks'], flush=True)
        return
    methods = list(dict.fromkeys(r['method'] for r in old if r['method']+'__'+str(r['index']) in saved['row_hashes']))
    # Ideal-calibration job returns its underlying method; caller restores its label.
    tasks = [(m, cases[0]) for m in methods]
    with ProcessPoolExecutor(8) as ex:
        fresh = list(ex.map(IC.closed_job, tasks, chunksize=1))
    checks = []
    for (m, c), r in zip(tasks, fresh):
        prior = next(x for x in old if x['method'] == m and x['index'] == c['index'])
        delta = {k: abs(r[k]-prior[k]) for k in ('L', 'force_rmse', 'peak', 'loss_time')}
        assert max(delta.values()) == 0., (m, delta)
        checks.append(dict(method=m, index=c['index'], max_abs_deltas=delta))
    saved['closed_loop_checks'] = checks
    saved['passed'] = True
    path.write_text(json.dumps(saved, ensure_ascii=False, indent=2), encoding='utf-8')
    print('PASS: 36 scene hashes, frozen configuration,', len(checks), 'method rows; all metric deltas zero', flush=True)


if __name__ == '__main__':
    main()
