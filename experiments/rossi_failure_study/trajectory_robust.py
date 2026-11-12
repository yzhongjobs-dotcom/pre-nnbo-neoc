"""Fixed round 3 robustness grid. Import round 2 arms, never alter their constants.
python -X utf8 trajectory_robust.py --project-root REPO --output-root MIRROR
"""
import os
import sys
sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
if '--project-root' in sys.argv:
    os.environ['TRAJECTORY_LOOP_PROJECT'] = sys.argv[sys.argv.index('--project-root')+1]
from pathlib import Path
PROJECT = Path(os.environ.get('TRAJECTORY_LOOP_PROJECT', Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(PROJECT/'experiments/rossi_failure_study'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import trajectory_loop as TL
import argparse
import hashlib
import inspect
import json
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
CC, SE, O, B, P3, M, D = TL.CC, TL.SE, TL.O, TL.B, TL.P3, TL.M, TL.D
ARMS = ('huber', 'delag')
PREDICTIONS = dict(
    R1='Fresh large: mean delta <0, max delta <=.02, >=9/12 delta <=.005.',
    R2='New seeds: every delta <=.02 and mean delta <0.',
    R3='Every eps large mean delta <0; delag spread for 36/37 <round2 gap.',
    R4='Expected positive event-window force RMS delta; full-episode mean L delta <=.05, all peaks <=15 N; away separately.',
    R5='Fresh nominal: every abs(delta) <=.01 and every used beta ==0.')


def grid():
    tasks = []
    cells = O.cells()
    for block, variant, offset in ((1, 'large_rib', 1000), (2, 'nominal', 2000), (5, 'large_rib', 3000)):
        cond = 'breath_rib_fast' if block == 5 else 'breath_rib'
        indices = (0, 2, 6) if block == 5 else (1, 3, 4, 5, 7, 8)
        running = 0
        for dy, ph in CC.SETTINGS:
            for i in indices:
                found = [c for c in cells if c[0] == cond and c[2:] == (dy, ph, i)]
                assert len(found) == 1, (cond, dy, ph, i, found)
                case = dict(block=block, scene=f'{dy}_{int(ph)}_i{i}', cell=found[0], variant=variant,
                            seed=CC.SENSOR_BASE+offset+running, index=None, eps=0., sign=1)
                for sign in ((1, -1) if block == 5 else (1,)):
                    for arm in ARMS:
                        tasks.append(dict(case, sign=sign, arm=arm))
                running += 1
    for index in range(36, 42):
        case = CC.cases()[index]
        common = dict(scene=f'case_{index}', index=index, cell=case['cell'], variant=case['variant'])
        for offset in (10000, 20000):
            for arm in ARMS:
                tasks.append(dict(common, block=3, seed=case['seed']+offset, eps=0., sign=1, arm=arm))
        for eps in (.005, -.005):
            for arm in ARMS:
                tasks.append(dict(common, block=4, seed=case['seed'], eps=eps, sign=1, arm=arm))
    for number, task in enumerate(tasks):
        task['run_id'] = number
    assert len(tasks) == 120
    return tasks


def shift_event(sc, flip):
    """Infer sign from curvature of stored u; remove original and shift to 6 s."""
    assert sc['source_cond'] == 'shift_fast' and sc['regime'] == 'nominal'
    assert len(sc['event_info']) == len(sc['events']) == len(sc['motion_components']) == 1
    ev = dict(sc['event_info'][0])
    assert ev['kind'] == 'fast' and ev['start'] == .2*sc['T']
    assert not np.any(sc['act']) and not np.any(sc['fext'])
    t = np.arange(len(sc['u']))*.001
    amp, start = ev['amplitude'], ev['start']
    dur, hold = 1.875*amp/ev['peak_limit'], .15
    assert abs(ev['end']-(start+2*dur+hold)) < 1e-12
    def wave(at):
        return amp*(B.S3._smooth((t-at)/dur)-B.S3._smooth((t-at-dur-hold)/dur))
    original = sc['u'].copy()
    w = wave(start)
    interval = (t[1:-1] >= start) & (t[1:-1] <= ev['end'])
    curvature = np.diff(w, n=2)[interval]
    sign = 1 if np.dot(np.diff(original, n=2)[interval], curvature) > 0 else -1
    w *= sign
    component_error = float(np.max(np.abs(w-sc['motion_components'][0])))
    assert component_error < 1e-14
    residual = original-w
    second = np.abs(np.diff(residual, n=2))
    elsewhere = ~((t[1:-1] >= start-.002) & (t[1:-1] <= ev['end']+.002))
    check = dict(original_sign=sign, component_max_abs_difference=component_error,
                 original_start=start, original_end=ev['end'],
                 residual_max_second_difference_interval=float(second[interval].max()),
                 residual_max_second_difference_elsewhere=float(second[elsewhere].max()))
    before = dict(sc)
    moved = sign*flip*wave(6.)
    sc['u'] = residual+moved
    ev.update(start=6., end=6.+2*dur+hold)
    assert ev['end']+.6 < sc['T']
    sc['event_info'] = [ev]
    sc['events'] = [('fast', ev['start'], ev['end'])]
    sc['motion_components'] = [moved]
    changed = {'u', 'event_info', 'events', 'motion_components'}
    assert all(O.equal(before[k], sc[k]) for k in sc if k not in changed)
    check.update(new_start=ev['start'], new_end=ev['end'], unchanged_other_fields=True, act_zero=True, fext_zero=True)
    return check


def job(task, output):
    started = time.perf_counter()
    O.setup()
    sc = SE.MC.scene_for(tuple(task['cell']), CC.BASE, task['variant'])
    assert sc['cond'] == task['cell'][0] and sc['T'] == 10.
    check = shift_event(sc, task['sign']) if task['block'] == 5 else None
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    recorder = TL.ReturnedEstimate(pol.estimator.estimate, task['arm'], task['index'])
    def returned(o, n, policy):
        k, sigma, chi = recorder(o, n, policy)
        k *= 1.+task['eps']
        recorder.k[-1] = k
        return k, sigma, chi
    pol.estimator.estimate = returned
    log = P3.episode(sc, pol, task['seed'])
    if log is None:
        raise RuntimeError(f'No episode: {task}')
    metric = M.metric(log, pol, sc)
    row = dict(task, L=float(np.log1p(metric['J'])),
               **{k:float(metric[k]) for k in ('rmse', 'peak', 'over_15', 'loss_time')},
               FD=float(sc['FD']), speed_mm_s=float(sc['v']*1000), builder_accepted=True,
               cell_exists=tuple(task['cell']) in O.cells(), event_info=sc['event_info'], filter_resets=recorder.filter.resets)
    te = np.asarray(recorder.t)
    ids = np.rint(te/TL.V3.TC).astype(int)
    assert np.max(np.abs(log[ids, 0]-te), initial=0.) < 1e-9
    if task['arm'] in ('delag', 'delag2'):
        beta = np.asarray(recorder.beta)
        assert np.all(np.isfinite(beta))
        row['beta'] = dict(median_after_3=float(np.median(beta[te > 3.])), max=float(beta.max()))
        if task['arm'] == 'delag2':
            row['beta_zero_identity'] = TL.beta_zero_identity(recorder)
    if task['block'] == 5:
        ev = sc['event_info'][0]
        kr, kt = np.asarray(recorder.k), log[ids, 26]
        valid = np.isfinite(kr) & (kr > 0) & (kt > 0)
        e = np.full(len(te), np.nan)
        e[valid] = np.log(kr[valid]/kt[valid])
        em = (te >= ev['start']) & (te <= ev['end']+.6) & (te > 1.5) & (kt > 0) & (log[ids, 18] > .5)
        good = em & np.isfinite(e)
        fm = (log[:, 0] >= ev['start']) & (log[:, 0] <= ev['end']+.6)
        before = (log[:, 0] >= 5.9) & (log[:, 0] < 6.)
        onset = (log[:, 0] >= 6.) & (log[:, 0] < 6.12)
        row['event'] = dict(mean_log_error=float(e[good].mean()) if good.any() else None,
                            min_log_error=float(e[good].min()) if good.any() else None,
                            n_scored=int(em.sum()), n_missing=int((em & ~good).sum()),
                            force_rms=D.rms(log[fm, 18]-sc['FD']),
                            onset_force_change=float(log[onset, 18].mean()-log[before, 18].mean()),
                            beta_before=float(beta[te < 6.][-1]) if task['arm'] in ('delag', 'delag2') else None)
        row['event_shift_check'] = check
    if task['block'] in (1, 5) or task.get('save_trace', False):
        filename = f"run_{task['run_id']:03d}_{task['arm']}.npz"
        arrays = dict(t=log[:, 0], fax=log[:, 18], k_eff=log[:, 26], y=log[:, 2], t_est=te,
                      k_ret=np.asarray(recorder.k), n_est=np.asarray(recorder.n), trace=np.asarray(pol.trace))
        if task['arm'] in ('delag', 'delag2'):
            arrays.update(k_w=np.asarray(recorder.kw), beta=np.asarray(recorder.beta),
                          force_ratio=np.asarray(recorder.force_ratio))
        with (output/filename).open('xb') as stream:
            np.savez_compressed(stream, **arrays)
        row['file'] = filename
    row['runtime_s'] = time.perf_counter()-started
    return row


def subprocess_job(task, output, spool):
    result, stdout = spool/f"{task['run_id']}.json", spool/f"{task['run_id']}.stdout"
    spec = spool/f"{task['run_id']}.task.json"
    spec.write_text(json.dumps(task), encoding='utf-8')
    cmd = [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()), '--project-root', str(PROJECT),
           '--worker-task', str(spec), '--worker-output', str(output), '--worker-result', str(result)]
    with stdout.open('x', encoding='utf-8') as stream:
        done = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=stream, stderr=stream,
                              creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
    if done.returncode:
        raise RuntimeError(stdout.read_text(encoding='utf-8'))
    return json.loads(result.read_text(encoding='utf-8'))


def pairs(rows, block):
    selected = [r for r in rows if r['block'] == block]
    keys = sorted({(r['scene'], r['seed'], r['eps'], r['sign']) for r in selected})
    return [{a:next(r for r in selected if (r['scene'], r['seed'], r['eps'], r['sign']) == k and r['arm'] == a)
             for a in ARMS} for k in keys]


def summarize(rows, output, round2):
    fresh = [r for r in rows if r['block'] == 1]
    summary = {a:dict(mean_L=float(np.mean([r['L'] for r in fresh if r['arm'] == a])),
                      median_L=float(np.median([r['L'] for r in fresh if r['arm'] == a]))) for a in ARMS}
    for arm in ARMS:
        errors, masks, forces = [], [], []
        for row in fresh:
            if row['arm'] != arm or row['speed_mm_s'] <= 15.:
                continue
            e, mm, _, force = TL.scored_arrays(row, output)
            errors.append(e)
            masks.append(mm)
            forces.append(force['fast'])
        e = np.concatenate(errors)
        summary[arm]['fast_scenes'] = [r['scene'] for r in fresh if r['arm'] == arm and r['speed_mm_s'] > 15.]
        summary[arm]['identification'] = {s:TL.error_stats(e, np.concatenate([m[s] for m in masks])) for s in ('rise', 'fall', 'reversal')}
        f = np.concatenate(forces)
        summary[arm]['crossing_force_rms'] = dict(n=len(f), rms=D.rms(f))
    for pair in pairs(rows, 5):
        direction = 'toward probe' if pair['huber']['event']['onset_force_change'] > 0 else 'away from probe'
        for r in pair.values():
            r['direction'] = direction
    b5 = {}
    for direction in ('all', 'toward probe', 'away from probe'):
        pp = [p for p in pairs(rows, 5) if direction == 'all' or p['huber']['direction'] == direction]
        b5[direction] = dict(n=len(pp), mean_L_delta=float(np.mean([p['delag']['L']-p['huber']['L'] for p in pp])),
            mean_event_force_rms_delta=float(np.mean([p['delag']['event']['force_rms']-p['huber']['event']['force_rms'] for p in pp])),
            max_peak=max(r['peak'] for p in pp for r in p.values()))
    ref = {(r['index'], r['arm']):r for r in round2['rows']}
    zero = [dict(ref[i, a], block=4, scene=f'case_{i}', eps=0., sign=1) for i in range(36, 42) for a in ARMS]
    perturb = [r for r in rows if r['block'] == 4]+zero
    means = {str(eps):{a:float(np.mean([r['L'] for r in perturb if r['eps'] == eps and r['arm'] == a])) for a in ARMS} for eps in (0., .005, -.005)}
    spreads = {str(i):dict(spread=float(np.ptp([r['L'] for r in perturb if r['index'] == i and r['arm'] == 'delag'])),
                           round2_gap=ref[i, 'huber']['L']-ref[i, 'delag']['L']) for i in (36, 37)}
    d1 = [p['delag']['L']-p['huber']['L'] for p in pairs(rows, 1)]
    d2 = [p['delag']['L']-p['huber']['L'] for p in pairs(rows, 3)]
    d5 = [abs(p['delag']['L']-p['huber']['L']) for p in pairs(rows, 2)]
    beta_max = max(r['beta']['max'] for r in rows if r['block'] == 2 and r['arm'] == 'delag')
    verdicts = dict(
        R1=dict(met=bool(np.mean(d1) < 0 and max(d1) <= .02 and sum(d <= .005 for d in d1) >= 9), mean_delta=float(np.mean(d1)), worst_delta=max(d1), improved_or_equal=sum(d <= .005 for d in d1)),
        R2=dict(met=bool(max(d2) <= .02 and np.mean(d2) < 0), mean_delta=float(np.mean(d2)), worst_delta=max(d2)),
        R3=dict(met=bool(all(m['delag'] < m['huber'] for m in means.values()) and all(s['spread'] < s['round2_gap'] for s in spreads.values())), means=means, spreads=spreads),
        R4=dict(met=bool(b5['all']['mean_event_force_rms_delta'] > 0 and b5['all']['mean_L_delta'] <= .05 and b5['all']['max_peak'] <= 15.), expected_weakness_met=b5['all']['mean_event_force_rms_delta'] > 0,
                bound_met=bool(b5['all']['mean_L_delta'] <= .05 and b5['all']['max_peak'] <= 15.), by_direction=b5),
        R5=dict(met=bool(max(d5) <= .01 and beta_max == 0.), max_abs_delta=max(d5), max_beta=beta_max))
    return dict(block1=summary, block5=b5), verdicts, [dict(index=i, **{a:ref[i, a]['L'] for a in ARMS}) for i in range(36, 42)]


def source_audit():
    source = inspect.getsource(B.S3.cell_scene3)
    required = ["start = .2*T+j*.17*T", "if kind in ('slow', 'fast'):", "act = np.maximum(act, level*np.sin", "components.append(wave)", "events.append((kind, start, end))", "amplitude=amp, peak_limit=peak"]
    assert all(s in source for s in required)
    return dict(sha256=hashlib.sha256(source.encode('utf-8')).hexdigest(),
                time_related_lines=[line.strip() for line in source.splitlines() if any(s in line for s in ('start', 'end', 'act =', 'components.append', 'events.append', 'event_info.append'))],
                finding='fast uses displacement only; activation belongs to else (flinch/cough). External-drive conversion inactive in nominal regime. Only u/events/event_info/motion_components depend on this fast event time; all other fields checked unchanged.')


ROUND4_PREDICTIONS = dict(
    S1='Block 2: every abs(delag2 L - huber L) <= .01, including V4.',
    S2='Original A_large mean within .01 of delag; block 1 mean within .02 of delag; every A_large/block1 delag2-huber delta <= .02.',
    S3='Original cases 0-23: every abs(delag2 L - huber L) <= .005.',
    S4='Block 5: mean delag2-huber L delta <= .02; away-from-probe separately.')


def diagnose_v4(rows, output):
    h, d = (next(r for r in rows if r['arm'] == a) for a in ('huber', 'delag'))
    with np.load(output/h['file']) as z:
        t, hf = z['t'], z['fax']
    with np.load(output/d['file']) as z:
        assert np.array_equal(t, z['t'])
        df, te, ratio, beta = z['fax'], z['t_est'], z['force_ratio'], z['beta']
    first = int(np.flatnonzero(np.abs(df-hf) > .2)[0])
    at = float(t[first])
    stats = {}
    for label, mask in (('before', te < at), ('after', te >= at)):
        values = ratio[mask & np.isfinite(ratio)]
        stats[label] = dict(n=len(values), minimum=float(values.min()), maximum=float(values.max()),
                            rms=D.rms(values), rms_deviation_from_one=D.rms(values-1.))
    error = df[first:]-d['FD']
    # Demean and Hann-window the complete post-divergence force error; exclude DC.
    spectrum = np.abs(np.fft.rfft((error-error.mean())*np.hanning(len(error))))**2
    frequency = np.fft.rfftfreq(len(error), d=float(t[1]-t[0]))
    dominant = int(np.argmax(spectrum[1:])+1)
    return dict(first_force_difference_over_02_s=at, ratio=stats,
                dominant_force_error_frequency_hz=float(frequency[dominant]),
                frequency_resolution_hz=float(frequency[1]),
                spectrum_convention='fax-FD after first >0.2 N divergence, demeaned Hann periodogram, largest non-DC bin',
                beta_zero_throughout=bool(np.all(beta == 0.)),
                peak={r['arm']:r['peak'] for r in rows}, L={r['arm']:r['L'] for r in rows},
                rows=rows)


def summarize_round4(original, rows, output, diagnosis):
    arms = ('huber', 'delag', 'delag2')
    lookup = {(r['index'], r['arm']):r for r in original}
    groups = {g:{a:dict(mean_L=float(np.mean([r['L'] for r in original if r['group'] == g and r['arm'] == a])),
                        median_L=float(np.median([r['L'] for r in original if r['group'] == g and r['arm'] == a])))
                 for a in arms} for g in dict.fromkeys(r['group'] for r in original)}
    def paired(block):
        selected = [r for r in rows if r['block'] == block]
        keys = sorted({(r['scene'], r['seed'], r['eps'], r['sign']) for r in selected})
        return [{a:next(r for r in selected if (r['scene'], r['seed'], r['eps'], r['sign']) == key and r['arm'] == a)
                 for a in arms} for key in keys]
    b1 = {}
    for arm in arms:
        selected = [r for r in rows if r['block'] == 1 and r['arm'] == arm]
        ee, mm, ff = [], [], []
        for row in selected:
            if row['speed_mm_s'] <= 15.:
                continue
            e, masks, _, force = TL.scored_arrays(row, output)
            ee.append(e)
            mm.append(masks)
            ff.append(force['fast'])
        b1[arm] = dict(mean_L=float(np.mean([r['L'] for r in selected])),
                       median_L=float(np.median([r['L'] for r in selected])),
                       identification={s:TL.error_stats(np.concatenate(ee), np.concatenate([m[s] for m in mm]))
                                       for s in ('rise', 'fall', 'reversal')},
                       crossing_force_rms=D.rms(np.concatenate(ff)))
    b5 = {}
    for pair in paired(5):
        for r in pair.values():
            r['direction'] = pair['huber']['direction']
    for direction in ('all', 'toward probe', 'away from probe'):
        pp = [p for p in paired(5) if direction == 'all' or p['huber']['direction'] == direction]
        b5[direction] = dict(n=len(pp),
            mean_L={a:float(np.mean([p[a]['L'] for p in pp])) for a in arms},
            mean_L_delta=float(np.mean([p['delag2']['L']-p['huber']['L'] for p in pp])),
            mean_event_force_rms={a:float(np.mean([p[a]['event']['force_rms'] for p in pp])) for a in arms},
            max_peak={a:max(p[a]['peak'] for p in pp) for a in arms})
    d1 = [p['delag2']['L']-p['huber']['L'] for p in paired(1)]
    d2 = [abs(p['delag2']['L']-p['huber']['L']) for p in paired(2)]
    d0 = [abs(lookup[i, 'delag2']['L']-lookup[i, 'huber']['L']) for i in range(24)]
    dl = [lookup[i, 'delag2']['L']-lookup[i, 'huber']['L'] for i in range(36,42)]
    gap_large = groups['A_large']['delag2']['mean_L']-groups['A_large']['delag']['mean_L']
    gap_b1 = b1['delag2']['mean_L']-b1['delag']['mean_L']
    predictions = dict(
        S1=dict(met=bool(max(d2) <= .01), max_abs_delta=max(d2)),
        S2=dict(met=bool(abs(gap_large) <= .01 and abs(gap_b1) <= .02 and max(d1+dl) <= .02),
                original_A_large_mean_gap=gap_large, block1_mean_gap=gap_b1,
                worst_original_A_large_delta=max(dl), worst_block1_delta=max(d1)),
        S3=dict(met=bool(max(d0) <= .005), max_abs_delta=max(d0), worst_case=int(np.argmax(d0))),
        S4=dict(met=bool(b5['all']['mean_L_delta'] <= .02), by_direction=b5))
    return dict(original_groups=groups, block1=b1, block5=b5, V4=diagnosis), predictions


def run_round4(base):
    loop_path, robust_path = base/'trajectory_loop.json', base/'trajectory_robust.json'
    loop = json.loads(loop_path.read_text(encoding='utf-8'))
    robust = json.loads(robust_path.read_text(encoding='utf-8'))
    if any(r['arm'] == 'delag2' for r in loop['rows']+robust['rows']):
        raise FileExistsError('Refusing duplicate round 4 rows')
    # Freeze predictions in the result before any diagnosis or candidate episodes.
    robust['round4'] = dict(status='running', prespecified_predictions=ROUND4_PREDICTIONS)
    robust_path.write_text(json.dumps(robust, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Round 4 predictions BEFORE RUN '+json.dumps(ROUND4_PREDICTIONS), flush=True)
    spool = Path(tempfile.mkdtemp(prefix='trajectory_round4_workers_'))
    output, loop_output = base/'trajectory_robust', base/'trajectory_loop'
    started = time.perf_counter()
    # A matched Huber rerun is required because block 2 had no saved force trace.
    diagnosis_rows = []
    with ThreadPoolExecutor(8) as pool:
        pending = []
        for i, arm in enumerate(('huber', 'delag')):
            stored = next(r for r in robust['rows'] if r['block'] == 2 and r['scene'] == 'viscous_14_i4' and r['arm'] == arm)
            task = {k:stored[k] for k in ('block','scene','cell','variant','seed','index','eps','sign','arm')}
            task.update(run_id=198+i, save_trace=True)
            pending.append(pool.submit(subprocess_job, task, output, spool))
        for future in as_completed(pending):
            row = future.result()
            stored = next(r for r in robust['rows'] if r['block'] == 2 and r['scene'] == row['scene'] and r['arm'] == row['arm'])
            assert all(row[k] == stored[k] for k in ('L','rmse','peak','over_15','loss_time'))
            diagnosis_rows.append(row)
    diagnosis = diagnose_v4(diagnosis_rows, output)
    print('V4 diagnosis '+json.dumps(diagnosis), flush=True)
    tasks = [dict(t, arm='delag2') for t in grid() if t['block'] in (1,2,3,5) and t['arm'] == 'huber']
    for i, task in enumerate(tasks):
        task['run_id'] = 120+i
        if task['block'] == 2 and task['scene'] == 'viscous_14_i4':
            task['save_trace'] = True
    assert len(tasks) == 48
    new_loop, new_robust = [], []
    with ThreadPoolExecutor(8) as pool:
        pending = {pool.submit(subprocess_job, t, output, spool):'robust' for t in tasks}
        pending.update({pool.submit(TL.subprocess_job, ('delag2', i, str(loop_output)), spool):'loop' for i in TL.CASES})
        for number, future in enumerate(as_completed(pending), 1):
            row = future.result()
            (new_loop if pending[future] == 'loop' else new_robust).append(row)
            print(f"{number}/78 {pending[future]} {row.get('scene', row['index'])} L={row['L']:.9f}", flush=True)
    new_loop.sort(key=lambda r:r['index'])
    new_robust.sort(key=lambda r:r['run_id'])
    # Score only new rows; retain all existing evidence verbatim as JSON objects.
    groups, identification = TL.summarize(new_loop, loop_output)
    for name, additions in groups.items():
        loop['groups'][name].update(additions)
    for name, additions in identification.items():
        loop['identification'][name].update(additions)
    loop['rows'].extend(new_loop)
    for c in loop['comparison']:
        r = next(r for r in new_loop if r['index'] == c['index'])
        c['L']['delag2'], c['peak']['delag2'] = r['L'], r['peak']
    robust['rows'].extend(new_robust)
    original_robust = [dict(r, block=0, scene=f"case_{r['index']}", eps=0., sign=1, run_id=168+i,
                            file='../../trajectory_loop/'+r['file']) for i,r in enumerate(new_loop)]
    # Relative to the robust trace directory, original traces live in its sibling.
    for r in original_robust:
        r['file'] = '../trajectory_loop/'+Path(r['file']).name
    robust['rows'].extend(original_robust)
    summary, predictions = summarize_round4(loop['rows'], robust['rows'], output, diagnosis)
    for obj, path, new_rows in ((loop, loop_path, new_loop), (robust, robust_path, new_robust+original_robust)):
        stored = json.loads((PROJECT/'experiments/rossi_failure_study/results'/path.name).read_text(encoding='utf-8'))
        assert obj['rows'][:len(stored['rows'])] == stored['rows']
        zero_checks = [r['beta_zero_identity'] for r in new_rows]
        obj['round4'] = dict(status='complete', workers=8, new_episodes=len(new_rows),
            diagnosis_episodes=2 if path == robust_path else 0,
            prespecified_predictions=ROUND4_PREDICTIONS, predictions=predictions,
            summary=summary, runtime_s=time.perf_counter()-started,
            beta_zero_identity=dict(n=sum(c['n'] for c in zero_checks), mismatches=sum(c['mismatches'] for c in zero_checks)),
            conventions=dict(correction='delag2: k_c=clip(k_w/r,200,3e5), r=clip(1+beta*x/Cw,.5,2); zero beta returns original k_w exactly; gain update unchanged.',
                             scoring='Unchanged round2 masks and round3 pooling/event/direction conventions.',
                             original_rows='Block 0 in robust JSON; primary traces in trajectory_loop; huber/delag original reference rows in trajectory_loop.json.',
                             provenance='All pre-round4 rows retained unchanged; stored huber/delag reused; V4 exact matched reruns saved separately in summary.'))
        path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Round 4 verdicts '+json.dumps(predictions), flush=True)
    print(f'Completed 78 delag2 episodes + 2 V4 matched reruns in {time.perf_counter()-started:.1f}s', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=PROJECT)
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--round4', action='store_true')
    parser.add_argument('--worker-task', type=Path)
    parser.add_argument('--worker-output', type=Path)
    parser.add_argument('--worker-result', type=Path)
    args = parser.parse_args()
    if args.worker_task:
        row = job(json.loads(args.worker_task.read_text(encoding='utf-8')), args.worker_output)
        with args.worker_result.open('x', encoding='utf-8') as stream:
            json.dump(row, stream, allow_nan=False)
        return
    base = (args.output_root or PROJECT)/'experiments/rossi_failure_study/results'
    if args.round4:
        run_round4(base)
        return
    output, destination = base/'trajectory_robust', base/'trajectory_robust.json'
    if destination.exists() or (output.exists() and any(output.iterdir())):
        raise FileExistsError('Refusing to overwrite round 3 outputs')
    output.mkdir(parents=True, exist_ok=True)
    tasks, audit = grid(), source_audit()
    round2 = json.loads((PROJECT/'experiments/rossi_failure_study/results/trajectory_loop.json').read_text(encoding='utf-8'))
    assert round2['round2']['status'] == 'complete'
    print('Predictions BEFORE RUN '+json.dumps(PREDICTIONS), flush=True)
    validation = []
    for task in tasks:
        if task['arm'] != 'huber':
            continue
        sc = SE.MC.scene_for(tuple(task['cell']), CC.BASE, task['variant'])
        check = shift_event(sc, task['sign']) if task['block'] == 5 else None
        validation.append(dict(run_id=task['run_id'], block=task['block'], cell=task['cell'], builder_accepted=True, cell_exists=tuple(task['cell']) in O.cells(), shift=check))
    print(f'Preflight accepted all {len(validation)} paired specifications; all overlay cells exist.', flush=True)
    spool = Path(tempfile.mkdtemp(prefix='trajectory_robust_workers_'))
    print(f'Worker IPC: {spool}', flush=True)
    rows, started = [], time.perf_counter()
    with ThreadPoolExecutor(8) as pool:
        futures = [pool.submit(subprocess_job, task, output, spool) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"{len(rows)}/120 block={row['block']} {row['scene']} {row['arm']} seed={row['seed']} eps={row['eps']} sign={row['sign']} L={row['L']:.9f}", flush=True)
    rows.sort(key=lambda r:r['run_id'])
    summary, verdicts, reference = summarize(rows, output, round2)
    result = dict(status='complete', workers=8, new_episodes=len(rows), rows=rows, prespecified_predictions=PREDICTIONS,
                  predictions=verdicts, summary=summary, validation=validation, event_source_audit=audit, round2_reference=reference,
                  round2_sha256=hashlib.sha256((PROJECT/'experiments/rossi_failure_study/results/trajectory_loop.json').read_bytes()).hexdigest(), runtime_s=time.perf_counter()-started,
                  conventions=dict(seeds='Running numbers zero-based per fresh block; settings outermost.',
                    perturbation='Only stiffness multiplied after imported arm; sigma/chi unchanged; k_ret actual controller return.',
                    identification='Round2 own-arm masks/validity; finite returns pooled, missing counted. Fast iff sc[v]*1000>15 mm/s.',
                    crossing_force='Round2 fast mask rise|fall|reversal, pooled within fast scenes; RMS fax-FD.',
                    beta='Gain used before update; median all t>3 ticks; pre-event last tick t<6.',
                    event_window='[6,event end+.6]; error round2 validity; force RMS all control ticks.',
                    direction='Huber mean fax [6,6.12) minus [5.9,6); positive=toward; paired delag same label.',
                    R4='Weakness: positive paired mean event RMS delta; bound: full-episode mean L delta and max peak both arms.',
                    smoothness='Raw second difference at 1ms; original closed event interval vs elsewhere excluding 2ms margins. Sign from stored-u curvature correlation.'))
    with destination.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print('Verdicts '+json.dumps(verdicts), flush=True)
    print(f'Written {destination}; elapsed {result["runtime_s"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
