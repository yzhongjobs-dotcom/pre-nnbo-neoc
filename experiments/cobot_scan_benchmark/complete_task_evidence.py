"""Finish the focused task evidence without changing frozen controllers.

Commands: protocol, check, run, ablate, time, reference, summarize.
Main comparisons use fresh paired scenes. Isolation tasks transfer the selected
T1/T3 settings without retuning; they are explicitly not task-retuned rankings.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
import plant as P
import scenarios as S
import report_physical_metrics as R
import refine_baselines as F

OUT = ROOT / 'completed_tasks'
TASKS = ('T1_static', 'T2_motion', 'T3_combined')
SEED = 7200
N = 40


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(name, value):
    OUT.mkdir(exist_ok=True)
    (OUT / name).write_text(json.dumps(F.clean(value), ensure_ascii=False,
                                     indent=1, allow_nan=False), encoding='utf-8')


def specs(task):
    source = 'T1' if task == 'T1_static' else 'T3'
    old = R.specifications(source)
    names = ('Best_fixed', 'EtR-A', 'EtR-A-Haninger', 'EtR-R', 'Candidate1_v2', 'Candidate2')
    out = [(name, 'frozen', (name, kind, cfg, seed), seed)
           for name, kind, cfg, seed in old if name in names]
    selected = read(F.OUT / 'selection.json')
    out += [(name, 'refined', selected[source + '_' + name]['spec'], -1)
            for name in F.FAMILIES]
    learned = read(F.OUT / 'learning.json')
    out += [('Learned_VIC', 'refined', learned[f'{source}_{seed}']['spec'], seed)
            for seed in range(3)]
    return out


def build(spec):
    _, kind, config, _ = spec
    return R.build(config) if kind == 'frozen' else F.make(config)


def scenes(task):
    if task == 'T1_static':
        return S.cohort('T1', N, SEED, dict(ratio=30., A=0.))
    if task == 'T3_combined':
        return S.cohort('T3', N, SEED, R.CELLS['T3'])
    # Four balanced cells: soft/stiff x initial push/retract. Surface returns
    # after a 2 s hold, giving two nonoverlapping recovery windows per scene.
    rng = np.random.default_rng(SEED)
    out = []
    for i in range(N):
        stiff = (i % 4) >= 2
        E = rng.uniform(5e5, 1.5e6) if stiff else rng.uniform(1e4, 4e4)
        p = S._base(rng, E, 0., .015, 8.)
        p[15] = 1
        p[16:20] = [2., .005 * (1 if i % 2 == 0 else -1), .5, 2.]
        out.append((p, 8.))
    return out


def verify():
    F.verify_frozen(read(F.OUT / 'protocol.json'))
    protocol = read(OUT / 'protocol.json')
    assert protocol['code_sha256'] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def protocol():
    F.verify_frozen(read(F.OUT / 'protocol.json'))
    value = dict(seed=SEED, n_scenes=N, tasks=TASKS, source_settings={
        'T1_static': 'frozen T1 settings, breathing removed',
        'T2_motion': 'frozen T3 settings, no retraining or per-scene tuning',
        'T3_combined': 'unchanged original main-task dose, fresh scenes'},
        force_target_N=4., warmup_s=1., logging_period_s=.008,
        recovery=dict(band_N=[3.5, 4.5], dwell_s=.2,
                      windows_s=[[2.5, 4.5], [5., 7.]],
                      right_censor_s=2., definition='time from ramp end to start of first 0.2 s in band; censored if absent'),
        timing='single process after warmup; observation assembly through output projection; excludes simulated plant and physical communications',
        ablations='frozen inference-time interventions; no retraining; same fresh scenes',
        reference='C1 final local OCP snapshots; 3 starts x 25/50/100 GN iterations; not global optimality or C2 optimality',
        new_scene_status='protocol fixed before new simulations; no adaptive stopping or test-based retuning',
        python=sys.version, numpy=np.__version__, processor=platform.processor(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        frozen_protocol_sha256=hashlib.sha256((F.OUT / 'protocol.json').read_bytes()).hexdigest())
    if (OUT / 'protocol.json').exists():
        assert read(OUT / 'protocol.json') == value
    else:
        save('protocol.json', value)
    for task in TASKS:
        ss = scenes(task)
        np.savez_compressed(OUT / f'scenes_{task}.npz', p=np.array([p for p, _ in ss]),
                            T=np.array([T for _, T in ss]), noise_seed=SEED + np.arange(N))
    print('PROTOCOL frozen', flush=True)


def episode(p, pol, T, seed, update, timed=False):
    """Same advance calls/noise/controller inputs as plant.episode; add logging."""
    pl = P.PLATFORMS['PC125']
    dt = .0005
    nc = round(pl['Tc'] / dt)
    nt = round(T / pl['Tc'])
    log = np.zeros((nt + 2, P.LOGC))
    s = P.init_state(p, dt)
    noise = np.random.default_rng(seed).standard_normal((int(T * 1000) + 10, 2))
    args = (nc, pl['ndel'], pl['tau_a'], pl['fc'], pl['VMAX'], dt, noise)
    is_policy = callable(pol)
    if is_policy:
        pol.reset(pl)
    per = max(round((update or pl['Tc']) / pl['Tc']), 1) if is_policy else nt
    n = 0
    tt, commands = [], []
    while n < nt:
        start = time.perf_counter_ns() if timed else 0
        obs = P.Obs(log, n)
        g = np.clip(P._pad(pol(obs) if is_policy else pol), P.GMIN, P.GMAX)
        if timed:
            tt.append((time.perf_counter_ns() - start) * 1e-6)
        commands.append([n * pl['Tc'], *g])
        n = P.advance(s, p, g, min(per, nt - n), *args, log, n)
        if n < 0:
            return None, np.array(tt), np.array(commands)
    return log[:n], np.array(tt), np.array(commands)


def recover(log):
    out = {}
    for label, (a, b) in zip(('out', 'return'), ((2.5, 4.5), (5., 7.))):
        out[label + '_recovered'] = 0
        out[label + '_restricted_recovery_s'] = 2.
        if log is None:
            continue
        q = log[(log[:, 0] >= a) & (log[:, 0] < b)]
        good = np.abs(q[:, 3] - 4.) <= .5
        dwell = int(np.ceil(.2 / .008))
        for k in range(len(good) - dwell + 1):
            if good[k:k + dwell].all():
                out[label + '_recovered'] = 1
                out[label + '_restricted_recovery_s'] = max(0., float(q[k, 0] - a))
                break
    return out


def metrics(log, T, task):
    out = R.physical(log, T)
    out.update(velocity_rms_mm_s=None, velocity_peak_mm_s=None)
    if log is not None:
        v = np.diff(log[:, 2]) / np.diff(log[:, 0]) * 1000
        v = v[log[1:, 0] >= 1.]
        out.update(velocity_rms_mm_s=float(np.sqrt(np.mean(v * v))),
                   velocity_peak_mm_s=float(np.max(np.abs(v))))
    if task == 'T2_motion':
        out.update(recover(log))
    return out


def bind(pol, p):
    if hasattr(pol, 'bind'):
        pol.bind(p)


def warm(spec):
    p, _ = S.cohort('T3', 1, 1000, R.CELLS['T3'])[0]
    pol, update = build(spec)
    bind(pol, p)
    P.episode(p, pol, T=.3, seed=1000, update=update)


def job(value):
    task, spec, indices, keep = value
    warm(spec)
    ss = scenes(task)
    rows = []
    for i in indices:
        p, T = ss[i]
        pol, update = build(spec)
        bind(pol, p)
        log, _, commands = episode(p, pol, T, SEED + i, update)
        row = dict(task=task, method=spec[0], training_seed=spec[3], scene=i,
                   **metrics(log, T, task))
        rows.append(row)
        if keep and i == 0:
            path = OUT / 'traces'
            path.mkdir(exist_ok=True)
            np.savez_compressed(path / f'{task}_{spec[0]}_{spec[3]}.npz',
                                log=np.zeros((0, P.LOGC)) if log is None else log,
                                failed=log is None, commands=commands, p=p)
    return rows


def csvwrite(name, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    R.write_csv(OUT / name, [{k: r.get(k) for k in fields} for r in rows])


def checks():
    verify()
    rows = []
    for task in TASKS:
        p, _ = scenes(task)[0]
        for spec in specs(task):
            warm(spec)
            a, update = build(spec)
            bind(a, p)
            original = P.episode(p, a, T=.6, seed=SEED, update=update)
            b, update = build(spec)
            bind(b, p)
            measured, _, _ = episode(p, b, .6, SEED, update)
            assert np.array_equal(original, measured), (task, spec[0], spec[3])
            rows.append(dict(task=task, method=spec[0], training_seed=spec[3], max_difference=0.))
    # Synthetic force traces check both recovery and right censoring.
    t = np.arange(1000) * .008
    lg = np.zeros((1000, P.LOGC)); lg[:, 0] = t; lg[:, 3] = 4.
    assert recover(lg)['out_recovered'] == 1
    lg[:, 3] = 0.
    assert recover(lg)['out_recovered'] == 0 and recover(lg)['out_restricted_recovery_s'] == 2.
    save('checks.json', dict(passed=True, exact_replays=rows, recovery_checks=True))
    print('CHECKS passed', len(rows), flush=True)


def run_jobs(jobs, filename):
    verify()
    assert read(OUT / 'checks.json')['passed']
    result = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(job, value) for value in jobs]
        for index, future in enumerate(as_completed(futures)):
            result += future.result()
            print(filename, index + 1, '/', len(jobs), 'episodes', len(result), flush=True)
    verify()
    csvwrite(filename, sorted(result, key=lambda r: (r['task'], r['method'], r['training_seed'], r['scene'])))


def main_runs():
    jobs = [(task, spec, list(range(i, min(i + 10, N))), True)
            for task in TASKS for spec in specs(task) for i in range(0, N, 10)]
    run_jobs(jobs, 'per_scene.csv')


def ablation_specs(task):
    out = []
    for spec in specs(task):
        if spec[0] == 'Candidate1_v2':
            for name, changes in (
                ('C1_no_surface_estimate', dict(drift=0)),
                ('C1_gains_only', dict(ffmode=0)),
                ('C1_information_delay100ms', dict(klat=.1)),
                ('C1_update40ms', dict(update=.04))):
                config = list(spec[2]); config[2] = {**config[2], **changes}
                out.append((name, 'frozen', config, -1))
        elif spec[0] == 'Candidate2':
            # Runtime interventions in exactly the trained full policy.
            for name, changes in (('C2_gains_only', dict(gains_only=True)),
                                  ('C2_no_force_features', dict(noforce=True)),
                                  ('C2_update96ms', dict(every=12))):
                out.append((name, 'c2_intervention', changes, spec[3]))
    return out


_original_build = build


def build(spec):
    if spec[1] != 'c2_intervention':
        return _original_build(spec)
    import candidate2_learned as C
    net, meta = C.load_net('full', spec[3])
    return C.Cand2Policy(net, **{**meta['cfg'], **spec[2]}), None


def ablate():
    jobs = [(task, spec, list(range(i, min(i + 10, N))), True)
            for task in ('T1_static', 'T3_combined') for spec in ablation_specs(task)
            for i in range(0, N, 10)]
    run_jobs(jobs, 'ablations.csv')


def timing():
    verify()
    rows = []
    raw = {}
    for task in ('T1_static', 'T3_combined'):
        for spec in specs(task):
            warm(spec)
            times = []
            for i in range(3):
                p, T = scenes(task)[i]
                pol, update = build(spec); bind(pol, p)
                _, tt, _ = episode(p, pol, T, SEED + i, update, timed=True)
                times.extend(tt)
            a = np.array(times)
            raw[f'{task}_{spec[0]}_{spec[3]}'] = a
            rows.append(dict(task=task, method=spec[0], training_seed=spec[3],
                             n_calls=len(a), median_ms=float(np.median(a)),
                             p99_ms=float(np.percentile(a, 99)), max_ms=float(a.max()),
                             over8ms=int((a > 8.).sum())))
            print('TIME', task, spec[0], spec[3], rows[-1]['p99_ms'], flush=True)
    csvwrite('timing.csv', rows)
    np.savez_compressed(OUT / 'timing_raw.npz', **raw)


def reference():
    verify()
    import candidate1_v2 as V
    real = V._solve_ord
    rows = []
    for task in ('T1_static', 'T3_combined'):
        spec = next(s for s in specs(task) if s[0] == 'Candidate1_v2')
        warm(spec)
        for scene, (p, T) in enumerate(scenes(task)[:6]):
            snapshots = []
            count = [0]
            def capture(*args):
                cost = real(*args)
                count[0] += 1
                if count[0] >= 125 and count[0] % 50 == 0:
                    snapshots.append((float(cost), [x.copy() if isinstance(x, np.ndarray) else x for x in args]))
                return cost
            V._solve_ord = capture
            try:
                pol, update = build(spec); bind(pol, p)
                P.episode(p, pol, T=T, seed=SEED + scene, update=update)
            finally:
                V._solve_ord = real
            for index, (candidate_cost, context) in enumerate(snapshots):
                lo, hi = context[12:14]
                starts = [context[0], np.clip(context[9], lo, hi), np.clip(context[10], lo, hi)]
                best = {}
                for budget in (25, 50, 100):
                    vals = []
                    for start in starts:
                        args = [x.copy() if isinstance(x, np.ndarray) else x for x in context]
                        args[0] = start.copy(); args[1] = budget
                        t0 = time.perf_counter_ns()
                        cost = real(*args)
                        vals.append((float(cost), (time.perf_counter_ns() - t0) * 1e-6))
                    best[budget] = min(x[0] for x in vals)
                ref = best[100]
                rows.append(dict(task=task, scene=scene, snapshot=index, candidate_cost=candidate_cost,
                                 reference25=best[25], reference50=best[50], reference100=ref,
                                 gap_percent=100 * (candidate_cost - ref) / max(abs(ref), 1e-12),
                                 convergence50_100_percent=100 * abs(best[50] - ref) / max(abs(ref), 1e-12)))
            print('REF', task, scene, len(snapshots), flush=True)
    csvwrite('reference.csv', rows)


def summarize():
    verify()
    import pandas as pd
    df = pd.read_csv(OUT / 'per_scene.csv')
    assert len(df) == 3 * 17 * N
    assert not df.duplicated(['task', 'method', 'training_seed', 'scene']).any()
    columns = [c for c in df if c not in ('task', 'method', 'training_seed', 'scene')]
    scene = df.groupby(['task', 'method', 'scene'])[columns].mean().reset_index()
    summary = scene.groupby(['task', 'method'])[columns].agg(['mean', 'std'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary = summary.reset_index()
    worst = df.groupby(['task', 'method']).peak_force_N.max().rename('worst_peak_force_N')
    summary = summary.merge(worst, on=['task', 'method'])
    failed = df.groupby(['task', 'method']).failed.sum().rename('failed_runs')
    summary = summary.merge(failed, on=['task', 'method'])
    # Never silently omit failures from finite-outcome group means.
    for ix in summary.index[summary.failed_runs > 0]:
        for c in summary:
            if c.endswith(('_mean', '_std')) and not c.startswith('failed'):
                summary.loc[ix, c] = np.nan
    summary.to_csv(OUT / 'summary.csv', index=False)
    # Simultaneous scene-paired bootstrap CIs for the two candidates versus all
    # eleven baselines in each task (66 comparisons; Bonferroni family level).
    rng = np.random.default_rng(20261022)
    pairs = []
    alpha = .05 / 66
    draws = rng.integers(0, N, size=(100000, N))
    for task in TASKS:
        wide = scene[scene.task == task].pivot(index='scene', columns='method', values='rmse_N')
        for candidate in ('Candidate1_v2', 'Candidate2'):
            for baseline in wide.columns:
                if baseline in ('Candidate1_v2', 'Candidate2'):
                    continue
                delta = (wide[candidate] - wide[baseline]).to_numpy()
                boot = delta[draws].mean(axis=1)
                lo, hi = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
                pairs.append(dict(task=task, candidate=candidate, baseline=baseline,
                                  delta_rmse_N=float(delta.mean()), simultaneous_low_N=float(lo),
                                  simultaneous_high_N=float(hi), family_comparisons=66))
    csvwrite('paired.csv', pairs)
    if (OUT / 'ablations.csv').exists():
        ab = pd.read_csv(OUT / 'ablations.csv')
        ab = ab.groupby(['task', 'method', 'scene']).mean(numeric_only=True).reset_index()
        ab.groupby(['task', 'method']).mean(numeric_only=True).to_csv(OUT / 'ablation_summary.csv')
    print(summary[['task', 'method', 'rmse_N_mean', 'peak_force_N_mean', 'no_contact_s_mean', 'failed_runs']].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['protocol', 'check', 'run', 'ablate', 'time', 'reference', 'summarize'])
    mode = parser.parse_args().mode
    {'protocol': protocol, 'check': checks, 'run': main_runs, 'ablate': ablate,
     'time': timing, 'reference': reference, 'summarize': summarize}[mode]()
