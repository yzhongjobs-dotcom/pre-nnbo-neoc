"""Fixed-seed H increment ablation. Writes only results/h_fix_eval.json.

python -X utf8 experiments/rossi_failure_study/h_fix_eval.py --workers 8
Phase 2 runs only after an explicit phase-1 effect check; no seed search.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import ast
import hashlib
import inspect
import json
import math
import textwrap
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import MethodType, SimpleNamespace
from pathlib import Path

import low_estimate_forensics as LF
import rib_transition_eval as R
import memory_segment_eval as SE
import contact_candidates as C
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'results/h_fix_eval.json'
CONFIGS = ('on_none', 'on_hfix', 'off_none', 'off_hfix',
           'off_hfix_forget', 'off_hfix_ewa')


def save(data):
    OUT.write_text(json.dumps(LF.clean(data), ensure_ascii=False,
                              allow_nan=False, separators=(',', ':')), encoding='utf-8')


def numerical_checks():
    # Deterministic algebraic fixtures, not simulation seeds.
    y = .018 + np.cumsum(np.r_[0., np.sin(np.arange(1, 1000)*.71)*6e-6])
    log = np.zeros((len(y), 22)); log[:, 2] = y
    o = SimpleNamespace(log=log)
    cases = []
    for tau, mag in ((0., 1.), (.002, 1.13), (.015, 2.7)):
        own = SimpleNamespace(tau_h=tau, mag_h=mag)
        old = C.Preprocess(); old._setup('H', own)
        new = C.Preprocess(); new.h_increment = True; new._setup('H', own)
        for n in range(1, len(y)+1):
            old._log(o, n); new._log(o, n)
        a, b = old.private[:len(y), 2], new.private[:len(y), 2]
        err = float(np.max(abs(np.diff(a)-np.diff(b))))
        assert err < 2e-17
        cases.append(dict(tau=tau, mag=mag, max_increment_difference_m=err,
                          max_level_difference_m=float(np.max(abs(a-b)))))
    # A phase step with no new movement must leave the fixed output still.
    log[:, 2] = np.linspace(0., .018, len(y))
    log[-2:, 2] = .018
    own = SimpleNamespace(tau_h=0., mag_h=1.)
    old = C.Preprocess(); old._setup('H', own)
    new = C.Preprocess(); new.h_increment=True; new._setup('H', own)
    old._log(o, len(y)-1); new._log(o, len(y)-1)
    a, b = old.y_out, new.y_out
    own.mag_h = 1.001
    old._log(o, len(y)); new._log(o, len(y))
    assert new.y_out == b
    return dict(constant_coefficients=cases, phase_step_old_m=old.y_out-a,
                phase_step_fixed_m=new.y_out-b,
                derivation='For constant tau,mag: Delta(y0+(y-y0+tau*dy/TC)/mag)=(dy+tau*(dy-dy_prev)/TC)/mag; dy_0=0. No correction needed. Variable coefficients intentionally act on current increments only.')


def instrument(pol, run):
    windows, commands, steps, contexts = [], [], [], {}
    LF.attach_instrumentation(pol, 'own', windows, commands, steps, contexts)
    if run.startswith('off'):
        SE.gates_off(pol)
        # Preserve the EXACT gates-off execution, adding return capture only.
        source = textwrap.dedent(inspect.getsource(SE.ES.SwapPolicy._stiffness_est))
        start = source.index("    if p['kgood_rule'] == 'filter':")
        end = source.index("    if not p['conservative']:", start)
        source = source[:start]+"    self.k_good = k_hat\n    self.f_set = max(float(o.log[n-1, 5]), .3)\n    self.kg_rel = sig/k_hat\n"+source[end:]
        start = source.index('    if precise and consistent:\n        return max(cand,')
        source = source[:start]+"    return max(cand, (1-p['rho_d'])*self.kbar), True\n"
        tree = ast.parse(source)
        LF.CaptureReturns().visit(tree); ast.fix_missing_locations(tree)
        ns = {}
        exec(compile(tree, '<captured gates-off instance>', 'exec'), SE.ES.__dict__, ns)
        pol._segment_swap = MethodType(ns['_stiffness_est'], pol)
        original_capture = pol._gate_capture

        def capture(v, result):
            out = original_capture(v, result)
            if 'k_hat' in v:
                windows[-1]['accepted'] = bool(np.isfinite(v['k_hat']) and
                    pol.estimator.fit.get('info', 0.) > 0 and
                    np.isfinite(pol.k_good))
                windows[-1]['path'] = 'gates_off'
            return out

        pol._gate_capture = capture
        stiffness = pol._stiffness_est

        def estimate(o, n):
            contexts[n] = dict(kgood=float(pol.k_good), kbar=float(pol.kbar),
                tau_h=float(pol.tau_h), mag_h=float(pol.mag_h), ndel=int(pol.ndel),
                eb=float(pol.eb), ns=int(pol.ns), afc=float(pol.afc))
            return stiffness(o, n)

        pol._stiffness_est = estimate
    return windows, commands, steps, contexts


def forensic_episode(case, run):
    LF.O.setup()
    sc = LF.MC.scene_for(tuple(case['cell']), case['base'], case['variant'])
    pol = LF.B.build('Huber-D4', LF.O.chosen()['Huber-D4'], sc)
    pol.estimator.h_increment = run != 'own'
    windows, commands, steps, contexts = instrument(pol, run)
    log = LF.P3.episode(sc, pol, case['seed'], diagnostic=True)
    assert log is not None
    metric = LF.M.metric(log, pol, sc)
    LF.finish_windows(windows, log, steps, pol.estimator.nw)
    peak_j = int(np.argmax(np.where(log[:, 0] >= 1., log[:, 18], -np.inf)))
    ev = LF.event_summary(windows, steps, log, float(log[-1, 0]))
    row = dict(case=case['label'], run=run, cell=case['cell'], seed=case['seed'],
               peak=float(metric['peak']), L=float(np.log1p(metric['J'])),
               peak_t=float(log[peak_j, 0]), first_low_event=ev['first_upper_exit'],
               first_either_exit=ev['first_either_exit'])
    return row, sc, log, windows, contexts


def replay_fixed(sc, log, contexts, targets):
    source = textwrap.dedent(inspect.getsource(LF.recorded_estimates))
    old = "pol=B.build('Huber-D4',O.chosen()['Huber-D4'],sc)"
    assert source.count(old) == 1
    source = source.replace(old, old+'\n    pol.estimator.h_increment=True')
    ns = {}
    exec(compile(source, '<fixed increment replay>', 'exec'), LF.__dict__, ns)
    return ns['recorded_estimates'](sc, log, contexts, targets)[0]


def case_job(task):
    case, run = task
    row, sc, log, windows, contexts = forensic_episode(case, run)
    stored = json.loads((HERE/'results/low_estimate_forensics.json').read_text(encoding='utf-8'))
    ref = next(r for r in stored['own'] if r['case'] == case['label'])
    wm = {w['n']: w for w in windows}
    row['analysed_windows'] = [dict(n=w['n'], t=w['window']['t'],
        original_k_hat=w['window']['k_hat'], original_k=w['window']['k'],
        closed_loop=wm.get(w['n'])) for w in ref['forensic_windows']]
    if run == 'own':
        row['baseline_peak_delta'] = row['peak']-ref['peak']
        row['baseline_L_delta'] = row['L']-ref['L']
        assert row['baseline_peak_delta'] == 0. and row['baseline_L_delta'] == 0.
        targets = [w['n'] for w in ref['forensic_windows']]
        replay = replay_fixed(sc, log, contexts, targets)
        for w in row['analysed_windows']:
            w['fixed_replay_original_trajectory'] = replay[w['n']]
        original_replay = LF.recorded_estimates(sc, log, contexts, targets)[0]
        row['original_replay_max_khat_delta'] = max(abs(original_replay[n]['k_hat']-wm[n]['k_hat']) for n in targets)
        assert row['original_replay_max_khat_delta'] == 0.
    return row


def stats_job(task):
    case, run = task
    row, sc, log, windows, contexts = forensic_episode(case, run)
    samples = [[r.get(k, np.nan) for k in LF.STAT_FIELDS] for r in windows
               if r['k'] > 0 and np.isfinite(r['k_hat']) and np.isfinite(r['sigma'])]
    return dict(case=case, run=run, peak=row['peak'], L=row['L'], rows=samples)


def phase2_job(task):
    config, cell, variant, q = task
    LF.O.setup()
    sc = LF.MC.scene_for(cell, LF.MC.DEV_BASE, variant)
    pol = LF.B.build('Huber-D4', LF.O.chosen()['Huber-D4'], sc)
    pol.estimator.h_increment = 'hfix' in config
    if config.startswith('off'):
        SE.gates_off(pol)
    kind = 'forget_time' if config.endswith('forget') else 'ewa' if config.endswith('ewa') else None
    sampling = R.make_filter(pol, kind) if kind else None
    records = []
    step = pol.step

    def capture(o):
        command = step(o)
        records.append((float(o.log[o.n-1, 0]), pol.k_good, pol.kbar))
        return command

    pol.step = capture
    log = LF.P3.episode(sc, pol, LF.MC.SEED+q)
    assert log is not None
    rec = np.asarray(records)
    assert len(log) == len(rec) and np.allclose(log[:, 0], rec[:, 0], atol=1e-9, rtol=0)
    metric = LF.M.metric(log, pol, sc)
    row = dict(config=config, cell=list(cell), scene=cell[0], variant=variant,
               dynamics=cell[2], phase=cell[3], i=cell[4], seed=LF.MC.SEED+q,
               L=float(np.log1p(metric['J'])), peak=float(metric['peak']),
               catastrophic=bool(metric['peak'] > 36.5))
    row.update(R.score(log, rec, sc))
    # R.score already retains uncorrected negative extremes; also retain the
    # matching high extreme (needed for both sides with/without local offset).
    with np.errstate(divide='ignore', invalid='ignore'):
        raw = np.log(rec[:, 1:3]/log[:, 21, None])
    for e in row['transitions']:
        a = e['start']; b = a+e['interval_ticks']
        for c, name in enumerate(('k_good', 'k_bar')):
            e['channels'][name]['raw_positive'] = float(np.nanmax(raw[a:b, c]))
    if sampling:
        gaps = np.diff(sampling['updates'])
        assert not len(gaps) or gaps.min() >= R.HOP
        row['sampling'] = dict(updates=len(sampling['updates']), invalid=sampling['invalid'],
            min_gap_ticks=int(gaps.min()) if len(gaps) else None,
            max_gap_ticks=int(gaps.max()) if len(gaps) else None,
            gaps_over_six=int((gaps > R.HOP).sum()))
    return row


def aggregate(rows):
    result = []
    for scene, variant in (('breath_rib', 'nominal'), ('breath_rib', 'slow_thin'), ('breath_fast', 'nominal'), ('all', 'all')):
        for config in CONFIGS:
            rr = [r for r in rows if r['config'] == config and
                  (scene == 'all' or (r['scene'], r['variant']) == (scene, variant))]
            if not rr:
                continue
            g = dict(scene=scene, variant=variant, config=config, n=len(rr),
                L_mean=float(np.mean([r['L'] for r in rr])), L_median=float(np.median([r['L'] for r in rr])),
                max_peak=max(r['peak'] for r in rr),
                above15=[dict(cell=r['cell'], seed=r['seed'], peak=r['peak']) for r in rr if r['peak'] > 15.],
                above36_5=[dict(cell=r['cell'], seed=r['seed'], peak=r['peak']) for r in rr if r['peak'] > 36.5],
                plateau={ch: {key: R.dist([r['plateau'][ch][key] for r in rr]) for key in ('median', 'robust_sd')}
                         for ch in ('k_good', 'k_bar')}, transitions={})
            g['plateau_noncatastrophic'] = {ch: {key: R.dist([r['plateau'][ch][key] for r in rr if not r['catastrophic']])
                for key in ('median', 'robust_sd')} for ch in ('k_good', 'k_bar')}
            for direction in ('rising', 'falling'):
                z = R.summarize_events(rr, 'large', direction)
                events = [e for r in rr if not r['catastrophic'] for e in r['transitions']
                          if e['size_class'] == 'large' and e['direction'] == direction]
                for ch in ('k_good', 'k_bar'):
                    z['channels'][ch]['raw_positive'] = R.dist([e['channels'][ch]['raw_positive'] for e in events])
                z['excluded_catastrophic'] = sum(e['size_class'] == 'large' and e['direction'] == direction
                    for r in rr if r['catastrophic'] for e in r['transitions'])
                z['excluded_prior_contact_loss'] = sum(e['size_class'] == 'large' and e['direction'] == direction
                    for r in rr for e in r['excluded_prior_contact_loss'])
                g['transitions'][direction] = z
            result.append(g)
    return result


def run(ex, fn, tasks):
    out = []
    for future in as_completed([ex.submit(fn, task) for task in tasks]):
        out.append(future.result())
        if len(out) % 12 == 0 or len(out) == len(tasks):
            print(fn.__name__, len(out), '/', len(tasks), flush=True)
    return out


def finalize(data):
    """Cross-check saved primary results only; no additional simulations."""
    rows = data['phase2']['episodes']
    assert len(rows) == 504
    assert all(sum(r['config'] == c for r in rows) == 84 for c in CONFIGS)
    by_key = {(r['config'], r['variant'], r['seed']): r for r in rows}
    assert len(by_key) == 504
    repeated = []
    for run, config in (('own', 'on_none'), ('hfix', 'on_hfix')):
        for s in data['phase1']['statistics'][run]['scene_samples']:
            r = by_key[config, s['case']['variant'], s['case']['seed']]
            assert s['peak'] == r['peak'] and s['L'] == r['L']
            repeated.append((config, s['case']['variant'], s['case']['seed']))
    for r in data['phase1']['cases']:
        if r['run'] == 'off_hfix' and r['case'] in ('dev1', 'dev4'):
            s = by_key['off_hfix', 'slow_thin', r['seed']]
            assert s['peak'] == r['peak'] and s['L'] == r['L']
    for config in CONFIGS:
        rr = [r for r in rows if r['config'] == config]
        assert sum(r['scene'] == 'breath_rib' and r['variant'] == 'nominal' for r in rr) == 36
        assert sum(r['scene'] == 'breath_rib' and r['variant'] == 'slow_thin' for r in rr) == 36
        assert sum(r['scene'] == 'breath_fast' for r in rr) == 12
    baseline = data['phase2']['baseline_check']
    assert baseline['max_abs_L_delta'] == 0. and baseline['max_abs_peak_delta'] == 0.
    filtered = [r for r in rows if 'sampling' in r]
    assert len(filtered) == 168
    assert all(r['sampling']['min_gap_ticks'] >= 6 for r in filtered)
    data['phase2']['summary'] = aggregate(rows)
    rows.sort(key=lambda r: (CONFIGS.index(r['config']), r['scene'], r['variant'], r['seed']))
    data['validation'] = dict(passed=True, own_cases_exact=7, original_target_windows_exact=23,
        population_baseline_episodes_exact=72, population_baseline_windows_exact=88582,
        phase2_baseline_episodes_exact=84, phase2_episodes=504,
        repeated_instrumented_vs_standard_episodes_exact=len(repeated)+2,
        memory_episodes_spacing_checked=len(filtered), min_memory_gap_ticks=6)
    data['code_path_audit'] = dict(
        memory='W/S/WF/WA/WT: y0+((y-y0)+sum(st))/mag_y uses online mag_y and cached rows; same normalization defect. Not changed or simulated here.',
        hc='H/HA/HT share this branch. Only Huber-D4 instances enable h_increment in these experiments.',
        hosts='id_baselines.build: Huber-D4/ContactRobust, D4/ContactEstimator, P5/ContactProjection; Rossi-ID/Rossi-R2/Kalman/MHE/Xue-const/Xue-chi through ContactAdapter._setup(H). Other hosts not rerun or modified.',
        other='CandidatePolicy._next_z HC inversion and MEMORY inversion update dc/Zd by increments, not episode-origin renormalization; sliding corrections likewise cumulative. No other mag_h/mag_y episode-origin transform found in this experiment directory.')
    for name in data['metadata']['source_sha256']:
        data['metadata']['source_sha256'][name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--phase1-only', action='store_true')
    args = ap.parse_args()
    start = time.perf_counter()
    paths = [Path(__file__), HERE/'contact_candidates.py', HERE/'low_estimate_forensics.py',
             HERE/'rib_transition_eval.py', HERE/'memory_segment_eval.py', HERE/'memory_filters2.py']
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if args.resume:
        data = json.loads(OUT.read_text(encoding='utf-8'))
        assert data['metadata']['source_sha256'] == hashes
    else:
        data = dict(metadata=dict(source_sha256=hashes, workers=args.workers, simulation_only=True,
            base=520000, seed_base=LF.MC.SEED, fresh_seeds=False, configs=CONFIGS,
            flag='estimator.h_increment (default False); only HC preprocessing changed',
            low_event='LF first_upper_exit: k/k_good or k/k_bar >1.3 after that channel first entered [0.75,1.3]; searched whole episode',
            replay='Fixed original measurements, calibrated sensing and phase schedule; chronological state, no feedback rerun. Closed-loop same-clock windows also retained.',
            low_tail='LF.population_statistics unchanged: warmup and overlapping windows included, bin-specific median c and existing sigma surrogate',
            phase2_metrics='rib_transition_eval.score/summarize_events unchanged, plus raw_positive; catastrophic episodes >36.5N excluded from transition aggregates; plateau both all/noncatastrophic',
            filters='R.make_filter: forgetting with floor / EWA, HOP=6, R=4*(sigma/k)^2, existing factory constants and likelihood correction unchanged',
            phase2_groups='A nominal 36, A slow_thin 36, B nominal 12; original O.cells indices and dev3 seeds'),
            numerical_checks=numerical_checks(), phase1={}, phase2=dict(episodes=[], completed=[]))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        if 'cases' not in data['phase1']:
            data['phase1']['cases'] = run(ex, case_job, [(c, r) for r in ('own', 'hfix', 'off_hfix') for c in LF.selected_cases()])
            save(data)
        cases = data['phase1']['cases']
        own = {r['case']: r for r in cases if r['run'] == 'own'}
        fixed = {r['case']: r for r in cases if r['run'] == 'hfix'}
        # Predeclared engineering checkpoint, not a winner selection: >=25%
        # peak reduction in >=2 forensic cases, and >=2x correction of at least
        # one originally sub-half-truth window on the unchanged trajectory.
        improved = [k for k in own if fixed[k]['peak'] <= .75*own[k]['peak']]
        corrected = [dict(case=r['case'], n=w['n']) for r in own.values() for w in r['analysed_windows']
                     if w['original_k_hat'] < .5*w['original_k'] and
                     w['fixed_replay_original_trajectory']['k_hat'] >= 2*w['original_k_hat']]
        effect = dict(clear=len(improved) >= 2 and bool(corrected), reduced_peak_cases=improved,
                      corrected_low_windows=corrected,
                      rule='>=25% peak reduction in >=2 cases plus >=2x correction of an originally <0.5*truth window on original trajectory')
        data['phase1']['effect_check'] = effect
        print('phase1_effect', effect, flush=True)
        save(data)
        if 'statistics' not in data['phase1']:
            grid = [c for c in LF.O.cells() if c[0] == 'breath_rib']
            tasks = [(dict(label=f'{v}_{q}', cell=c, base=520000, variant=v, seed=LF.MC.SEED+q), r)
                     for r in ('own', 'hfix') for v in ('nominal', 'slow_thin') for q, c in enumerate(grid)]
            samples = run(ex, stats_job, tasks)
            stats = {}
            for r in ('own', 'hfix'):
                ss = [s for s in samples if s['run'] == r]
                stats[r] = LF.population_statistics(ss)
                stats[r]['by_variant'] = {}
                for v in ('nominal', 'slow_thin'):
                    vs = LF.population_statistics([s for s in ss if s['case']['variant'] == v])
                    vs.pop('scene_samples')
                    stats[r]['by_variant'][v] = vs
            data['phase1']['statistics'] = stats
            ref = json.loads((HERE/'results/low_estimate_forensics.json').read_text(encoding='utf-8'))['statistics']
            refs = {(s['case']['variant'], s['case']['seed']):s for s in ref['scene_samples']}
            checks = [dict(variant=s['case']['variant'], seed=s['case']['seed'],
                peak_delta=s['peak']-refs[(s['case']['variant'],s['case']['seed'])]['peak'],
                L_delta=s['L']-refs[(s['case']['variant'],s['case']['seed'])]['L'],
                samples_exact=LF.clean(s['rows'])==refs[(s['case']['variant'],s['case']['seed'])]['rows'])
                for s in stats['own']['scene_samples']]
            assert all(c['peak_delta']==0. and c['L_delta']==0. and c['samples_exact'] for c in checks)
            data['phase1']['population_baseline_checks'] = checks
            save(data)
        if args.phase1_only or not effect['clear']:
            return
        episodes = [(c, 'nominal', q) for q, c in enumerate(LF.O.cells()) if c[0] in ('breath_rib', 'breath_fast')]
        episodes += [(c, 'slow_thin', q) for q, c in enumerate(LF.O.cells()) if c[0] == 'breath_rib']
        assert len(episodes) == 84
        for config in CONFIGS:
            if config in data['phase2']['completed']:
                continue
            rr = run(ex, phase2_job, [(config, c, v, q) for c, v, q in episodes])
            data['phase2']['episodes'] += rr
            if config == 'on_none':
                check = SE.check_baseline(rr)
                data['phase2']['baseline_check'] = check
                assert check['passed']
            data['phase2']['completed'].append(config)
            save(data)
        finalize(data)
    data['seconds_this_invocation'] = time.perf_counter()-start
    save(data)
    print('saved', OUT, flush=True)


if __name__ == '__main__':
    main()
