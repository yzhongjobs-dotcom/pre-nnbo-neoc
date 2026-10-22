"""Physical-outcome tables after baseline repair, with frozen-candidate pairing."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import report_physical_metrics as R
import refine_baselines as F

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'refined_baselines'
NAMES = {'Xue': 'Xue_A1', 'Xue_nominal': 'Xue_nominal_ys0',
         'Duan': 'Duan2018', 'Song': 'Song2025',
         'Jin_OS': 'Jin_OS', 'Jin_MS': 'Jin_MS', 'Learned_VIC': 'Learned_VIC'}
MAIN = ['Best_fixed', 'EtR-A', 'EtR-A-Haninger', 'EtR-R', 'Duan2018', 'Song2025',
        'Jin_OS', 'Jin_MS', 'Xue_A1', 'Xue_nominal_ys0', 'Learned_VIC',
        'Candidate1_v1', 'Candidate1_v2', 'Candidate2']
DIAGNOSTICS = ['Fixed_low', 'Fixed_high', 'Jin_MS_literal',
               'TruthK_0.0s', 'TruthK_0.1s', 'TruthK_0.2s']
LABELS = dict(R.LABELS)
LABELS.update(Xue_A1='Xue, repaired (initial-surface prior)',
              Xue_nominal_ys0='Xue, repaired (nominal surface)',
              Learned_VIC='Learned VIC, repaired/retrained (3 seeds)')


def load_csv(path):
    rows = list(csv.DictReader(path.open(encoding='utf-8', newline='')))
    for row in rows:
        for k, v in row.items():
            if k in ('task', 'method'):
                continue
            if not v:
                row[k] = None
            elif k in ('cohort', 'scene', 'noise_seed', 'training_seed', 'failed'):
                row[k] = int(v)
            else:
                row[k] = float(v)
    return rows


def validate(rows, old, protocol):
    F.verify_frozen(protocol)
    # The old primary data are reusable only while every recorded source,
    # configuration and weight remains byte-identical.
    old_protocol = json.loads((ROOT / 'physical_results/protocol.json').read_text())
    mismatches = [n for n, h in old_protocol['source_sha256'].items()
                  if hashlib.sha256((ROOT / n).read_bytes()).hexdigest() != h]
    assert not mismatches, mismatches
    expected = {(r['cohort'], r['task'], r['method'], r['training_seed'], r['scene']) for r in old}
    keys = [(r['cohort'], r['task'], r['method'], r['training_seed'], r['scene']) for r in rows]
    assert len(keys) == len(set(keys)) and set(keys) == expected
    for row in rows:
        if row['failed']:
            assert row['rmse_N'] is None
            continue
        assert all(np.isfinite(row[m]) for m in R.METRICS[:5])
        assert row['max_abs_error_N'] + 1e-10 >= row['rmse_N'] >= 0
        assert row['max_abs_error_N'] + 1e-10 >= abs(row['peak_force_N'] - 4)
        assert 0 <= row['longest_no_contact_s'] <= row['no_contact_s'] + 1e-10
        assert row['no_contact_s'] <= row['low_contact_s'] + 1e-10 <= row['sampled_seconds'] + 1e-10
        assert row['noise_seed'] == row['cohort'] + row['scene']
    return dict(passed=True, records=len(rows), unique_records=len(set(keys)),
                frozen_records_reused=sum(r['method'] not in NAMES.values() for r in rows),
                repaired_records=sum(r['method'] in NAMES.values() for r in rows),
                failed_runs=sum(r['failed'] for r in rows),
                original_sources_and_weights_unchanged=True,
                reporting_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def fmt(x, digits=4):
    return 'NA' if x is None else f'{x:.{digits}f}'


def table(lines, summary, cohort, task, methods):
    lookup = {r['method']: r for r in summary if r['cohort'] == cohort and r['task'] == task}
    lines += ['| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for name in methods:
        row = lookup[name]
        cols = [row[m + '_mean'] for m in R.METRICS[:4]] + [row['worst_peak_force_N']]
        lines.append('| ' + LABELS[name] + ' | ' + ' | '.join(fmt(x) for x in cols) + f" | {row['failed_runs']} |")
    lines += ['']


def report(summary, old_summary, pairs, budget, learning):
    lines = ['# Fair physical comparison after baseline repair and additional tuning', '',
             '## Material Passport', '',
             'Experiment: PC125 baseline refinement. Status: completed numerical comparison; '
             'source and parameter hashes in protocol.json. Scope: the existing T1/T3 main tasks, '
             'with frozen proposed controllers. Primary data: combined_per_scene.csv; '
             'selection.json and learning.json retain train/dev decisions. This is a task adaptation '
             'study, not a reproduction of the original papers’ robot experiments.', '',
             'The position-controlled robot regulates 4 N at 125 Hz with a one-cycle transport delay, '
             '19 ms servo lag, noisy 1 kHz force sampling and a 30 Hz filter. T1 scans a stiffness '
             'ladder with modulus ratio 30; T3 scans bone-like ridges with 0.5 s patient push/retraction '
             'ramps. Both use 1 mm respiratory motion. Every controller receives the same scene and '
             'noise realization, with 40 scenes per task in each of the seed-3000 and seed-4000 '
             'cohorts. Parameters were selected only on train-1000 and dev-2000. Both test cohorts '
             'had been inspected historically; this is not an untouched-test claim.', '',
             'Force RMSE, per-episode maximum absolute error, peak force and zero-force duration '
             '(F <= 0, an unloaded-contact proxy) use the true simulated force after the original '
             '1 s warm-up at 8 ms intervals. Columns report means across scenes; the worst peak '
             'additionally reports the largest sampled peak. For neural controllers, three training '
             'seeds are averaged within each scene, keeping independent n = 40. Their worst peaks '
             'cover 120 runs and should not be treated as an equal-run-count extreme-value comparison. '
             'Failed integrations remain NA, with the failure count shown. The 7 N threshold is '
             'the inherited tuning constraint, not a clinical safety limit.', '',
             'Frozen proposed/unchanged-reference rows reuse the prior primary episode data after '
             'verifying all recorded source/configuration/weight hashes. Repaired controllers were '
             'rerun at the exact same seeds with unchanged plant/scenario code. Candidate 1 v2 and '
             'Candidate 2 were not retuned. The recently added Candidate 1 v3 uses the same 1 mm '
             'task configurations as v2 here; its separate shared-configuration variant is not '
             'silently substituted. Candidate 1 v1 remains a historical diagnostic and its runaway '
             'is retained.', '',
             '## What changed and what was verified', '',
             '- Xue: restore positive variable stiffness, initialize the constrained stiffness at '
             'kmax, project the full HC parameter vector, and account for the bounded discrete '
             'stiffness change in the energy ledger. Search measured-state synchronization, '
             'bounded servo lead, causal contact reanchoring, model/measured force linearization, '
             'and the initial-surface and nominal-surface variants separately. The selected '
             'parameters, including whether synchronization is enabled, are explicit below. '
             'The original continuous-time passivity claim is not transferred to this delayed '
             'position interface.',
             '- Jin: the old 8 ms one-step endpoint preceded the new command reaching the servo, '
             'making its force term action-blind. Retain a single held decision and extend its '
             'prediction endpoint by train-selected ticks. A finite-action regression gives old '
             'response 0 and corrected response 6.7916 in the normalized force residual. Multistep '
             'equations are retained. Extend objective weights, fixed gains, rate bounds, update '
             'period, horizon and solver iterations; then refine the selected gains and weights.',
             '- Duan/Song: retain the algebraic impedance/compensator adaptations and extend joint '
             'gain searches plus local refinement. Enforce Duan’s sigma bound at Tc = 8 ms and '
             'the common gain limits; enforce Song’s positive native-parameter realizability and '
             'M ki < D. These are mappings to the available discrete interface, not claims of '
             'exactly identical sampling implementations.',
             '- Learned VIC: restore the reference-position output, retain the B1 spring term, '
             'use 8 ms updates with filtered/rate-bounded residuals, and initialize at the same '
             'task-specific B1 controller. Retrain three independent seeds per task with causal '
             'history only. Zero residual reproduces B1 exactly; the compiled training rollout '
             'and the timed Python interface produce identical logs in the regression check. '
             'No privileged teacher is used by this baseline. Candidate 2 retains its original '
             'privileged-teacher preparation, so offline training information/budgets are not '
             'identical across the methods.', '',
             'This closes the specific inactive-branch, prediction-endpoint, omitted-output and '
             'narrow-search defects. Finite random/coordinate searches and retained learning curves '
             'do not prove global tuning optimality or successful transfer of every original method. '
             'In particular, a still-poor Xue/Duan task result is reported as such.', '']
    # Put the actual answer before the complete audit tables.
    at = lines.index('## What changed and what was verified')
    outcome = ['## Does the proposed controller retain an advantage?', '',
               'The reference below has the lowest observed mean RMSE among the competitive '
               'baselines in MAIN. This descriptive selection does not identify the population-optimal '
               'baseline. RMSE reduction does not establish superiority on peak force or every scene.', '',
               '| Cohort | Task | Lowest-RMSE baseline | Baseline RMSE (N) | C1 v2 RMSE (N) | C1 v2 reduction | C2 RMSE (N) | C2 reduction |',
               '|---|---|---|---:|---:|---:|---:|---:|']
    lookup = {(r['cohort'], r['task'], r['method']): r for r in summary}
    for cohort in (4000, 3000):
        for task in ('T1', 'T3'):
            refs = [lookup[cohort, task, n] for n in MAIN if not n.startswith('Candidate')]
            best = min(refs, key=lambda r: r['rmse_N_mean'])
            b = best['rmse_N_mean']
            c1 = lookup[cohort, task, 'Candidate1_v2']['rmse_N_mean']
            c2 = lookup[cohort, task, 'Candidate2']['rmse_N_mean']
            outcome.append(f"| {cohort} | {task} | {LABELS[best['method']]} | {b:.4f} | "
                           f"{c1:.4f} | {100*(b-c1)/b:.1f}% | {c2:.4f} | {100*(b-c2)/b:.1f}% |")
    outcome += ['', 'Positive percentages mean a lower mean RMSE; negative percentages mean a higher '
                'mean RMSE. The complete physical table and paired uncertainty below determine the '
                'scope of this conclusion.', '']
    jm = lookup[4000, 'T3', 'Jin_MS']
    c1 = lookup[4000, 'T3', 'Candidate1_v2']
    c2 = lookup[4000, 'T3', 'Candidate2']
    outcome += [f"The accuracy advantage is not uniform peak-force dominance. On T3/cohort 4000, "
                f"Jin MS has a worst sampled peak of {jm['worst_peak_force_N']:.4f} N versus "
                f"{c1['worst_peak_force_N']:.4f} N for C1 v2 and {c2['worst_peak_force_N']:.4f} N "
                "for C2 (C2's extreme additionally covers three fitted models). Jin MS and C1 v2 "
                "both have zero sampled zero-force time in this cohort; C2 has a small nonzero "
                "mean. C1 v1's 46.9164 N runaway remains in the full table and is not evidence "
                "of a robust proposed controller.", '']
    lines[at:at] = outcome
    for cohort in (4000, 3000):
        for task in ('T1', 'T3'):
            lines += [f'## {task}: main comparisons, cohort {cohort}', '']
            table(lines, summary, cohort, task, MAIN)
    lines += ['## Before/after repair on the same test scenes', '',
              'These test outcomes were computed after parameter selection. They did not feed the search.', '',
              '| Cohort | Task | Baseline | RMSE before → after (N) | Mean max error before → after (N) | Mean peak before → after (N) | Zero-force time before → after (s) |',
              '|---|---|---|---:|---:|---:|---:|']
    for new in summary:
        if new['method'] not in NAMES.values():
            continue
        old = next(r for r in old_summary if all(r[k] == new[k] for k in ('task', 'cohort', 'method')))
        vals = [fmt(old[m + '_mean']) + ' → ' + fmt(new[m + '_mean']) for m in R.METRICS[:4]]
        lines.append(f"| {new['cohort']} | {new['task']} | {LABELS[new['method']]} | " + ' | '.join(vals) + ' |')
    lines += ['', '## Paired uncertainty for proposed controllers', '',
              'Scene-paired bootstrap, 10,000 resamples, three neural seeds averaged within scene. '
              'Intervals condition on the saved fitted models and are marginal/exploratory, without '
              'multiple-comparison adjustment. Training-seed variation is in training_seed_summary.csv. Negative '
              'differences favor the proposed controller. Full metric/reference pairs are in paired_differences.csv.', '',
              '| Cohort | Task | Proposed | Reference | RMSE difference (N) | Marginal 95% interval (N) |',
              '|---|---|---|---|---:|---:|']
    refs = set(MAIN) - {'Candidate1_v1', 'Candidate1_v2', 'Candidate2'}
    for r in pairs:
        if r['candidate'] == 'Candidate1_v1' or r['reference'] not in refs or r['metric'] != 'rmse_N':
            continue
        lines.append(f"| {r['cohort']} | {r['task']} | {LABELS[r['candidate']]} | {LABELS[r['reference']]} | "
                     f"{r['difference_mean']:.4f} | [{r['ci95_low']:.4f}, {r['ci95_high']:.4f}] |")
    lines += ['', '## Search and training budget', '',
              'Screen on 8 training scenes, reevaluate shortlisted configurations on all 20 training '
              'scenes, then select using 20 development scenes. Two subsequent coordinate-refinement '
              'rounds include the incumbent and reevaluate all configurations on 20 training scenes; '
              'development selects among the top three and the incumbent. Tests are not used. '
              'The inherited scalar J selects parameters; the tables above expose its physical '
              'components rather than presenting J as a physical measurement.', '',
              '| Task | Family | Coarse configurations | Local configurations evaluated | Train J | Dev J | Train peak criterion | Dev peak criterion |',
              '|---|---|---:|---:|---:|---:|---|---|']
    for r in budget:
        lines.append(f"| {r['task']} | {r['family']} | {r['coarse_configs']} | {r['local_configs']} | "
                     f"{r['train_J']:.5f} | {r['dev_J']:.5f} | {r['train_feasible']} | {r['dev_feasible']} |")
    lines += ['', '| Task | Training seed | Generations | Training/evaluation episodes | Best generation | Train J | Dev J |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for key, r in learning.items():
        task, seed = key.split('_')
        lines.append(f"| {task} | {seed} | {r['gens']} | {r['episodes']} | {r['best_gen']} | "
                     f"{r['train']['J']:.5f} | {r['dev']['J']:.5f} |")
    lines += ['', 'Learning: population 16; six common training scenes per generation; initial ES '
              'sigma 0.025; network 12–8–3, with learned constant gain offsets; 100 generations '
              'per task/seed. Train/dev checkpoint histories and all final parameter vectors are '
              'in learning.json. A best generation of zero means the finite training budget did '
              'not improve on B1 under the declared selection rule. This is not reported as '
              'learning convergence. After training, each selected network receives a 15-setting '
              'train/dev calibration of output filtering (tau = 0.016/0.04/0.08 s) and reference '
              'span (0.25/0.5/1/3/10 mm), including the incumbent. These evaluation episodes are '
              'included above; weights are unchanged by calibration.', '', '## Selected parameter vectors', '']
    selection = F.read_refined('selection.json', {})
    for key, r in selection.items():
        lines += [f"- {key}: `{json.dumps(r['spec'][1], sort_keys=True)}`"]
    lines += ['', '## Fixed extremes, unanchored diagnostic and truth references', '',
              'These rows are unchanged historical controls/reference bounds; they were not '
              'substituted for repaired competitive baselines. Truth references use privileged '
              'online stiffness and are not deployable baselines.', '']
    for cohort in (4000, 3000):
        for task in ('T1', 'T3'):
            lines += [f'### {task}, cohort {cohort}', '']
            table(lines, summary, cohort, task, DIAGNOSTICS)
    lines += ['## Reproduction', '', '```powershell',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py checks',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py tune --workers 6',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py polish --workers 6',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py learn --workers 6',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py calibrate --workers 6',
              'python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py test --workers 6',
              'python -X utf8 experiments/cobot_scan_benchmark/report_refined_baselines.py',
              '```', '',
              'Completed search/training entries are reused. protocol.json freezes selection, '
              'weights, source files and cached reference data before testing; changed hashes '
              'stop the test. Runtime samples are diagnostic on this workstation and are not '
              'real-time certification; old and new concurrent-process timing conditions differ. '
              'No task/surface/controller values in the frozen B0 snapshot were edited.', '']
    (OUT / 'results.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    protocol = json.loads((OUT / 'protocol.json').read_text())
    old = load_csv(ROOT / 'physical_results/per_scene.csv')
    new = json.loads((OUT / 'test_records.json').read_text())
    for r in new:
        r['method'] = NAMES[r['method']]
    rows = [r for r in old if r['method'] not in NAMES.values()] + new
    fields = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k) for k in fields} for r in rows]
    validation = validate(rows, old, protocol)
    R.write_csv(OUT / 'combined_per_scene.csv', rows)
    summary = R.aggregate(rows)
    old_summary = R.aggregate(old)
    pairs = R.paired(rows)
    R.write_csv(OUT / 'summary.csv', summary)
    R.write_csv(OUT / 'paired_differences.csv', pairs)
    seeds = []
    for name in ('Learned_VIC', 'Candidate2'):
        for seed in range(3):
            seeds += [dict(r, training_seed=seed) for r in R.aggregate(
                [v for v in rows if v['method'] == name and v['training_seed'] == seed])]
    R.write_csv(OUT / 'training_seed_summary.csv', seeds)
    selection = F.read_refined('selection.json', {})
    learning = F.read_refined('learning.json', {})
    budget = []
    for key, r in selection.items():
        task, family = key.split('_', 1)
        coarse = r['coarse']
        budget.append(dict(task=task, family=family, coarse_configs=coarse['n_screen'],
                           local_configs=sum(v['n_full'] for v in r['local_rounds']),
                           train_J=r['train']['J'], dev_J=r['dev']['J'],
                           train_feasible=r['train']['feasible'], dev_feasible=r['dev']['feasible']))
    R.write_csv(OUT / 'search_budget.csv', budget)
    report(summary, old_summary, pairs, budget, learning)
    F.save('validation.json', validation)
    print(json.dumps(validation))


if __name__ == '__main__':
    main()
