# Completed task evidence

This is the focused completion of the October 21 benchmark: isolate stiffness,
isolate surface motion, and retain the original combined scan. The Chinese
[experiment section](../../../research/paper_experiments_final_zh.md) presents
the setup, results, interpretation and original-paper mapping.

## Experimental contract

- Protocol and source hashes were saved before simulation. Seed 7200 supplies
  40 paired scenes per task; each method receives identical scene/noise inputs.
- Thirteen method variants, including two Xue prior variants, are evaluated.
  Learned VIC and C2 each have three frozen training seeds: 17 model instances,
  2,040 primary episodes. Neural models are averaged within each scene before
  computing scene statistics. Training seeds do not multiply independent n.
- T1 removes breathing from the original T1 distribution. T2 uses 20 soft and
  20 stiff homogeneous scenes, balanced by initial push/retract, with a 5 mm
  displacement and return. T3 retains the original main distribution.
- T1 uses frozen T1 settings; T2/T3 use frozen T3 settings. No new task-specific
  tuning or training was performed. Thus T1/T2 are isolation/transfer tests,
  rather than rankings of separately retuned algorithms. T3 is the new primary
  in-distribution comparison. Original baseline corrections remain unchanged.
- Recovery is evaluated after both T2 ramps. Band: 4 +/- 0.5 N, 25 consecutive
  8 ms samples; horizon: 2 s. Unrecovered events are right-censored and included
  as 2 s in the restricted mean, alongside recovery fractions.
- Four C1 and three C2 frozen runtime interventions provide 1,040 additional
  episodes. Input masking in the frozen C2 network is an inference intervention,
  not a replacement for retraining a missing-feature architecture.
- Final-C1 local OCP snapshots use three starts and 25/50/100 GN iterations.
  Numerical cost stability is checked; this is not global optimality, a C2
  certificate, or a reference policy with identical action classes for all methods.
- Online timing is serial and warmed: observation object construction through
  estimation/decision/output projection. Simulated plant integration and physical
  sensor/robot communications are excluded. Fixed/Duan initialization samples
  are not equivalent to periodic decision calls.

## Results and scientific scope

All 3,080 primary/ablation episodes are finite. C1 and C2 have the lowest two
mean force RMSEs on the stiffness and combined tasks. Neither dominates every
task: Haninger's damping-rule adaptation is best on the balanced homogeneous
push/retract cohort. C1 transfers well on the soft subset but poorly on sustained
stiff contact. These are frozen-configuration results.

C1's reference/feedforward action is important: suppressing it raises combined
task RMSE from 0.1471 to 0.4888 N. Its advantages cannot all be attributed to
impedance-gain updates. Normal velocity RMS is reported descriptively; lower
velocity alone is not a tracking-quality claim. Lateral motion is independent of
the normal controller, so lateral tracking improvement cannot be inferred.

All measured online updates are below 8 ms. Optimization is not demonstrated to
be the main bottleneck in this carrier. Final-C1 local-reference gaps are usually
small but have large outliers, so uniform near-optimality is unsupported. These
results do not establish implementation of complete B0 or robot certification.

## Reproduction

From the repository root, use Python with NumPy, Numba, pandas and matplotlib.
All commands use `python -X utf8` on Windows. Existing selected-controller
sources/configurations/weights must retain their recorded hashes.

```powershell
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py protocol
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py check
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py run
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py ablate
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py time
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py reference
python -X utf8 experiments/cobot_scan_benchmark/complete_task_evidence.py summarize
python -X utf8 experiments/cobot_scan_benchmark/write_completed_report.py
```

Figure rendering uses `render_completed_tasks.py --qa-tools <figure-QA scripts>`.
The figure shows the first soft/stiff push scenes (indices 0/2) and model seed 0;
all methods/scenes remain in the result tables. It uses no smoothing or axis
clipping. Alignment, exported PDF font sizes/collisions, and the final image are
checked separately. The figure-QA scripts are supplied by the nature-figure skill;
the simulation and numerical tables do not require that skill.

## Retained evidence

- `protocol.json`, `checks.json`: frozen contract and 51 exact full-log comparisons
  against the original episode implementation, plus recovery/censoring checks.
- `scenes_*.npz`, `per_scene.csv`, `summary.csv`, `paired.csv`: scene parameters,
  primary outcomes, scene statistics, paired bootstrap differences. The 66
  candidate/baseline/task comparisons use Bonferroni-adjusted percentile bootstrap
  intervals (100,000 draws), conditional on the fixed trained models.
- `ablations.csv`, `ablation_summary.csv`: every frozen-intervention outcome.
- `timing.csv`, `timing_raw.npz`: per-model summaries and all raw update durations.
- `reference.csv`: every sampled local OCP cost and iteration-stability comparison.
- `traces/`, `figures/`: first-scene observations/applied gains and selected raw
  illustration traces. `report_manifest.json` hashes the inputs to the draft.

Older seed-3000/4000 results, seed-5000/6000 pressure tests and frozen B0 sources
are preserved. They are not silently pooled into this fresh comparison.
