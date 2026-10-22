# Physical outcomes of the cobot scanning benchmark

The PC125 position-controlled plant runs at 125 Hz with a 4 N force target, one control-cycle transport delay, a 19 ms servo lag, noisy 1 kHz force sensing and a 30 Hz filter. T1 scans a stiffness ladder (modulus ratio 30); T3 scans bone-like ridges with 0.5 s push/retraction ramps. Both main tasks use 1 mm respiration. Each cohort has 40 paired scenes per task. All controller parameters and neural weights are frozen from the existing training records; no retuning is performed.

Seed 3000 is the original, already inspected test set. Seed 4000 is the previously used v2 evaluation set, extended here to all controllers. Neither is described as a new untouched test set. Three neural training seeds are averaged within each scene; the independent sample size remains 40, not 120. Worst sampled peaks for neural methods cover all 120 scene/seed runs and are descriptive, not an equal-run-count comparison against a single fitted controller.

Force outcomes use the true simulated force after the original 1 s warm-up. RMSE, maximum absolute error and peak force are computed per episode, then averaged over scenes. The worst-peak column additionally gives the largest sampled peak across all scenes/seeds. Zero-force duration counts F <= 0 as an unloaded-contact proxy, without the old 0.5 N operational cutoff. Failed integrations remain in the data and give NA physical means, rather than silently dropping scenes. Samples are 8 ms apart, so all peaks and durations are sampled quantities. No clinical safety threshold is asserted. Units are N and s.

## T3, cohort seed 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | NA | NA | NA | NA | NA | 1 |
| Fixed high damping | 5.8100 | 23.0524 | 27.0524 | 4.1216 | 66.4453 | 0 |
| Best tuned fixed (B1) | 0.3110 | 1.5365 | 5.4618 | 0.0042 | 12.0202 | 0 |
| Estimate + stiffness rule | 0.2710 | 1.6107 | 5.5094 | 0.0020 | 10.7190 | 0 |
| Estimate + Haninger damping rule | 0.2976 | 1.5755 | 5.4808 | 0.0000 | 7.2436 | 0 |
| RLS estimate + rule | 0.5556 | 5.2244 | 9.0959 | 0.0350 | 43.1995 | 0 |
| Duan (RAS 2018), adapted | 1.0546 | 4.5967 | 8.5478 | 0.0998 | 22.5533 | 0 |
| Song (RA-L 2025), adapted | 0.4396 | 1.9211 | 5.8644 | 0.0066 | 11.5774 | 0 |
| Jin (TMECH 2023), one-step | 0.3156 | 1.5985 | 5.4971 | 0.0032 | 11.8837 | 0 |
| Jin (TMECH 2023), multistep | 0.3490 | 1.7046 | 5.5597 | 0.0030 | 7.7716 | 0 |
| Jin multistep, unanchored diagnostic | 1.6606 | 12.5960 | 16.5430 | 0.3668 | 58.0744 | 0 |
| Xue MPVIC, initial-surface prior | 1.2762 | 4.3995 | 8.2126 | 0.2082 | 13.3817 | 0 |
| Xue MPVIC, nominal initial surface | 1.2593 | 4.6191 | 8.4276 | 0.1350 | 19.9934 | 0 |
| E27-inspired learned VIC (3 seeds) | 0.3179 | 1.4615 | 5.4004 | 0.0011 | 13.9957 | 0 |
| Candidate 1 v1 | 0.1633 | 1.1804 | 5.0847 | 0.0000 | 7.1215 | 0 |
| Candidate 1 v2 | 0.1627 | 1.1865 | 5.0912 | 0.0000 | 7.1222 | 0 |
| Candidate 2 (3 seeds) | 0.2260 | 1.3582 | 5.2411 | 0.0000 | 7.9067 | 0 |
| Truth-stiffness reference, 0 s delay | 0.1452 | 0.8919 | 4.7769 | 0.0000 | 6.5007 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1538 | 0.8966 | 4.7809 | 0.0000 | 6.4518 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.3041 | 1.4630 | 5.3045 | 0.0004 | 7.8143 | 0 |

## T1, cohort seed 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 15.4310 | 152.2074 | 156.2074 | 5.3912 | 261.1714 | 0 |
| Fixed high damping | 3.1345 | 11.0831 | 15.0244 | 1.1708 | 30.1981 | 0 |
| Best tuned fixed (B1) | 0.1554 | 0.7883 | 4.7809 | 0.0034 | 9.3978 | 0 |
| Estimate + stiffness rule | 0.1249 | 0.6839 | 4.6411 | 0.0000 | 6.2447 | 0 |
| Estimate + Haninger damping rule | 0.1591 | 1.1194 | 5.0796 | 0.0018 | 10.6807 | 0 |
| RLS estimate + rule | 0.2542 | 1.4526 | 5.4009 | 0.0000 | 7.5865 | 0 |
| Duan (RAS 2018), adapted | 0.6236 | 2.8784 | 6.8702 | 0.0344 | 14.6680 | 0 |
| Song (RA-L 2025), adapted | 0.3219 | 1.2124 | 5.1954 | 0.0024 | 10.0896 | 0 |
| Jin (TMECH 2023), one-step | 0.1489 | 0.8888 | 4.8820 | 0.0024 | 13.2501 | 0 |
| Jin (TMECH 2023), multistep | 0.1396 | 0.7251 | 4.7154 | 0.0000 | 7.8094 | 0 |
| Jin multistep, unanchored diagnostic | 0.7074 | 5.3513 | 9.2668 | 0.0570 | 22.8351 | 0 |
| Xue MPVIC, initial-surface prior | 0.5472 | 2.1467 | 5.9855 | 0.0000 | 8.8504 | 0 |
| Xue MPVIC, nominal initial surface | 0.5561 | 2.2232 | 6.0358 | 0.0000 | 9.1001 | 0 |
| E27-inspired learned VIC (3 seeds) | 0.1837 | 1.0675 | 5.0532 | 0.0086 | 32.7635 | 0 |
| Candidate 1 v1 | 0.0672 | 0.3462 | 4.3270 | 0.0000 | 5.0026 | 0 |
| Candidate 1 v2 | 0.0674 | 0.3454 | 4.3249 | 0.0000 | 4.9782 | 0 |
| Candidate 2 (3 seeds) | 0.0934 | 0.4759 | 4.4463 | 0.0000 | 5.4009 | 0 |
| Truth-stiffness reference, 0 s delay | 0.0674 | 0.3434 | 4.3349 | 0.0000 | 5.1022 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1065 | 0.5170 | 4.4413 | 0.0000 | 5.8104 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.1273 | 0.7465 | 4.7042 | 0.0000 | 7.6770 | 0 |

## T3, cohort seed 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 18.0291 | 149.2457 | 153.2457 | 8.2284 | 306.6901 | 0 |
| Fixed high damping | 6.2277 | 25.0353 | 29.0353 | 4.4138 | 52.3818 | 0 |
| Best tuned fixed (B1) | 0.3236 | 1.6481 | 5.5519 | 0.0062 | 9.8619 | 0 |
| Estimate + stiffness rule | 0.2726 | 1.6205 | 5.4244 | 0.0032 | 10.3161 | 0 |
| Estimate + Haninger damping rule | 0.2949 | 1.5992 | 5.4440 | 0.0000 | 6.9814 | 0 |
| RLS estimate + rule | 0.5352 | 4.6341 | 8.5816 | 0.0256 | 15.1715 | 0 |
| Duan (RAS 2018), adapted | 1.1371 | 4.8019 | 8.7045 | 0.1726 | 19.0302 | 0 |
| Song (RA-L 2025), adapted | 0.4374 | 2.0736 | 5.9655 | 0.0052 | 9.7068 | 0 |
| Jin (TMECH 2023), one-step | 0.3297 | 1.8181 | 5.6944 | 0.0038 | 9.7644 | 0 |
| Jin (TMECH 2023), multistep | 0.3523 | 1.7903 | 5.6591 | 0.0004 | 9.5561 | 0 |
| Jin multistep, unanchored diagnostic | 1.8048 | 14.3591 | 18.2691 | 0.4532 | 55.3832 | 0 |
| Xue MPVIC, initial-surface prior | 1.2699 | 4.5187 | 8.3413 | 0.2142 | 12.5215 | 0 |
| Xue MPVIC, nominal initial surface | 1.2359 | 4.4236 | 8.2952 | 0.1842 | 12.2522 | 0 |
| E27-inspired learned VIC (3 seeds) | 0.3185 | 1.5964 | 5.4604 | 0.0005 | 10.8890 | 0 |
| Candidate 1 v1 | 0.2903 | 2.3068 | 6.1326 | 0.0550 | 46.9164 | 0 |
| Candidate 1 v2 | 0.1726 | 1.3044 | 5.1007 | 0.0000 | 7.1640 | 0 |
| Candidate 2 (3 seeds) | 0.2308 | 1.4369 | 5.3028 | 0.0001 | 7.9063 | 0 |
| Truth-stiffness reference, 0 s delay | 0.1524 | 0.9506 | 4.8061 | 0.0000 | 6.1974 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1654 | 0.9869 | 4.8685 | 0.0000 | 6.4316 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.3327 | 1.6563 | 5.4878 | 0.0018 | 8.2387 | 0 |

## T1, cohort seed 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 15.3557 | 147.0219 | 151.0219 | 5.4516 | 269.7349 | 0 |
| Fixed high damping | 2.9783 | 9.9165 | 13.8014 | 1.1018 | 27.2236 | 0 |
| Best tuned fixed (B1) | 0.1472 | 0.7565 | 4.7375 | 0.0026 | 9.3777 | 0 |
| Estimate + stiffness rule | 0.1204 | 0.6855 | 4.6230 | 0.0000 | 5.9249 | 0 |
| Estimate + Haninger damping rule | 0.1402 | 0.9353 | 4.8861 | 0.0000 | 8.2923 | 0 |
| RLS estimate + rule | 0.2447 | 1.4674 | 5.4174 | 0.0000 | 7.6144 | 0 |
| Duan (RAS 2018), adapted | 0.5927 | 2.8776 | 6.8664 | 0.0344 | 13.8899 | 0 |
| Song (RA-L 2025), adapted | 0.2892 | 1.0662 | 5.0397 | 0.0000 | 7.1594 | 0 |
| Jin (TMECH 2023), one-step | 0.1230 | 0.6492 | 4.6346 | 0.0000 | 6.8259 | 0 |
| Jin (TMECH 2023), multistep | 0.1302 | 0.6307 | 4.6132 | 0.0000 | 7.2333 | 0 |
| Jin multistep, unanchored diagnostic | 0.5888 | 4.1478 | 8.1270 | 0.0464 | 16.3158 | 0 |
| Xue MPVIC, initial-surface prior | 0.5147 | 2.0301 | 5.8833 | 0.0000 | 9.9906 | 0 |
| Xue MPVIC, nominal initial surface | 0.5165 | 2.0010 | 5.7827 | 0.0000 | 9.8201 | 0 |
| E27-inspired learned VIC (3 seeds) | 0.1878 | 1.1144 | 5.0963 | 0.0093 | 24.5411 | 0 |
| Candidate 1 v1 | 0.0669 | 0.3547 | 4.3154 | 0.0000 | 5.5974 | 0 |
| Candidate 1 v2 | 0.0673 | 0.3604 | 4.3123 | 0.0000 | 5.6032 | 0 |
| Candidate 2 (3 seeds) | 0.0896 | 0.4739 | 4.4294 | 0.0000 | 5.9867 | 0 |
| Truth-stiffness reference, 0 s delay | 0.0658 | 0.3196 | 4.3052 | 0.0000 | 5.4954 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.0996 | 0.4753 | 4.3873 | 0.0000 | 6.1857 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.1188 | 0.7081 | 4.6483 | 0.0000 | 8.5450 | 0 |

## Policy-update timing

Timing includes each policy's estimator, feature extraction and decision call. It excludes shared plant/servo execution, communication, and JIT compilation. Warm-up uses training scenes. Measurements are host prototype timings during a two-process replay, not robot real-time certification. Each table entry averages per-episode medians or 95th percentiles; it is not the pooled call-distribution percentile. A dash means a fixed parameter vector, not a claim of zero control computation.

| Controller | Decision period (ms) | T1 mean episode-median (ms) | T1 mean episode-p95 (ms) | T3 mean episode-median (ms) | T3 mean episode-p95 (ms) |
|---|---:|---:|---:|---:|---:|
| Fixed low damping | — | — | — | — | — |
| Fixed high damping | — | — | — | — | — |
| Best tuned fixed (B1) | — | — | — | — | — |
| Estimate + stiffness rule | 8 | 0.02482 | 0.02757 | 0.03181 | 0.04077 |
| Estimate + Haninger damping rule | 8 | 0.02478 | 0.02723 | 0.03213 | 0.04097 |
| RLS estimate + rule | 8 | 0.00477 | 0.00532 | 0.00515 | 0.00621 |
| Duan (RAS 2018), adapted | — | — | — | — | — |
| Song (RA-L 2025), adapted | 8 | 0.00133 | 0.00149 | 0.00133 | 0.00152 |
| Jin (TMECH 2023), one-step | 16 | 0.05480 | 0.06063 | 0.05657 | 0.06582 |
| Jin (TMECH 2023), multistep | 16 | 0.24762 | 0.28938 | 0.27602 | 0.32918 |
| Jin multistep, unanchored diagnostic | 16 | 0.24592 | 0.29808 | 0.24160 | 0.28835 |
| Xue MPVIC, initial-surface prior | 8 | 0.10101 | 0.11719 | 0.09619 | 0.11052 |
| Xue MPVIC, nominal initial surface | 8 | 0.09905 | 0.11373 | 0.09497 | 0.10759 |
| E27-inspired learned VIC (3 seeds) | 48 | 0.08150 | 0.09323 | 0.08122 | 0.09469 |
| Candidate 1 v1 | 8 | 0.08742 | 0.09671 | 0.08627 | 0.09934 |
| Candidate 1 v2 | 8 | 0.09014 | 0.10589 | 0.08847 | 0.10103 |
| Candidate 2 (3 seeds) | 8 | 0.01511 | 0.01660 | 0.01481 | 0.01572 |
| Truth-stiffness reference, 0 s delay | 8 | 0.00224 | 0.00256 | 0.00223 | 0.00248 |
| Truth-stiffness reference, 0.1 s delay | 8 | 0.00229 | 0.00256 | 0.00217 | 0.00240 |
| Truth-stiffness reference, 0.2 s delay | 8 | 0.00221 | 0.00254 | 0.00210 | 0.00248 |

## Replay agreement

Original files are unchanged. Physical outcomes above are newly replayed values. The comparison tolerance is abs(error) <= 1e-8 + 1e-5 abs(archived); timing is not compared across machines. See replay_checks.csv for every check. 19 of 232 archived-metric checks differ beyond this strict tolerance, including unstable Jin_MS_literal trajectories and the failed low-damping extreme. This is not a claim of exact cross-environment reproduction. Small differences also occur in selected Jin, Xue and v1 outputs. The cause of the unstable variant's larger discrepancy has not been isolated.

Flagged cohort/task/controller groups: 3000/T1/Jin_MS_literal; 3000/T1/Xue_A1; 3000/T3/Fixed_low; 3000/T3/Jin_MS; 3000/T3/Jin_MS_literal; 3000/T3/Jin_OS; 4000/T3/Candidate1_v1.

## Observations on the supplementary seed-4000 cohort

- T1: best-fixed mean force RMSE 0.1472 N, Candidate1_v2 0.0673 N and Candidate2 0.0896 N. The corresponding mean episode-max absolute errors are 0.7565, 0.3604 and 0.4739 N.
- T3: best-fixed mean force RMSE 0.3236 N, Candidate1_v2 0.1726 N and Candidate2 0.2308 N. The corresponding mean episode-max absolute errors are 1.6481, 1.3044 and 1.4369 N.
- T3 v1's worst sampled force is 46.92 N; v2 reduces it to 7.16 N. The Haninger-rule comparator has a lower worst sampled peak (6.98 N) despite its larger force RMSE. Mean tracking gains do not imply dominance on every outcome.
- The 0.25 s push and 4 mm respiration stress cells are outside these two main-cell replays. Their existing negative findings remain in candidate1_v2_log.md and are not replaced by these positive main-task results.

## Scope and reproducibility

These are physical performance tables, not the old weighted tuning J. The original J is retained only in per-scene source data for replay checks. The exact scenario, configuration, source and weight hashes are in protocol.json. Source data and scene-level SDs are in per_scene.csv and summary.csv. Seed-specific neural summaries are in training_seed_summary.csv. Pairwise scene-bootstrap intervals are in paired_differences.csv, without multiplicity correction; no significance stars or universal superiority claim is attached.

Jin and Xue rows retain the existing parameter searches and adaptation choices. In particular, Xue has kmax=0 and no train-feasible grid setting; those rows are not evidence that a fully tuned Xue family is inferior. Learned_VIC uses the existing E27-inspired ES adaptation, with unchanged weights. TruthK rows use unavailable true stiffness and are information references, not deployable baselines or globally optimal controllers. The per-scene fixed-grid B2 and isolated preview/feedforward oracle probes remain historical gate diagnostics; their missing physical outputs are not fabricated in this report.

The lateral axis is position-servoed and not adapted by these policies. No claim of improved lateral tracking, complete B0 validation or certified near-optimality follows from these normal-force tables.

Reproduce: `python -X utf8 experiments/cobot_scan_benchmark/report_physical_metrics.py`.
