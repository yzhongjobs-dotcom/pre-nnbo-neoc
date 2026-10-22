# Fair physical comparison after baseline repair and additional tuning

## Material Passport

Experiment: PC125 baseline refinement. Status: completed numerical comparison; source and parameter hashes in protocol.json. Scope: the existing T1/T3 main tasks, with frozen proposed controllers. Primary data: combined_per_scene.csv; selection.json and learning.json retain train/dev decisions. This is a task adaptation study, not a reproduction of the original papers’ robot experiments.

The position-controlled robot regulates 4 N at 125 Hz with a one-cycle transport delay, 19 ms servo lag, noisy 1 kHz force sampling and a 30 Hz filter. T1 scans a stiffness ladder with modulus ratio 30; T3 scans bone-like ridges with 0.5 s patient push/retraction ramps. Both use 1 mm respiratory motion. Every controller receives the same scene and noise realization, with 40 scenes per task in each of the seed-3000 and seed-4000 cohorts. Parameters were selected only on train-1000 and dev-2000. Both test cohorts had been inspected historically; this is not an untouched-test claim.

Force RMSE, per-episode maximum absolute error, peak force and zero-force duration (F <= 0, an unloaded-contact proxy) use the true simulated force after the original 1 s warm-up at 8 ms intervals. Columns report means across scenes; the worst peak additionally reports the largest sampled peak. For neural controllers, three training seeds are averaged within each scene, keeping independent n = 40. Their worst peaks cover 120 runs and should not be treated as an equal-run-count extreme-value comparison. Failed integrations remain NA, with the failure count shown. The 7 N threshold is the inherited tuning constraint, not a clinical safety limit.

Frozen proposed/unchanged-reference rows reuse the prior primary episode data after verifying all recorded source/configuration/weight hashes. Repaired controllers were rerun at the exact same seeds with unchanged plant/scenario code. Candidate 1 v2 and Candidate 2 were not retuned. The recently added Candidate 1 v3 uses the same 1 mm task configurations as v2 here; its separate shared-configuration variant is not silently substituted. Candidate 1 v1 remains a historical diagnostic and its runaway is retained.

## Does the proposed controller retain an advantage?

The reference below has the lowest observed mean RMSE among the competitive baselines in MAIN. This descriptive selection does not identify the population-optimal baseline. RMSE reduction does not establish superiority on peak force or every scene.

| Cohort | Task | Lowest-RMSE baseline | Baseline RMSE (N) | C1 v2 RMSE (N) | C1 v2 reduction | C2 RMSE (N) | C2 reduction |
|---|---|---|---:|---:|---:|---:|---:|
| 4000 | T1 | Jin (TMECH 2023), one-step | 0.1161 | 0.0673 | 42.1% | 0.0896 | 22.8% |
| 4000 | T3 | Estimate + stiffness rule | 0.2726 | 0.1726 | 36.7% | 0.2308 | 15.3% |
| 3000 | T1 | Learned VIC, repaired/retrained (3 seeds) | 0.1235 | 0.0674 | 45.4% | 0.0934 | 24.3% |
| 3000 | T3 | Estimate + stiffness rule | 0.2710 | 0.1627 | 40.0% | 0.2260 | 16.6% |

Positive percentages mean a lower mean RMSE; negative percentages mean a higher mean RMSE. The complete physical table and paired uncertainty below determine the scope of this conclusion.

The accuracy advantage is not uniform peak-force dominance. On T3/cohort 4000, Jin MS has a worst sampled peak of 6.4988 N versus 7.1640 N for C1 v2 and 7.9063 N for C2 (C2's extreme additionally covers three fitted models). Jin MS and C1 v2 both have zero sampled zero-force time in this cohort; C2 has a small nonzero mean. C1 v1's 46.9164 N runaway remains in the full table and is not evidence of a robust proposed controller.

## What changed and what was verified

- Xue: restore positive variable stiffness, initialize the constrained stiffness at kmax, project the full HC parameter vector, and account for the bounded discrete stiffness change in the energy ledger. Search measured-state synchronization, bounded servo lead, causal contact reanchoring, model/measured force linearization, and the initial-surface and nominal-surface variants separately. The selected parameters, including whether synchronization is enabled, are explicit below. The original continuous-time passivity claim is not transferred to this delayed position interface.
- Jin: the old 8 ms one-step endpoint preceded the new command reaching the servo, making its force term action-blind. Retain a single held decision and extend its prediction endpoint by train-selected ticks. A finite-action regression gives old response 0 and corrected response 6.7916 in the normalized force residual. Multistep equations are retained. Extend objective weights, fixed gains, rate bounds, update period, horizon and solver iterations; then refine the selected gains and weights.
- Duan/Song: retain the algebraic impedance/compensator adaptations and extend joint gain searches plus local refinement. Enforce Duan’s sigma bound at Tc = 8 ms and the common gain limits; enforce Song’s positive native-parameter realizability and M ki < D. These are mappings to the available discrete interface, not claims of exactly identical sampling implementations.
- Learned VIC: restore the reference-position output, retain the B1 spring term, use 8 ms updates with filtered/rate-bounded residuals, and initialize at the same task-specific B1 controller. Retrain three independent seeds per task with causal history only. Zero residual reproduces B1 exactly; the compiled training rollout and the timed Python interface produce identical logs in the regression check. No privileged teacher is used by this baseline. Candidate 2 retains its original privileged-teacher preparation, so offline training information/budgets are not identical across the methods.

This closes the specific inactive-branch, prediction-endpoint, omitted-output and narrow-search defects. Finite random/coordinate searches and retained learning curves do not prove global tuning optimality or successful transfer of every original method. In particular, a still-poor Xue/Duan task result is reported as such.

## T1: main comparisons, cohort 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Best tuned fixed (B1) | 0.1472 | 0.7565 | 4.7375 | 0.0026 | 9.3777 | 0 |
| Estimate + stiffness rule | 0.1204 | 0.6855 | 4.6230 | 0.0000 | 5.9249 | 0 |
| Estimate + Haninger damping rule | 0.1402 | 0.9353 | 4.8861 | 0.0000 | 8.2923 | 0 |
| RLS estimate + rule | 0.2447 | 1.4674 | 5.4174 | 0.0000 | 7.6144 | 0 |
| Duan (RAS 2018), adapted | 0.5654 | 2.5971 | 6.5930 | 0.0122 | 14.4272 | 0 |
| Song (RA-L 2025), adapted | 0.1963 | 0.7649 | 4.7538 | 0.0000 | 6.5617 | 0 |
| Jin (TMECH 2023), one-step | 0.1161 | 0.5920 | 4.5765 | 0.0000 | 6.8467 | 0 |
| Jin (TMECH 2023), multistep | 0.1247 | 0.6057 | 4.5911 | 0.0000 | 6.9013 | 0 |
| Xue, repaired (initial-surface prior) | 0.3281 | 1.2853 | 5.2338 | 0.0000 | 7.6957 | 0 |
| Xue, repaired (nominal surface) | 0.3314 | 1.2431 | 5.1973 | 0.0000 | 7.3440 | 0 |
| Learned VIC, repaired/retrained (3 seeds) | 0.1179 | 0.5937 | 4.5783 | 0.0007 | 9.3282 | 0 |
| Candidate 1 v1 | 0.0669 | 0.3547 | 4.3154 | 0.0000 | 5.5974 | 0 |
| Candidate 1 v2 | 0.0673 | 0.3604 | 4.3123 | 0.0000 | 5.6032 | 0 |
| Candidate 2 (3 seeds) | 0.0896 | 0.4739 | 4.4294 | 0.0000 | 5.9867 | 0 |

## T3: main comparisons, cohort 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Best tuned fixed (B1) | 0.3236 | 1.6481 | 5.5519 | 0.0062 | 9.8619 | 0 |
| Estimate + stiffness rule | 0.2726 | 1.6205 | 5.4244 | 0.0032 | 10.3161 | 0 |
| Estimate + Haninger damping rule | 0.2949 | 1.5992 | 5.4440 | 0.0000 | 6.9814 | 0 |
| RLS estimate + rule | 0.5352 | 4.6341 | 8.5816 | 0.0256 | 15.1715 | 0 |
| Duan (RAS 2018), adapted | 1.1370 | 4.8073 | 8.7098 | 0.1730 | 19.0161 | 0 |
| Song (RA-L 2025), adapted | 0.4317 | 2.0551 | 5.9507 | 0.0030 | 9.6768 | 0 |
| Jin (TMECH 2023), one-step | 0.2915 | 1.5274 | 5.4442 | 0.0004 | 7.4090 | 0 |
| Jin (TMECH 2023), multistep | 0.2819 | 1.4494 | 5.3714 | 0.0000 | 6.4988 | 0 |
| Xue, repaired (initial-surface prior) | 0.6939 | 2.5474 | 6.3557 | 0.1142 | 13.6064 | 0 |
| Xue, repaired (nominal surface) | 0.5997 | 2.3490 | 6.2597 | 0.0484 | 13.8924 | 0 |
| Learned VIC, repaired/retrained (3 seeds) | 0.3096 | 1.4300 | 5.3698 | 0.0158 | 10.5038 | 0 |
| Candidate 1 v1 | 0.2903 | 2.3068 | 6.1326 | 0.0550 | 46.9164 | 0 |
| Candidate 1 v2 | 0.1726 | 1.3044 | 5.1007 | 0.0000 | 7.1640 | 0 |
| Candidate 2 (3 seeds) | 0.2308 | 1.4369 | 5.3028 | 0.0001 | 7.9063 | 0 |

## T1: main comparisons, cohort 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Best tuned fixed (B1) | 0.1554 | 0.7883 | 4.7809 | 0.0034 | 9.3978 | 0 |
| Estimate + stiffness rule | 0.1249 | 0.6839 | 4.6411 | 0.0000 | 6.2447 | 0 |
| Estimate + Haninger damping rule | 0.1591 | 1.1194 | 5.0796 | 0.0018 | 10.6807 | 0 |
| RLS estimate + rule | 0.2542 | 1.4526 | 5.4009 | 0.0000 | 7.5865 | 0 |
| Duan (RAS 2018), adapted | 0.6337 | 2.9808 | 6.9771 | 0.0144 | 13.9570 | 0 |
| Song (RA-L 2025), adapted | 0.2110 | 0.8574 | 4.8489 | 0.0000 | 6.5456 | 0 |
| Jin (TMECH 2023), one-step | 0.1264 | 0.6668 | 4.6598 | 0.0000 | 6.4357 | 0 |
| Jin (TMECH 2023), multistep | 0.1359 | 0.6864 | 4.6790 | 0.0000 | 6.4277 | 0 |
| Xue, repaired (initial-surface prior) | 0.3514 | 1.3875 | 5.3568 | 0.0000 | 7.7038 | 0 |
| Xue, repaired (nominal surface) | 0.3419 | 1.3561 | 5.3151 | 0.0000 | 7.4130 | 0 |
| Learned VIC, repaired/retrained (3 seeds) | 0.1235 | 0.6176 | 4.6093 | 0.0013 | 9.4256 | 0 |
| Candidate 1 v1 | 0.0672 | 0.3462 | 4.3270 | 0.0000 | 5.0026 | 0 |
| Candidate 1 v2 | 0.0674 | 0.3454 | 4.3249 | 0.0000 | 4.9782 | 0 |
| Candidate 2 (3 seeds) | 0.0934 | 0.4759 | 4.4463 | 0.0000 | 5.4009 | 0 |

## T3: main comparisons, cohort 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Best tuned fixed (B1) | 0.3110 | 1.5365 | 5.4618 | 0.0042 | 12.0202 | 0 |
| Estimate + stiffness rule | 0.2710 | 1.6107 | 5.5094 | 0.0020 | 10.7190 | 0 |
| Estimate + Haninger damping rule | 0.2976 | 1.5755 | 5.4808 | 0.0000 | 7.2436 | 0 |
| RLS estimate + rule | 0.5556 | 5.2244 | 9.0959 | 0.0350 | 43.1995 | 0 |
| Duan (RAS 2018), adapted | 1.0542 | 4.5938 | 8.5449 | 0.0996 | 22.5535 | 0 |
| Song (RA-L 2025), adapted | 0.4385 | 1.9457 | 5.8905 | 0.0060 | 12.9494 | 0 |
| Jin (TMECH 2023), one-step | 0.3099 | 1.5735 | 5.5014 | 0.0088 | 13.8536 | 0 |
| Jin (TMECH 2023), multistep | 0.2848 | 1.4015 | 5.3295 | 0.0000 | 7.1223 | 0 |
| Xue, repaired (initial-surface prior) | 0.5806 | 2.2962 | 6.0719 | 0.0618 | 18.0851 | 0 |
| Xue, repaired (nominal surface) | 0.5276 | 2.0669 | 6.0137 | 0.0220 | 16.3626 | 0 |
| Learned VIC, repaired/retrained (3 seeds) | 0.2783 | 1.3578 | 5.2923 | 0.0115 | 20.7932 | 0 |
| Candidate 1 v1 | 0.1633 | 1.1804 | 5.0847 | 0.0000 | 7.1215 | 0 |
| Candidate 1 v2 | 0.1627 | 1.1865 | 5.0912 | 0.0000 | 7.1222 | 0 |
| Candidate 2 (3 seeds) | 0.2260 | 1.3582 | 5.2411 | 0.0000 | 7.9067 | 0 |

## Before/after repair on the same test scenes

These test outcomes were computed after parameter selection. They did not feed the search.

| Cohort | Task | Baseline | RMSE before → after (N) | Mean max error before → after (N) | Mean peak before → after (N) | Zero-force time before → after (s) |
|---|---|---|---:|---:|---:|---:|
| 3000 | T1 | Duan (RAS 2018), adapted | 0.6236 → 0.6337 | 2.8784 → 2.9808 | 6.8702 → 6.9771 | 0.0344 → 0.0144 |
| 3000 | T1 | Jin (TMECH 2023), multistep | 0.1396 → 0.1359 | 0.7251 → 0.6864 | 4.7154 → 4.6790 | 0.0000 → 0.0000 |
| 3000 | T1 | Jin (TMECH 2023), one-step | 0.1489 → 0.1264 | 0.8888 → 0.6668 | 4.8820 → 4.6598 | 0.0024 → 0.0000 |
| 3000 | T1 | Learned VIC, repaired/retrained (3 seeds) | 0.1837 → 0.1235 | 1.0675 → 0.6176 | 5.0532 → 4.6093 | 0.0086 → 0.0013 |
| 3000 | T1 | Song (RA-L 2025), adapted | 0.3219 → 0.2110 | 1.2124 → 0.8574 | 5.1954 → 4.8489 | 0.0024 → 0.0000 |
| 3000 | T1 | Xue, repaired (initial-surface prior) | 0.5472 → 0.3514 | 2.1467 → 1.3875 | 5.9855 → 5.3568 | 0.0000 → 0.0000 |
| 3000 | T1 | Xue, repaired (nominal surface) | 0.5561 → 0.3419 | 2.2232 → 1.3561 | 6.0358 → 5.3151 | 0.0000 → 0.0000 |
| 3000 | T3 | Duan (RAS 2018), adapted | 1.0546 → 1.0542 | 4.5967 → 4.5938 | 8.5478 → 8.5449 | 0.0998 → 0.0996 |
| 3000 | T3 | Jin (TMECH 2023), multistep | 0.3490 → 0.2848 | 1.7046 → 1.4015 | 5.5597 → 5.3295 | 0.0030 → 0.0000 |
| 3000 | T3 | Jin (TMECH 2023), one-step | 0.3156 → 0.3099 | 1.5985 → 1.5735 | 5.4971 → 5.5014 | 0.0032 → 0.0088 |
| 3000 | T3 | Learned VIC, repaired/retrained (3 seeds) | 0.3179 → 0.2783 | 1.4615 → 1.3578 | 5.4004 → 5.2923 | 0.0011 → 0.0115 |
| 3000 | T3 | Song (RA-L 2025), adapted | 0.4396 → 0.4385 | 1.9211 → 1.9457 | 5.8644 → 5.8905 | 0.0066 → 0.0060 |
| 3000 | T3 | Xue, repaired (initial-surface prior) | 1.2762 → 0.5806 | 4.3995 → 2.2962 | 8.2126 → 6.0719 | 0.2082 → 0.0618 |
| 3000 | T3 | Xue, repaired (nominal surface) | 1.2593 → 0.5276 | 4.6191 → 2.0669 | 8.4276 → 6.0137 | 0.1350 → 0.0220 |
| 4000 | T1 | Duan (RAS 2018), adapted | 0.5927 → 0.5654 | 2.8776 → 2.5971 | 6.8664 → 6.5930 | 0.0344 → 0.0122 |
| 4000 | T1 | Jin (TMECH 2023), multistep | 0.1302 → 0.1247 | 0.6307 → 0.6057 | 4.6132 → 4.5911 | 0.0000 → 0.0000 |
| 4000 | T1 | Jin (TMECH 2023), one-step | 0.1230 → 0.1161 | 0.6492 → 0.5920 | 4.6346 → 4.5765 | 0.0000 → 0.0000 |
| 4000 | T1 | Learned VIC, repaired/retrained (3 seeds) | 0.1878 → 0.1179 | 1.1144 → 0.5937 | 5.0963 → 4.5783 | 0.0093 → 0.0007 |
| 4000 | T1 | Song (RA-L 2025), adapted | 0.2892 → 0.1963 | 1.0662 → 0.7649 | 5.0397 → 4.7538 | 0.0000 → 0.0000 |
| 4000 | T1 | Xue, repaired (initial-surface prior) | 0.5147 → 0.3281 | 2.0301 → 1.2853 | 5.8833 → 5.2338 | 0.0000 → 0.0000 |
| 4000 | T1 | Xue, repaired (nominal surface) | 0.5165 → 0.3314 | 2.0010 → 1.2431 | 5.7827 → 5.1973 | 0.0000 → 0.0000 |
| 4000 | T3 | Duan (RAS 2018), adapted | 1.1371 → 1.1370 | 4.8019 → 4.8073 | 8.7045 → 8.7098 | 0.1726 → 0.1730 |
| 4000 | T3 | Jin (TMECH 2023), multistep | 0.3523 → 0.2819 | 1.7903 → 1.4494 | 5.6591 → 5.3714 | 0.0004 → 0.0000 |
| 4000 | T3 | Jin (TMECH 2023), one-step | 0.3297 → 0.2915 | 1.8181 → 1.5274 | 5.6944 → 5.4442 | 0.0038 → 0.0004 |
| 4000 | T3 | Learned VIC, repaired/retrained (3 seeds) | 0.3185 → 0.3096 | 1.5964 → 1.4300 | 5.4604 → 5.3698 | 0.0005 → 0.0158 |
| 4000 | T3 | Song (RA-L 2025), adapted | 0.4374 → 0.4317 | 2.0736 → 2.0551 | 5.9655 → 5.9507 | 0.0052 → 0.0030 |
| 4000 | T3 | Xue, repaired (initial-surface prior) | 1.2699 → 0.6939 | 4.5187 → 2.5474 | 8.3413 → 6.3557 | 0.2142 → 0.1142 |
| 4000 | T3 | Xue, repaired (nominal surface) | 1.2359 → 0.5997 | 4.4236 → 2.3490 | 8.2952 → 6.2597 | 0.1842 → 0.0484 |

## Paired uncertainty for proposed controllers

Scene-paired bootstrap, 10,000 resamples, three neural seeds averaged within scene. Intervals condition on the saved fitted models and are marginal/exploratory, without multiple-comparison adjustment. Training-seed variation is in training_seed_summary.csv. Negative differences favor the proposed controller. Full metric/reference pairs are in paired_differences.csv.

| Cohort | Task | Proposed | Reference | RMSE difference (N) | Marginal 95% interval (N) |
|---|---|---|---|---:|---:|
| 3000 | T1 | Candidate 1 v2 | Best tuned fixed (B1) | -0.0880 | [-0.1501, -0.0473] |
| 3000 | T1 | Candidate 1 v2 | Duan (RAS 2018), adapted | -0.5663 | [-0.6984, -0.4431] |
| 3000 | T1 | Candidate 1 v2 | Estimate + stiffness rule | -0.0575 | [-0.0722, -0.0458] |
| 3000 | T1 | Candidate 1 v2 | Estimate + Haninger damping rule | -0.0917 | [-0.1313, -0.0637] |
| 3000 | T1 | Candidate 1 v2 | RLS estimate + rule | -0.1868 | [-0.2110, -0.1640] |
| 3000 | T1 | Candidate 1 v2 | Jin (TMECH 2023), multistep | -0.0685 | [-0.0929, -0.0465] |
| 3000 | T1 | Candidate 1 v2 | Jin (TMECH 2023), one-step | -0.0590 | [-0.0775, -0.0416] |
| 3000 | T1 | Candidate 1 v2 | Learned VIC, repaired/retrained (3 seeds) | -0.0561 | [-0.0841, -0.0357] |
| 3000 | T1 | Candidate 1 v2 | Song (RA-L 2025), adapted | -0.1436 | [-0.1733, -0.1145] |
| 3000 | T1 | Candidate 1 v2 | Xue, repaired (initial-surface prior) | -0.2840 | [-0.3267, -0.2423] |
| 3000 | T1 | Candidate 1 v2 | Xue, repaired (nominal surface) | -0.2745 | [-0.3186, -0.2329] |
| 3000 | T1 | Candidate 2 (3 seeds) | Best tuned fixed (B1) | -0.0620 | [-0.1249, -0.0232] |
| 3000 | T1 | Candidate 2 (3 seeds) | Duan (RAS 2018), adapted | -0.5403 | [-0.6710, -0.4159] |
| 3000 | T1 | Candidate 2 (3 seeds) | Estimate + stiffness rule | -0.0315 | [-0.0453, -0.0209] |
| 3000 | T1 | Candidate 2 (3 seeds) | Estimate + Haninger damping rule | -0.0657 | [-0.1037, -0.0396] |
| 3000 | T1 | Candidate 2 (3 seeds) | RLS estimate + rule | -0.1608 | [-0.1834, -0.1392] |
| 3000 | T1 | Candidate 2 (3 seeds) | Jin (TMECH 2023), multistep | -0.0425 | [-0.0638, -0.0225] |
| 3000 | T1 | Candidate 2 (3 seeds) | Jin (TMECH 2023), one-step | -0.0330 | [-0.0491, -0.0179] |
| 3000 | T1 | Candidate 2 (3 seeds) | Learned VIC, repaired/retrained (3 seeds) | -0.0300 | [-0.0573, -0.0117] |
| 3000 | T1 | Candidate 2 (3 seeds) | Song (RA-L 2025), adapted | -0.1175 | [-0.1456, -0.0900] |
| 3000 | T1 | Candidate 2 (3 seeds) | Xue, repaired (initial-surface prior) | -0.2580 | [-0.2996, -0.2178] |
| 3000 | T1 | Candidate 2 (3 seeds) | Xue, repaired (nominal surface) | -0.2485 | [-0.2896, -0.2090] |
| 3000 | T3 | Candidate 1 v2 | Best tuned fixed (B1) | -0.1482 | [-0.2011, -0.1099] |
| 3000 | T3 | Candidate 1 v2 | Duan (RAS 2018), adapted | -0.8915 | [-1.0590, -0.7572] |
| 3000 | T3 | Candidate 1 v2 | Estimate + stiffness rule | -0.1083 | [-0.1346, -0.0853] |
| 3000 | T3 | Candidate 1 v2 | Estimate + Haninger damping rule | -0.1349 | [-0.1484, -0.1217] |
| 3000 | T3 | Candidate 1 v2 | RLS estimate + rule | -0.3928 | [-0.4923, -0.3080] |
| 3000 | T3 | Candidate 1 v2 | Jin (TMECH 2023), multistep | -0.1221 | [-0.1441, -0.1017] |
| 3000 | T3 | Candidate 1 v2 | Jin (TMECH 2023), one-step | -0.1472 | [-0.2021, -0.1096] |
| 3000 | T3 | Candidate 1 v2 | Learned VIC, repaired/retrained (3 seeds) | -0.1155 | [-0.2010, -0.0640] |
| 3000 | T3 | Candidate 1 v2 | Song (RA-L 2025), adapted | -0.2758 | [-0.3311, -0.2341] |
| 3000 | T3 | Candidate 1 v2 | Xue, repaired (initial-surface prior) | -0.4179 | [-0.5880, -0.2861] |
| 3000 | T3 | Candidate 1 v2 | Xue, repaired (nominal surface) | -0.3649 | [-0.4621, -0.2999] |
| 3000 | T3 | Candidate 2 (3 seeds) | Best tuned fixed (B1) | -0.0849 | [-0.1366, -0.0473] |
| 3000 | T3 | Candidate 2 (3 seeds) | Duan (RAS 2018), adapted | -0.8281 | [-0.9915, -0.6983] |
| 3000 | T3 | Candidate 2 (3 seeds) | Estimate + stiffness rule | -0.0450 | [-0.0660, -0.0262] |
| 3000 | T3 | Candidate 2 (3 seeds) | Estimate + Haninger damping rule | -0.0716 | [-0.0810, -0.0619] |
| 3000 | T3 | Candidate 2 (3 seeds) | RLS estimate + rule | -0.3295 | [-0.4235, -0.2490] |
| 3000 | T3 | Candidate 2 (3 seeds) | Jin (TMECH 2023), multistep | -0.0587 | [-0.0801, -0.0366] |
| 3000 | T3 | Candidate 2 (3 seeds) | Jin (TMECH 2023), one-step | -0.0839 | [-0.1374, -0.0464] |
| 3000 | T3 | Candidate 2 (3 seeds) | Learned VIC, repaired/retrained (3 seeds) | -0.0522 | [-0.1351, -0.0017] |
| 3000 | T3 | Candidate 2 (3 seeds) | Song (RA-L 2025), adapted | -0.2125 | [-0.2645, -0.1725] |
| 3000 | T3 | Candidate 2 (3 seeds) | Xue, repaired (initial-surface prior) | -0.3546 | [-0.5214, -0.2280] |
| 3000 | T3 | Candidate 2 (3 seeds) | Xue, repaired (nominal surface) | -0.3016 | [-0.3969, -0.2415] |
| 4000 | T1 | Candidate 1 v2 | Best tuned fixed (B1) | -0.0800 | [-0.1343, -0.0419] |
| 4000 | T1 | Candidate 1 v2 | Duan (RAS 2018), adapted | -0.4982 | [-0.6291, -0.3813] |
| 4000 | T1 | Candidate 1 v2 | Estimate + stiffness rule | -0.0531 | [-0.0634, -0.0441] |
| 4000 | T1 | Candidate 1 v2 | Estimate + Haninger damping rule | -0.0729 | [-0.0918, -0.0571] |
| 4000 | T1 | Candidate 1 v2 | RLS estimate + rule | -0.1774 | [-0.2020, -0.1541] |
| 4000 | T1 | Candidate 1 v2 | Jin (TMECH 2023), multistep | -0.0575 | [-0.0780, -0.0393] |
| 4000 | T1 | Candidate 1 v2 | Jin (TMECH 2023), one-step | -0.0488 | [-0.0635, -0.0349] |
| 4000 | T1 | Candidate 1 v2 | Learned VIC, repaired/retrained (3 seeds) | -0.0506 | [-0.0739, -0.0333] |
| 4000 | T1 | Candidate 1 v2 | Song (RA-L 2025), adapted | -0.1291 | [-0.1562, -0.1038] |
| 4000 | T1 | Candidate 1 v2 | Xue, repaired (initial-surface prior) | -0.2608 | [-0.3027, -0.2207] |
| 4000 | T1 | Candidate 1 v2 | Xue, repaired (nominal surface) | -0.2641 | [-0.3085, -0.2208] |
| 4000 | T1 | Candidate 2 (3 seeds) | Best tuned fixed (B1) | -0.0576 | [-0.1104, -0.0206] |
| 4000 | T1 | Candidate 2 (3 seeds) | Duan (RAS 2018), adapted | -0.4758 | [-0.5959, -0.3576] |
| 4000 | T1 | Candidate 2 (3 seeds) | Estimate + stiffness rule | -0.0307 | [-0.0393, -0.0232] |
| 4000 | T1 | Candidate 2 (3 seeds) | Estimate + Haninger damping rule | -0.0505 | [-0.0678, -0.0361] |
| 4000 | T1 | Candidate 2 (3 seeds) | RLS estimate + rule | -0.1551 | [-0.1775, -0.1342] |
| 4000 | T1 | Candidate 2 (3 seeds) | Jin (TMECH 2023), multistep | -0.0351 | [-0.0539, -0.0186] |
| 4000 | T1 | Candidate 2 (3 seeds) | Jin (TMECH 2023), one-step | -0.0265 | [-0.0391, -0.0141] |
| 4000 | T1 | Candidate 2 (3 seeds) | Learned VIC, repaired/retrained (3 seeds) | -0.0282 | [-0.0503, -0.0120] |
| 4000 | T1 | Candidate 2 (3 seeds) | Song (RA-L 2025), adapted | -0.1067 | [-0.1323, -0.0828] |
| 4000 | T1 | Candidate 2 (3 seeds) | Xue, repaired (initial-surface prior) | -0.2384 | [-0.2796, -0.1992] |
| 4000 | T1 | Candidate 2 (3 seeds) | Xue, repaired (nominal surface) | -0.2418 | [-0.2850, -0.2003] |
| 4000 | T3 | Candidate 1 v2 | Best tuned fixed (B1) | -0.1510 | [-0.2235, -0.0971] |
| 4000 | T3 | Candidate 1 v2 | Duan (RAS 2018), adapted | -0.9644 | [-1.1807, -0.7840] |
| 4000 | T3 | Candidate 1 v2 | Estimate + stiffness rule | -0.1000 | [-0.1248, -0.0798] |
| 4000 | T3 | Candidate 1 v2 | Estimate + Haninger damping rule | -0.1223 | [-0.1353, -0.1094] |
| 4000 | T3 | Candidate 1 v2 | RLS estimate + rule | -0.3626 | [-0.4192, -0.3081] |
| 4000 | T3 | Candidate 1 v2 | Jin (TMECH 2023), multistep | -0.1093 | [-0.1304, -0.0869] |
| 4000 | T3 | Candidate 1 v2 | Jin (TMECH 2023), one-step | -0.1190 | [-0.1423, -0.0950] |
| 4000 | T3 | Candidate 1 v2 | Learned VIC, repaired/retrained (3 seeds) | -0.1371 | [-0.2421, -0.0624] |
| 4000 | T3 | Candidate 1 v2 | Song (RA-L 2025), adapted | -0.2591 | [-0.3128, -0.2128] |
| 4000 | T3 | Candidate 1 v2 | Xue, repaired (initial-surface prior) | -0.5214 | [-0.7188, -0.3563] |
| 4000 | T3 | Candidate 1 v2 | Xue, repaired (nominal surface) | -0.4272 | [-0.5642, -0.3167] |
| 4000 | T3 | Candidate 2 (3 seeds) | Best tuned fixed (B1) | -0.0928 | [-0.1614, -0.0415] |
| 4000 | T3 | Candidate 2 (3 seeds) | Duan (RAS 2018), adapted | -0.9062 | [-1.1195, -0.7348] |
| 4000 | T3 | Candidate 2 (3 seeds) | Estimate + stiffness rule | -0.0418 | [-0.0677, -0.0216] |
| 4000 | T3 | Candidate 2 (3 seeds) | Estimate + Haninger damping rule | -0.0641 | [-0.0737, -0.0538] |
| 4000 | T3 | Candidate 2 (3 seeds) | RLS estimate + rule | -0.3044 | [-0.3557, -0.2538] |
| 4000 | T3 | Candidate 2 (3 seeds) | Jin (TMECH 2023), multistep | -0.0511 | [-0.0728, -0.0294] |
| 4000 | T3 | Candidate 2 (3 seeds) | Jin (TMECH 2023), one-step | -0.0608 | [-0.0827, -0.0386] |
| 4000 | T3 | Candidate 2 (3 seeds) | Learned VIC, repaired/retrained (3 seeds) | -0.0789 | [-0.1816, -0.0064] |
| 4000 | T3 | Candidate 2 (3 seeds) | Song (RA-L 2025), adapted | -0.2009 | [-0.2513, -0.1589] |
| 4000 | T3 | Candidate 2 (3 seeds) | Xue, repaired (initial-surface prior) | -0.4632 | [-0.6554, -0.2984] |
| 4000 | T3 | Candidate 2 (3 seeds) | Xue, repaired (nominal surface) | -0.3690 | [-0.5050, -0.2636] |

## Search and training budget

Screen on 8 training scenes, reevaluate shortlisted configurations on all 20 training scenes, then select using 20 development scenes. Two subsequent coordinate-refinement rounds include the incumbent and reevaluate all configurations on 20 training scenes; development selects among the top three and the incumbent. Tests are not used. The inherited scalar J selects parameters; the tables above expose its physical components rather than presenting J as a physical measurement.

| Task | Family | Coarse configurations | Local configurations evaluated | Train J | Dev J | Train peak criterion | Dev peak criterion |
|---|---|---:|---:|---:|---:|---|---|
| T3 | Xue | 96 | 28 | 2.40233 | 1.49340 | False | True |
| T3 | Xue_nominal | 96 | 32 | 1.58233 | 1.07879 | True | True |
| T3 | Jin_OS | 84 | 34 | 0.36283 | 0.34793 | True | True |
| T3 | Jin_MS | 84 | 32 | 0.38030 | 0.34233 | True | True |
| T3 | Duan | 96 | 13 | 5.64603 | 5.06574 | False | False |
| T3 | Song | 96 | 16 | 0.79032 | 0.72362 | True | True |
| T1 | Xue | 96 | 32 | 0.72257 | 0.56919 | True | True |
| T1 | Xue_nominal | 96 | 32 | 0.70684 | 0.55003 | True | True |
| T1 | Jin_OS | 84 | 34 | 0.10948 | 0.10669 | True | True |
| T1 | Jin_MS | 84 | 32 | 0.11807 | 0.14477 | True | True |
| T1 | Duan | 96 | 14 | 2.82747 | 2.54873 | False | False |
| T1 | Song | 96 | 18 | 0.28874 | 0.21227 | True | True |

| Task | Training seed | Generations | Training/evaluation episodes | Best generation | Train J | Dev J |
|---|---:|---:|---:|---:|---:|---:|
| T3 | 0 | 100 | 10820 | 80 | 0.60519 | 0.21017 |
| T3 | 1 | 100 | 10800 | 80 | 0.47545 | 0.23750 |
| T3 | 2 | 100 | 10800 | 90 | 0.38857 | 0.25816 |
| T1 | 0 | 100 | 10800 | 80 | 0.08447 | 0.08588 |
| T1 | 1 | 100 | 10800 | 80 | 0.07340 | 0.07304 |
| T1 | 2 | 100 | 10800 | 100 | 0.08261 | 0.09353 |

Learning: population 16; six common training scenes per generation; initial ES sigma 0.025; network 12–8–3, with learned constant gain offsets; 100 generations per task/seed. Train/dev checkpoint histories and all final parameter vectors are in learning.json. A best generation of zero means the finite training budget did not improve on B1 under the declared selection rule. This is not reported as learning convergence. After training, each selected network receives a 15-setting train/dev calibration of output filtering (tau = 0.016/0.04/0.08 s) and reference span (0.25/0.5/1/3/10 mm), including the incumbent. These evaluation episodes are included above; weights are unchanged by calibration.

## Selected parameter vectors

- T3_Xue: `{"Delta": 0.03, "KPF": 0.0003, "M": 3.0, "dmax": 50.0, "dmin": 50.0, "efb": 40.0, "est_gain": 0.2, "f0src": "model", "hp": 10, "ki": 30.0, "kmax": 1000.0, "n0": 1.3, "reanchor": false, "servo_lead": 0.5, "state_reset": true, "ys0_mode": "A1"}`
- T3_Xue_nominal: `{"Delta": 0.1, "KPF": 0.0001, "M": 2.0, "dmax": 100.0, "dmin": 100.0, "efb": 40.0, "est_gain": 0.1, "f0src": "measured", "hp": 10, "ki": 10.0, "kmax": 15.0, "n0": 1.0, "reanchor": true, "servo_lead": 1.0, "state_reset": true, "ys0_mode": "nominal"}`
- T3_Jin_OS: `{"D0": 180.5354225770494, "K": 320.43248977694043, "M": 1.0004041935912389, "dist": 1, "kgain": 1.0, "ki": 29.999999999999996, "kpf0": 0.00010492401101884874, "lookahead": 1, "mode": "os", "qv": 0.0, "rc": 10.0, "rd": 0.025, "update": 0.008}`
- T3_Jin_MS: `{"D0": 180.5354225770494, "K": 320.43248977694043, "M": 1.0004041935912389, "dist": 1, "kgain": 1.0, "ki": 29.999999999999996, "kpf0": 0.00010492401101884874, "mode": "ms", "qv": 0.0, "rc": 100.0, "rd": 0.025, "tau_k": 0.3, "update": 0.008}`
- T3_Duan: `{"b": 288.4068840731621, "m": 0.25, "sigma": 0.08212337307347552}`
- T3_Song: `{"D": 197.52320565324726, "KPF": 6.815353786161588e-05, "M": 0.25, "ki": 18.49663884218113}`
- T1_Xue: `{"M": 1.0, "dmax": 200.0, "dmin": 100.0, "efb": 10.0, "f0src": "measured", "ki": 6.0, "kmax": 30.0, "servo_lead": 0.5, "state_reset": false, "ys0_mode": "A1"}`
- T1_Xue_nominal: `{"M": 1.0, "dmax": 200.0, "dmin": 100.0, "efb": 10.0, "f0src": "measured", "ki": 6.0, "kmax": 30.0, "servo_lead": 0.5, "state_reset": false, "ys0_mode": "nominal"}`
- T1_Jin_OS: `{"D0": 104.27236953871011, "K": 0.00293892907502169, "M": 1.6882664267023355, "dist": 1, "kgain": 1.7, "ki": 11.324934016090703, "kpf0": 0.00011638706231124378, "lookahead": 1, "mode": "os", "qv": 0.1, "rc": 10.0, "rd": 0.1, "update": 0.008}`
- T1_Jin_MS: `{"D0": 104.27236953871011, "K": 0.0007347322687554225, "M": 1.6882664267023355, "dist": 1, "kgain": 1.0, "ki": 11.324934016090703, "kpf0": 0.00011638706231124378, "mode": "ms", "qv": 0.1, "rc": 100.0, "rd": 0.25, "tau_k": 0.1, "update": 0.008}`
- T1_Duan: `{"b": 490.85779472977504, "m": 0.25, "sigma": 0.09208401731464654}`
- T1_Song: `{"D": 137.55146364701568, "KPF": 0.00018749470517734468, "M": 1.4587871947392728, "ki": 6.59471354166848}`

## Fixed extremes, unanchored diagnostic and truth references

These rows are unchanged historical controls/reference bounds; they were not substituted for repaired competitive baselines. Truth references use privileged online stiffness and are not deployable baselines.

### T1, cohort 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 15.3557 | 147.0219 | 151.0219 | 5.4516 | 269.7349 | 0 |
| Fixed high damping | 2.9783 | 9.9165 | 13.8014 | 1.1018 | 27.2236 | 0 |
| Jin multistep, unanchored diagnostic | 0.5888 | 4.1478 | 8.1270 | 0.0464 | 16.3158 | 0 |
| Truth-stiffness reference, 0 s delay | 0.0658 | 0.3196 | 4.3052 | 0.0000 | 5.4954 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.0996 | 0.4753 | 4.3873 | 0.0000 | 6.1857 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.1188 | 0.7081 | 4.6483 | 0.0000 | 8.5450 | 0 |

### T3, cohort 4000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 18.0291 | 149.2457 | 153.2457 | 8.2284 | 306.6901 | 0 |
| Fixed high damping | 6.2277 | 25.0353 | 29.0353 | 4.4138 | 52.3818 | 0 |
| Jin multistep, unanchored diagnostic | 1.8048 | 14.3591 | 18.2691 | 0.4532 | 55.3832 | 0 |
| Truth-stiffness reference, 0 s delay | 0.1524 | 0.9506 | 4.8061 | 0.0000 | 6.1974 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1654 | 0.9869 | 4.8685 | 0.0000 | 6.4316 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.3327 | 1.6563 | 5.4878 | 0.0018 | 8.2387 | 0 |

### T1, cohort 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | 15.4310 | 152.2074 | 156.2074 | 5.3912 | 261.1714 | 0 |
| Fixed high damping | 3.1345 | 11.0831 | 15.0244 | 1.1708 | 30.1981 | 0 |
| Jin multistep, unanchored diagnostic | 0.7074 | 5.3513 | 9.2668 | 0.0570 | 22.8351 | 0 |
| Truth-stiffness reference, 0 s delay | 0.0674 | 0.3434 | 4.3349 | 0.0000 | 5.1022 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1065 | 0.5170 | 4.4413 | 0.0000 | 5.8104 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.1273 | 0.7465 | 4.7042 | 0.0000 | 7.6770 | 0 |

### T3, cohort 3000

| Controller | Force RMSE (N) | Mean max absolute error (N) | Mean peak force (N) | Mean zero-force time (s) | Worst sampled peak (N) | Failed runs |
|---|---:|---:|---:|---:|---:|---:|
| Fixed low damping | NA | NA | NA | NA | NA | 1 |
| Fixed high damping | 5.8100 | 23.0524 | 27.0524 | 4.1216 | 66.4453 | 0 |
| Jin multistep, unanchored diagnostic | 1.6606 | 12.5960 | 16.5430 | 0.3668 | 58.0744 | 0 |
| Truth-stiffness reference, 0 s delay | 0.1452 | 0.8919 | 4.7769 | 0.0000 | 6.5007 | 0 |
| Truth-stiffness reference, 0.1 s delay | 0.1538 | 0.8966 | 4.7809 | 0.0000 | 6.4518 | 0 |
| Truth-stiffness reference, 0.2 s delay | 0.3041 | 1.4630 | 5.3045 | 0.0004 | 7.8143 | 0 |

## Reproduction

```powershell
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py checks
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py tune --workers 6
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py polish --workers 6
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py learn --workers 6
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py calibrate --workers 6
python -X utf8 experiments/cobot_scan_benchmark/refine_baselines.py test --workers 6
python -X utf8 experiments/cobot_scan_benchmark/report_refined_baselines.py
```

Completed search/training entries are reused. protocol.json freezes selection, weights, source files and cached reference data before testing; changed hashes stop the test. Runtime samples are diagnostic on this workstation and are not real-time certification; old and new concurrent-process timing conditions differ. No task/surface/controller values in the frozen B0 snapshot were edited.
