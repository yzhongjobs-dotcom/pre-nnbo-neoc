# C1 repair on the original three scanning tasks

The user stopped additional simple-controller tuning and scenario expansion. Only the selected C1 repair and one estimator ablation were evaluated on the original seed-7200 tasks (40 paired scenes each). Historical controllers, weights and results remain unchanged. No neural training was run.

## Method and capability tested

The repair profiles a cubic force-drift polynomial out of the stiffness regression before updating stiffness. It retains the original 128 ms predictive controller, uses three local solve iterations, and replaces the task-tuned fixed recovery with stiffness-based recovery and command-state compensation. The estimated drift is an equivalent contact response, not established anatomical motion.

The tested capability is maintaining contact and limiting impact during sustained push/retract on stiff tissue, while retaining useful combined-scan performance. This is not a first-joint-estimation or near-optimal-control claim.

## Three-task outcomes

| Method | T1 stiffness | T2 soft motion | T2 stiff motion | T3 combined |
|---|---:|---:|---:|---:|
| Jin_MS | 0.117 | 0.265 | 1.268 | 0.285 |
| Xue_nominal | 0.276 | 0.383 | 3.414 | 0.511 |
| Learned_VIC | 0.113 | 0.188 | 6.132 | 0.269 |
| Candidate1_v2 | 0.060 | 0.057 | 1.372 | 0.147 |
| C1_updated | 0.087 | 0.068 | 0.291 | 0.142 |
| Candidate2 | 0.083 | 0.143 | 0.493 | 0.210 |

RMSE in N. Literature rows are the previously evaluated task adaptations. T2 was used to develop this repair, whereas those historical rows retain their T3-selected settings; this is not an equally retuned ranking or a claim about the original papers' performance limits.

## Stiff push/retract and estimator ablation

| Configuration | RMSE N | Worst sampled peak N | Scenes below 0.5 N | 2 s event recovery |
|---|---:|---:|---:|---:|
| Candidate1_v2 | 1.372 | 22.228 | 6/20 | 80.0% |
| C1_old_model_new_recovery | 0.531 | 13.230 | 6/20 | 100.0% |
| C1_updated | 0.291 | 6.779 | 0/20 | 100.0% |

The ablation retains the new recovery and three solve iterations but restores the original joint nonlinear estimator. Its remaining stiff-motion error and high-force tail support an additional role for the estimator repair. The experiment does not separately identify every contribution of iteration count and recovery.

The stiff-motion repair costs static accuracy: T1 RMSE increases from 0.060 to 0.087 N; the worst sampled T1 peak increases from 5.327 to 7.584 N. T3 mean RMSE remains similar (0.147 to 0.142 N). The new method is not uniformly better and has no global 7 N guarantee.

## Paired evidence

| Task/cell | Reference | Mean difference N | Family 95% interval N |
|---|---|---:|---|
| T2_motion/hard | Candidate1_v2 | -1.0811 | [-2.0543, -0.2470] |
| T2_motion/hard | EtR-A-Haninger | -0.1424 | [-0.1782, -0.1136] |
| T1_static/all | Candidate1_v2 | +0.0270 | [+0.0178, +0.0423] |
| T3_combined/all | Candidate1_v2 | -0.0047 | [-0.0139, +0.0026] |
| T2_motion/hard | C1_old_model_new_recovery | -0.2405 | [-0.4216, -0.0887] |

Intervals refer to new minus reference, using paired scene resampling and Bonferroni coverage for the five comparisons above. They do not remove the different development exposure noted above.

## Update time and reproducibility

| Task | P99 ms | Maximum measured ms |
|---|---:|---:|
| T1_static | 0.501 | 2.450 |
| T2_motion | 0.628 | 2.009 |
| T3_combined | 0.368 | 2.137 |

Single-process policy call including observation assembly and output projection, excluding separate JIT warmup, the first 10 calls of each scene, the simulated plant and physical communications. The control period is 8 ms; this is not robot real-time certification.

Reproduce checks and the frozen cohort from the repository root:
```powershell
python -X utf8 experiments/cobot_scan_benchmark/evaluate_contact_update.py check
python -X utf8 experiments/cobot_scan_benchmark/evaluate_contact_update.py evaluate
python -X utf8 experiments/cobot_scan_benchmark/analyze_contact_update.py
```

Frozen configuration and source hashes: [frozen.json](frozen.json). Per-scene outcomes: [test_rows.json](test_rows.json). Raw trajectories, commands and update times: [traces](traces/). Estimator and truth-input checks: [checks.json](checks.json). Development/validation data refer only to C1 and predate test evaluation.

The closest DOSMAC, Cao and Wang comparisons remain unfinished; exploratory unverified adaptations are not included as publication evidence. Closed-loop stability, physical motion recovery and image-quality benefits are not established by this repair.
