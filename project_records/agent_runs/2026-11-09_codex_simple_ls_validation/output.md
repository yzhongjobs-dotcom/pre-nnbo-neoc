# Final output

Implemented two ordinary least-squares estimators using the existing compensated/H-preprocessed regression signals and unchanged 24-tick H host, phase estimator, gates, and controller. Simple-2's intercept uses the `common._build()` filtered constant basis. Exact synthetic coefficient checks passed.

The retained result reuses the original corrected-comparison rows for all 17 methods and appends six new closed-loop episodes (A, B, combo; `neural_white`, 28°, `i=0`). The stiffness diagnostics compare raw OLS solutions with current `k_eff` in `log[:,26]`.

Across all three fixed cases, both simple baselines had higher (L), force RMSE, and peak force than the original methods. Simple-2 was less poor than Simple-1 on force error, but raw OLS (k) underestimated `k_eff` in 90.1–95.6% of valid ticks and still produced large task error. This is limited evidence for these three paired instances.

Execution deviation: the six cases were run three times (18 simulation attempts total). The first pass failed during JSON serialization after all six episodes; the second pass exposed that `log[:,21]` was not current `k_eff`; the third pass used `log[:,26]` and is the retained result. The duplicate passes used no new cases, seeds, or parameter changes. This should be considered when assessing compliance with the requested six-run budget.

Artifacts: `experiments/rossi_failure_study/simple_ls_validation.py`, `experiments/rossi_failure_study/results/simple_ls_validation.json`.
