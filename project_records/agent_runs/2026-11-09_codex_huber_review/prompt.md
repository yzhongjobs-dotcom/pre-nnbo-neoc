Independent review. You did not write this work. Find what is wrong; do not defend it.

Project: robot ultrasound constant-force scanning; window estimator z = kY + d + e with a prior on d. Question studied: Huber-D4 was expected to beat D4 on sudden disturbances, but the tables show no stable gain. Claimed root causes are in `research/huber_d4_theory_vs_experiment_20261109.md`.

Read that note, then these scripts in `experiments/rossi_failure_study/`: `huber_windows.py`, `huber_prior_check.py`, `huber_synthetic_map.py`, `huber_same_config_loop.py`, `huber_loop_sensitivity.py`, `huber_case_trace.py`; the estimator `d_robust.py`, `ridge_estimator.py`; results `results/huber_prior_check.json`, `huber_synthetic_map.json`, `huber_same_config_loop.json`, `huber_loop_sensitivity.json`. Window data is already in `results/huber_windows/`; do not regenerate it.

Check exactly these, nothing else:

1. Derivation. Section 1 of the note: the stationarity conditions and the influence-bound formula, with its factor rho*s = sigma*sqrt(rho) and the objective's factor 2 rho s^2. Is the "doesn't depend on step height" claim valid, and under what condition? Is the scale claim (bound proportional to c*sigma*sqrt(rho)/(kA)) right? Is "c near .72 gives about .30 median and 1.4 worst" right? Is the criterion delta*a_band >= c*s/G derived correctly, and why does activation occur near Q .3-.5 rather than 1?

2. Budget. In `huber_prior_check.py`, the telescoping error budget. Does k's part (weighted mean of k(t) with weights M Y) follow from the model? Are the row indices in `analyse` right (window entry to log row), including when a window is shorter than 24? Is the rest term fairly described?

3. Attribution claims. Do these follow from the data? (a) The plateau difference in the tables is entirely window length. (b) B's error is mainly calibration, not disturbance. (c) Case 26 cough is a real prior effect and the larger closed-loop gaps are not. Does the perturbation test (multiply k by 1+eps) really separate these? (d) In windows where Huber is much worse, the objective's own minimum is at k<0. Reproduce this for one window with your own short solver (for example at fixed k, minimize over d by a method other than the script's IRLS), not by rerunning the script.

4. Alternative explanations the note does not consider, if any are serious.

Constraints: no commit, no push. Do not edit research documents or existing scripts. Leave `00_基线快照_2026-10-14/` alone. Do not run closed-loop experiments beyond one or two episodes if a check needs them; no sweeps. Put any short check script in the scratchpad path below, not in the repo.
Scratchpad: D:/pip_env/temp/claude/8fbcb2674095e6be/4c5dc125-7da2-41f3-aa59-d2d3f2f2c904/scratchpad

Report: at most one and a half pages, in Chinese. For each of the four items: confirmed, or what is wrong, with the number or line that shows it. Say what you did not check. No suggestions for method improvements.
