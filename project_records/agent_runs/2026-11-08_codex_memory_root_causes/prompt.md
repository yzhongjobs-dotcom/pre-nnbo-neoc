# Task: test whether each memory method's weak result is a theoretical limit or a fixable defect

Simulation only. Do not read project_records/HANDOFF.md. Read `research/identification_module_scheme_20261108.md` section 2, `experiments/rossi_failure_study/memory_filters.py`, `memory_segment_eval.py` (its gates-off host, segment masks and scoring are to be reused unchanged), and the output of the last run in `project_records/agent_runs/2026-11-08_codex_gate_free_segments/output.md`.

A method may only be excluded for a theoretical limit that cannot be improved, not because one implementation did badly. For each method below the suspected root cause and a theory-based fix are given. Implement the fix as stated, check the predicted property, then measure. If you find a derivation wrong, say so and use the correct form.

## Common setting

ell = ln k. Window estimates have variance s_t^2 = (sigma/k)^2. N = s^2*T0 with T0 = 0.192 s. Parameters come from theory and are fixed, not tuned: tau_long = 2 s (plateau requirement 2 sd <= 0.09 at s = 0.2), tau_short = 0.2 s (optimum for ramps of 0.3-0.7 /s), so q = N/tau^2.

Innovations from overlapping windows are correlated, which violates the white-innovation assumption of every likelihood and of ARD. Use non-overlapping updates (hop = 24 ticks, R = s_t^2) for all methods here, holding the output between updates, so that this defect is removed from all of them. Note this adds up to 0.192 s latency and state it in the results.

## Methods

1. **Kalman, fixed q** (reference). q = N/tau^2 for tau_long and for tau_short. No Huber. Root cause: one memory length cannot serve plateau and transition. This is the theoretical limit of the fixed form.
2. **Two-mode IMM on q** (fix for 1). Interacting multiple model with q_lo = N/tau_long^2, q_hi = N/tau_short^2, Markov stay probability 0.98 per update, Gaussian likelihood on the correctly normalized innovation. The Method's motion tracker is already an IMM, so this is a mature mechanism in the same family.
3. **Lepski** Root causes: the old code intersects confidence intervals using each scale's own variance, which is conservative; selection is a hard switch. Fix: pairwise test between nested means with the correct difference variance, Var(mean_m - mean_N) = v_m - v_N for m < N with v the variance of the precision-weighted mean; choose the longest N with |mean_N - mean_m| <= c*sqrt(v_m - v_N) for all shorter m; c = 2.5. Scales 1, 2, 4, 8, 16 windows.
4. **Exponentially weighted aggregation** (soft version of 3, Vovk/Herbster-Warmuth fixed share). Run Kalman filters at tau in {0.2, 0.5, 1, 2, 4} s; weight each by exp(-sum of one-step log loss) with fixed share 0.02 per update; output the weighted mean and mixture variance.
5. **Innovation forgetting** Root causes found: stationary P becomes very small, so recovery after a change is slow because P grows only geometrically (delay ~ ln(P_target/P)/ln(1/lambda)); and the old normalization used the inflated R. Fix: use S = P + s_t^2 and keep a floor through a process noise q_floor = N/tau_long^2 added each update. Keep lam_max .998, lam_min .8, gain 1. Also run the same with the floor from distance, q_s*ds with q_s = q_floor/(0.192*v_ref) and v_ref = 20 mm/s, so that the same q_s applies at all speeds; this is the test of the scan-distance prior.
6. **Huang ARD** Root causes found: at single-window SNR the closed form alpha* = (nu^2 - S)_+ fires about 32% of the time, so memory stays short; the old precision learning had a scale error (see `project_records/agent_runs/2026-11-08_codex_memory_scheme_check/output.md`). Fixes: (a) implement eq (15)/(17) correctly in scalar form from that report; (b) block version: apply ARD to the mean of M = 4 consecutive windows, whose noise is s^2/M. Prediction: stationary memory grows about in proportion to M. The state filter stays replaced by the D4 window, the reference layer stays off. Say whether (a) or (b) changes the conclusion.

## Checks and measurements

A. **Synthetic check.** In `memory_theory_check.py` style (white window noise at hop 24, s = 0.2): stationary sd, ramp 0.29 and 0.73 /s lag and RMSE, step 0.3 and 1.1 delay. Compare to the predictions stated above, and to the exponential-memory values in the note's section 2.

B. **Closed loop, gates off, dev seeds only.** Reuse the gates-off host of `memory_segment_eval.py`, and its baseline check for the gates-on reference. Episodes as in that script: scene A nominal 36, scene A slow_thin 36, scene B 12. Report the same segment table: e_g median/robust sd, fraction k_good/k < 1 and < 0.85, force RMSE and P99, plus whole-episode L mean/median, max peak and number above 15 N. Include `off_none` as the reference. Add a lag measure for transitions: median of e_g on rising k versus falling k.

## Report

Under 2 pages in Chinese, final message only: one table per check; then for each method one line: root cause as found (theory, assumption, parameter or implementation), whether the fix worked, what remains, and whether any exclusion is justified on theoretical grounds or the evidence does not support exclusion. Report derivations you found wrong. State only what was measured. Do not select a winner on fresh seeds; none are used here.

Files: `experiments/rossi_failure_study/memory_filters2.py`, `memory_root_causes.py`, results `results/memory_root_causes.json`. Up to 8 workers, `python -X utf8`. Edit no other repo file. No commit, no push. Leave `00_基线快照_2026-10-14/` alone.
