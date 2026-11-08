# Task: make the preprocessing correct, then compare the methods on a few representative cases

Simulation only. Do not read project_records/HANDOFF.md. Read `project_records/agent_runs/2026-11-08_codex_h_fix/output.md`, `experiments/rossi_failure_study/h_fix_eval.py`, `rib_transition_eval.py`.

## 1. Correct the code

- In `contact_candidates.py`, make the increment form the default for the HC branch. Fix the MEMORY branch the same way: today's coefficients act only on today's increment, so a coefficient change cannot move earlier rows. Keep a way to get the old behaviour (a flag defaulting to the corrected form) only so old results can be reproduced.
- `ridge_estimator.py` line ~62 has the same pattern (`y_ref+((y-y_ref)+state)/mag`). Search the estimator and preprocessing code for any other place where a time-varying coefficient multiplies an accumulated quantity in cached rows. Fix those on the live paths of the methods below, and list every place, fixed or not.
- Checks: at constant coefficients old and new increments agree to round-off; a coefficient step at rest moves the output by 0; with the old flag, one stored baseline episode reproduces exactly. If any comparison method's own code bypasses the fix, say so.

## 2. Compare, corrected code only

Fresh, unused base 874500; sensor seeds from 910000 + index. Probe mass present.

Cases: two dynamics/phase settings (`neural_white` 28 deg, `viscous` 14 deg), three cells each:
- scene A `breath_rib`, i = 0, 2, 6 (nominal)
- scene A with `scene_for(..., 'slow_thin')`, i = 1, 4, 5
- scene B `breath_fast`, i = 0, 1, 2
- `breath_rib_fast`, i = 0, 2, 6
- `cough` and `combo` via `B.cells()`, i = 0, 3, 7

That is 36 episodes per method.

Methods, parameters as they stand, no tuning:
- Identifiers in the Method with gates: Huber-D4, D4, Rossi-ID, Kalman, MHE.
- Huber-D4 gates off (as in `memory_segment_eval.py`), no memory.
- Huber-D4 gates off plus memory, 6-tick updates as in `rib_transition_eval.py`: Kalman tau=2, Kalman tau=0.2, IMM, Lepski, EWA, forgetting with floor, Huang ARD fixed precision.
- Complete controllers: Rossi-control, DO-MPAC, Xing-DAHIC.
- Oracle: window estimate = 1.2*truth through the gates-off host.

## Report

In Chinese, final message only, as short as the content allows:
- the code changes and checks;
- one table per scene, method by: L mean/median, force RMSE median, max peak, episodes above 15 N;
- one table of the transition metrics of `rib_transition_eval.py` over large transitions in the two scene A sets and A+B: lag and settling time, largest under-estimate, force peak, rising and falling separately, with the count of transitions (say when too few);
- plateau bias and sd per method.
Report what was measured. No winner selection, no method suggestions.

Script `experiments/rossi_failure_study/corrected_comparison.py`, results `results/corrected_comparison.json`. Up to 8 workers, `python -X utf8`. No commit, no push.
