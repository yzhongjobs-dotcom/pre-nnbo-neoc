# Experiment index

## Current (2026-11-09)

| Directory | Use | Boundary |
|---|---|---|
| [rossi_failure_study](rossi_failure_study/) | Current physical plant (`physics_plant3.py`, `scenes3.py`), Method entry (`method3.py`), Huber-D⁴ estimator, comparison code, real-data derivations | Current results: [three comparison tables](rossi_failure_study/comparison_tables.md) from `identification_completion.py` and `corrected_comparison.py`. The [baseline](rossi_failure_study/baseline_interpretation.md) and [overlay](rossi_failure_study/overlay_interpretation.md) interpretations predate the H correction and are invalid as conclusions. TVK and jump GLR results also predate it and have not been redone |
| [motion_aware_scan](motion_aware_scan/README.md) | IMM, admittance and feed-forward code the current Method uses; `Proposed-IMMb` evidence for the paper draft; feed-forward loop checks | The draft describes the 2026-10-31 method, not the current default |
| [problem_mechanism_check](problem_mechanism_check/README.md) | Linear failure demonstration for the paper draft | Not a benchmark |
| `cobot_scan_benchmark/` | Only `plant.py` and `scenarios.py`, which current code imports | Older controllers removed |
| `method_redesign/` | Only `v1_servo_id.py` (online servo identification), which current code imports | Other prototypes removed |

## Stale intermediates removed from the working tree (2026-11-09)

Removed from `rossi_failure_study`: pre-H-correction baseline and overlay logs, replays, module results and summaries; `memory_closed_loop_*`, `memory_root_causes`, `memory_segment_eval`, `rib_transition_eval`, `h_fix_eval`, `low_estimate_forensics` and `id_layer_root_cause` results; and the scripts `h_fix_eval.py`, `memory_root_causes.py`, `overlay_report.py`, `overlay_plots.py`, `baseline_report_r2.py`, `baseline_diagnostics.py`. Recover from Git `e5b2de9`.

Kept because current code uses them: `baseline_test_r2.json` (frozen configurations), `baseline_dev_r2.json`, `baseline_test.json`, `baseline_calibration.json`, `baseline_report.py`. The old-form H switches are checked against `results/legacy_h_fixture.*`, one saved episode. The `main` stages of `memory_segment_eval.py`, `rib_transition_eval.py` and `low_estimate_forensics.py` read removed pre-correction references; those modules are used as libraries now.

## Discarded methods removed from the working tree (2026-11-08)

Removed: `contact_information`, `contact_model`, `conditional_linear_contact`, `evidence_memory_contact`, the rest of `method_redesign` (estimators E0–E4, constrained MPC, predictors, V3/V4) and the rest of `cobot_scan_benchmark` (C1/C2 controllers, repaired baselines, seed-7200 evidence, paper tables). Recover with `git show ca398be:experiments/<dir>/<file>`, or check out `ca398be` into a copy. The reasoning stays in `research/` and `project_records/handoff/03_方法演进审计.md`; their links to these paths resolve only at that commit.

## Early experiments removed from the working tree

On 2026-10-26 the early B0 and 2D directories were removed from the working tree. Recover them from Git `b7a3f91`, for example with `git show b7a3f91:experiments/<dir>/README.md`, or by checking that commit out into a copy.

The removed directories are:
- NNBO source reproduction, finite-horizon contact and continuous slow variation;
- dynamic strong comparators, observation/local roots, parameter/observation boundaries and constraint-switching value roots;
- B0 evidence audit and contact-history audit;
- optimal impedance update and constrained contact update;
- within-scan prototype, spatial-estimation audit and realistic scan task;
- MHE pilot, separate-candidate validation and the 2D dual-candidate task simulation.

Lessons retained from them:
- **NEOC speed:** local NEOC/sensitivity computation does not stably beat adequate ordinary solvers; even one Newton step may be better in finite-horizon cases. This limits specific speed claims, not the original question.
- **Networks versus RBF:** steady-state networks/RBF approximate numerical references; RBF is faster, and tuned fixed gains approach the references in some tasks.
- **Background model improvement:** release-aligned background improvement gives local gains of about 0.67%–0.78%, so old waiting negatives alone cannot reject the model route.
- **Within-scan CEM feedback:**
  - the initial 48-case joint CEM feedback improves about 0.433% over analytic control;
  - 24 pure-K/D cases do not beat continuously tuned fixed control;
  - a corrected-physics 24-case joint comparison gives about 0.201%, with an interval crossing zero.
- **2D dual-candidate simulation:** both candidates trailed the Jin MS adaptation on force and tangential error. Low totals driven by gain preferences did not establish a win.
- **Scope of these instances:** none of them was the complete B0, the original candidate 2 or a robot validation. Detailed mistakes are in the [failure review](../project_records/research_failures_and_corrections.md).

## Usage notes

Result hashes refer to the versions stated in each report. Reproduce in a copy so original data are not overwritten. File checks are not scientific reruns or physical validation.
