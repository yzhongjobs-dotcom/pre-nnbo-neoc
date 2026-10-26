# Experiment index

## Current

[Cobot scan benchmark](cobot_scan_benchmark/README.md): synthetic PC125 normal-axis plant (8 ms outer loop, 8 ms transport delay, 19 ms servo, 1 kHz force with 30 Hz low-pass and 0.08 N noise). It holds the current C1 and C2 implementations, repaired field baselines and the seed-7200 three-task evidence.

| Entry | Use | Boundary |
|---|---|---|
| [C1 repair results](cobot_scan_benchmark/contact_update_results/README.md), [controller](cobot_scan_benchmark/contact_update.py), [evaluation](cobot_scan_benchmark/evaluate_contact_update.py) | Current C1 on the original three tasks, estimator ablation, peaks, recovery, timing | T2 development data used for the repair; baselines keep earlier settings |
| [Completed tasks](cobot_scan_benchmark/completed_tasks/README.md) | Seed-7200 protocol, 40 paired scenes per task, fixed extremes, figures | T1/T2 transfer frozen T3 settings; not separately retuned rankings |
| [Repaired baselines](cobot_scan_benchmark/refined_baselines/results.md) | Xue, Jin, learned VIC and others after adaptation repair | Task adaptations, not reproductions of the original robot experiments |
| [Paper tables](cobot_scan_benchmark/paper_tables/), [builder](cobot_scan_benchmark/build_paper_experiments.py) | Original-paper task/metric mapping used by the experiment section | No new control simulations |

The Chinese simulation section is [paper_experiments_final_zh.md](../research/paper_experiments_final_zh.md).

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
