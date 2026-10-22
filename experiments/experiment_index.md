# Experimental resources for subsequent research

Retained materials serve environment estimation, both candidate designs, reasonable simple/strong comparators, and physical/numerical checks. They use different carriers, protocols and implementations; they are not a unified ranking of complete candidates 1/2 and need not be rerun in table order.

## Purpose and reproduction entries

| Resource | Subsequent use | Required boundary |
|---|---|---|
| [Cobot benchmark and repaired baseline comparison](cobot_scan_benchmark/README.md) | October 21 task redesign; same-scene physical metrics after Xue/Jin/learned adaptation repair and additional train/dev tuning | Synthetic PC125 normal-axis plant; frozen C1 v2/C2; historically inspected test cohorts; accuracy gains do not establish peak-force dominance or full B0/journal readiness |
| [Realistic scan task](realistic_scan_task/README.md) | Planar 3R arm, viscoelastic B0 footprint contact, F/T noise/delay, three representative tasks, tolerance-based score, headroom check against fixed impedance and force PI | First benign tasks only; no candidate ported yet; synthetic parameters |
| [NNBO source reproduction](nnbo_source_reproduction/README.md) | Critic equations, exploration and Newton position inversion | Minimal mechanism, without complete identification/robot experiments; not B0 four gains |
| [Finite-horizon contact](finite_horizon_contact/README.md) | Cylindrical contact, Bellman carrier and independent trajectories; dependency of slow-variation tests | Single K rate/fixed D, not the full-dimensional B0 critic |
| [Continuous slow variation](continuous_slow_variation_updates/README.md) | Sensitivity/re-anchoring, cost accounting and slow/small-step instances | Retain both original failures and exploratory groups, not favorable selection |
| [Dynamic strong comparators](dynamic_contact_strong_comparators/README.md) and [Q2](dynamic_contact_strong_comparators/Q2/README.md) | Actual K/D, energy execution, feedback/cheap prediction and independent integration | Strong conceptual controls, not whole-paper reproduction; Q2 assumes a known motion program |
| [Observation/local roots](observation_local_roots/test_report.md) | Separate estimation, geometry and critic approximation errors | Known geometry/ideal normal moment, not complete MHE |
| [Parameter/observation boundaries](parameter_observation_boundaries/README.md) | Descriptor revision, physical slow variation and full-wrench conditions | Rechecks B0 boundaries, not a new field failure |
| [Constraint switching/value roots](constraint_switching_value_roots/candidate_description.md) | Regularity, reduced-root/projection tools candidate 1 may borrow | Conditional revision, not a third selected question; multinode actions may equal simple prediction |
| [B0 evidence audit](b0_problem_evidence_audit/README.md) | Leading-footprint information, geometry timing and model scales | Known structure/synthetic parameters, not useful-preview, real-load or MHE validation |
| [Contact-history audit](contact_history_mechanism_audit/README.md) | Realizable histories, ordinary observers and prediction versus control benefits | Two known Maxwell histories, not unknown-tissue/scan-control validation |
| [Optimal impedance updates](optimal_impedance_update/README.md) | Network/RBF/simple solvers, adjoints, finite-rate optimization and delay accounting | Known local model, not complete candidate 2 |
| [Constrained contact updates](contact_constrained_impedance_update/README.md) | Four-gain execution, task-cost learning, prediction release and background improvement | 0.3s known-model transition, not spatial MHE/original PPO |
| [Within-scan prototype](within_scan_impedance_feedback/README.md) | Local responses, fixed/analytic/reference baselines, CEM feedback and physics corrections | Small offline formula, not NN/PPO; original/corrected physics separate, not complete candidates |
| [Spatial estimation audit](spatial_contact_estimation_audit/README.md) | Extracted material/displacement estimation, noise propagation and scale tools | Synchronous sampling/affine conditions; no closed-loop/novelty evidence |

Added [MHE–impedance engineering precheck, paused at that time](mhe_guided_impedance_pilot/README.md). Reduced rolling estimation, gain execution and local-gradient tools are reusable. It stopped after the user's baseline correction, without strong NNBO-type methods, completed NN training or independent tests; it cannot establish candidate superiority or journal readiness.

Work then resumed under the user's latest instruction: [separate-candidate iteration](separate_candidate_validation/README.md), [grouped independent results](separate_candidate_validation/independent_results.md). Includes actual rate actions, state coverage/sampling alignment, separate model/learning experiments and paper-mechanism adaptations. Nonzero-force NNBO inversion failures are isolated from contribution comparisons. This carrier omits full B0 objectives/constraints and does not establish complete candidates/journal readiness. The old paused directory is a historical record, not the present stopping point.

**Latest: [dual-candidate task simulations](dual_candidate_task_simulation/README.md), [results](dual_candidate_task_simulation/results.md).** Implemented both routes from the research question, then selected Jin multistep predictive VIC task adaptation. Restored B0 two-region Hertz, positive-definite gain/rate costs, terminal terms, four actual gains and trajectory checks. Candidate 1 uses local quadratic HJB feedback/directional Riccati updates, not original B0 collocation roots. Candidate 2 is an actual task-gradient NN with three training seeds. Keep 24 independent main cases, preview/neighbor ablations, eight softening/faster-surface cases and isolated compute/release experiments separate. Both candidates still trail comparator force/tangential performance; low totals due to gain preferences do not establish a win.

## How to use key results

- Local NEOC/sensitivity computation does not stably beat adequate ordinary solvers; even one Newton step may be better in finite-horizon cases. This limits specific speed claims, not the original complete question.
- Steady-state networks/RBF approximate numerical references; RBF is faster, and tuned fixed gains approach references in some tasks. Weights, training configuration and results support reproducibility of this comparison, not original candidate 2.
- Release-aligned background model improvement has local gains about 0.67%–0.78%; old waiting negatives alone cannot reject the model route. Difficult contact/energy constraints were not active.
- Initial 48-case joint CEM feedback improves about 0.433% over analytic control; 24 new pure-K/D cases do not beat continuously tuned fixed control; a separate corrected-physics 24-case joint comparison gives about 0.201%, interval crossing zero. Protocols differ; no pooled sample or attribution of joint feedforward benefit to gain learning.

## Materials removed from the workspace

Removed misread candidate-2 direct PPO training, candidate-A scan/coverage/public-acquisition branches, unselected policy-library rollout and two-stage posterior planning pipelines. No-advantage outcomes/lessons appear in the [failure review](../project_records/research_failures_and_corrections.md); original versions remain in Git. Independent spatial-estimation tools were extracted and retained. Removal does not revoke measurements or invalidate original candidate 2.

## Usage notes

Source, parameters/protocols, necessary models, per-case results and independent checks correspond to their historical versions. Old result source hashes may refer to the commits stated in reports; never silently relabel them as current hashes. Select the relevant Git version and reproduce in a copy to avoid overwriting results. Historical paths retain their spelling in the corresponding Git version.

Root requirements cover older numerical resources; newer environments are specified by local requirements/environment.json. Preservation checks are not scientific reruns or validation of synthetic physics on real tissue. See [HANDOFF](../project_records/HANDOFF.md) for continuation.
