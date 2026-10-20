> Existing local validation. Obsolete topics/plans were removed; original configurations/results remain. See[HANDOFF](../../../project_records/HANDOFF.md).

# Continuous slow-variation contact loop: step 2

2026-10-13. **Small-model theory/contact simulation support directNEOC in slowly varying environments; large-jump failures do not invalidate the idea.** A smaller-neighborhood, denser-query nonlinear test preserves accuracy/reduces total computation; simple fixed control is also sufficient. Results concern these carriers, not full R15 or paper-level contribution.

## 1. Target mechanism

Physical environments vary continuously; methods receive identical accurate **current** parameters at common periods, without future parameters. OrdinaryNEOC uses no exact-root queries, candidate residual acceptance or rollout checks. Maintain anchors by preset neighborhoods, or sensitivities along predicted points; count the latter's real cost.

First use finite-horizon single-value nonlinear normal cylindricalHertz projection. State: indentation, velocity, actual stiffness; action: stiffness rate, integrated with energy records. Damping8, mass0.2, target≈0.1414N; E(t)30→36kPa with a smooth cyclic test. Predict the same1.5s frozen-current problem while physical parameters vary between samples. No full MHE, preview, tangential/imaging goals or compressedNN critic.

Compare local-reuseNEOC, maintainedNEOC, cached-J chord, oneNewton, full structured solves and fixed stiffness. Full solving is one backwardBellman sweep; warm starts cannot reduce sweep count. Do not substitute inefficient genericNewton as strongest full solve. Share problem, grid, terminal cost/actions.

## 2. Numerical-carrier diagnosis/repair

Euler-grid refinement did not remove control-cost bias; retain[initial record](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98NEOC%E9%AA%8C%E8%AF%81_2026-10-13/%E9%98%B6%E6%AE%B52_%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98%E6%8E%A5%E8%A7%A6%E9%97%AD%E7%8E%AF/00_%E5%88%9D%E5%A7%8BEuler%E8%BD%BD%E4%BD%93%E8%AE%B0%E5%BD%95/pilot_results.json). Switching contact prediction/stage costs toRK4 gives0.333%,0.294%,0.168% gaps across three grids versus independent nonlinear trajectory optimization; see[resolution audit](resolution_audit.json). This reference is feasible/local, not global; directional gradient error≈2.4e−8.

Reuse old minimizing actions inF_η to avoid repeated minimization. At two anchors, relative derivative differences from original differences≈4.4e−10/6.0e−10, and from full-root differences≈0.0122%/0.000705%. Preparation falls≈3.3–3.8s→0.5s; see[derivative audit](01_reused_action_sensitivity/derivative_audit.json). This optimizes implementation, not new control modules.

Piecewise node interpolation does not establish original whole-neighborhoodC² roots. Parameter radii below are diagnostic presets, not strictR15 deployment certificates.

## 3. Four preset tests and one exploratory test

Original radius0.025, period0.25/0.5s. All pass retrospective sampled contact, gain/rate/energy checks. Cost gaps use full solves with equal information. Original joint diagnostic thresholds remain: cost increase≤1%, force RMSE increase≤5%, total computation reduction≥30%.

|Path, duration/parameterperiod|Localreusecostgap|Maintainedcostgap|Localanchorsincludinginitial|Jointcheck|
|---|---:|---:|---:|---|
|ramp,3s/0.25s|1.612%|0.000%|6|Fail|
|ramp,6s/0.25s|1.559%|-0.158%|6|Fail|
|ramp,6s/0.5s|1.588%|0.000%|6|Fail|
|cycle,6s/0.25s|0.270%|0.387%|13|Fail|

Maintenance matches control quality but global update costs do not beat one structured solve. Local monotone reuse costs≈1.6% more; small parameter steps alone do not validate neighborhoods, which also involve finite approximations/value-gradient sensitivity.

**Open a separate exploratory test:**same material path, physics, costs and thresholds, halve radius to0.0125 and set0.1s updates. Retain original failures and criteria; see[configuration/results](03_small_neighborhood_dense_updates_exploratory/results.json).

|Method|Accumulatedcost|ForceRMSE/N|Totalcontrolcomputation/s|
|---|---:|---:|---:|
| local_NEOC | 0.04280730 | 0.00440104 | 30.9501 |
| one_Newton | 0.04259169 | 0.00437526 | 87.7647 |
| full_solve | 0.04263196 | 0.00437943 | 72.8112 |
| fixed_gain | 0.04259169 | 0.00437526 | 0.0015 |

NEOC increases cost**0.411%**/force RMSE**0.494%**.60 queries use initial+14 anchors and45 predictions. Total solving, sensitivity, maintenance, prediction/action recovery falls**72.81s→30.95s**,**57.5%** reduction. RecheckNewton/full/fixed control, not chord; this is not all-optimizer superiority.

Quality simulation spans6s virtual time without injecting solve latency.14NEOC maintenance events exceed0.1s; full solving has59 overruns. Total cost exceeds task duration. This supports relative amortization,**not real-time robots**. Timings are localPython, not compressedNN critics.

## 4. Smooth-root test: prove range before directNEOC

Independent local linear contact uses exact three-coefficientRiccati critics for normal PD. Actual gains interpolate causally; effective stiffness varies300→450N/m, and all receive identical current surface motion. Different physics/costs from gain-rate optimization cannot form completeR15 certification.

[Analytical conditions](analytic_conditions_for_smooth_comparator.md)give roots, sensitivities/gain-curvature bounds. Error budgets≤0.1N/m,≤0.002Ns/m yield radius**13.027N/m**. Stability includes variation, sampling/interpolation; all four decay margins are positive. Predict/release within range without candidate-dependent acceptance. Floating-point checks are not interval certification.

|Path, duration/parameterperiod|LocalNEOCcostgap|Initialfixedimpedance|Midpointfixedimpedance|
|---|---:|---:|---:|
| ramp，3s / 0.25s | -0.0100% | 3.70% | -2.27% |
| ramp，6s / 0.25s | -0.0073% | 3.36% | -1.58% |
| ramp，6s / 0.5s | -0.0189% | 6.73% | -4.35% |
| cycle，6s / 0.25s | -0.0018% | 5.06% | -1.17% |

The reference is the full analytic control solution with identical information and frozen parameters. Local NEOC differences are below 0.03% in all four groups; maximum stiffness error is about 0.0902 N/m, within the preset bound, and all sampled constraints pass. Along-path maintenance also stays close to the reference, but cyclic paths accumulate more error; see [complete data](02_smooth_quadratic_critic_comparator/results.json).

Midpointfixed control beats the frozen reference, which is not globally optimal for moving/time-varying tasks. Negative gaps are not beyond-optimality, and weak initial-fixed control does not prove adaptation necessary. Cheap closed forms prevent universal speed claims from genericCARE comparisons.

## 5. Local conclusions and limits

**Mechanism supports further work:**bounded neighborhoods/proper sensitivity maintenance permit directNEOC under slow changes; few maintenance events/dense queries can amortize computation. No mandatoryNewton after every prediction or online candidate reacceptance is needed.

Fixed stiffness matches optimization in nonlinear simplified tasks; midpointfixed is better in smooth cases. These bound carrier benefits, not all B0. Preview/noisy estimation are untested; this report does not prescribe later research.

Reproduction:16 nonlinear baseline trajectories match across original/optimized code;32 smooth trajectories match across reruns; denseNEOC rerun matches. Four original/optimizedNEOC cost gaps<1e−6 relative. Timing varies; see[record](reproducibility.json). Not all nonlinearNEOC trajectories were independently rerun.

Entries:[carrier](continuous_contact.py),[fourtests](run_contact_validation.py),[structure](run_structured_check.py),[denseexploration](run_dense_exploration.py),[smoothcontrol](smooth_contact_check.py),[reproduction](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98NEOC%E9%AA%8C%E8%AF%81_2026-10-13/%E9%98%B6%E6%AE%B52_%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98%E6%8E%A5%E8%A7%A6%E9%97%AD%E7%8E%AF/verify_stage2.py).Retain original results, failures and three core-theory originals at this checkpoint.
