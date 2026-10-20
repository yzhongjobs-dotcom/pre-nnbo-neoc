> Existing local validation. Obsolete topics/plans removed; original configurations/results remain. See [HANDOFF](../../../project_records/HANDOFF.md).

# Q2 continuous phases/causal-model updates: results and decisions

2026-10-13.Stage validation executed; complete objective unfinished.
The main revision concerned question interpretation; see [question/method choices](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E9%97%AE%E9%A2%98%E6%9C%BA%E5%88%B6%E9%87%8D%E6%96%B0%E5%AF%B9%E9%BD%90_2026-10-13/00_%E7%A0%94%E7%A9%B6%E9%97%AE%E9%A2%98%E4%B8%8E%E6%96%B9%E6%B3%95%E9%80%89%E6%8B%A9_%E8%BF%AD%E4%BB%A3%E4%B8%AD.md).

## Material Passport

- Skill: academic-research-suite / experiment-agent，run + independent validation
- Type: continuous local contact / same-information rival discrimination
- Status: executed, with independent equations/integration, gradients, storage/static-target checks
- Scope: local soft contact, two actual gains/rates and energy, not full R15/real tissue
- Numerical references are converged local trajectories, not global optima.

## Scenarios and same-problem comparisons

4s near-preload contact, material300→550N/m, two smooth surface motions;
compute actual states/force from dynamics. See [protocol](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E9%97%AE%E9%A2%98%E6%9C%BA%E5%88%B6%E9%87%8D%E6%96%B0%E5%AF%B9%E9%BD%90_2026-10-13/Q2_%E8%BF%9E%E7%BB%AD%E9%98%B6%E6%AE%B5/%E7%A0%94%E7%A9%B6%E5%A5%91%E7%BA%A6.md) and configurations.
Phases are not preset stiffness paths. Fix position, force, compliance preferences, action costs/bounds together.

|Comparator|Accumulated cost|Gap from repaired reference|
|---|---:|---:|
|Tuned fixed targets|1.464683|+25.63%|
|Equal-information anticipatory schedule|1.283045|+10.05%|
|Strict accuracy/whole-stage-feasible one-step|1.200463|+2.96%|
|Low-level analytical constraints/one-step|1.191440|+2.19%|
|Numerically guarded trajectory|1.165898|—|

Tune fixed targets by multistart Powell; search schedule preview/damping thoroughly,
not merely three gains/few dampings. Static targets solve common physical-cost stationarity,
with maximum K difference about 9e-6 from independent scalar optimization.

Analytical low-level projection replaces expensive stage checks as baseline repair, not our new mechanism.
It enforces state-dependent rate/energy actions each low-level cycle; full prediction remains timed.
Its actions exceed the flexibility of125ms fixed-rate reference segments, so negative gaps would
first indicate restricted reference parameterization, not beyond-optimality.

## Timing limits

Repeated whole-stage checking costs about 11.02s total, max decision about 1.77s;
analytical interfaces reduce this to about 1.78s/max0.339s.
Full prototype trajectory optimization about 28–32s; separate preparation/JIT, control service/physics.

These do not prove necessary NEOC speedups: oracle future materials permit baseline pre-solving;
prototype timings differ from structured RTI/hardware, and latency is not injected into loops.

## Continuous-time numerical repairs

[Independent audit](audit.json) finds sampled-feasible coarse roots violate continuous power limits by about 1.18e-8W.
Retain old data and add common normalized planning guard3e-5.
[Repaired results](guarded_validation/results.json) converge from two starts to
about 1.165898, with independent minimum power margin about 9.17e-9W and energy above original bounds.

Independent DOP853 relative cost difference about 6.65e-8, storage identity≤1.39e-17;
three directional objective-gradient errors≤1.42e-9, constraints about 2.05e-8.
Analytical-action projection independent relative cost difference about 1.65e-6.
These are numerical/sample checks, not reliable full-horizon interval certification.

## Removing future material information

Controllers receive current/past material observations; all still know future surface programs as task inputs.
Material continues changing; low-level controls see only current actual states, not future material.
300→1000 is a holdout range; tail policies tuned on300→550 remain frozen.

|Material range|Frozen initial model|Current updates/frozen future|History-slope prediction|
|---|---:|---:|---:|
| 300→550 | 1.167293 | 1.170405 | 1.191166 |
| 300→1000 | 1.589625 | 1.620307 | 1.584444 |

Updates yield no stable significant gains; wide-range trend prediction improves frozen models only about 0.326%.
This does not prove that accurate models lack value: the comparator uses approximate tail policies, local action optimization and approximate predictions,
model, representation and solution errors can compensate. It shows these results cannot
establish model updating as a necessary task capability.

Strengthen cheap controls with maxiter=8, ftol=1e-5 and prior-action warm starts;
some retain quality with maximum 77–146ms decisions, still prototype observations.
Do not presume all strong simple solvers miss real time, or call virtual-time quality deployment.
See [causal controls](material_update_discriminator/results.json) and
[trend repair](trend_discriminator/results.json).

## Structural constraints on NEOC

For fixed nodes, weights/references, material enters only contact dynamics/cost;
the rate Hamiltonian allowed set/action recovery does not directly depend on theta.
Thus finite value residual F(Wa, theta)=A0(Wa)+A1(Wa)theta+A2theta² exactly.
Cached-J chord can store two response columns and evaluate linear/quadratic terms online;
include it as a strong comparator with equal preparation opportunities.

[Structure code](prepared_residual_structure.py) checks this at arbitrary weights,
with error about 2.2e-16; cached/explicit chord differs about 1.7e-13.
These weights are not qualified HJB roots; no critic-quality/control-acceleration test is implied.
This structure is not standalone novelty or automatically valid for all η, changing references/state boundaries.

## Scope

Local continuous optimization improves fixed control; strong simple methods approximate it; causal-model gains are not stable/significant. Limited model, cost/disturbance evidence neither invalidates B0 nor creates an algorithmic subquestion.

Position changes in known fixed spatial fields are state changes; model-estimate revisions alter descriptors.

## Running and recovery

- [Physics](contact_stages.py),[strong schedules/prediction](strong_rivals.py)
- [Analytical interface](projected_rival.py),[independent equations](independent_audit.py)
- [Guard repair](guarded_validation.py),[causal ablation](material_update_discriminator.py)
- [Trend prediction](trend_discriminator.py),[prepared structure](prepared_residual_structure.py)

strong_rivals/, projected_rival/ retain preguard results; guarded_validation/ retains repaired results. First-pilot intermediates left the working tree but remain in Git.

Use new outputs/copies; order: physics→strong controls→independent checks→
guard repair→information ablation→trend prediction. Directories retain parameters/raw states/actions.

## Checked precedents and boundaries

- [Beber ICRA2024](https://arxiv.org/html/2309.14893v2): tissue characterization, impedance QP/energy mechanisms already exist.
- [Rai NEOC](https://arxiv.org/html/2312.04733v1): closed-loop parameter sensitivity/multistep treatment exist.
- [ENE](https://arxiv.org/html/2306.04830v1): state/preview and piecewise active changes exist.
- [Parameter prediction/correction, SIAM2012](https://doi.org/10.1137/110844349): parameter tracking/contraction analysis exist.
- [Multistage MPC sensitivity, LCSS2018 author version](https://jaschke.folk.ntnu.no/preprints/2018/Krishnamoorthy_LCSS_Path-following/PrimalDecompQP_L-CSS.pdf): neighboring parameters already avoid repeated full solves.

This is a bounded mechanism check, not final novelty review.
