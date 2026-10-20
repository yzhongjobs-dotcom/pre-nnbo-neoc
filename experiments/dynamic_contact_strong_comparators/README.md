> Existing local validation. Obsolete topics/plans removed; original configurations/results remain. See [HANDOFF](../../project_records/HANDOFF.md).

# Realigning problems, strong comparators and mechanisms

Retain local contact dynamics, strong feedback/short-prediction code/results as references, without replacing B0.

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + independent numerical validation
- Version Label: physical_transient_alignment_v1
- Verification Status: executed; optimization actions/trajectories rerun identically from the same start.
  Independent integration, gradient/energy checks follow, not full R15/global certification.
- Type: local physical contact simulation / rival-mechanism discrimination
- Status: local-instance review complete; complete question–method matching remains unfinished.

## Local test question

Q1: disturbed operation needs precision; recovered operation needs compliance, while actual stiffness/damping face
rate/energy limits. Does this dynamic tradeoff require continued optimization and naturally support original NEOC?

Unlike information proxies, contact dynamics determines errors; no target stiffness trajectory is imposed.
Model: local linear soft contact near preload; states are position error, velocity, actual K/D and energy.
All methods share physical costs, force/gain/rate/energy bounds.
Robot incremental control/contact forces and actual gain integration are included, not original four-dimensional gain-rate HJB.
Ranges distinguish mechanisms without tissue/platform calibration; material holdouts freeze local parameters,
not a full continuously changing environment. They do not invalidate continuous scanning.

## Judgments changed by results

|Order|Hypothesis|Result|Decision|
|---|---|---|---|
|1|Recovery/compliance defeats fixed targets|Full rate cost 0.07374682 versus fixed0.08377740, about 11.97%better; stronger preselected-initial-gain diagnostic0.08220324|State-dependent impedance remains useful, complex optimization not necessary yet|
|2|Dynamic feedback closes fixed-control gaps|Global search/Powell six-parameter schedule is only0.07586%above reference|Tuned feedback approximates local optimization|
|3|One schedule covers materials/states|Six holdouts: worst frozen gap9.60%; independent retuning≤0.363%|Coverage differs from mechanism inability; retuning is not free online|
| 4 | Can material-conditioned lookup repair coverage? | Four material anchors with the same initial disturbance are accurate on new materials, but some new-state costs are up to 12.30% higher | Material lookup alone cannot ensure state coverage; this local failure does not refute all conditioned policies |
|5|Cheap one-step prediction repairs gaps|26-candidate rollout: five ample-energy groups0.061%–0.176%; low energy 1.42%–2.41%|Check implementations before claiming full HJB necessary|
|6|Low-energy gaps: coupling or grids?|Refine current two actions: seven gaps0.00240%–0.06478%, sample constraints pass|Action resolution dominates, not intrinsic class differences|

Original simple-schedule L-BFGS-B line-search failure was not treated as incapability.
Independent global search/Powell found feasible successful solutions; retain originals.
Fixed targets start from common actual gains; preselected fixed-initial gains are separate diagnostics.

## What the strong simple predictor does

For each current gain-rate action, predict one step in the same contact model, then use cheap fixed
feedback to the endpoint. Choose feasible low-cost current actions; execute
only that action. Refined versions optimize current two actions continuously and retain cheap future feedback,
without querying optimal trajectories, true values or NEOC answers.

This is established rollout/one-step policy improvement, not a new contribution.
It is not whole Beber/Jin/Xue reproduction; it demonstrates a simpler same-problem mechanism
repairs this failure, not exhaustive comparators or theoretical dominance.

Seven-group maximum observed decisions12.92–21.22ms versus50ms periods.
These are local single-run maxima, not WCET, platform real time or end-to-end deployment.
Record JIT separately; exclude physical simulation/independent references from controller time.
Cheap tail policies need offline tuning, not free preparation, but preparation was not fully timed,
so report online timings only, not total/amortized advantage.

## Verification evidence

- Four optimization starts reach similar feasible costs; same-start actions/trajectories
  match bitwise; see [reproduction](pilot_compiled/reproducibility.json).
- Compiled/original Python RK4 real trajectories, constraints/costs match; derivative differences are floating-point scale.
  Compilation reduces overhead, not a methodological contribution.
- Five directions: maximum cost-gradient error4.4e-11, constraint error1.36e-9.
-50-state robot-storage identity error≤1.04e-17.
- Main DOP853 relative cost difference0.00251%; halved-step 0.00234%.
- Refined rollout on one ample/two low-energy groups: independent relative cost differences≤
  0.00056%; quarter steps retain low-energy bounds.
  These are numerical/sample checks, not continuous interval certificates.
- Cost uses RK4 stage integration; RMSE uses common sample grids,
  so differing sample-density RMSE is not interpreted as changed physics.

See [main audit](pilot_compiled/audit.json),[predictor audit](rollout_rival_refined/audit.json).
References are converged local nonlinear optimizations, not true global optima; tiny gaps are not
strict absolute-optimality certificates.

## Scope

Q1 shows fixed-target loss and near-reference tuned feedback/refined prediction in this local task. It limits mechanism benefits, not B0 validity or future tasks.

See [Q2](Q2/README.md) for continuous phases/causal-information ablations.

## Code and data

- [Physics](dynamic_tradeoff.py),[equivalent compiled flow](compiled_flow.py)
- [Independent checks/holdouts](audit_and_family.py),[energy boundary](energy_boundary.py)
- [Conditioned control](conditioned_rival.py),[one-step control](rollout_rival.py)
- pilot_compiled/: full main case, audits/reproduction.
- material_and_state_holdout/: six frozen/retuned groups.
- energy_boundary/: four actual energy budgets; tightest main case uses about 98%initial margin.
- conditioned_rival/: eight lookup/state-coverage groups.
- rollout_rival/, rollout_rival_refined/: retain both coarse/refined versions.

Run order in this directory:

```powershell
python -X utf8 dynamic_tradeoff.py --name pilot_compiled
python -X utf8 audit_and_family.py audit
python -X utf8 audit_and_family.py family
python -X utf8 energy_boundary.py
python -X utf8 conditioned_rival.py
python -X utf8 rollout_rival.py
python -X utf8 rollout_rival.py --refine
```

Use copies/new output directories rather than overwriting results for reproduction.
