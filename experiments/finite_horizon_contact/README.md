> Existing local validation. Obsolete topics/plans were removed; original configurations/results remain. Current plans:[HANDOFF](../../project_records/HANDOFF.md).

# R15 rapid-validation repairs and actual closed-loop results

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: carrier repair + numerical validation + reproducibility
- Origin Date:2026-10-13(work began2026-10-12)
- Verification Status: REPRODUCED; two results, candidate hashes and9 original CSVs match exactly; timing retained separately
- Status: LOCAL_FINITE_HORIZON_CONTACT_COMPARISON
- Main text, appendix, eight proofs and existing failures were unchanged at this experiment checkpoint.
- Main runs≈188.1/188.9s, both exit normally; scope follows.

## 1. Local results

The user authorized autonomous fixes to implementations, settings/carriers, discussing core changes only where evidence supports them. Work **continued beyond initial root failure**, repairing finite-horizon carriers, checking policy quality independently and executingNEOC contact loops.

Three findings:

1. At50% modulus increase, directNEOC preserves constraints but costs≈29% more than reliable full solving; oneNewton gap≈0.8%.
2. At5% increase, direct prediction is accurate, but fixed stiffness matches full solving.
3. Preset stiffness150 is only≈0.15% above independent optimization. Meaningful necessity of exact optimization is unestablished, and global checks consume most generation speed gains.

These limit this one-change/one-stiffness carrier, not continuous slow variation, B0 or project investment. Spatial-preview benefits are not compared.

## 2. Repaired causes

|Old issue|Repair|Limit|
|---|---|---|
|Infinite-horizon NNBO meaning under sustained contact incomplete|Use existingR15 1.5s horizon/terminal cost|No whole-paper repair/rejection|
|Three-term critic interpolates nodes but not domain values|Fixed-grid single value/semi-LagrangianBellman|Validation carrier, not frozen main collocation implementation|
|No qualified common root|Backward triangularBellman gives unique discrete root/sensitivity|Discrete uniqueness is not exact continuous optimality|
|Uniform-grid diffusion|Concentrate normal-velocity/indentation resolution, refine time|Retain coarse failures; no method-specific weights|
|Self-evaluation by one solver|Independent constrained shooting/analytical gradients|Feasible local reference, not global proof|

Final root residual0.0; sensitivity versus perturbed roots differs≈0.228%.NEOC uses zero onlineNewton. Exact terminal parameter dependence is retained, not predicted linearly.

## 3. Physical setting and simplifications

Original30→45kPa soft-to-hard cylinder scan: R=30mm, L=20mm, initial half-width3mm, speed8mm/s, path−6→6mm, duration1.5s, target≈0.141372N. Synthetic scales, not clinical advice; see configuration/source.

Retain normal inertia, actual stiffness/rate, cylindrical spatial contact, gain/rate bounds and energy. Damping fixed8Ns/m; no tangential/orientation optimization or full four-action validation.

Root state:(indentation, normal velocity, actual stiffness). Eliminate energy's value dimension only in the inactive envelope: even ignoring recovery, energy≥0.019678J>0.001J. Check actual energy separately; no resets/instant gain jumps.

Use fixed-nodeV=Φ(η)+W, terminalW=0, handling directΦ changes and critic sensitivity together. Next-value interpolation is piecewise linear; enumerate quadratic segment minima to obtain global action minima. Actions are stiffness rates, not handwritten gain paths.

Grid boundaries truncate **value queries**, not physics. ActualRK4 paths remain in-domain; refinement/independent references check truncation/discretization without all-state proofs. Piecewise equations and one soft-to-hard sensitivity are not whole-domainC²/PublishOK certification.

## 4. Common-root policy quality

Independent optimization uses30 actual stiffness nodes, explicit gain/rate constraints and physicalRK4.Analytical directional error2.22e-08, successful solve, maximum rate127.777<500.

|Grid/time|Layout|Physicalclosed-loopcost|Gapfromdirectreference|
|---|---|---:|---:|
|17×17/100steps|Uniform|0.135257|29.194%|
|33×33/200steps|Operating-domainfocus|0.105810|1.066%|
|49×49/200steps|Operating-domainfocus|0.104691|-0.002%|

Directcost0.104693408 versusBellman≈0.104691 differs≈0.0021%. A tiny negative gap does not prove better optimality; action-time representations/integration differ slightly. Use agreement as rapid quality checking, not global certification.

## 5. Common causal history and forward estimates

Both use actual-model history under old policy. After normal-force deviation>0.5%Fd, collect12 raw normal-force/pitch-moment samples5ms apart. Moments include known tangential lever effects; fitting reads no generating truth.

At t=0.455s, estimate45000.0Pa/boundary≈0m, with≈0.295s until center arrival. Noiseless known geometry, surface/tangential parameters make this finite-history fitting, not full MHE/real tissue.

## 6. Executed updates and loops

The common anchor is homogeneous; new descriptors share fitted outputs. All update identical discrete roots, costs, actions, integrators, gain/energy checks.

|Method|Loopcost|Increaseoverfullsolve|ForceRMSE/N|Samplechecks|
|---|---:|---:|---:|---|
|NEOC|0.135109|29.06%|0.011932|Pass|
|cached_chord|0.119520|14.16%|0.008887|Pass|
|one_Newton|0.105545|0.82%|0.011280|Pass|
|full_backward_solve|0.104691|0.00%|0.010916|Pass|
|Fixedstiffness150|0.104846|0.15%|0.010825|Pass|

NEOC retains contact/constraints but29.06% cost increase fails the 5%loss budget. Estimation is nearly floating-point exact, and derivative checks support local tangents;50% exceeds useful first-order performance range. NoNewton repair is relabeledNEOC success.

All quality comparisons assume **ideal zero update latency** to isolate policy error. ActualPython times follow separately; no artificial delay weakens baselines and offline trajectories are not real-time deployment.

## 7. Costs beyond prediction multiplication

|Method|Two candidategenerationtimes/s|Two modelprep/generation/globalcheck/sampledrollouttimes/s|
|---|---:|---:|
| NEOC | 0.072 / 0.065 | 10.910 / 10.500 |
| cached_chord | 11.964 / 10.916 | 22.743 / 23.234 |
| one_Newton | 10.408 / 12.301 | 22.270 / 23.785 |
| full_backward_solve | 10.574 / 9.992 | 22.573 / 21.236 |

NEOC generation is fast, but global residual checks dominate. Backward solving reuses layer minima; do not force duplicate residual sweeps. With valid fusion, fullsolves11.931/11.277s andNEOC end-to-end reductions≈8.6%/6.9%, below the 30%diagnostic criterion.

ExtraNEOC sensitivity preparation≈27.444s is separate. One descriptor event does not justify amortization over imaginary infinite queries. Future reuse must also allow baseline caching/precomputation.

About7.78million critic coefficients do not represent the specified low-dimensional main-text instance. All globally checked updates exceed remaining scan time:**real time unvalidated**. Unimplemented offline/local certificates have no credit in current results.

Nonexpansive discreteBellman givesE≤Σ_k||F_k||∞ and greedyMDPgap≤2E.NEOC global bound≈189.622 cannot certify small budgets. This is projected-grid, not physicalODE near-optimality; sample checks are notPublishOK.

## 8. Small changes

At5% modulus increase, NEOCcost0.001685543 versusfull0.001685127, gap≈0.0247%, constraints pass. This supports local prediction, not arbitrary amplitudes.

Fixed stiffness also gives0.001685127. At50% it is only0.146%above independent reference. This task does not establish sufficient practical necessity for real-time optimal impedance.

## 9. Reuse and conclusion levels

Reuse critic/policy checks, Newton/sensitivity comparisons and complete cost accounting. Large-change single-step near-optimality is unsupported; fixed gains suffice in small cases and important optimization benefits remain unestablished.

These are instance limits. Earlier necessity-before-preview ordering is withdrawn; this report sets neither future questions nor mandatory falsification of every simple control. Retain tools/valid outcomes; user requests/HANDOFF govern.

## 10. Files

- [Driver](run_repaired_validation.py),[Bellmancarrier](bellman_carrier.py),[independentreference](independent_reference.py).
- [Configuration](config.json),[results](results.json),[first run](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/results_run1.json),[reproduction](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/reproducibility.json).
- [Shared sensor history](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/common_wrench_history.csv),[independenttrajectory](independent_shooting_trace.csv),[gain knots](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/independent_shooting_knots.csv).
- `trace_NEOC.csv`,`trace_cached_chord.csv`,`trace_one_Newton.csv`,`trace_full_backward_solve.csv`,`trace_fixed_gain.csv`retain full trajectories.
- [First log](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/run_log.txt),[reproduction log](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/run_log_reproduction.txt),[summary script](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81%E4%BF%AE%E5%A4%8D%E4%B8%8E%E9%97%AD%E7%8E%AF_2026-10-12/summarize_repaired_validation.py).

Run `python -X utf8 run_repaired_validation.py`here with NumPy/SciPy. Configuration declares scale/costs. Old failure directories/three theory originals were unchanged at this experiment checkpoint.
