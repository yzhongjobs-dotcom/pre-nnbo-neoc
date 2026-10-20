> Existing local validation. Obsolete topics/plans were removed; original configurations/results remain. See[HANDOFF](../../project_records/HANDOFF.md).

# R15 small feasibility tests

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + reproducibility check
- Origin Date: 2026-10-12
- Verification Status: VERIFIED means identical independently rerun fixed-seed numerical results, not identical timing
- Version Label: quick_probe_v1
- Type: simulation / algebraic feedback check
- Status: completed
- Python: 3.10.11；NumPy: 1.26.4；SciPy: 1.13.1
- Run `python -X utf8 run_quick_checks.py`here. This retained report describes results; old summary code is outside the working tree.
- Latest run5.65seconds, exit0.

## Purpose and scope

The user authorized a simple feasibility test using simulation or theory. These tests examine two local mechanisms without rewriting the paper or replacing the R15 cylindrical-robot main line with a new toy model.

**Results support ideal spatial-moment information and small-change direct acceptance, not completeR15 acceleration, real tissue or platformPublishOK.**

Tests are independent; estimates do not feed control roots. Contact uses ideal normal first moments; control uses appendixM.3's scalar example. No full robot, sensors, surface/bias states, energy orMHE/control loop.

## TestA: can moments distinguish forward materials earlier?

Use main(9)–(11)/appendixD.1 closed half-ellipse integrals forF, m1: R=30mm, L=20mm, a=3mm, E0*=30kPa, Ef*=45kPa, b=0.5mm.12 positions−2→−0.2mm all precede center crossing.

FitF-history andF+m1-history by weighted bounded nonlinear least squares with common history/force noise, not weak single-threshold baselines. EstimateEf, b; other geometry, current modulus, surface/bias are known.100 paired fixed-seed trials per condition.

Noise: nominalσF=0.003N,σm=2×10^-6N·m; higherσF=0.01N,σm=8×10^-6N·m. Width mismatch uses truea=3mm, fitted1.05a. Artificial assumptions, not real-sensor calibration.

|Condition|Method|BoundaryRMSE/mm|ForwardmodulusRMSE/%|
|---|---|---:|---:|
| nominal_noise | F_history | 0.2105 | 5.801 |
| nominal_noise | F+M_history | 0.0484 | 1.068 |
| higher_noise | F_history | 0.9456 | 25.804 |
| higher_noise | F+M_history | 0.1823 | 4.644 |
| width_5percent_mismatch | F_history | 1.6006 | 41.073 |
| width_5percent_mismatch | F+M_history | 0.3230 | 2.790 |

All 600 fits terminated successfully; this does not mean all estimates were accurate. In the calibrated-noise scenario, the F+m1 boundary RMSE fell from about 0.210 mm to 0.048 mm. A width error of only 5% still raised F+m1 boundary error to about 0.323 mm, mainly systematic bias.

**Moment information helps this ideal model, but model/geometry mismatch may dominate random noise.** Wrist moments also include traction/transform effects; reliable m1 separation with unknown friction and joint identifiability are untested. Matched simulation/fitting checks only internal approximation consistency.

## TestB: must predicted candidates receiveNewton correction?

Use appendixM.3: z_dot=u, u≥−0.7, T=1, runningη²z²+u², criticV=(1-t)Wz², node(t, z)=(0,1), finite root

```text
F(Z,η)=(-W+η²+u²+2Wu, 2u+2W-λ, λ(u+0.7))=0
Z=(W,u,λ)
```

Anchorηa=0.8, inactiveW*=(-1+sqrt(1+4η²))/2, u*=−W*,λ*=0.ANEOC prepares sensitivities; chord cachesLU; oneNewton uses new residuals.

HereF_Z at the old root does not directly depend on newη, so chord/Newton coincide, not universally. Add eliminated`structured_chord`and`structured_ANEOC`to avoid small-matrix call-overhead advantages.

### Common analytical screen

Accept only inactive parameters,0≤W≤0.7, u=−W,λ=0.Recoveredu(t)=−(1-t)Wz(t), z(t)=exp[-W(t-t²/2)]fromz0=1 always respects the lower bound.

ScalarB=η²−W−W² gives

```text
|W-W*| = |B|/(1+W+W*) ≤ |B|/(1+W)
|J(W)-J(W*)| ≤ L_J |B|/(1+W)
L_J=2U/3+2η²/3+2U²/15，U=0.7
```

The final cost bound explicitly integrates derivatives inW. All share0.001budget; exact roots are retrospective only. Actual root/cost errors remain below bounds.

**This screens scalar finite-root-feedback errors/action feasibility, not full continuousHJB near-optimality.**

### Small and larger changes

Each radius uses identical2,000uniform descriptor queries.

|Anchor-relative range|Method|Directpasses|Passrate|
|---|---|---:|---:|
| ±0.02 | ANEOC | 2000/2000 | 100.00% |
| ±0.02 | cached_chord | 2000/2000 | 100.00% |
| ±0.02 | one_Newton | 2000/2000 | 100.00% |
| ±0.08 | ANEOC | 1842/2000 | 92.10% |
| ±0.08 | cached_chord | 1137/2000 | 56.85% |
| ±0.08 | one_Newton | 1137/2000 | 56.85% |
| ±0.20 | ANEOC | 687/2000 | 34.35% |
| ±0.20 | cached_chord | 431/2000 | 21.55% |
| ±0.20 | one_Newton | 431/2000 | 21.55% |

AtΔη=0.02, rooterror≈8.27×10^-5, costgap≈4.47×10^-6, bound≈7.54×10^-5: passes0.001without correction. At0.08, bound≈0.001247 rejects despite small actual error. Conservative rejection does not prove poor actual policies.

Δη=0.4 crosses active branches and is rejected; old guarantees cannot continue. Recovery is not implemented/guaranteed.

### Equal-screen timing

Generate/screen3,000identical small queries in7interleaved batches; report median batchmeanμs/query for two runs.

|Implementation|Firstrun|Reproduction|
|---|---:|---:|
| ANEOC | 6.886 | 7.493 |
| cached_chord | 19.335 | 20.912 |
| one_Newton | 22.082 | 23.550 |
| structured_ANEOC | 6.704 | 6.843 |
| structured_chord | 7.240 | 7.558 |

Generic prediction is faster; structured chord removes most differences. Remaining microseconds reflectPython/calls/scheduling, without significance/WCET analysis, not established strong-baseline acceleration. `results.json`separately records incremental sensitivity preparation; common roots, full certification, state recovery, MHE, communication/execution are excluded.

### More important than prediction error

Atη=0.82, continuous scalar optimum is independentlyηtanh(η), with feasible optimal actions. Finitecritic feedback differs by**0.005133**, far above ANEOC extra error**0.00000447**.

Passing0.001prediction budgets does not make total near-optimality error<0.001.Exact roots retain critic/configuration gaps. At larger parameters unconstrained optima violate actions and are no longer compared as constrained optima.

## Implementation/reproduction checks

- Analytical sensitivities agree with closed-root central differences.
-96-pointGauss costs agree with independent adaptive integration.
- Actual candidate errors stay below analytical bounds.
- Structured/matrix counterparts agree.
- Fixed-seed errors, acceptance, sensing summaries/settings repeat exactly; retain varying timing records.
- Initial checks found inactive collocation roots need not imply unconstrained continuous-optimal feasibility; comparison ranges were corrected and failed outputs excluded.

## Limits

Ideal sensing, local predictions and finitecritic accuracy are separate levels. FullMHE, multidimensional impedance and complete route remain unvalidated. This report sets no future priorities.

## Files

- [Driver](run_quick_checks.py),[summary/reproduction](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/summarize_results.py).
- [Latest results](results.json),[first results](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/results_run1.json),[log](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/run_log.txt).
- [Root steps](root_steps.csv),[acceptance](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/acceptance.csv).
- [Sensing trials](sensing_trials.csv),[sensing summary](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/sensing_summary.csv).
- [Reproduction record](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/reproducibility.json).
