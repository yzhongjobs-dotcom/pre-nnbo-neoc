> Existing local validation record. Obsolete topics/continuation plans were removed; original configurations/numerical conclusions remain. See [HANDOFF](../../project_records/HANDOFF.md) for current plans.

# Continuous slow-variation NEOC: first-step results

2026-10-13. Continuous-root mechanism tests of appendix M.3. The subsequent [contact closed-loop step](contact_closed_loop/README.md) was executed with a reproduction subset, supporting direct NEOC and amortization. Complete R15, real-time deployment and paper-level contribution remain unvalidated; preview increments were not tested.

## Setup

Use the existing M.3 inactive branch: W+W²=η², W*=(√(1+4η²)−1)/2, u*=−W*,λ*=0. For η∈[0.65,1.05], exact roots/all candidates satisfy 0≤W<0.7. This is a finite collocation root, not an exact continuous-HJB solution.

The monotone path increases 0.65→1.05; the oscillatory path traverses four cycles in the same interval. Each uses 100,200,400,800 updates. Slow variation here refers to parameter steps, not physical time, scan speed or control period. Compare frozen sensitivity, budgeted local reuse, maintained sensitivity, maintenance with budgeted anchor updates, one Newton step and fixed old policy.

Prediction queries no exact errors, computes no candidate residuals and performs no candidate rollout/global acceptance. Exact roots serve only initial/budget-triggered anchor solves and retrospective evaluation. Local reuse triggers compare only parameter distance with a preset radius; this is not candidate-quality rechecking.

## Checked bounds

W*''=2/(1+4η²)^(3/2), with global M=2/(1+4×0.65²)^(3/2). Local error ≤M|η−η_a|²/2; a normalized W budget 0.001 permits radius √(2×0.001/M)≈0.0664.

Maintained C(W,η)=2η/(1+2W) comes from implicit differentiation, without unknown exact roots. For W≥0, its W derivative magnitude ≤4×1.05=: a, giving E_{k+1}=(1+a d_k)E_k+M d_k²/2. Long oscillatory paths make this conservative; large bounds do not imply equally large actual errors.

These certify only example-root error/inactive branches, not physical-loop performance. Floating-point errors lack reliable interval enclosure; validation tolerance is 1e−12, not a machine-certified proof.

## Results

| Path/implementation | Updates | Maximum W error | Anchor solves including initialization |
|---|---:|---:|---:|
| Monotone, initial sensitivity frozen | 100 / 800 | 0.02586 / 0.02586 | 1 / 1 |
| Monotone, maintained sensitivity | 100 / 800 | 0.0001822 / 0.00002272 | 1 / 1 |
| Monotone, budgeted local reuse | 800 | 0.0009276 | 7 |
| Four cycles, maintained sensitivity | 800 | 0.002334 | 1 |
| Four cycles, budgeted local reuse | 800 | 0.0008003 | 33 |
| Four cycles, one Newton step | 800 | 0.00001492 | 1, plus 800 residual/correction steps |

All tested analytical bounds enclose observed errors. Over 800 monotone updates, local reuse needs initialization plus 6 new anchors and directly predicts other queries, computing sensitivity once per anchor. Maintenance improves with smaller monotone steps, but accumulates oscillatory error. Newton is more accurate; retain this strong comparator.

**Continuous slow-variation direct NEOC has testable support; earlier large-jump failure did not invalidate the core idea.** Sensitivity cannot remain frozen indefinitely or be maintained for free. This conclusion concerns scalar roots; physical-loop results are in the subdirectory.

Scalar formulas for roots/sensitivities do not establish high-dimensional NNBO speed advantages. Fewer anchors need not reduce total time. Fixed-W/root differences are not evidence of worse physical fixed impedance.

## Reproduction

Run `python -X utf8 run_continuous_check.py` here. Outputs: [summary](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98NEOC%E9%AA%8C%E8%AF%81_2026-10-13/summary.csv),[complete trajectories](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98NEOC%E9%AA%8C%E8%AF%81_2026-10-13/trace.csv) and [results/configuration](results.json). See [reproduction record](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E6%8C%81%E7%BB%AD%E7%BC%93%E5%8F%98NEOC%E9%AA%8C%E8%AF%81_2026-10-13/reproducibility.json) for two-run hash checks. No random sampling or new dependency.
