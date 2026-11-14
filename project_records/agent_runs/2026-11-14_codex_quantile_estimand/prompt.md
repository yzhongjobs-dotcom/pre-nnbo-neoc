# Task: closed-loop test of returning an upper quantile of the stiffness posterior

Work in `experiments/rossi_failure_study/`. Start from `v2_trajectory_prior_loop.py` (arm `v2t`) and reuse the runner, scenes, metrics and stored rows of that script and of `trajectory_loop.py` / `trajectory_robust.py` (30 original cases, 12 fresh A_large scenes, 12 fresh A scenes). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `quantile_estimand_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

Background and predictions: `research/v2_representation_adjudication_20261114.md`, section 4 and the predictions at the end (read them first).

## Arms (ungated host, as `v2t`)

`v2t` returns `(k, sig, chi)`, with `sig` its posterior standard deviation. New arms change only the returned stiffness:

- `v2q05`: return `(k*exp(0.5*sig/k), sig, chi)`.
- `v2q10`: return `(k*exp(1.0*sig/k), sig, chi)`.
- `v2c`: return `(k*C, sig, chi)` with one constant `C` for all runs, equal to the mean of `exp(sig/k)` over all estimator ticks of the stored or rerun `v2t` runs on the 30 original cases (compute it first and report it; it makes the average margin of `v2c` equal to that of `v2q10`).

Apply the factor after all of `v2t`'s own processing and clipping, then clip to [200, 3e5]. The gain learner inside `v2t` keeps using the window output as before.

Check: `v2q` with factor 0 must reproduce the stored `v2t` L exactly on cases 0, 7, 36 (report max |dL|).

## Runs

`v2q05`, `v2q10`, `v2c` on the 30 original cases, the 12 fresh A_large scenes and the 12 fresh A scenes. One closed loop each, 8 workers. Reuse stored rows for `huber`, `delag2`, `v2t`.

## Output

`results/quantile_estimand_loop.json`, one row per run (scene identity, arm, L, rmse, peak, over_15, loss_time, ripple). No traces.

Final message (at most 1.5 pages):
1. The constant `C`, and the median and 90th percentile of `exp(sig/k)` in the `v2t` runs.
2. Per group (A, A_slow_thin, A_slow_thin without case 7, B, A+B, A_large original, A_large fresh, fresh A, all 42 earlier scenes): mean/median L for `huber`, `v2t`, `v2q05`, `v2q10`, `v2c`; mean ripple rms for `v2t` and `v2q10`.
3. Per case for A_slow_thin, A_large original and fresh: L for the five arms.
4. Identification with the returned stiffness, `e = ln(k_ret/k_eff)`, for `v2t`, `v2q10`, `v2c`, per group: mean e on plateau, rise, fall; sustained under-estimate median/worst (masks and conventions of the earlier rounds).
5. Verdicts on predictions Q1–Q5 of the research note with numbers; the reproduction check.

State plainly what failed.
