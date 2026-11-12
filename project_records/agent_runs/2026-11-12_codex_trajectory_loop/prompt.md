# Task: closed-loop test of the trajectory-corrected stiffness estimate

Work in `experiments/rossi_failure_study/`. Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option.** If the sandbox refuses to write inside the repository, write the new files to a temporary folder that mirrors the repository paths and say where; that is acceptable.

## Scientific question

Replay showed that correcting the Huber-D4 window output with the encoder trajectory removes most of the window's lag bias on large, fast rib crossings (`research/indentation_identity_20261112.md`). Replay cannot tell whether this survives in closed loop, where the estimate changes the trajectory itself. Question: **with the corrected estimate driving the unchanged controller, does force control on rib crossings improve, and do the normal scenes stay about the same?**

Predictions written before the run (report each as met / not met with numbers):
- P1. A_large fast cases (36, 37, 39, 40): mean error of the stiffness the controller uses, `ln(k_hat/k_eff)`, on rise ticks moves from about −0.13 to within ±0.06, and on fall ticks from about +0.2 to below +0.12; median |error| on fall and reversal ticks drops by at least 25%.
- P2. A_large (36–41): mean L below the ungated Huber-D4 mean and at most 0.22; no case worse by more than 0.05.
- P3. A (0–5) and A_slow_thin (6–11): group mean L within ±0.02 of ungated Huber-D4.
- P4 (expected weakness). B (12–17) and A+B (18–23): L gets worse; quantify.
- P5. On A_large the privileged-gain arm is at least as good as the learned-gain arm.

## Template and harness

Follow `probe_band_loop.py` exactly for building and running an arm: `CC.cases()[index]`, `O.setup()`, `SE.MC.scene_for`, `B.build('Huber-D4', O.chosen()['Huber-D4'], sc)`, `CC.own_mechanism(pol)`, then wrap `pol.estimator.estimate`, `P3.episode(sc, pol, case['seed'])`, `M.metric(log, pol, sc)`. Stored references are in `results/completion_closed_loop.json` (`Huber-D4_nogate`, `Huber-D4_imm`, `oracle_keff_off`).

## Arms

**huber.** Unchanged ungated Huber-D4. Must reproduce the stored `Huber-D4_nogate` L for every case run (report max |dL|; it must be 0 or the harness is wrong, stop and report).

**traj (the candidate).** Wrap the estimator: call the original `estimate(o, n, policy)` to get `(k_w, sigma, chi)`, then pass them through the online filter below and return `(k_c, sigma*k_c/k_w, chi)`. When the original returns a non-finite `k_w`, return it unchanged.

The online filter is a tick-by-tick port of `dual_channel_v3_replay.run` with parameters `(qs, qo, qbeta) = (1.0, 1.0, 1e-7)` and `f0 = 0.25`. Keep its arithmetic exactly; only the data source changes:
- `y[j]` = `o.log[n-1, 2]` (the encoder row the policy sees, not the estimator's preprocessed log); `yw[j]` = mean of the last 24 such values (available prefix at start, like `D.trailing`).
- `fs[j]` = mean of the last 6 values of `o.log[n-1, 5]` (available prefix at start).
- `m[j] = ln(clip(k_w, 200, 3e5))`; `R[j] = 24*(sigma/k_w)^2`; `valid[j]` = finite `k_w` and `sigma/k_w <= 1`; `finite[j]` = finite `k_w`.
- The estimator is called once per control tick with increasing `n`; keep the histories the filter needs (`fs[j-12]`) in the wrapper. If a tick is skipped or `n` does not advance by one, reset the filter.
- Output `k_c = clip(exp(out[j]), 200, 3e5)`; if `out[j]` is not finite return `k_w`.

**traj_priv (diagnosis, A_large 36–40 only).** No learning: `k_c = F_now / (C_w + beta*(y_now - yw))` with `C_w = fs[j-12]/k_w`, `F_now = fs[j]`, `y_now = y[j]`, `beta = 1/g`, `g` = 12, 8, 3, 1.5, 3 for cases 36, 37, 38, 39, 40; if the denominator is below `0.3*C_w` use `0.3*C_w`. These gains were chosen from truth in replay, so this arm is only a ceiling for the mechanism.

## Required equivalence check (before any closed loop)

Feed the online filter class, tick by tick, with the stored arrays of cases 0, 36 and 39 (`results/huber_windows`: `y = log[:,2]`, force `log[:,4]`, `k_w = raw`, and `sigma = rel*k_w` with `rel` from `dual_channel_replay.preprocess`), and compare its output with `dual_channel_v3_replay.run(case, (1.0, 1.0, 1e-7), 0.25)[0]`. Report the max absolute difference; it must be below 1e-9.

## Runs

Cases 0–23 and 36–41 (groups A, A_slow_thin, B, A+B, A_large), arms `huber` and `traj`; `traj_priv` on 36–40. One closed loop per (arm, case), 8 workers. No tuning, no extra arms, no cough or combo.

## Output

Per run save `results/trajectory_loop/case_{index:02d}_{arm}.npz` with per-tick arrays: `t` (col 0), `fax` (col 18), `k_eff` (col 26), `y` (col 2), the stiffness the estimator returned to the controller each tick it was called (`t_est`, `k_ret`, and for `traj` also `k_w` and the filter's `g = 1/beta`), and `pol.trace` as an array.

Per run report: L = ln(1+J), rmse, peak, `over_15`, `loss_time`, `osc`, `ripple` from `M.metric`.

Identification error in closed loop, `e = ln(k_ret/k_eff)` at the estimator ticks (truth interpolated to `t_est`), scored for `t > 1.5`, `k_eff > 0`, `fax > 0.5`, with the segment definitions of `dual_channel_replay.masks_for` applied to this run's own `k_eff` (plateau, slow, rise, fall, reversal; event ticks separate). Per group and arm: median |e| and mean e per segment, and median/worst sustained under-estimate (`dual_channel_replay.sustained`).

Force control on crossings: per group and arm, rms of `fax − FD` over the ticks of the rise, fall and reversal segments together, and over plateau ticks.

Tables in the final message (at most 2 pages):
1. Per case: group, L for `huber`, `traj`, `traj_priv` (where run), stored `Huber-D4_imm` and `oracle_keff_off`; peak force for `huber` and `traj`.
2. Per group: mean and median L per arm, and the two force-rms numbers.
3. Per group: the identification table above (A_large split into fast 36, 37, 39, 40 and slow 38, 41).
4. P1–P5 verdicts with numbers, the harness regression and the equivalence check.

Deliverables: `experiments/rossi_failure_study/trajectory_loop.py` (one file), `results/trajectory_loop.json` with all rows, the npz traces. State plainly what failed; do not tune or add analyses to rescue a prediction.
