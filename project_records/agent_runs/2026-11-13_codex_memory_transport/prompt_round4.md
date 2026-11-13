# Task (round 4): hold on inconsistent windows only

Reuse the runner, scenes and metrics of the earlier rounds (`gated_delag_loop.py`, `memory_transport_loop.py`, `trajectory_loop.py`, `trajectory_robust.py`; latest copies, from the repository or your temporary mirrors). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `chi_hold_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Why

On the ungated host the estimate collapses in case 7 and fresh V4. The full gates avoid that but cost a lot on rib crossings. In an ungated run of case 7, the windows that lead into the collapse (error −0.5 to −1) have a normal reported precision (`sigma/k_w` about 0.15) but a residual statistic `chi` of about 6–10, against about 0.2 in normal windows; on a large rib crossing (case 36) `chi` stays near 0.2–0.6 even when the error is large. So the window's own misfit marks the collapse and not the crossing. Test: keep the ungated host and hold the last accepted value only when the window is inconsistent with its model.

## Arms (ungated host: `CC.own_mechanism(pol)` as in `trajectory_loop.py`)

Read `chi_max = pol.params['chi2_max']` **before** calling `CC.own_mechanism(pol)` (it is about 2.68; `own_mechanism` then sets the parameter to infinity). A window is *accepted* when `k_w`, `sigma`, `chi` are finite, `chi <= chi_max` and `sigma/k_w <= 1`.

- `huber_chi`: wrap the estimator. Accepted window: return the original triple and remember it. Rejected window with a remembered value: return the remembered triple. Rejected with nothing remembered: return the original triple. Forget the remembered value when `n` does not advance by one.
- `delag2_chi`: the same rule around `delag2`: an accepted window goes through the unchanged `delag2` correction and its output triple is remembered; a rejected window returns the remembered triple and is stored as missing (NaN `C_w`, `rel`) for the gain learner, so it takes no part in learning `beta`.

Reuse stored rows for `huber`, `delag2`, `gated`, `gated_delag2`.

## Runs

The 30 original cases, the 12 fresh A_large scenes and the 12 fresh A scenes of round-3 block 2 of `trajectory_robust.py` (same scenes and seeds; run `huber` on the fresh A scenes if no stored row exists), both new arms once each, 8 workers.

## Predictions written before the run (report met / not met with numbers)

- C1. `delag2_chi` removes the collapses: L below 0.5 on case 7 and on fresh V4.
- C2. It keeps the rib-crossing benefit: original A_large mean L at most 0.20, fresh A_large mean at most 0.27 with V4 excluded from this mean; no A_large scene worse than `delag2` by more than 0.05.
- C3. Original cases 0–23 except 7 and the 12 fresh A scenes: `delag2_chi` within ±0.02 of `huber`, case by case.
- C4. Over the 42 scenes of the earlier rounds, mean L of `delag2_chi` is below that of `delag2` (0.220) and of `gated_delag2` (0.228).
- C5. Holds are rare and short outside collapses: in cases 0–23 except 7 the held fraction of estimator ticks is below 5% and the longest hold is below 0.5 s.

## Output

`results/chi_hold_loop.json`, one row per run (arm, scene identity, L, rmse, peak, over_15, loss_time, held fraction, longest hold in seconds). Traces for cases 7, 36, 37 and fresh V4, N4, N5 (`delag2_chi` only) to `results/chi_hold/` with `t`, `fax`, `k_eff`, `t_est`, `k_w`, `sigma`, `chi`, `k_ret`, `held`, `beta`.

Final message (at most 1.5 pages):
1. Per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh, fresh A, all 42 earlier scenes): mean/median L for `huber`, `delag2`, `gated`, `gated_delag2`, `huber_chi`, `delag2_chi` (stored where available).
2. Per case for A_slow_thin, A_large original and fresh: L for the six arms; held fraction and longest hold for `delag2_chi`.
3. Identification for `delag2` and `delag2_chi` on original fast A_large (36, 37, 39, 40) and fresh A_large: median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst (returned stiffness).
4. C1–C5 verdicts with numbers.

State plainly what failed.
