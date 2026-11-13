# Task (round 5): recursive memory with a coherent transport error

Reuse the runner, scenes and metrics of the earlier rounds (`memory_transport_loop.py` with its `kf` arm, `chi_hold_loop.py` for the fresh A scenes; latest copies, from the repository or your temporary mirrors). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `kf2_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Why

`kf` (round 1) was the best memory overall but failed on case 37: the per-tick transport carried the state to 18 times the truth because the learned gain was about three times too large early in the run, and the filter added the transport uncertainty in quadrature per tick although a gain error is common to all ticks. `kf2` carries the relative gain error as a second state, so the transport uncertainty grows coherently and the window measurements can correct it.

## Arm `kf2` (ungated host, Huber-D4 window, `delag2` histories and gain learner unchanged)

State `x = [lam, eps]` (ln stiffness now; relative error of the transport gain), covariance `P` (2×2). Constants: `Q0 = 5e-4`, `QE = 1e-5`, `RHO = 0.5`, `NOV = 24`. With `beta` as learned up to the previous tick, at each call:

1. Transport (skip on the first call after a reset or when `F24 < 1`):

       u = -ln( clip(1 + beta*(y6 - y6_prev) / (F24/exp(lam)), 0.8, 1.25) )
       s = clip(1 + eps, 0, 2)
       lam = lam + s*u
       Fm = [[1, u], [0, 1]];   P = Fm P Fm' + diag(Q0, QE)

2. Measurement, when `k_w`, `sigma` are finite and `rel = sigma/k_w <= 3`:

       d  = ln( clip(1 + beta*(y6 - y24)/C_w, 0.5, 2.0) )
       nu = ln(k_w) - lam - s*d
       H  = [1, d];   R = NOV*rel^2
       S = H P H' + R;  K = P H'/S;  x = x + K*nu;  P = (I - K H) P (I - K H)' + R K K'
       eps = clip(eps, -1, 1)

3. Initialisation or re-initialisation (first finite `k_w` with `rel <= 1`; `n` not advancing by one; more than 25 consecutive calls without a finite `k_w`): `lam = ln(k_w)`, `eps = 0`, `P = diag((2*rel)^2, RHO^2)`.

4. Return `(k, k*sqrt(P[0,0]), chi)` with `k = clip(exp(lam), 200, 3e5)`; before initialisation return the original triple. Then update `beta` as in `delag2`.

## Runs

- The 30 original cases, the 12 fresh A_large scenes and the 12 fresh A scenes, once each.
- Perturbation: cases 7, 36, 37 and fresh N4, N5, V4, V5 with the returned stiffness multiplied by `(1+e)`, `e` = +0.005 and −0.005 (after the arm's processing).

Reuse stored rows for `huber`, `delag2`, `gated`, `kf`.

## Predictions written before the run (report met / not met with numbers)

- D1. Collapses removed: case 7 and fresh V4 below 0.5.
- D2. Rib-crossing benefit kept: original A_large mean L at most 0.20; fresh A_large mean at most 0.25; no A_large scene worse than `delag2` by more than 0.05.
- D3. Original cases 0–23 except 7 and the 12 fresh A scenes: within ±0.02 of `huber`, case by case.
- D4. Mean L over the 42 earlier scenes below 0.18.
- D5. No lock-up: no run has `k_ret/k_eff > 3` for more than 0.5 s in total after t = 0.5 s.
- D6. Under ±0.5% perturbation each of the seven scenes stays within 0.05 of its unperturbed `kf2` L.

## Output

`results/kf2_loop.json`, one row per run (arm, scene identity, perturbation, L, rmse, peak, over_15, loss_time, seconds with ratio above 3 and max ratio). Traces for cases 7, 36, 37 and fresh N4, N5, V4, V5 (unperturbed) to `results/kf2/` with `t`, `fax`, `k_eff`, `t_est`, `k_w`, `k_ret`, `P00`, `eps`, `beta`.

Final message (at most 1.5 pages):
1. Per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh, fresh A, all 42): mean/median L for `huber`, `delag2`, `gated`, `kf`, `kf2`.
2. Per case for A_slow_thin, A_large original and fresh: L for the five arms; for `kf2` the median of `eps` after 3 s.
3. Identification for `delag2` and `kf2` per group (A, A_slow_thin, original fast A_large, fresh A_large): median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst.
4. Perturbation table.
5. D1–D6 verdicts with numbers.

State plainly what failed.
