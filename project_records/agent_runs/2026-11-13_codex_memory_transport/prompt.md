# Task: closed-loop test of a log-domain memory with trajectory transport

Work in `experiments/rossi_failure_study/`. Start from `trajectory_loop.py` and `trajectory_robust.py` (arms `huber`, `delag2`, the fresh A_large scenes of round-3 block 1) and `probe_band_loop.py` (`BandEstimate`). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning.** New code goes in a new file `memory_transport_loop.py` that imports what it needs. If the sandbox refuses repository writes, write to a temporary mirror that mirrors the repository paths and say where.

## Scientific question

The stateless trajectory correction (`delag2`) removes the window's lag on large rib crossings but does nothing about noise or about the collapse of the estimate (case 7). A memory would reduce noise and ride through uninformative windows, but a plain memory lags on stiffness changes, and an earlier memory that kept its own level locked up on an over-estimate. Question: **does a scalar log-domain Kalman memory, whose state is moved by the trajectory and whose uncertainty grows when it is moved, reduce noise and remove the collapse without losing the rib-crossing benefit and without locking up?**

## Estimator (per estimator call; all arms use the ungated host as in `trajectory_loop.py`)

Inputs from the wrapped window estimator: `(k_w, sigma, chi)`. Histories as in `delag2`: `y6`, `y24`, `F24`, `C_w = F24/k_w`, and the gain `beta` learned exactly as in `delag2` (same update, same constants).

State: `lam` (ln stiffness now) and its variance `P`. Constants: `Q0` (per-tick process variance), `RHO = 0.5`, `NOV = 24`.

At each call, in this order:

1. Transport of the state by the trajectory since the previous call (skip on the first call after a reset, or when `F24 < 1`):

       u = -ln( clip(1 + beta*(y6 - y6_prev) / (F24/exp(lam)), 0.8, 1.25) )
       lam = lam + u
       P   = P + Q0 + (RHO*u)^2

2. Measurement, when `k_w` and `sigma` are finite and `rel = sigma/k_w <= 3`:

       d = ln( clip(1 + beta*(y6 - y24)/C_w, 0.5, 2.0) )      # the delag2 correction, in log form
       z = ln(k_w) - d
       R = NOV*rel^2 + (RHO*d)^2
       K = P/(P+R);  lam = lam + K*(z - lam);  P = (1-K)*P

   Use the `beta` learned up to the previous tick for both steps, then update `beta` as in `delag2`.

3. Initialisation or re-initialisation (first finite `k_w` with `rel <= 1`, or after `n` did not advance by one, or after more than 25 consecutive calls without a finite `k_w`): `lam = ln(k_w)`, `P = (2*rel)^2`.

4. Return `(k, k*sqrt(P), chi)` with `k = clip(exp(lam), 200, 3e5)`. Before initialisation return the original triple.

## Arms

- `kf`: the estimator above on the unchanged Huber-D4 window output, `Q0 = 5e-4`.
- `kf_long`: same, `Q0 = 1e-4` (a planned second memory length, not a tuning sweep).
- `kf_noT`: `kf` with `beta` forced to 0 (plain memory, no transport).
- `band_kf`: `kf` fed by the probe-band estimate instead of the Huber-D4 window: use `probe_band_loop.BandEstimate(huber=False)` but **without** its `REL_MAX` rejection (return the value and its sigma whenever `a > 0` and `k` is finite and positive; the Kalman step weights it by its sigma). Keep everything else of that class.

Reuse stored rows for `huber` and `delag2`.

## Runs

The 30 original cases (0–23, 36–41) and the 12 fresh A_large scenes of round-3 block 1 (same scenes and seeds), each arm once, 8 workers.

## Predictions written before the run (report met / not met with numbers)

- K1. `kf` keeps the rib-crossing benefit: original A_large mean L at most 0.195; fresh A_large mean at most 0.265.
- K2. `kf` removes the collapse of case 7: L below 0.5.
- K3. `kf` leaves the other original normal cases (0–23 except 7) within ±0.02 of `huber`, case by case.
- K4. Transport matters: on original and fresh A_large together, `kf_noT` has a mean L at least 0.03 above `kf`.
- K5. No lock-up: in no run of any arm is `k_ret/k_eff > 3` for more than 0.5 s in total after t = 0.5 s.
- K6. `band_kf` is at least as good as `kf` on case 7 and within ±0.02 of `kf` in group means elsewhere.
- K7. Identification with `kf`: plateau median |e| at least 15% below `huber` (pooled over the 30 original cases); with `kf_long` at least 35% below; median sustained under-estimate of `kf` below that of `huber` in every group.

## Output

`results/memory_transport_loop.json`, one row per run (arm, case or scene identity, L, rmse, peak, over_15, loss_time, osc, ripple, and the K5 statistic: seconds with `k_ret/k_eff > 3` after 0.5 s, and max ratio). Save traces only for cases 7, 36, 37 and fresh scenes N5, V4, V5 (all new arms) to `results/memory_transport/`, with `t`, `fax`, `k_eff`, `t_est`, `k_w`, `k_ret`, `P`, `beta`.

Final message (at most 2 pages):
1. Per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh): mean/median L for `huber`, `delag2`, `kf`, `kf_long`, `kf_noT`, `band_kf`.
2. Per case for A_slow_thin, A_large original and A_large fresh: L for the six arms.
3. Identification, `e = ln(k_ret/k_eff)`, with the segment masks of the earlier rounds: per group and arm, median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst. A_large split into fast and slow as before; fresh A_large separate.
4. Ripple and force rms on crossings per group and arm.
5. K1–K7 verdicts with numbers; the K5 statistic for the worst run of each arm.

State plainly what failed; no tuning, no extra arms or analyses.
