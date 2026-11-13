# Task: closed-loop test of V2 with a trajectory slope prior

Work in `experiments/rossi_failure_study/`. Reuse the runner, scenes, metrics and stored rows of `trajectory_loop.py` and `trajectory_robust.py` (arms `huber`, `delag2`; blocks: 30 original cases, 12 fresh A_large, 12 fresh A, 12 late push/pull runs), latest copies in the repository. Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no extra arms.** New code in a new file `v2_trajectory_prior_loop.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

The derivation and the predictions are in `research/v2_trajectory_prior_20261113.md`; read it first.

## Arm `v2t` (ungated host, `CC.own_mechanism(pol)`, as `delag2`)

Wrap the Huber-D4 estimator. At each call:

1. Call the original `estimate(o, n, policy)` → `(k_w, sigma, chi)`. If `k_w` is not finite return the triple unchanged.

2. Rebuild the window exactly as `d_robust.RobustRidge.estimate` does (same preprocessed log `est._log(o, n)`, same `j0`, `m`, same call of `common._build`, same trimming after a broken contact), giving `Y` (= h0) and `z`. Build the time-weighted regressor `h1` with `huber_windows.time_weighted` on the same log, `j0`, `m` and trimming (time relative to the window end `te`). Final Huber weights `w = est.weights` (length `m-4`), `D` = 4th difference, `M = I - inv(I + est.rho * D' diag(w) D)`.

       a = Y'MY;  b = Y'M h1;  c1 = h1'M h1;  p0 = Y'Mz;  p1 = h1'Mz
       m0 = p0/a;  Jv = c1 - b*b/a;  vt = (p1 - b*m0)/Jv;  tau_w = -b/a
       s_eff = est.sigma_scale * policy.params['noise'];  sv = s_eff/sqrt(Jv)

   Check on every call that `m0` equals the unclipped stiffness of the original estimate (compare with `k_w` when `k_w` is strictly inside the clip range) to 1e-8 relative; count and report violations. If `Jv <= 0`, `tau_w <= 0` or anything is not finite, fall back to the `delag2` output for this tick.

3. Trajectory prior, with the `delag2` histories and gain learner unchanged (`y6`, `y24`, `F24`, `C_w = F24/k_w`, `beta` as learned up to the previous tick):

       r  = clip(1 + beta*(y6 - y24)/C_w, 0.5, 2.0)
       vp = k_w*(1/r - 1)/tau_w
       Sp = 0.5*abs(vp)
       g  = Sp^2/(Sp^2 + sv^2)            # 0 when vp = 0
       vh = vp + g*(vt - vp)
       k  = clip(k_w + tau_w*vh, 0.5*k_w, 2.0*k_w), then clip to [200, 3e5]
       sig = sqrt(sigma^2 + tau_w^2*(1 - g)*Sp^2)
       return (k, sig, chi)

   Skip the correction (return the original triple) when `F24 < 1` or fewer than 24 ticks of history exist. Then update `beta` exactly as in `delag2`. With `beta = 0` the arm must return `k_w` exactly; check on one run with `beta = 0` throughout.

## Runs

`v2t` on: the 30 original cases; the 12 fresh A_large scenes; the 12 fresh A scenes; the 12 late push/pull runs of round-3 block 5 (same scenes, seeds, event shift and both signs). One closed loop each, 8 workers. Reuse stored `huber` and `delag2` rows.

## Output

`results/v2_trajectory_prior_loop.json`: one row per run (block, scene identity, sign, arm, L, rmse, peak, over_15, loss_time). Traces for cases 36, 37, fresh N5, V5 and all 12 late push/pull runs to `results/v2_trajectory_prior/` with `t`, `fax`, `k_eff`, `t_est`, `k_w`, `k_ret`, `g`, `vt`, `vp`, `sv`, `tau_w`, `beta`.

Final message (at most 2 pages):
1. Per group (A, A_slow_thin, B, A+B, A_large original, A_large fresh, fresh A): mean/median L for `huber`, `delag2`, `v2t`. Per case for A_large original and fresh.
2. Identification for `delag2` and `v2t` on original fast A_large (36, 37, 39, 40) and fresh A_large, returned stiffness: median |e| / mean e on plateau, rise, fall, reversal; sustained under-estimate median/worst.
3. Late push/pull: per run L for the three arms, peak force, `beta` before the event, and over the event window (start to end + 0.6 s): force rms for the three arms, mean and minimum of `ln(k_ret/k_eff)` for `delag2` and `v2t`, mean `g`. Split toward / away from the probe.
4. Mechanism: mean of `g`, median of `tau_w`, median of `sv/k_w` on (a) change ticks (rise, fall, reversal) of the traced A_large runs and (b) event-window ticks of the late push/pull runs with `beta > 0` before the event.
5. Verdicts on predictions V1–V4 of the research note with numbers; the `m0` check count; the `beta = 0` check.

State plainly what failed.
