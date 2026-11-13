# Task: does a small scan-direction dither expose the spatial stiffness gradient? (phase-1 physics check)

Work in `experiments/rossi_failure_study/`. Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no closed-loop use of the new signal.** New code in a new file `lateral_dither_check.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Scientific question

V2 (time-varying stiffness in the window) fails under noise because the window holds almost no information about the stiffness rate: the rate only modulates the 0.03 N probe response, giving a rate noise of 2.5–4.3 /s per window against true rates of 1–8 /s. Physics offers another route. With the probe at fixed depth, moving it by `dx` along the scan changes the contact force by `dF = (∂f/∂x) dx ≈ κ F (∂ln k/∂x) dx` (κ between about 0.5 and 1). A small known dither along the scan, superimposed on the scan motion, should therefore produce a force response proportional to the local spatial stiffness gradient, and `d ln k/dt = v_scan · ∂ln k/∂x`. It should be zero on uniform tissue and insensitive to body motion, which is not correlated with the dither.

This phase only asks whether the plant's physics shows that signal, how large and how clean it is, and what it costs. The controller and estimator stay unchanged.

## Setup

Use the harness of `trajectory_loop.py` (ungated Huber-D4 host, `CC.own_mechanism(pol)`, unchanged estimator and controller). Before running an episode, add to the scene's commanded scan trajectory

    sc['xcmd'] = sc['xcmd'] + a_x * sin(2π f_x t) * ramp(t),     f_x = 15.625 Hz (8 control ticks per cycle)

with `ramp` rising smoothly from 0 to 1 over the first 0.3 s, on the same time grid as `sc['xcmd']`. Amplitudes `a_x` = 0, 0.05, 0.1, 0.2 mm (commanded). Check in `physics_plant3.py` how `xcmd` drives the probe (`s[1]`, first-order lag) and how the contact point follows it (the tangential spring, mass, damping and friction around `xc`), and report in two or three sentences what the plant assumes there, including whether the dither velocity can reverse the sliding direction at the scan speeds used (it should not: `2π f_x a_x` must stay below the scan speed; report the ratio per case).

Cases: 36, 39 (A_large, fast), 38 (A_large, slow), 0 (A), 12 (B, push/pull, no ribs). One closed loop per (case, amplitude): 20 runs.

## Analysis (offline, from each run's own log; truth only where stated)

Record per tick what `huber_windows.job` records so that the estimator's own regression windows can be rebuilt: `Y`, `z` for the last 24 ticks via `common._build` with the estimator's preprocessed log. Build the lateral regressor `X` the same way by calling `common._build` on a copy of that log in which column 2 (y) is replaced by the measured scan encoder column 1 (x), without the H preprocessing. For each tick fit, with `M = I − (I + ρ D'D)^-1` (4th difference, `ρ = (2 sin(π·14·Tc))^-8`, no Huber weights):

    [k_hat, c_hat] = argmin (z − k Y − c X)' M (z − k Y − c X)

and the single-regressor `k_only` fit (`c = 0`). Noise covariance of `(k_hat, c_hat)`: `σ² (R'MR)^-1` with `R = [Y X]`, `σ = 0.0245`.

Measured-only quantities: `c_hat` (N/m), `F24` (24-tick mean of the measured force), scan speed `v24` (24-tick mean of the measured x velocity).

Truth for scoring: `lam = ln k_eff` (column 26), smoothed by an 8-tick mean (one dither cycle), `G_true = (d lam/dt) / v24` wherever `v24 > 2 mm/s`.

Report per case and amplitude:
1. Achieved dither amplitude at 15.625 Hz of the probe (`x` encoder) and of the contact point (integrate column 16, the contact velocity), in mm; ratio of dither velocity amplitude to scan speed.
2. Force ripple at 15.625 Hz (rms, true axial force column 18) on plateau ticks and on changing ticks (`|d lam/dt| > 1 /s`), and at 20.833 Hz; L of the run; change of L from `a_x = 0`.
3. Gradient signal: correlation of `c_hat/F24` with `G_true` over ticks outside events with `t > 1.5 s`, the regression gain `κ` (slope of `c_hat/F24` on `G_true`), and the same with an optimal lead/lag of ±0.2 s (report the lag).
4. Noise: sd of `c_hat` on plateau ticks (`|d lam/dt| < 0.2 /s`), the sd predicted by the covariance formula, and the implied rate noise `v24 · sd(c_hat) / (κ F24)` in 1/s. For comparison compute V2's own trend statistic noise on the same windows: with `h1` = the time-weighted regressor of `huber_windows.time_weighted`, `g = h1 − (b0/a0) h0`, `J_v = g'Mg`, rate noise `σ / (k_eff sqrt(J_v))` in 1/s (median).
5. Stiffness channel: median |ln(k_hat/k_eff)| and median |ln(k_only/k_eff)| on plateau and on changing ticks, for each amplitude (does adding the lateral regressor or the dither disturb the stiffness estimate?).
6. Body motion (case 12 only): mean and sd of `c_hat/F24` inside the event window (start to end + 0.3 s) against outside.

Predictions written before the run (report met / not met with numbers):
- L1. Cases 36 and 39, `a_x = 0.1 mm`: correlation of `c_hat/F24` with `G_true` at least 0.8 (allowing the best lag), gain κ between 0.3 and 1.5.
- L2. Implied rate noise at `a_x = 0.1 mm` at most 1 /s in cases 36 and 39, and at least 3 times smaller than V2's trend noise on the same windows.
- L3. Cost at `a_x = 0.1 mm`: added force ripple at 15.625 Hz at most 0.08 N rms on changing ticks and at most 0.02 N rms on plateau ticks; L within +0.03 of the run without dither in every case.
- L4. The stiffness estimate is not disturbed: median |ln(k_hat/k_eff)| within 0.01 of the no-dither `k_only` value in every case at `a_x = 0.1 mm`.
- L5. Case 12: the mean of `c_hat/F24` inside the push/pull window stays within two plateau standard deviations of its mean outside.

## Output

`results/lateral_dither_check.json` with all numbers; a figure per case 36, 39, 12 at `a_x = 0.1 mm` to `research/figures/lateral_dither/`: `G_true` and `c_hat/(κ F24)` against time (top), `ln k_eff`, `ln k_hat` (bottom).

Final message (at most 2 pages): the plant description asked above, tables for items 1–6, L1–L5 verdicts, anything that looked wrong. State plainly what failed; no tuning, no extra arms.
