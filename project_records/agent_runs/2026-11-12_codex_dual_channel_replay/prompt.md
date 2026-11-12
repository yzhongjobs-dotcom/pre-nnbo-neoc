# Task: replay test of a two-channel stiffness estimator (probe + indentation identity)

Work in `experiments/rossi_failure_study/`. Use only the stored closed-loop recordings `results/huber_windows/case_00.npz … case_41.npz` (42 cases; see `huber_windows.py` for how they were made). **No new plant simulation, no closed loop, no change to any existing file, no commit or push, no edits to research documents.** Leave `00_基线快照_2026-10-14/` alone. Windows: use `python -X utf8`.

## Scientific question

In constant-force scanning the encoder position satisfies the identity

    y(t) = b(t) + p(t),   p(t) = g F(t) / k(t)

- `b`: position of the tissue surface along the probe axis (breathing, surface relief, support). Unknown, smooth.
- `p`: indentation. `g` is a dimensionless shape-times-dispersion factor of the tissue, unknown, slowly varying (expected range 1–10).
- `k`: the probe-frequency stiffness the existing window estimator measures with noise.

So a stiffness change shows up twice: weakly in the 20.8 Hz probe response (existing channel, sd about 0.17 in ln k per 24-tick window), and strongly in the low-frequency encoder trajectory (millimetres). Question: **does fusing both channels in one causal filter track rib-crossing stiffness changes (rise, fall, reversal) better than the probe channel alone, without hurting plateaus?**

Predictions written before the run (report each as met / not met with numbers):
- P1. Fast large crossings (A_large cases 36, 37, 39, 40): rms error on rise/fall/reversal ticks at most half of the probe-only filter; sustained under-estimate (definition below) at most 0.12.
- P2. Plateau rms error not larger than the probe-only filter (within 0.01).
- P3. Slow crossings (A_large 38, 41; A_slow_thin) and A: within ±0.02 of the probe-only filter.
- P4. The learned `g = exp(c)` after two crossings is within a factor 1.5 of the per-episode fit of `F/k_eff` on `y + u`.
- P5 (expected weakness). During push/pull events (groups B, A+B) the fused filter is worse than probe-only; quantify.

## Data (per case)

`log` columns in order `[t, x, y, fn, F_sensor, cmd, u_true, fax_true, |kt|_true, k_eff_true, phase_true]`, one row per 8 ms tick. `raw`: Huber-D4 window output per tick (N/m, NaN if none). `Y`, `z`: the regression windows per tick (32 wide; use the last 24 columns). `meta` (JSON): group, `event_info`, etc.

The estimators may read only `t`, `y`, `F_sensor`, `raw`, `Y`. Truth columns are for scoring, and for arm A3 below.

## Common preprocessing (causal)

- `Tc = 0.008`. `ys[j]` = mean of `y[j-5..j]`, `Fs[j]` = mean of `F_sensor[j-5..j]` (six ticks = one probe cycle, removes the ripple; delay 2.5 ticks).
- `m[j] = ln(clip(raw[j], 200, 3e5))` when `raw[j]` is finite.
- Probe information: `M = I - (I + rho D'D)^-1` with `D` the 4th difference on 24 points, `rho = (2 sin(pi*14*Tc))^-8`; `rel[j] = 0.0245 / (raw[j] * sqrt(Yw' M Yw))` with `Yw` = last 24 columns of `Y[j]`. Measurement variance of `m[j]` used by every filter: `R[j] = 24 * rel[j]^2` (each tick adds 1/24 of a window). Skip the probe update when `raw` or `Yw` is not finite, `Yw' M Yw <= 0`, or `rel[j] > 1`.
- The window output is an average over its 24 ticks, and `ys` is delayed 2.5 ticks; use a single lag `LAG = 9` ticks between the probe observation and the trajectory state.
- Contact: if `Fs[j] < 1` or no finite `raw` for more than 25 ticks, stop updating; re-initialise at the next valid tick.

## Arms

**A0 raw.** `exp(m[j])`, held when missing.

**A1 probe-only IMM (the fair no-coupling comparison).** Scalar state `lam = ln k`, random walk, two modes with per-tick variances `q0` (quiet) and `q1` (transit), Markov `p_stay`; observation `m[j] = lam + noise(R[j])`. Standard IMM (mixing, mode-conditioned Kalman update, Gaussian likelihood, combination), as in `../motion_aware_scan/imm.py`.

**A2 fused filter (the candidate).** IMM of two EKFs differing only in the process variance of `ell`.
- State `x = [b, bd, bdd, ell, c]`: `b` background position with velocity and acceleration; `ell = ln p`; `c = ln g`.
- Dynamics: `b` constant-acceleration with white jerk of spectral density `q_b` (exact discretisation, as `imm.py`); `ell` random walk with per-tick variance `q0` or `q1`; `c` random walk with per-tick variance `q_c`.
- Observation 1, every valid tick: `ys[j] = b + exp(ell) + v`, `var(v) = (5e-6)^2`.
- Observation 2, when the probe update is allowed: `m[j] = c + ln Fs[j-LAG] - ln p_lag + noise(R[j])`, with `p_lag = exp(ell) - (ys[j] - ys[j-LAG]) + bd*tau - bdd*tau^2/2`, `tau = LAG*Tc` (the indentation LAG ticks ago, from the identity and the motion of the background itself). If `p_lag <= 0.05*exp(ell)` skip the update. Use the analytic Jacobian.
- Output: `lam_hat[j] = c + ln Fs[j] - ell`, `k_hat = exp(lam_hat)`, and its standard deviation from the covariance.
- Initialisation at the first valid tick: `g0 = 3`; `c = ln g0` with variance `(ln 3)^2`; `ell = c + ln Fs - m`, `b = ys - exp(ell)`; `bd` sd 0.03 m/s, `bdd` sd 0.3 m/s^2; set the covariance of `(b, ell, c)` consistently with these definitions (the errors of `ell` and `c` share the `c` error; the `b` error is `-p` times the `ell` error).
- IMM: mode likelihood uses both observations of the tick.

**A3 fused filter, background removed (privileged ceiling, diagnosis only).** Same as A2 with `ys` replaced by the six-tick mean of `y + u_true`. Shows what is left when breathing is known; surface relief is still unknown.

## Selection and scoring

- Truth: `lam_true = ln k_eff_true`. Error `e = lam_hat - lam_true` on ticks with `t > 1.5`, `k_eff_true > 0`, `fax_true > 0.5`.
- Segments from truth: `dl` = six-tick mean of `lam_true`, then `np.gradient(.., Tc)`, then five-tick mean. Plateau `|dl| < 0.2`; slow `0.4 < |dl| <= 1`; rise `dl > 1`; fall `dl < -1`; reversal = ticks within 0.1 s of a sign change of `dl` that has `|dl| > 1` on both sides within 0.15 s. Event ticks (`event_info` start minus 0.05 s to end plus 0.6 s) are a separate segment and are excluded from the others.
- Sustained under-estimate of a case: `max(0, -min over time of the 0.1 s moving average of e)` outside events.
- Development cases: even indices in groups A, A_slow_thin, A_large, B, A+B. Test cases: odd indices of the same groups. Cough and combo: report only, never used for selection.
- Selection objective (same for A1 and A2, each with its own grid): mean over development cases of `rms(e on plateau) + rms(e on rise, fall, reversal and slow ticks together)` (a missing segment contributes 0).
- Grids: A1: `q0 in {0, 1e-7, 1e-6, 1e-5}`, `q1 in {1e-4, 1e-3, 3e-3, 1e-2}`, `p_stay in {0.95, 0.99}`. A2: the same three, times `q_b in {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}` (m^2/s^5), `q_c in {1e-8, 1e-6}`. A3 uses the values selected for A2.
- Report on test cases, per group and segment, for A0–A3: rms, median |e|, mean e, 5th percentile of e; per group: median and worst sustained under-estimate. Also the same table for cough and combo (all cases), and a per-case table for A_large (36–41, mark which are test).
- P4: per case, final `exp(c)` versus the slope fit: regress `F/k_eff` (six-tick means, truth) on `y + u_true` over the episode after removing a 2 s centred moving average from both; `g_fit = 1/slope`.

## Checks you must run and report

1. Synthetic unit test of A2 before touching the recordings: `b = 0.004 sin(2 pi 0.25 t)`, `F = 8`, `g = 4`, `lam` = 7.6 with ramps of +0.7 and -0.7 over 0.4 s every 1.2 s, `y = b + g F / exp(lam)`, probe observation = `lam` delayed 9 ticks plus white noise with `R = 24 * 0.17^2`. Report rms error of A1 and A2 on ramp and flat ticks, and the learned `g`.
2. A0 reproduces `raw` exactly.
3. No truth column is read by A0–A2 (state how you enforce it).

## Deliverables

- `experiments/rossi_failure_study/dual_channel_replay.py` (one file, runnable with `python -X utf8 dual_channel_replay.py`), results in `results/dual_channel_replay.json`.
- Figures `research/figures/dual_channel/case_{36,39,0,6,12}.png`: `lam_true`, A0, A1, A2 (with the one-sd band of A2) against time.
- Final message (at most 2 pages): selected parameters, the test tables, P1–P5 verdicts with numbers, the three checks, anything that looked wrong. State plainly if a prediction fails; do not tune beyond the grids to rescue it, and do not add other analyses.
