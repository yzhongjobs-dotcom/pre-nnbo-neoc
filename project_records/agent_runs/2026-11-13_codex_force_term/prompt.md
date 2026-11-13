# Task: does the force term of the indentation identity carry the timing that the trajectory term misses? (replay diagnosis)

Work in `experiments/rossi_failure_study/`. Use only the stored closed-loop recordings `results/huber_windows/case_00.npz … case_41.npz` (see `huber_windows.py`; `log` columns in order `[t, x, y, fn, F_sensor, cmd, u_true, fax_true, |kt|_true, k_eff_true, phase_true]`, one row per 8 ms tick; `raw` = window stiffness output). Windows: `python -X utf8`. **No plant simulation, no closed loop, do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone.** New code in a new file `force_term_replay.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Question

The identity `y = b + g F / k` gives, for the change of log stiffness over a short span during which the background `b` barely moves,

    Δ ln k = Δ ln F − Δ ln(y − b)  ≈  θ · Δ ln F − β · Δy / C,      C = F/k .

The trajectory correction used so far keeps only the second term. Hypothesis: when the stiffness changes, the force deviates first and the trajectory follows only after the controller has restored the force, so the trajectory term alone is late; that would explain why memory with trajectory transport lags on stiffness falls and why the stateless correction under-corrects rises. Test on recorded data whether the force term adds explanatory power for the stiffness change across the spans that matter (window lag and memory).

## Definitions (per case, per tick `j`, span `h` ticks, `h ∈ {9, 24, 48}`)

- `y6`, `F6`: trailing 6-tick means of `y` and `F_sensor` (removes the 20.8 Hz ripple). `lam = ln k_eff_true`, trailing 6-tick mean `lam6`.
- Target: `T = lam6[j] − lam6[j−h]`.
- Regressors: `x1 = −(y6[j] − y6[j−h]) / C[j−h]` with `C = F6 / exp(lam6)` (truth stiffness in `C`, to isolate the physics); `x2 = ln(F6[j] / F6[j−h])`.
- Measured-only variants: `x1m` with `C = F6 / k_s`, where `k_s` is a causal smooth of `raw` (exponential, time constant 0.3 s, held through gaps).
- Score on ticks with `t > 1.5 s`, `k_eff_true > 0`, `fax_true > 0.5`, both ends of the span in the same contact run; events (`meta['event_info']`, start − 0.05 s to end + 0.6 s) reported separately. Segments from truth as in `dual_channel_replay.masks_for` (plateau, rise, fall, reversal).

## Report

Groups: A_large fast (36, 37, 39, 40), A_large slow (38, 41), A (0–5), A_slow_thin (6–11); events from B (12–17) and A+B (18–23).

1. For each group and span, ordinary least squares without intercept, fitted per case and then pooled as medians across cases: R² and coefficients for `[x1]`, `[x2]`, `[x1, x2]`; the same with the measured-only `x1m`. Also the plain correlation between `x1` and `x2`.
2. Residual bias by segment for `h = 9` and `h = 24`, per group: mean of `T − fit` on rise ticks and on fall ticks, for `[x1]` and for `[x1, x2]` (per-case fits).
3. Timing: per group, the cross-correlation lag (in ms, within ±200 ms) at which `d x1/dt` (5-tick difference) best matches `d lam6/dt`, and the same for `d (x1 + x2)/dt` with the fitted coefficients, and for `x2` alone.
4. Events (B, A+B; the stiffness barely changes): rms of `x1`, `x2`, and of the fitted `[x1, x2]` combination using the median A_large-fast coefficients, at `h = 9`, inside the event windows; and the correlation between `x1` and `x2` there. (Does the force term cancel or add to the false change claimed by the trajectory term during body motion?)

Predictions written before the run (report met / not met with numbers):
- F1. A_large fast, `h = 9`: adding `x2` raises R² by at least 0.15 over `[x1]` alone, with a coefficient on `x2` between 0.3 and 2.
- F2. A_large fast, `h = 9`: with `[x1]` alone the residual is negative on rises or positive on falls by more than 0.03 in magnitude; with `[x1, x2]` both are within ±0.03.
- F3. Timing: `x1` lags the true change by at least 24 ms in A_large fast; the fitted combination lags by at most 8 ms.
- F4. Events: inside push/pull windows the combination's rms is not larger than that of `x1` alone.

Output: `results/force_term_replay.json` and a final message of at most 1.5 pages with the tables for items 1–4 and the F1–F4 verdicts. State plainly what failed; no extra analyses.
