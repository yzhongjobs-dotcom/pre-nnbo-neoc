# Task (round 2): why the lateral-dither signal did not match the stiffness rate — measurement or physics?

Continue from round 1 (`lateral_dither_check.py`; latest copy in the repository or in your temporary mirror `D:/pip_env/temp/lateral_dither_20261113_144358/`). Windows: `python -X utf8`. **Do not change any existing file, do not commit or push, do not edit research documents, leave `00_基线快照_2026-10-14/` alone, no dynamic arm or other plant option, no tuning, no new estimator or controller arm.** New code in a new file `lateral_dither_truth.py`. If the sandbox refuses repository writes, write to a temporary mirror of the repository paths and say where.

## Question

Round 1: the dither coefficient `c_hat` (force per unit scan displacement at the dither frequency) correlated with the stiffness rate per unit scan distance in some cases (case 38: 0.87) and not in others (case 39: 0.1–0.3), with a gain between 0.2 and 4 instead of the predicted 0.5–1. Two different things could be wrong:

1. **Measurement:** `c_hat` does not measure the plant's lateral force derivative at fixed probe height (tangential contact dynamics, the force loop reacting at the dither frequency, in-band disturbances).
2. **Physics:** it measures that derivative well, but the derivative is not the stiffness gradient: it also contains the slope of the surface relief times the stiffness, and its link to the rate of the probe-frequency stiffness depends on the tissue.

Separate the two with the plant's own static contact functions.

## Work

Rerun cases 36, 39, 38, 0 with `a_x = 0.1 mm` exactly as in round 1 (4 closed loops), keeping the full plant log. For each logged tick compute, from the scene and the plant's static contact routine (`physics_plant3.contact` / `_contact` as used for the reference force; read the code to call it correctly with the logged depth, contact position and pose; say which columns and which function you used):

- `f_x`: derivative of the static contact force with respect to scan position at fixed depth, by central difference with ±0.2 mm.
- `k_s`: static tangent stiffness at that state (derivative with respect to depth, same routine).
- `h_x`: slope of the surface relief field at the contact position (the field added to the depth as `base = y + body + h`; read the sign from the code).
- `c_static = f_x + k_s * h_x` (lateral force derivative at fixed probe height; state the sign convention and check it on one tick by perturbing the position directly).
- `dlnk_x`: derivative of `ln k_s` with respect to scan position at fixed depth.

Use the same `c_hat`, `F24`, `v24`, `G_true` as round 1. Score on ticks with `t > 1.5 s`, in contact, outside events.

Report per case:
1. Measurement: correlation of `c_hat` with `c_static` (zero lag and best lag within ±0.2 s), regression gain of `c_hat` on `c_static`, and the rms of the residual `c_hat − gain*c_static` against the noise sd predicted in round 1. Also the ratio of the dynamic to the static stiffness at 15.625 Hz if the log provides it (`|K|` column 21 against `k_s`); otherwise say so.
2. Physics: variance shares of `f_x` and `k_s*h_x` in `c_static` (and their correlation); correlation and gain of `f_x/F` with `dlnk_x`; correlation and gain of `c_static/F` with `G_true`; correlation of `dlnk_x * v24` with `d ln k_eff/dt` (how much of the stiffness rate is the spatial gradient at fixed depth, how much is depth and other effects).
3. Signal size: rms of `c_static * a_contact` (N) with the achieved contact dither amplitude, against 0.0245 N.

Predictions written before the run (report met / not met with numbers):
- T1. The measurement works: correlation of `c_hat` with `c_static` at least 0.85 (best lag) in cases 36, 38 and 0.
- T2. The gain mismatch is physics: in cases 36 and 0 the relief term carries at least 30% of the variance of `c_static`.
- T3. Case 39: the signal is small: rms of `c_static * a_contact` below 0.02 N.
- T4. `f_x/F` and `dlnk_x` agree: correlation at least 0.9 and gain between 0.4 and 1.1 in every case.

Final message (at most 1.5 pages): functions and columns used, tables for items 1–3, T1–T4 verdicts, anything that looked wrong. State plainly what failed.
