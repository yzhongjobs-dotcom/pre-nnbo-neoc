Continuing. Two jobs, in order.

## 1. Physics of the arm

The plant uses a quasi-static arm: depth from x = y - (F + m*y_ddot)/k_r. Add a dynamic option where the tool is a mass on the arm spring: m x_ddot = k_r(y-x) + c_r(y_dot-x_dot) - F, with damping ratio zeta_r. No measured value exists for zeta_r: take 0.1 as nominal and also run 0.03 and 0.3 as sensitivity, and mark it an assumption. The sensor reading is then F + m*x_ddot, using the probe's actual acceleration. ForceSensing stays as is: it only knows the encoder and its calibration, and any new mismatch is real and should be reported. Keep the old equation behind a flag.

Checks before anything else:
- Without the arm option change, stored results reproduce exactly.
- Open-loop frequency response of the dynamic option against the analytic x/y = (k_r + j w c_r)/(k_r + K - m w^2 + j w c_r), from 5 to 150 Hz including near the arm's natural frequency. Integration step stable and converged at 0.5 ms: compare to 0.25 ms.
- Low-frequency limit agrees with the old equation.

Then show whether the approximation matters. On your 6 episodes' conditions, old versus dynamic (zeta 0.03, 0.1, 0.3): the gain layers 3, 5, 6, 7 at 20.8 Hz; noise-free k_hat relative to layer 3; and closed loop with gated Huber-D4 and gates-off with Kalman tau=0.2: L, peak, force RMSE. Also say whether the 2 ms fast pole matters under the dynamic arm (1 and 4 ms, gated Huber-D4 only).

Decide from this. If the dynamic arm is the correct physics and it changes results, make it the default, with zeta 0.1, and say so; if it changes nothing that matters, say so and still report which is default. Do not retune anything.

## 2. Rescore against the effective stiffness

Truth k_eff = |G_3(j w0)| at 20.8 Hz, the continuous axial force to encoder gain under the plant in use. Log or compute it per tick from the plant's own state, using the formula matching the arm option, including geometry and support. Check it against a noise-free lock-in at several operating points, and report the error. Keep |K| logged too.

Using the plant chosen in part 1 and the replay method of your identification-layer study, identical observations, generate 12 episodes: scene A base 874500 with nominal i = 0, 2, 6 and slow_thin i = 1, 4, 5, each with (`neural_white`, 28) and (`viscous`, 14). Estimators: gated, raw, Kalman tau 2 and 0.2, IMM, Lepski, EWA, forgetting with floor, Huang ARD fixed.

Report per method with e = ln(k_est/k_eff):
- plateau: bias, robust sd, RMSE;
- transitions, rising and falling, large and small, with counts: lag, settling, mean error, largest under-estimate, fraction of ticks below 1 and below 0.85;
- the same plateau and transition under-estimate numbers against the old truth col 21, so the change is visible;
- the gated estimator's bias split into the 1.054 centering and the rest.
Update the noise-versus-tracking figure and the curve figures for two slow_thin episodes, in `research/figures/id_layer_keff/`.

Report in Chinese, under 2 pages, final message only. State what was measured, and any check that failed. No method suggestions.

Scripts: `experiments/rossi_failure_study/plant_arm_check.py`, `keff_rescore.py`, results in `results/`. No commit, no push.
