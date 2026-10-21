# Candidate 1, v2: contact-loss supervisor (design log and results)

Code: `candidate1_v2.py`, which imports from `candidate1_model.py` (none of that file's code changed).
Numbers: `candidate1_v2_results.json`.

## Protocol

- **Splits.** Development used train (seed 1000) and dev (seed 2000) only.
  - The fast-push dev split was enlarged to 60 scenes and the fast-push train split to 40.
- **Fresh final test.** Seed 4000, 40 scenes per cell. This seed had never been used before.
  - B1 (gate gains), B3b at 0.1 s (gate rules) and EtR-A were recomputed on the same scenes and noise seeds.
- **Original test (seed 3000).** Re-run for continuity. It had been seen before, so treat those numbers
  as already known.
- **Diagnosis of the failure mechanism.** I used the two diverging test traces the brief referred to
  (seed 3000, scenes 6 and 36) and two dev scenes that diverged in the same way (dev 6 and 53).
  - I never tuned on the test traces.

## Failure mechanism (T3, 0.25 s push ramp)

Traces: test 6/36 and dev 6/53/45/43, with v1 and the T3fast-tuned configuration. The sequence is:

1. **Fast push or retraction.** It is 7 mm at up to 42-44 mm/s, often on bone (k about 20-26 kN/m).
   - The window fit trades stiffness against drift: k_hat drops from 26 to 3 kN/m and vs hits its
     0.2 m/s bound.
   - The OCP then over-retracts: on dev 6 the command goes 15 mm back for a 7 mm push.
   - The result is contact loss.
2. **During the loss, three things go wrong at once:**
   - (i) The estimator keeps fitting the few pre-loss samples left in its window, so vs swings between
     +200 and -200 mm/s.
   - (ii) The local contact model is anchored at the current robot position with a force level of about
     0. The predictor therefore believes the surface touches the probe while it is really 10-40 mm away.
   - (iii) The plant force integral winds up, because ef = FD for the whole loss and ki = 30.
3. **The OCP fights the integrator.** It uses u_ff against the integrator-driven approach, which it
   predicts will hit the surface.
   - u_ff runs to its -20 mm bound and stays there.
   - The loss lasts 0.6-1.2 s.
4. **Re-contact.** It happens with the integral wound, the feedforward saturated and a window that mixes
   pre-loss and post-loss data.
   - Result: 20-100 N impacts and a bounce limit cycle (J 75-731).

Why `lossfix` (v1) could not work: it changed the OCP's loss barrier, but the OCP's information is invalid
during a loss. The integrator windup was not addressed at all.

## Development on dev (about 38 supervisor variants; no test data)

Mean J on fast-push dev (60 scenes) unless stated otherwise. For reference: v1 2.88, B1 1.84, B3b 1.08.

| Step | Result | Kept? |
|---|---|---|
| Loss mode with an approach-speed damping law (D = drive / v_app), u_ff held, estimator frozen | worse on the 6 hard scenes (40 vs 18.7) | no |
| + bumpless integrator (conditional integration through ki) | 6 hard scenes 10.2; without it 49 | **yes** (essential) |
| Estimator reset with n_settle 12 | 2.48 | yes |
| Feedforward rate clamp (0.05-0.1 m/s) | 2.48 → 2.8-3.0 | no |
| Approach relative to the pre-loss drift estimate ("chase"), v_app 0.04 | 2.26 | superseded |
| Stiffness-scaled approach speed c/k_hold (limits the impact force) | 3.2-6.1 | no |
| Bumpless gain handback | 2.15 | superseded |
| Robustness probe speed 0.03/0.05, wI 20/80 | 2.2-12.4 | no |
| **Simplex fallback to the cell's B1 admittance during loss** (until 12 ticks in a [1, 7] N band) | **1.60-1.63** | **yes** |
| Fallback also entered at F > 8-10 N | 1.75 | no |

Why the approach laws fail:
- J weighs loss heavily. Each loss tick costs 64 + 20 against about 1 for normal tracking.
- So a slow, safe approach loses more J through long losses than it saves on impacts.
- No single approach speed works on bone, where k·v·T_react is about 30 N at 35 mm/s.

Why the fallback works:
- The B1 admittance has low D, high KPF and M ≈ 1.6. It keeps or regains contact quickly and tolerates
  oscillation on bone instead of losing contact.
- C1's bone gains (D ≈ 500) are well damped in contact but follow neither a push nor a loss.

## Declared tuning

- **Supervisor.**
  - 49 configurations: the default + 32 random draws + 16 one-parameter perturbations from `SUP_SPACE`.
  - Run on fast-push train, 40 scenes. The v1 T3fast OCP configuration was held fixed.
  - Train J: v1 9.29, B1 3.44, best 2.14 (F_off 1.0, F_fb_hi 10). No setting met the peak constraint.
- **Non-degradation constraint (declared after seeing it was needed, train only).**
  - The best setting raised T3 (0.5 s) train J from 0.144 to 0.216, because F_off = 1.0 enters the
    fallback on normal force dips.
  - Rule: take the lowest fast-train J among the 15 best settings whose T3 and T1 train J do not exceed
    v1 by more than 1e-3.
  - Selected: F_off 0.3, band [1, 6] N, n_back 12, n_settle 12, aw 1, vs_max 0.1, u_rate 0, no entry on
    high force. Fast-train J 2.73.
  - Dev check: fast 1.534 (unconstrained setting 1.546; B1 1.84). T3 and T1 dev and train are
    J-identical to v1.
- **OCP configurations.**
  - T3, T1 and T3fast: v1's tuned configurations, with k0 = 3e4 as the start prior.
  - T3fast4, T3_4 and T1_4 had no v1 tuning. For these I ran the declared 49-configuration C1 tuning
    (20 train scenes) with the final supervisor. Train J: 1.96, 0.319 and 0.039.
- **EtR-A.**
  - T3 and T1: from `adaptive_results.json`.
  - The 4 other cells: the `tune_etr` procedure, 120 configurations, 20 train scenes.

## Fresh test (seed 4000, 40 scenes)

Columns:
- **J / RMSE N / >7 N:** mean J, mean RMSE in N, and the fraction of scenes with peak above 7 N.
- **loss:** the time fraction in contact loss.
- **J>20:** the number of scenes with J above 20.
- No method had a plant divergence (div 0 throughout).

| Cell | C1v2 | C1v1 (as tested) | B1 | EtR-A | B3b 0.1 s |
|---|---|---|---|---|---|
| T3 0.5 s, 1 mm | 0.153 / 0.173 / 0.025, loss 0, J>20 0 | 3.25 / 0.290 / 0.05, J>20 1 | 0.638 / 0.324 / 0.05 | 0.397 / 0.273 / 0.025 | 0.141 / 0.165 / 0 |
| T1 r30, 1 mm | 0.022 / 0.067 / 0 | 0.022 / 0.067 / 0 | 0.205 / 0.147 / 0.025 | 0.071 / 0.120 / 0 | 0.056 / 0.100 / 0 |
| T3 0.25 s, 1 mm | 2.13 / 0.596 / 0.55, loss 0.0061, J>20 0 | 5.65 / 0.713 / 0.525, J>20 2 | 2.26 / 0.670 / 0.725, loss 0.0060 | 6.49 / 0.948 / 0.675, J>20 2 | 1.41 / 0.472 / 0.425 |
| T3 0.25 s, 4 mm | 2.95 / 0.675 / 0.525, loss 0.0095, J>20 0 | - | 3.41 / 0.774 / 0.775, J>20 1 | 4.78 / 0.882 / 0.775 | 1.44 / 0.477 / 0.425 |
| T3 0.5 s, 4 mm | 0.739 / 0.327 / 0.10, loss 0.0021 | - | 0.666 / 0.334 / 0.05 | 0.390 / 0.281 / 0.05 | 0.153 / 0.173 / 0 |
| T1 r30, 4 mm | 0.038 / 0.090 / 0 | - | 0.215 / 0.163 / 0.025 | 0.130 / 0.166 / 0 | 0.089 / 0.127 / 0 |

Paired 95 % CIs for C1v2 minus the reference (`tune.paired_ci`). Wins are out of 40 scenes.

| Cell | ΔJ vs B1 | ΔRMSE vs B1 | Δfrac>7N vs B1 | ΔJ vs EtR-A | ΔJ vs B3b 0.1 s |
|---|---|---|---|---|---|
| T3 0.5/1 | -0.48 (-0.98, -0.16), 38 wins | -0.15 (-0.22, -0.10) | -0.025 (-0.10, +0.05) | -0.24 (-0.38, -0.14) | +0.012 (-0.001, +0.029) |
| T1/1 | -0.18 (-0.41, -0.04), 37 wins | -0.08 (-0.13, -0.04) | -0.025 (-0.075, 0) | -0.048 (-0.066, -0.033) | -0.034 (-0.052, -0.018) |
| T3 0.25/1 | -0.13 (-0.87, +0.67), 30 wins | -0.07 (-0.16, +0.02) | **-0.175 (-0.30, -0.075)** | -4.36 (-7.76, -2.01) | +0.72 (+0.28, +1.23) |
| T3 0.25/4 | -0.46 (-1.75, +0.73), 27 wins | -0.10 (-0.21, +0.01) | **-0.25 (-0.40, -0.125)** | -1.84 (-3.10, -0.90) | +1.51 (+0.83, +2.28) |
| T3 0.5/4 | **+0.07 (-0.57, +0.69)**, 28 wins | -0.007 (-0.10, +0.09) | +0.05 (-0.075, +0.175) | **+0.35 (+0.03, +0.79)** | +0.59 (+0.24, +1.05) |
| T1/4 | -0.18 (-0.41, -0.05), 34 wins | -0.07 (-0.12, -0.04) | -0.025 (-0.075, 0) | -0.092 (-0.118, -0.068) | -0.051 (-0.079, -0.027) |

C1v2 vs C1v1, fresh:
- T3 0.5: -3.09 (-9.28, +0.00).
- T3 0.25: -3.52 (-10.05, +0.11).

Diagnostic ablation, run after the fresh test and not used for any selection:

| Fresh cell | v2 without the supervisor (k0 = 3e4) | v2 with k0 = 5000 |
|---|---|---|
| T3 0.5/1 | 0.153 | 0.180 |
| T3 0.25/1 | 16.99, 4 runaways | 2.01 |
| T3 0.25/4 | 7.70 | 3.03 |
| T3 0.5/4 | 2.03 | 0.70 |

- v1's fresh T3 0.5 runaway is a start-prior failure (k0 = 5000), and k0 = 3e4 removes it.
- In the fast cells the supervisor is what removes the runaways.

## Original test (seed 3000, seen before)

- **T3 0.5:** C1v2 0.139 vs C1v1 0.142 vs B1 0.570.
- **T1:** 0.021 vs 0.021 vs 0.247.
- **T3 0.25:**
  - C1v2: J 1.86, RMSE 0.538, >7N 0.475, 1 scene with J > 20 (old scene 6: 23.9, B1 18.0).
  - The re-run of C1v1 gave 17.2. The originally reported value was 25.5, and these chaotic runaway
    scenes vary between runs.
  - B1: 2.27 / 0.673 / 0.65. B3b: 1.08.
  - ΔJ vs B1: -0.42 (-0.81, +0.05), 34 wins. Δfrac>7N: -0.175 (-0.30, -0.075).
  - Old scene 36 (previously 731) is now 5.2.

## N1 known-behaviour check (final setting; homogeneous tissue, 1 mm breathing, 3 scenes each)

- **T3 configuration.**
  - Soft 15 kPa: k_hat/k 1.1, J 0.0086 vs B1 0.0028, RMSE 0.046 vs 0.026 N.
  - Stiff 1 MPa: J 0.014 vs B1 50.2, peak 4.16 N.
  - No loss entries.
- **T1 configuration.** Soft: RMSE 0.039 vs 0.023. Stiff: J 0.025 vs 0.069.
- **T3fast configuration.** J is fine on both tissues: soft 0.0054 vs B1 0.014, stiff 0.009 vs 52. The
  stiffness estimate is wrong, though:
  - k_hat stays near the 3e4 start prior on soft tissue (median 20 kN/m vs 1.06 true).
  - It is ×1.6 high on stiff tissue.
  - This configuration has no dither (Fdith = 0), so homogeneous tissue gives no excitation. The
    estimate is unidentified and errs on the stiff, safe side.
- The soft homogeneous accuracy gap of open negative 3 is **not** addressed: RMSE is still 0.04-0.05 vs
  B1 0.023-0.026 N.

## What failed, and the remaining negatives

- **Every approach-speed law failed:** constant, drift-relative, stiffness-scaled, and the feedforward
  clamp. Integrator compensation alone was not enough without the fallback.
- **Fast push (0.25 s).** The runaways are gone: 0 of 40 fresh scenes in both amplitudes, against 2 for v1.
  Peak violations drop by 17.5-25 points against B1, and those CIs exclude 0.
  - J is not significantly better than B1: the CIs include 0.
  - C1v2 stays well behind the B3b oracle.
  - A residual bad scene type remains, seen on fresh scene 23 in three cells: a push on bone right after a
    short loss. The B1 fallback gains chatter at 0-21 N on bone, and the band exit (12 consecutive ticks
    in [1, 6] N) keeps the fallback active through the push.
- **New negative, T3 0.5 s with 4 mm breathing** (cell newly tuned for C1).
  - C1v2 does not beat B1: ΔJ +0.07 (-0.57, +0.69), with a lower median (0.185 vs 0.349).
  - It is significantly worse than EtR-A: +0.35 (+0.03, +0.79).
  - It has more scenes with loss (17 % vs 5 %).
  - The same scene-23 mechanism contributes, and the supervisor halves the damage (2.03 → 0.74). The OCP
    configuration tuned on 20 train scenes for this cell may also be weaker than the 1 mm one.
- **The fix borrows B1's gains as the safety controller.** C1v2 is therefore C1 plus the cell's best
  fixed controller, and B1's tuning is part of its cost.
