# Candidate 1, v3: the T3 0.5 s / 4 mm negative, the fixes tried, and a fresh seed-5000 test

Code: `candidate1_v3.py`. It imports `candidate1_model.py` and `candidate1_v2.py`, and neither file was changed.
Numbers are in `candidate1_v3_results.json`:
- `shared_tuning`
- `nondegradation`
- `final_seed5000`
- `known_behaviour`
- `development_notes`, which holds the diagnostic numbers quoted below
- `protocol`

## Protocol

- **Development.** Every development step and every selection used train (seed 1000+) and dev (seed 2000+) only.
  - The non-degradation check used train 40 and dev 40 scenes; the fast-push cells used dev 60.
  - The shared search used 40 train scenes per cell. Train was enlarged from v2's 20 scenes.
- **Seed 4000.** Its traces were re-simulated once to confirm the mechanism, read only. Nothing was selected on them.
- **Final test.** A new cohort with seed 5000, 40 scenes per cell, all six cells, run once after every choice was
  fixed. These were recomputed on the same scenes and noise seeds:
  - B1, with gains from `gate_PC125.json`;
  - EtR-A, tuned per cell as in candidate1_v2 (`adaptive_results.json` for T3 and T1, v2 `etrtune` for the other 4);
  - B3b at 0.1 s;
  - the previous C1v2.
- **Paired CIs.** All paired 95 % CIs come from `tune.paired_ci`.
- **Compute.** At most 2 worker processes; numpy/BLAS/numba threads set to 1.

## Mechanism of negative (1): T3 0.5 s ramp with 4 mm breathing

**This is a configuration-selection failure, not a breathing-model failure.**

1. **Selection.** The 4 mm cell had no v1 tuning, so v2 ran the 49-configuration C1 search on 20 train scenes.
   - Only 8 of the 49 were peak-feasible (at most 1 of 20 scenes above 7 N). One or two scenes therefore decided the
     choice.
   - The search space did not contain the transferable 1 mm configuration. That configuration is the hand-set default
     with `probe_v = 0.003`; the default itself, with `probe_v = 0.01`, had 2 of 20 train scenes above 7 N.
   - The selected setting was: Fdith 0 (no excitation), sk 0.3 (twice the stiffness random walk), ki 10, tu 0.08,
     lam 0.9, Nw 12, wI 40, cP 4, pD 0.1.
2. **Estimation.**
   - Without the dither, the window cannot separate stiffness from surface motion while the probe follows the
     breathing surface. The robot motion is nearly collinear with the drift basis {1, τ, τ²}.
   - With the faster random walk on top, k_hat drifts high. In quiet periods k_hat/k_true has a median of:

     | Setting | Tissue | dev | train | seed 4000 |
     |---|---|---|---|---|
     | v2 setting | soft | 2.53 | 2.12 | 2.25 |
     | v2 setting | bone | 1.40 | 1.37 | 1.36 |
     | 1 mm configuration | soft | 1.11 | 1.14 | 1.13 |
     | 1 mm configuration | bone | 1.05 | 1.02 | 1.04 |

   - The p90 of the v2 setting on soft tissue is 4-5.
3. **Control, traced on dev scene 6.**
   - On bone (21.6 kN/m) k_hat is 39 kN/m.
   - At push onset, k_used jumps to 184 kN/m through the spike plus the conservative release. D goes to 1465 N s/m and
     KPF to 1.2e-5.
   - The probe cannot follow, and F overshoots to 10 N.
   - The long feedforward ramp then retracts to u_ff -4.4 mm while the surface moves at 25 mm/s. The result is
     contact loss, then the fallback.
   - With the 1 mm configuration on the same scene, k_hat is 22-24 kN/m and F peaks at 7.1 N.
4. **The cause is the selected setting, not breathing.**
   - The same v2 setting is equally bad on the 1 mm cell: dev J 0.58, against 0.12 for the 1 mm configuration.
   - One-factor swaps (dev J; the base, the 1 mm configuration, is 0.165) show which parameters hurt:

     | Swapped parameter | dev J |
     |---|---|
     | sk 0.3 | 0.57 |
     | tu 0.08 | 0.45 |
     | ki 10 | 0.45 |
     | lam 0.9 | 0.34 |
     | Fdith 0 | 0.33 |
     | others | 0.10-0.23 |

5. **The drift model needs no periodic term.**
   - Hypothesis checked: the constant-velocity plus acceleration drift model mis-fits sinusoidal breathing. **It does
     not.**
   - The cubic remainder of 4 mm, 0.33 Hz breathing over the 0.13 s window is about 1e-8 m.
   - The drift-velocity error on bone is ±0.3 mm/s (p10/p90) at both 1 and 4 mm.
   - So no periodic drift term or adaptive window is needed.
6. **Seed-4000 diagnostic (read only).** The 1 mm configuration gives J 0.171 and loss in 2.5 % of scenes. The v2
   setting gives 0.739 and 17.5 %, and B1 gives 0.666.

## Fix

**C1v3 drops per-cell re-tuning for breathing amplitude.** Each 4 mm cell uses the v1 configuration tuned on its 1 mm
counterpart:
- T3_4 takes the T3 configuration.
- T1_4 takes the T1 configuration.
- T3fast4 takes the T3fast configuration, which v2 had already selected.

Nothing else changes: estimator, OCP, v2 supervisor, k0 = 3e4, and the cell's B1 gains as fallback. The class defaults
reproduce v2 exactly; this was checked scene by scene.

Why this fix: breathing amplitude is a patient dose variable, not a design variable. The 1 mm configurations were
tuned under the same push regime. The failure came from the small-sample selection, not from a model that is missing
4 mm physics.

## Non-degradation check (train / dev mean J; the number of scenes with J > 20 in brackets)

| Cell | split | C1v2 | C1v3 | C1v3-shared | B1 |
|---|---|---|---|---|---|
| T3 0.5 s, 1 mm | train 40 | 0.166 | 0.166 | 0.171 | 0.733 |
| T3 0.5 s, 1 mm | dev 40 | 0.116 | 0.116 | 0.283 | 0.343 |
| T3 0.5 s, 4 mm | train 40 | 0.497 | 0.677 (1) | 0.187 | 0.766 |
| T3 0.5 s, 4 mm | dev 40 | 0.755 | 0.165 | 0.348 | 0.367 |
| T3 0.25 s, 1 mm | train 40 | 2.729 (1) | 2.729 (1) | 5.053 (2) | 3.435 (1) |
| T3 0.25 s, 1 mm | dev 60 | 1.534 | 1.534 | 3.894 (2) | 1.838 |
| T3 0.25 s, 4 mm | train 40 | 2.752 (1) | 2.752 (1) | 4.024 (1) | 4.568 (2) |
| T3 0.25 s, 4 mm | dev 60 | 2.478 (1) | 2.478 (1) | 3.668 (1) | 2.275 |
| T1 r30, 1 mm | train 40 | 0.029 | 0.029 | 0.030 | 0.139 |
| T1 r30, 1 mm | dev 40 | 0.031 | 0.031 | 0.032 | 0.139 |
| T1 r30, 4 mm | train 40 | 0.044 | 0.041 | 0.038 | 0.154 |
| T1 r30, 4 mm | dev 40 | 0.049 | 0.044 | 0.040 | 0.151 |

**Rule met.**
- T3 0.5/1, T1/1, T3 0.25/1 and T3 0.25/4 have identical configurations to C1v2.
- T1/4 improves: train −0.0035 (−0.006, −0.001), dev −0.0045 (−0.015, +0.004).

**Caveat found with the enlarged train set.** On T3_4 train scenes 21-40, C1v3 has one runaway (scene 36, J 21.1).
- That scene is the residual failure type: a loss on bone, followed by a 0-12 N, 15.6 Hz bounce under the B1 fallback
  for about 0.6 s.
- Its median is 0.092 against 0.212 for C1v2.

## Single shared configuration (reported as a separate row)

**The search** (`shared_tune`) evaluated 51 configurations:
- 3 anchors: the v1 T3, T1 and T3fast configurations;
- the 32 C.SPACE random draws;
- 16 one-parameter perturbations of the best.

**Scoring.**
- Score: the mean over the six cells of log(J / J_B1), on 40 train scenes per cell.
- Constraint: peak-feasible in the four non-fast cells.
- Fallback: the cell's B1 gains.

**Selected: `{probe_v 0.003, Fdith 0.05}`.** This is the T3 configuration with the T1 dither level. Train J:
T3 0.171, T3_4 0.187, T3fast 5.05, T3fast4 4.02, T1 0.030, T1_4 0.038.

**Anchor scores.**

| Anchor | Score |
|---|---|
| T3 configuration | −0.70 |
| T3fast configuration | −0.39 |
| T1 configuration | infeasible on T3 |
| Selected configuration | −0.93 |

No single configuration works in both push regimes:
- the slow-push configurations lose contact in fast pushes;
- the fast configuration (Fdith 0, cP 1, wI 40) gives J 0.35 on T3 and 0.13 on T1.

## Fresh test, seed 5000, 40 scenes (J / RMSE N / fraction of scenes > 7 N / fraction of scenes with contact loss / J > 20)

No method had a plant divergence.

| Cell | C1v3 | C1v3-shared | C1v2 | B1 | EtR-A | B3b 0.1 s |
|---|---|---|---|---|---|---|
| T3 0.5 s, 1 mm | 0.184 / 0.185 / 0.025 / 0.000 / 0 | 0.213 / 0.190 / 0.025 / 0.050 / 0 | 0.184 / 0.185 / 0.025 / 0.000 / 0 | 0.379 / 0.290 / 0.025 / 0.025 / 0 | 0.489 / 0.301 / 0.150 / 0.100 / 0 | 0.187 / 0.185 / 0.000 / 0.025 / 0 |
| T3 0.5 s, 4 mm | 0.262 / 0.203 / 0.100 / 0.025 / 0 | 0.196 / 0.193 / 0.025 / 0.025 / 0 | 1.022 / 0.344 / 0.125 / 0.125 / 0 | 0.398 / 0.299 / 0.025 / 0.025 / 0 | 0.450 / 0.295 / 0.075 / 0.050 / 0 | 0.200 / 0.192 / 0.025 / 0.000 / 0 |
| T3 0.25 s, 1 mm | 2.468 / 0.610 / 0.600 / 0.400 / 1 | 6.146 / 0.859 / 0.725 / 0.550 / 4 | 2.468 / 0.610 / 0.600 / 0.400 / 1 | 2.288 / 0.682 / 0.725 / 0.550 / 0 | 7.976 / 0.998 / 0.725 / 0.600 / 4 | 1.995 / 0.529 / 0.400 / 0.375 / 0 |
| T3 0.25 s, 4 mm | 3.545 / 0.694 / 0.625 / 0.400 / 2 | 4.716 / 0.794 / 0.650 / 0.550 / 2 | 3.545 / 0.694 / 0.625 / 0.400 / 2 | 2.926 / 0.744 / 0.775 / 0.575 / 0 | 7.487 / 1.059 / 0.750 / 0.675 / 4 | 2.053 / 0.535 / 0.400 / 0.375 / 0 |
| T1 r30, 1 mm | 0.033 / 0.080 / 0 / 0 / 0 | 0.032 / 0.078 / 0 / 0 / 0 | 0.033 / 0.080 / 0 / 0 / 0 | 0.224 / 0.179 / 0.075 / 0.025 / 0 | 0.128 / 0.153 / 0.025 / 0 / 0 | 0.098 / 0.129 / 0 / 0 / 0 |
| T1 r30, 4 mm | 0.047 / 0.098 / 0 / 0 / 0 | 0.040 / 0.092 / 0 / 0 / 0 | 0.055 / 0.105 / 0 / 0.025 / 0 | 0.274 / 0.199 / 0.125 / 0.025 / 0 | 0.202 / 0.202 / 0 / 0 / 0 | 0.156 / 0.163 / 0.025 / 0 / 0 |

### C1v3 paired 95 % CIs (method minus reference)

| Cell | ΔJ vs B1 | ΔJ vs EtR-A | ΔJ vs B3b 0.1 s | ΔJ vs C1v2 | Δfrac>7N vs B1 | Δfrac>7N vs EtR-A | Δloss-scenes vs B1 | wins vs B1 |
|---|---|---|---|---|---|---|---|---|
| T3 0.5/1 | **−0.195 (−0.250, −0.145)** | **−0.305 (−0.455, −0.184)** | −0.003 (−0.035, +0.024) | 0 (identical) | 0 | −0.125 (−0.225, −0.025) | −0.025 (−0.075, 0) | 38/40 |
| T3 0.5/4 | −0.136 (−0.264, +0.053) | **−0.188 (−0.248, −0.127)** | +0.062 (−0.041, +0.227) | **−0.760 (−1.555, −0.213)** | +0.075 (0, +0.175) | +0.025 (0, +0.075) | 0 (−0.075, +0.075) | 38/40 |
| T3 0.25/1 | +0.180 (−0.635, +1.215) | −5.51 (−10.53, −2.02) | +0.474 (−0.562, +1.646) | 0 (identical) | −0.125 (−0.275, +0.025) | −0.125 (−0.225, −0.025) | −0.150 (−0.275, −0.050) | 29/40 |
| T3 0.25/4 | +0.619 (−0.622, +2.109) | −3.94 (−7.26, −1.40) | +1.49 (+0.10, +3.23) | 0 (identical) | −0.150 (−0.325, 0) | −0.125 (−0.250, 0) | −0.175 (−0.300, −0.075) | 27/40 |
| T1/1 | **−0.191 (−0.335, −0.083)** | −0.096 (−0.147, −0.057) | −0.065 (−0.098, −0.036) | 0 (identical) | −0.075 (−0.175, 0) | −0.025 (−0.075, 0) | −0.025 (−0.075, 0) | 38/40 |
| T1/4 | **−0.228 (−0.414, −0.090)** | −0.155 (−0.207, −0.109) | −0.109 (−0.167, −0.059) | −0.008 (−0.018, +0.001) | −0.125 (−0.225, −0.025) | 0 | −0.025 (−0.075, 0) | 38/40 |

The old negative reproduces on these fresh scenes. On T3 0.5/4, C1v2 has:
- ΔJ vs B1 +0.624 (+0.007, +1.573);
- ΔJ vs EtR-A +0.572 (+0.033, +1.386);
- contact loss in 12.5 % of scenes.

C1v3 on T3 0.5/4, scenes above 7 N: 3 of its 4 have peaks of 7.03-7.08 N, and one has 10.9 N. That last scene is
also 9.95 N for EtR-A and 16.1 N for C1v2.

### C1v3-shared paired 95 % CIs

| Cell | ΔJ vs B1 | ΔJ vs EtR-A | ΔJ vs B3b 0.1 s | ΔJ vs C1v2 | Δfrac>7N vs B1 | Δfrac>7N vs EtR-A | wins vs B1 |
|---|---|---|---|---|---|---|---|
| T3 0.5/1 | −0.166 (−0.230, −0.103) | −0.277 (−0.433, −0.136) | +0.026 (−0.020, +0.097) | +0.029 (−0.004, +0.084) | 0 | −0.125 (−0.225, −0.025) | 37/40 |
| T3 0.5/4 | −0.202 (−0.267, −0.144) | −0.254 (−0.386, −0.163) | −0.005 (−0.039, +0.025) | −0.827 (−1.752, −0.217) | 0 | −0.050 (−0.125, 0) | 37/40 |
| T3 0.25/1 | **+3.86 (+1.06, +7.60)** | −1.83 (−6.94, +2.55) | +4.15 (+1.34, +7.87) | +3.68 (+1.53, +6.54) | 0 (−0.15, +0.15) | 0 | 21/40 |
| T3 0.25/4 | **+1.79 (+0.43, +3.34)** | −2.77 (−6.21, −0.06) | +2.66 (+1.12, +4.58) | +1.17 (+0.21, +2.21) | −0.125 (−0.25, 0) | −0.10 (−0.25, +0.05) | 23/40 |
| T1/1 | −0.192 (−0.335, −0.084) | −0.096 (−0.148, −0.057) | −0.066 (−0.099, −0.038) | −0.001 (−0.002, +0.001) | −0.075 (−0.175, 0) | −0.025 (−0.075, 0) | 40/40 |
| T1/4 | −0.234 (−0.424, −0.094) | −0.162 (−0.216, −0.113) | −0.116 (−0.176, −0.064) | −0.014 (−0.025, −0.007) | −0.125 (−0.225, −0.025) | 0 | 40/40 |

## Negative (2): soft homogeneous tissue

Tested on N1 dev scenes (rng 70, 8 × 15 kPa); the declared check below uses rng 7.

**Decomposition** (RMSE in N):

| Setting | RMSE |
|---|---|
| C1 with dither | 0.046 |
| without dither | 0.040 |
| B1 | 0.024 |

So the dither is about a quarter of the gap. Most of the rest comes from the soft-side gain regime (D at the
20 N s/m floor, KPF at its 1e-3 cap). Less aggressive soft-side settings recover at most 0.004 N and do not help T3.
Single changes, with soft RMSE and T3 dev-20 J (base 0.135):

| Change | Soft RMSE | T3 dev-20 J |
|---|---|---|
| Dmin 50 | 0.045 | 0.136 |
| Dmin 100 | 0.044 | 0.136 |
| Dmin 200 | 0.044 | 0.147 |
| KPF cap 3e-4 | 0.043 | - |
| KPF cap 1e-4 | 0.042 | - |
| feedforward prior pu 1 | 0.040 | 0.24 |
| M 1 | 0.055 | 0.28 |
| no dither | 0.040 | 0.205 |

**Excitation-aware dither: not adopted.**
- The posterior std of log k (arrival prior included) relative to the prior-only std is 0.92-0.99 everywhere, so
  scaling the dither by it does nothing useful. Closed loop it was harmful: T3 dev-12 J 0.31 vs 0.14.
- The data-only std of log k per window is 0.34-0.47 with the dither (0.47 on soft N1) and 0.6-1.0 without it.
- **The window never identifies k confidently.** An on-demand dither would stay on, and switching it off removes the
  information that bone-entry detection relies on (T3 dev J 0.135 → 0.205 without the dither).
- No fix was adopted. Negative (2) remains, inside the band by a factor of about 10.

### N1 known-behaviour check (rng 7, 3 + 3 scenes, 1 mm breathing; final setting)

| Configuration | Soft 15 kPa: k_hat/k | J vs B1 | RMSE N vs B1 | Stiff 1 MPa: k_hat/k | J vs B1 | RMSE N vs B1 |
|---|---|---|---|---|---|---|
| C1v3 T3 (= T3_4) | 1.10 | 0.0086 vs 0.0028 | 0.046 vs 0.026 | 1.07 | 0.014 vs 50.2 | 0.060 vs 3.34 |
| C1v3 T1 (= T1_4) | 1.25 | 0.0060 vs 0.0022 | 0.039 vs 0.023 | 1.10 | 0.025 vs 0.069 | 0.078 vs 0.130 |
| C1v3 T3fast (= T3fast4) | 18.8 (unidentified, no dither) | 0.0054 vs 0.0140 | 0.035 vs 0.059 | 1.58 | 0.009 vs 52.1 | 0.046 vs 3.40 |
| C1v3-shared | 1.22 | 0.0064 vs 0.0028 | 0.040 vs 0.026 | 1.15 | 0.009 vs 50.2 | 0.048 vs 3.34 |

All runs had peaks of at most 4.34 N and no loss entries. The table reproduces v2 for the unchanged configurations.

## Residual failure: push on bone just after a short loss, B1 fallback bounce

**Trace.** The B1 fallback gains bounce at 0-12 N (15.6 Hz) on bone at 25-29 kN/m. The band exit (12 ticks in
[1, 6] N) cannot be met, so the bounce persists. This was seen on T3_4 train 36 and on fast-push dev scenes 6 and 43.

**Tried on fast-push train 40 / dev 60.** Values are ΔJ vs the v2 configuration as (train, dev); all were rejected.

| Variant | T3fast (train, dev) | T3fast4 (train, dev) | Other cells |
|---|---|---|---|
| Stiffness-aware fallback D_fb = max(D_B1, k_hold/wI), KPF_fb = min(KPF_B1, cP/k_hold), always | +7.6, +0.96 | +2.8, +1.7 | |
| Same, with k memory 0.2 / 0.5 s | +2.7 to +7.9 | worse | |
| KPF only | +20, +16 | | |
| Contact-gated | +2.7, +0.54 | +0.26, +0.46 | |
| Latched after first re-contact | +1.7, +0.59 | +0.35, +0.77 | T3_4: fixes train 36 (train −0.39) but dev +0.017 |
| Exit band without upper bound / band_hi 10 / n_back 6 | mixed | mixed | already part of the v2 supervisor search |

**Why they fail.**
- A stiffer fallback slows the approach during the loss: D 466 against 73 gives about 12 instead of about 60 mm/s.
  Each loss tick costs more than 20 in J.
- A lower KPF removes the fast re-contact path.

**Unresolved.**

## Configuration budget (all on train/dev or N1 dev scenes; seed 5000 run once)

| Item | Configurations | Evaluated on |
|---|---|---|
| T3_4 configuration variants | 16 | train 20 / dev 40: the v2 setting, the 1 mm configuration, probe_v 0.01, the v2 setting with dither, 12 one-factor swaps |
| T1_4 transfer | 1 | |
| N1 dither and gain variants | 17 | N1 dev, plus T3 dev 12/20 |
| Fallback and exit variants | 11 | |
| Shared search | 51 | 6 cells × 40 train |
| **Total** | **96** | |

## Judgement

- **Negative (1) is explained and fixed.**
  - It was a small-sample selection of a non-excited configuration, not a missing breathing model.
  - On fresh scenes, C1v3 on T3 0.5/4 is better than EtR-A (−0.19 (−0.25, −0.13)) and than C1v2 (−0.76).
  - It is not significantly better than B1 in J (−0.14 (−0.26, +0.05), 38 of 40 wins).
  - It has more scenes above 7 N than B1: 4 against 1, three of them at 7.03-7.08 N.
  - It is within noise of the 0.1 s oracle.
  - The enlarged train set shows one residual-type runaway, so the 4 mm cell is not fully safe.
- **The shared configuration** is as good as or better than the regime configurations on the four slow and ladder
  cells: T3_4 −0.20 vs B1, with the CI excluding 0. It is significantly worse than B1 in both fast-push cells. A
  single configuration does not span both push regimes.
- **The sourced cells reproduce on a third fresh cohort.** T3 0.5/1: −0.195 vs B1 and −0.305 vs EtR-A; T1/1: −0.19 vs
  B1.
- **The fast-push cells remain a robustness-only result.** J is not better than B1 (CIs include 0), but there are
  fewer contact-loss scenes (−0.15 and −0.175, CIs excluding 0) and fewer scenes above 7 N (CIs touching 0).
  1-2 runaways appear in 40 scenes, and there were none on seed 4000. B1 had none here.
- **Negative (2) and the fallback-bounce failure remain open.**
