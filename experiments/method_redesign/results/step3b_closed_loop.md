# Step 3b: closed-loop replacement

Candidates were frozen using only validation prediction MSE before ID/OOD or closed-loop evaluation. All other controller parameters remain FINAL. No confidence-dependent feed-forward attenuation was introduced. All errors are true force RMSE in N.

| Family | Model | Median / mean | Peak>50 N (incl divergence) | Long loss >0.5s | Paired mean Δ [95% CI] |
|---|---|---:|---:|---:|---|
| S2 | Proposed | 0.124/0.281 | 0/50 | 2/50 |  |
| S2 | P3 | 1.690/1.537 | 0/50 | 0/50 | +1.255 [+1.060, +1.443] |
| S2 | P1 | 0.131/0.252 | 0/50 | 2/50 | -0.030 [-0.086, +0.008] |
| S2 | TrueK | 0.109/0.111 | 0/50 | 0/50 | -0.170 [-0.305, -0.059] |
| S2 | Oracle+exc | 0.562/0.603 | 0/50 | 0/50 | +0.322 [+0.210, +0.418] |
| T3m | Proposed | 0.151/0.157 | 0/50 | 0/50 |  |
| T3m | P3 | 0.665/0.735 | 0/50 | 1/50 | +0.578 [+0.487, +0.671] |
| T3m | P1 | 0.148/0.160 | 0/50 | 0/50 | +0.003 [-0.000, +0.007] |
| T3m | TrueK | 0.140/0.142 | 0/50 | 0/50 | -0.015 [-0.019, -0.011] |
| T3m | Oracle+exc | 0.611/0.631 | 0/50 | 0/50 | +0.474 [+0.425, +0.523] |
| T2m-soft | Proposed | 0.077/0.078 | 0/50 | 0/50 |  |
| T2m-soft | P3 | 0.104/0.108 | 0/50 | 0/50 | +0.030 [+0.024, +0.037] |
| T2m-soft | P1 | 0.079/0.081 | 0/50 | 0/50 | +0.003 [+0.001, +0.005] |
| T2m-soft | TrueK | 0.101/0.101 | 0/50 | 0/50 | +0.024 [+0.022, +0.025] |
| T2m-soft | Oracle+exc | 0.072/0.072 | 0/50 | 0/50 | -0.006 [-0.008, -0.004] |
| T2m-stiff | Proposed | 0.283/2.305 | 2/50 | 1/50 |  |
| T2m-stiff | P3 | 1.796/2.806 | 1/50 | 27/50 | +0.502 [-1.490, +1.597] |
| T2m-stiff | P1 | 0.280/1.378 | 1/50 | 1/50 | -0.926 [-2.780, +0.002] |
| T2m-stiff | TrueK | 0.263/0.284 | 0/50 | 0/50 | -2.021 [-5.066, -0.027] |
| T2m-stiff | Oracle+exc | 0.086/0.091 | 0/50 | 0/50 | -2.214 [-5.268, -0.210] |
| T2f-stiff | Proposed | 0.597/1.085 | 0/50 | 16/50 |  |
| T2f-stiff | P3 | 2.600/2.807 | 3/50 | 44/50 | +1.722 [+1.526, +1.916] |
| T2f-stiff | P1 | 0.538/3.227 | 1/50 | 17/50 | +2.143 [-0.096, +6.524] |
| T2f-stiff | TrueK | 0.493/3.331 | 6/50 | 3/50 | +2.246 [-0.007, +5.727] |
| T2f-stiff | Oracle+exc | 0.098/0.102 | 0/50 | 0/50 | -0.982 [-1.230, -0.752] |

## Dose response

| Kind / dose | Proposed | P3 | P1 | TrueK | Oracle+exc |
|---|---|---|---|---|---|
| ramp / 0.6 | 0.243/0.336 (0/20) | 0.990/1.216 (0/20) | 0.242/0.333 (0/20) | 0.228/1.216 (1/20) | 0.087/0.091 (0/20) |
| ramp / 0.45 | 0.232/0.286 (0/20) | 1.008/1.206 (0/20) | 0.260/8.087 (1/20) | 0.213/0.224 (0/20) | 0.080/0.084 (0/20) |
| ramp / 0.3 | 0.302/0.474 (0/20) | 1.882/2.096 (0/20) | 0.368/0.533 (0/20) | 0.262/0.307 (0/20) | 0.082/0.090 (0/20) |
| ramp / 0.2 | 0.996/1.479 (1/20) | 3.001/3.110 (2/20) | 1.102/1.576 (1/20) | 0.613/4.467 (3/20) | 0.102/0.108 (0/20) |
| ramp / 0.15 | 2.421/14.497 (5/20) | 3.504/3.745 (5/20) | 2.547/5.696 (4/20) | 3.598/9.333 (8/20) | 0.118/0.119 (0/20) |
| ramp / 0.1 | 2.497/6.120 (5/20) | 3.750/5.346 (7/20) | 2.061/5.879 (4/20) | 4.787/6.786 (10/20) | 0.102/0.106 (0/20) |
| ratio / 3.0 | 0.080/0.082 (0/20) | 0.064/0.096 (0/20) | 0.079/0.077 (0/20) | 0.095/0.096 (0/20) | 0.174/0.176 (0/20) |
| ratio / 10.0 | 0.105/0.112 (0/20) | 1.473/1.320 (0/20) | 0.110/0.130 (0/20) | 0.101/0.100 (0/20) | 0.418/0.425 (0/20) |
| ratio / 20.0 | 0.136/0.275 (0/20) | 1.466/1.499 (0/20) | 0.134/0.174 (0/20) | 0.108/0.110 (0/20) | 0.621/0.594 (0/20) |
| ratio / 30.0 | 0.135/1.122 (2/20) | 1.728/1.477 (0/20) | 0.151/3.521 (2/20) | 0.115/0.118 (0/20) | 0.647/0.723 (0/20) |
| ratio / 45.0 | 0.158/1.299 (2/20) | 1.899/1.652 (0/20) | 0.159/0.575 (0/20) | 0.114/0.116 (0/20) | 0.827/0.897 (0/20) |
| ratio / 60.0 | 0.254/4.468 (2/20) | 2.046/1.907 (0/20) | 0.769/1.237 (1/20) | 0.126/0.126 (0/20) | 1.078/1.135 (0/20) |

## Measurement-gain tolerance

Linear analysis reproduces ff_loop.py substep dynamics, inverse force filter, encoder interpolation and one-command delay, with FINAL force gains. P1 freezes each IMM mode and recomputes its steady Kalman gain for each κ_g; arbitrary switching is outside this analysis. P2 uses its complete 63-tap FIR realization. The accepted log-gain band b≤ln1.2 corresponds to κ_g∈[1/1.2,1.2]. Intervals are grid observations, not analytic stability guarantees.

| Predictor | Common grid interval | Entire accepted band stable? |
|---|---|---|
| Proposed | [0.9250000000000004, 1.1250000000000004] | False |
| P1 | [0.9500000000000004, 1.1000000000000005] | False |
| P2 | [0.9500000000000004, 1.0750000000000006] | False |

Gain stress scales inverse measurement stiffness ONLY inside drift prediction; feedback stiffness, excitation normalization and gain acceptance stay FINAL. It multiplies the existing gain error rather than imposing known true stiffness. P3 is assessed empirically at 0.75/0.82/1.18/1.25; this cannot certify a continuum or convert an estimator bias bound into a stability proof.

| Gain factor | Model | Runaway | Long loss>0.5s | Median / mean |
|---|---|---:|---:|---:|
| 0.75 | Proposed | 3/250 | 23/250 | 0.179/0.727 |
| 0.75 | P3 | 5/250 | 81/250 | 1.548/1.832 |
| 0.75 | P1 | 2/250 | 23/250 | 0.188/0.692 |
| 0.82 | Proposed | 3/250 | 20/250 | 0.162/0.623 |
| 0.82 | P3 | 7/250 | 77/250 | 1.437/1.755 |
| 0.82 | P1 | 3/250 | 21/250 | 0.167/0.631 |
| 1.18 | Proposed | 20/250 | 42/250 | 0.478/7.225 |
| 1.18 | P3 | 4/250 | 66/250 | 1.128/1.482 |
| 1.18 | P1 | 32/250 | 33/250 | 0.416/10.334 |
| 1.25 | Proposed | 51/250 | 89/250 | 3.176/15.250 |
| 1.25 | P3 | 3/250 | 61/250 | 1.023/1.417 |
| 1.25 | P1 | 92/250 | 92/250 | 3.224/28.183 |

## Computation

Whole-controller hot serial timing (ms), one diagnostic scene/family, includes stiffness estimation, predictor and control, excludes plant and cold compilation:

| Model | Median P50 | P99 of scene P99 |
|---|---:|---:|
| Proposed | 0.338 | 0.551 |
| P3 | 0.421 | 0.502 |
| P1 | 0.379 | 0.512 |
| TrueK | 0.341 | 0.438 |
| Oracle+exc | 0.012 | 0.018 |

## Statistical boundaries

Diagnostic seeds110000+,50/family; sweeps exactly reuse sweep.py ramp doses0.6/0.45/0.3/0.2/0.15/0.1s and ridge ratios3/10/20/30/45/60,20/dose, seeds160000+/170000+. No final cohort accessed. Paired bootstrap resamples scenes within family,10000 replicates, seed219901; CI is percentile95%, exploratory and without multiplicity correction. Means/medians include finite runaway trajectories. Divergent pairs would be excluded only from RMSE/CI and explicitly counted, never from safety counts. TrueK/Oracle+exc are privileged references, not feasible competitors. Long loss means maximum contiguous true force<0.5N duration>0.5s. Runtime in parallel workers is retained for auditing but serial hot timing is used for interpretation. No evaluated outcome feeds selection or retuning.

Pooled individual-call hot serial timing, no concurrent batch jobs (ms):

| Model | P50 | P99 | n calls |
|---|---:|---:|---:|
| Proposed | 0.324 | 0.404 | 4750 |
| P3 | 0.385 | 0.482 | 4750 |
| P1 | 0.403 | 0.467 | 4750 |
| TrueK | 0.326 | 0.442 | 4750 |
| Oracle+exc | 0.012 | 0.024 | 4750 |

The legacy FINAL IMM remains only to preserve recovery velocity; new FF mean comes entirely from the replacement. Whole-controller timings include this legacy update and the new predictor. Standalone predictor timing excludes the legacy recovery tracker. P0 is NumPy, P1/P0r filtering is Numba; speed comparisons concern these implementations, not intrinsic algorithm complexity.

## Feedback-sector dependence

| Predictor | κ=0.5/1 common interval | Accepted band covered? | κ=0.5/1/2 common interval |
|---|---|---|---|
| Proposed | [0.7250000000000002, 1.2250000000000005] | True | [0.9250000000000004, 1.1250000000000004] |
| P1 | [0.8750000000000003, 1.2250000000000005] | False | [0.9500000000000004, 1.1000000000000005] |
| P2 | [0.9500000000000004, 1.18] | False | [0.9500000000000004, 1.0750000000000006] |

The legacy approximate[0.75,1.2] tolerance is consistent with conservative feedback κ≤1, not the full grid includingκ=2. The existing ff_loop.json already has a coarse all-grid common interval[0.95,1.1]; our finer grid gives[0.925,1.125]. This is an assumption boundary, not a newly discovered baseline failure. P1 fails full b≤ln1.2 coverage even in the conservative subset; P3 has no analytic certificate and suffers contact losses in gain stress. The tested gain multipliers act on the existing estimated gain, so they are not the actual κ_g of every tick.
