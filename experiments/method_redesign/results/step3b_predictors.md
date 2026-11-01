# Step 3b: causal motion prediction

All prediction errors below are in μm. Predictions are evaluated on the SAME FINAL closed-loop measurements; no truth enters predictor inputs. The target sign is eta=surface+reset constant=-s_eq+constant. Oracle uses future physical surface at the current observed x, not future x.

| Model | ID RMSE | ID P95 | OOD RMSE | OOD P95 | ID 1σ/2σ | P50/P99 ms |
|---|---:|---:|---:|---:|---|---|
| P0 | 2390.81 | 5347.91 | 2253.50 | 4881.87 | 0.053/0.105 | 0.198/0.258 |
| P1 | 2389.90 | 5347.71 | 2226.23 | 4874.66 | 0.051/0.104 | 0.051/0.075 |
| P2 | 2389.17 | 5339.95 | 2175.30 | 4869.85 | 0.171/0.304 | 0.036/0.057 |
| P3 | 2134.09 | 4563.11 | 1911.02 | 4114.49 | 0.710/0.896 | 0.059/0.084 |
| P0r | 2387.52 | 5341.83 | 2164.16 | 4871.56 | 0.029/0.059 | 0.036/0.050 |

## id: family and stage errors

| Group | P0 | P1 | P2 | P3 | P0r |
|---|---:|---:|---:|---:|---:|
| all / all | 2390.81/5347.91 | 2389.90/5347.71 | 2389.17/5339.95 | 2134.09/4563.11 | 2387.52/5341.83 |
| S2 / all | 3399.19/6739.91 | 3395.63/6763.03 | 3385.99/6745.23 | 2967.55/5757.78 | 3388.45/6764.05 |
| all / breathing | 2675.00/5619.05 | 2673.62/5620.20 | 2670.62/5614.26 | 2344.00/4776.43 | 2670.60/5618.22 |
| T3m / all | 3618.26/5746.07 | 3618.40/5742.94 | 3618.39/5742.80 | 2904.41/4937.99 | 3618.04/5738.88 |
| all / pre50ms | 1767.22/4505.29 | 1767.03/4525.78 | 1770.51/4555.50 | 1663.72/3477.34 | 1765.89/4572.54 |
| all / ramp | 1925.07/4859.74 | 1926.21/4860.54 | 1928.96/4831.48 | 1797.85/4026.01 | 1924.27/4861.25 |
| all / hold | 1798.54/4638.56 | 1798.47/4627.20 | 1803.64/4621.39 | 1718.34/3846.24 | 1798.02/4625.73 |
| all / retract | 1817.96/4348.89 | 1817.71/4354.10 | 1823.53/4354.00 | 1716.31/3592.87 | 1816.45/4347.38 |
| T2m-soft / all | 964.61/1775.80 | 963.97/1776.09 | 971.78/1794.75 | 1279.25/2093.48 | 960.99/1774.56 |
| T2m-stiff / all | 544.57/449.27 | 545.88/448.84 | 562.16/459.73 | 1143.95/1472.80 | 545.01/444.67 |
| T2f-stiff / all | 208.69/398.27 | 208.90/400.39 | 242.64/444.56 | 1070.95/1487.42 | 207.71/396.68 |

## ood: family and stage errors

| Group | P0 | P1 | P2 | P3 | P0r |
|---|---:|---:|---:|---:|---:|
| all / all | 2253.50/4881.87 | 2226.23/4874.66 | 2175.30/4869.85 | 1911.02/4114.49 | 2164.16/4871.56 |
| S2 / all | 3300.75/6633.64 | 3229.23/6607.79 | 3196.70/6608.33 | 2646.27/5633.42 | 3139.43/6609.90 |
| all / breathing | 2449.10/5014.57 | 2402.93/5009.46 | 2366.82/5000.69 | 2049.25/4324.38 | 2345.49/5005.75 |
| T3m / all | 3287.06/5419.55 | 3243.96/5415.75 | 3217.00/5408.91 | 2611.04/4655.97 | 3211.23/5409.40 |
| all / pre50ms | 1550.55/3855.71 | 1551.59/3866.19 | 1558.33/3802.12 | 1490.95/3168.69 | 1541.01/3795.14 |
| all / ramp | 1723.23/4289.33 | 1761.97/4285.78 | 1731.38/4277.98 | 1601.16/3487.09 | 1699.19/4274.57 |
| all / hold | 1797.80/4332.68 | 1826.86/4300.92 | 1667.81/4275.45 | 1549.31/3434.96 | 1722.51/4294.90 |
| all / retract | 1735.02/4490.90 | 1744.04/4484.25 | 1759.18/4509.02 | 1641.51/3677.98 | 1727.93/4471.43 |
| T2m-soft / all | 687.88/1292.76 | 688.66/1298.22 | 701.94/1329.36 | 1057.21/2010.13 | 684.60/1291.56 |
| T2m-stiff / all | 187.95/389.52 | 188.71/388.91 | 191.92/406.86 | 994.26/1425.98 | 177.49/388.38 |
| T2f-stiff / all | 1090.57/494.32 | 1186.41/497.16 | 753.84/623.34 | 1094.52/1513.95 | 885.13/494.99 |

## Protocol and boundaries

Training seeds 200000+ (100/family); selection 190000+ (exactly10/family); ID110000+ (20/family), OOD115000+/120000+/125000+ (20/family/profile, max seed129019). No final130000+ cohort accessed. The OOD125000 block is separate from FINAL and stops before130000. Each new model has at most61 validation looks (GP one); P0r61, P1 61, NN≤61epochs with12epoch patience and16384 sampled windows/epoch. All models see the same causal eta and measurement-variance sequences; GP/NN use63 samples, IMM recursively summarizes history. No evaluation tuning.

Selection is frozen BEFORE diagnostic evaluation. P0r is a budget-matched retuned two-CA sensitivity reference (not a new architecture); P0 stays FINAL. P1 has stationary / CA / constant-jerk modes. P2 fits SE+periodic covariance on training observations, profiles the constant mean, and predicts the exact linear functional with derivative covariance. P3 is a 63-sample temporal convolution and 1x1 Gaussian head, 3170 parameters, trained with Gaussian NLL. Early stopping is by validation NLL; closed-loop candidate selection is validation scene MSE, never diagnostic/OOD outcomes.

Calibration is tick coverage, descriptive only. P0/P1 covariances model state uncertainty but omit contact-model and stiffness bias; none is calibrated on evaluation. Reset offsets and gain error can accumulate: absolute physical-target error includes these, so scene increment RMSE is also retained in JSON. Equal causal information does not imply equal supervision: P3 uses training physical-target labels, P2 uses measurement likelihood, and IMM uses validation labels for selection. Computational budgets are explicit, rather than pretending those training operations are equivalent.

Training includes 100 scenes/family with 20 each of smoothstep, cosine, linear ramp, minimum jerk, and smoothstep with varied breathing; nonzero profiles include a second respiratory harmonic. ID is the original cohort, 20/family. OOD has trapezoidal velocity, two-stage pushes and phase-continuous frequency jumps, 20/family/profile, with matching noise/control structure. Stage labels are temporal event intervals; breathing means outside event intervals and is not necessarily free of stiffness transitions in S2/T3m. All ticks after 1 s are retained, including post-loss resets; incomplete episodes are counted in metadata.

No variance scaling of feed-forward is applied. Multiplying an absolute prediction by a time-varying confidence gain can create offset-dependent command jumps; a risk tradeoff and state-consistent increment rule would need separate design. Hot serial timings include prediction/update, exclude stiffness estimation, plant and cold compilation, and are not real-time certification.

## Offset-free diagnostics (post-selection only)

Local horizon error = (prediction-current measured eta) - (surface(t+12ms)+tau*surface_dot(t+12ms)-surface(t)). Oracle increment error = diff(prediction)-diff(physical target), excluding tracker resets. Segment-centered error removes one best-fit constant per reset segment for diagnosis only (uses future labels, never control or selection). None changes the original absolute target or validation criterion.

| Split / family | Model | Local horizon RMSE μm | Oracle increment RMSE μm | Centered RMSE μm |
|---|---|---:|---:|---:|
| id / all | P0 | 226.66 | 275.92 | 1035.72 |
| id / all | P1 | 214.76 | 262.15 | 1033.62 |
| id / all | P2 | 195.19 | 203.08 | 1032.08 |
| id / all | P3 | 897.54 | 50.69 | 1009.20 |
| id / all | P0r | 189.38 | 203.39 | 1028.08 |
| id / S2 | P0 | 353.88 | 497.12 | 1895.06 |
| id / S2 | P1 | 310.77 | 454.33 | 1888.40 |
| id / S2 | P2 | 199.52 | 219.07 | 1871.69 |
| id / S2 | P3 | 899.26 | 58.23 | 1831.35 |
| id / S2 | P0r | 228.08 | 308.04 | 1875.59 |
| id / T3m | P0 | 102.91 | 99.96 | 1182.98 |
| id / T3m | P1 | 101.34 | 94.02 | 1183.25 |
| id / T3m | P2 | 112.70 | 65.78 | 1183.61 |
| id / T3m | P3 | 901.27 | 27.58 | 1147.74 |
| id / T3m | P0r | 93.75 | 54.76 | 1182.25 |
| id / T2m-soft | P0 | 92.35 | 154.38 | 112.29 |
| id / T2m-soft | P1 | 93.36 | 150.84 | 113.11 |
| id / T2m-soft | P2 | 149.59 | 180.48 | 164.40 |
| id / T2m-soft | P3 | 885.13 | 48.41 | 236.68 |
| id / T2m-soft | P0r | 56.80 | 81.30 | 77.79 |
| id / T2m-stiff | P0 | 341.60 | 329.20 | 409.62 |
| id / T2m-stiff | P1 | 343.69 | 335.97 | 411.41 |
| id / T2m-stiff | P2 | 319.92 | 350.97 | 432.50 |
| id / T2m-stiff | P3 | 889.94 | 71.73 | 356.03 |
| id / T2m-stiff | P0r | 342.13 | 328.18 | 410.22 |
| id / T2f-stiff | P0 | 65.73 | 70.98 | 72.85 |
| id / T2f-stiff | P1 | 66.48 | 72.62 | 73.47 |
| id / T2f-stiff | P2 | 139.75 | 79.97 | 140.89 |
| id / T2f-stiff | P3 | 912.06 | 40.89 | 263.27 |
| id / T2f-stiff | P0r | 63.43 | 55.75 | 70.14 |
| ood / all | P0 | 760.99 | 1078.09 | 1292.08 |
| ood / all | P1 | 675.00 | 959.04 | 1243.69 |
| ood / all | P2 | 501.64 | 628.79 | 1149.03 |
| ood / all | P3 | 901.47 | 120.67 | 1022.20 |
| ood / all | P0r | 457.19 | 649.60 | 1130.02 |
| ood / S2 | P0 | 1147.59 | 1574.29 | 2195.50 |
| ood / S2 | P1 | 932.68 | 1266.28 | 2085.95 |
| ood / S2 | P2 | 840.20 | 1092.55 | 2033.15 |
| ood / S2 | P3 | 901.68 | 186.69 | 1795.80 |
| ood / S2 | P0r | 616.00 | 828.20 | 1947.40 |
| ood / T3m | P0 | 769.77 | 1150.65 | 1483.75 |
| ood / T3m | P1 | 562.45 | 846.75 | 1384.37 |
| ood / T3m | P2 | 415.51 | 624.19 | 1320.24 |
| ood / T3m | P3 | 898.85 | 126.22 | 1219.47 |
| ood / T3m | P0r | 374.53 | 575.98 | 1307.07 |
| ood / T2m-soft | P0 | 104.14 | 155.25 | 126.34 |
| ood / T2m-soft | P1 | 107.33 | 155.11 | 130.07 |
| ood / T2m-soft | P2 | 169.85 | 192.11 | 188.21 |
| ood / T2m-soft | P3 | 900.39 | 58.39 | 264.32 |
| ood / T2m-soft | P0r | 88.11 | 84.97 | 106.05 |
| ood / T2m-stiff | P0 | 83.75 | 126.39 | 91.50 |
| ood / T2m-stiff | P1 | 84.73 | 120.76 | 91.69 |
| ood / T2m-stiff | P2 | 91.44 | 63.27 | 99.77 |
| ood / T2m-stiff | P3 | 891.93 | 32.24 | 265.81 |
| ood / T2m-stiff | P0r | 58.22 | 78.55 | 67.99 |
| ood / T2f-stiff | P0 | 1004.85 | 1420.20 | 1039.91 |
| ood / T2f-stiff | P1 | 1091.22 | 1570.31 | 1142.07 |
| ood / T2f-stiff | P2 | 612.84 | 587.67 | 677.94 |
| ood / T2f-stiff | P3 | 916.69 | 132.06 | 381.07 |
| ood / T2f-stiff | P0r | 756.59 | 1090.19 | 822.27 |

This separates cumulative measurement/contact bias from short-horizon extrapolation, but does not provide privileged bias correction to any predictor. A lower pooled absolute error need not imply better motion compensation. GP fitting reached its61-evaluation cap without convergence; its result is a bounded prototype, not a GP performance upper bound.

## Protocol correction

An initial implementation used50 native+50 augmented validation scenes. After discovering the mismatch, selection was replayed with ONLY50 native190000+ scenes. Architecture, training inputs, all61 prespecified IMM candidates, GP optimization and network update trajectory were unchanged; no pilot ID/OOD/closed-loop result was used to add or alter candidates. Initial labels also had an implementation sign error, corrected and independently checked against Oracle before fitting the delivered models. Final numbers use corrected labels and the exact50-scene validation set. This is exploratory diagnostic work, not a preregistered confirmatory evaluation. The correction repeated computation (same candidate sets), and that cost is not hidden as a single executed pass.
