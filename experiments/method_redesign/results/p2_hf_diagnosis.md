# P2 high-frequency diagnosis: why original C1 wins the static span (S2) median

Diagnosis cohort (110000+), 20 scenes per family, nominal platform; error of the true contact force, RMS per band
from Welch PSD (median over scenes). Scripts: `p2_c1_decompose.py`, `p2_s2_spectrum.py`, `p2_hf_attribution.py`
(argument selects the variant set), `p2_z_noise.py`, `p2_calibrated_imm.py`.

## 1. Where the gap is (S2, median RMSE N)

| | total | slow (<2.5 Hz) | 20.8 Hz ripple | rest |
|---|---|---|---|---|
| Final method | 0.132 | 0.035 | 0.037 | 0.123 |
| Final, excitation 0.10 N | 0.124 | 0.036 | 0.023 | 0.119 |
| C1-GN (excitation 0.10 N) | 0.088 | 0.031 | 0.037 | 0.073 |
| C1-GN, excitation 0.17 N | 0.102 | 0.030 | 0.052 | 0.078 |

Bands (Hz) 0-5 / 5-10 / 10-23 / 23-62.5:

| S2 | 0-5 | 5-10 | 10-23 | 23-62.5 |
|---|---|---|---|---|
| Final | 0.038 | 0.069 | 0.077 | 0.076 |
| True measurement gain | - | 0.026 (5-10 Hz of 2-band split) | - | 0.045-0.054 |
| C1-GN | 0.027-0.031 | 0.032 | 0.051 | 0.005-0.017 |
| Sensor noise off (diagnosis) | 0.075 | 0.059 | 0.068 | 0.012 |
| Proportional term K_PF1 0.5 | 0.079 | 0.106 | 0.067 | 0.079 |
| Feed-forward without servo-lag lead | 0.119 | 0.013 | 0.049 | 0.025 |

T3m: the >23 Hz bands carry ~44% of the final method's error variance (0.075/0.071 of 0.156).

## 2. Lead filtering and its scheduling (median/mean RMSE N; r runaways, L long losses)

| Variant | S2 | T3m | stiff push | fast push |
|---|---|---|---|---|
| Final | 0.133/0.339 | 0.156/0.163 | 0.302/2.958 | 1.119/1.230 L10 |
| Lead low-pass 3 Hz | 0.083/0.126 | 0.292/0.311 | 0.740/3.514 | 2.410/2.447 |
| Lead low-pass 10 Hz | 0.104/0.149 | 0.151/0.163 | 0.384/3.040 | 1.986/1.918 |
| 3 Hz while true kappa_g outside [0.87,1.15] (oracle) | 0.105/0.128 | 0.224/0.260 | 0.427/3.182 | 2.026/2.060 |
| 6 Hz while footprint spans a modulus boundary (oracle) | 0.123/0.142 | 0.179/0.195 | = final | = final |
| Posterior shrinkage v^2/(v^2+4P_vv) | 0.129/0.442 | 0.151/0.161 | 0.294/2.959 | 0.823/1.169 |
| Two-band Wiener gain, minimum-statistics floor (B=5, 60 ms) | 0.109/0.222 | 0.152/0.159 | 0.304/2.961 | 1.180/1.910 |
| Selection {raw,10,6,3 Hz} by z prediction at L+tau, 0.3 s | 0.078/0.137 | 0.163/0.169 | 0.345/2.998 | 1.568/1.773 |
| Selection, servo-exact criterion, 0.3 s | 0.093/0.115 | 0.161/0.168 | 0.306/2.959 | 1.477/1.515 |
| Selection, servo-exact, 1.0 s | 0.112/0.148 | 0.148/0.157 | 0.305/2.963 | 1.476/1.492 |
| Selection, counterfactual force error | 0.133/0.338 | 0.247/0.271 | 0.738/3.455 | 2.496/2.449 |
| Notch 20.8 Hz on z (r=0.8) | 0.113/0.119 | 0.176/0.183 | 0.479/3.159 | 1.935/1.886 |

## 3. Noise model of the motion measurement (10 per family)

| family | var(HP z)/(sigma_F/k)^2 | var/R_tracker | 20.8 Hz share of >10 Hz | NIS median/mean |
|---|---|---|---|---|
| S2 | 2.10 | 15.6 | 0.05 | 1.45/3.6 |
| T3m | 2.92 | 21.6 | 0.04 | 1.25/3.4 |
| T2m-soft | 0.76 | 5.6 | 0.06 | 1.61/3.5 |
| T2m-stiff | 7.01 | 51.9 | 0.03 | 0.80/2.2 |

Online-estimated R (minimum statistics of the second difference) with a 3x3 process-noise grid on the tuning seeds
(190000+, 10 per family): best mean log(1+J) 0.737 vs final 0.605; larger maneuver noise destabilises the
feed-forward loop (S2 22-28 N).
