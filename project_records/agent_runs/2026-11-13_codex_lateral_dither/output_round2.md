**Both the measurement and the stiffness-rate interpretation fail in these runs.** Relief does **not** explain the large gain mismatch. Case 38 tracks the static derivative well, but with gain 2.40; cases 36, 39 and 0 fail the measurement-correlation threshold.

Four unchanged closed loops were rerun at **0.1 mm**. Their original 28 log columns exactly equal round 1. Repository writes were refused; [new code](D:/pip_env/temp/lateral_truth_20261113/experiments/rossi_failure_study/lateral_dither_truth.py), full logs and [results](D:/pip_env/temp/lateral_truth_20261113/experiments/rossi_failure_study/results/lateral_truth_20261113_150000/summary.json) are in the temporary mirror. No existing files, research documents or baseline files were changed; no commit or push.

Used `physics_plant3._contact(xc,d,theta,vc,geometry,q_field,E_field,ribs,activation)[0]`, the static force used inside the plant. Original zero-based columns: time **0**, depth **6**, contact velocity **16**, pose **17**, dynamic stiffness magnitude **21**, effective stiffness magnitude **26**. Logging-only instrumentation appended exact contact position **28**, activation **29**, static force **30**, and hereditary states **31–32**. Recomputed static force matched column 30 exactly.

Derivatives used ±0.2 mm in position and ±0.1 µm in depth, holding pose, velocity and activation fixed. Positive relief increases depth (`base=y+body+h`), giving the requested \(c_{\rm static}=f_x+k_sh_x\). At case 36, \(t=5.752\) s, this gave **−223.624 N/m**, versus **−223.584 N/m** from direct position-plus-relief perturbation. Strict fixed-height geometry includes a \(\cos\theta\) projection; its correction was small: **0.266 N/m RMS** in case 36.

Scoring used **1,062 ticks per case**: \(t>1.5\) s, contact, outside events. `c_hat`, `F24`, `v24` and `G_true` are unchanged from round 1. Gains below are OLS with intercept.

**1. Measurement**

Residuals are the requested \(c_{\hat{}}-g\,c_{\rm static}\), using best lag, without subtracting the intercept. Positive lag means the estimate follows truth. Residual and predicted noise units are N/m; prediction is the RMS of round 1’s per-tick noise SD over matched ticks.

| Case | Correlation zero / best (lag) | Gain zero / best | Residual RMS / predicted SD | \(|K|/k_s\), 15.625* / 20.8 Hz |
|---|---|---|---|---|
| 36 | .588 / .610 (+64 ms) | 1.754 / 1.826 | 2503 / 742 | 4.890 / 5.284 |
| 39 | .145 / .169 (+80 ms) | .217 / .255 | 640 / 414 | 1.216 / 1.255 |
| 38 | .956 / .967 (+72 ms) | 2.406 / 2.397 | 394 / 327 | 3.813 / 4.133 |
| 0 | .638 / .689 (+112 ms) | 3.725 / 3.953 | 353 / 362 | 5.123 / 5.554 |

*Column 21 actually uses **20.8 Hz**, not 15.625 Hz. The starred ratio reconstructs the same plant expression at 15.625 Hz using appended hereditary states. Ratios are medians. Round 1’s **plateau-median** predicted SDs were respectively **344, 344, 296, 348 N/m**.

**2. Physics**

Variance entries are \(\operatorname{Var}(f_x)\), \(\operatorname{Var}(k_sh_x)\), and \(2\operatorname{Cov}\), each divided by \(\operatorname{Var}(c_{\rm static})\). They sum to 100%; cancellation permits individual shares above 100%. Here \(F=F24\).

| Case | Variance: \(f_x\) / relief / cross (%) | Component correlation | \(f_x/F\) vs \(d\ln k_x\): r / gain | \(c_{\rm static}/F\) vs \(G_{\rm true}\): r / gain |
|---|---|---|---|---|
| 36 | 100.49 / 0.25 / −0.74 | −.074 | .992 / .425 | .939 / 1.401 |
| 39 | 101.51 / 0.87 / −2.38 | −.126 | .993 / .418 | .935 / .723 |
| 38 | 101.02 / 0.10 / −1.11 | −.178 | .993 / .450 | .906 / .935 |
| 0 | 124.12 / 10.76 / −34.88 | −.477 | .982 / .449 | .817 / 1.029 |

Let \(s=d\ln k_x\,v24\) and \(r=d\ln k_{\rm eff}/dt\).

| Case | corr(\(s,r\)) | RMS spatial \(s\) / remainder \(r-s\) / total \(r\) (s⁻¹) |
|---|---|---|
| 36 | .958 | 7.135 / 5.209 / 2.045 |
| 39 | .955 | 2.960 / 1.510 / 1.600 |
| 38 | .930 | 1.213 / .745 / .531 |
| 0 | .829 | .595 / .448 / .194 |

Depth and other effects **strongly oppose** the fixed-depth spatial contribution. The remainder also includes effective-versus-static stiffness and filtering differences; it cannot be attributed entirely to depth.

**3. Signal size**

Used round 1’s achieved-contact harmonic amplitude at each tick.

| Case | Median contact amplitude (mm) | RMS \(c_{\rm static}a_{\rm contact}\) (N) | Relative to 0.0245 N |
|---|---|---|---|
| 36 | .03833 | .03961 | 1.62× |
| 39 | .03829 | .01630 | .67× |
| 38 | .06729 | .03439 | 1.40× |
| 0 | .03848 | .00335 | .14× |

**Predictions:** T1 **not met** (.610, .967, .689 for 36/38/0); T2 **not met** (relief shares 0.25%, 10.76%); T3 **met** (.01630 < .020 N); T4 **met in every case** (r=.982–.993, gain=.418–.450).

The failed assumption is that `c_hat` reliably measures the static lateral derivative across these cases. Independently, that static derivative is not a calibrated stiffness-rate proxy. Dynamic/static stiffness differences are substantial, but these four runs do not uniquely identify the remaining measurement error’s cause.