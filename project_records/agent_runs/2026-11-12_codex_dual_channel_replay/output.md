**A2 did not reliably improve tracking: P1–P4 failed; P5 was met.** The replay and figures are complete, but Windows denied writing into the repository. Artifacts are staged in temporary storage:

[Runnable script](D:/pip_env/temp/dual_channel_replay_asodag0c/experiments/rossi_failure_study/dual_channel_replay.py) · [Results JSON](D:/pip_env/temp/dual_channel_replay_asodag0c/experiments/rossi_failure_study/results/dual_channel_replay.json) · Figures: [36](D:/pip_env/temp/dual_channel_replay_asodag0c/research/figures/dual_channel/case_36.png), [39](D:/pip_env/temp/dual_channel_replay_asodag0c/research/figures/dual_channel/case_39.png), [0](D:/pip_env/temp/dual_channel_replay_asodag0c/research/figures/dual_channel/case_0.png), [6](D:/pip_env/temp/dual_channel_replay_asodag0c/research/figures/dual_channel/case_6.png), [12](D:/pip_env/temp/dual_channel_replay_asodag0c/research/figures/dual_channel/case_12.png).

The JSON contains the **complete requested tables**: every group and segment, all four error statistics, sustained-error summaries, A_large per-case results, and per-case \(g\) fits.

Selected parameters:

| Arm | \(q_0\) | \(q_1\) | \(p_{\rm stay}\) | \(q_b\) | \(q_c\) | Dev objective |
|---|---:|---:|---:|---:|---:|---:|
| A1 | 0 | 0.01 | 0.99 | — | — | 0.31873 |
| A2/A3 | 1e−5 | 1e−4 | 0.95 | 1e−6 | 1e−6 | 0.35551 |

All 32 A1 and 320 A2 grid points were evaluated using only the 15 development cases.

Pooled RMS log-stiffness error over all scored ticks:

| Group | A0 | A1 | A2 | A3 |
|---|---:|---:|---:|---:|
| A | .334 | .271 | .273 | .234 |
| A_slow_thin | .153 | .094 | .146 | .146 |
| B | .385 | .373 | .376 | .392 |
| A+B | .165 | .126 | .114 | .109 |
| A_large | .292 | .239 | .304 | .401 |
| cough—all cases | .256 | .201 | .207 | .202 |
| combo—all cases | .813 | .318 | .645 | .439 |

The first five rows use odd-index test cases. Cough and combo were never used for selection.

A_large per-case RMS over all scored ticks:

| Case | Split | A0 | A1 | A2 | A3 |
|---|---|---:|---:|---:|---:|
| 36 | Dev | .257 | .324 | .284 | .258 |
| 37 | Test | .424 | .239 | .341 | .409 |
| 38 | Dev | .122 | .181 | .290 | .429 |
| 39 | Test | .227 | .299 | .297 | .491 |
| 40 | Dev | .227 | .247 | .222 | .217 |
| 41 | Test | .155 | .160 | .269 | .272 |

- **P1—not met.** For cases 36, 37, 39, 40, fast-segment A2/A1 RMS ratios were **.762, 1.175, .913, .779**, respectively; none reached .5. A2 sustained under-estimates were **.711, .360, .487, .503**, all above .12.
- **P2—not met.** Plateau RMS increases were **.0152 in A, .0334 in A_slow_thin, .1434 in A_large**, exceeding .01. B and A+B passed.
- **P3—not met**, evaluated using all-scored-tick RMS. A2−A1 differences were **+.0022 in A** (passed), **+.0526 in A_slow_thin**, **+.1089 in case 38**, and **+.1090 in case 41**.
- **P4—not met.** Final learned \(g\) was within factor 1.5 in **5/15 test cases**, and **11/42 overall**. Case 39 learned **6.517**, versus \(g_{\rm fit}=.974\): factor **6.69**.
- **P5—met.** Event RMS increased from **.3174 to .3873 in B** (**22.0%**) and **.0819 to .1407 in A+B** (**71.8%**).

The three checks:

1. **Synthetic check ran before recordings were opened**, with fixed seed and fixed grid-valued parameters:

   | Arm | Ramp RMS | Flat RMS | Final learned \(g\) |
   |---|---:|---:|---:|
   | A1 | .2659 | .1808 | — |
   | A2 | .8293 | .8628 | 6.365; truth 4 |

   A2 performed poorly here too. Its analytic Jacobian matched finite differences within \(6.1\times10^{-9}\).

2. **A0 reproduced every finite raw value exactly across all 42 cases.**
3. **Truth isolation:** A0 receives only `raw`; A1/A2 receive only arrays derived from copied `t`, `y`, `F_sensor`, `raw`, and `Y`. Their functions receive no log, metadata, or truth arrays; static inspection confirmed no truth references.

Other concerns: background removal did not establish a ceiling—A3 was worse on A_large. Learned \(g\) ranged outside the expected 1–10 range in combo (.525–18.912). No tuning beyond the grids was performed.

All recordings passed unchanged-hash checks. No existing project file, research document, or baseline file was changed; no plant simulation, closed loop, commit, or push was run. **Installing the staged artifacts at the requested repository paths remains blocked by write permissions.**