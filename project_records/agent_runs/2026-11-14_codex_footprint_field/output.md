Completed the offline replay from the repository root with `python -X utf8`. The only implementation correction was making input paths relative to the script, so the documented root-level command works. The saved traces have time at compact column 0, scan position at 1, force at 3, and `k_eff` at 9 (source log column 26). All 52,369 nonzero position steps are positive, so the field filter’s high-edge-only cell-entry assumption holds for these traces.

The saved [JSON evidence](</E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/footprint_field_filter.json>) contains per-case results for all 55 method/grid settings, plus held-out and all-case summaries for all seven groups. The [experiment note](</E:/github_project/scan-impedance-update/research/footprint_field_replay_20261114.md>) records the results and verdicts. Global settings selected on the even cases were time RW `q=0.0003`, position RW `q=3`, position LLT `qs=10000`, and field `qe=10, c=0.35`.

Each vector below lists RMS for **24-tick window / time RW / position RW / position LLT / field**.

| Group | Odd held-out cases | All cases |
|---|---|---|
| A | .169 / .100 / .103 / .122 / .119 | .173 / .106 / .110 / .127 / .126 |
| A_slow_thin | .165 / .101 / .099 / .118 / .102 | .167 / .109 / .112 / .134 / .108 |
| B | .169 / .089 / .093 / .104 / .110 | .169 / .089 / .094 / .102 / .108 |
| A+B | .172 / .090 / .093 / .098 / .102 | .169 / .087 / .089 / .093 / .097 |
| A_large | .198 / .205 / .202 / .257 / .183 | .208 / .221 / .217 / .272 / .195 |
| combo | .222 / .211 / .211 / .247 / .213 | .247 / .225 / .223 / .253 / .229 |
| cough | .170 / .107 / .111 / .122 / .121 | .169 / .109 / .111 / .119 / .120 |

**Prediction verdicts:**

- **P1: Not supported.** Held-out field RMS for A, A_slow_thin, B and A+B was `.119, .102, .110, .102`; all-case RMS was `.126, .108, .108, .097`. Only all-case A+B met `.10`.
- **P2: Partly supported.** On A_large, field RMS was 7.6% below window on held-out cases (`.183` vs `.198`) and 6.3% below across all cases (`.195` vs `.208`). Rise bias did not improve by a third: its absolute value increased from `.0747` to `.0766` on held-out cases, and from `.0802` to `.0857` across all cases.
- **P3: Field claim not supported; time-memory claim depends on metric.** Field RMS penalty over the retrospective group-best exceeded 15% in B and A+B: `54.5%` and `42.4%` held out, `43.7%` and `41.0%` all-case. Time RW’s largest RMS penalty over its group-best was `24.2%` held out and `23.4%` all-case, below 30%. By MSE, its penalty exceeded 30% in B and A+B held out (`39.7%`, `54.3%`) and in B and A+B all-case (`34.6%`, `52.2%`). Thus the time-memory prediction is supported only under the MSE reading.
- **P4: Not supported overall.** Five-group mean MSE/window was `.475` for position RW and `.470` for time RW held out; across all cases it was `.506` and `.494`. Position RW was slightly better in held-out A_slow_thin and A_large, and all-case A_large; pooled normalized ratios also favored time RW (`.512` vs `.514` held out; `.553` vs `.560` all-case).

This is an offline comparison of simulated `k_eff` truth with synthetic white observation noise. It is not a fresh closed-loop evaluation; group-best settings are retrospective, and the field-direction check covers these monotonically increasing traces only. The JSON structure was checked: 42 cases, 55 settings per case, and all seven groups in both summaries. No baseline files were edited, and nothing was committed or pushed.