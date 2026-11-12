**A_large force control improved, but normal-scene preservation failed. P2 and P4 met; P1, P3 and P5 did not.** No tuning or extra arms were run.

65 episodes completed with 8 subprocess workers in 142.7 s. Sandbox restrictions required a temporary mirror:

[Runner](D:/pip_env/temp/trajectory_loop_wyxvbgas/experiments/rossi_failure_study/trajectory_loop.py) · [JSON—all per-run metrics and summaries](D:/pip_env/temp/trajectory_loop_wyxvbgas/experiments/rossi_failure_study/results/trajectory_loop.json) · [65 NPZ traces](D:/pip_env/temp/trajectory_loop_wyxvbgas/experiments/rossi_failure_study/results/trajectory_loop)

H = huber; T = traj; P = traj_priv; S = A_slow_thin; G = A_large. IMM and Oracle are stored references.

| Case | Group | L H | L T | L P | IMM | Oracle | Peak H/T, N |
|---:|:---|---:|---:|---:|---:|---:|:---|
| 0 | A | .056 | .050 | — | .046 | .052 | 8.53/8.31 |
| 1 | A | .226 | 1.732 | — | .210 | .339 | 9.16/15.69 |
| 2 | A | .169 | .168 | — | .160 | .172 | 8.80/8.78 |
| 3 | A | .031 | .030 | — | .030 | .031 | 8.22/8.19 |
| 4 | A | .064 | .060 | — | .058 | .076 | 8.47/8.49 |
| 5 | A | .092 | .088 | — | .084 | .086 | 8.94/8.58 |
| 6 | S | .212 | .202 | — | .195 | .196 | 9.15/9.21 |
| 7 | S | 2.881 | .226 | — | .218 | .229 | 36.83/9.20 |
| 8 | S | .196 | .225 | — | .250 | .181 | 9.02/9.09 |
| 9 | S | .119 | .125 | — | .119 | .143 | 8.69/8.66 |
| 10 | S | .035 | .035 | — | .034 | .039 | 8.33/8.30 |
| 11 | S | .104 | .393 | — | .101 | .109 | 8.64/9.07 |
| 12 | B | .095 | .107 | — | .093 | .087 | 8.94/8.96 |
| 13 | B | .038 | .045 | — | .032 | 2.362 | 8.63/8.66 |
| 14 | B | .037 | 3.493 | — | .030 | .046 | 8.72/11.55 |
| 15 | B | .080 | 2.533 | — | .071 | .052 | 8.65/9.73 |
| 16 | B | .022 | .028 | — | .026 | 1.153 | 8.32/8.36 |
| 17 | B | .026 | .028 | — | .026 | .202 | 8.41/8.44 |
| 18 | A+B | .127 | .130 | — | .126 | .134 | 8.74/8.74 |
| 19 | A+B | .050 | .088 | — | .052 | .082 | 8.48/8.99 |
| 20 | A+B | .104 | 1.633 | — | .109 | .096 | 9.20/18.23 |
| 21 | A+B | .081 | .080 | — | .081 | .083 | 8.48/8.48 |
| 22 | A+B | .052 | .053 | — | .049 | .057 | 8.66/8.84 |
| 23 | A+B | .050 | .046 | — | .045 | .044 | 8.72/8.63 |
| 36 | G | .790 | .476 | .491 | 2.827 | .458 | 12.97/9.29 |
| 37 | G | .316 | .195 | .245 | .400 | .218 | 9.93/9.14 |
| 38 | G | .106 | .104 | .110 | .122 | .124 | 8.65/8.63 |
| 39 | G | .082 | .053 | .075 | .058 | .048 | 8.94/8.51 |
| 40 | G | .142 | .106 | .124 | .181 | .070 | 9.31/9.01 |
| 41 | G | .085 | .085 | — | .085 | .073 | 8.90/9.16 |

Group L cells are **mean/median**; force RMS cells are **crossings/plateau**, in N. Crossings combine rise, fall and reversal outside events; — means no qualifying ticks. P aggregates only cases 36–40.

| Group | L H | L T | L P | Force RMS H | Force RMS T | Force RMS P |
|:---|:---|:---|:---|:---|:---|:---|
| A | .106/.078 | .355/.074 | — | .245/.150 | .336/.432 | — |
| S | .591/.158 | .201/.214 | — | 4.421/.272 | .256/.207 | — |
| B | .050/.038 | 1.039/.076 | — | —/.096 | —/1.218 | — |
| A+B | .077/.067 | .338/.084 | — | —/.116 | —/.116 | — |
| G | .253/.124 | .170/.105 | .209/.124 | .391/.194 | .286/.181 | .311/.186 |

Identification cells are **median |e| / mean e**. “Under” is median/worst sustained underestimate. G-fast = 36, 37, 39, 40; G-slow = 38, 41; G-slow/P contains only 38. Statistics pool finite returns: eight H/S rise/reversal returns and eight T/B event returns were undefined; missing counts are retained in JSON.

| Group/arm | Plateau | Slow | Rise | Fall | Reversal | Event | Under |
|:---|:---|:---|:---|:---|:---|:---|:---|
| A/H | .10/+.04 | .12/+.12 | .13/+.14 | .24/+.29 | .21/+.21 | — | .41/.69 |
| A/T | .11/+.10 | .14/+.24 | .16/+.20 | .27/+.69 | .21/+.72 | — | .27/.39 |
| S/H | .09/+.03 | .10/+.06 | .12/−.04 | .23/−.12 | .10/−.08 | — | .38/.96 |
| S/T | .10/+.11 | .12/+.26 | .15/+.01 | .22/+.21 | .25/+.24 | — | .27/.43 |
| B/H | .28/+.20 | — | — | — | — | .25/+.18 | .25/.48 |
| B/T | .41/+.95 | .50/+.34 | — | — | — | .47/+1.27 | .03/.49 |
| A+B/H | .10/+.07 | — | — | — | — | .07/+.05 | .24/.42 |
| A+B/T | .11/+.14 | — | — | — | — | .11/+.20 | .18/.26 |
| G-fast/H | .10/+.05 | .12/+.03 | .14/−.12 | .22/+.22 | .18/+.12 | — | .36/.54 |
| G-fast/T | .16/+.15 | .11/+.09 | .13/+.09 | .10/+.05 | .13/+.11 | — | .24/.36 |
| G-fast/P | .13/+.07 | .10/+.04 | .10/−.03 | .12/+.08 | .13/+.07 | — | .30/.57 |
| G-slow/H | .09/+.03 | .09/+.00 | .12/−.10 | .11/+.08 | — | — | .27/.33 |
| G-slow/T | .08/+.01 | .09/+.04 | .06/+.03 | .03/+.03 | — | — | .25/.33 |
| G-slow/P | .05/+.02 | .07/+.05 | .20/+.20 | .15/+.16 | — | — | .24/.24 |

| Prediction/check | Verdict | Numbers |
|:---|:---|:---|
| P1 | **Not met** | Rise mean −.116 → +.092 misses ±.06. Fall +.224 → +.047 passes; median absolute fall/reversal errors drop 53.9%/26.1%, both passing. |
| P2 | **Met** | G mean L .25346 → .16989 ≤ .22. All six improve; largest case ΔL = −.000485. |
| P3 | **Not met** | Mean ΔL: A **+.24817**, S **−.39018**; both exceed ±.02. |
| P4 | **Met** | Mean ΔL: B **+.98897**, A+B **+.26114**. |
| P5 | **Not met** | P has higher L in all five pairs: ΔL +.00598 to +.04957. Paired mean T/P = .18695/.20881. |
| Harness | **Passed** | All 30 H cases: max \|ΔL\| = **0**. |
| Equivalence | **Passed** | Max absolute differences: case 0 **5.33e−15**, 36 **1.26e−13**, 39 **7.11e−15**; finite patterns identical. |

All requested metrics and trace fields were verified. Existing repository files and the baseline snapshot were untouched; no commit or push was performed.