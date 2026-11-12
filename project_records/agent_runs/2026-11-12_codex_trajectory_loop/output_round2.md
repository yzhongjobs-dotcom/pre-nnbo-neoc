**Q1–Q4 met. Case 7 remains a substantial failure:** delag L = **2.475**, peak force **31.56 N**. Delag also trails traj on A_large mean L: **.189 versus .170**. No constants were tuned.

30 episodes completed with 8 workers in **76.5 s**. The sandbox refused repository writes, so outputs are in this temporary mirror:

[Extended runner](D:/pip_env/temp/trajectory_loop_round2_zct2hvku/experiments/rossi_failure_study/trajectory_loop.py) · [Combined JSON](D:/pip_env/temp/trajectory_loop_round2_zct2hvku/experiments/rossi_failure_study/results/trajectory_loop.json) · [NPZ traces](D:/pip_env/temp/trajectory_loop_round2_zct2hvku/experiments/rossi_failure_study/results/trajectory_loop)

H = huber; T = traj; D = delag; S = A_slow_thin; G = A_large. IMM and Oracle are stored references. The last column gives delag’s **median/final β for t > 3 s**, using the gain applied before each tick’s update.

| Case | Group | L H | L T | L D | IMM | Oracle | Peak H/D, N | β med/final |
|---:|:---|---:|---:|---:|---:|---:|:---|:---|
| 0 | A | .056 | .050 | .056 | .046 | .052 | 8.53/8.52 | 0/0 |
| 1 | A | .226 | 1.732 | .229 | .210 | .339 | 9.16/9.19 | 0/0 |
| 2 | A | .169 | .168 | .170 | .160 | .172 | 8.80/8.80 | 0/0 |
| 3 | A | .031 | .030 | .031 | .030 | .031 | 8.22/8.22 | 0/0 |
| 4 | A | .064 | .060 | .064 | .058 | .076 | 8.47/8.47 | 0/0 |
| 5 | A | .092 | .088 | .092 | .084 | .086 | 8.94/8.95 | 0/0 |
| 6 | S | .212 | .202 | .213 | .195 | .196 | 9.16/9.16 | 0/0 |
| 7 | S | 2.881 | .226 | 2.475 | .218 | .229 | 36.83/31.56 | 0/0 |
| 8 | S | .196 | .225 | .193 | .250 | .181 | 9.02/9.01 | .032/0 |
| 9 | S | .119 | .125 | .119 | .119 | .143 | 8.69/8.71 | 0/0 |
| 10 | S | .035 | .035 | .035 | .034 | .039 | 8.34/8.34 | 0/0 |
| 11 | S | .104 | .393 | .104 | .101 | .109 | 8.64/8.64 | 0/0 |
| 12 | B | .095 | .107 | .096 | .093 | .087 | 8.94/8.94 | 0/0 |
| 13 | B | .038 | .045 | .039 | .032 | 2.362 | 8.63/8.63 | 0/0 |
| 14 | B | .037 | 3.493 | .037 | .030 | .046 | 8.72/8.71 | 0/0 |
| 15 | B | .080 | 2.533 | .079 | .071 | .052 | 8.65/8.65 | 0/0 |
| 16 | B | .022 | .028 | .022 | .026 | 1.153 | 8.32/8.32 | 0/0 |
| 17 | B | .026 | .028 | .026 | .026 | .202 | 8.41/8.41 | 0/0 |
| 18 | A+B | .127 | .130 | .127 | .126 | .134 | 8.74/8.73 | 0/0 |
| 19 | A+B | .050 | .088 | .049 | .052 | .082 | 8.48/8.49 | 0/0 |
| 20 | A+B | .104 | 1.633 | .103 | .109 | .096 | 9.20/9.20 | 0/0 |
| 21 | A+B | .081 | .080 | .081 | .081 | .083 | 8.48/8.49 | 0/0 |
| 22 | A+B | .052 | .053 | .052 | .049 | .057 | 8.66/8.66 | 0/0 |
| 23 | A+B | .050 | .046 | .048 | .045 | .044 | 8.72/8.69 | 0/0 |
| 36 | G | .790 | .476 | .497 | 2.827 | .458 | 12.98/9.98 | .104/.105 |
| 37 | G | .316 | .195 | .250 | .400 | .218 | 9.93/9.66 | .132/.114 |
| 38 | G | .106 | .104 | .106 | .122 | .124 | 8.65/8.65 | 0/.091 |
| 39 | G | .082 | .053 | .076 | .058 | .048 | 8.94/8.86 | .529/.403 |
| 40 | G | .142 | .106 | .122 | .181 | .070 | 9.31/8.93 | .343/.265 |
| 41 | G | .085 | .085 | .085 | .085 | .073 | 8.90/8.88 | 0/0 |

Group L cells are **mean/median**; force RMS cells are **crossings/plateau**, in N. Crossings combine rise, fall and reversal outside events; — means no qualifying ticks.

| Group | L H | L T | L D | Force RMS H | Force RMS T | Force RMS D |
|:---|:---|:---|:---|:---|:---|:---|
| A | .106/.078 | .355/.074 | .107/.078 | .245/.150 | .336/.432 | .247/.150 |
| S | .591/.158 | .201/.214 | .523/.156 | 4.421/.272 | .256/.207 | 3.413/.307 |
| B | .050/.038 | 1.039/.076 | .050/.038 | —/.096 | —/1.218 | —/.096 |
| A+B | .077/.067 | .338/.084 | .077/.066 | —/.116 | —/.116 | —/.116 |
| G | .253/.124 | .170/.105 | .189/.114 | .391/.194 | .286/.181 | .314/.180 |

Identification cells are **median |e| / mean e**, with \(e=\ln(k_{\rm ret}/k_{\rm eff})\). “Under” is median/worst sustained underestimate using round-1’s contiguous 13-tick windows. G-fast = 36, 37, 39, 40; G-slow = 38, 41.

| Group/arm | Plateau | Slow | Rise | Fall | Reversal | Event | Under |
|:---|:---|:---|:---|:---|:---|:---|:---|
| A/H | .10/+.04 | .12/+.12 | .13/+.14 | .24/+.29 | .21/+.21 | — | .41/.69 |
| A/T | .11/+.10 | .14/+.24 | .16/+.20 | .27/+.69 | .21/+.72 | — | .27/.39 |
| A/D | .10/+.04 | .12/+.12 | .14/+.15 | .23/+.25 | .21/+.20 | — | .41/.69 |
| S/H | .09/+.03 | .10/+.06 | .12/−.04 | .23/−.12 | .10/−.08 | — | .38/.96 |
| S/T | .10/+.11 | .12/+.26 | .15/+.01 | .22/+.21 | .25/+.24 | — | .27/.43 |
| S/D | .09/+.03 | .10/+.06 | .14/+.02 | .22/−.09 | .12/−.07 | — | .38/1.00 |
| B/H | .28/+.20 | — | — | — | — | .25/+.18 | .25/.48 |
| B/T | .41/+.95 | .50/+.34 | — | — | — | .47/+1.27 | .03/.49 |
| B/D | .28/+.20 | — | — | — | — | .25/+.18 | .25/.48 |
| A+B/H | .10/+.07 | — | — | — | — | .07/+.05 | .24/.42 |
| A+B/T | .11/+.14 | — | — | — | — | .11/+.20 | .18/.26 |
| A+B/D | .10/+.07 | — | — | — | — | .07/+.05 | .24/.42 |
| G-fast/H | .10/+.05 | .12/+.03 | .14/−.12 | .22/+.22 | .18/+.12 | — | .36/.54 |
| G-fast/T | .16/+.15 | .11/+.09 | .13/+.09 | .10/+.05 | .13/+.11 | — | .24/.36 |
| G-fast/D | .13/+.08 | .10/+.04 | .09/−.04 | .13/+.11 | .13/+.06 | — | .30/.56 |
| G-slow/H | .09/+.03 | .09/+.00 | .12/−.10 | .11/+.08 | — | — | .27/.33 |
| G-slow/T | .08/+.01 | .09/+.04 | .06/+.03 | .03/+.03 | — | — | .25/.33 |
| G-slow/D | .09/+.03 | .09/+.00 | .12/−.10 | .11/+.08 | — | — | .27/.33 |

| Prediction | Verdict | Numbers |
|:---|:---|:---|
| Q1 | **Met** | G mean L **.189388 ≤ .22**; largest case ΔL versus H **+.000484**, case 38, below +.02. |
| Q2 | **Met, using stated exception** | Largest case ΔL **+.002644**, case 1, below +.05. Group mean ΔL: A **+.000730**, S **−.068016**, B **−.000093**, A+B **−.000400**. S exceeds ±.02 through improvement: case 7 ΔL **−.406213**. |
| Q3 | **Met** | G-fast pooled rise mean **−.041008**, within ±.06; fall mean **+.105898 < +.12**. |
| Q4 | **Met** | Every G-fast case had β > .05 on **875/875 ticks (100%)** after 3 s. Median β was **0 in every B case**. |

All delag scored returns were finite. Stored H/S retains eight missing rise/reversal returns; stored T/B retains eight missing event returns. All **65 stored rows remained identical**, and their **65 traces remained byte-identical**. Repository files, research documents and the baseline snapshot were untouched; no commit or push was performed.