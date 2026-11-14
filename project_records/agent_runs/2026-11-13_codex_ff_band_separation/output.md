Completed **236 specified new runs + 3 reproduction runs**, reusing 96 stored rows. **B1–B5 all failed as complete predictions**, although averaging rescued case 7 and V4.

Repository writes were refused; files are in the temporary mirror: [script](D:/pip_env/temp/ff_band_separation_loop_20261113/experiments/rossi_failure_study/ff_band_separation_loop.py), [results JSON](D:/pip_env/temp/ff_band_separation_loop_20261113/experiments/rossi_failure_study/results/ff_band_separation_loop.json). The JSON contains 335 unique rows, without traces. Existing files were unchanged; no commit or push.

**Block S.** Each cell is `(1,0) → (6,0)`, expressed as **L / peak force (N)**.

| f | Case 0 | Case 1 | Case 3 | Case 36 |
|---|---|---|---|---|
| .5 | 5.608/38.0 → 5.998/28.7 | 5.640/74.1 → 6.431/39.6 | 2.272/11.7 → 3.728/16.6 | 5.630/50.6 → 6.382/45.5 |
| .7 | 4.430/17.1 → .100/8.4 | 5.309/35.6 → .241/8.9 | 1.960/10.9 → .071/8.3 | 5.405/45.4 → 1.045/10.7 |
| .85 | .073/8.3 → .067/8.3 | 4.962/26.2 → .164/8.8 | .046/8.2 → .043/8.2 | .650/10.1 → 1.214/11.2 |
| 1 | .052/8.3 → .053/8.3 | .339/9.3 → .139/8.7 | .031/8.2 → .030/8.2 | .458/9.2 → 1.401/11.6 |
| 1.2 | .040/8.3 → .046/8.3 | .221/9.2 → .129/8.6 | .021/8.2 → .021/8.1 | .461/9.4 → 1.653/12.2 |
| 1.5 | .038/8.4 → .045/8.3 | .196/9.1 → .128/8.6 | .015/8.1 → .014/8.1 | .531/9.9 → 2.017/13.1 |
| 2 | .060/8.5 → .055/8.3 | .654/9.9 → .131/8.5 | 3.575/13.9 → .011/8.1 | .756/10.6 → 2.568/14.9 |

**Block M.** Columns throughout: **H** = stored huber; **H0/H1/H2.5** = huber `(6,extra)`; **D** = stored delag2; **D0** = delag2 `(6,0)`. Group cells give **mean L / median L / mean rmse (N)**.

| Group | H | H0 | H1 | H2.5 | D | D0 |
|---|---|---|---|---|---|---|
| A | .106/.078/.160 | .409/.069/.328 | .462/.075/.368 | .602/.095/.483 | .106/.078/.160 | .452/.069/.361 |
| A_slow_thin | .591/.158/.476 | .132/.127/.182 | .144/.140/.191 | .201/.191/.228 | .591/.156/.476 | .132/.127/.182 |
| B | .050/.038/.109 | .096/.079/.153 | .072/.068/.131 | .116/.094/.170 | .050/.038/.109 | .096/.079/.153 |
| A+B | .077/.067/.139 | .151/.084/.191 | .124/.084/.174 | .158/.127/.196 | .077/.067/.139 | .151/.084/.191 |
| A_large original | .253/.124/.252 | .498/.118/.388 | .403/.115/.330 | .315/.185/.291 | .191/.115/.218 | .443/.112/.352 |
| A_large fresh | .439/.199/.352 | .727/.296/.525 | .732/.220/.604 | .768/.235/.749 | .261/.192/.258 | .481/.273/.373 |
| All 42 | .280/.104/.263 | .392/.103/.328 | .381/.107/.343 | .418/.158/.409 | .220/.103/.231 | .319/.103/.284 |

Per-case **L** below. Fresh **N** = `neural_white_28`, **V** = `viscous_14`; suffix is scene index.

| Case | H | H0 | H1 | H2.5 | D | D0 |
|---|---:|---:|---:|---:|---:|---:|
| 6 | .212 | .191 | .213 | .303 | .212 | .191 |
| 7 | 2.881 | .149 | .163 | .225 | 2.882 | .149 |
| 8 | .196 | .210 | .223 | .338 | .193 | .209 |
| 9 | .119 | .101 | .111 | .136 | .119 | .101 |
| 10 | .035 | .035 | .037 | .043 | .035 | .035 |
| 11 | .104 | .104 | .117 | .157 | .104 | .104 |
| 36 | .790 | 1.909 | 1.561 | .958 | .494 | 1.657 |
| 37 | .316 | .696 | .485 | .319 | .260 | .638 |
| 38 | .106 | .117 | .118 | .152 | .106 | .117 |
| 39 | .082 | .074 | .066 | .164 | .075 | .066 |
| 40 | .142 | .119 | .111 | .205 | .124 | .106 |
| 41 | .085 | .072 | .078 | .090 | .085 | .072 |
| N1 | .711 | 1.743 | 1.286 | .906 | .484 | 1.572 |
| N3 | .081 | .077 | .079 | .098 | .081 | .077 |
| N4 | .209 | .979 | .613 | .303 | .196 | .908 |
| N5 | 1.551 | 1.418 | 1.464 | 1.648 | .355 | 1.246 |
| N7 | .189 | .252 | .214 | .223 | .189 | .247 |
| N8 | .435 | .854 | .601 | .450 | .281 | .830 |
| V1 | .130 | .340 | .225 | .248 | .102 | .299 |
| V3 | .048 | .042 | .045 | .052 | .048 | .042 |
| V4 | .915 | .062 | .055 | .106 | .915 | .062 |
| V5 | .829 | 2.790 | 4.011 | 4.885 | .307 | .327 |
| V7 | .116 | .095 | .103 | .129 | .116 | .095 |
| V8 | .060 | .076 | .095 | .169 | .058 | .074 |

**Late push/pull.** Each cell gives **L / event-window force rms (N)**. Directions use stored Huber labels; scene and direction uniquely identify each run.

| Scene | Toward H | Toward H0 | Away H | Away H0 |
|---|---:|---:|---:|---:|
| N0 | .296/.542 | .746/1.190 | .258/.465 | .767/1.218 |
| N2 | .088/.238 | .239/.560 | .087/.235 | .254/.585 |
| N6 | .379/.662 | 1.534/1.908 | .333/.549 | 1.418/1.627 |
| V0 | .102/.210 | .112/.266 | .101/.203 | .114/.271 |
| V2 | .054/.159 | .060/.191 | .055/.164 | .062/.199 |
| V6 | .126/.394 | .319/.756 | .111/.347 | .282/.686 |
| Mean | .174/.368 | .502/.812 | .157/.327 | .483/.764 |

**Predictions and checks:**

- **B1 failed.** In case order 0/1/3/36, averaged `L(.7)/L(1)` was **1.905/1.728/2.379/.746**: case 3 exceeded 2. At `.5`, averaging worsened L in **all four cases**, rather than improving it threefold. At `f=2`, case 36’s ratio was **1.832**, exceeding 1.5; all `f=1.5` ratios passed.
- **B2 failed overall:** case 7 **.149 < .5** and V4 **.062 < .5** passed; N5 **1.418 > .6** failed.
- **B3 failed overall:** H0 mean-L changes were A **+.303**, slow-thin excluding 7 **−.005**, B **+.046**, A+B **+.073**. A and A+B exceeded their .03 bounds.
- **B4 failed:** D0 mean L was **.443** original and **.481** fresh, exceeding **.19/.26**.
- **B5 failed:** case 7 with `(6,2.5)` had **L=.225**, below the predicted collapse threshold of 1.
- **Both reproduction checks passed exactly:** Huber cases **0,7,36**, max `|dL|=0`; oracle cases **0,1,3,36**, max `|dL|=0`.

The sweep supports substantial tolerance improvement near `f=.7`, but the predicted wider robustness and acceptable Method-level costs did **not** hold across these closed-loop tests. Late-event force rms increased in every paired run.