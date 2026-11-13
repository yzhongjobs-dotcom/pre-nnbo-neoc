**C1, C2 and C4 failed; C3 and C5 passed.** The hold reduced the collapse scores in case 7 and V4, but neither fell below 0.5. It also caused a severe regression in case 36.

The repository write was denied by the sandbox. Outputs are in the temporary mirror:

- [Runner](D:/pip_env/temp/chi_hold_mirror/experiments/rossi_failure_study/chi_hold_loop.py)
- [Results JSON](D:/pip_env/temp/chi_hold_mirror/experiments/rossi_failure_study/results/chi_hold_loop.json)
- [Requested traces](D:/pip_env/temp/chi_hold_mirror/experiments/rossi_failure_study/results/chi_hold/)

Stored earlier rows were reused. The fresh A Huber reruns matched their stored rows exactly.

Group entries are **mean/median L**. Columns are Huber, delag2, gated, gated_delag2, Huber-χ and delag2-χ, respectively; `—` means no stored row.

| Group | Huber | delag2 | gated | gated_delag2 | Huber-χ | delag2-χ |
|---|---:|---:|---:|---:|---:|---:|
| A | .106/.078 | .106/.078 | .129/.122 | .128/.115 | .106/.078 | .106/.078 |
| A_slow_thin | .591/.158 | .591/.156 | .164/.173 | .164/.172 | .481/.158 | .480/.156 |
| B | .050/.038 | .050/.038 | .086/.100 | .086/.100 | .050/.038 | .050/.038 |
| A+B | .077/.067 | .077/.067 | .076/.059 | .076/.059 | .077/.067 | .077/.067 |
| A_large original | .253/.124 | .191/.115 | .616/.261 | .269/.235 | .588/.178 | .552/.115 |
| A_large fresh | .439/.199 | .261/.192 | .546/.399 | .437/.305 | 1.015/.199 | .238/.192 |
| Fresh A | .068/.057 | .068/.057 | — | — | .067/.057 | .067/.057 |
| All 42 earlier | .280/.104 | .220/.103 | .309/.123 | .228/.123 | .476/.104 | .249/.103 |

Per-case columns follow the same arm order; the final column is delag2-χ held fraction / longest hold in seconds.

| Case | Huber | delag2 | gated | gated_delag2 | Huber-χ | delag2-χ | Hold |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | .212 | .212 | .200 | .200 | .212 | .212 | .001/.008 |
| 7 | 2.881 | 2.882 | .223 | .223 | 2.217 | 2.217 | .107/.656 |
| 8 | .196 | .193 | .194 | .191 | .196 | .193 | .000/.000 |
| 9 | .119 | .119 | .110 | .110 | .119 | .119 | .000/.000 |
| 10 | .035 | .035 | .153 | .153 | .035 | .035 | .000/.000 |
| 11 | .104 | .104 | .106 | .106 | .104 | .104 | .013/.080 |
| 36 | .790 | .494 | 2.444 | .710 | 2.679 | 2.655 | .041/.400 |
| 37 | .316 | .260 | .556 | .263 | .316 | .260 | .000/.000 |
| 38 | .106 | .106 | .102 | .102 | .106 | .106 | .000/.000 |
| 39 | .082 | .075 | .316 | .265 | .082 | .075 | .000/.000 |
| 40 | .142 | .124 | .206 | .208 | .251 | .124 | .021/.160 |
| 41 | .085 | .085 | .069 | .069 | .092 | .092 | .019/.184 |
| N1 | .711 | .484 | .928 | 1.142 | .711 | .484 | .000/.000 |
| N3 | .081 | .081 | .066 | .069 | .081 | .081 | .015/.152 |
| N4 | .209 | .196 | 1.403 | .490 | .209 | .196 | .004/.032 |
| N5 | 1.551 | .355 | .968 | 1.329 | 6.286 | .342 | .002/.024 |
| N7 | .189 | .189 | .171 | .172 | .189 | .189 | .000/.000 |
| N8 | .435 | .281 | .745 | .481 | 2.336 | .281 | .000/.000 |
| V1 | .130 | .102 | .429 | .419 | .130 | .102 | .000/.000 |
| V3 | .048 | .048 | .045 | .045 | .048 | .048 | .000/.000 |
| V4 | .915 | .915 | .369 | .141 | .658 | .656 | .234/.904 |
| V5 | .829 | .307 | 1.130 | .666 | 1.358 | .307 | .000/.000 |
| V7 | .116 | .116 | .101 | .101 | .116 | .116 | .000/.000 |
| V8 | .060 | .058 | .192 | .192 | .060 | .058 | .000/.000 |

Identification entries are pooled **median |e| / mean e**, with \(e=\log(k_{\rm ret}/k_{\rm eff})\). Sustained entries are median/worst under-estimate.

| Scenes, arm | Plateau | Rise | Fall | Reversal | Sustained |
|---|---:|---:|---:|---:|---:|
| Original fast, delag2 | .128/.078 | .094/−.037 | .129/.103 | .117/.057 | .296/.560 |
| Original fast, delag2-χ | .128/.072 | .092/−.040 | .129/.106 | .125/.070 | .300/.560 |
| Fresh large, delag2 | .127/.098 | .148/−.021 | .188/.155 | .174/.111 | .317/1.526 |
| Fresh large, delag2-χ | .127/.100 | .151/−.020 | .186/.161 | .176/.115 | .281/.707 |

- **C1 not met:** case 7 L = 2.217; V4 L = .656.
- **C2 not met:** original large mean = .552 versus the .20 limit; fresh mean excluding V4 = .200 versus .27. Case 36 worsened by 2.161 versus delag2, exceeding .05.
- **C3 met:** largest absolute difference from Huber was .0081, on fresh A8.
- **C4 not met:** 42-scene mean = .249, above delag2’s .220 and gated_delag2’s .228.
- **C5 met:** outside case 7, the largest held fraction was 2.4% and the longest hold was .136 s, both in case 14.