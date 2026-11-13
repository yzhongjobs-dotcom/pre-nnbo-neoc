**Finite transported memory did not retain the large-rib benefit or remove the fresh V4 collapse.** It passed M3 and M4, but failed M1, M2, M5 and M6. The H = 0 trace comparison also failed the requested tolerance.

The sandbox refused to modify the existing runner, so the updated [runner](D:/pip_env/temp/finite_memory_round2/experiments/rossi_failure_study/memory_transport_loop.py), [results JSON](D:/pip_env/temp/finite_memory_round2/experiments/rossi_failure_study/results/memory_transport_loop.json), and [new traces](D:/pip_env/temp/finite_memory_round2/experiments/rossi_failure_study/results/memory_transport/) are in the temporary mirror. All 126 requested new episodes completed with eight workers; the three comparison arms reuse their stored rows and pooled summaries.

All tables use the column order **huber, delag2, kf, tmem24, tmem48, mem24**. Values are rounded to three decimals; “—” means the group had no eligible ticks.

### 1. Group L: mean / median

| Group | huber | delag2 | kf | tmem24 | tmem48 | mem24 |
|---|---|---|---|---|---|---|
| A | .106/.078 | .106/.078 | .097/.070 | .103/.077 | .099/.075 | .103/.077 |
| A_slow_thin | .591/.158 | .591/.156 | .146/.155 | .147/.161 | .155/.159 | .155/.161 |
| B | .050/.038 | .050/.038 | .049/.032 | .048/.035 | .046/.032 | .048/.035 |
| A+B | .077/.067 | .077/.067 | .076/.063 | .076/.069 | .075/.065 | .076/.069 |
| A_large original | .253/.124 | .191/.115 | .306/.107 | .239/.111 | .242/.199 | .841/.235 |
| A_large fresh | .439/.199 | .261/.192 | .222/.175 | .884/.276 | .888/.314 | 1.224/.774 |

### 2. Per-case L

| Case | huber | delag2 | kf | tmem24 | tmem48 | mem24 |
|---|---:|---:|---:|---:|---:|---:|
| 6 | .212 | .212 | .192 | .199 | .194 | .199 |
| 7 | 2.881 | 2.882 | .211 | .212 | .220 | .212 |
| 8 | .196 | .193 | .225 | .211 | .255 | .262 |
| 9 | .119 | .119 | .119 | .122 | .123 | .122 |
| 10 | .035 | .035 | .033 | .035 | .034 | .035 |
| 11 | .104 | .104 | .097 | .101 | .101 | .101 |
| 36 | .790 | .494 | .494 | .572 | .507 | 3.675 |
| 37 | .316 | .260 | .994 | .471 | .260 | .697 |
| 38 | .106 | .106 | .120 | .117 | .139 | .117 |
| 39 | .082 | .075 | .057 | .081 | .356 | .344 |
| 40 | .142 | .124 | .075 | .105 | .092 | .126 |
| 41 | .085 | .085 | .094 | .087 | .096 | .087 |
| N1 | .711 | .484 | .478 | .662 | .519 | 1.501 |
| N3 | .081 | .081 | .077 | .090 | .084 | .086 |
| N4 | .209 | .196 | .269 | 4.192 | 4.904 | 4.014 |
| N5 | 1.551 | .355 | .335 | .822 | .817 | 2.210 |
| N7 | .189 | .189 | .178 | .198 | .191 | .198 |
| N8 | .435 | .281 | .415 | .354 | .437 | 1.280 |
| V1 | .130 | .102 | .100 | .138 | .116 | .269 |
| V3 | .048 | .048 | .052 | .056 | .054 | .056 |
| V4 | .915 | .915 | .172 | 3.360 | 2.505 | 3.375 |
| V5 | .829 | .307 | .402 | .523 | .807 | 1.435 |
| V7 | .116 | .116 | .106 | .124 | .115 | .124 |
| V8 | .060 | .058 | .085 | .086 | .105 | .141 |

### 3. Identification: median |e| / mean e

| Group; segment | huber | delag2 | kf | tmem24 | tmem48 | mem24 |
|---|---|---|---|---|---|---|
| A; plateau | .099/.041 | .098/.041 | .071/.079 | .084/.054 | .076/.066 | .084/.054 |
| A; rise | .130/.144 | .130/.150 | .208/.220 | .171/.161 | .155/.169 | .171/.161 |
| A; fall | .237/.286 | .237/.288 | .362/.373 | .318/.341 | .353/.407 | .318/.341 |
| A; reversal | .209/.212 | .209/.214 | .184/.248 | .182/.206 | .207/.297 | .182/.206 |
| Slow thin; plateau | .092/.033 | .091/.032 | .065/.056 | .067/.048 | .062/.054 | .066/.047 |
| Slow thin; rise | .119/−.043 | .119/−.038 | .111/−.076 | .106/−.065 | .108/−.085 | .126/−.109 |
| Slow thin; fall | .235/−.124 | .235/−.124 | .083/.100 | .222/.166 | .060/.088 | .222/.166 |
| Slow thin; reversal | .103/−.078 | .102/−.077 | — | — | — | — |
| B; plateau | .275/.202 | .275/.202 | .296/.240 | .287/.213 | .291/.227 | .287/.213 |
| A+B; plateau | .101/.066 | .101/.066 | .090/.088 | .087/.076 | .083/.083 | .087/.076 |
| Large original; plateau | .096/.041 | .104/.058 | .096/.102 | .097/.069 | .094/.068 | .082/.028 |
| Large original; rise | .139/−.116 | .095/−.038 | .154/.099 | .117/−.039 | .135/−.011 | .173/−.102 |
| Large original; fall | .213/.220 | .128/.102 | .252/.257 | .122/.105 | .186/.188 | .336/.353 |
| Large original; reversal | .177/.122 | .117/.057 | .235/.293 | .142/.100 | .173/.168 | .273/.252 |
| Large fresh; plateau | .123/.093 | .127/.098 | .129/.113 | .116/.090 | .110/.088 | .105/.078 |
| Large fresh; rise | .156/−.105 | .148/−.021 | .183/.109 | .162/−.045 | .175/−.014 | .210/−.080 |
| Large fresh; fall | .319/.289 | .188/.155 | .329/.368 | .272/.271 | .319/.344 | .465/.495 |
| Large fresh; reversal | .194/.129 | .174/.111 | .334/.376 | .252/.213 | .290/.249 | .308/.257 |

B and A+B had no eligible rise, fall or reversal ticks.

**Sustained under-estimate: median / worst**

| Group | huber | delag2 | kf | tmem24 | tmem48 | mem24 |
|---|---|---|---|---|---|---|
| A | .408/.694 | .408/.661 | .115/.207 | .231/.366 | .155/.285 | .231/.366 |
| A_slow_thin | .380/.955 | .380/.955 | .091/.177 | .212/.316 | .120/.181 | .212/.316 |
| B | .251/.479 | .251/.479 | .000/.326 | .108/.457 | .000/.425 | .108/.457 |
| A+B | .237/.422 | .237/.422 | .018/.132 | .136/.207 | .063/.134 | .136/.207 |
| Large original | .345/.539 | .296/.560 | .231/.353 | .302/.593 | .253/.271 | .292/.548 |
| Large fresh | .377/1.526 | .317/1.526 | .135/.425 | .266/1.576 | .252/1.934 | .322/1.960 |

### 4. Force: mean ripple RMS / pooled crossing-force RMS

| Group | huber | delag2 | kf | tmem24 | tmem48 | mem24 |
|---|---|---|---|---|---|---|
| A | .080/.245 | .080/.245 | .076/.223 | .080/.224 | .078/.234 | .080/.224 |
| A_slow_thin | .078/4.421 | .078/4.421 | .090/.369 | .090/.383 | .093/.382 | .092/.465 |
| B | .071/— | .071/— | .067/— | .071/— | .070/— | .071/— |
| A+B | .083/— | .083/— | .081/— | .083/— | .083/— | .083/— |
| Large original | .082/.391 | .083/.315 | .069/.451 | .084/.365 | .079/.315 | .089/.596 |
| Large fresh | .077/.604 | .076/.380 | .067/.306 | .086/1.358 | .107/2.305 | .089/1.225 |

### 5. Prespecified predictions

- **M1 — not met.** `tmem24` mean L was .239 original and .884 fresh, versus limits .195 and .265. Its largest deterioration versus `delag2` was **+3.996 in N4**, above the .05 per-scene limit.
- **M2 — not met.** Case 7 was .212, but fresh V4 was **3.360**, above .5.
- **M3 — met.** Across original cases 0–23 except 7, the largest absolute `tmem24 − huber` difference was .0158.
- **M4 — met.** Combined large-rib mean L for `mem24` exceeded `tmem24` by .428, above .03. This average includes severe failures.
- **M5 — not met.** Worst accumulated time above the ratio threshold: `tmem24` V4 **.216 s** (peak ratio 8.187); `tmem48` V4 **.256 s** (3.583); `mem24` N5 **.520 s** (4.027). The last exceeds .5 s.
- **M6 — not met.** Pooled original plateau median |e| fell **12.7%** for `tmem24` and **15.3%** for `tmem48`, below the respective 15% and 25% targets. Fast original large-rib fall mean e was **+.104**, meeting its +.15 limit.

**H = 0 check — not met.** On case 0, the returned stiffness trace had matching NaN positions but a maximum absolute difference of **3.59×10⁻⁷** from the stored `delag2` trace, exceeding 10⁻¹². This comparison checked the stored stiffness trace; it does not establish equality of the full returned triple.