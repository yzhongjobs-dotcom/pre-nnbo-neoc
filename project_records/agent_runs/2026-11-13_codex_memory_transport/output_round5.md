`kf2` passed D1, D4 and D5, but failed D2, D3 and D6. It reduced case 37’s loss from `kf`’s 0.994 to 0.289, while case 36 and fresh N5 remained problematic. Small output perturbations also caused substantial loss changes in three scenes.

The repository denied the new-file write. The script, 68 run rows and seven unperturbed traces are in the temporary mirror:

- [kf2_loop.py](</D:/pip_env/temp/kf2_round5/experiments/rossi_failure_study/kf2_loop.py>)
- [kf2_loop.json](</D:/pip_env/temp/kf2_round5/experiments/rossi_failure_study/results/kf2_loop.json>)
- [Trace directory](</D:/pip_env/temp/kf2_round5/experiments/rossi_failure_study/results/kf2/>)

Existing files were not changed. Fresh A has no stored `gated` or `kf` rows, so those entries are unavailable.

| Group | huber | delag2 | gated | kf | kf2 |
|---|---:|---:|---:|---:|---:|
| A | .106/.078 | .106/.078 | .129/.122 | .097/.070 | .097/.070 |
| A_slow_thin | .591/.158 | .591/.156 | .164/.173 | .146/.155 | .146/.155 |
| B | .050/.038 | .050/.038 | .086/.100 | .049/.032 | .049/.032 |
| A+B | .077/.067 | .077/.067 | .076/.059 | .076/.063 | .076/.063 |
| A_large original | .253/.124 | .191/.115 | .616/.261 | .306/.107 | .203/.107 |
| A_large fresh | .439/.199 | .261/.192 | .546/.399 | .222/.175 | .245/.175 |
| Fresh A | .068/.057 | .068/.057 | — | — | .056/.043 |
| All 42 | .280/.104 | .220/.103 | .309/.123 | .160/.095 | .151/.095 |

Cells show mean/median L.

| Scene | huber | delag2 | gated | kf | kf2 | kf2 median ε after 3 s |
|---|---:|---:|---:|---:|---:|---:|
| 6 | .212 | .212 | .200 | .192 | .192 | .000 |
| 7 | 2.881 | 2.882 | .223 | .211 | .211 | .000 |
| 8 | .196 | .193 | .194 | .225 | .224 | −.008 |
| 9 | .119 | .119 | .110 | .119 | .119 | .000 |
| 10 | .035 | .035 | .153 | .033 | .033 | .104 |
| 11 | .104 | .104 | .106 | .097 | .097 | .000 |
| 36 | .790 | .494 | 2.444 | .494 | .589 | −.306 |
| 37 | .316 | .260 | .556 | .994 | .289 | −.685 |
| 38 | .106 | .106 | .102 | .120 | .120 | .000 |
| 39 | .082 | .075 | .316 | .057 | .049 | −.347 |
| 40 | .142 | .124 | .206 | .075 | .076 | −.065 |
| 41 | .085 | .085 | .069 | .094 | .094 | .000 |
| N1 | .711 | .484 | .928 | .478 | .510 | −.430 |
| N3 | .081 | .081 | .066 | .077 | .077 | −.062 |
| N4 | .209 | .196 | 1.403 | .269 | .270 | .164 |
| N5 | 1.551 | .355 | .968 | .335 | .590 | −.249 |
| N7 | .189 | .189 | .171 | .178 | .178 | .050 |
| N8 | .435 | .281 | .745 | .415 | .377 | .335 |
| V1 | .130 | .102 | .429 | .100 | .096 | −.414 |
| V3 | .048 | .048 | .045 | .052 | .052 | .000 |
| V4 | .915 | .915 | .369 | .172 | .172 | .000 |
| V5 | .829 | .307 | 1.130 | .402 | .422 | −.044 |
| V7 | .116 | .116 | .101 | .106 | .106 | .000 |
| V8 | .060 | .058 | .192 | .085 | .085 | .000 |

Identification values below are pooled median \(|e|\)/mean \(e\), with \(e=\ln(k_\mathrm{ret}/k_\mathrm{eff})\). Sustained values are the median/worst per-scene sustained underestimate from the existing metric.

| Group, arm | Plateau | Rise | Fall | Reversal | Sustained median/worst |
|---|---:|---:|---:|---:|---:|
| A, delag2 | .098/.041 | .130/.150 | .237/.288 | .209/.214 | .408/.661 |
| A, kf2 | .071/.079 | .208/.220 | .362/.372 | .184/.248 | .115/.207 |
| A_slow_thin, delag2 | .091/.032 | .119/−.038 | .235/−.124 | .102/−.077 | .380/.955 |
| A_slow_thin, kf2 | .065/.056 | .110/−.074 | .083/.100 | — | .091/.184 |
| Original fast A_large, delag2 | .128/.078 | .094/−.037 | .129/.103 | .117/.057 | .296/.560 |
| Original fast A_large, kf2 | .104/.086 | .133/.075 | .283/.286 | .299/.312 | .232/.256 |
| Fresh A_large, delag2 | .127/.098 | .148/−.021 | .188/.155 | .174/.111 | .317/1.526 |
| Fresh A_large, kf2 | .115/.103 | .188/.122 | .405/.447 | .395/.420 | .169/.465 |

| Scene | Unperturbed L | +0.5% | −0.5% |
|---|---:|---:|---:|
| 7 | .211 | .210 | .212 |
| 36 | .589 | .572 | .650 |
| 37 | .289 | .277 | .296 |
| N4 | .270 | .271 | .270 |
| N5 | .590 | .631 | .828 |
| V4 | .172 | .204 | .167 |
| V5 | .422 | .563 | .396 |

- **D1 met:** case 7 = .211; V4 = .172.
- **D2 not met:** original A_large mean = .203, above .20; fresh mean = .245, below .25. The largest excess over `delag2` was N5, +.235, above .05.
- **D3 not met:** case 8 was +.0287 versus `huber`; fresh A9 and A11 were −.0292 and −.0244.
- **D4 met:** the 42-scene mean was .151, below .18.
- **D5 met:** the longest accumulated time above ratio 3 was .336 s, for V5 at +0.5%; the limit was .5 s.
- **D6 not met:** absolute L changes exceeded .05 for case 36 at −0.5% (.060), N5 at −0.5% (.238), and V5 at +0.5% (.141).