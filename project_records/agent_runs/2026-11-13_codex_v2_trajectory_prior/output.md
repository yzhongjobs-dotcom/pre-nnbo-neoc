Completed **66 `v2t` runs plus one fixed-beta-zero validation**, using 8 workers and 132 stored reference rows. Repository writes were denied, so outputs are in the [temporary mirror](D:/pip_env/temp/v2_trajectory_prior_2809ef305ed5489995ad98f62ab5880d): [runner](D:/pip_env/temp/v2_trajectory_prior_2809ef305ed5489995ad98f62ab5880d/experiments/rossi_failure_study/v2_trajectory_prior_loop.py), [JSON results](D:/pip_env/temp/v2_trajectory_prior_2809ef305ed5489995ad98f62ab5880d/experiments/rossi_failure_study/results/v2_trajectory_prior_loop.json), and [traces](D:/pip_env/temp/v2_trajectory_prior_2809ef305ed5489995ad98f62ab5880d/experiments/rossi_failure_study/results/v2_trajectory_prior). No existing repository files changed; no commit or push.

Throughout, **H/D/V = `huber` / `delag2` / `v2t`**. N/V scene labels use the actual scene index.

Mean / median \(L\):

| Group | H | D | V |
|---|---:|---:|---:|
| A | .106356/.077943 | .105942/.077943 | .105961/.077943 |
| A_slow_thin | .591190/.157503 | .590797/.156138 | .590721/.156113 |
| B | .049903/.037889 | .049903/.037889 | .049903/.037889 |
| A+B | .077288/.066649 | .077288/.066649 | .077288/.066649 |
| A_large original | .253460/.123797 | .190634/.114892 | .190202/.113443 |
| A_large fresh | .439460/.198818 | .261112/.192332 | .255390/.192299 |
| Fresh A | .067998/.057384 | .067998/.057384 | .067998/.057384 |

Per-case A_large \(L\), each cell **H/D/V**:

| Original | Fresh neural | Fresh viscous |
|---|---|---|
| 36: .7903/.4943/.4996 | N1: .7106/.4839/.4916 | V1: .1298/.1015/.1014 |
| 37: .3157/.2598/.2562 | N3: .0808/.0813/.0812 | V3: .0483/.0483/.0483 |
| 38: .1059/.1058/.1058 | N4: .2087/.1959/.1959 | V4: .9155/.9152/.9152 |
| 39: .0821/.0749/.0734 | N5: 1.5510/.3551/.3590 | V5: .8287/.3074/.2283 |
| 40: .1417/.1240/.1211 | N7: .1889/.1887/.1887 | V7: .1165/.1164/.1164 |
| 41: .0851/.0851/.0851 | N8: .4346/.2811/.2801 | V8: .0599/.0584/.0584 |

Identification uses returned stiffness and \(e=\ln(k_{\rm ret}/k_{\rm eff})\), with unchanged masks and validity rules. Cells show **median \(|e|\) / mean \(e\)**.

| Segment | Original fast D | Original fast V | Fresh A_large D | Fresh A_large V |
|---|---:|---:|---:|---:|
| Plateau | .1281/+.0776 | .1292/+.0783 | .1266/+.0984 | .1279/+.1006 |
| Rise | .0941/−.0367 | .0936/−.0303 | .1477/−.0205 | .1434/−.0153 |
| Fall | .1288/+.1026 | .1307/+.1048 | .1876/+.1552 | .1950/+.1586 |
| Reversal | .1166/+.0572 | .1216/+.0612 | .1737/+.1106 | .1796/+.1172 |

Original fast means cases **36, 37, 39, 40**. Sustained under-estimate median/worst: original D **.2960/.5596**, V **.2938/.5637**; fresh D **.3171/1.5259**, V **.3211/1.5259**. Per-case identification is also in the JSON.

Late push/pull results follow. Peak and RMS are in N; RMS and stiffness errors use **event start through end + 0.6 s**. Beta is **D/V before the event**; error cells are **mean/minimum**, D then V.

| Direction/run (sign) | L H/D/V | Peak H/D/V | Beta D/V | RMS H/D/V | Error D; V | Mean g |
|---|---|---|---|---|---|---:|
| Toward N0 (−) | .296/.298/.297 | 9.372/9.372/9.372 | 0/0 | .542/.546/.545 | −.094/−.351; −.091/−.351 | .028 |
| Toward N2 (+) | .088/.087/.087 | 9.246/9.240/9.240 | .0515/.0515 | .238/.237/.236 | .174/−.114; .174/−.114 | .013 |
| Toward N6 (+) | .379/.338/.321 | 9.641/9.739/9.741 | .1502/.1502 | .662/.682/.644 | −.008/−.453; −.030/−.480 | .217 |
| Toward V0 (−) | .102/.102/.102 | 8.607/8.607/8.607 | 0/0 | .210/.210/.210 | −.066/−.314; −.066/−.314 | 0 |
| Toward V2 (+) | .054/.055/.055 | 8.454/8.466/8.465 | 0/0 | .159/.164/.163 | .057/−.159; .058/−.159 | .002 |
| Toward V6 (+) | .126/.115/.114 | 8.959/8.825/8.822 | .2595/.2604 | .395/.376/.374 | .105/−.127; .104/−.128 | .032 |
| Away N0 (+) | .258/.258/.258 | 9.313/9.313/9.313 | 0/0 | .465/.465/.465 | −.100/−.352; −.100/−.352 | 0 |
| Away N2 (−) | .0867/.0865/.0864 | 8.699/8.702/8.700 | .0515/.0515 | .235/.236/.235 | .158/−.114; .158/−.114 | .005 |
| Away N6 (−) | .333/.450/.267 | 9.282/11.059/9.281 | .1502/.1502 | .549/.903/.530 | .003/−.525; −.014/−.360 | .159 |
| Away V0 (+) | .101/.101/.101 | 8.581/8.581/8.581 | 0/0 | .203/.203/.203 | −.063/−.350; −.063/−.350 | 0 |
| Away V2 (−) | .055/.055/.055 | 8.496/8.496/8.496 | 0/0 | .163/.163/.163 | .058/−.121; .058/−.121 | 0 |
| Away V6 (−) | .111/.114/.124 | 8.731/9.045/8.962 | .2595/.2604 | .347/.374/.404 | .081/−.392; .099/−.380 | .050 |

Mechanism statistics, pooled over finite ticks:

| Tick subset | Mean g | Median τw (s) | Median sv/kw |
|---|---:|---:|---:|
| Change ticks: 36, 37, fresh N5, V5 | .10609 | .09275 | 3.84654 |
| Event ticks: six runs with pre-event β > 0 | .07184 | .09740 | 2.41378 |

- **V1 passes.** Original/fresh A_large mean-\(L\) differences from D are **−.000432 / −.005722**. Largest rise/fall mean-error change is **.00634**; largest median-\(|e|\) increase is **.00741**.
- **V2 fails overall.** N6 improves substantially: away \(L\) **.450→.267**, RMS **.903→.530**; both N6 losses are below .40. But away V6 worsens: \(L\) **.1144→.1244**, RMS **.3741→.4038**, exceeding the prescribed H/D-average RMS bound by **.04340 N**. Away N2 also exceeds that bound by **.000199 N**.
- **V3 passes.** Maximum \(|\Delta L|\) versus H across the original 24 and fresh 12 A runs is **.002781**, below .003.
- **V4 fails.** Crossing mean \(g=.106\) meets ≤.3, but learned-event mean \(g=.0718\) is far below ≥.5. The predicted stronger event weighting did not occur.

The `m0` audit found **0 violations in 81,718 unclipped comparisons**, maximum relative difference **\(2.014\times10^{-13}\)**; 102 calls were clipped. There were 20 specified fallbacks. All **59,004 beta-zero ticks** returned `k_w` bitwise. With beta fixed at zero throughout case 36, all five metrics matched H exactly, and force/truth traces matched bitwise.

The method preserves crossing performance and helps N6, but **does not reliably suppress push/pull errors or deliver the predicted mechanism separation**. Large sustained identification errors also remain. Missing repository reference traces were read from the prior round-4 mirror after verifying both reference JSON files were byte-identical to the repository copies.