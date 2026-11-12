**S1–S4 all met. V4 was fixed; the late pull-away weakness remains.** Completed 78 `delag2` episodes plus two matched V4 reruns, without tuning.

The sandbox refused repository writes. Deliverables are in [the round-4 mirror](D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study): [loop code](D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/trajectory_loop.py), [robust code](D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/trajectory_robust.py), [loop JSON](D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/results/trajectory_loop.json), [robust JSON](D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/results/trajectory_robust.json). All prior rows and traces were preserved. Repository files remained unchanged; no commit or push.

Throughout, triples mean **huber / delag / delag2**. N = neural_white, 28°; V = viscous, 14°.

**V4 diagnosis.** Both reruns reproduced the stored metrics exactly. Forces first differed by more than 0.2 N at **0.224 s**. Before/after that time, `F6/F24` ranges were **0.92614–1.12756 / 0.74148–1.16733**, and RMS values were **1.010033 / 1.000352**. Post-divergence force error had its dominant frequency at **20.765 Hz**, using a demeaned Hann periodogram, excluding DC. Peak forces were **8.474 / 10.040 / 8.474 N**. `delag` β stayed zero throughout.

Removing the ratio restored V4 L to **0.028616570**, supporting the proposed cause. On this `delag2` run, all **1,250** returned stiffness values matched `k_w` **bit for bit**, with β = 0 throughout.

**Original 30 cases.** Group statistics:

| Group | Mean L, H/D/D2 | Median L, H/D/D2 |
|---|---|---|
| A | .10636/.10709/.10594 | .07794/.07813/.07794 |
| A_slow_thin | .59119/.52318/.59080 | .15750/.15596/.15614 |
| B | .04990/.04981/.04990 | .03789/.03782/.03789 |
| A+B | .07729/.07689/.07729 | .06665/.06643/.06665 |
| A_large | .25346/.18939/.19063 | .12380/.11419/.11489 |

| A_large case | L, H/D/D2 |
|---|---|
| 36 | .79032/.49673/.49429 |
| 37 | .31572/.25011/.25978 |
| 38 | .10588/.10636/.10579 |
| 39 | .08205/.07581/.07487 |
| 40 | .14172/.12201/.12399 |
| 41 | .08508/.08531/.08508 |

**Blocks 1 and 2.** β is the `delag2` median over estimator ticks after 3 s.

| Scene | Block 1 L, H/D/D2 | β | Block 2 L, H/D/D2 |
|---|---|---:|---|
| N1 | .71061/.49574/.48392 | .12933 | .07893/.07886/.07893 |
| N3 | .08084/.08163/.08126 | 0 | .05626/.05621/.05626 |
| N4 | .20870/.19669/.19592 | .06155 | .05305/.05207/.05305 |
| N5 | 1.55105/.35765/.35511 | .17389 | .05616/.05581/.05616 |
| N7 | .18894/.18900/.18875 | .03308 | .12286/.12281/.12286 |
| N8 | .43458/.28433/.28109 | .09386 | .08995/.09053/.08995 |
| V1 | .12984/.10140/.10152 | .33621 | .04531/.04535/.04531 |
| V3 | .04834/.04832/.04834 | 0 | .04464/.04460/.04464 |
| V4 | .91549/.90325/.91522 | 0 | .02862/.18994/.02862 |
| V5 | .82874/.29135/.30743 | .95546 | .05850/.05821/.05850 |
| V7 | .11648/.11603/.11642 | .01843 | .11436/.11465/.11436 |
| V8 | .05990/.05847/.05835 | .28053 | .06733/.06618/.06733 |

Block-1 mean L: **.439460/.260321/.261112**; median: **.198818/.192844/.192332**. All 12 scenes exceeded 15 mm/s. Identification errors below use the unchanged masks and \(e=\ln(k_{\rm ret}/k_{\rm eff})\); entries are **median |e| / mean e**.

| Segment | Huber | Delag | Delag2 |
|---|---|---|---|
| Rise | .15633/−.10530 | .15378/−.01609 | .14768/−.02051 |
| Fall | .31866/.28886 | .19409/.16605 | .18758/.15515 |
| Reversal | .19385/.12949 | .19010/.12700 | .17372/.11063 |

Crossing force RMS: **.603818/.386209/.380231 N**.

**Block 3.** Seeds are the column base plus case number.

| Case | Seed 920000+case: L, H/D/D2 | Seed 930000+case: L, H/D/D2 |
|---|---|---|
| 36 | .74344/.47468/.46916 | .70810/.50020/.49470 |
| 37 | .28240/.23992/.24017 | .30101/.23435/.23120 |
| 38 | .10196/.10213/.10187 | .10633/.10708/.10628 |
| 39 | .06131/.05304/.05319 | .07129/.05557/.05598 |
| 40 | .09323/.08232/.08142 | .11341/.09976/.09957 |
| 41 | .09108/.09199/.09108 | .08285/.08337/.08284 |

**Block 5.** “+” preserves the original event sign; “−” flips it. T = toward probe; A = away. Peak and event-window force-error RMS are in N. β lists **delag/delag2**, last used before 6 s; Huber has no learned β. Event window remains [6 s, event end + 0.6 s].

| Run | L, H/D/D2 | Peak, H/D/D2 | β, D/D2 | RMS, H/D/D2 |
|---|---|---|---|---|
| T N0− | .29588/.29710/.29760 | 9.372/9.392/9.372 | 0/0 | .542/.545/.546 |
| T N2+ | .08779/.08614/.08684 | 9.246/9.244/9.240 | .05166/.05150 | .238/.234/.237 |
| T N6+ | .37874/.32645/.33790 | 9.641/9.742/9.739 | .14531/.15022 | .662/.644/.682 |
| T V0− | .10243/.10209/.10244 | 8.607/8.603/8.607 | 0/0 | .210/.210/.210 |
| T V2+ | .05404/.05518/.05504 | 8.454/8.458/8.466 | 0/0 | .159/.165/.164 |
| T V6+ | .12636/.11975/.11530 | 8.959/8.824/8.825 | .25704/.25949 | .394/.388/.376 |
| A N0+ | .25764/.25852/.25764 | 9.313/9.303/9.313 | 0/0 | .465/.467/.465 |
| A N2− | .08672/.08624/.08649 | 8.699/8.695/8.702 | .05166/.05150 | .235/.235/.236 |
| A N6− | .33334/.38913/.44986 | 9.282/10.619/11.059 | .14531/.15022 | .549/.784/.903 |
| A V0+ | .10071/.10054/.10071 | 8.581/8.590/8.581 | 0/0 | .203/.203/.203 |
| A V2− | .05490/.05479/.05491 | 8.496/8.490/8.496 | 0/0 | .164/.163/.164 |
| A V6− | .11090/.10711/.11438 | 8.731/8.909/9.045 | .25704/.25949 | .347/.349/.374 |

**Prespecified verdicts**, with ΔL = `delag2 − huber`:

- **S1 met:** block-2 maximum |ΔL| = **2.4654×10⁻⁷**, below .01.
- **S2 met:** original A_large mean **.190634**, only **+.001246** versus `delag`; block-1 mean **.261112**, only **+.000791**. Worst ΔL was **0** across original A_large and **+.000416** across block 1.
- **S3 met:** cases 0–23 maximum |ΔL| = **.002731**, at case 8, below .005.
- **S4 met:** overall mean ΔL **+.005805**. Toward: **−.008354**; away: **+.019964**. Mean event RMS, H/D/D2: toward **.36777/.36427/.36923 N**; away **.32707/.36684/.39071 N**.

**Remaining failure:** away N6 deteriorated by **+.116525 L** versus Huber, with peak **11.059 N** and event RMS **.902704 N**. Removing the ratio fixes V4 and preserves the specified crossing benefit, but does not resolve this learned-gain pull-away failure.