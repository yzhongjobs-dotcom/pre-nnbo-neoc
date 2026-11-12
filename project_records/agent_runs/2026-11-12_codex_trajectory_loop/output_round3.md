**R1–R4 met; R5 failed.** Completed 120 episodes with eight workers. Repository files remained unchanged; no commit or push.

Sandbox fallback: [code](D:/pip_env/temp/trajectory_robust_round3/experiments/rossi_failure_study/trajectory_robust.py), [JSON results](D:/pip_env/temp/trajectory_robust_round3/experiments/rossi_failure_study/results/trajectory_robust.json), [48 traces](D:/pip_env/temp/trajectory_robust_round3/experiments/rossi_failure_study/results/trajectory_robust). Arms and constants were imported unchanged.

Throughout, pairs are **huber/delag**. N = neural_white, 28°; V = viscous, 14°.

**Blocks 1–2: fresh scenes.** All requested overlay cells existed and the scene builder accepted them. β columns describe `delag`.

| Scene | Block 1 L | Median β, t>3 s | Block 2 L | Maximum β |
|---|---:|---:|---:|---:|
| N1 | .71061/.49574 | .12797 | .07893/.07886 | 0 |
| N3 | .08084/.08163 | 0 | .05626/.05621 | 0 |
| N4 | .20870/.19669 | .06163 | .05305/.05207 | 0 |
| N5 | 1.55105/.35765 | .17035 | .05616/.05581 | 0 |
| N7 | .18894/.18900 | .03598 | .12286/.12281 | 0 |
| N8 | .43458/.28433 | .09276 | .08995/.09053 | 0 |
| V1 | .12984/.10140 | .33384 | .04531/.04535 | .02218 |
| V3 | .04834/.04832 | 0 | .04464/.04460 | 0 |
| V4 | .91549/.90325 | 0 | .02862/.18994 | 0 |
| V5 | .82874/.29135 | .95234 | .05850/.05821 | 0 |
| V7 | .11648/.11603 | .01704 | .11436/.11465 | 0 |
| V8 | .05990/.05847 | .27992 | .06733/.06618 | .04284 |

Block 1 mean L: **.439460/.260321**; median: **.198818/.192844**. All 12 scenes were fast: 15.846–25.454 mm/s. Pooled identification errors, \(e=\ln(k_{\rm ret}/k_{\rm eff})\):

| Segment | Huber median \|e\| / mean e | Delag median \|e\| / mean e |
|---|---:|---:|
| Rise | .15633/−.10530 | .15378/−.01609 |
| Fall | .31866/.28886 | .19409/.16605 |
| Reversal | .19385/.12949 | .19010/.12700 |

Crossing force RMS: **.603818/.386209 N**, using the round 2 masks.

**Block 3: new sensor seeds.** Column seeds are the indicated base plus case number.

| Case | Seed 920000+case: L | Seed 930000+case: L |
|---|---:|---:|
| 36 | .74344/.47468 | .70810/.50020 |
| 37 | .28240/.23992 | .30101/.23435 |
| 38 | .10196/.10213 | .10633/.10708 |
| 39 | .06131/.05304 | .07129/.05557 |
| 40 | .09323/.08232 | .11341/.09976 |
| 41 | .09108/.09199 | .08285/.08337 |

**Block 4: perturbation after each arm’s processing.** ε=0 values are the stored round 2 results.

| Case | ε=0 | ε=+.005 | ε=−.005 |
|---|---:|---:|---:|
| 36 | .79032/.49673 | .76209/.49674 | .82654/.49698 |
| 37 | .31572/.25011 | .30506/.24903 | .32552/.25359 |
| 38 | .10588/.10636 | .10448/.10500 | .10737/.10782 |
| 39 | .08205/.07581 | .07850/.07427 | .08555/.07758 |
| 40 | .14172/.12201 | .13896/.11829 | .14537/.12540 |
| 41 | .08508/.08531 | .08410/.08430 | .08610/.08636 |
| Mean | .253460/.189388 | .245534/.187939 | .262743/.191289 |

**Block 5: late events.** “+” preserves the original wave sign; “−” flips it. T = toward probe, A = away, classified from Huber’s onset force change. β is the last used gain before 6 s. Error entries are **mean/minimum** during [start, end+0.6 s]; RMS is force error in N.

| Run | L H/D | Peak H/D, N | β | e H | e D | RMS H/D, N |
|---|---:|---:|---:|---:|---:|---:|
| T N0− | .29588/.29710 | 9.372/9.392 | 0 | −.087/−.351 | −.094/−.352 | .542/.545 |
| T N2+ | .08779/.08614 | 9.246/9.244 | .05166 | .158/−.126 | .177/−.113 | .238/.234 |
| T N6+ | .37874/.32645 | 9.641/9.742 | .14531 | −.051/−.369 | −.046/−.400 | .662/.644 |
| T V0− | .10243/.10209 | 8.607/8.603 | 0 | −.066/−.313 | −.067/−.324 | .210/.210 |
| T V2+ | .05404/.05518 | 8.454/8.458 | 0 | .061/−.161 | .057/−.167 | .159/.165 |
| T V6+ | .12636/.11975 | 8.959/8.824 | .25704 | .082/−.188 | .105/−.175 | .394/.388 |
| A N0+ | .25764/.25852 | 9.313/9.303 | 0 | −.100/−.352 | −.099/−.353 | .465/.467 |
| A N2− | .08672/.08624 | 8.699/8.695 | .05166 | .147/−.126 | .157/−.113 | .235/.235 |
| A N6− | .33334/.38913 | 9.282/10.619 | .14531 | −.066/−.338 | .005/−.371 | .549/.784 |
| A V0+ | .10071/.10054 | 8.581/8.590 | 0 | −.063/−.350 | −.063/−.348 | .203/.203 |
| A V2− | .05490/.05479 | 8.496/8.490 | 0 | .058/−.121 | .058/−.120 | .164/.163 |
| A V6− | .11090/.10711 | 8.731/8.909 | .25704 | .088/−.248 | .084/−.355 | .347/.349 |

Reconstructed waves matched stored event components exactly: maximum difference **0**. Residual maximum second differences inside/elsewhere were, for both settings: i0 **3.932e−6/4.323e−6**, i2 **7.641e−7/1.236e−6**, i6 **1.037e−6/1.822e−6 m**, at 1 ms spacing. Source inspection confirmed fast events have no activation coupling; `act` and `fext` stayed zero. Only displacement and event/component metadata changed.

**Prediction verdicts** (\(\Delta L=L_{\rm delag}-L_{\rm huber}\)):

- **R1 met:** mean ΔL **−.179139**; worst **+.000784**; **12/12** improved or equal within .005.
- **R2 met:** mean ΔL **−.052669**; worst **+.000908**, below .02.
- **R3 met:** all three mean comparisons favored `delag`. Delag spreads: case 36 **.000252 < .293595** round 2 gap; case 37 **.004564 < .065612**.
- **R4 met overall:** event RMS increased **.018132 N**; mean ΔL **−.000534**; maximum peak **10.619 N**. Toward: ΔL **−.009756**, event RMS change **−.003504 N**. Away: **+.008688**, **+.039768 N**. Away N6 deteriorated individually by **.055796 L**. Only three of six scenes had nonzero pre-event β.
- **R5 failed both conditions:** V4’s ΔL was **+.161321 despite β=0**; maximum β elsewhere reached **.042839**. Fresh nominal scenes therefore did not preserve the round 2 behavior.