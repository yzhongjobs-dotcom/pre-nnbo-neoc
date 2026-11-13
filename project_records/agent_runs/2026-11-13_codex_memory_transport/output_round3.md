The gated Method passed the specified collapse threshold in both named scenes: **L = 0.223 for case 7 and 0.369 for fresh V4**. The correction reduced mean L in both large-rib groups, but its effects were uneven: it worsened fresh N1 by 0.215 and N5 by 0.360. Ungated `delag2` had the lowest mean L across all 42 scenes.

The original-case harness matched exactly: **max |ΔL| = 0** over all 30 cases. The sandbox refused repository writes, so the new [runner](<D:/pip_env/temp/gated_delag_round3/experiments/rossi_failure_study/gated_delag_loop.py>) and [results JSON](<D:/pip_env/temp/gated_delag_round3/experiments/rossi_failure_study/results/gated_delag_loop.json>) are in the temporary mirror. The JSON contains one row per arm and scene, with no traces.

Mean/median L, in arm order **huber / delag2 / gated / gated_delag2**:

| Group | huber | delag2 | gated | gated_delag2 |
|---|---:|---:|---:|---:|
| A | .106/.078 | .106/.078 | .129/.122 | .128/.115 |
| A_slow_thin | .591/.158 | .591/.156 | .164/.173 | .164/.172 |
| B | .050/.038 | .050/.038 | .086/.100 | .086/.100 |
| A+B | .077/.067 | .077/.067 | .076/.059 | .076/.059 |
| A_large original | .254/.124 | .191/.115 | .616/.261 | .269/.235 |
| A_large fresh | .440/.199 | .261/.192 | .546/.399 | .437/.305 |
| All 42 | .280/.104 | .220/.103 | .309/.123 | .228/.123 |

Per-case L and gated window acceptance rates (`gated / gated_delag2`):

| Case | huber | delag2 | gated | gated_delag2 | Acceptance |
|---|---:|---:|---:|---:|---:|
| 6 | .212 | .212 | .200 | .200 | .956/.956 |
| 7 | 2.881 | 2.882 | .223 | .223 | .974/.974 |
| 8 | .196 | .193 | .194 | .191 | .920/.926 |
| 9 | .119 | .119 | .110 | .110 | .992/.992 |
| 10 | .035 | .035 | .153 | .153 | .976/.976 |
| 11 | .104 | .104 | .106 | .106 | .765/.765 |
| 36 | .790 | .494 | 2.444 | .710 | .903/.958 |
| 37 | .316 | .260 | .556 | .263 | .960/.935 |
| 38 | .106 | .106 | .102 | .102 | .982/.982 |
| 39 | .082 | .075 | .316 | .265 | .931/.925 |
| 40 | .142 | .124 | .206 | .208 | .948/.918 |
| 41 | .085 | .085 | .069 | .069 | .942/.942 |
| N1 | .711 | .484 | .928 | 1.142 | .928/.950 |
| N3 | .081 | .081 | .066 | .069 | .983/.982 |
| N4 | .209 | .196 | 1.403 | .490 | .782/.883 |
| N5 | 1.551 | .355 | .968 | 1.329 | .837/.836 |
| N7 | .189 | .189 | .171 | .172 | .979/.990 |
| N8 | .435 | .281 | .745 | .481 | .857/.868 |
| V1 | .130 | .101 | .429 | .419 | .885/.928 |
| V3 | .048 | .048 | .045 | .045 | .977/.977 |
| V4 | .915 | .915 | .369 | .141 | .497/.566 |
| V5 | .829 | .307 | 1.130 | .666 | .570/.565 |
| V7 | .116 | .116 | .101 | .101 | .954/.954 |
| V8 | .060 | .058 | .192 | .192 | .578/.578 |

Identification uses the per-tick feed-forward `k_good`, with \(e=\log(k_\mathrm{good}/k_\mathrm{eff})\). Entries are pooled median |e| / mean e; sustained under-estimate is the median/worst across scenes using the existing 0.104 s convention.

| Scenes and arm | Plateau | Rise | Fall | Reversal | Sustained median/worst |
|---|---:|---:|---:|---:|---:|
| Original fast, gated | .143/.136 | .140/.059 | .416/.436 | .346/.326 | .303/.553 |
| Original fast, gated_delag2 | .178/.150 | .121/.061 | .292/.314 | .272/.228 | .328/.648 |
| Fresh, gated | .153/.170 | .179/.122 | .501/.560 | .297/.409 | .182/.560 |
| Fresh, gated_delag2 | .174/.189 | .205/.147 | .412/.467 | .318/.394 | .189/1.023 |

Predictions were recorded before the runs:

- **G1 met:** case 7 = .223 and V4 = .369, both below .5.
- **G2 not met:** mean L fell by .346 original and .108 fresh, but N1 and N5 worsened by more than .05.
- **G3 not met:** cases 1 and 4 changed by +.0108 and −.0141 respectively.
- **G4 not met:** all-scene mean was .2283 for `gated_delag2` versus .2195 for ungated `delag2`.

The gates avoided the two specified collapses, but did not prevent high L in other large-rib scenes. The correction improved pooled fall identification and group mean L, while worsening some individual crossings and the worst sustained under-estimate.