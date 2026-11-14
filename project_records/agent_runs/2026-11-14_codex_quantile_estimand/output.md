The z = 1 quantile **reduced sustained underestimation but worsened mean control loss in most groups**. The fixed-margin comparison was dominated by extreme uncertainty values and drove the returned stiffness to its upper clip.

The requested files are in the temporary mirror because repository writes were denied: [runner](D:/pip_env/temp/quantile_estimand_mirror/experiments/rossi_failure_study/quantile_estimand_loop.py) and [results JSON](D:/pip_env/temp/quantile_estimand_mirror/experiments/rossi_failure_study/results/quantile_estimand_loop.json). The JSON contains 324 rows and no traces.

The prescribed **C = 2.01446 × 10¹⁹**, computed over 37,180 finite estimator ticks in matched original V2 runs. The median and 90th percentile of `exp(sig/k)` were **1.188** and **1.313**. Thirty margins exceeded one million; the largest was 7.44 × 10²³. Thus the nominal mean margin was determined by a small extreme tail.

Each loss entry below is **mean / median L**. Ripple is mean RMS for V2 / q10.

| Group | Huber | V2 | q05 | q10 | Fixed | Ripple V2 / q10 |
|---|---:|---:|---:|---:|---:|---:|
| A | .106/.078 | .106/.078 | .101/.064 | .156/.112 | 6.347/6.909 | .080/.045 |
| A slow thin | .591/.158 | .591/.156 | .466/.148 | .408/.239 | 6.899/6.909 | .078/.052 |
| Slow thin, excluding 7 | .133/.119 | .133/.119 | .125/.108 | .436/.211 | 6.909/6.909 | .085/.049 |
| B | .050/.038 | .050/.038 | .051/.050 | .372/.292 | 3.899/4.111 | .071/.028 |
| A+B | .077/.067 | .077/.067 | .073/.056 | .293/.076 | 5.352/5.491 | .083/.043 |
| A_large original | .254/.124 | .190/.113 | .187/.120 | .230/.092 | 6.183/6.909 | .083/.041 |
| A_large fresh | .440/.199 | .255/.192 | .193/.135 | .991/.516 | 6.176/6.701 | .075/.040 |
| Fresh A | .068/.057 | .068/.057 | .190/.057 | 1.019/.264 | 5.323/5.044 | .072/.028 |
| All 42 earlier scenes | .280/.104 | .218/.103 | .180/.092 | .492/.236 | 5.862/6.670 | .078/.041 |

Per-case L values below use the column order **Huber, V2, q05, q10, fixed**. N and V abbreviate `neural_white_28` and `viscous_14`.

| Scene | L values |
|---|---|
| Slow: 6 | .212, .212, .201, .278, 6.909 |
| 7 | 2.881, 2.881, 2.170, .267, 6.848 |
| 8 | .196, .193, .188, .211, 6.909 |
| 9 | .119, .119, .108, .103, 6.909 |
| 10 | .035, .035, .033, .047, 6.909 |
| 11 | .104, .104, .097, 1.541, 6.909 |
| Large original: 36 | .790, .500, .549, .715, 6.909 |
| 37 | .316, .256, .223, .373, 6.909 |
| 38 | .106, .106, .091, .086, 5.734 |
| 39 | .082, .073, .044, .052, 6.909 |
| 40 | .142, .121, .148, .098, 6.909 |
| 41 | .085, .085, .066, .054, 3.730 |
| Fresh N1 / N3 | .711, .492, .568, .685, 6.909 / .081, .081, .065, .056, 6.165 |
| N4 / N5 | .209, .196, .220, 1.307, 6.909 / 1.551, .359, .436, 2.006, 6.909 |
| N7 / N8 | .189, .189, .165, .153, 6.493 / .435, .280, .295, .346, 6.909 |
| Fresh V1 / V3 | .130, .101, .111, .270, 6.909 / .048, .048, .042, .040, 5.320 |
| V4 / V5 | .915, .915, .128, 2.731, 3.840 / .829, .228, .141, 3.176, 6.909 |
| V7 / V8 | .116, .116, .092, .083, 5.451 / .060, .058, .048, 1.035, 5.392 |

Identification entries are **plateau, rise, fall mean log error; sustained-underestimate median/worst**. Segment means pool valid ticks using each arm’s earlier-round masks; sustained values use the earlier 13-tick convention. “—” means no scored ticks.

| Group | V2 | q10 | Fixed |
|---|---|---|---|
| A | .041, .151, .289; .408/.660 | .396, .449, .805; .050/.235 | 3.573, 3.536, 3.597; 0/0 |
| Slow thin | .032, −.038, −.124; .380/.955 | .325, 1.351, .360; .054/.287 | 3.423, 3.741, 3.530; 0/0 |
| Slow excluding 7 | .032, −.022, .210; .344/.535 | .344, 1.351, .360; .032/.111 | 3.463, 3.772, 3.672; 0/0 |
| B | .202, —, —; .251/.479 | 1.539, —, —; 0/.149 | 5.529, —, —; 0/0 |
| A+B | .066, —, —; .237/.422 | .910, .434, —; 0/.025 | 4.859, 3.447, 3.605; 0/0 |
| Large original | .058, −.032, .104; .294/.564 | .365, .284, .596; .170/.364 | 3.524, 3.624, 3.987; 0/0 |
| Large fresh | .101, −.015, .159; .321/1.526 | .420, .515, 1.119; .066/.602 | 3.753, 3.712, 4.002; 0/0 |
| Fresh A | .101, −.038, .026; .292/.606 | .584, 1.683, 1.641; 0/1.732 | 4.540, 4.374, 4.572; 0/0 |
| All 42 | .096, −.015, .144; .295/1.526 | .772, .467, .941; .024/.602 | 4.402, 3.670, 3.976; 0/0 |

- **Q1 failed:** q10 mean L increased in all four specified groups; for fresh large it rose from .255 to .991.
- **Q2 failed:** B increased by .322 and A+B by .216, both beyond .02.
- **Q3 passed:** sustained-underestimate medians fell by at least one-third in every reported group.
- **Q4 failed as stated:** case 7 improved from 2.881 to **.267**, below the predicted .5 threshold.
- **Q5 passed numerically:** q10’s 42-scene mean was .492 versus 5.862 for fixed. It is a weak test of adaptive margins: extreme ticks made C enormous, and clipping destroyed the intended practical margin match.

The factor-zero check reproduced stored V2 L exactly on cases 0, 7 and 36: **max |ΔL| = 0**. All 54 V2 reruns matched stored metrics exactly; 23 matched reference reruns recovered missing ripple fields.

A cleanup command for temporary worker files was rejected by automatic policy review as “blocked by policy.” Those intermediate files remain under `D:\pip_env\temp\quantile_estimand_srfvmaal` and `D:\pip_env\temp\quantile_ripple_2hd9a3jx`; the delivered JSON itself contains no traces.