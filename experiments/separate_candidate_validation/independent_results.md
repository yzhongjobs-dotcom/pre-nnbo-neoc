# Independent numerical results for the separate candidates

This is a reduced test model, not complete B0, an original-paper reproduction or a judgment of leading-journal readiness. There are 24 independent physical scenarios reused at each weight, not 72 independent samples. Percentages are scenario-paired means with bootstrap 95% intervals. Positive values mean lower aggregate candidate cost, not improvement in every force/motion component.

Main tables use zero release wait and each prototype's design frequency (model 12.5Hz, inexpensive feedback/network 100Hz); they do not establish that ordinary MPC cannot run faster. Actual computational release and matched-frequency checks are separate.

## Normal-pose weight 1

| Method | Mean aggregate cost | Mean force RMS/N | Mean tangential RMS/mm | Peak over all scenarios/N |
|---|---:|---:|---:|---:|
| nnbo_adapted | 0.67127 | 0.10045 | 0.007912 | 3.3684 |
| mpvic_os_adapted | 0.565822 | 0.017158 | 0.0077505 | 3.2392 |
| rate_neighbor | 0.598393 | 0.052729 | 0.015717 | 3.4492 |
| rate_mpc | 0.57656 | 0.017845 | 0.0083287 | 3.1491 |
| candidate2 | 0.563017 | 0.01476 | 0.01617 | 3.1115 |
| candidate2_no_front | 0.560947 | 0.018414 | 0.01549 | 3.1595 |

The NNBO adaptation has substantial inversion failures: this row is diagnostic only, and related improvement percentages are not contribution evidence. Original-paper/nonzero-target-force adaptation still needs repair and verification. The E41 row is limited to this directory's declared OS scanning adaptation, not a conclusion about the whole original paper.

| Paired comparison | Mean improvement/% | 95% interval/% | Improved scenarios |
|---|---:|---|---:|
| rate_neighbor_vs_mpvic_os_adapted | -3.9588 | [-6.8831, -1.3974] | 8/24 |
| rate_neighbor_vs_rate_mpc | -2.6136 | [-5.2427, -0.4915] | 10/24 |
| rate_neighbor_vs_candidate2_no_front | -4.8503 | [-7.8039, -2.2702] | 6/24 |
| candidate2_vs_mpvic_os_adapted | 0.7072 | [0.2748, 1.1972] | 19/24 |
| candidate2_vs_rate_mpc | 1.8761 | [0.6178, 3.5212] | 20/24 |
| candidate2_vs_candidate2_no_front | -0.1392 | [-0.2919, -0.0129] | 6/24 |

## Normal-pose weight 12.5

| Method | Mean aggregate cost | Mean force RMS/N | Mean tangential RMS/mm | Peak over all scenarios/N |
|---|---:|---:|---:|---:|
| nnbo_adapted | 5.60969 | 0.1911 | 0.0080907 | 4.3534 |
| mpvic_os_adapted | 5.77167 | 0.11356 | 0.008028 | 4.5034 |
| rate_neighbor | 5.68205 | 0.32927 | 0.019939 | 5.0097 |
| rate_mpc | 5.61019 | 0.32619 | 0.02469 | 5.0182 |
| candidate2 | 5.46106 | 0.34048 | 0.019954 | 4.5515 |
| candidate2_no_front | 5.68037 | 0.37022 | 0.012982 | 5.0109 |

The NNBO adaptation has substantial inversion failures: this row is diagnostic only, and related improvement percentages are not contribution evidence. Original-paper/nonzero-target-force adaptation still needs repair and verification. The E41 row is limited to this directory's declared OS scanning adaptation, not a conclusion about the whole original paper.

| Paired comparison | Mean improvement/% | 95% interval/% | Improved scenarios |
|---|---:|---|---:|
| rate_neighbor_vs_mpvic_os_adapted | 3.6806 | [-3.4795, 10.5466] | 14/24 |
| rate_neighbor_vs_rate_mpc | -3.2586 | [-7.4936, 0.1167] | 12/24 |
| rate_neighbor_vs_candidate2_no_front | -0.2648 | [-5.5539, 4.5002] | 13/24 |
| candidate2_vs_mpvic_os_adapted | 6.6283 | [1.2891, 12.4636] | 16/24 |
| candidate2_vs_rate_mpc | -0.2422 | [-2.6565, 2.1458] | 14/24 |
| candidate2_vs_candidate2_no_front | 3.1086 | [1.5765, 4.8420] | 23/24 |

## Normal-pose weight 25

| Method | Mean aggregate cost | Mean force RMS/N | Mean tangential RMS/mm | Peak over all scenarios/N |
|---|---:|---:|---:|---:|
| nnbo_adapted | 10.0499 | 0.22764 | 0.0081534 | 4.6065 |
| mpvic_os_adapted | 8.84358 | 0.43512 | 0.0092317 | 5.0166 |
| rate_neighbor | 8.73555 | 0.39248 | 0.029267 | 5.0097 |
| rate_mpc | 8.73773 | 0.39255 | 0.01468 | 5.0166 |
| candidate2 | 8.81032 | 0.42267 | 0.037336 | 5.011 |
| candidate2_no_front | 8.81919 | 0.42222 | 0.026019 | 5.0112 |

The NNBO adaptation has substantial inversion failures: this row is diagnostic only, and related improvement percentages are not contribution evidence. Original-paper/nonzero-target-force adaptation still needs repair and verification. The E41 row is limited to this directory's declared OS scanning adaptation, not a conclusion about the whole original paper.

| Paired comparison | Mean improvement/% | 95% interval/% | Improved scenarios |
|---|---:|---|---:|
| rate_neighbor_vs_mpvic_os_adapted | 2.2711 | [0.4380, 4.3739] | 16/24 |
| rate_neighbor_vs_rate_mpc | -0.0390 | [-0.6245, 0.5175] | 11/24 |
| rate_neighbor_vs_candidate2_no_front | 1.6697 | [0.1507, 3.3183] | 16/24 |
| candidate2_vs_mpvic_os_adapted | 0.7449 | [0.3020, 1.2355] | 18/24 |
| candidate2_vs_rate_mpc | -1.8226 | [-3.8158, -0.0938] | 12/24 |
| candidate2_vs_candidate2_no_front | 0.0843 | [-0.5134, 0.5584] | 23/24 |

## Matched-frequency and computation-wait checks

6 seed9201 validation scenarios use 100Hz requests and release according to measured local computation time; they are separate from the 24 independent zero-wait scenarios above. Front-end costs, preparation/tail time and unreleased updates must be included. This is not robot real-time certification. NNBO remains an adaptation whose inversion has not passed verification.

| Weight | Method | Median decision/ms | p95/ms | Requests skipped before release/3600 |
|---|---|---:|---:|---:|
| 1 | rate_neighbor | 1.399 | 2.872 | 41 |
| 1 | rate_mpc | 2.575 | 3.524 | 12 |
| 1 | candidate2 | 0.142 | 0.290 | 2 |
| 12.5 | rate_neighbor | 1.400 | 3.101 | 33 |
| 12.5 | rate_mpc | 2.539 | 5.012 | 61 |
| 12.5 | candidate2 | 0.136 | 0.241 | 0 |
| 25 | rate_neighbor | 1.411 | 3.310 | 33 |
| 25 | rate_mpc | 2.256 | 9.386 | 242 |
| 25 | candidate2 | 0.139 | 0.263 | 0 |
