The force term adds explanatory power in these recordings, but **F1–F3 were not met**. For A_large fast at \(h=9\), it raised median \(R^2\) from 0.623 to 0.770, narrowly missing the required increase of 0.15. The timing results do not support the predicted trajectory lag for that group.

The sandbox denied repository writes. The new [replay script](D:/pip_env/temp/force_term_replay_mirror/experiments/rossi_failure_study/force_term_replay.py) and [JSON results](D:/pip_env/temp/force_term_replay_mirror/experiments/rossi_failure_study/results/force_term_replay.json) are in a temporary mirror. Existing files were not changed.

Fits below are per-case, no-intercept OLS; entries are medians across cases. Each fit cell gives **\(R^2\); coefficients**. \(M\) denotes measured-only \(x_{1m}\); \(r\) and \(r_m\) are its correlations with \(x_2\), using truth-based and measured-only trajectory terms respectively.

| Group | \(h\) | \(x_1\) | \(x_2\) | \(x_1,x_2\) | \(M\) | \(M,x_2\) | \(r,r_m\) |
|---|---:|---|---|---|---|---|---|
| A_large fast | 9 | .623; .293 | .059; .945 | .770; .296, 1.394 | .670; .252 | .790; .250, 1.309 | −.118, −.124 |
| | 24 | .674; .288 | .040; .839 | .765; .279, 1.209 | .732; .256 | .800; .246, 1.144 | −.062, −.067 |
| | 48 | .636; .253 | .050; .680 | .682; .241, .970 | .696; .239 | .751; .229, .935 | .009, −.012 |
| A_large slow | 9 | .169; −.028 | .010; .599 | .213; −.006, .724 | .194; −.049 | .242; −.028, .720 | −.260, −.273 |
| | 24 | .227; −.012 | −.008; .361 | .231; −.012, .359 | .262; −.043 | .267; −.044, .326 | −.090, −.095 |
| | 48 | .238; −.031 | −.009; .175 | .241; −.037, −.424 | .273; −.063 | .277; −.070, −.534 | −.055, −.060 |
| A | 9 | .079; .028 | .119; .745 | .257; .036, .862 | .077; .027 | .254; .034, .859 | −.189, −.185 |
| | 24 | .153; .036 | .027; .738 | .244; .037, .899 | .151; .035 | .236; .036, .890 | −.077, −.071 |
| | 48 | .173; .036 | .031; .758 | .250; .036, 1.021 | .162; .034 | .234; .034, 1.009 | −.046, −.029 |
| A_slow_thin | 9 | −.054; −.018 | .091; .522 | .145; −.012, .534 | −.054; −.017 | .144; −.011, .529 | −.183, −.180 |
| | 24 | −.117; −.010 | −.047; .438 | −.024; −.009, .476 | −.115; −.009 | −.024; −.008, .468 | −.061, −.059 |
| | 48 | −.154; −.003 | −.139; .356 | −.088; −.003, .324 | −.151; −.004 | −.086; −.004, .315 | −.050, −.048 |

Median residuals \(T-\text{fit}\), shown as **rise / fall**:

| Group | \(h=9\), \(x_1\) | \(h=9\), both | \(h=24\), \(x_1\) | \(h=24\), both |
|---|---|---|---|---|
| A_large fast | .065 / −.030 | .046 / −.018 | .145 / −.030 | .105 / −.007 |
| A_large slow | .071 / −.068 | .066 / −.055 | .120 / −.107 | .114 / −.099 |
| A | .058 / −.052 | .049 / −.043 | .079 / −.082 | .072 / −.076 |
| A_slow_thin | .059 / −.018 | .052 / −.011 | .114 / −.010 | .112 / −.005 |

Timing is the median lag in ms maximizing correlation of five-tick differences; **positive means the predictor follows truth**. Columns are \(x_1\) / fitted combination / \(x_2\).

| Group | \(h=9\) | \(h=24\) | \(h=48\) |
|---|---|---|---|
| A_large fast | −80 / −104 / −124 | −56 / −44 / −152 | −4 / 0 / −128 |
| A_large slow | 60 / 12 / 12 | −124 / −100 / 4 | −80 / −64 / −4 |
| A | 64 / 8 / 8 | 60 / 8 / 0 | −12 / 4 / 0 |
| A_slow_thin | 56 / 8 / 8 | 184 / 0 / 0 | −4 / 0 / 0 |

Event-window results at \(h=9\), using A_large-fast median coefficients \((0.296,1.394)\):

| Group | RMS \(x_1\) | RMS \(x_2\) | RMS combination | Correlation |
|---|---:|---:|---:|---:|
| B | .192 | .022 | .073 | .244 |
| A+B | .236 | .010 | .067 | −.254 |

The combination is smaller than **unscaled** \(x_1\) in both event groups. Relative to the fitted trajectory component, however, the force term tends to add in B and slightly cancel in A+B.

- **F1 — not met:** \(R^2\) gain was **0.147**, below 0.15; the force coefficient was **1.394**, within range.
- **F2 — not met:** trajectory-only residuals were **+0.065 on rises / −0.030 on falls**, opposite the predicted signs. The combined rise residual remained **+0.046**, outside ±0.03.
- **F3 — not met:** at \(h=9\), \(x_1\) led by **80 ms** and the combination led by **104 ms**, under the stated lag convention.
- **F4 — met:** combination RMS was **0.073 vs 0.192** in B and **0.067 vs 0.236** in A+B.

Fits and timing exclude event windows; event statistics are reported separately. The contact-run definition, smoothing, aggregation, and per-case results are recorded in the JSON.