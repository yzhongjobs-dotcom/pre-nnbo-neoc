# Separate simulations and comparison results for the two candidates

24 physical scenarios, 6 per category. Three learning seeds reuse the same scenarios; these are not 72 independent physical trials. The main table assumes zero computation wait.

| Method | Force RMSE/N | Tangential RMSE/mm | Normal RMSE/mm | Shared physical score (lower is better) | Sampled constraint violations |
|---|---:|---:|---:|---:|---:|
| Candidate 1: local value feedback + neighboring updates | 0.008543 | 0.2370 | 0.1061 | 0.07080 | 0 |
| Candidate 1 without neighboring updates | 0.008487 | 0.2375 | 0.1062 | 0.07072 | 0 |
| More thorough same-family value solve | 0.008479 | 0.2373 | 0.1062 | 0.07065 | 0 |
| Candidate 2: seed16111 | 0.006577 | 0.2178 | 0.1074 | 0.06112 | 0 |
| Candidate 2: seed16112 | 0.006060 | 0.2156 | 0.1076 | 0.05921 | 0 |
| Candidate 2: seed16113 | 0.006753 | 0.2114 | 0.1073 | 0.06054 | 0 |
| Jin multistep scanning adaptation | 0.003986 | 0.1647 | 0.1087 | 0.05646 | 0 |

The shared physical score contains only motion, velocity and normal force, integrated with the same sampling; it excludes controller-family gain-preference/rate regularization. A large Jin cost under this candidate's complete OCP objective does not establish a win. Component errors show the tradeoff between normal position and force/tangential motion.

| Measured-wait release: 8 scenarios | Force RMSE/N | Tangential RMSE/mm | Normal RMSE/mm | Mean of per-scenario median decision times/ms | Skipped updates |
|---|---:|---:|---:|---:|---:|
| Candidate 1 | 0.009282 | 0.2366 | 0.1043 | 52.867 | 13 |
| Without neighboring prediction | 0.008887 | 0.2376 | 0.1045 | 32.280 | 0 |
| More thorough value solve | 0.010998 | 0.2370 | 0.1033 | 77.283 | 37 |
| Candidate 2 | 0.006748 | 0.2198 | 0.1058 | 0.099 | 1 |
| Jin scanning adaptation | 0.003820 | 0.1329 | 0.1074 | 0.720 | 0 |

Timing uses CPU prototypes with a shared compiled scalar contact kernel and a 2ms release grid; it is not robot real-time certification. Old timed_* files diagnose the pointwise NumPy version and do not establish current speed. Model updates and individual network inference have different granularity, so their two median times are not an end-to-end speedup ratio. Estimation overhead is reported separately and included in release delay.

See [summary.json](summary.json) for mechanisms, stratified data, confidence intervals, distribution shifts and training-seed scope. Narrow intervals reflect scenario variation for fixed implementations, not attainment of leading-journal standards.

Per-case components and example trajectories: [PNG figure](results.png), [PDF figure](results.pdf).

## What these results establish

The preview ablation was checked and redone: the old version also removed new-material information after it entered the current region, so the original 26.86%/6.17% is not contribution evidence. The correction retains a uniform model matched to current-contact Fn/Ft and removes spatial variation only; the learning version is independently retrained. Comparison uses 24 new seed16501 scenarios without selecting parameters from the original test.
This comparison matches current Fn/Ft, not every moment or local derivative. It tests the overall role of spatially parameterized information; benefits cannot be attributed solely to temporal anticipation.

Candidate 1 corrected preview ablation: mean paired improvement in complete task cost is 2.428%, paired 95% interval [1.445%,3.433%]. NEOC on/off still gives only about 0.007% improvement, with interval [-0.260%,0.207%]; neighboring prediction is currently slower than disabling it under the same budget. Information and neighboring-update value must be judged separately. More thorough same-family solving is not ground-truth globally optimal control; tiny measured wins cannot mean surpassing the optimum.

Candidate 2 corrected preview ablation (same training seed16111, separately trained): paired complete-cost improvement is 3.740%, paired 95% interval [-0.168%,8.922%]; this is conditional on a fixed seed. The three main seeds have mean shared physical score 0.06029 versus 0.05646 for the Jin adaptation, about 6.78% higher. Fast individual inference or a one-seed ablation does not establish overall superiority to in-field methods.

Both candidates have slightly smaller normal-displacement errors but larger force/tangential-motion errors; this is an actual tradeoff. No leading-journal-level method advantage is established. Across 8 softer/faster-surface scenarios, candidate 1/2 mean force RMSE is about 0.01792/0.02482N. This extra cohort has no external baseline and does not establish superiority to prior work. Constraint checks concern sampled trajectories; energy constraints are inactive in the main scenarios, so these results do not establish active-constraint behavior or real-contact safety.

Next explain and improve why candidate 1's extra prediction/correction does not yield better cost than ordinary feedback updates, and why candidate 2 uses preview yet trails strong force/tangential performance; check training objectives against actual task tradeoffs. Existing strong impedance structures or learned corrections may be improved, with explicit changes verified on new scenarios. Changing scores, inventing computational bottlenecks or counting added modules cannot establish success. The actual platform/calibration still needs verification.
