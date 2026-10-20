# Reduced MHE–impedance execution engineering precheck (paused)

2026-10-19. The user identified that tuned fixed impedance had been used as the main comparison without implementing NNBO and strong field methods; NN training was stopped immediately. Retain estimation/execution tools and actual positive/negative results here, **not as validation of candidate effectiveness, field advantage or leading-journal contribution**. The recorded resumption order prioritizes [field comparators](../../research/candidate_improvements_targeted_literature_and_design.md#4-user-correction-main-comparisons-must-use-field-methods), rather than resuming training first.

[Protocol and boundaries](design_and_protocol.md): local seven-parameter constrained joint receding-window estimation, not full B0 MHE; ordinary two-segment shooting, not HJB/NEOC; planned task-gradient-guided NN, not PPO. Ground-truth contact and gain/energy execution reuse the corrected earlier model read-only; old outputs were not overwritten.

Completed: two training-distribution pilots, fixed auxiliary tuning on 12 cases, one fixed/analytical trajectory each on 24 cases, and 3600 correlated training contexts. Training was interrupted; no actor weights, validation or independent test results exist. Model-gradient diagnostics on 28 correlated states are not 28 independent closed-loop trials. Zero-release-wait scores do not establish deployment real-time execution.

Core files: `core.py`, `frontend.py`, `run.py`; `tune.py` only reproduces fixed auxiliary tuning; `diagnose.py` reproduces the two local gradient diagnostics. `first_pilot*`, `tuned_fixed.json`, `train_contexts*` and `contexts.npz` retain original data; do not overwrite them with new defaults. First-pilot trajectories have 29 columns, later collected trajectories 34 (adding anchor/sensed). Indices are at the end of `run.py`; diagnostics use the first-pilot format. Dependencies: NumPy/SciPy/Numba/Torch; see `environment.json`.

Numerical checks (not full control validation):

```bash
python experiments/mhe_guided_impedance_pilot/checks.py
```

Checks cover homogeneous-field contact consistency, independent directional differences for task gradients, and causality: undelivered delayed samples cannot affect estimates. See `checks.json`. NN deployment does not read material truth; simulation training cost may use privileged information. The recorded future work still includes field baselines, measured timed release, integration/footprint convergence and strong task tests.
