# Separate task simulations for both candidates

Experiments follow the original question: how can arrived contact information become timely, effective impedance changes during continuous scanning? Implement candidates/tests first, then suitable paper methods, rather than choosing an MPC to beat. Parameters are checkable synthetic engineering settings, not tissue calibration/clinical standards.

## Methods and code

| Object | Code and role | Boundaries |
|---|---|---|
| Contact/execution | [physics.py](physics.py): B0 cylindrical Hertz two-region forces/moments, planar two-axis force/torque interface, four actual gains, energy ledger, gain/rate/energy constraints; costs include tracking, velocity, force, positive-definite gain/rate penalties and terminal terms | Synthetic2DOF, not full arm/calibrated tissue. True normal contact is quasistatic; surface motion is generated separately |
| Estimation | [estimator.py](estimator.py): joint two-region materials, boundary, constant-velocity surface and three biases in arrived windows, with arrival regularization | B0 constant-velocity/window-fixed-bias instance, not full stochastic per-time-state MHE/statistical certificates. Weak parameters depend on priors |
| **Candidate 1** | [local_hjb.py](local_hjb.py): finite-horizon local quadratic values/feedback, directional state/model neighboring Riccati updates, correction/nominal checks;[run.py](run.py)executes | Earlier66-term global collocation found no qualified root; explicitly changed to trajectory-local HJB/iLQR. Retains feedback/NEOC, not original B0 roots or inherited certificates; uses Gauss–Newton neighboring approximations |
| **Candidate 2** | [learning.py](learning.py): MHE/preview, motion/gains/energy inputs;64–64 network outputs four gain rates; offline task gradients through contact/execution; rolling feedback | Contact adaptation of YOPO training roles, not PPO/original YOPO. Truth evaluates costs only, not actor inputs; short-time true-field training uses local constant-velocity surfaces, not actual future sensors |

Candidate 1 is not merely open-loop: backward recursion gives timewise gradients, Hessians/feedback matrices, and current states modify gain rates. Optimality remains discrete/local. Candidate 2 is a trained network, not a hand formula; data, seeds and weights are retained.

All methods share gain/rate/actual-energy projection without extra candidate force servos. Force, pose, velocity, drive, contact, gain and energy constraints enter nominal checks/results. Sample checks are not continuous safety proofs. Model updates do not reset energy.

## Task-first experiments

- **Tracking:** stiffness, friction, surface-motion and noise/delay categories,6 each,24 untuned physical scenarios.1.6s trajectories,2ms integration,100Hz force,50Hz requests,10Hz estimation/model updates,20/40ms sensing delay.
- **Candidate 1 mechanisms:** disable neighboring prediction with identical two corrections; compare a more extensive same-value-function solve for numerical quality. These are mechanism tests, not external baselines.
- **Candidate 2 mechanisms:** train with/without material preview independently; three main seeds. Ablations retain current-region/surface estimates, removing forward differences/boundary, not all MHE or a reproduction ofE27.
- **Generalization:** separate unseen softening/faster-surface cases, not merged into the 24 in-distribution cohort.
- **Computation/execution:** zero wait isolates control. Separately, without concurrent jobs, add measured estimation/decision costs to release queues; use prior available feedback while pending, and record skips/expiry. Training-concurrent timings cannot prove real time.

Pointwise NumPy overhead motivated [fast_physics.py](fast_physics.py), compiling identical scalar formulas for all methods.1000 states differ at most7.1e-15. Four replay cases show maximum full-cost difference1.84e-6/force RMSE1.86e-7N, not bitwise trajectories; see [compiled_replay_check.json](compiled_replay_check.json). Preserve original 24-case zero-wait tables. Current timing uses only new single-task`compiled_timed_*`; old`timed_*`diagnoses overhead. Timing versions do not add independent scenarios. Distinguish small reproduction variation from these recorded differences.

Report force, normal/tangential RMSE, peaks and violations separately. Different controller families share physical tracking scores, not candidate gain/rate preferences. Use fixed scales/common sampling and report provenance/components; full family-specific OCP costs are separate.

48 training scenarios/3840 correlated contexts;600 fixed gradient steps, batch48. Developmentseed16201, independentseed16301, shiftseed16401. Actorseeds16111–16113, without best-test-seed selection. Initial descriptors lagged queries by up to80ms; propagate using their own velocities via [align_training.py](align_training.py), without future observations.

**Spatial-information ablation correction:** the first `*_without_preview_test` always replaced forward-region parameters with old-region values, even after new material became current contact, so it cannot isolate preview. Its 26.86%/6.17% figures were withdrawn to diagnostics and are not contribution evidence. `physics.current_equivalent` now converts the causally estimated current footprint to a uniform contact model, strictly preserving current Fn/Ft (checked to about 2.2e-16) while removing explicit spatial material variation. Candidate 2 was retrained as `candidate2_no_spatial_preview_seed16111.pt`; new seed 16501, 24-scenario `*_preview_audit_*` results supply the formal information ablation. This comparator does not match every moment/local derivative: it tests spatial parameterization as a whole, so benefit cannot be attributed solely to temporal anticipation. Old networks load with legacy features; `--legacy-old-region-mask` only rechecks the old error and is not a valid baseline.

Training encountered0×∞ contamination in inactive energy branches; all branches are now bounded with immediate nonfinite failure. Failed runs produced no main model. [checks.py](checks.py)independently checks contact integrals, energy QP, costs and directional gradients. Initial collocation failures remain in[critic.py](critic.py)/[collocation_diagnosis.py](collocation_diagnosis.py), not success evidence.

## External comparison selected after candidates

Choose Jin et al.,*Model Predictive Variable Impedance Control of Manipulators for Adaptive Precision–Compliance Tradeoff*(TMECH2023, E41), for directK/D and contact tradeoffs. [baseline.py](baseline.py)reconstructs multistep dynamics, tanh parameterization, soft state penalties and tank conditions, using analytical-gradientL-BFGS on the same problem;[baseline_check.py](baseline_check.py)checks source costs/gradients. This is an explicit scan adaptation, not original Panda or old OS as whole-paper reproduction.

Retain both-axisK/D freedom/tanks. Common nonzero normal target is fixed feedforward; predicted external force subtracts it and common known damping from measurements. All share actual gain/rate/energy limits. No unusedMHE overhead is charged to the baseline; retain native sensing. TuneQ/R on training, allowing separate axisR as the paper permits instead of one binding shared weight.

Train8 scenarios with seven logarithmic candidates per axis; selectR=(1e-10,1e-5), check jointly, freeze, then test24 cases. In8 extra execution-audit cases,5/640 optimizations report non-success but yield bounded executable local outputs; tank minimum 0.1J. Neither baseline nor candidates are exact global solves; source/cost/gradient checks are not all original performance reproduction.

NNBO's target-force/reference inversion is not successfully adapted to nonzero sustained contact; exclude failed adaptation as the main comparator without rejecting the paper. E27 informs design, but full joint-reference/impedance output/training is not reconstructed or replaced by a generic network. One suitable comparator does not exhaust the field.

## Experimental practices from sources

Our task defines objectives; papers inform evidence: E41§VI-A1/A2 separates tradeoffs/constraints; E27§IV-B/C, TableIII separates training, typical/unseen terrain; YOPO§III-A/B, TableI checks cost, deployment information/timing; NNBO§V/TablesII, V compare models/position–force tradeoffs. These passages were read, not newly audited full proofs. See[protocol.json](protocol.json). Source publication does not make all project experiments publication-ready.

## Reproduction

Use prior separate-candidate NumPy/SciPy/CasADi/PyTorch versions; oneBLAS/torch thread per process, totalCPU quota2 cores.

```bash
OPENBLAS_NUM_THREADS=1 python checks.py
OPENBLAS_NUM_THREADS=1 python run.py --method candidate1 --n 24 --seed 16301 --output replay_candidate1.json
OPENBLAS_NUM_THREADS=1 python learning.py --data training_aligned_contexts.npz --seed 16111 --iterations 600 --output replay_actor.pt
OPENBLAS_NUM_THREADS=1 python run.py --method candidate2 --model replay_actor.pt --n 24 --seed 16301 --output replay_candidate2.json
```

Scripts reject overwrites; timing reruns vary. Older experiments, papers, three B0 core documents and frozen snapshots retain original bytes at the recorded experiment checkpoint. Current numerical results do not replace hardware, material, cross-model or clinical evidence.
