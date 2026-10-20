# Separate candidate design, rapid validation and iteration

2026-10-19; start`a25d887`. The user requests independent candidate design/tests toward leading-journal methods, with improvements where needed. Comparisons serve design, not an all-baselines-first approval gate. See[design/iteration](../../research/separate_candidate_design_and_iteration.md).

**Reduced prototypes with actual numerical iteration; neither leading-journal contribution nor full B0 is established.** NNBO/Jin predictiveVIC comparators are explicit scan adaptations, not complete reproductions. Do not generalize to wholeNNBO, E41, B0 or candidate success/failure.

## Current routes

- Candidate 1: shared reducedMHE→contact-predictiveOCP→neighboring/curvature reuse and correction→actual gain rates. DirectOCP temporarily replaces unimplementedHJB roots; compare sufficient ordinary warmMPC separately, without B0 certificates. Compatible-model weighting was tested, slower with insufficient benefit, and not selected.
- Candidate 2: shared reducedMHE plus motion/actual gains/energy/task inputs→network→actual gain rates. Simulation task-gradient training; deployment reads no material truth/future. Independently train state/arrived-force-only ablations.

Shared ground truth is the previously corrected finite-footprint foundation model, not B0 cylindrical Hertz two-region contact or calibrated tissue; moments are not a complete 6D stress model. Shared MHE is a local seven-parameter moving-window fit, not B0 full-state/bias MHE. The initial surface reference is known; no ground truth is read thereafter.

## Iterations and reasons

|Stage|Finding/revision|Boundary|
|---|---|---|
|Target gains/sensitivity|Expensive curvature/cross-derivative refresh tails without stable warm-start speed benefits; use directional residuals/cached curvature|Established quasi-Newton/ENE implementation, not new general solver|
|Initial joint NN|Task-gradientNN trained, but stronger comparator tuning reduced benefits|Joint compensation gains are not K/D gains|
|Pure gains/task weights|Retainw=1, addNNBOQ/R-relatedw=12.5,25 and retune each comparator|Numerical tradeoffs, not clinical calibration or72 independent cases|
|Actual gain rates|Far-targetKn gradients zero on32 states; rate gradients nonzero; return to B0 rate actions|Repair of B0 deviation, not new discovery; common bounds/rates/energy|
|State coverage|5400 contexts had zero range for all four gains; add varying-gain trajectories on same72 fields, total10800|72 physical scenarios, not144 or10800 independent trials; no truth inputs|
|Feedback/derivatives|Align100Hz sampling/smoothing,20ms delay, and compiledAD with independent checks|First20ms future queue still initialized from arrived force, not exact full-history prediction/real-time certification|
|Projection|Two nonnegative stiffness-power coefficients allow closed-form rate projection; check 1000 cases,393 active|Maximum rate difference6.6e-6; plant retains independent original projection; equivalence is not novelty|

Final`candidate2_sampled_ad*.pt`models train with/without front ends independently for2000 fixed steps on mixed contexts using sampling-consistent bisectionAD. Later closed-form projection is equivalent;`compile_model.py --bisection`rebuilds training kernels. Retain intermediate positive/negative results, not test-based selection.

## Comparisons and splits

Training: seed 9101,12 fixed-tuning/24 early-context cases; newseed9102,72 fixed/varying-gain scenarios. Comparator tuning uses6 trainingseed9101 cases; seed 9201 supports design changes. Independentseed9401 has24 scenarios, reported perw with shared noise. See`independent_protocol.json`. Do not select per-test winners or deletew=1.

`nnbo_adapted`retainsHC preidentification, single critic and Newton references, explicitly shifts drift/control forFd=3N, uses localARE initialization/declared position loops. Original zero-force/full identification/robot execution are not exactly reproduced; scope is this adaptation.

**The NNBO adaptation fails verification:** independentw=1 median 88/600 inversion failures, maximum 572/600, and failures at higher weights. Local equation/root agreement does not ensure full-trajectory nonzero-target feasibility. Isolate this baseline from contribution comparisons; large apparent gains cannot establish superiority. Repair task/policy/execution matching first.

`nnbo_adapter_diagnosis.py`replays all critics and derives inverse ranges:4597/4774 failures lack positive-contact roots under this shift/quadratic critic, not insufficient iterations. See design§9. This does not apply to the original zero-force range or prove candidate novelty.

`mpvic_os_adapted` implements E41 constant acceleration, impedance gains as inputs and an OS box-constrained QP, with exact two-dimensional active-set enumeration on the normal axis. The tangential axis uses a shared strong servo; execution uses the shared rate/energy interface. This is not its MS controller, original tank or full experimental reproduction. `rate_mpc` is an ordinary warm-started, long-iteration reference on this same task, not a reproduction of all E41.

Main tests use zero wait/design frequencies: model12.5Hz, cheap feedback/network100Hz, not inherentMPC speed limits. Separately check 100Hz requests/measured releases, without certification. Comparators/ablations not usingMHE pay no useless overhead, though diagnostics generate estimates. Timing requires single-task runs; concurrent-training measurements are developmental only.

Intervals reflect24-case variation for fixed training/implementation, not seed/hardware/out-of-family uncertainty; exploratory intervals lack multiple-comparison adjustment. Since seed 9401 informed this iteration, new revisions need unseen scenarios rather than relabeling this cohort blind.

## Explicit B0 differences

Current costs include force, tangential/normal motion, velocity/effort;**B0(33)positive-definiteQg/Rg and terminal terms are absent**. Full hard trajectory constraints, robot dynamics, strict positive gain bounds and two-region Hertz fronts are incomplete. Approximate/saturated outcomes neither decide B0 nor invalidate the question. Restore these relationships; cross-interface comparisons separate common physical objectives from family-specific gain preferences.

These are iteration/carrier-adaptation results requiring real-interface/task validation. Platform information was requested; no reply does not imply no platform or justify invented clinical thresholds.

## Reproduction

Versions are in`environment.json`. Rebuild CPU prediction kernels into designated temporary directories, not binaries as source:

```bash
python experiments/separate_candidate_validation/compile_model.py --rates --path /tmp/impedance_rate_fast
python experiments/separate_candidate_validation/compile_model.py --rates --true-plant --bisection --path /tmp/impedance_training_sampled
python experiments/separate_candidate_validation/checks.py
```

`experiment.py --help`offers weights, methods, seeds, frequencies, timing and new output names, rejecting overwrites. Separate model/truth training kernels; manifests check purpose and prohibit truth kernels in online optimization. Full second-order graph compilation exhausted memory; online value/gradient compilation plus separately charged occasional curvature replaced it. This is not runtime/hardware evidence.

Training example uses identical data/fixed budget; timing/logs need not be bitwise:

```bash
python experiments/separate_candidate_validation/learning.py --data mixed72_contexts.npz --iterations 2000 --rates --compiled /tmp/impedance_training_sampled/prediction.so --output replay_new.pt
```

Add`--no-front`for ablations. `summary.py`only summarizes existing outputs, without training/new experiments. Earlier experiments/frozen snapshots remain unchanged.
