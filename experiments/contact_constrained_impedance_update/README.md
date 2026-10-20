# Are fast impedance updates worthwhile with actual contact and energy execution?

2026-10-18. Comparison `ed3a3299f5baf855685665bcf5d70226e759f0c4`. Two exploratory validation rounds under the original question, not full B0, hardware results or an established paper contribution.

## Question

Earlier linear-preload/expected-cost examples had small optimization gains erased by release delay. Add tangential/normal motion, unilateral nonlinear normal contact, four actual gains, rate limits and energy accounting to compare model optimization/learned impedance against tuned simple schedules. Conditions follow original interfaces, without increasing disturbances/tightening bounds after observing outcomes.

Neither round supports a leading-journal contribution. Learning improves timeliness, but ordinary RBF is similar/better and task increments are small. Inactive contact/force limits do not test benefits under difficult contact.

## Model, information and objective

`model.py` independently implements two translation axes and tangential/normal K/D. Normal force is

\[
F_n=3[\max(\delta,0)/\delta_d]^{3/2},\quad \delta_d=4.5/k,
\qquad F_t=\mu F_n\tanh(v_t/0.002).
\]

This 3/2 unilateral spring **is not B0 cylindrical Hertz finite-footprint contact**; no normal moment, spatial material field, MHE or full robot. Local slope k=300–1500N/m, friction0.1–0.6, target3N, tangential speed8mm/s. State contains two position errors, two velocity errors, four gains and energy; masses0.8/0.5kg. Old slope determines indentation; old/new ratio0.75–1.25. Gains/errors/ledger are jointly sampled. Only the postupdate0.3-second transition is modeled, with fixed model/no future randomness.

**All methods share exact current state, updated model and energy.** This ideal information does not solve environmental understanding; old/new coefficient sampling is not continuous scanning/deployment estimation. Reference/energy initial conditions are given, excluding preceding identification, reference changes and contact establishment.

Feedforward plus impedance feedback saturates actual force at±12N. Cost contains position, velocity, contact-force errors and changes in actual control, with scales/weights in `flow()`, not gain preferences, information rewards or imaging proxies. Optimize two constant target segments through shared low-level gain dynamics, not instantaneous K/D jumps.

K=5–1600N/m, D=2–120Ns/m; rates800N/m/s and80Ns/m/s. Energy filtering simplifies original(31):

\[
\dot E=\chi P_D-\tfrac12(e_t^2\dot K_t+e_n^2\dot K_n),\quad
\tfrac12(e_t^2\dot K_t+e_n^2\dot K_n)\le\chi P_D+\alpha(E-E_{\min}).
\]

\(\chi=0.5,\alpha=40/s,E_{\min}=10^{-5}J\). Euclidean projection onto a rate box/halfspace has feasible zero action. Exact continuous execution gives \(\dot{(E-E_{\min})}\ge-\alpha(E-E_{\min})\), preserving ledger lower bounds, with analogous gain bounds. Check numerical errors separately. This uses established constraints/energy accounting, not a new theorem or inherited full B0 proof.

Include power from moving references, feedforward and environmental changes. Ledger bounds do not establish whole-robot passivity, asymptotic stability or hard force safety.0.2–6N is a **simulation diagnostic range**, not a clinical or guaranteed trajectory limit.

## Round1: optimization, imitation and strong simple controls

`screen.py` writes `protocol.json` first:192 training,48 validation,96 test cases, seeds7201/7202/7203. Optimize fixed targets/critical-damping schedules on training; the latter changes normal K/D with slope. The30-target library includes holding old values, two tuned simple targets and27 stiffness/damping-ratio combinations, evaluated through identical complete execution.

Two-segment local L-BFGS-B starts from the library best under common physical objectives/bounds, without global-optimum claims.333/336 solves report convergence; retain failures/start-retention rules. Train three state-conditioned networks and an RBF on common optimizer labels; old PPO is not retrained.

Mean per-test-case relative improvement,**positive is better**, against tuned fixed targets:

| Method | Immediate release | Including locally measured delay | Median single-call time |
|---|---:|---:|---:|
| Tuned critical-damping schedule | 0.587% | 0.579% | 0.040ms |
| 30-target prediction library | 1.515% | −2.005% | 15.44ms |
| Two-segment local optimization | 1.837% | −10.313% | 74.53ms |
| Three target-imitation networks | 0.231%–0.309% | 0.218%–0.298% | about 0.047–0.051ms |
| Target-imitation RBF | −1.789% | −1.795% | 0.039ms |

Optimization timing includes library initialization; label production/training are separate. Library/network reselection adds tiny immediate benefits and remains worse with delay. These limits concern this implementation, not all RTI, sparse solvers or full B0.

## Revision: learn task costs rather than target gains

Different targets may produce identical actual gain trajectories. If energy is inactive and a component satisfies
\(g_d-g_0\ge\bar r(T+\tau_g)\), low-level \((g_d-g)/\tau_g\) stays saturated at\(\bar r\). Two sufficiently large targets then execute identically. A constructed example differs by72N/m in K and7.2Ns/m in D yet has zero cost/final-state difference. Label error need not equal decision loss; this does not attribute all first-round network error to this mechanism.

`ranking.py` reuses **30 actual execution costs** on the same192/48 training/validation cases. The network minimizes expected relative task cost over candidate distributions and selects candidates at deployment; quadratic regression/RBF predict candidate relative costs from identical data. Choose all hyperparameters/checkpoints on validation task loss, then generate96 new cases(seed 7503), without test retuning. This transfers mature task-oriented learning, not novelty from a named combination.

New results against tuned critical-damping scheduling:

| Method | Immediate improvement | With measured delay | Median single-call time |
|---|---:|---:|---:|
| 30-target prediction library | 1.168% | −3.092% | 15.58ms |
| Two-segment local optimization | 1.456% | −14.547% | 110.36ms |
| Quadratic cost regression | 0.050% | 0.003% | 0.232ms |
| RBF cost regression | 0.411% | 0.372% | 0.198ms |
| Task-cost network, three seeds | 0.191%–0.400% | 0.141%–0.352% | 0.232–0.252ms |

95 of 96 new local numerical optimizations report convergence. Networks have no reliable advantage over RBF: paired differences for seeds11/33 span zero, and seed 22 is worse. Against simple scheduling, the case-resampled 95% improvement interval is about [0.150%,0.615%] for RBF and [0.154%,0.583%], [−0.001%,0.312%], [0.136%,0.560%] for the three networks. These intervals resample this synthetic case family only, excluding robot, model-error, task and material-transfer uncertainty. Three training seeds are not three groups of physical trials.

![Benefits and timing on96 new cases](comparison.png)

Single-threaded, warm-cache, random-order local single-call timing includes features, candidates and selection, not throughput, WCET or robot guarantees. Hold old actual gains during calculation, then execute the original plan with exact submillisecond release handling. This is neither delay-compensated RTI nor receding replanning; negative delay outcomes do not condemn all online optimization.

## Stronger online execution: start simple, improve in background

Old-impedance waiting is diagnostic and insufficient as an online comparator. `anytime.py` adds established execution on the same96 cases: send critical-damping targets immediately, continue during computing, predict20ms release state/gains/energy, and optimize only remaining0.28 seconds. Include retaining the simple plan, choose by common cost, and retain it on timeout.20ms is an explicit engineering diagnostic, not a literature-mandated deadline or permanent acceptance rule.

Local optimization allows at most four L-BFGS-B iterations, reserving return time from measured single-evaluation costs and retaining the best evaluated candidate at budget expiry without requiring full convergence. RBF warm starts include measured feature, prediction and candidate-generation time; offline outputs are not treated as free. All methods share actual gain/energy filtering.

| Background method | Mean improvement over immediate simple execution | Median worker time | Over20ms |
|---|---:|---:|---:|
| 30-target prediction library | 0.780% | 14.17ms | 0/96 |
| Budgeted local optimization, at most four steps | 0.674% | 19.27ms | 0/96 |
| RBF warm start + budgeted local optimization | 0.758% | 19.20ms | 0/96 |

Independent DOP853 on the first two preselected cases gives maximum relative error2.88e−7. Including the baseline and predicting release correctly avoids old-plan waiting degradation under the same model. This is established anytime candidate improvement, not our theorem or real-system nondegradation/WCET guarantee. Added after observing prior results on the same96 cases,**not a third independent test**.

This positive result corrects the possible inference that delay always erases online benefits: execution design/RBF warm starts help. Yet all gains are below 1%, library prediction is best, and networks are unnecessary; leading-journal contribution remains unestablished. Retain this strengthening alongside earlier negative numbers.

## Independent checks and boundaries

`verify.py` independently reconstructs physics/Brent projection and compares zero-start SLSQP convex QP.64 random states plus16 energy-active algebraic examples give maximum field/cost difference2.16e−10, projection1.32e−13 and explicit-power storage identity1.39e−17. Algebraic activations are not performance samples or additional benefits.

Independent DOP853 on the first two preselected tests gives first-round relative difference8.63e−7, new 96 cases3.27e−7. Also check doubled resolution, weight reload and exact release times. Background checks are in `anytime_results.json`; JSON retains errors.

No performance cases lose contact, exceed6N or saturate actuators. Energy filtering activates only1/96 fixed-target cases in round 1 and none of the listed round 2 methods. **These do not test contact switching, limited energy or hard-force-bound advantages.** Algebraic energy validation does not demonstrate task importance.

Strong simple controls here stay within one impedance-execution family. Full Xue, Haninger, Fu, Roveda and capable direct force feedback with free reference/feedforward are not compared. Small mean gains do not establish a new field capability or leading-journal increment.

## Retained uses and status

Reuse four-gain execution, task-cost training, RBF and release prediction/background improvement. This known-local-model exploration does not integrate candidate 2's MHE/preview front end or prescribe future training, solving or contribution forms. Old continuation plans are withdrawn; see HANDOFF.

## Reproduction and material

Dependencies: Python3.12, NumPy1.26.4, SciPy1.13.1, Numba0.68.0, Torch2.7.1+cpu; plots Matplotlib3.9.4.

```bash
cd experiments/contact_constrained_impedance_update
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python screen.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ranking.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_ranking.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python anytime.py
python plot.py
```

Do not run other heavy jobs during timing. Retain source, protocols, results, two paired-array sets, six networks, regression weights and independent checks together. Reruns overwrite outputs; preserve comparison commits before new research changes. Historical hashes do not describe new results.
