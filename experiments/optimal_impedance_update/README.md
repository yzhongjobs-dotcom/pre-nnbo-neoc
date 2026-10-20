# Optimal impedance updates after change: value, speed and a method iteration

2026-10-18. Return to the confirmed [rapid control-beneficial near-optimal impedance question](../../research/optimal_impedance_question_and_methods.md). Completed literature-grounded derivation, model/learning comparison and finite-rate transient optimization, not full B0, hardware or leading-journal results.

## Model and comparisons

Independent two-coordinate linear contact near preload with four diagonalK/D. The second axis is lever-normalized rotation, not B0 tangential/normal Hertz sliding. With fixed known environment,

\[
M\ddot q+(B_e+D)\dot q+(K_e+K)q=w,
\quad u=-K(q+n_q)-D(\dot q+n_v).
\]

Use n=0 without noise; control-force noise also enters dynamics otherwise. Contact disturbances are \(K_eq+B_e\dot q\). Costs include displacement, velocity, contact force and actual control,**without artificial target gains, information reward or imaging proxies**. All start at common preload equilibrium with known models, not earned identification benefits.

Test white forces, band-limited correlated forces, and correlated forces plus band-limited observation noise. Spectra/amplitudes illustrate mechanisms, not calibrated tissue. See `protocol.json` for coordinates, ranges/training.

## 1. Stationary theory and mapping comparisons

For fixedK/D, Lyapunov equations give covariance/expected cost. Single-axis white-force interior optima are

\[
K^*=-k+\sqrt{k^2+(q_x+q_fk^2)/r},\qquad
D^*=-b+\sqrt{b^2+mA/r},
\]

\[
A=\frac{q_x+q_fk^2+r(K^*)^2}{k+K^*}+\frac{q_v+q_fb^2}{m}.
\]

This standard quadratic/stationary computation is not new theory but a required cheap comparator. Common noise-intensity scaling changes cost, not automatically optimal gains. Coupling, band-limited disturbance and observation error need separate calculations.

Each family has320 training,64 validation,96 independent tests. Three networks train2500 supervised updates on optimizer labels; RBF uses identical320 labels. This is **known-model supervised exploration**, not old PPO or direct unknown-tissue optimal inference. Teachers are bounded local solutions; two extra test starts do not establish globality.

| Test condition | Tuned-fixed mean cost gap from numerical reference | Analytical schedule gap | RBF gap | Three-network gaps |
|---|---:|---:|---:|---:|
| White forces | 0.1108% | 0.01124% | 0.0000128% | 0.0000300%–0.0000401% |
| Band-limited correlated forces | 0.2575% | 48.05% | 0.002451% | 0.002327%–0.002601% |
| Band-limited forces/observation noise | 0.9544% | 20.33% | 0.000348% | 0.000264%–0.000328% |

Analytical rules failing their assumptions do not defeat all analytical control; tuned fixed gains are stronger counterfactuals. Network/RBF approximate references; single-call RBF≈0.03ms, networks≈0.09–0.11ms. Randomized order reduces first-network cache effects; these are not batch-time divisions. Teacher/training preparation is separately counted.

White-noise full-state LQR provides an exact broader-feedback reference, unrestricted by four diagonal gains. Under colored/observation noise it is not this problem's optimal lower bound; beating it does not beat LQR theory.

**Static mappings are fast approximations, but strong fixed impedance is already near-optimal here, without major control increments.** ChangingK/D values do not establish a need for frequent exact updates; stationary results do not invalidate B0 finite-horizon problems.

## 2. Ordinary Newton, chord and sensitivity updates

`update_benchmark.py`:1152 updates per method, smooth sequences plus two model revisions, shared\(10^{-6}\)projected-gradient stopping. Refresh anchors every8 updates; derive sensitivities from common Hessian/mixed derivatives and count preparation, refresh, gradient checks, line search/fallbacks.

- Ordinary warm Newton≈1.49–2.23ms; fixed-Hessian chord≈0.69–0.89ms.
- Sensitivity prediction/chord correction≈0.76–1.01ms, no advantage after preparation.
- All updates reach common tolerance without fallback. Cheap gradients/Hessian reuse explain chord benefits in this small problem, not all NEOC, B0 finite roots or robot-chain timing.

Further RBF/network correction costs≈1–3ms with almost unchanged task cost. Distinguish root/control accuracy; stopping-threshold design is not a new research question here.

## 3. Newly changed environments with finite actual gain rates

`transient.py`: eight exploratory cases,0.3 seconds,12 gain-rate segments, actualK/D states. Initial gains/covariance use old stationary solutions; new environment/current distribution are known, while future noise is unknown to all. No terminal-gain targets/preferences. Compare holding old gains, tracking new stationary/analytical gains, and1/2/12-segment optimization.

| Four-case mean gap from best current12-segment reference | Hold old | Track new stationary | One-segment optimization | Two-segment optimization |
|---|---:|---:|---:|---:|
| White disturbance | 0.641% | 0.673% | 0.356% | 0.221% |
| Band-limited disturbance/observation noise | 1.734% | 1.936% | 0.413% | 0.148% |

One noisy case: holdingJ=0.121377, tracking new stationaryJ=0.123486, short-horizon reference J=0.117292. New stationary targets need not give best transitions. **B0/PPO already include states; this only rejects treating this static network as a complete dynamic solution.**

Two-segment worst gap is0.3551%; two12-segment solves hit150 iterations. Call it best current numerical trajectory, not a certified global gap. Eight mechanism examples do not establish population/clinical effects.

## 4. Solver revision driven by results

Original complex-step expansion integrates whole trajectories for every variable. `adjoint.py` uses forward second moments/backward adjoints under identical tasks, starts and budgets, borrowing established optimal-control structure.

|12-segment prototype, four-case mean | Complex-step solve | Adjoint solve |
|---|---:|---:|
| White disturbance |12.60 seconds|0.919 seconds|
| Band-limited disturbance/observation noise |15.87 seconds|0.779 seconds|

Same-integrator reevaluation gives maximum relative cost difference6.3e−7; the same two iteration-limited cases remain unconverged. Slow prototypes establish neither a field bottleneck, new solver nor robot real-time control.

`bounded_iterations.py` keeps two segments and uses adjoints for at most4 line-searched sequential LP steps, without tiny residual targets. Same eight cases end0.110%–0.483% above the 12-segment reference in47–82ms: seven use4 steps, one3. Include initialization, gradients, LP/search, excludeJIT/shared mechanics. This follows full-solve inspection, not independent testing.

**Usable improvements are blocking, structured derivatives and limited iterations, not exclusive NEOC/network advantages.** Test independent cases/full execution; single-state trajectory timing cannot replace B0 multistate feedback-reuse accounting.

`optimization_ms=0` in `transient_results.json` means no trajectory optimization, not free stationary targets/zero total delay. Their production costs are elsewhere; compare like1/2/12-segment solves only. All timings are local measurements, not WCET/certified real time.

## 5. Include computation delay in control outcomes

`delay_audit.py` retains the eight cases, costs/plans: old K/D controls during computation, then execute stored rates truncated at the same0.3-second horizon. Independent DOP853 gives:

| Mean improvement over holding old impedance | Plan immediately available | Release after measured47–82ms |
|---|---:|---:|
| White disturbance |0.331%|0.111%|
| Band-limited disturbance/observation noise |1.525%|0.192%|

Near-reference scores therefore are not deployment gains. Timing is part of the original rapid-update question, not a new question. **This diagnostic executes old-state plans, not delay-aware replanning or B0 postrelease feedback**; it does not show equivalent losses for those methods. Preparation/feedback separation, alignment and fast policies remain possible, untested in their full strength.

## Independent verification and reproduction

`verify.py` rebuilds physical matrices, checks costs by Kronecker covariance/frequency integration, gradients/Hessians by differences, reloads9 networks, and checks preselected transients/all original/budgeted methods with DOP853 plus independent gain/rate bounds. See `independent_checks.json` for nonzero errors.

Use root NumPy/SciPy/Numba plusTorch/Matplotlib(here2.7.1+cpu/3.9.4; see `environment.json`). Run from root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/screen.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/update_benchmark.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/transient.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/adjoint.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/bounded_iterations.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/optimal_impedance_update/delay_audit.py
python experiments/optimal_impedance_update/summarize.py
```

**Not covered: B0 energy tanks, hard force constraints, MHE, Hertz sliding, full robots or imaging.** Gaussian expected costs do not ensure samplewise force safety; fixed-gain stability does not ensure arbitrary time-varying stability. Retain solver/mapping/cost references, not a final topic or universal percentage threshold. Old reports do not prescribe current work.
