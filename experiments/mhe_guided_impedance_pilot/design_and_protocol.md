# MHE front end/task-cost-guided impedance: first executable design

2026-10-19; comparison`9cd2672`. The user authorized search, design, quick tests and iteration, not a new approval gate. FrozenB0/old data remain unchanged. This reduces both candidates, not full B0 or leading-journal completion.

**Status after user correction: paused.** Fixed/analytical tuning dominated validation without implementedNNBO/strongVIC. Training stopped with no completedactor/independent tests. Retain engineering prechecks, not effectiveness/contribution evidence. Recorded resumption prioritizes[field methods](../../research/candidate_improvements_targeted_literature_and_design.md#4-user-correction-main-comparisons-must-use-field-methods), with fixed control auxiliary. Confirmed requirements were neither lowered nor replaced.

## Adaptations from mature work

YOPO backpropagates numerical task gradients rather than only imitating labels; Roveda2022 already approximates predictive impedance with actors. Replacing optimization by networks is not separately new. Adaptation: network sees delayed/noisy joint estimates/current execution states; offline gradients pass through contact, actual gain transitions, energy projection/feedback, training executable rolling actions. Truth may score training, but deployment cannot read materials/future contact.

Jointly fit finite-footprintFn/Ft/My and measured trajectories for local material/surface/motion. Both share this front. Candidate 1 first uses ordinary warm shooting as a capable model route, notHJB/NEOC. Candidate 2 uses task-gradient networks, with four gain targets/bounded normal compensation and pure-gain variants for attribution. Experiments determine joint-action value, without requiring gains to win.

This may reduce losses from model error, slow optimization and delayed gain execution; compare actual increments against nearest prior/strong simple control. No new theorem or first joint estimation/gradient learning is claimed.

## Executable interfaces

- Truth: reuse corrected finite cylindrical foundation, heterogeneous materials, moving surfaces/material convection; fixed attitude, not6D/calibrated tissue. Old files unchanged.
- Front end:31 samples at 100Hz,20ms wrench delay,20Hz constrained nonlinear window fitting. Seven parameters: log stiffness, gradient, surface height/slope/velocity, viscous time, friction; propagated previous estimates provide arrival regularization. Known pose inputs, not joint robot-state/full B0MHE. Known initial surface reference, no later truth. Linear local fields mismatch true sine/hyperbolic fields.
- Execution:100Hz force feedback,12.5Hz decisions, shared rate/energy projections;1ms/65-node truth. Prediction240ms,2ms/33nodes, two120ms segments. Roll future sensing; initialize first20ms queue with arrived force, not future-oracle data.
- Objective: previous engineering force, motion, velocity/effort scales, normal-pose weight1. Preserve prior distributions without artificially increasing difficulty. Not clinical thresholds.
- Auxiliary controls only: delayed force feedback/tuned fixedK/D, shared-front analytical force/reference compensation, ordinary warm finite-horizon optimization. NotNNBO/strong fieldVIC and insufficient for main contribution comparisons; the user corrected that arrangement.
- NN:64–64tanh, joint estimates/actual robot state/gains/energy/task inputs,10 bounded outputs; brief feasible analytical initialization then numerical task gradients. NotPPO/YOPO reproduction.

## Splits and decisions

Trainingseed9101; first two pilots from training. Fixed control tunes12 cases, not tests. Collect24 fixed/analytical trajectories; validationseed9201 supports design; newseed9401 reserved after freezing. Repeated trajectories are not independent cases; share noise within scenarios. Viewed/reused batches become validation, not blind tests.

First check local-model consistency with homogeneous ground truth, task-gradient direction and causal delay. Then compare actual closed-loop cost/force/motion, energy/gain execution and estimation/decision costs. Finally check relevant ablations, integration resolution and delay sensitivity. Current zero-release-wait results show algorithm/closed-loop potential only; actual computation or explicit sensitivity must be accounted for, without claiming deployment real-time guarantees. Complete MHE, nonlinear tissue memory, hardware and full strong-paper reproductions remain outside this prototype scope.

Completed: two trainingfixed/analytic/mpc pilots,12 fixed-tuning cases,24 paired trajectory collections/3600 correlated contexts. User-pausedNN has no completed model; candidate 2 training success is not established. Validationseed9201/testseed9401 were not run.28 correlated-state gradients diagnose possible local action errors, not28 independent loops.
