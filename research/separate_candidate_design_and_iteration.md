# Separate design and iteration of the two candidates

2026-10-19, starting from `a25d887`. The user's latest clarification is to develop methods that solve the original question with innovation and evidence sufficient for leading journals; comparisons support design. Experiment with and judge each candidate separately, and continue improving deficiencies. Completing baselines does not complete the research. Work resumed on this basis, without reinstating the old pause or introducing approval gates.

## 1. Required capabilities and design sources

The task remains continuous soft-tissue scanning: as new contact information arrives, rapidly update actual impedance to suit the current environment, robot state and task tradeoff, improve combined force and motion performance, and account for execution and computation. The task is neither abstract whole-acquisition planning nor accelerating a root detached from control.

Design sources were checked in the [targeted literature](candidate_improvements_targeted_literature_and_design.md): NNBO supplies the nonlinear-contact → performance-tradeoff → feedback-decision → reference-execution route; Jin E41 supplies direct impedance decisions, constraints and energy handling; preview ENE supplies solution reuse after state/model changes; YOPO supplies task-cost gradient guidance and the division between privileged training information and restricted deployment observations; Li E27 is an in-field learning neighbor combining random-terrain pretraining with variable-impedance feedback. Adding these module names together does not establish a contribution.

The specific adjustment is to include **delayed local environmental knowledge, current actual gains/stored energy, future contact consequences and executable gain transitions** in one decision object. The two candidates divide computation differently; each requires its own control/cost and mechanism evidence.

## 2. Shared test model and boundaries

Reuse the previously corrected finite-footprint contact model and actual execution, without changing old data. The shared front end is a reduced seven-parameter nonlinear constrained moving-window estimator. It takes 31 measured motion/delayed-wrench points and fits local log stiffness, material gradient, surface height/slope/velocity, viscous time and friction; the propagated previous estimate provides arrival regularization. This is neither complete B0 two-region MHE nor a calibrated tissue model. Current pose is known; force sensing has a 20ms delay at 100Hz, and estimation runs at 20Hz. The true field and local prediction field differ.

Target gains pass through actual rate limits, bounds and stored-energy projection; neither network nor optimizer can make K/D jump instantly. Robot integration uses 1ms and the true contact footprint has 65 points; optimization uses a 240ms horizon, 4ms/33 points and two action segments; learning cost uses 240ms and 2ms/33 points. Passing numerical-consistency, directional-gradient and energy-projection checks does not establish tissue-model correctness, passivity/force safety or global optimality.

This iteration uses two complementary instances: joint outputs of four gains plus normal compensation, and a learning version that outputs gains only. The latter avoids attributing reference/feedforward benefits to K/D. Paper baselines retain their original interface freedoms rather than being restricted to the same number of outputs.

## 3. Candidate 1: executable optimization and neighboring updates with contact prediction

### 3.1 Relation between the numerical implementation and the original candidate

This iteration tests control and reuse in the model route using a parameterized finite-horizon OCP, rather than first solving a high-dimensional critic nodal root. **This changes candidate 1's numerical implementation: a direct OCP optimality system temporarily replaces the unimplemented HJB critic root.** The original research question and estimation → optimal impedance → rapid update relationship remain. This is not a complete NNBO/HJB reproduction and does not automatically inherit B0 root-certification proofs. Reconnecting the critic depends on its practical control/computation role, rather than preserving an algorithm name.

Let $c$ contain current robot state, actual gains/stored energy, the estimated local contact field, delayed-sensing initialization and task parameters; $a$ contains target gains/compensation for two action segments. The optimization is

$$\min_{a\in[-1,1]^{10}}J(a;c),$$

Here $J$ integrates actual gain-transition dynamics, predicted contact and delayed force feedback over the horizon; bounded actions and the actual execution layer handle limits. This is a local finite-horizon problem, with no claim of global optimality or hard force safety.

### 3.2 First version and redesign after failure

The first version constructs sensitivity from $F=\nabla_aJ$, conditional on a locally fixed active set,

$$H\,\delta a=-F_c\,\delta c$$

with projection/model-cost comparisons and two correction steps. Nonoptimal nominal points and regularization make this a numerical approximation; exact-solution IFT guarantees do not apply. Full curvature/cross derivatives are refreshed every four updates, with all costs counted.

Measurements showed a fast median but expensive refresh tails, making the first version slower on average than ordinary warm starts. The claim that sensitivity necessarily accelerates computation was therefore dropped. The revision computes a directional residual from the parameter change that has actually occurred,

$$\Delta F=F(a_k;c_{k+1})-F(a_k;c_k),\qquad
a_{\rm pred}=\Pi(a_k-B_k\Delta F),$$

where $B_k$ is cached inverse curvature, corrected by BFGS using accepted action/gradient differences under the same new model. Positive-curvature checks, step limits and cost rejection prevent blind reuse of an old anchor. A few projected corrections follow. This reuses established quasi-Newton/neighboring-update mechanisms and **does not claim a new general solver**. Its value must appear as more timely scanning decisions with comparable control quality, rather than fewer iterations alone.

Cost decrease holds only for the current estimated model. The first iteration already showed that low current-force prediction error does not guarantee correct action gradients. Model mismatch, actual waiting and strong optimization references remain unresolved boundaries for this candidate.

## 4. Candidate 2: rapid impedance feedback guided by task cost

A 64–64 tanh network takes the shared estimate, current motion/actual gains/stored energy and task weights, and produces impedance targets. Simulation ground truth is used only for offline rollout cost and numerical gradients, never as a deployment input. Gradients pass through contact, delayed feedback, gain rates and stored-energy projection to train actual task consequences, rather than fit ideal-stiffness labels. This adapts YOPO's training division to contact; it is neither PPO nor a YOPO reproduction.

The first version jointly outputs four K/D gains plus compensation. It uses 3600 correlated contexts from fixed/analytic trajectories in 24 training scenarios and 600 gradient steps. It improved over the first paper-adapted comparator, but the advantage shrank after comparator tuning. The separate contributions of learning, preview and impedance have not yet been established.

The second version outputs gains only, retaining the same actual execution. Feasible gains provide initialization, followed by task-gradient training. Training expands to 72 scenarios with new seed9102 and 5400 causal contexts; repeated samples are not independent scenarios. Task-tradeoff inputs are added. A separate ablation uses the same training mechanism with robot state and available delayed force only, omitting MHE parameters/predictions. It tests the information role in this design and does not stand in for E27 or the old PPO.

## 5. Why study task tradeoffs, and how to avoid favorable-setting selection

The original force-error scale is 0.3N, position scale 1.5mm and normal-pose weight $w=1$, giving a squared position/force weight ratio $Q/R=40000$. This emphasizes steady force tracking, where adequate direct force control can approach good performance. Original NNBO simulations use $Q/R=500000,1000000$. Accordingly, $w=1$ is retained in advance and $w=12.5,25$ added. These are three tradeoffs on the same physical test model, reported separately without removing the unfavorable original group.

This is a numerical mechanism range from a specific paper, **not clinical scanning weight calibration or a new task requirement**. Changing the original NNBO zero target force to 3N must also be explicit: the scanning adaptation uses $u=f-F_d$ and adds $F_d$ to the drift, retaining a single critic and Newton reference structure. It is not a literal reproduction of the original task. In-field comparators are retuned for each weight rather than using original-group parameters to manufacture an advantage. Lower aggregate cost may accompany greater force error and cannot be described as simultaneous force/motion improvement.

## 6. Is a leading-journal contribution currently established?

This is method iteration with specific sources, executable relationships and separate experiments, rather than a finished leading-journal manuscript. The substantive claims to judge are whether new contact information improves impedance decisions; whether the two divisions of computation advance actual execution relative to strong existing capabilities; and under which task/model conditions these advances are reliable.

A temporary win over an adapted baseline, fast network inference, mature sources, module counts or mathematical presentation cannot establish the target standard. Original-paper reproduction/fair task adaptation, front-end ablations, strong rolling optimization, computational release, generalization and platform evidence each support their corresponding judgments. Missing evidence calls for further design/quick verification. The comparator must not become the research endpoint, and local negatives must not directly decide complete B0.

## 7. Subsequent iteration and current judgment

### Actual gain rates and physical projection

The long saturation region in target-gain execution was checked on 32 correlated training states: all task gradients for distant Kn targets were zero, whereas direct rate actions had nonzero gradients. The implementation returns to B0's existing actual gain-rate action and uses a bounded-target inverse map to realize the same hardware gain/rate bounds and original energy projection, without allowing unbounded candidate requests. Both stiffness channels enter stored-energy power with nonnegative coefficients; rate requests are already box-feasible, so energy projection can only touch the lower bound and admits a closed form. Across 1000 random states, including 393 energy-active cases, the maximum rate difference from the independent original bisection projection is about 6.6e-6. The actual plant retains its original implementation.

This adjustment improves action representation for both optimization and training. It faithfully implements the existing B0 design rather than inventing gain rates. The closed form has practical computational value; mathematical equivalence alone does not establish leading-journal innovation.

### Training coverage and prediction/execution consistency

In the original 5400 contexts from fixed trajectories in 72 scenarios, all four actual gains had zero range. Gain-state inputs alone cannot teach transitions. Actual constrained variable-gain trajectories from those same 72 physical scenarios were added, giving 10800 contexts and approximate four-gain ranges of 1215, 767, 142 and 89. Models with and without front-end information are trained independently; deleting inputs from one model at test time would be an out-of-distribution test, not this ablation.

The first prediction model also used continuous delayed-force feedback, while actual control included 100Hz sampling and smoothing. Prediction/training now align sampling, smoothing and the 20ms delay, while explicitly retaining the approximation that initializes the first 20ms sensing queue with the force currently received. Compiled AD was checked against independent Numba contact rollout costs and directional differences; interpreter slowness is not treated as the optimization bottleneck. Compiling the full large second-order graph exceeded memory; online value/gradient computation is now compiled, with occasional curvature costs counted separately.

### Compatible-model tests were not selected as the main route

Using the received window, fit low/high-stiffness local models still compatible with the observed wrench, then optimize weighted task cost. Also compare a single model with minimum observation residual to avoid attributing all benefits of multistart fitting to uncertainty handling. Changes were small and computation increased in six validation scenarios, so this was not selected as the main improvement. It does not plan future observations, restore the withdrawn posterior-planning question or claim a credible statistical posterior.

### How accumulated evidence changes the next step

Independent tests still show that fast prediction/small residuals do not guarantee a near-optimal actual policy; sensitivity reuse does not consistently beat adequate ordinary solving; and learning differences with/without the front end are insufficient to support a main contribution from spatial information. The roughly 3ms front end and compiled ordinary MPC can already run quickly on this test model. A general computational bottleneck sufficient to support a paper cannot be assumed.

Rechecking B0 Eqs.31–34 also confirmed that this test model omitted positive-definite gain preference Qg, gain-rate cost Rg, terminal terms and full trajectory constraints; its reference was fixed and its contact model was not B0 two-region Hertz. Thus this iteration cannot decide the original question or complete candidate designs, but “only local” must not shield the current implementation indefinitely. The next iteration must **restore the original task → estimation → decision → execution relationships and verify the test model against actual interfaces/public originals**, rather than append small modules to this short-horizon integral objective or select favorable weights.

Comparisons between NNBO position references and four-gain control should use shared physical task metrics. Within-family Qg/Rg near-optimality should be checked against strong optimization in the same family; penalties for one parameterization must not be imposed on another interface to manufacture an advantage. Actual platform/control frequency and target scales have been requested; no reply does not mean none exist. Independent formalization and public-benchmark work can continue without inventing clinical requirements. A complete method meeting leading-journal standards has not yet been established.

## 8. Move from results to the next design, rather than repackage this prototype

For the static linear approximation of the current simplified test model, let $k$ be the environmental slope, $e$ normal displacement error and $k_p$ the shared direct-force-feedback coefficient. At equilibrium, $(1+k_p)(F-F_d)=-K e$. Thus, for current quadratic force/position criteria alone, defining $K_c=K/(1+k_p)$ gives

$$J_{\rm static}=\Delta^2\frac{q_FK_c^2+q_z}{(k+K_c)^2},\qquad
K_{\rm static}^{*}=\operatorname{clip}\!\left((1+k_p)\frac{q_z}{q_F k}\right).$$

This is neither a new discovery, the B0 optimum nor a new topic. It only explains why a static structured method is already capable on this engineering test model and why constraint saturation under extreme tradeoffs weakens the action role of parameter information. This analytic variable-impedance method is now a principled validation reference, rather than considering only fixed gains or arbitrary networks. Compatible-model/optimal-single-model tests at intermediate tradeoffs are also checked separately, avoiding an overall verdict on uncertainty design from one insensitive range.

In six intermediate-tradeoff validation cases, the analytic structured method is close to ordinary optimization. Compatible-set improvements are concentrated mainly in one case and add computation; there is no evidence supporting this module on this test model. Mature structures should be retained and studied for adaptation rather than excluded because they are simple. Their existing capabilities must also not be presented as this project's new contribution.

The next iteration implements the two candidates separately:

- **Model route:** restore full B0 task criteria, reference and feasible-trajectory relationships, then decide between a critic root and high-performance direct OCP as the main numerical implementation. Retain ordinary optimization when it is already fast enough and do not invent a general bottleneck. Study how new contact models change dynamic impedance and timely execution. Solution reuse must demonstrate practical value against ordinary solving of the same complete problem. Do not omit B0's positive-definite Rg and then use singular/saturated local results to decide its validity.
- **Learning route:** anchor the method in capable physical/analytic impedance structures and study actual corrections under time-varying contact, references and constrained gain transitions. Include complete task cost and boundary-state coverage in training. Mature residual or value structures may be adapted; being established does not disqualify them, and adding a residual does not predetermine a contribution. The existing direct-rate NN can provide an implementation starting point, but this iteration's small advantage or front-end benefit under one tradeoff cannot be presented as universal capability.

Failed NNBO inversion adaptations have been isolated from contribution comparisons. Further analysis concerns how the original zero target force, policy representation/inversion range and actual inner loop transfer to this task. An algorithm name cannot establish coverage, and one failed adaptation cannot establish incapability. These next designs have not yet been implemented; revise them using subsequent evidence rather than present plans as achieved innovation.

## 9. Adaptation failure located by returning to NNBO Eq.52

This is not a case of Newton iterations being too slow. For this iteration's declared nonzero-force offset, quadratic critic and $x_e$, let $d=x_r-x_e>0$, $A=W_1$ and $b=W_1(x_e-r)+W_2r/2$. The inversion equation becomes

$$S(d)=f-F_d+\frac{Ad+b}{B_ed^n}=0.$$

When $A>0,b<0,n\ge1$, the second term has a finite positive upper bound: for $n=1$ its supremum is $A/B_e$; for $n>1$ the extremum occurs at $d_p=-nb/[(n-1)A]$. Thus, if $f-F_d$ is below the negative of this upper bound, **no root exists in the positive-contact domain**. The original $F_d=0,f\ge0$ has a different inversion range; this adaptation's missing roots cannot be attributed to original NNBO.

Stepwise critic replay over the independent 24 scenarios at w=1 exactly matched final weights and failure counts. Of 4774 inversion failures, 4597 satisfy the no-root condition above; see [diagnostic code](../experiments/separate_candidate_validation/nnbo_adapter_diagnosis.py) and [case results](../experiments/separate_candidate_validation/nnbo_inverse_range_diagnosis.json). This explains the current adaptation's failure; it neither refutes the original paper nor constitutes a new general theorem. Entering main tests after passing local root examples, without first covering the actual force-error range, was an execution error to correct.

The concrete design requirement is that a nonzero sustained-contact task cannot reuse the original policy representation/inversion interface merely by shifting a force target. The finite-horizon objective, policy representation and actual execution must be designed together. B0 already has a direct actual-gain-rate channel and a Hamiltonian with positive-definite Rg, so it does not depend on this failed force–position inversion. These relationships must actually be implemented. Predictive VIC also already makes direct impedance decisions, so removing inversion alone is not innovation. The contribution must come from actual progress across the new-information → dynamic-impedance → execution chain.

## 10. Task simulations after the user clarified the sequence

The user specified that experiments follow this research question and the two candidate designs: first implement and verify each method, then select suitable external comparisons. Compare fairly when suitable methods exist and do not force their inclusion otherwise. Latest implementations and actual results are in [dual-candidate task simulation](../experiments/dual_candidate_task_simulation/README.md). This section does not rename prototype work as complete B0 success.

The shared test model returns to B0 cylindrical Hertz two-region resultant forces/moments, a two-dimensional force/moment interface, four actual gains and an energy state. Positive-definite Qg/Rg, terminal terms and trajectory-constraint checks are added, and the previous extra direct-force servo is removed. Environmental estimation uses a joint moving-window instance with constant surface velocity and a fixed bias within the window. It is neither complete per-sample stochastic-state MHE nor real-tissue calibration.

Candidate 1 first tried a global quadratic basis with 66 terms; none of four horizon settings yielded a qualified root. This does not establish nonexistence or invalidate B0. It then **explicitly switched to a quadratic HJB/iLQR value representation near a trajectory**: backward recursion gives the value gradient/Hessian and state feedback; directional residuals between old/new models and Riccati recursion supply neighboring prediction, followed by a few corrections and nominal-trajectory checks. This retains value feedback and a NEOC role, but is no longer the B0 fixed-node root; the original implicit-function/certification proofs do not apply automatically. Preparing full-parameter differences with 24 extra solves was expensive and was replaced by directional recursion. The present Gauss–Newton approximate neighboring update is not a new general solver.

Candidate 2 is a 64–64 network conditioned on the shared estimate/preview and actual robot state, trained offline for 600 task-gradient steps. Training uses 48 physical scenarios and 3840 correlated contexts, with three seeds trained separately. Simulation truth is used only for cost; deployment inputs are available estimates. Models with/without ahead-material information are trained separately. Query-timestamp propagation and nonfinite gradients in inactive energy branches have been fixed.

Separate task tests used 24 untuned scenarios, six in each of four categories, plus preview ablations, candidate 1 neighboring-update ablations and a more thoroughly solved same-family reference, and three seeds for candidate 2. Eight additional scenarios test softer material and faster surfaces. Computation-wait tests ran independently as single tasks. Measured software time includes the front end, updates, caching and release; expired horizons are not queried further. All methods were retimed after compiling the shared scalar Hertz kernel, so old Python scalar overhead is not treated as algorithmic difficulty.

After completing the candidates, Jin TMECH2023 multistep predictive variable impedance was selected as a matching comparator. The source model, parameterization, soft-state penalties and energy tank were reconstructed, with numerical optimization changed to L-BFGS. Nonzero target force, sampling and hardware bounds are explicit scanning adaptations, rather than reproduction of the original Panda experiments. Q/R was tuned separately for both axes on training scenarios, retaining its two-axis K/D freedoms. Without MHE, it is not charged unused MHE overhead. The failed NNBO adaptation is no longer used as contribution evidence.

Results do not support a win: candidate 1 force RMSE is about 0.00854N, the candidate 2 three-seed mean 0.00646N and the Jin adaptation about 0.00399N. Jin also has smaller tangential error; both candidates have slightly smaller normal-displacement error. The first ablation incorrectly removed information already belonging to current material, so its 26.86%/6.17% claims were withdrawn. Replacing it with a uniform model matched to current Fn/Ft and retraining the ablation network, a new seed16501 cohort of 24 scenarios gives a candidate 1 full-cost paired improvement of 2.43% (interval [1.45%,3.43%]); candidate 2 at a fixed seed gives 3.74% (interval [-0.17%,8.92%], uncertain). The overall spatial-modeling effect cannot all be attributed to temporal anticipation. NEOC on/off differences span zero, and current updates are slower. Full original data, components, costs, intervals and limitations are in the [results](../experiments/dual_candidate_task_simulation/results.md). Continue improving methods from these failures; extra comparisons, low gain-preference cost or completing the workflow cannot establish the target standard.
