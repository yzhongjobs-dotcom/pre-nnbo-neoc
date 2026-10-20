# Understanding method design, improvement scale and contribution presentation through specific papers

2026-10-19. At the user's request, this analysis returns to original papers, understanding each paper first and only then comparing project failures. Sections 1–4 analyze papers independently; Section 5 connects them to the project. This report reconstructs publicly presented arguments, does not infer the authors' actual order of discovery, and does not infer that all equations are correct from publication status.

## 1. Reading scope and analyzed materials

The source is the [pinned sample repository](https://github.com/yzhongjobs-dotcom/academic-writing-skill/tree/c88d48eb08a876eb9b81a30d37b62d70fa028581/文献资料), commit `c88d48eb08a876eb9b81a30d37b62d70fa028581`. Earlier identity/abstract screening of 22 papers is not counted as new reading here. This round rechecked methods and evidence in ten local originals: all five T-RO papers, plus R01, R03, R15, R17 and R19. PDF pages and equation/figure numbers are given per case; the new scope is recorded in the [source list](../literature/top_venue_sample_sources_and_reading_scope.json).

The ten papers span T-RO, TMECH, TIE, TNNLS, Mechatronics and ICRA and cannot all be called leading-journal papers. Their limited years, author groups and tasks do not represent the entire 2026 field. This round did not reprove every theorem line by line, rerun paper source code or investigate every preceding work for each paper. “Existing foundations” primarily means explicitly cited, reused or compared material; “contribution interpretation” is this analysis. “Size of change” means structural change relative to those foundations, not invented code-line counts, development effort or discovery history.

## 2. Six detailed cases: how problems shape design and how large the contributions actually are

### R01: a small extension along the DMP route enables synchronized motion/stiffness generalization

Yang et al., TMECH2018, *A DMPs-Based Framework for Robot Learning and Generalization of Humanlike Variable Impedance Skills*, [original](../literature/papers/reverse_analysis_samples/R01.pdf), DOI10.1109/TMECH.2018.2817589. Checked PDF2–8, especially §II-C, Eqs.(13)–(18), Figs.4,7 and §III-B.

**Research question.** Demonstration should convey suitable stiffness evolution for cutting and pick-and-place as well as motion. The authors acknowledge existing DMP stiffness scheduling and ask how stiffness profiles can be obtained directly from human demonstration and generalized with target and duration, reducing repeated trial and error on new tasks.

**Inheritance and changes.** DMPs, EMG stiffness estimation, dual-arm demonstration and virtual-spring feedback all have predecessors. The central change puts motion and stiffness indices on equal footing in two transformation systems:

\[
\dot s=h(s),\quad \dot x=g_1(x,s,\gamma),\quad \dot p=g_2(p,s,\omega).
\]

Both share phase and temporal scaling; slowing motion changes stiffness according to the same progress variable. EMG provides a directly learnable index instead of starting trial and error from constant stiffness. The original explicitly says implementation “only needs minor modifications,” including replacing x with p and choosing corresponding parameters.

**Technical substance and scale.** The local computational change is small. The design work selects the represented object and synchronization relationship, connecting demonstration, generalization and execution. Unified representation makes it possible to provide appropriate stiffness when needed. The introduction acknowledges existing stiffness scheduling before identifying limitations in direct demonstration/generalization. Contributions first describe combining EMG and DMPs and then emphasize unified motion/stiffness encoding with “more importantly.”

**Evidence supporting the contribution.** Cutting and pick-and-place each include changing motion with fixed stiffness, changing stiffness with fixed motion, and changing both. With a heavier object, pick-and-place succeeds 1/7 without stiffness generalization and 6/7 with it; changing the target for the original object succeeds 7/7. Specific skill capabilities support the small structural extension, beyond fitted curves alone.

**Boundaries.** p is an EMG-derived stiffness index, not an exact recovery of all joint stiffnesses. Humans set the generalization parameters; this does not establish automatic recognition and optimal stiffness for arbitrary new objects. Success rates are small samples from selected tasks, not universal effects.

### R13: redesigning the certificate representation that limits stable learning performance

Jin et al., T-RO2023, *Learning a Flexible Neural Energy Function With a Unique Minimum…*, [original](../literature/papers/reverse_analysis_samples/R13.pdf), DOI10.1109/TRO.2023.3303011. Checked PDF2–7,9–17, including the original figure on PDF5; especially Eqs.(5),(18),(45)–(49) and Figs.10–19.

**Research question.** Stable motion learning can already guarantee trajectory convergence, but simple energy certificates force large corrections to complex demonstrated trajectories. Flexible networks may introduce extra minima that stop motion at the wrong location. The paper specifically compares quadratic forms, WSAQF, FDM, NILC and ICNN; it does not treat learning Lyapunov functions as unexplored.

**Inheritance and changes.** The route retains energy certificates, regression of the original vector field and minimum CLF correction. The central redesign is the energy structure and its trainable conditions:

\[
V=\omega^Tf(z(x))-\omega^Tf(z(0))+\alpha\|x\|^2,
\quad z(x)=(\|x\|^{1+\epsilon},\|x\|^\epsilon x^T)^T.
\]

Writing x as rh(θ), each hidden unit is constrained by \(a_{k,1}>\|a_{k,2:}\|\) and \(\omega_k>0\), making radial derivatives positive on nonzero rays. Thus the absence of additional stationary points throughout the space is expressed through finitely many parameter conditions. \(\epsilon>0\) and the quadratic term contribute to smoothness at the origin and radial unboundedness, respectively. Training distinguishes a convex problem with fixed features from a nonconvex problem jointly learning features.

Execution still uses an established minimum-correction QP to enforce \(\nabla V^T(o+u)\le-\rho\). The same radial property also helps establish compatibility between added ellipsoidal position constraints and descent. **Representation, parameter conditions, training and execution correspond to each other**, giving the design its technical substance.

**Scale and presentation.** This is a substantive technical redesign focused on representation and properties, beyond swapping networks, while retaining the broader learning/execution framework. The introduction progressively explains existing-certificate limitations; the method develops structure, learning and execution from required properties. Experiments first compare energy-representation compatibility, then final-trajectory SEA/DTWD, and finally starting points, goals, disturbances, assembly and reproduction of ultrasound trajectories in 7-dimensional joint space.

**Boundaries.** Some LASA cases are close to CLF-DM; JShape_2 and Leaf_2 are even identical, as reported. The network does not represent every function with a unique minimum; compatibility metrics do not establish quality on all out-of-domain tasks, and motion-field stability does not automatically ensure contact-force safety. This analysis learns the construction relationships without claiming independent reproof of every expressivity/stability assertion.

### R12: dividing mechanisms around the conflict between ease of manipulation and stable contact

Xing et al., T-RO2023, *Impedance Learning for Human-Guided Robots in Contact With Unknown Environments*, [original](../literature/papers/field_papers/E28.pdf), DOI10.1109/TRO.2023.3281483. Checked PDF1–6,8–13; especially Eqs.(5),(11)–(12),(19)–(21) and Figs.6–11.

**Research question.** Human-guided tools should move easily while remaining stable and accurate under environmental disturbances. Higher impedance can require more human effort. Importance comes from this conflict; an unknown environment is a condition.

**How the method forms a complete relationship.** Start with desired easy-manipulation behavior \(M\ddot x_r+B\dot x_r=f_h\), allowing human force to determine reference motion. Then design compensation/learning for environmental disturbance \(d=Kx+K_v\dot x+k_p\) to recover that behavior. Control includes \(-\hat Kx-\hat K_v\dot x-\hat k_p-K_e\dot e+g\); learning laws correspond to errors/regressors and handle parameter-error cross terms in the Lyapunov derivative.

Reference generation and impedance learning therefore have explicit roles and a shared objective. Nonrepetitive tasks use constant-coefficient adaptation; repetitive tasks exploit repeatability of parameters at the same position, using spatial learning and velocity-dependent factors for different durations, with regularization near zero velocity. The technical substance is **retaining human guidance, compensating environmental disturbance and matching learning to task repetition**.

**Scale and presentation.** Adaptation, iterative learning, admittance references and Lyapunov analysis are mature tools. The authors organize them for three-party interaction, deriving specific updates and closed-loop analysis. The contribution reduces human burden while improving task performance. Known virtual disturbances test learning/compensation; real carving/grinding tests path error and human work. H-MAN has six participants; the two Panda tasks have twelve in total, six per task. A dedicated handle also separates human force from tool-contact force, making the information interface serve the method.

**Boundaries.** Constant coefficients and spatial repetition cannot be merged into guarantees for arbitrary time-varying environments. The paper reports weak advantages under small disturbances. The velocity-based comparator also considers humans and the environment, so it cannot be said to ignore the environment entirely. Negative-semidefinite Lyapunov derivatives alone do not automatically establish that all parameters converge to their true values.

### R11: transferring human adaptation patterns into joint robot update laws

Yang et al., T-RO2011, *Human-Like Adaptation of Force and Impedance in Stable and Unstable Interactions*, [original](../literature/papers/reverse_analysis_samples/R11.pdf), DOI10.1109/TRO.2011.2158251. Checked PDF1–10; especially Eqs.(3),(14)–(16) and Figs.4–5,7–10.

**Research question.** Repeatable external forces can be compensated by feedforward, whereas unstable interactions amplify deviations. Why do humans adjust both applied force and mechanical impedance, and can robots acquire similar behavior? Human experimental patterns provide a positive design source, without first listing failures of every algorithm.

**Inheritance and changes.** Control retains dynamic compensation, feedforward, stiffness and damping feedback. Changes focus on joint learning and forgetting: feedforward updates from composite motion error; K and D update from products of that error with position and velocity errors, respectively. Forgetting increases at small errors to reduce unnecessary effort. A representative structure is

\[
\delta\tau=Q_\tau(\varepsilon-\gamma\tau),\quad
\delta K=Q_K(\varepsilon e^T-\gamma K),\quad
\delta D=Q_D(\varepsilon\dot e^T-\gamma D).
\]

The equations are concise, but each term corresponds to compensating repeatable forces, resisting unstable deviations or reducing unnecessary impedance. Human mechanisms, the control model, costs and stability analysis jointly constrain the transfer.

**Contribution and evidence.** The contribution is adaptive behavior that automatically allocates feedforward and impedance. Constant-bias, position-divergent and velocity-divergent disturbances separately test feedforward, K and D; comparison with feedforward-only learning shows the unstable-case distinction. VIA, single-joint and 7-DOF platforms further show force compensation for low-frequency disturbances, increased impedance for high-frequency disturbances and predominantly feedforward changes for added load. Such experiments explain mechanisms rather than merely rank total error.

**Boundaries.** Conditions include periodic references, dynamics compensation, noise and compensable interactions. “Optimal” in the paper cannot directly mean global optimality for arbitrary finite-horizon tasks. The lesson concerns positive transfer and behavioral roles; the 2011 comparator scope is not today's complete standard.

### R19: an additive structure can support substantive robotics research

Johannink et al., ICRA2019, *Residual Reinforcement Learning for Robot Control*, [original](../literature/papers/reverse_analysis_samples/R19.pdf). Checked PDF1–6; especially Eqs.(4)–(5), Algorithm1, Figs.2–5 and §VII. This is a leading-conference supplement, not a T-RO paper.

**Research question.** Feedback can readily handle geometric motion near a goal, but friction, small couplings and object toppling after contact make effective feedback harder to design directly. The authors retain a hand-designed controller and let RL learn corrections only:

\[
u=\pi_H(s_m)+\pi_\theta(s_m,s_o).
\]

The learner uses existing TD3. The paper explicitly does not depend on a particular RL algorithm and acknowledges concurrent independent residual policy learning and earlier neural-compensation ideas.

**Technical substance and scale.** The core structure is small. Known motion structure reduces exploration burden, while learning handles contact behavior, with training/execution correctly defined at that interface. The paper notes that dynamics couple both components, requiring optimization of total return rather than two strictly independent objectives. Baseline stability is not automatically inherited by arbitrary residuals.

**How the contribution is supported.** Practical training efficiency and contact adaptation organize the contribution. In aligned physical tasks both methods succeed 20/20; under misalignment the residual succeeds 15/20 and the hand-designed controller 2/20, with learned nudging and correction. Pure RL, environmental changes, control noise and simulation initialization explain when priors suffice and when learning matters. Reinventing TD3 is unnecessary, but the paper does not establish superiority over every strong contact-feedback method. “First control plus learning” is inaccurate.

### R21: physical modeling and mature optimization produce executable criteria

Chen et al., T-RO2025, *Stability Criterion and Stability Enhancement for a Thruster-Assisted Underwater Hexapod Robot*, [original](../literature/papers/reverse_analysis_samples/R21.pdf), DOI10.1109/TRO.2024.3492374. Checked PDF1–4,7–10,15–18, rechecked metric definitions on PDF13 and the original figure on PDF10.

**Research question.** The robot moves using both legs and thrusters while contacts alternate. How can the effects of fluid, inertia, friction, joint torque and thrust on walking become quantities that can be assessed and adjusted in practice? The paper also compares the authors' earlier flat-ground-only method, continuing an existing route.

**Inheritance and changes.** Force balance, contact feasibility and SLSQP/LP/QP are established tools. Work centers on a model combining C-shaped legs with fluid/thruster effects, friction/actuation limits and approximate alternating support, then projects these conditions into executable thrust. It first derives a vertical thrust interval, extends to a six-dimensional force/moment set and finally selects thrust within the set or projects the preceding thrust.

**Scale and organization.** This is a coordinated modeling, criterion and execution construction, beyond changing the application and without reinventing general solvers. Physical conditions become usable execution choices through explicit models/approximations. Fig.2 presents two levels; methods establish ranges; simulation separates friction/torque/mass effects; pool ramps and steps test practical capability. Ramp tests use twelve groups of five trials. The paper explains the 37kg physical robot versus a principal 5kg simulation setup. It also acknowledges that suitable constant thrust can climb steps; reducing trial and error in selecting constant thrust is part of the value.

**Boundaries of independent judgment.** “Stability probability” is the fraction of counts judged stable during crawling, not successful trials divided by trials. The approximately 50% on PDF17 does not mean half the climbs succeeded or provide an independent full-system stability certificate. Eqs.(34)–(35) on PDF10 construct a convex hull from extrema in finitely many directions. Under the stated convex feasible-set setting, this is an inner approximation; being outside does not automatically mean true infeasibility. Eq.(37) minimizes the sum of squared distances to vertices, whose solution is the vertex mean; this is generally not equivalent to maximizing minimum distance to the boundary.

The existing [general mathematical counterexamples](../literature/top_venue_sample_math_checks.json) remain. This round rechecked original equations without extending general counterexamples into rejection of all experiments. The physical → computational → execution relationship can be learned while preserving approximation direction and claim strength accurately.

## 3. Four additional cases: specific supplements rather than one imposed formation path

| Sample and pages checked this round | Actual design and changes | Technical substance and contribution boundaries |
|---|---|---|
| R05, T-RO2018, [original](../literature/papers/reverse_analysis_samples/R05.pdf), PDF1–6; Eqs.(1)–(2),(6)–(19),(23)–(33), TableII | A relative Jacobian describes one arm holding and another operating. Force/position estimates produce impedance, including changing stiffness, and references, followed by constrained adaptive control and a force channel | Estimation actually changes the tracked reference; modules serve the same cooperative task. Tracing with five participants against conventional impedance supports overall capability. Force residual D is not an independent measurement of true intention; nonsingularity, initial conditions and human-model assumptions remain |
| R03, TIE2018, [original](../literature/papers/reverse_analysis_samples/R03.pdf), PDF1–4,10; Eqs.(10)–(13) | Filtered history constructs D,T so computable residual P=T−Dθ̂ carries parameter-error information. Newton–Euler parameterization is extended by finite-time updates and SVD reduction | Faster convergence requires specific information/update design rather than tuning learning rates. The paper explicitly draws on error-driven identification. SVD removes redundancy but does not create excitation; kinematic IE and dynamic PE/reduction conditions remain distinct. Not all convergence arguments were reproved here |
| R17, TNNLS2021, [original](../literature/papers/reverse_analysis_samples/R17.pdf), PDF1–5,10; Eqs.(8)–(11),(18)–(20) | BLF control addresses saturation through feedback defined on the whole axis, with bounded derivative and strongest near boundaries; a changed saturation approximation lets auxiliary states decay without saturation; RBFNN is simplified | Functions, auxiliary dynamics and approximation are coordinated around execution capability. Error limits are explicitly soft constraints; infeasible cases concern boundedness/recovery, tested physically on PDF10. This is not an always-within-bounds guarantee; assumptions such as known inertia remain |
| R15, Mechatronics2020, [original](../literature/papers/reverse_analysis_samples/R15.pdf), PDF2–5,8,11–12; Eqs.(19),(23)–(24),(30)–(36) | K/D and references update iteratively. Minimal null-space coordinates and an inertia-weighted inverse satisfy J#ᵀMZ=0, dividing task tracking from joint compliance under pushes, with task-channel disturbance compensation | Combination is realized through dynamic compatibility. Twenty learning rounds of Sawyer/cardboard contact and body-push experiments serve distinct capabilities. References also change, so reduced force cannot all be attributed to K/D. This is neither a leading-journal identity nor an arbitrary-time-varying-disturbance guarantee |

## 4. Contribution scale and paper organization learned from these designs

Both formation routes proposed by the user have concrete support: R01, R13 and R17 modify key structures along mature routes; R05, R12 and R15 organize mechanisms around tasks; R11 transfers human adaptation patterns into joint updates. R19 can mean adding a residual to feedback or reallocating control/learning roles. Classification must not constrain design.

Technical substance appears in different relationships: R01's shared phase, R13's radial properties/trainable conditions, R12's human-force-reference/environment-compensation division, R15's dynamic compatibility and R21's executable thrust set. These relationships let components change capabilities. Some are concise, others require extensive derivation; no universal new-component count, algorithm-change size or improvement percentage follows.

Paper presentation matters, and specific writing shows how:

- **R01 expresses a small change at the appropriate capability level.** The equation-level change is two transformation systems; the contribution is directly demonstrable, synchronously generalizable impedance skills; experiments separate motion/stiffness. Merely saying “apply DMP to p” understates its value.
- **R12 lets the conflict determine narrative and metrics.** The introduction explains effort versus contact performance, Eq.(5) gives target behavior, methods assign roles and experiments measure both error and human work, making the design's importance traceable.
- **R13 supports related contributions with one structure.** It first explains certificate-induced trajectory limitations, then derives the structure. Representation compatibility, convergence, position limits and generalization surround the same design. It does not require large wins in every case.
- **R11 and R19 explain mechanisms through behavioral comparisons.** Constant/divergent disturbances and aligned/misaligned tasks show when each adjustment is needed, giving algorithm rankings concrete meaning.
- **R21 organizes multiphysics analysis into execution choices.** Moving from one to six dimensions makes the coordinated work understandable; actual tasks test its role, while approximation/metric limits remain.

Such presentation clarifies the level, meaning and evidence of actual improvements and can substantially strengthen persuasion. It neither requires portraying small changes as a new paradigm nor lets publication status erase evidence gaps. Limited comparators or questionable theorems in an old paper cannot automatically become universal standards today.

## 5. Comparison with project failures: confirmed errors and unresolved causes

Sources are the [failure record](../project_records/research_failures_and_corrections.md), retained experiments, and old decisions/posterior-planning reports at Git `db28365f33dc314ce8bfe83a23796dd18cf27f5b`. This round reran no control experiments and does not assign one cause to all method failures.

| Actual project situation | Specific paper comparison | Supported explanation and prevention |
|---|---|---|
| Old records repeatedly correct mandatory NEOC, algorithm exclusivity and requiring comprehensive exclusion of repairs first | R01 explicitly uses small modifications; R19 reuses TD3; R21 uses mature optimizers | Extra gates exclude normal method sources. Study how existing methods can be adapted/improved. Do not invent a history in which some already excellent specific method was rejected for this reason |
| Premature shift to Q3 followed by drift toward candidate A; old records contain withdrawal and user correction | R12 stays focused on human burden/contact performance; R13 on accuracy/convergence | Using local algorithm comparisons to shift the question's focus is a documented level error. Strong prior work calls for deeper inquiry into worthwhile improvements within the original question, not automatic narrowing to a task where one tool can win |
| Added policy rolling and posterior planning lack stable advantages over adequate simple methods | R12 divides roles from target dynamics; R15 realizes combinations through decoupling | Adaptation and prototypes are legitimate and informative, but these instances have not established how connections improve important task tradeoffs. Workflow completeness does not form a paper method; this also does not prohibit natural fusion |
| Short-history direct PPO omitted complete MHE/preview; small CEM formulas also inherited candidate 2's name | R01 clearly defines demonstration → representation → generalization → execution; R19 defines a fixed controller/residual training object | Retaining the complete candidate name after changing the experimental object accumulates mistaken judgments. This is a validation-object mismatch, not falsification of original candidate 2. Preserve specific implementation names and omitted conditions |
| Two histories produce up to about 0.2806N future-force difference, yet the same K/D has only about 0.4334% mean loss over two cases relative to each case's best candidate | R11 uses bias, position-divergent and velocity-divergent disturbances to test different action roles | Information differs, but important control improvement remains unestablished; the knowledge-to-action-value link lacks evidence. These two cases neither refute history information nor establish an important impedance contribution |
| Joint scanning outputs have small positive results; a new gains-only cohort does not consistently beat fixed gains; corrected-model joint output improves about 0.201% over analytic control with an interval spanning zero | R01 separates motion/stiffness effects; R19's baseline suffices when aligned, while learning matters under misalignment | Explain whether benefit comes from references, feedforward, gains or their coordination. Small numbers alone do not decide contribution; a clear capability improvement like the samples has not yet appeared |
| Nodal roots, policy quality and complete guarantees were confused; counting prediction only or waiting for old plans affects timeliness conclusions | R13 separates energy compatibility/final trajectories; R21 requires distinguishing feasible approximations/true stability | Numerical objects and execution interfaces must match claims. Later release-aligned/background improvements in constrained-contact tests yield about 0.67%–0.78% local benefit; weaker execution cannot represent all online methods |

Posterior planning adds direct evidence: across 64 new cases, additional mean paired improvement over a frozen planner already able to update knowledge is 0.08363%; median planning is about 668ms plus about 162ms construction, while reported performance omits waiting. This supports the conclusion that this instance cannot yet serve as the main rapid-control method. No evidence attributes all small benefits to absence of a network, an overly simple model or poor understanding of leading journals. These require specific analysis rather than choosing any explanation to protect the proposal.

Failure is understood at two levels: **judgment departed from normal research design illustrated by the samples**, through excessive defensiveness, narrowing and object confusion; **the specific prototypes also still lack sufficient task progress**, which lower standards or renaming cannot remove. Correct design/judgment at the first level; substantive method improvements and evidence remain necessary at the second. These local instances have not decided complete B0 as a whole.

## 6. Subsequent confirmation and implementation in requirements

Learnable design actions have original sources: follow a strong route to ask how representation/conditions limit capability as in R13; allocate roles from required behavior as in R12; permit small changes and identify the added skill as in R01; translate mature-component combinations into identifiable task behavior as in R19. Subsequent project work can use modification, transfer and natural fusion without prescribing one formation path.

In subsequent discussion the user further confirmed that the goal is to learn how to construct methods and what leading-journal methods actually look like. A finished paper's narrative of methods arising naturally from the problem does not reveal the real discovery sequence. Method formation may improve a mature route or naturally combine/adapt mechanisms around a problem; problems and methods can co-evolve through adaptation, derivation, implementation and verification. Improvements may be small; technically sound adaptation and actual progress confer research value. Organizing real work into clear, persuasive contributions is an important part of paper preparation.

The user then explicitly required these understandings and reverse-analysis lessons to guide practices and judgments at each research stage; they are now in the [long-term requirements](../project_records/long_term_requirements.md). This report retains examples, source locations and analysis boundaries; the requirements retain confirmed stable understanding. Sample patterns do not become a uniform acceptance checklist. This update selects no new main method and adds no control experiment. Continuation follows the [handoff](../project_records/HANDOFF.md) and latest user task.
