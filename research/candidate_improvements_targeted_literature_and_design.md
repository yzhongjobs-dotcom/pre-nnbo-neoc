# Targeted literature,design and quick tests for both candidates

2026-10-19.The user authorizes search,design,rapid tests and iteration within the original question. This follows [candidate definitions](optimal_impedance_question_and_methods.md),without changing requirements or adding journal-acceptance gates. Revise conclusions with evidence.

## 1. Design-specific search and reading

Searches on2026-10-19 include MHE variable impedance control,moving parameters optimization adaptive impedance,Sobolev approximate nonlinear MPC,impedance guidance learning,model predictive impedance gradient neural,and specific titles/author sources. Ranking/nonretrieval does not prove absence. Abstract screening is not full-text auditing.

|Paper/source|Reading scope|Established capability,adaptation and implications|
|---|---|---|
|Roveda et al.,*Q-Learning-based model predictive variable impedance control for physical human-robot collaboration*,Artificial Intelligence312(2022)103771,[author PDF](https://ipg.idsia.ch/preprints/Roveda2022a.pdf),[DOI](https://doi.org/10.1016/j.artint.2022.103771)|Introduction,control/network/update passages,especially§2.6–2.8,Algorithm3,experiments;not all stability proofs|Predictive control/actor–critic already combine,with reference/damping actions and periodic model/critic/actor updates. Learn problem→optimization→feedback approximation,not novelty from replacing optimization with networks. Different tasks,outputs/timing prevent abstract-only claims of missing contact capabilities.|
|Lu et al.,*You Only Plan Once: A Learning-Based One-Stage Planner With Guidance Learning*,RA-L2024,[local paper](../literature/papers/method_transfer/YOPO.pdf),[author PDF](https://tju-air-lab.github.io/projects/YOPO/YOPO.pdf),[code](https://github.com/TJU-Aerial-Robotics/YOPO)|Objectives,numerical gradients,privileged maps/deployment observations;repository/loss entry checked,not full reproduction|Trajectory optimization/primitives guide cost-gradient learning and candidate/scoring generation. Transfer task-consequence training rather than gain-label fitting,with gradients through contact/actual gains. Do not copy flight geometry or call it PPO.|
|Vahidi-Moghaddam et al.,*Extended Neighboring Extremal Optimal Control with State and Preview Perturbations*,[arXiv2306.04830v1](https://arxiv.org/html/2306.04830v1)|§III-B/C,(11)–(36),Algorithm1,constraint singularities;not full proofs or journal status|Two feedback matrices extend neighboring extremals to state/preview disturbances and imperfect nominal solutions/larger changes. Candidate1 may borrow this;parameter/preview sensitivity is not exclusive. Check positivity,active changes/state singularities. Not yet implemented here as contact control.|
|Zhang et al.,*Observer-Based Variable Impedance Control Using Moving Horizon Estimation for Robot Machining Thin-Walled Workpieces*,TIE71(6):5972–5982(2024),[publisher](https://doi.org/10.1109/TIE.2023.3292856)|Publication identity/searchable abstract;no readable full text obtained|MHE-driven VIC directly exists;the combination is not new. Missing full text prevents attributing absent capabilities or trusting automated secondary summaries.|
|Wang et al.,*Moving Parameters Optimization-Based Model Predictive Adaptive Impedance Control*,[publisher](https://doi.org/10.1109/TIE.2025.3587119)|Identity/index abstract,not full text|Grinding-related parameter optimization/predictive adaptive impedance is close prior. Obtain methods before delimiting rapid-update differences.|
|Lüken et al.,*Sobolev Training for Data-efficient Approximate Nonlinear MPC*,IFAC2023,[publisher](https://www.sciencedirect.com/science/article/pii/S2405896323009126)|Publication/abstract index;direct access restricted,derivations unread|Sensitivity-based derivative supervision is a traceable cross-field lead. Do not transfer unread derivations/stability claims or impose it as a new mandatory training module.|

## 2. Form methods rather than rewriting the story

Start with shared estimation/actual execution. Build a comparable warm-start predictive route and actual task-gradient candidate2 network. Train against shared estimation error/actual gain transitions,not ideal labels/instant jumps. Retain legitimate reference/force compensation and assess it separately from pure K/D.

Borrow mature MHE,contact,rolling optimization,neural feedback/energy constraints. The specific relationship is **delayed finite-footprint information→local joint estimation→execution/estimate-conditioned decisions→force/motion costs after actual transitions**. This is testable adaptation,not established leading-journal contribution. Combination does not prove value,nor mature components negate potential.

See [protocol](../experiments/mhe_guided_impedance_pilot/design_and_protocol.md). At this research checkpoint B0 core/frozen files were unchanged;ordinary shooting is not HJB/NEOC,and reduced estimation is not full B0 MHE.

## 3. Iteration evidence

Two training pilots connect reduced estimation/execution:current Fn RMSE≈0.023–0.026N,MHE≈3ms,shooting≈31ms locally with zero release wait,not certification. One model case improves,the other worsens. Among28 correlated states,25%gradient-direction cosines are negative and35.7%of local model optimizations increase true short-time costs. This diagnoses this prototype,not an established field problem.

Tune fixed control on12 training cases,then collect24 fixed/analytical paired trajectories. Training started but stopped after user correction,without completed actor,validation or independent tests. These cannot judge complete candidate2 or journal contribution.

## 4. User correction:main comparisons must use field methods

The user identified weak comparison level:at least original NNBO/related methods were needed rather than tuned fixed impedance as main innovation evidence. **Correction accepted;training paused.** Fixed control remains engineering/ablation support,while existing variable/optimal impedance is the contribution comparator. Better than fixed is not progress over existing VIC. Implement strong field methods before more fixed tuning/training.

The recorded capability-based implementation order follows;these are planned,not completed comparisons:

|Main comparator|Original capability/check entry|Required comparison and safeguards|
|---|---|---|
|**Kong NNBO,TSMC2025**,[paper](../literature/papers/NNBO.pdf),[minimal reproduction](../experiments/nnbo_source_reproduction/README.md)|HC identification,single critic,Newton Downhill references;(22)position/force-square tradeoff,Remark6 zero desired force|Check original identification/robot reference execution before scan adaptation. Retain position-reference outputs,not forced four K/D. Derive/record nonzero-force/moving-contact adaptations;target mismatch is not innovation. Minimal reproduction has no inner position loop.|
|**Jin predictive VIC,TMECH2023,E41**,[paper](../literature/papers/field_papers/E41.pdf),DOI10.1109/TMECH.2022.3204350|Gain controls,OS QP,MS prediction,energy correction,precision/compliance|Read/implement representative gain optimization/interfaces/constraints and compare loops/full cost;arbitrary shooting scripts are not E41 reproductions.|
|**Li unknown-terrain impedance learning,T-RO2025,E27**,[paper](../literature/papers/field_papers/E27.pdf)|Randomized-terrain prelearning,VIC feedback,sim-to-real;full text local|Direct candidate2 neighbor:verify outputs,feedback,training/transfer before implementation. Old short-history PPO or reused names cannot substitute.|

Zhao TIE2022(E38,DOI10.1109/TIE.2021.3134082)and Roveda2022 also inform model-learning design/comparison. Select by task fit/achievable scope,not method counts.

Validate original algorithms/interfaces first,then explicit shared-scan adaptations. State task,tuning budget,training/identification,sensing,action freedom,execution constraints and compute/release cost. Preserve native configurations and use information/action alignment to identify benefit sources;removing original advantages is not fairness. Until implementations/checks are complete,these are engineering prechecks,not proven NNBO/strong-VIC superiority or journal readiness.
