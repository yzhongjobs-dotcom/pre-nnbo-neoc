# Handoff: rapid impedance updates under environmental change

Research status updated: 2026-10-22. Check the current user request and Git state first; historical plans are not today's instructions.

**Current numerical entry, superseding the October 19 status below:** the October 21
[36-paper setup analysis](../research/simulation_setup_reverse_analysis.md) led to
the implemented [cobot benchmark](../experiments/cobot_scan_benchmark/README.md),
paper-method adaptations, two candidates and Candidate 1 v2. On October 22,
[physical-outcome reporting](../experiments/cobot_scan_benchmark/physical_results/results.md)
replayed all frozen deployable baseline variants and candidates on the same two
main tasks, keeping seed-3000 and seed-4000 cohorts separate. It contains 3,840
episode records; neural training seeds do not increase the 40 independent scenes
per task/cohort. The tables separate force RMSE, maximum error, peak force,
no-contact duration and update timing. Parameters and original results were not
changed. Replay differences and one failed low-damping run are explicit.

The new candidates change representation and interface: Candidate 1 is local
stiffness/surface-motion estimation plus direct OCP and feedforward, and Candidate
2 uses history/lock-in estimation and a distilled network. They do not implement
the complete original B0/HJB/spatial-preview design. Candidate 1 v2 and Candidate 2
improve main-task mean errors, but not every physical outcome; the v1 seed-4000
runaway, the fast-push limitations and the 4 mm respiration negative remain.
Xue's disabled stiffness branch and unsuccessful feasibility search, and limited
Jin/learning tuning, remain comparison limitations. The latest user request was
to produce physical-outcome tables with missing quantities computed, not to start
another method redesign. All replay work stops at its completed checkpoint.

The following sections preserve the earlier design and evidence context. Their
old "latest" experiment references describe the October 19 checkpoint, not the
current numerical entry above. Scientific asset consolidation is recoverable at
`368a29539a723f960f27fb7e28362a65f51428c8`.


## 1. Understand these six points first

1. **Goal:** one paper with a substantive scientific contribution and evidence strong enough to compete at a leading journal. This has not yet been achieved.
2. **Core question:** when the environment changes during continuous robotic scanning of soft tissue, how can new contact information rapidly update near-optimal impedance and improve actual contact-force and motion behavior? This includes environment understanding, impedance decisions and timely execution.
3. **Candidates:** the shared front end is MHE with spatial preview. Candidate 1 explores NNBO/critic–HJB feedback and fast updates such as NEOC. Candidate 2 draws on YOPO to train a network offline in simulation and produce impedance quickly at deployment. Both can be improved, combined or retired.
4. **Actual progress:** both candidates have been implemented and run under the research question: two-region Hertz contact/joint rolling estimation, four actual gains and energy, positive-definite gain/rate costs, terminal and trajectory checks. Candidate 1 is local quadratic HJB feedback with directional neighboring Riccati updates; candidate 2 is a task-gradient neural network with three training seeds. Work includes 24 main scenarios, mechanism ablations, eight distribution-shift scenarios and independent compute/release checks. These tests have not established a journal-level contribution; they neither implement the complete B0 fixed-node-root method nor constitute physical-robot validation.
5. **Largest deviations:** excessive defensiveness, retreat/narrowing around strong prior work and misunderstanding candidate 2. An earlier reduced carrier omitted objectives/constraints; the latest implementation must also avoid concealing physical performance deficits with controller-preference costs. Specific failures are recorded. Do not resume mistaken judgments or indefinitely shield the implementation as “only local.”
6. **Current stopping point:** the user required experiments to follow the research question/method: implement and simulate each candidate first, then fairly compare suitable baselines if any exist, without forcing unsuitable comparisons. This order was followed. Jin's adapted multistep variable impedance scan controller has lower force/tangential errors. Corrected ablation gives a small spatial-information benefit for candidate 1; candidate 2's interval crosses zero. NEOC on/off benefit is near zero and updates are slower. Next improve the relation between physical performance and cost objectives, update efficiency, and the front-end-to-action pathway. Do not manufacture wins through gain-preference terms or add baselines without purpose.

**Latest task and primary objective.** The user requests methods that solve the original question with leading-journal innovation and evidence under the confirmed workflow. Implement/test each candidate separately and improve its own deficits. Experiments follow the question/method, followed by suitable comparisons; none need be forced. These requirements are in the [long-term requirements](long_term_requirements.md). Do not return to fixed control as the main comparison, repeatedly audit the standards or silently change the question.

**Latest entry: [dual-candidate task simulations](../experiments/dual_candidate_task_simulation/README.md) and [results](../experiments/dual_candidate_task_simulation/results.md).** Candidate 1 changes the value representation and cannot be called the original B0 root or inherit its certificate. E41 is a source-equation reconstruction adapted to MS scanning, with a different numerical optimizer, not a reproduction of the original Panda experiments. Common physical scoring removes controller-family gain-preference/rate penalties and reports component errors. Use `compiled_timed_*` for independent compiled timing; old `timed_*` files diagnose software overhead only. Failed NNBO adaptations remain isolated. The original question and the candidates' computational division of labor remain.

## 2. Information priority within this project

| Level | How to use it |
|---|---|
| Current user task and confirmed corrections | Establish this turn's scope first; the latest explicit user request overrides old arrangements. Briefly restate understanding and scope for each new request. |
| [Long-term requirements](long_term_requirements.md) and [AGENTS](../AGENTS.md) | Confirmed goals, stage practices/judgment, method-formation understanding and file/synchronization constraints. The research checkpoint updated requirements as explicitly requested; it did not change AGENTS. |
| Current research files and HANDOFF | Candidate definitions, actual evidence and stopping point; all remain checkable and correctable. Being committed does not make a statement true. |
| Older reports, decisions and reference workflows | Consult for explanations or resource recovery. Old “next steps,” withdrawn topics and early gates do not automatically become instructions. |

This is the order for project materials, not a change to system/developer instruction priority. Scientific judgment must return to the relevant source version, equations, protocol, code and results. Publication status, agent synthesis and file recency cannot establish scientific correctness by themselves.

## 3. Read in layers, not by rereading everything

| Layer | Materials | What the successor should be able to explain |
|---|---|---|
| **1: goals and avoiding drift** | This file, current request, [requirements](long_term_requirements.md), [failure review](research_failures_and_corrections.md) | What the user wants, what to do now and which mistakes must not be continued |
| **2: research and judgment foundations** | [Question](../research/research_question.md), [field state](../research/field_state.md), [candidates](../research/optimal_impedance_question_and_methods.md), [independent paper analysis](../research/independent_reverse_analysis_of_top_journal_methods.md) | Existing capabilities, candidate inputs/training/outputs, verified scope, how sample papers form methods, and where the analysis needs calibration |
| **3: evidence for the current decision** | [Literature](../literature/literature_index.md), [experiments](../experiments/experiment_index.md), [transferable mechanisms](../research/transferable_mechanisms.md) and their original papers/data/code | Grounds for key judgments; distinguish source conclusions, agent synthesis, hypotheses and measurements; know what is reusable and missing |

For B0 models/theory, consult the relevant sections of the [method](../research/method_design.md), [appendix](../research/technical_appendix.txt) and [conditional proofs](../research/theoretical_proofs.txt). Starting work does not require reading all history, all 65 papers or rerunning every experiment.

## 4. Research definitions: prevent another misunderstanding

**Shared front-end concept:** contact force/moment and motion observations → MHE and finite-footprint spatial preview → environment/contact-response estimates → impedance decisions → actual execution. Complete noisy MHE and usable preview control benefits remain unverified.

| Object | Correct definition and key distinction |
|---|---|
| B0 | Revisable complete written starting point: cylindrical contact → MHE → descriptor η → finite-horizon HJB for four actual gain rates → fixed-node root → ANEOC and conditional release of control updates. Written conditional proofs are not implemented robot guarantees. |
| Candidate 1 | Build an optimal-control problem from new environmental information, obtain feedback through critic/HJB ideas and update quickly using NEOC or other mechanisms. Original NNBO generates position references using HC pre-identification, critic and Newton. Four K/D gains and gain-rate control are B0 extensions, not the original paper interface. |
| Candidate 2 | Retain estimation/spatial preview and draw on RA-L2024 YOPO to move decision burden into simulation/offline training. PPO is optional; original YOPO uses guidance learning, not PPO. Adaptation/fusion remains a design question; do not freeze the architecture or training algorithm. |

“Optimal” is relative to an explicit task, model, information and constraints. Candidate outputs, node roots and local numerical references are not automatically globally optimal. Earlier stochastic-expectation, joint feedforward/reference actions and simplified contact extensions do not automatically inherit B0 proofs.

## 5. Current evidence and reusable assets

| Material/fact | Use and boundaries |
|---|---|
| NNBO reproduction, finite-horizon/slow variation and dynamic strong controls | Check critic, ordinary Newton/chord, NEOC, errors and cost. These tests neither establish a general NEOC advantage nor invalidate B0 as a whole. |
| Known-model networks/RBF and constrained contact updates | Reuse mappings, actual gain/energy execution, adjoint gradients, task-cost learning and release alignment. Networks do not reliably beat RBF; background model improvement has local positive results, so old delay negatives alone cannot reject online methods. |
| Within-scan prototype | Local estimation with model prediction or small CEM feedback formulas. The initial 48-case joint-output cohort improves roughly 0.433% over the analytic comparator; later 24-case pure-K/D and separate 24-case corrected-physics cohorts show no stable advantage over adequate simple controls. Not complete MHE/HJB/NEOC or PPO/NN; do not merge the three cohorts. |
| Old direct PPO | Omits complete MHE/spatial preview and does not validate original candidate 2. Unneeded training pipelines were removed; conclusions and mistakes remain in the failure review. Do not restore “complete candidate 2 failed.” |
| Candidate A, policy-library rollout and posterior planning | Removed from the workspace/default main line; useful results, sources and lessons were extracted. Independent tools remain in [spatial contact estimation](../experiments/spatial_contact_estimation_audit/README.md); identification is not closed-loop benefit. |
| Local literature | 65 deduplicated PDFs/texts and 68 usage records, including relevant material from the two requested repositories, supplemental originals and YOPO. Abstract, passage, full-text and reproduction scopes are separate; downloading does not expand reading. |
| Paper-by-paper reverse analysis | Earlier screening covered 22 identities/abstracts. Ten original methods/evidence sections were reread, including all five T-RO papers: six detailed cases and four supplements, then comparison with project failures. Not 22 full-text or complete-proof audits. Confirmed lessons are in the requirements. |
| [Earlier reduced MHE precheck, paused at that time](../experiments/mhe_guided_impedance_pilot/README.md) | Two pilots and old trajectories retained. The pause was historical. Latest work is in separate-candidate directories; do not restore a present stop instruction or claim the NN remains untrained. |
| [Previous separate design/tests](../experiments/separate_candidate_validation/README.md) | Model/learning iteration, rate-action/100Hz feedback alignment, 72 training scenarios and 10800 correlated contexts, 24 independent scenarios repeated at three tradeoffs. Model reuse does not stably beat ordinary optimization; learning has positive results relative to the declared E41 adaptation, but information effects depend on the tradeoff. Not complete B0, original-paper reproduction or journal-level judgment. |
| [Current dual-candidate task simulation](../experiments/dual_candidate_task_simulation/README.md) | Original question → separate implementations/experiments → Jin MS selection. New 48 training scenarios/3840 contexts, 24 independent scenarios in four categories, three learning seeds, eight distribution-shift and eight independently timed cases. Hertz/objective relationships restored; synthetic 2DOF model, CV/fixed-bias MHE instance and local value feedback rather than B0 root. Comparator force/tangential behavior is better; the candidates' normal-position errors are slightly smaller. No completion claim. |

**Overall: the original question remains; a specific important advance and complete method advantage are not yet established.** Platform/control interface, material/sensor calibration and references remain to be verified. These are candidate concerns, not selected new topics or fixed gates.

## 6. What to do next, and to what extent

The table records the agreed continuation order. **A new user task takes precedence; standards calibration is not an approval gate to repeat every turn.**

| Priority | Work/deliverable | Completion/boundary |
|---|---|---|
| **P0: requirements implementation complete** | Confirmed method-formation understanding and paper-analysis lessons are in the requirements; source evidence remains in the analysis. | Continue at the appropriate stage without repeating understanding audits. Actual design/evidence must demonstrate implementation; a completion record is not method success. |
| **P1: improve each candidate's actual task effect** | Hertz/Qg/Rg/terminal/trajectory checks were added while retaining separate mechanisms. Determine why preview helps yet force/tangential performance trails; improve information-to-impedance relationships and update efficiency, and verify platform/task scale. A low gain-preference cost is not a physical win. | Comparisons follow the question/method; do not impose a fixed MPC-based pass/fail gate or force an unsuitable NNBO adaptation. Deepen Jin/E27/NNBO where useful and reproduce/adapt fairly only if suitable. Explicitly state representation changes and evidence limits; do not retreat from the question because strong prior work exists. |
| **P2: iterate and converge from evidence** | Decide whether implementation, carrier, mechanism or question needs change; progressively form the paper's contribution and complete evidence. | Do not overgeneralize local evidence or shield the idea indefinitely. Explain substantive question changes and their effects. Judge actual contribution/evidence, not workflow completion or wording. |

NN weights, training data and grouped results are available for inspection. They belong to reduced implementations, which are not automatically selected as final methods. PPO has not been selected; candidate A/future posterior planning has not been restored. Platform, interface, frequency and task scale were requested and remain pending; silence does not mean no platform exists. Complete-method formalization and public-original-paper work can continue without inventing clinical thresholds. Stop all training/experiments at completed checkpoints; do not claim background work continues.

## 7. Recovery, reproduction and update rules

- Preserve `00_基线快照_2026-10-14/` byte-for-byte, including names. The three working B0 core files may be revised with research evidence and explicit effects on assumptions/proofs. The research checkpoint added independent precheck code/data without changing B0, older experiments or the snapshot.
- Paper/file check: `python -X utf8 literature/verify_materials.py`. Select experiments by their reports and Git versions; do not silently refresh historical result hashes or overwrite original data. File checks are not scientific reruns or physical validation. Consult Git at the reported version for historical paths.
- Consolidation checkpoint: `368a295`; preceding materials: `db28365f33dc314ce8bfe83a23796dd18cf27f5b`; earlier complete materials: `b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c`. Use `git show <revision>:<historical-path>` for necessary individual files; do not restore withdrawn branches wholesale.
- Reports must explain understanding of the question/candidates, actual work and results, comparator choice, relevant numbers and limits. File counts, “24 more cases,” abstract plans or commit IDs alone are insufficient.
- Git history tracks changes. Keep scientific reasoning with the corresponding method, experiment or failure analysis; do not maintain separate operational change logs. After checking a completed checkpoint, automatically synchronize using `python -X utf8 scripts/research_sync.py --message "research: specific change and purpose"` under AGENTS, without repeated upload approval.

## 8. Concise first handoff report

Explain in this order: **goal/core requirements → research question/two candidates → completed/pending work and evidence boundaries → largest past deviations → what to do now and why.** Point each item to its source. State unread material and unconfirmed judgments honestly; do not claim every file has been read.
