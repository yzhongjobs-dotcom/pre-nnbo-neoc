# Handoff: rapid impedance updates under environmental change

Research status updated: 2026-10-22. Check the current user request and Git state first; historical plans are not today's instructions.

**Paper-style overview (Chinese):** `research/paper_draft_zh.md` lays out field status, problem, C1, C2, contributions and results in journal structure, with 【待补】 placeholders for missing baselines, fixes and theory. Keep it in sync when results change.

**Confirmed by the user (2026-10-22): refined question, contribution framing and baseline set.**
The refined question is in `research/research_question.md`; contributions, candidate fixes
and the baseline set are in `research/optimal_impedance_question_and_methods.md` §2.6.
The current final Chinese section reflects the old baseline set and will be revised after
the plan below. Ordered plan for the next (local) session; state the quota cost before each
compute-heavy step and run nothing large without the user's go-ahead:
1. Short derivation: identifiability of stiffness vs equivalent drift under constant-force
   regulation with/without dither; why a lumped disturbance observer cannot schedule gains.
2. Implement and tune the three new baselines (DOSMAC-style DOB + sliding-mode admittance,
   Cao-style vision compensation with noisy/delayed true surface motion, Wang TIE2026 G02)
   on the existing train/dev protocol.
3. Improve C1: T2 stiff push/retract and fallback bounce; replace the tuned B1 fallback with
   a rule-based one; support or drop "near-optimal".
4. Improve C2: drift estimate as input plus a structured reference output.
5. Stability/passivity condition for the variable gains.
6. Re-run the seed-7200 cohort for the final method set and revise the Chinese section.

Background of the decision (the analysis the user accepted):
Novelty assessment, non-defensive: the idea survives at the problem/mechanism level, not as
"first joint motion-stiffness estimation". Socratic reading of strong priors (inferences marked):
Cao G13 needs external vision/markers and does not handle stiffness heterogeneity; DOSMAC lumps
both changes into one disturbance, so it cannot schedule gains by stiffness (inference: the same
motion gives ~30x force disturbance on bone, so one observer bandwidth cannot suit both);
Xue/Wang G02 attribute motion-induced force error to stiffness (supported by our adaptations:
C1 without drift 1.502 N, Xue T2-stiff 3.414 N); CARE/Chalasani rely on static palpation with
deliberate excitation, while constant-force regulation removes the excitation needed for
identification (C1 v3 log: k unidentified without dither) - a tracking-vs-identification
conflict to be checked theoretically. Packaging proposal: problem = force-only, vision-free
constant-force scanning where motion and stiffness are confounded; method = joint estimate of
local stiffness and equivalent contact drift, source-based allocation (stiffness -> D/K_PF,
drift -> u_ff) in one short-horizon OCP, small dither with an identifiability condition
(transfer of persistent excitation/dual control). C2: possible transfer of RMA-style
privileged teacher-student adaptation (in-field precedent unchecked); likely needs the drift
channel and a structured reference output. Missing baselines, by priority: DOSMAC (standard
DOB + sliding-mode admittance if full text unavailable), Cao-style vision compensation using
noisy/delayed true surface motion, Wang TIE2026 G02, and full Haninger GP (or label the current
row as a damping-rule adaptation). Proposed main set: two fixed, Jin, Xue, DOSMAC, Cao-style,
Wang; EtR-A as ablation; Li and Haninger rule to appendix. User to decide: refine the question
wording (and the K=0 vs "stiffness" wording) and the baseline set.

**Latest completed checkpoint (2026-10-22): original-source field review for the user's motion/stiffness focus.**
The user asked to select relevant papers from both local libraries and update the
field assessment, then explicitly rejected abstract reframing and extra research
questions. [Field state](../research/field_state.md) was rewritten in Chinese as a
related-work analysis of adaptive/optimal force control, online predictive models,
learned impedance, respiratory compensation, joint geometry/compliance estimation
and stable execution. [Research question](../research/research_question.md) preserves
the remote's parent question/current-implementation descriptions and records the
user's force/position-based motion-and-stiffness formulation as this review's focus.
Neither candidate was deleted, and no new method was selected.

Twenty-two selected originals/author manuscripts and texts were added: 87 unique
PDFs and 90 usage records now. [Review sources](../literature/motion_stiffness_review_sources.json)
and the unified manifest retain identities, roles, hashes, versions and actual
reading pages. Wang TIE2026 (G02), Roveda AI2022 (G08), VICES (G24) and CARE (G21)
are now local. The two user-named breathing papers are G13 (Cao TIM2025:
vision-predicted reference plus force admittance) and G14 (Velikova ICRA2024:
constant-force acquisition and INR volume reconstruction). They are different
contributions. Zhang TIE2024 MHE-VIC and DOSMAC remain abstract/excerpt-level
sources; their equations and some joint-estimation papers still need full text.

Scientific conclusion: do not claim all VIC attributes force changes to stiffness,
assumes stationary/known surfaces, or fails when motion and stiffness vary together.
CARE already analyzes joint compliance/registration observability; Haninger and
learned approaches already coordinate reference/impedance; DOSMAC explicitly
handles both changes as a lumped disturbance. These precedents do not establish
that the exact new candidate is fully covered, but the current broad novelty claim
is not established. C1's implemented drift is an equivalent force-prediction state
that also absorbs material/relaxation mismatch, not demonstrated recovery of true
surface position/velocity. Ablation numbers establish the channel's utility, not
physical source separation. The benchmark Haninger row is a damping-rule task
adaptation, not its original full GP predictive controller.

The user wants concrete scholarly related-work prose, faithful to original papers;
do not resume the rejected abstract reframing or protect C1 through a new question.
No controllers, weights, experimental results, B0 equations or frozen-snapshot files
were changed, and no simulations/training were run. Next work follows the user's
next request and the confirmed stage-1-to-3 workflow below.

**Previous checkpoint (2026-10-22, cloud Claude session): final Chinese simulation
section and an open question/method decision. Read this block first; it supersedes
the presentation choices in the blocks below.**

What the user wants (stated repeatedly in this session, treat as standing):
- Goal: one top-journal paper (robotics: TMECH/RA-L/T-RO level). Real-robot
  experiments are deferred for now; simulation must already be top-journal level.
- The simulation section must follow how top papers actually write experiments:
  2-4 scenario experiments tied to the research question, each "purpose -> setup ->
  result table/figure -> conclusion", a few representative baselines, plus one
  ablation and computation time. No special additions (no repeated fresh cohorts,
  no seed-pressure tests, no long caveat lists in the main text). When unsure,
  check the papers in `literature/extracted_text/` (Jin E41, Xue, Li E27, Haninger, YOPO).
- Fixed impedance follows field practice (Jin TMECH2023: fixed small/large impedance
  at the method's own bounds). The best-tuned fixed controller B1 must NOT be the
  main baseline (the user rejected it as non-standard). Both candidates C1 and C2
  must appear.
- Report concrete numbers, tables and conclusions directly in the chat, in Chinese;
  the user reads Chinese (via translation). Do not answer with only a file link.
  This overrides the English-response line in AGENTS.md for chat replies; repo
  documents keep their existing language.
- Be careful with quota: say what a large run will cost before dispatching many
  background agents; do not auto-resume large work after rate-limit resets.
- The user has a large local impedance-control literature library on their PC
  (path not yet given; ask for it and request folder access). Use it first for
  the novelty check.

Final section: `research/paper_experiments_final_zh.md` (Chinese). Data: seed-7200
`completed_tasks/` (40 scenes per task, identical scenes for all methods) plus
`completed_tasks/fixed_extremes.json` from `fixed_extremes_eval.py` (compliant fixed:
C1's own rule at k=1.75 kN/m, D=22, KPF=1e-3; stiff fixed: k=25 kN/m, D=313,
KPF=8e-5; M=0.3, K=0, ki=30). Figures: `final_fig_T3.png`, `final_fig_T2.png`
(`render_final_figures.py`, traces from `final_figure_traces.py`).
`research/paper_experiments_merged_zh.md` is an intermediate version; prefer the final one.

Force RMSE (N), 40 scenes: T1 stiffness change / T2 push soft, stiff / T3 combined.
Compliant fixed 7.524 / 0.049, 22.72 / 7.099; stiff fixed 0.170 / 0.775, 1.261 / 0.502;
EtR-A 0.107 / 0.147, 0.551 / 0.230; Haninger 0.122 / 0.178, 0.433 / 0.273;
Jin MS 0.117 / 0.265, 1.268 / 0.285; Xue nominal 0.276 / 0.383, 3.414 / 0.511;
Li learned VIC 0.113 / 0.188, 6.132 / 0.269; C1 0.060 / 0.057, 1.372 / 0.147;
C2 0.083 / 0.143, 0.493 / 0.210. Ablation T3: C1 without surface-motion estimate
1.502, C1 gains only 0.489; C2 gains only 0.217.

Assessment given to the user (honest, agreed by the user as the current state):
- Neither C1 nor C2 yet meets a top-journal method standard. C1 wins T1/T3
  (56-64% of the best existing method's error, significant) but fails sustained
  stiff push/retract (T2 stiff 1.372 vs Haninger 0.433), which contradicts its own
  "separate tissue motion from stiffness" selling point. Its parts (moving-window
  estimation, MPC-VIC, fallback) are mature; joint estimation of environment
  position/motion and stiffness has close relatives (online HC identification,
  Xue, physiological-motion compensation) and no systematic novelty check exists.
  "Near-optimal" is not supported (local reference gaps up to 199%).
- C2 is standard teacher-student distillation and is not significantly better than
  EtR-A on T3 (0.210 vs 0.230).
- Question/method mismatch: the question stresses near-optimal impedance updates,
  while the evidence shows the gain comes from decoupling motion from stiffness plus
  reference compensation (gains-only T3 0.489 vs 0.147).

Status of question and candidates (corrected 2026-10-22 after the user asked why they
had "changed"): nothing has been discarded or changed. The research question in
`research/research_question.md` is unchanged and remains the active question. C1 and C2
both remain candidates; per `long_term_requirements.md` each is improved against its own
specific failures, not dropped. My earlier "demote/drop C2" and "reframe the question"
lines were unconfirmed suggestions made without checking `research/field_state.md`, and
are withdrawn. Note that `field_state.md` already records that unknown/changing
environment position (Liu TCST2021, E45) and joint reference-plus-impedance planning
(Haninger RAS2023) exist, so "decouple motion from stiffness and adjust reference with
impedance" is not by itself a contribution. The ablation finding (gain comes from the
motion/stiffness separation plus reference compensation) is a mechanism result inside
the current question, not grounds to replace it.
Question/method wording gap: the question says "update stiffness and damping", but in both current implementations the admittance stiffness K is fixed at 0 and the actions are D, K_PF and u_ff; most of C1's gain comes from motion estimation plus u_ff. Resolve with the user (reword the question or make the method adapt K).
Known gaps: C1 fails sustained stiff push/retract (T2 1.372 vs Haninger 0.433 N) and its
"near-optimal" claim is unsupported (local reference gaps up to 199%); C2 is not
significantly better than EtR-A on T3 (0.210 vs 0.230 N). The evaluated C1/C2 are reduced
implementations, not the full B0/HJB/spatial-preview/NEOC design.
Next work: build on `research/field_state.md` (stage 1 already exists, not a full
systematic review) and supplement it with the user's local literature library where it
is thin; then, per stages 2-5, decide with the user how each candidate is improved and
whether any question revision is warranted (a substantive revision is the user's call).
Start the new conversation by restating, in Chinese: goal, standing requirements
(especially stages 1-3), current results, this plan, and the first concrete step with
its quota cost; get the user's confirmation before running anything.

**Latest completed checkpoint: focused simulation completion, superseding the
missing-task presentation below.** The user requested finished simulations and
direct Chinese task/purpose/results/conclusions in chat. The
[completed task evidence](../experiments/cobot_scan_benchmark/completed_tasks/README.md)
now contains a pre-simulation protocol, fresh seed-7200 paired scenes (40 per
task), 2,040 primary episodes, 1,040 final-controller inference interventions,
51 exact full-log replay checks, serial raw update timing, and 234 final-C1 local
OCP snapshot comparisons. No frozen controller source/configuration/weights were
changed. The [Chinese experiment section](../research/paper_experiments_draft.md)
was replaced in place with completed evidence; its original-paper appendix stays.

Tasks are T1 with breathing removed, T2 homogeneous soft/stiff push/retract
(balanced material x sign, 5 mm, 0.5 s ramps), and original-distribution T3.
All settings remain frozen: T1 uses old T1 settings, T2/T3 old T3 settings.
Thus T1/T2 are isolation/transfer experiments, not independently retuned rankings;
T3 is the fresh in-distribution main comparison. All 3,080 episodes are finite.
RMSEs C1/C2/best baseline: T1 0.0599/0.0829/0.1070 N; T2
0.7144/0.3177/0.3054 N (Haninger); T3 0.1471/0.2096/0.2298 N.
C1 T2 soft/stiff RMSE is 0.0572/1.3717 N: sustained stiff motion is a clear
transfer deficit. C1/C2 are not uniformly superior. Adjusted paired intervals
support C1 versus all baselines on T1/T3; C2 versus EtR-A/learned VIC on T3
cross zero, despite lower means. Neural seeds are averaged within scenes.

Worst sampled T2 peaks C1/C2/Haninger are 22.228/12.116/8.668 N; T3
7.012/13.067/6.926 N. C2 extrema cover 120 model/scene runs versus 40 classical
runs; finite integration does not imply safety or tracking success.

C1 gains-only T3 RMSE is 0.4888 N versus full 0.1471 N; reference compensation
is important and the full gain cannot all be attributed to impedance updates.
No-surface-state C1 also deteriorates on static T1, so this state may compensate
model mismatch, not only identify physical surface motion. Frozen feature
masking in C2 is an input-distribution intervention, not a retrained architecture.
T3 normal velocity RMS is similar for C1/C2/B1; motion superiority is unproven.
The force/actual-motion figure uses first soft/stiff push scenes, with all data
retained in tables; figure alignment, PDF font size/collision and visual checks
passed. No spectral/jerk/energy/image or new cross-platform sweep was added.

Final-C1 P99 updates: T1/T3 0.1523/0.1379 ms; C2 0.0272/0.0376 ms.
All 101,262 measured software updates are below 8 ms. This excludes physical
communications and low-level plant integration; fixed/Duan samples initialize
constant gains, not periodic decisions. C1 three-start 100-iteration reference
costs are stable from 50 to 100 iterations; mean gaps T1/T3 0.720/2.480%,
maxima 34.946/198.960%. These are local-OCP cost differences, not global
optimality or physical error percentages. Uniform near-optimality and a primary
compute bottleneck are not established; no corresponding C2 optimality claim.
The requested focused simulation evaluation is complete; do not automatically
start retuning, redesign or larger experiments absent the next user request.

**Previous checkpoint: presentation and evidence interpretation.** The user requested original-paper
reading followed by a Chinese, paper-style experimental section with existing
results and explicit missing experiments. The [experiment section](../research/paper_experiments_draft.md)
was updated in place to Chinese at the user's explicit request. It links the
questions, tasks and metrics of 12 directly relevant original papers; this is
focused experimental-section reading, not a fresh full-text reading of all 36.
[Re-extraction code](../experiments/cobot_scan_benchmark/build_paper_experiments.py)
and [tables/source manifest](../experiments/cobot_scan_benchmark/paper_tables/)
keep the comparisons reproducible. No new control simulations were run.

The section includes the repaired-baseline physical results, historical policy
call timing, old C1 numerical-reference gaps, information/frequency ablations,
and newly pulled seed-5000/6000 results. Ordinary one-step GN is already about
0.11-0.13 ms; the old neighbouring update has no demonstrated unique advantage.
C2 fails under fast pushes, and retuned B1 outperforms retuned C1 on PC500.
PC125dead C2 has nonfinite scene outcomes (T3: 1/40, T1: 2/40); the Chinese table
shows NA instead of the finite means obtained by omitting failed integrations.
Cross-period C2 results remain direct, uncalibrated transfers, not fair rankings
of fully adapted algorithms. Raw source results were preserved.

The remaining evidence directly needed for the full claim is final-method normal
motion/vibration and matched end-to-end timing/deadlines, isolated-event recovery
and information/cadence curves, and final-method quality-versus-time against an
explicitly converged reference in a declared action class. Existing C1 reference
gaps use two initializations with 25 GN iterations, not a global-optimum certificate.
T2, the lag2 sweep and an additional two-cycle delay test are not complete; the
newer cohorts also lack all repaired paper baselines. Image quality/actuator energy
are not measured by this kinematic contact model and are only required if those
application claims are made. This presentation task is complete; no automatic
continuation of the old gap-experiment plan is implied.

**Earlier numerical entry, superseding the October 19 status below:** the October 21
[36-paper setup analysis](../research/simulation_setup_reverse_analysis.md) led to
the implemented [cobot benchmark](../experiments/cobot_scan_benchmark/README.md).
The latest user request was to repair Xue/other baseline adaptations, extend weak
tuning, and compare the frozen proposed controllers using physical metrics on
identical scenes. The completed [refined comparison](../experiments/cobot_scan_benchmark/refined_baselines/results.md)
is now the numerical entry. The earlier [frozen replay](../experiments/cobot_scan_benchmark/physical_results/results.md)
remains the before-repair evidence, including its replay differences.

The refined entry restores Xue's variable-stiffness branch and HC projections,
checks stiffness/energy initialization and measured-state/servo alternatives,
fixes Jin OS's action-blind 8 ms prediction endpoint, expands gain/solver tuning,
and restores the learned baseline's reference output with per-task training.
It retains 1,104 coarse configurations, 317 local evaluations, three neural
seeds per task (64,820 added training/selection episodes including calibration),
selected parameter/native-constraint checks and hashes frozen before testing.
The 1,440 repaired episode runs combine with 2,400 unchanged primary records;
all original source/configuration/weight hashes were verified. Neural seeds do
not increase the 40 independent scenes per task/cohort. Tests 3000 and 4000 were
historically inspected and are not called untouched tests. One historical
low-damping integration failure remains visible in the diagnostic rows.

On cohort 4000, C1 v2/C2 mean force RMSE is 0.0673/0.0896 N on T1, versus
0.1161 N for the lowest-RMSE baseline (repaired Jin OS), and 0.1726/0.2308 N on T3,
versus 0.2726 N for EtR-A. Both candidates retain lower mean RMSE than all
competitive baselines in both cohorts; marginal scene-paired intervals are
reported without multiple-comparison correction. This is not dominance on every
physical outcome: T3/4000 Jin MS has a worst sampled peak of 6.499 N versus
C1 v2's 7.164 N; C2 also has a small nonzero zero-force duration. Xue improves
substantially, but its A1 T3 selection still fails the train peak criterion;
Duan remains weak and additional tuning does not consistently improve test
outcomes. Finite search/training does not certify global tuning optimality or
reproduction of the original papers' robot experiments.

The candidates change representation and interface: Candidate 1 is local
stiffness/surface-motion estimation plus direct OCP and feedforward; Candidate 2
uses history/lock-in estimation and a privileged-teacher distilled network.
They do not implement the complete original B0/HJB/spatial-preview design.
Candidate 1 v3 was added in the preceding checkpoint; its main 1 mm configurations
are the same as v2 here. The shared-configuration variant and 4 mm/fast-push cells
were not silently substituted or reevaluated against the repaired baselines.
C1 v1's 46.916 N runaway is retained. No core-question revision was made, and the
baseline-repair request is complete at this checkpoint; further method/robustness
work should follow the next user instruction.

The following sections preserve the earlier design and evidence context. Their
old "latest" experiment references describe the October 19 checkpoint, not the
current numerical entry above. Scientific asset consolidation is recoverable at
`368a29539a723f960f27fb7e28362a65f51428c8`.


## 1. Understand these six points first

1. **Goal:** one paper with a substantive scientific contribution and evidence strong enough to compete at a leading journal. This has not yet been achieved.
2. **Core question:** when the environment changes during continuous robotic scanning of soft tissue, how can new contact information rapidly update near-optimal impedance and improve actual contact-force and motion behavior? This includes environment understanding, impedance decisions and timely execution.
3. **Candidates (original concepts):** candidate 1 builds an optimal-control problem from newly estimated environment information and updates it quickly (originally NNBO critic/HJB plus NEOC); candidate 2 moves the decision burden into offline simulation training in the spirit of YOPO. Both can be improved, combined or retired.
4. **Current implementations (cobot scan benchmark, see `research/optimal_impedance_question_and_methods.md` §2.5):** C1 v2 = moving-window joint estimate of local stiffness and surface motion + direct 128 ms OCP over D, K_PF and u_ff solved by warm-started Gauss-Newton every 8 ms + contact-loss supervisor falling back to B1 gains; critic/HJB and the NEOC neighbouring update were dropped (no benefit; GN already ~0.14 ms). C2 = residual MLP on 0.4 s causal history plus the lock-in stiffness estimate, DAgger-distilled from a truth-stiffness teacher with expectile loss and CMA-ES-tuned; its u_ff head was dropped as unstable. Earlier 2D dual-candidate implementations (local quadratic HJB/Riccati, task-gradient network) are history. None is the complete B0 or a robot validation.
5. **Largest deviations:** excessive defensiveness, retreat/narrowing around strong prior work and misunderstanding candidate 2. An earlier reduced carrier omitted objectives/constraints; the latest implementation must also avoid concealing physical performance deficits with controller-preference costs. Specific failures are recorded. Do not resume mistaken judgments or indefinitely shield the implementation as “only local.”
6. **Current stopping point:** see the 2026-10-22 checkpoint at the top of this file. The final Chinese simulation section is done; the question and both candidates are unchanged; next work supplements `research/field_state.md` with the user's local literature and then improves each candidate against its own failures.

**Latest task and primary objective.** The user requests methods that solve the original question with leading-journal innovation and evidence under the confirmed workflow. Implement/test each candidate separately and improve its own deficits. Experiments follow the question/method, followed by suitable comparisons; none need be forced. These requirements are in the [long-term requirements](long_term_requirements.md). Do not return to fixed control as the main comparison, repeatedly audit the standards or silently change the question.

**Earlier entry (2D, superseded by the cobot benchmark): [dual-candidate task simulations](../experiments/dual_candidate_task_simulation/README.md) and [results](../experiments/dual_candidate_task_simulation/results.md).** Candidate 1 changes the value representation and cannot be called the original B0 root or inherit its certificate. E41 is a source-equation reconstruction adapted to MS scanning, with a different numerical optimizer, not a reproduction of the original Panda experiments. Common physical scoring removes controller-family gain-preference/rate penalties and reports component errors. Use `compiled_timed_*` for independent compiled timing; old `timed_*` files diagnose software overhead only. Failed NNBO adaptations remain isolated. The original question and the candidates' computational division of labor remain.

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

Rows below are the original concepts; current implementations are in point 4 of section 1.

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
| [Current cobot scan benchmark](../experiments/cobot_scan_benchmark/README.md) | Current C1 v2/v3 and C2 implementations, repaired field baselines, seed-7200 completed tasks and the final Chinese section. C1 leads on stiffness change and combined scans but fails sustained stiff push/retract; C2 is not significantly better than EtR-A on the combined scan. |
| [Earlier dual-candidate task simulation (2D, historical)](../experiments/dual_candidate_task_simulation/README.md) | Original question → separate implementations/experiments → Jin MS selection. New 48 training scenarios/3840 contexts, 24 independent scenarios in four categories, three learning seeds, eight distribution-shift and eight independently timed cases. Hertz/objective relationships restored; synthetic 2DOF model, CV/fixed-bias MHE instance and local value feedback rather than B0 root. Comparator force/tangential behavior is better; the candidates' normal-position errors are slightly smaller. No completion claim. |

**Overall: the original question remains; a specific important advance and complete method advantage are not yet established.** Platform/control interface, material/sensor calibration and references remain to be verified. These are candidate concerns, not selected new topics or fixed gates.

## 6. What to do next, and to what extent

The table records the agreed continuation order. **A new user task takes precedence; standards calibration is not an approval gate to repeat every turn.**

| Priority | Work/deliverable | Completion/boundary |
|---|---|---|
| **P0: requirements implementation complete** | Confirmed method-formation understanding and paper-analysis lessons are in the requirements; source evidence remains in the analysis. | Continue at the appropriate stage without repeating understanding audits. Actual design/evidence must demonstrate implementation; a completion record is not method success. |
| **P1: improve each candidate's actual task effect** | On the cobot benchmark: fix C1's sustained stiff push/retract failure and fallback bounce without hurting stiffness-change/combined results, and support or drop the "near-optimal" claim; improve C2 beyond EtR-A on the combined scan. A low gain-preference cost is not a physical win. | Comparisons follow the question/method; do not impose a fixed MPC-based pass/fail gate or force an unsuitable NNBO adaptation. Deepen Jin/E27/NNBO where useful and reproduce/adapt fairly only if suitable. Explicitly state representation changes and evidence limits; do not retreat from the question because strong prior work exists. |
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
