# Related-work scope and current coverage

Assessment date: 2026-10-30. This is a targeted related-work audit, not a systematic review, a novelty verdict, or a revision of the approved research question. Local source records are in [coverage sources](../literature/related_work_coverage_sources.json); the working Chinese related-work text remains [field state](field_state.md).

## What the paper needs to explain

The task is continuous force-regulated ultrasound scanning while tissue moves and local contact mechanics vary. The related work must explain existing ways of maintaining useful contact under those changes, how observations affect feedback/reference decisions, and the conditions and costs of doing so. It must not define relevance as “uses the same stiffness/drift decomposition as C1.” Direct force adaptation, lumped-disturbance compensation, external motion sensing, and learned/reference policies are legitimate alternatives.

The parent research question still concerns using online environmental information to promptly improve interaction control. The current C1 is one implementation. Neither its chosen inputs nor its polynomial model defines the entire field, and the present audit makes no new method selection.

## How original-paper introductions calibrate the scope

| Inspected source | What its introduction actually organizes | Consequence for this project |
|---|---|---|
| Jin, TMECH2023 E41, PDF1–2, Introduction | Precision/compliance tradeoffs; demonstration/RL alternatives; predictive optimization; objective and constraint limitations | Organize work around control capability and conditions, rather than listing every impedance or MPC paper. General optimization references explain tools only when they support a concrete method statement. |
| Li, T-RO2025 E27, PDF1–2, Introduction | Unknown stiffness/position in force tracking; reference and impedance adjustments; adaptive and learned approaches | Unknown environment position and reference correction are already within the relevant field. An estimator is not the only valid way to address force tracking. |
| Xue, RAS2025, PDF1–2, Introduction | Unknown soft interaction, safety/performance objectives, online contact estimation, predictive constraints and execution | RAS is retained because of direct scientific proximity, not excluded for being outside the project's target venue list. Its claims about other methods are not automatically adopted as facts. |
| Independent reverse analysis, Sections 2–5 | R12 connects desired interaction behavior to mechanism roles; R11 separates feedforward and impedance roles; R01/R19 show that concise changes can support capabilities | Use these papers to calibrate argument and evidence design. They are not automatically direct ultrasound comparators or a universal publication checklist. |

The inspected original introductions and the reverse analysis support a compact, capability-led related-work section. They do not prescribe a citation quota, require every source to solve every condition, or prove that C1 has a new contribution.

### Depth required for this project

The user's research-experience collection (`D:/桌面/科研经验总结`, README, workflow/01_literature, references/stages/01_literature and its Li, Xue, TUM-ultrasound and Jin comparison examples, inspected 2026-10-29) supplies the reasoning standard, rather than a new permanent project rule. Li's T-RO2025 introduction organizes uncertain-contact force tracking around reference/gain adaptation and learning; Jin's TMECH2023 introduction connects the objective, admissible control variables and constraints. The TUM sequence illustrates a different, task-led question: which part of ultrasound acquisition can actually be completed under motion. These are complementary calibrations; an imaging-system result is not automatically a variable-impedance result.

For the present **contact-control method study**, the account must explain the strongest relevant alternatives through a common relation: **available online information → model, reference, gain or force-target update → actual force/motion benefit → preparation, computation and execution conditions**. This connects several routes without pretending that they share one universal failure cause. Closest control neighbors and closest ultrasound systems need separate, explicit comparisons. Read their methods and experiments, acknowledge supported strengths first, and locate a remaining boundary in an assumption or mechanism. Lack of a reported test is an evidence limit, not demonstrated failure.

The deliverable is therefore a compact route account with a capability/condition/cost comparison and a precise statement of this project's current position. It must distinguish environment stiffness from controller stiffness, nominal trajectory from physical surface motion, gain learning from target-force learning, and quasi-static geometry variation from independent tissue motion. It must also explain why parameter estimation, lumped disturbance compensation, predictive decisions and learned policies coexist. A stiffness/drift decomposition is a design choice whose additional control value is still to be tested; it cannot define the inclusion boundary or manufacture the literature gap.

Completion means that original evidence supports these comparisons and a targeted refinement plus backward/forward citation check reveals no new direct neighbor that materially changes them. It does not require every downloaded paper in the prose, every adjacent paper reread, or a proof that all existing methods fail. Remaining non-decisive identity/version limits can be disclosed. A missing closest-neighbor method or unresolved decisive timing claim must be resolved before claiming the current scope is complete. This literature task does not certify C1's novelty, superiority or publication readiness, and does not authorize method changes or experiments.

## Inclusion by scientific role

| Role | What belongs | What must remain explicit |
|---|---|---|
| Direct task prior | Ultrasound/medical scanning contact control with relevant force, compliance, motion, or contact-retention results | Geometry, vision, US-image input, pre-calibration, and motion conditions may differ from C1. These differences are comparison dimensions, not automatic exclusion grounds. |
| Direct control mechanism | Unknown-contact force tracking, online environment identification, predictive/adaptive impedance/admittance, joint reference/gain decisions, motion or disturbance compensation | Manufacturing and human-interaction tasks qualify when the control relationship is directly relevant. Do not equate robot-model identification with tissue identification. |
| Supporting foundation | Identifiability/excitation and joint geometry/compliance estimation; stable variable-gain execution; contact transitions and model uncertainty | A condition for an ideal impedance model is not a theorem for C1's delayed position-servo loop. A static estimation result is not dynamic scanning validation. |
| Method transfer/background | General optimization differentiation, RTI/NEOC, offline planning/learning, generic uncertainty methods, and task context | Keep separate from already demonstrated target-field capabilities. Their presence in the library is not a reason to cite them in the main related-work argument. |
| Outside the core | Pure diagnosis, reconstruction, segmentation, path planning, unrelated locomotion/SLAM, or intent inference without a material connection to the declared control question | A paper may be valuable elsewhere. Exclusion here is about scope, not scientific quality. |

The practical groups in the Chinese field assessment remain: adaptive/optimal force interaction; predictive contact-model updates; learned/reference-and-impedance control; physiological motion and disturbance compensation; joint contact estimation and stable execution. Adjacent image-reconstruction work is retained only where needed to distinguish the user's named sources, not as evidence of contact-control capability.

## Local findings that change the assessment

| Source added from the supplied library | Inspected source locations | Why it belongs and what it changes |
|---|---|---|
| H01: Fu et al., TIM2024, DOI 10.1109/TIM.2024.3372209 | PDF1–4, Sections II–III; selected PDF5–7 experiment passages | Direct medical-contact/scanning prior: force-error feedback, online QP robot-stiffness selection, N=1 implementation, and energy-tank execution. Tests include soft/rigid material transitions and uneven/sloped surfaces. “Online optimized impedance for medical scanning” is therefore not a sufficient standalone C1 contribution. The inspected formulation does not establish C1's particular local stiffness/equivalent-drift predictor. |
| H02: Dyck et al., RA-L2022, DOI 10.1109/LRA.2022.3184800 | PDF1–3, II-A; PDF7–8, III-B/Fig.8 and IV | Direct ultrasound contact prior: surface coordinates and negative-distance impedance are compared with feedforward force strategies, including observed contact loss and recontact impact. Contact retention is not an unprecedented outcome. The given mesh/silicone experiments and the future-work paragraph on stiffness/mesh updates must not be generalized into arbitrary moving-tissue capability. |
| H03: Duan et al., arXiv:2203.02126v1 | PDF1, PDF3–6, II-B/C, III and selected IV | Task-specific optimization and variable-impedance precedent. The QP coordinates robot/contact tasks; gains come from expert demonstration covariance and US images guide scanning. This differs from online contact-parameter inference. The retained file is an author manuscript; its final venue/version has not been independently verified in this audit. |
| H04: Kronander and Billard, T-RO2016, DOI 10.1109/TRO.2016.2593492 | PDF1–3, II–III/Theorem1; PDF6, V | Supporting stability foundation: sufficient conditions link varying stiffness/damping under constant target inertia and associated assumptions. It prevents presenting generic variable-impedance stability as new; it does not certify C1 with K=0, integral force control, changing reference, delay, and an identified position servo. |

No source-paper RMSE is ranked against the project's synthetic benchmark: platforms, loads, materials, sensing, references, and metrics differ. No new controller or simulation was run in this audit.

## What was not imported into the core

- Caldarelli et al., *Perturbation-Based Stiffness Inference in Variable Impedance Control*, DOI 10.1109/LRA.2022.3187866, infers desired **robot-control stiffness** through demonstrations and human intervention. It is not evidence of online **environmental stiffness** identification in moving tissue. Its older context-specific screening decision is left intact; this audit records the present role separately.
- Riccati-ZORO, DDP/implicit differentiation, decision-focused learning and generic safety frameworks remain possible method resources, not direct ultrasound prior work merely because C1 uses optimization.
- Robot-dynamics/uncertainty compensation, EMG intent inference and grasp/peg-in-hole policies are not relabelled as physiological motion/stiffness separation. Representative learned contact-control precedents already exist in the project; additional similar sources are not mandatory main-text additions.
- Alternative or annotated versions of Anand/E40, Pollayil, Wijayarathne and R01 are not counted as additional works. The alternate local H02 manuscript and journal PDF share one work identity.
- The scanned *Direct adaptive impedance control* file exposes only licence text through ordinary extraction; its content is unresolved here. A failed text extraction is not a scientific exclusion verdict.

## Current coverage and remaining work

The account now includes the three supplied originals, retaining their actual control roles:

| Source and inspected passages | What it establishes | Boundary affecting this project |
|---|---|---|
| Piedra IJRR2026 C10, PDF6–8,15–18,24 | Torque-executed explicit force/moment control, including independent moving-phantom comparisons with retuned IAAC | Position-interface admittance cannot inherit the execution mechanism unchanged; motion handling is already demonstrated, and coupling confidence is not diagnostic utility |
| Zhang RA-L2025 C11, PDF3–7 | Pre-calibrated force–strain and geometry predict desired force; image confidence updates the force target; fixed feedback mapping executes it | Target-force adjustment is not gain selection. The published Eq11–12 direction ambiguity is preserved rather than silently corrected or used to dismiss the paper |
| Li TAFAC CEP2026 S04, PDF1,3–8 | Known reference velocity/acceleration and control-cycle compensation with fixed admittance | Online2025-12-06, DOI10.1016/j.conengprac.2025.106674. Unknown patient motion differs from commanded motion; the delay proof uses a first-order approximation; AVAC has slightly smaller RMSE in one reported condition |

The archive has106 unique PDFs and109 usage records, with27 accepted selected-passage cards.32 query attempts include27 completed events and five retained historical failures;221 candidate/local/alias records have118 OUT_OF_SCOPE,96 IN_SCOPE,3 DUPLICATE and4 HOLD_IDENTITY decisions. These counts do not establish saturation. Original DOI/version/pages/hashes are in the materials manifest. Preexisting original and text bytes remain unchanged.

The fresh [independent review](../literature/aris_coverage_review.json) returned **CONTINUE**, with one specific material doubt: Santos and Cortesão TASE2018, *Computed-Torque Control for Robotic-Assisted Tele-Echography Based on Perceived Stiffness Estimation*, DOI10.1109/TASE.2018.2790900. Its actual abstract links pre-contact anticipation, in-contact perceived stiffness, effective robot mass and force adaptation, with task validation. The prior broad-family retention reason does not resolve whether this interface changes the capability or bottleneck account. No local original is retained. Request it from the user rather than inferring its detailed mechanism from the abstract.

The review supports the four bounded development traces and the coherence of the inspected branches. It explicitly does not clear all221 candidate dispositions. The candidate canonical map had declared SUFFICIENT for review; that declaration was not accepted as a completed landscape. Controller state is QUERY_PLANNING with the review gap. The mirrored map remains the exact reviewed candidate, not an edited acceptance record. A minor problem–method crosswalk omission for the three new cards is recorded for the next lawful map revision; the family/assumption rows already contain them.

See [manual requests](../literature/manual_fulltext_requests.md) for the single current original request. The previous C10/C11/S04 download issues are closed. Missing Jiang2024 identity and other non-decisive version limits remain disclosed and carry no fine-grained claims. No blanket expansion or exhaustive-coverage assertion follows.

## Connection to the present scientific assessment

The user's current task additionally authorizes assessing and refining the question and proposed method. This does not turn the ARIS literature run into a new idea pipeline. The section-by-section standards, reference-paper calibration, question wording proposal and C1 mechanism derivation are in [methods §7](optimal_impedance_question_and_methods.md#7-现状问题方法的阶段评估与完善2026-10-30). Current experiments are unchanged. Resolve the material original, revise the same account/map, and reassess; do not advance automatically to experiments.
