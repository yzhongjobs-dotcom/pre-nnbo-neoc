# Related-work scope and current coverage

Assessment date: 2026-10-29. This is a targeted related-work audit, not a systematic review, a novelty verdict, or a revision of the approved research question. Local source records are in [coverage sources](../literature/related_work_coverage_sources.json); the working Chinese related-work text remains [field state](field_state.md).

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

## Coverage judgment and next search

The earlier set was broad but not complete even within the supplied library: H01 and H02 supply directly relevant capabilities that were missing from the working account; H03 supplies a task-specific gain/optimization precedent and H04 a foundational execution reference. Adding them repairs identifiable omissions, not the entire coverage problem.

At audit start the library contained 87 unique PDF hashes and 90 usage records. The field assessment linked 37 distinct local PDFs plus online sources; these numbers measure different things. The four additions yield 91 unique PDFs and 94 usage records. None is a count of “all relevant studies.” The local source folder has 145 PDF-named files, 141 with a PDF signature and 138 unique valid-file hashes. At file level, 71 valid files already matched project hashes and 70 did not; title/version checks are needed before calling any of those 70 an omitted work.

Remaining coverage gaps are scientific, not a target paper count:

1. Obtain decision-grade method passages for DOSMAC and Zhang's MHE-VIC; their current abstract/excerpt access cannot settle all distinctions from C1.
2. Complete the direct joint geometry/compliance and dynamic force-regulation citation chains, distinguishing palpation/static registration from moving-tissue control.
3. Search recent observer/predictive-admittance ultrasound work, including the existing *Autonomous robotic ultrasound scanning on complex surfaces* lead, and verify identity and publication timing.
4. Use backward/forward citations and mechanism-specific queries to look for counterexamples to a claimed gap, rather than searching only C1 terminology.

## ARIS supplementation delivered for acceptance: 2026-10-29

This run supplements an existing project literature account. Its source policy remains bound to the approved 2026-10-26 publication cutoff; retrieval and user-file intake occurred on October28–29. It does not restart ideas, select a method, change C1/C2, or run experiments.

### Available-material result and self-corrections

Of the 11 user-downloaded PDF-named files, ten have valid PDF signatures and matching original identities. They were retained byte-for-byte with extracted texts. S03 is an Elsevier security-verification HTML page, not a paper; it remains in the ignored inbox and contributes no original-method evidence. The library is now 101 unique PDFs /104 usage records. Existing original and extracted-text hashes and the frozen B0 snapshot were preserved.

The Chinese field account links51 different local original PDFs, plus explicitly graded online leads; this is separate from the101-PDF library count. It now integrates five additional task/control or joint-estimation sources (C02,C03,C04,C06,C07), deepens two previously abstract-level closest priors (DOSMAC S01 and MHE-VIC S02), and uses C05 briefly for execution foundations. C01 provides review calibration. C08 is retained as a cutting/CBF execution resource with unresolved earliest publication date, outside verified pre-cutoff core claims. Library retention does not determine main-text inclusion.

| Source | Original capability now supported | Condition that matters for this project |
|---|---|---|
| C02 Roveda/Piga2021 | Dual EKF force/stiffness estimates inform Cartesian coupling compensation and impedance rules | Sensorless/torque-model assumptions and estimator/controller bandwidth separation |
| C03 Wu T-RO2026 | Passive constant-torque buffering plus active deformation control and modular force-range changes | Additional hardware; axial/normal force distinction under inclination; original p1 verifies publication on2026-09-15 |
| C04 Huang TIE2025 | Structured adaptive feedforward plus robust feedback commands motion for force tracking | Model K has no special physical meaning; bounded parameters and uncertainties, not physical tissue decomposition |
| C06 Duan TASE2025 | Ultrasound-image confidence drives learned target-force increments | Not impedance-gain learning; expected cumulative constraints are not pointwise hard safety |
| C07 Chalasani ICRA2016 | Continuous varying-depth/force palpation estimates geometry and stiffness through GP models | Quasi-static physical data; adaptive-path simulation; dense-data reference is not independent material truth |
| S01 DOSMAC CEP2026 | Lumped observer/sliding compensation in zero-stiffness admittance improves its robot-scan force errors | Known direction/low-frequency conditional model; fixed curved surfaces in hardware and true nominal surface position in simulation |
| S02 MHE-VIC TIE2024 | MHE estimates environment stiffness/damping and robot disturbance, feeding reference/model compensation and gain rules | Surface measured beforehand, static normal position, local slow-change and gain-rate conditions; torque execution |
| C05 Zeng TIE2026 | Smooth force/admittance activation plus tank injection/freezing in assembly | Supporting execution mechanism; assembly is not ultrasound evidence |

Three already-local originals were additionally admitted and read through ARIS to repair a mismatch between the increment map and the project account: E27 supplies actual learned gain/reference increments, prior-trial regression and pretraining conditions; G23 supplies the joint position/gain action interface and simulated constant-force manipulation; E41 supplies direct predictive gain optimization with distinct OS-QP and MS nonlinear/soft-state-constraint behavior. No new PDF was added for these reused identities.

The self-check corrected four issues: (1) the selected download set was not treated as a mandatory main-text bibliography; (2) missing formal-registry entries were repaired by reusing existing representatives instead of restarting the project search; (3) source mechanisms and evidence replaced repetitive defensive language; (4) simulated/physical motion, force-policy/gain-policy actions, and model parameters/physical stiffness were kept distinct. Experimental run history was removed from field-state prose; it remains in its original experiment reports and Git history. No new permanent requirement was invented.

### Actual search and evidence boundary

The [record-level data](../literature/aris_literature_search.json) preserves 25 query attempts: 20 completed, five historical zero-hit adapter failures. The latest 158 candidate/local/alias records comprise 103 OUT_OF_SCOPE, two DUPLICATE, two HOLD_IDENTITY and 51 IN_SCOPE; these are not unique-work or core-reference counts. Twenty accepted cards explicitly cover original selected passages, not full-paper reading or independent theorem verification.

The last two queries were a no-lower-era mechanism/backbone refinement and an ACT-US respiratory-policy identity refinement. The former returned two actual abstracts about location/stiffness-informed reference control; their citation counts do not meet the approved high-citation thresholds and they are not relabelled high-citation Evidence. The latter returned off-scope health/legal/index noise, which was screened out. Repeated ACT-US metadata verification still returned no actual abstract. Empty or noisy responses are not proof of absence.

An independent fresh-context review verified all six artifact hashes and returned **CONTINUE**; its [exact verdict](../literature/aris_coverage_review.json) was accepted by Controller using the actual native-child transcript, with no Hook provenance claim. It identified learning-representative, historical-backbone and decisive-access/metadata gaps. The E27/G23/E41 registry omission was corrected afterward. Historical high-citation coverage and inaccessible direct-task evidence remain limited; no second sufficient-coverage verdict is claimed. The archived [formal map](../literature/aris_active_field_map.md) is the last 17-card PARTIAL version; the three further cards await the next allowed map revision after the pending reading batch. The project field account already incorporates their findings.

### Delivery stop and remaining originals

The supplied eight-paper batch was completed. After the bounded correction, Controller now requests one original: ACT-US, DOI10.1109/IJCNN69854.2026.11728409. The same [manual checklist](../literature/manual_fulltext_requests.md) also identifies the invalid S03 DOB-MPC file needing replacement. Neither source is described as having or lacking a method capability from its title/snippet. Their earliest public dates matter because issue/DOI year alone does not establish eligibility under the approved cutoff.

Current stage: **HUMAN_SEARCH_REQUIRED / human_fulltext_batch**. Deliver the updated field account for user acceptance and stop, as requested. Do not automatically resume search, create an idea, change methods or run simulation. Remaining gaps are disclosed rather than converted into a false completion or novelty verdict.
