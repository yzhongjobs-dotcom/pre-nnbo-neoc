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

## ARIS evidence and coverage checkpoint: 2026-10-29

The approved task is to supplement the existing research status. The policy and its 2026-10-26 publication cutoff remain unchanged; actual retrieval occurred on October 28–29. No idea selection, method revision or simulation is started by this checkpoint.

The Controller completed the initial review/primary reading and one map-guided expansion. Across the two plans and two fresh retries there are 20 query attempts: 15 completed and five preserved historical zero-result adapter failures. The fixed adapter treats a successful empty Scholar response as zero hits; two gap queries were retried with new budget and completed empty. This does not retroactively erase failures or establish saturation. The follow-up included a small Lin forward-citation response; actual citation arguments, queried dates/pages and reasons are retained in the [search data](../literature/aris_literature_search.json). Planned venue lists are not claimed as executed venue filters or an exhaustive venue sweep.

The corpus now contains 130 records: 124 gateway candidate identities plus six local registrations, including explicit aliases. Latest decisions are 83 OUT_OF_SCOPE (excluded or resource-only), two DUPLICATE, two HOLD_IDENTITY and 43 IN_SCOPE. Screening bases are 39 TITLE_ABSTRACT, 64 TITLE_ONLY, 18 TITLE_ONLY_ABSTRACT_UNAVAILABLE, seven FULL_TEXT and two IDENTITY_UNRESOLVED. These are record-level administrative counts, not 130 unique works or 43 established capabilities. The two identity holds remain unread and cannot support inclusion or irrelevance. Two Scholar2026/metadata2027 records retain their date conflict.

### What the accepted originals change

Seven accepted partial-text cards preserve exact read-event/source-hash bindings. The reviews calibrate concepts; the originals establish specific capabilities, rather than using titles or publication prestige as proof.

| Source and inspected packet | Source-grounded finding | Placement and limit |
|---|---|---|
| E32, G22 reviews | Robot desired impedance, environment parameters, reference adaptation and learned impedance are distinct objects; stability and learning conditions remain conditional. | Terminology and family cross-check, not a substitute for primary capability evidence. [E32](../literature/aris_evidence/E32.json), [G22](../literature/aris_evidence/G22.json). |
| Lin TMECH2021 E57, selected PDF1–8 | RLSE estimates stiffness/damping; Remark4 obtains equilibrium position at first contact. A reference-force/motion QP is separate from experiment-specific preset target-matrix scheduling; the model assumes stable single contact and decoupled Kelvin–Voigt interaction. | Direct control mechanism in field-state Section2. Neither continuously recovered tissue motion nor online task-optimal gain selection is established by these passages. [Card](../literature/aris_evidence/E57.json). |
| Fu TIM2024 H01, selected PDF3–5 | Medical force tracking already uses robot-stiffness QP, energy-tank execution, N=1, soft/rigid transitions and uneven-surface scanning. | Predictive/optimized contact-control route, with declared low-speed approximation. The new reader packet verifies mechanism/setup, not a new performance-table audit; earlier reading scope remains intact. [Card](../literature/aris_evidence/H01.json). |
| Wang TIE2026 G02, selected PDF5–6,8–9 | Historical force-prediction QP updates both impedance damping Bd and environment stiffness Ke; a separate MPC chooses future velocity corrections. The force-direction equilibrium velocity is set to zero in the inspected prediction approximation. | Correct the earlier compressed description: damping is optimized, while the two QPs have different direct objectives. Maximum control frequency is not a worst-case solver-time certificate. [Card](../literature/aris_evidence/G02.json). |
| Cao TIM2025 G13, selected PDF4–5,7–9 | External visual sinusoidal prediction and triggered refitting supply the admittance reference; the one-cycle position prediction is not evidence of finite-horizon gain optimization. The 5N/5mm/0.2Hz phantom case reports tracking-error mean±SD of .29±.21N against 1.28±.83N and 2.10±1.00N for its two comparators. | Physiological-motion compensation route, with fitting transients and measured signal latency. These are paper-specific errors, not C1 RMSE comparisons or clinical certification. Do not conflate nominal signal-update time and induced-delay thresholds. [Card](../literature/aris_evidence/G13.json). |
| CARE RSS2014 G21, selected PDF2–4 | Static IEKF jointly estimates registration and relative stiffness with prior geometry, force and reliable tip position. Its rank-nine example uses three known, differently oriented planes, each palpated at two force levels. | Joint-estimation foundation in field-state Section5. No general dynamic-tissue identifiability or gain-update guarantee follows from this example. [Card](../literature/aris_evidence/G21.json). |

These papers are placed by their control/inference role, not collected in an unrelated new section. No source-paper metric is ranked against the project's simulations. The library remains 91 unique PDFs / 94 usage records; no PDF was added during this ARIS reading checkpoint.

### Discovery leads that still need originals

The new TASE2025 *Safe Learning by Constraint-Aware Policy Optimization for Robotic Ultrasound Imaging* is a directly relevant task lead: its actual abstract describes force-constrained policy learning and spinal-phantom evidence. It must remain outside new source-grounded capability claims until its actions, constraints and experiments are read. Do not assume that policy optimization means impedance-gain optimization.

Other selected leads cover sensorless environment estimation, adaptive robust force control, energy-aware force admittance, continuous adaptive palpation and soft-tissue passive/CBF execution. The lightweight T-RO system and latest TMECH identity also require publication-timing checks. None is included as demonstrated dynamic-tissue capability merely from a title or abstract. Pure reconstruction, planning, needle-position tracking, soft-robot actuator deformation and generic optimization remain outside core contact-control evidence unless an explicit control relationship is verified.

### Current stop and exact continuation

`finish-reading` accepted all seven available-original cards and entered **HUMAN_SEARCH_REQUIRED / human_fulltext_batch**, with [one eight-paper request](../literature/aris_human_fulltext_request.json). The [manual-acquisition checklist](../literature/manual_fulltext_requests.md) adds the three pre-existing DOSMAC, MHE-VIC and DOB-MPC ultrasound gaps in the same user-facing list, clearly separate from the formal eight IDs. The DOB-MPC lead has an existing January2027 issue record and an unverified earliest online date; no pre-cutoff capability is attributed from the DOI year.

The two unresolved identities no longer block reading other verified sources: they are explicitly held rather than falsely excluded. Cwd and Hooks are not prerequisites. The canonical Active Field Map is still PARTIAL and currently synthesizes the first E32/G22/E57 evidence; the four new cards await the next allowed map revision after the human batch resumes reading. No independent final coverage acceptance has occurred.

After the user supplies originals, verify and register the batch, read methods/experiments, revise the same map and research status, and continue discriminating citation-chain/gap checks. A small forward-citation return and failed exact-title searches do not establish that no further direct neighbor exists. Completion and novelty remain unresolved until the decisive originals and coverage gaps are addressed.
