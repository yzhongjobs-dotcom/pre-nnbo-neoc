# Literature index

This directory retains original papers, extracted texts, sources and actual reading scopes for subsequent research. Use field-state sources, method-transfer references and independent top-venue samples separately. Downloading/extraction is neither full reading nor theoretical reproduction.

## Usage and sources

- Scientific entries: [field state](../research/field_state.md), [candidate methods](../research/optimal_impedance_question_and_methods.md), [transferable mechanisms](../research/transferable_mechanisms.md).
- Local files/SHA256: [materials manifest](local_materials_manifest.json). Each record identifies PDF/text paths, purpose, source URL, known version and reading scope. Deduplicated files are shared rather than copied.
- Requested repository 1: [academic-writing-skill literature](https://github.com/yzhongjobs-dotcom/academic-writing-skill/tree/c88d48eb08a876eb9b81a30d37b62d70fa028581/文献资料). Screened 22 identities/abstracts and examined key passages in ten papers; those ten are local. Metadata for twelve others remains in [sample sources](top_venue_sample_sources_and_reading_scope.json), without implying full-text reading.
- Requested repository 2: [e2e/User Batch V2](https://github.com/yzhongjobs-dotcom/e2e/tree/8fff5f3cbb1d47405f54888c0e92c61d08e31ced/source-materials/user-batch-v2) and [idea-stage](https://github.com/yzhongjobs-dotcom/e2e/tree/8fff5f3cbb1d47405f54888c0e92c61d08e31ced/idea-stage). Of the original 61 screened papers, local material includes 25 directly relevant, six partly relevant and five background papers; 25 currently irrelevant papers retain exclusion reasons only. Ten originals added from existing leads are also local.
- Complete screening, original pages and neighbors without obtained full text appear in [sources/screening](impedance_literature_sources_and_screening.json). Do not copy all 2653 discovery records, 105 old evidence cards or 3822 search events, or treat them as confirmed conclusions.

There are **65** deduplicated local PDFs with corresponding texts. Verify selected IDs with `python -X utf8 literature/verify_materials.py --ids E27 F08 YOPO`; verify all with `python -X utf8 literature/verify_materials.py`. `--download` restores missing selected files and checks hashes; it does not search/screen again. Extraction uses system `pdftotext`. Record version-dependent text differences explicitly; do not overwrite original identity hashes.

## Original methods and direct evidence

| ID | PDF and text |
|---|---|
| beber | [beber](papers/beber.pdf) · [text](extracted_text/beber.txt) |
| pollayil | [pollayil](papers/pollayil.pdf) · [text](extracted_text/pollayil.txt) |
| wijayarathne | [wijayarathne](papers/wijayarathne.pdf) · [text](extracted_text/wijayarathne.txt) |
| haninger | [haninger](papers/haninger.pdf) · [text](extracted_text/haninger.txt) |
| xue | [xue](papers/xue.pdf) · [text](extracted_text/xue.txt) |
| risk_sensitive | [risk sensitive](papers/risk_sensitive.pdf) · [text](extracted_text/risk_sensitive.txt) |
| fu | [fu](papers/fu.pdf) · [text](extracted_text/fu.txt) |
| piedra_volumetric | [piedra volumetric](papers/piedra_volumetric.pdf) · [text](extracted_text/piedra_volumetric.txt) |
| NEOC | [NEOC](papers/NEOC.pdf) · [text](extracted_text/neoc.txt) |
| implicit_differentiation | [implicit differentiation](papers/implicit_differentiation.pdf) · [text](extracted_text/idoc.txt) |
| NNBO | [NNBO](papers/NNBO.pdf) · [text](extracted_text/nnbo.txt) |

## Field-related and background material

Inclusion grounds and reading pages are governed by source JSON. Direct relevance does not imply completed full-text reading. E30 reuses the existing risk_sensitive version; external E30 hashes and the local version are recorded separately.

| ID | Role | PDF and text |
|---|---|---|
| E01 | Background review | [A Tutorial Survey and Comparison of Impedance Control on Robotic Manipulation](papers/field_papers/E01.pdf) · [text](extracted_text/field_papers/E01.txt) |
| E04 | Directly relevant | [A Novel Robust Adaptive Impedance Control Scheme for Contact Force Tracking With Unknown Environment](papers/field_papers/E04.pdf) · [text](extracted_text/field_papers/E04.txt) |
| E06 | Directly relevant | [Adaptive Energy-Bounding Approach for Robustly Stable Interaction Control of Impedance-Controlled Industrial Robot With Uncertain Environments](papers/field_papers/E06.pdf) · [text](extracted_text/field_papers/E06.txt) |
| E09 | Directly relevant | [Adaptive Impedance Control of HumanRobot Cooperation Using Reinforcement Learning](papers/field_papers/E09.pdf) · [text](extracted_text/field_papers/E09.txt) |
| E11 | Partly relevant | [An Optimal Variable Impedance Control With Consideration of the Stability](papers/field_papers/E11.pdf) · [text](extracted_text/field_papers/E11.txt) |
| E12 | Background review | [Applied Bionics and Biomechanics - 2018 - Al-Shuka - Active Impedance Control of Bioinspired Motion Robotic Manipulators ](papers/field_papers/E12.pdf) · [text](extracted_text/field_papers/E12.txt) |
| E14 | Directly relevant | [Constrained Variable Impedance Control using Quadratic Programming](papers/field_papers/E14.pdf) · [text](extracted_text/field_papers/E14.txt) |
| E15 | Partly relevant | [Deflection-Domain Passivity Control of Variable Stiffnesses Based on Potential Energy Reference](papers/field_papers/E15.pdf) · [text](extracted_text/field_papers/E15.txt) |
| E16 | Directly relevant | [Development of a Safety- and Energy-Aware Impedance Controller for Collaborative Robots](papers/field_papers/E16.pdf) · [text](extracted_text/field_papers/E16.txt) |
| E17 | Partly relevant | [Diffusion-Based Impedance Learning for Contact-Rich Manipulation Tasks](papers/field_papers/E17.pdf) · [text](extracted_text/field_papers/E17.txt) |
| E18 | Directly relevant | [Dynamic Adaptive Hybrid Impedance Control for Dynamic Contact Force Tracking in Uncertain Environments](papers/field_papers/E18.pdf) · [text](extracted_text/field_papers/E18.txt) |
| E23 | Directly relevant | [Force tracking impedance control of robot manipulators under unknown environment](papers/field_papers/E23.pdf) · [text](extracted_text/field_papers/E23.txt) |
| E24 | Directly relevant | [Human-Inspired Adaptive Optimal Control Framework for Robot-Environment Interaction](papers/field_papers/E24.pdf) · [text](extracted_text/field_papers/E24.txt) |
| E25 | Background review | [Impedance and Interaction Control](papers/field_papers/E25.pdf) · [text](extracted_text/field_papers/E25.txt) |
| E27 | Directly relevant | [Impedance Learning-Based Adaptive Force Tracking for Robot on Unknown Terrains](papers/field_papers/E27.pdf) · [text](extracted_text/field_papers/E27.txt) |
| E28 | Directly relevant | [Impedance Learning for Human-Guided Robots in Contact With Unknown Environments](papers/field_papers/E28.pdf) · [text](extracted_text/field_papers/E28.txt) |
| E29 | Directly relevant | [Impedance Learning for Robotic Contact Tasks Using Natural Actor-Critic Algorithm](papers/field_papers/E29.pdf) · [text](extracted_text/field_papers/E29.txt) |
| E30 | Directly relevant | [Impedance Optimization for Uncertain Contact Interactions Through Risk Sensitive Optimal Control](papers/risk_sensitive.pdf) · [text](extracted_text/risk_sensitive.txt) |
| E32 | Background review | [Impedance Variation and Learning Strategies in HumanRobot Interaction](papers/field_papers/E32.pdf) · [text](extracted_text/field_papers/E32.txt) |
| E34 | Directly relevant | [Learning-Based Predictive Impedance Control Towards Safe Predefined-Time Physical Robotic Interaction](papers/field_papers/E34.pdf) · [text](extracted_text/field_papers/E34.txt) |
| E37 | Directly relevant | [Model-Based Reinforcement Learning Variable Impedance Control for Human-Robot Collaboration](papers/field_papers/E37.pdf) · [text](extracted_text/field_papers/E37.txt) |
| E38 | Directly relevant | [Model-Based ActorCritic Learning of Robotic Impedance Control in Complex Interactive Environment](papers/field_papers/E38.pdf) · [text](extracted_text/field_papers/E38.txt) |
| E39 | Partly relevant | [Model-Based Predictive Impedance Variation for Obstacle Avoidance in Safe HumanRobot Collaboration](papers/field_papers/E39.pdf) · [text](extracted_text/field_papers/E39.txt) |
| E40 | Directly relevant | [Model-based variable impedance learning control for robotic manipulation](papers/field_papers/E40.pdf) · [text](extracted_text/field_papers/E40.txt) |
| E41 | Directly relevant | [Model Predictive Variable Impedance Control of Manipulators for Adaptive Precision-Compliance Tradeoff](papers/field_papers/E41.pdf) · [text](extracted_text/field_papers/E41.txt) |
| E43 | Directly relevant | [Neural-Network-Based Optimal Impedance Control for Robots in Physical Interaction With Soft Environments](papers/NNBO.pdf) · [text](extracted_text/nnbo.txt) |
| E44 | Directly relevant | [Optimal Impedance Force-Tracking Control Design With Impact Formulation for Interaction Tasks](papers/field_papers/E44.pdf) · [text](extracted_text/field_papers/E44.txt) |
| E45 | Directly relevant | [Optimized Impedance Adaptation of Robot Manipulator Interacting With Unknown Environment](papers/field_papers/E45.pdf) · [text](extracted_text/field_papers/E45.txt) |
| E46 | Partly relevant | [Passive Interaction Control With Dynamical Systems](papers/field_papers/E46.pdf) · [text](extracted_text/field_papers/E46.txt) |
| E47 | Directly relevant | [Passive Model-Predictive Impedance Control for Safe Physical HumanRobot Interaction](papers/field_papers/E47.pdf) · [text](extracted_text/field_papers/E47.txt) |
| E48 | Partly relevant | [Passivity-Based Skill Motion Learning in Stiffness-Adaptive Unified Force-Impedance Control](papers/field_papers/E48.pdf) · [text](extracted_text/field_papers/E48.txt) |
| E49 | Directly relevant | [Passivity Filter for Variable Impedance Control](papers/field_papers/E49.pdf) · [text](extracted_text/field_papers/E49.txt) |
| E54 | Directly relevant | [Robot Impedance Iterative Learning with Sparse Online Gaussian Process(1)](papers/field_papers/E54.pdf) · [text](extracted_text/field_papers/E54.txt) |
| E57 | Directly relevant | [Unified Motion Force Impedance Control for Manipulators in Unknown Contact Environments Based on Robust Model-Reaching Approach](papers/field_papers/E57.pdf) · [text](extracted_text/field_papers/E57.txt) |
| E58 | Directly relevant | [Variable Impedance Control - A Reinforcement Learning Approach](papers/field_papers/E58.pdf) · [text](extracted_text/field_papers/E58.txt) |
| E59 | Background review | [Variable impedance models including fuzzy fractional order for control of human–robot interaction a systematic review](papers/field_papers/E59.pdf) · [text](extracted_text/field_papers/E59.txt) |

## Originals added from existing leads

| ID | PDF and text |
|---|---|
| F01 | [C-GMS, arXiv:2511.16330v2, 2026-03-01](papers/field_supplements/F01.pdf) · [text](extracted_text/field_supplements/F01.txt) |
| F02 | [Duan et al., RAS2018, 10.1016/j.robot.2018.01.009](papers/field_supplements/F02.pdf) · [text](extracted_text/field_supplements/F02.txt) |
| F03 | [Ge et al., IJC2014, 10.1080/00207179.2013.827799](papers/field_supplements/F03.pdf) · [text](extracted_text/field_supplements/F03.txt) |
| F04 | [Matschek et al., TCST2023, 10.1109/TCST.2023.3262987](papers/field_supplements/F04.pdf) · [text](extracted_text/field_supplements/F04.txt) |
| F05 | [Tiezzi and Kao, TRO2007, 10.1109/TRO.2006.889494](papers/field_supplements/F05.pdf) · [text](extracted_text/field_supplements/F05.txt) |
| F06 | [Jin et al., TMECH2025, 10.1109/TMECH.2024.3466950](papers/field_supplements/F06.pdf) · [text](extracted_text/field_supplements/F06.txt) |
| F07 | [Zhang et al., LRA2025, 10.1109/LRA.2025.3604746](papers/field_supplements/F07.pdf) · [text](extracted_text/field_supplements/F07.txt) |
| F08 | [Duan et al., AuSoScan, TMECH2025, 10.1109/TMECH.2025.3583041](papers/field_supplements/F08.pdf) · [text](extracted_text/field_supplements/F08.txt) |
| F09 | [Wang et al., TASE2024, 10.1109/TASE.2023.3282974](papers/field_supplements/F09.pdf) · [text](extracted_text/field_supplements/F09.txt) |
| F10 | [Moreira et al., CMPB2014, 10.1016/j.cmpb.2014.01.017](papers/field_supplements/F10.pdf) · [text](extracted_text/field_supplements/F10.txt) |

## Independent reverse-analysis samples

Samples were not selected by relevance to this idea. The ten below had relevant passages examined; they are neither all journal papers nor complete theoretical audits. They inform question/method/evidence design; see [independent reverse analysis](../research/independent_reverse_analysis_of_top_journal_methods.md).

| ID | PDF and text |
|---|---|
| R01 | [A DMPs-Based Framework for Robot Learning and Generalization of Humanlike Variable Impedance Skills](papers/reverse_analysis_samples/R01.pdf) · [text](extracted_text/reverse_analysis_samples/R01.txt) |
| R03 | [Adaptive Parameter Estimation and Control Design for Robot Manipulators With Finite-Time Convergence](papers/reverse_analysis_samples/R03.pdf) · [text](extracted_text/reverse_analysis_samples/R03.txt) |
| R05 | [Asymmetric Bimanual Control of Dual-Arm Exoskeletons for Human-Cooperative Manipulations](papers/reverse_analysis_samples/R05.pdf) · [text](extracted_text/reverse_analysis_samples/R05.txt) |
| R11 | [Human-Like Adaptation of Force and Impedance in Stable and Unstable Interactions](papers/reverse_analysis_samples/R11.pdf) · [text](extracted_text/reverse_analysis_samples/R11.txt) |
| R12 | [Impedance Learning for Human-Guided Robots in Contact With Unknown Environments](papers/field_papers/E28.pdf) · [text](extracted_text/field_papers/E28.txt) |
| R13 | [Learning a Flexible Neural Energy Function With a Unique Minimum for Globally Stable and Accurate Demonstration Learning](papers/reverse_analysis_samples/R13.pdf) · [text](extracted_text/reverse_analysis_samples/R13.txt) |
| R15 | [Multi-hierarchy interaction control of a redundant robot using impedance learning](papers/reverse_analysis_samples/R15.pdf) · [text](extracted_text/reverse_analysis_samples/R15.txt) |
| R17 | [Neural Control of Robot Manipulators With Trajectory Tracking Constraints and Input Saturation](papers/reverse_analysis_samples/R17.pdf) · [text](extracted_text/reverse_analysis_samples/R17.txt) |
| R19 | [Residual Reinforcement Learning for Robot Control](papers/reverse_analysis_samples/R19.pdf) · [text](extracted_text/reverse_analysis_samples/R19.txt) |
| R21 | [Stability Criterion and Stability Enhancement for a Thruster-Assisted Underwater Hexapod Robot](papers/reverse_analysis_samples/R21.pdf) · [text](extracted_text/reverse_analysis_samples/R21.txt) |

## Candidate 2 source and access limits

[You Only Plan Once: A Learning-Based One-Stage Planner With Guidance Learning](papers/method_transfer/YOPO.pdf) · [text](extracted_text/method_transfer/YOPO.txt). The original uses guidance learning, not PPO. Text is local and no longer depends on the removed PPO report.

VICES (arXiv:1906.08880), SRL-VIC (arXiv:2406.13744) and Zanella et al. (arXiv:2301.12759) have prior records of relevant original-passage access, but no local PDF was obtained in that consolidation. Their roles/limits are in transferable mechanisms/field state. Yang et al., IFAC2025 (DOI10.1016/j.ifacol.2025.10.191) was checked at abstract level only. Recover earlier online-reading details from the old candidate-2 report at `db28365`; do not register them as local full texts.

Two TIE2026 papers, the original Roveda2022 journal article and recent admittance/joint-estimation neighbors still have full-text or local-acquisition gaps; see field state/source JSON. Missing full text does not imply that a relevant difficulty remains unsolved.

The unchanged [sample mathematical-boundary results](top_venue_sample_math_checks.json) and [historical verification script](verify_top_venue_samples.py) examine generic constructions only and cannot refute entire papers. Those scientific checks were not run in the consolidation.
