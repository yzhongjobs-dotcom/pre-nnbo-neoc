# 文献索引

本目录保存可服务后续研究的原文、检索文本、来源及实际阅读范围。领域现状、方法借鉴、独立顶刊样本分开使用；下载和文本提取不等于读完全文或复证理论。

## 使用与来源

- 科学判断入口：[研究现状](../research/field_state.md)、[候选方法](../research/optimal_impedance_question_and_methods.md)、[可借鉴机制](../research/transferable_mechanisms.md)。
- 本地文件与SHA256：[本地资料清单](local_materials_manifest.json)。每项有原文/文本路径、用途、来源URL、已知版本及阅读范围；同一文件去重共享，不重复复制。
- 指定仓库一：[academic-writing-skill/文献资料](https://github.com/yzhongjobs-dotcom/academic-writing-skill/tree/c88d48eb08a876eb9b81a30d37b62d70fa028581/文献资料)，22篇身份/摘要、其中10篇关键段落分析；本地保存这10篇样本。其余12篇元数据留在[样本来源](top_venue_sample_sources_and_reading_scope.json)，不冒充已读全文。
- 指定仓库二：[e2e/User Batch V2](https://github.com/yzhongjobs-dotcom/e2e/tree/8fff5f3cbb1d47405f54888c0e92c61d08e31ced/source-materials/user-batch-v2)及[idea-stage](https://github.com/yzhongjobs-dotcom/e2e/tree/8fff5f3cbb1d47405f54888c0e92c61d08e31ced/idea-stage)。原61篇筛查中的25直接相关、6相关部分、5背景资料落地；25篇当前不相关的只保留排除理由。由idea-stage等既有线索补取的10份原文也落地。
- 完整筛选、原文页码和未获取全文的近邻见[来源与筛选](impedance_literature_sources_and_screening.json)。2653条发现记录、105张旧证据卡和3822条检索事件不整库复制，也不作为已确认科学结论。

本地共 **65** 份去重PDF，另有相应文本。资料ID可用于 `python literature/verify_materials.py --ids E27 F08 YOPO`；检查全库用 `python literature/verify_materials.py`。`--download`只恢复缺失的已选文件并核哈希，不重新搜索/筛选论文。文本提取依赖系统`pdftotext`；不同版本输出有差异时需明确记录，不覆盖原文身份哈希。

## 原有方法与直接依据

| ID | 原文与文本 |
|---|---|

| beber | [beber](papers/beber.pdf) · [文本](extracted_text/beber.txt) |
| pollayil | [pollayil](papers/pollayil.pdf) · [文本](extracted_text/pollayil.txt) |
| wijayarathne | [wijayarathne](papers/wijayarathne.pdf) · [文本](extracted_text/wijayarathne.txt) |
| haninger | [haninger](papers/haninger.pdf) · [文本](extracted_text/haninger.txt) |
| xue | [xue](papers/xue.pdf) · [文本](extracted_text/xue.txt) |
| risk_sensitive | [risk sensitive](papers/risk_sensitive.pdf) · [文本](extracted_text/risk_sensitive.txt) |
| fu | [fu](papers/fu.pdf) · [文本](extracted_text/fu.txt) |
| piedra_volumetric | [piedra volumetric](papers/piedra_volumetric.pdf) · [文本](extracted_text/piedra_volumetric.txt) |
| NEOC | [NEOC](papers/NEOC.pdf) · [文本](extracted_text/neoc.txt) |
| implicit_differentiation | [implicit differentiation](papers/implicit_differentiation.pdf) · [文本](extracted_text/idoc.txt) |
| NNBO | [NNBO](papers/NNBO.pdf) · [文本](extracted_text/nnbo.txt) |

## 领域相关与背景资料

纳入理由和核读页码以来源JSON为准；“直接相关”不等于已完成全文核读。E30复用已有risk_sensitive版本，外部E30原文件哈希与本地版本分别标记。

| ID | 用途 | 原文与文本 |
|---|---|---|
| E01 | 背景综述 | [A Tutorial Survey and Comparison of Impedance Control on Robotic Manipulation](papers/field_papers/E01.pdf) · [文本](extracted_text/field_papers/E01.txt) |
| E04 | 直接相关 | [A Novel Robust Adaptive Impedance Control Scheme for Contact Force Tracking With Unknown Environment](papers/field_papers/E04.pdf) · [文本](extracted_text/field_papers/E04.txt) |
| E06 | 直接相关 | [Adaptive Energy-Bounding Approach for Robustly Stable Interaction Control of Impedance-Controlled Industrial Robot With Uncertain Environments](papers/field_papers/E06.pdf) · [文本](extracted_text/field_papers/E06.txt) |
| E09 | 直接相关 | [Adaptive Impedance Control of HumanRobot Cooperation Using Reinforcement Learning](papers/field_papers/E09.pdf) · [文本](extracted_text/field_papers/E09.txt) |
| E11 | 相关部分 | [An Optimal Variable Impedance Control With Consideration of the Stability](papers/field_papers/E11.pdf) · [文本](extracted_text/field_papers/E11.txt) |
| E12 | 背景综述 | [Applied Bionics and Biomechanics - 2018 - Al-Shuka - Active Impedance Control of Bioinspired Motion Robotic Manipulators ](papers/field_papers/E12.pdf) · [文本](extracted_text/field_papers/E12.txt) |
| E14 | 直接相关 | [Constrained Variable Impedance Control using Quadratic Programming](papers/field_papers/E14.pdf) · [文本](extracted_text/field_papers/E14.txt) |
| E15 | 相关部分 | [Deflection-Domain Passivity Control of Variable Stiffnesses Based on Potential Energy Reference](papers/field_papers/E15.pdf) · [文本](extracted_text/field_papers/E15.txt) |
| E16 | 直接相关 | [Development of a Safety- and Energy-Aware Impedance Controller for Collaborative Robots](papers/field_papers/E16.pdf) · [文本](extracted_text/field_papers/E16.txt) |
| E17 | 相关部分 | [Diffusion-Based Impedance Learning for Contact-Rich Manipulation Tasks](papers/field_papers/E17.pdf) · [文本](extracted_text/field_papers/E17.txt) |
| E18 | 直接相关 | [Dynamic Adaptive Hybrid Impedance Control for Dynamic Contact Force Tracking in Uncertain Environments](papers/field_papers/E18.pdf) · [文本](extracted_text/field_papers/E18.txt) |
| E23 | 直接相关 | [Force tracking impedance control of robot manipulators under unknown environment](papers/field_papers/E23.pdf) · [文本](extracted_text/field_papers/E23.txt) |
| E24 | 直接相关 | [Human-Inspired Adaptive Optimal Control Framework for Robot-Environment Interaction](papers/field_papers/E24.pdf) · [文本](extracted_text/field_papers/E24.txt) |
| E25 | 背景综述 | [Impedance and Interaction Control](papers/field_papers/E25.pdf) · [文本](extracted_text/field_papers/E25.txt) |
| E27 | 直接相关 | [Impedance Learning-Based Adaptive Force Tracking for Robot on Unknown Terrains](papers/field_papers/E27.pdf) · [文本](extracted_text/field_papers/E27.txt) |
| E28 | 直接相关 | [Impedance Learning for Human-Guided Robots in Contact With Unknown Environments](papers/field_papers/E28.pdf) · [文本](extracted_text/field_papers/E28.txt) |
| E29 | 直接相关 | [Impedance Learning for Robotic Contact Tasks Using Natural Actor-Critic Algorithm](papers/field_papers/E29.pdf) · [文本](extracted_text/field_papers/E29.txt) |
| E30 | 直接相关 | [Impedance Optimization for Uncertain Contact Interactions Through Risk Sensitive Optimal Control](papers/risk_sensitive.pdf) · [文本](extracted_text/risk_sensitive.txt) |
| E32 | 背景综述 | [Impedance Variation and Learning Strategies in HumanRobot Interaction](papers/field_papers/E32.pdf) · [文本](extracted_text/field_papers/E32.txt) |
| E34 | 直接相关 | [Learning-Based Predictive Impedance Control Towards Safe Predefined-Time Physical Robotic Interaction](papers/field_papers/E34.pdf) · [文本](extracted_text/field_papers/E34.txt) |
| E37 | 直接相关 | [Model-Based Reinforcement Learning Variable Impedance Control for Human-Robot Collaboration](papers/field_papers/E37.pdf) · [文本](extracted_text/field_papers/E37.txt) |
| E38 | 直接相关 | [Model-Based ActorCritic Learning of Robotic Impedance Control in Complex Interactive Environment](papers/field_papers/E38.pdf) · [文本](extracted_text/field_papers/E38.txt) |
| E39 | 相关部分 | [Model-Based Predictive Impedance Variation for Obstacle Avoidance in Safe HumanRobot Collaboration](papers/field_papers/E39.pdf) · [文本](extracted_text/field_papers/E39.txt) |
| E40 | 直接相关 | [Model-based variable impedance learning control for robotic manipulation](papers/field_papers/E40.pdf) · [文本](extracted_text/field_papers/E40.txt) |
| E41 | 直接相关 | [Model Predictive Variable Impedance Control of Manipulators for Adaptive Precision-Compliance Tradeoff](papers/field_papers/E41.pdf) · [文本](extracted_text/field_papers/E41.txt) |
| E43 | 直接相关 | [Neural-Network-Based Optimal Impedance Control for Robots in Physical Interaction With Soft Environments](papers/NNBO.pdf) · [文本](extracted_text/nnbo.txt) |
| E44 | 直接相关 | [Optimal Impedance Force-Tracking Control Design With Impact Formulation for Interaction Tasks](papers/field_papers/E44.pdf) · [文本](extracted_text/field_papers/E44.txt) |
| E45 | 直接相关 | [Optimized Impedance Adaptation of Robot Manipulator Interacting With Unknown Environment](papers/field_papers/E45.pdf) · [文本](extracted_text/field_papers/E45.txt) |
| E46 | 相关部分 | [Passive Interaction Control With Dynamical Systems](papers/field_papers/E46.pdf) · [文本](extracted_text/field_papers/E46.txt) |
| E47 | 直接相关 | [Passive Model-Predictive Impedance Control for Safe Physical HumanRobot Interaction](papers/field_papers/E47.pdf) · [文本](extracted_text/field_papers/E47.txt) |
| E48 | 相关部分 | [Passivity-Based Skill Motion Learning in Stiffness-Adaptive Unified Force-Impedance Control](papers/field_papers/E48.pdf) · [文本](extracted_text/field_papers/E48.txt) |
| E49 | 直接相关 | [Passivity Filter for Variable Impedance Control](papers/field_papers/E49.pdf) · [文本](extracted_text/field_papers/E49.txt) |
| E54 | 直接相关 | [Robot Impedance Iterative Learning with Sparse Online Gaussian Process(1)](papers/field_papers/E54.pdf) · [文本](extracted_text/field_papers/E54.txt) |
| E57 | 直接相关 | [Unified Motion Force Impedance Control for Manipulators in Unknown Contact Environments Based on Robust Model-Reaching Approach](papers/field_papers/E57.pdf) · [文本](extracted_text/field_papers/E57.txt) |
| E58 | 直接相关 | [Variable Impedance Control - A Reinforcement Learning Approach](papers/field_papers/E58.pdf) · [文本](extracted_text/field_papers/E58.txt) |
| E59 | 背景综述 | [Variable impedance models including fuzzy fractional order for control of human–robot interaction a systematic review](papers/field_papers/E59.pdf) · [文本](extracted_text/field_papers/E59.txt) |

## 从既有线索补入的原文

| ID | 原文与文本 |
|---|---|
| F01 | [C-GMS, arXiv:2511.16330v2, 2026-03-01](papers/field_supplements/F01.pdf) · [文本](extracted_text/field_supplements/F01.txt) |
| F02 | [Duan et al., RAS2018, 10.1016/j.robot.2018.01.009](papers/field_supplements/F02.pdf) · [文本](extracted_text/field_supplements/F02.txt) |
| F03 | [Ge et al., IJC2014, 10.1080/00207179.2013.827799](papers/field_supplements/F03.pdf) · [文本](extracted_text/field_supplements/F03.txt) |
| F04 | [Matschek et al., TCST2023, 10.1109/TCST.2023.3262987](papers/field_supplements/F04.pdf) · [文本](extracted_text/field_supplements/F04.txt) |
| F05 | [Tiezzi and Kao, TRO2007, 10.1109/TRO.2006.889494](papers/field_supplements/F05.pdf) · [文本](extracted_text/field_supplements/F05.txt) |
| F06 | [Jin et al., TMECH2025, 10.1109/TMECH.2024.3466950](papers/field_supplements/F06.pdf) · [文本](extracted_text/field_supplements/F06.txt) |
| F07 | [Zhang et al., LRA2025, 10.1109/LRA.2025.3604746](papers/field_supplements/F07.pdf) · [文本](extracted_text/field_supplements/F07.txt) |
| F08 | [Duan et al., AuSoScan, TMECH2025, 10.1109/TMECH.2025.3583041](papers/field_supplements/F08.pdf) · [文本](extracted_text/field_supplements/F08.txt) |
| F09 | [Wang et al., TASE2024, 10.1109/TASE.2023.3282974](papers/field_supplements/F09.pdf) · [文本](extracted_text/field_supplements/F09.txt) |
| F10 | [Moreira et al., CMPB2014, 10.1016/j.cmpb.2014.01.017](papers/field_supplements/F10.pdf) · [文本](extracted_text/field_supplements/F10.txt) |

## 独立逆向分析样本

样本不按本idea相关性挑选。下列10篇此前核读过相关段落，不代表全部都是顶刊，也不代表完整理论审计。其作用是理解问题/方法/证据设计，分析见[独立逆向分析](../research/independent_reverse_analysis_of_top_journal_methods.md)。

| ID | 原文与文本 |
|---|---|
| R01 | [A DMPs-Based Framework for Robot Learning and Generalization of Humanlike Variable Impedance Skills](papers/reverse_analysis_samples/R01.pdf) · [文本](extracted_text/reverse_analysis_samples/R01.txt) |
| R03 | [Adaptive Parameter Estimation and Control Design for Robot Manipulators With Finite-Time Convergence](papers/reverse_analysis_samples/R03.pdf) · [文本](extracted_text/reverse_analysis_samples/R03.txt) |
| R05 | [Asymmetric Bimanual Control of Dual-Arm Exoskeletons for Human-Cooperative Manipulations](papers/reverse_analysis_samples/R05.pdf) · [文本](extracted_text/reverse_analysis_samples/R05.txt) |
| R11 | [Human-Like Adaptation of Force and Impedance in Stable and Unstable Interactions](papers/reverse_analysis_samples/R11.pdf) · [文本](extracted_text/reverse_analysis_samples/R11.txt) |
| R12 | [Impedance Learning for Human-Guided Robots in Contact With Unknown Environments](papers/field_papers/E28.pdf) · [文本](extracted_text/field_papers/E28.txt) |
| R13 | [Learning a Flexible Neural Energy Function With a Unique Minimum for Globally Stable and Accurate Demonstration Learning](papers/reverse_analysis_samples/R13.pdf) · [文本](extracted_text/reverse_analysis_samples/R13.txt) |
| R15 | [Multi-hierarchy interaction control of a redundant robot using impedance learning](papers/reverse_analysis_samples/R15.pdf) · [文本](extracted_text/reverse_analysis_samples/R15.txt) |
| R17 | [Neural Control of Robot Manipulators With Trajectory Tracking Constraints and Input Saturation](papers/reverse_analysis_samples/R17.pdf) · [文本](extracted_text/reverse_analysis_samples/R17.txt) |
| R19 | [Residual Reinforcement Learning for Robot Control](papers/reverse_analysis_samples/R19.pdf) · [文本](extracted_text/reverse_analysis_samples/R19.txt) |
| R21 | [Stability Criterion and Stability Enhancement for a Thruster-Assisted Underwater Hexapod Robot](papers/reverse_analysis_samples/R21.pdf) · [文本](extracted_text/reverse_analysis_samples/R21.txt) |

## 候选2的借鉴原文及访问边界

[You Only Plan Once: A Learning-Based One-Stage Planner With Guidance Learning](papers/method_transfer/YOPO.pdf) · [文本](extracted_text/method_transfer/YOPO.txt)。原文采用guidance learning；不是PPO。本文本已落地，不再依赖被删除的旧PPO报告。

VICES（arXiv:1906.08880）、SRL-VIC（arXiv:2406.13744）、Zanella等（arXiv:2301.12759）此前记录过原文相关段落访问，本次未取得本地PDF；其作用和边界已提炼至可借鉴机制及研究现状。Yang等IFAC2025（DOI10.1016/j.ifacol.2025.10.191）仅核摘要。原在线阅读细节可由清理前提交`db28365`中的旧候选2报告恢复；不把它们登记为本地全文。

两篇TIE2026、Roveda2022原刊、近期导纳/联合估计近邻仍有全文或本地获取缺口，具体记录见研究现状与来源JSON。没有全文不能推定没有解决相关困难。

原样保留的[样本数学边界结果](top_venue_sample_math_checks.json)与[历史复核脚本](verify_top_venue_samples.py)只检查一般构造，不能据此否定整篇论文；本轮未运行这些科学核查。
