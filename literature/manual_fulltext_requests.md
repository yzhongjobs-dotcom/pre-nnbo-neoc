# 人工原文获取：8篇待下载（ARIS近邻补查，2026-10-31）

更新：2026-10-31。ARIS `research-lit`＋`novelty-check`补查（OpenAlex/S2，24条查询，35条候选）发现：“在线刚度＋力反馈补偿未知生理运动”“无视觉的力预测补偿”“IMM生理运动预测”均已有先例，研究现状的缺口表述必须收窄到这些工作之外的具体增量。以下原文用于逐篇界定增量，请下载到 `literature/manual_inbox/`（F10 Moreira CMPB2014 已在本地并已核读）：

| 优先 | 编号 | 文献 | DOI | 需要核对 |
|---|---|---|---|---|
| 1 | L02 | Cortesão, Dominici, "Robot Force Control on a Beating Heart," IEEE/ASME TMECH, 2017 | 10.1109/TMECH.2017.2696259 | 双AOB观测对象（位置/速度还是等效扰动）、刚度是否更新、执行接口 |
| 1配套 | L03 | Correction to "Robot Force Control on a Beating Heart," IEEE/ASME TMECH, 2017 | 10.1109/TMECH.2017.2765798 | 勘误是否影响机制或稳定性 |
| 2 | L04 | Dominici, Cortesão, "Cascade force control for autonomous beating heart motion compensation," Control Eng. Practice, 2015 | 10.1016/j.conengprac.2014.12.012 | 由力到运动参考的状态、对刚度的依赖 |
| 3 | L06 | Yuan, Lee, "Robot-Assisted Ultrasound Probe Force Control Under Respiration-Induced Motion," IEEE AIM, 2024 | 10.1109/AIM55361.2024.10637095 | 输入信号、准周期预测、刚度与接口 |
| 4 | L35 | Beber et al., "Towards Robotised Palpation for Cancer Detection through Online Tissue Viscoelastic Characterisation with a Collaborative Robotic Arm," IROS, 2024（arXiv:2404.09542） | 10.1109/IROS58592.2024.10802582 | 压入量与黏弹参数如何联合可观、是否允许组织运动 |
| 5 | L05 | Moreira, Abayazid, Misra, "Towards physiological motion compensation for flexible needle interventions," IROS, 2015 | 10.1109/IROS.2015.7353468 | 由力换算运动的方法、标定、额外触点 |
| 6 | L07 | Liang, Xie, Yu, "Beating Heart Motion Accurate Prediction Method Based on Interactive Multiple Model: An Information Fusion Approach," BioMed Res. Int., 2017 | 10.1155/2017/1279486 | IMM的观测来源、模式与噪声结构 |
| 7 | L16 | Yao et al., "Safe Variable-Force Control for Robotic Ultrasound Scanning via ℓ2-Regularized Fuzzy Broad Offline Reinforcement Learning," IEEE TMRB, 2026 | 10.1109/TMRB.2026.3699233 | 在体呼吸条件、状态输入与控制接口 |

第二批（按论证需要再下）：L08 Wang等MSSP2021（10.1016/j.ymssp.2021.107819，联合阻抗与接触位置估计）；L31 Tsuji IROS2022（10.1109/IROS47612.2022.9981146，概率可靠度门控刚度辨识）。调研报告与全部35条候选见[ARIS补查报告](aris_concurrent_survey_2026-10-31.md)，要点将并入[研究现状](../research/field_state.md)。
