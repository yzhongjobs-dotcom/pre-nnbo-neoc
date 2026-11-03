# 交接入口：运动且软硬不均组织上的位置接口恒力超声扫描

更新于 2026-11-03（`git log -1 -- project_records/HANDOFF.md` 给出提交）。**用户最新消息优先于这里的一切。** 旧版交接内容：`git show c64be05:project_records/HANDOFF.md`。

## 现在的情况

**用户最新交流状态：先讨论模块，不继续自行跑模型变体。** 用户要求记录“模块锚定”设计偏差以供后续优化skill，随后明确先讲清微正弦激励、漂移投影与刚度回归，以及Rossi等已有路线。事实、执行/skill初步归因见[失败与纠偏](research_failures_and_corrections.md)；本轮只澄清与记录，未修改skill、控制器或仿真。后续依据用户讨论结果继续，不以旧待办自动启动动态逆实现。

**2026-11-03 最新：保留精简候选，继续完成了接触曲率/记忆的实现与小规模闭环。** 先读[第一模块 §8–10](../research/first_module_information.md)。端点曲率混入整体运动后已修移动参考、当前增益功能量和正曲率信息筛选；硬组织 .266780→.265463 N但仅两次启用，T3m完全回退，无充分采纳证据。正约束局部SLS改善硬组织条件预测 .1544→.1257 N，但只接探测实部增益后闭环退化到.28631 N；指数漂移消融排除了一部分替代解释。较慢 ke 的条件不确定性明显大于20.8Hz响应，不可混称同一刚度。全部12个闭环、同输入正负和数值检查已留档。**下一项仍在接触模块内：** 对动态柔顺观测的逆功能量、初值和参数更新一致性作物理/信息闭合；正文已给冻结参数逆式但尚未实现。不能把本次有限标量接口失败说成SLS/接触路线无用，不能把弱曲率差值包装成完成顶刊升级。精简候选和原正式方法均保留，仿真未改。

**2026-11-03 当前用户授权与实质进展：第一个模块联合优化（主动激励、接触表示、估计）。** 已实现滤波对齐的窗口末端有效增益＋变化率估计、条件信息探测设计及物理力波动预算。小规模定向闭环：S2 .13425→.11498 N，T3m .12114→.11483 N，硬推拉 .26805→.26678 N（近似持平）。只延长常数窗反而变差；较强精简候选保留原20.8 Hz探测和全部上层控制参数，不强求速率项。噪声实现与未知伺服滞后各一组配对也未见明显退化。细节/正负结果/理论边界先读 [第一模块研究](../research/first_module_information.md)，代码与复现 [contact_information](../experiments/contact_information/README.md)。已推检查点，不是只有计划；尚未正式采纳或完成顶刊贡献判断。下一项在同一候选上深入模型剩余、闭环信息有效性和最接近工作区别，不默认扩大任务矩阵或另开 MPC 专项。

**2026-11-03 用户规则纠正（持续有效）：** 后续仿真按“每场景的科学问题—每组比较排除的替代解释—整条证据链”控制规模，取消统一每族 50＋剂量扫描、全面参数/分布外矩阵的默认要求；solid 与包装同时重视。转贴分析核查见接触专项 §8。该核查之后，用户已明确授权并开始上述第一个模块联合优化；不预定三个部件全部替换或必须选 MHE/GLS。

**2026-11-02 接触专题阶段（现已扩展为上述联合优化）：** 已精读用户给定 AuSoScan、T-RO 2023 原文，完成独立有限足迹/柱面候选的因果回放、同记忆消融与局部干扰检查。结果、假设、修正和未解项见 [接触模型专项研究](../research/contact_model_study.md)，复现见 [实验说明](../experiments/contact_model/README.md)。这些复杂接触候选未有一致优势，未采纳，也未据此否定接触路线。用户六项追问仍适用：依据、合理性、物理数学、仿真真实性、因果验证和失败修正。

**后续明确与诊断：** 借鉴前先确认论文真实痛点是否在本项目存在；仿真是否调整由物理需求决定。专项研究 §6 已补上论文—任务—仿真对应：当前 Ft 不反作用于位置/法向力，不能检验摩擦负载补偿的那部分收益；静态逆压入原型还混用了总力与动态增益，该项在 §7 已修正，尚无一致任务优势。现阶段保留对象，不预定最终模型或另开模块。原始诊断在 `experiments/contact_model/results/paper_painpoint_sensitivity.json` 与 `normal_contact_state_diagnostic.json`。

**2026-11-03 接触模型续研：** 已修动态状态/目标力关系、前向力求解、独立回退和到达信息；完成五族静态/动态/柱面闭环，定位并修正剖面导数和累计观测参数方差。完整正负结果见专项研究 §7.5–7.10；汇总 `experiments/contact_model/results/dynamic_contact_five_family_summary.json`。没有一致改进，未采纳候选；低松弛假设在两场景有收益但参数未辨识，拟合仍数百毫秒至秒。**下一项：** 同输入回放比较松弛/空间假设的条件预测和控制功能量，区分结构偏差、可辨识性与反馈分布；若缺独立信息，再具体设计标定/压力矩接触候选，不继续靠调 R/门限补结果。只研究接触模型，不改仿真求胜、不换其他专项；本轮 Fn-only 负面不裁决未完成的动态力矩/空间预瞄。

- **任务：** 位置接口协作臂只用力和编码器，在运动、软硬不均的组织上保持恒定接触力。
- **目标：** 把已胜过多个基线的 Proposed Method 优化到 T-RO/IJRR 水准，work 第一，方法与 story 同样重要。
- **可靠主体：** 最终方法 `Proposed-IMMb` 加 V1 伺服在线辨识（V1 的失配恢复已验证）。
- **未完成：** 用户要求的接触模型、数学投影/MHE、GN/MPC 架构和整体方法升级。V1 之后的候选（`3b81615`、`c64be05`）都未并入方法，论文初稿也未修改。
- **最近的教训：** 反复出现“诊断出原因却不修原候选，转向下一机制”。开始前务必看第 3 步的偏差模式。

## 阅读顺序

分层读，每步读完能回答对应问题再往下走。不要一开始就读 methods 全文（634 行）或 Git 全历史。

| 步 | 读什么 | 要能回答 | 约用时 |
|---|---|---|---|
| 1 | [handoff/01_要求与标准](handoff/01_要求与标准.md) | 用户的五条要求是什么？哪些历史方案已失效？有哪些待确认冲突？ | 短 |
| 2 | Skill `top-researcher`：`README.md`、`references/stages/README.md`，再按当前阶段读 `workflow/03_method.md` 和 `references/stages/03_method.md` | 判据是什么？跨阶段闭环怎样检查？ | 中 |
| 3 | [handoff/02_研究主线](handoff/02_研究主线.md)，并按其指引读 [论文初稿](../research/paper_draft_zh.md) §一至§四、§五.1、§五.3、§五.9 | 现状 → 问题 → 方法 → 证据这条线是什么？代码在哪？ | 中 |
| 4 | [handoff/03_方法演进审计](handoff/03_方法演进审计.md) | 每项尝试实际做到哪一层？失败原因是什么？反复犯了哪些错？ | 中 |
| 5 | [handoff/04_当前状态与待解决问题](handoff/04_当前状态与待解决问题.md) | 什么已验证、什么没研究透、什么可放下、当前要解决什么？ | 短 |
| 6 | 只在需要时：[methods 记录](../research/optimal_impedance_question_and_methods.md) §16、§18–21；[B0 设计](../research/method_design.md) §3；[决策记录](decision_log.md)；[agent_runs](agent_runs/)；[实验索引](../experiments/experiment_index.md) | 原始推理与证据 | 按需 |

读完后先向用户说明你的理解和下一步打算，再开始实质工作。本交接**不规定最终技术路线**。

## 下一项决定属于谁

- **题内技术探索：** 接手者依据 skill 主动研究，大规模实现前先报告设计及理由。
- **研究问题的实质修订、长期约束的变更、[01 第 6 节](handoff/01_要求与标准.md)的冲突：** 由用户决定。

## 运行与规则

- **环境：** Windows，`python -X utf8`，Python 3.10，numpy、scipy、numba。
- **同步：** 工作前 `git pull --rebase origin main`；每个完成的检查点运行 `python -X utf8 scripts/research_sync.py --message "research: <具体改动与目的>"`。不改写历史，不强推。
- **记录：** Codex/ARIS 的任务书与最终输出存 [agent_runs](agent_runs/)，关键决定写入 [decision_log](decision_log.md)，科学推理放研究文档。不另建变更日志。
- **保护：** `00_基线快照_2026-10-14/` 逐字节不动（`python -X utf8 literature/verify_materials.py`）。不翻译已有材料。
- **实验目录：** [motion_aware_scan](../experiments/motion_aware_scan/README.md)（最终方法）、[method_redesign](../experiments/method_redesign/)（11-01 的 P1、P2）、[contact_model](../experiments/contact_model/)（11-01 傍晚起）、[cobot_scan_benchmark](../experiments/cobot_scan_benchmark/README.md)（对象与 C1）。
- **文献：** 索引 [literature_index](../literature/literature_index.md)；ARIS 运行 `recent-related-work-20261031`（RUN_ROOT `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run`，机器本地）；当前没有待下载文献（[清单](../literature/manual_fulltext_requests.md)）；用户下载件放 `literature/manual_inbox/`。
- **更早的实验和文档：** 早期 B0/二维实验已移出工作区，用 `git show b7a3f91:<路径>` 恢复。
