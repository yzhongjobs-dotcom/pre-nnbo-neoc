# 交接入口：运动且软硬不均组织上的位置接口恒力超声扫描

更新于 2026-11-02（`git log -1 -- project_records/HANDOFF.md` 给出提交）。**用户最新消息优先于这里的一切。** 旧版交接内容：`git show c64be05:project_records/HANDOFF.md`。

## 现在的情况

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

- **环境：** Windows，`python -X utf8`，Python 3.10，numpy、scipy、numba，8 线程以上为宜。每族 50 的 5 个配置约 30 分钟本机 CPU。
- **同步：** 工作前 `git pull --rebase origin main`；每个完成的检查点运行 `python -X utf8 scripts/research_sync.py --message "research: <具体改动与目的>"`。不改写历史，不强推。
- **记录：** Codex/ARIS 的任务书与最终输出存 [agent_runs](agent_runs/)，关键决定写入 [decision_log](decision_log.md)，科学推理放研究文档。不另建变更日志。
- **保护：** `00_基线快照_2026-10-14/` 逐字节不动（`python -X utf8 literature/verify_materials.py`）。不翻译已有材料。
- **实验目录：** [motion_aware_scan](../experiments/motion_aware_scan/README.md)（最终方法）、[method_redesign](../experiments/method_redesign/)（11-01 的 P1、P2）、[contact_model](../experiments/contact_model/)（11-01 傍晚起）、[cobot_scan_benchmark](../experiments/cobot_scan_benchmark/README.md)（对象与 C1）。
- **文献：** 索引 [literature_index](../literature/literature_index.md)；ARIS 运行 `recent-related-work-20261031`（RUN_ROOT `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run`，机器本地）；当前没有待下载文献（[清单](../literature/manual_fulltext_requests.md)）；用户下载件放 `literature/manual_inbox/`。
- **更早的实验和文档：** 早期 B0/二维实验已移出工作区，用 `git show b7a3f91:<路径>` 恢复。
