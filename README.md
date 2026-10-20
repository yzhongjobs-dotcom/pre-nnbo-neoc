# 环境变化时的快速阻抗更新
<!-- project-update:start -->

**2026-10-20：关键内容已恢复中文原文。** 长期要求、交接、决策、失败复盘、科研流程参考及索引以英文化前的Git中文版本为基础，目录和文件名继续使用英文。Claude Code通过[CLAUDE.md](CLAUDE.md)及[AGENTS.md](AGENTS.md)直接阅读同一套中文要求，不另维护英文摘要约束。历史路径见[迁移记录](project_records/english_migration.json)。本次只恢复文件并审查一致性，不启动科研或实验。

可运行 `python -X utf8 scripts/verify_english_migration.py` 核对中文原文、路径、代码语法、链接和证据保全；技术检查不等于证明所有科研结论或翻译语义正确。

<!-- project-update:end -->

研究目标：机器人连续扫描软组织时，利用新接触信息，及时生成有实际力与运动控制收益的近最优阻抗。目标论文须有实质科学贡献和足以参与顶刊竞争的证据；当前尚未形成这样的完整方法与证据。

先读[接手指南](project_records/HANDOFF.md)。用户确认的方法形成与贡献尺度已落实[长期要求](project_records/long_term_requirements.md)，具体依据见[论文逐例分析](research/independent_reverse_analysis_of_top_journal_methods.md)。最新已按原问题分别实现[两候选任务仿真](experiments/dual_candidate_task_simulation/README.md)，完成[24场景、消融、泛化及发布延迟实验](experiments/dual_candidate_task_simulation/results.md)，随后与Jin多步预测变阻抗的明确任务适配比较。目前力/切向误差仍落后，两候选均未建立顶刊级贡献；候选1采用局部二次价值反馈，并非B0固定配置点根。

| 入口 | 用途 |
|---|---|
| [长期要求](project_records/long_term_requirements.md) | 已确认的目标、流程、实验顺序与具体汇报要求 |
| [研究问题](research/research_question.md) | 核心任务、B0起点及尚未验证的前提 |
| [研究现状](research/field_state.md) | 相关工作的已有能力、条件与全文访问边界 |
| [两套候选方法](research/optimal_impedance_question_and_methods.md) | 正确的候选定义、具体实验与完整构想的区别 |
| [可借鉴机制](research/transferable_mechanisms.md) | 已查阅的成熟机制及适配条件，不是已选新方案 |
| [文献索引](literature/literature_index.md) | 两个外部仓库及项目原有文献的本地原文、文本、来源与阅读范围 |
| [实验索引](experiments/experiment_index.md) | 后续可用的模型、求解器、对照、结果与复现入口 |
| [科研失败与纠偏](project_records/research_failures_and_corrections.md) | 本项目实际做偏的原因、后果及纠正办法 |
| [独立顶刊逆向分析](research/independent_reverse_analysis_of_top_journal_methods.md) | 样本论文的问题、方法和证据组织；与失败复盘互补 |
| [研究决策与变更](project_records/decisions_and_changes.md) | 必要历史结论、版本与本轮清理依据 |

B0的[方法](research/method_design.md)、[附录](research/technical_appendix.txt)、[条件性证明](research/theoretical_proofs.txt)是可修订的研究起点。冻结快照 `00_基线快照_2026-10-14/` 保持原样；完整旧材料由Git保存，不另建历史副本目录。

仓库：[yzhongjobs-dotcom/pre-nnbo-neoc](https://github.com/yzhongjobs-dotcom/pre-nnbo-neoc)。完成授权改动后，按[AGENTS](AGENTS.md)检查、记录并运行 `python -X utf8 scripts/research_sync.py --message "research: 具体改动和目的"`。
