R15 顶会顶刊研究问题与 Method 收敛交接包
版本：2026-09-08

一、用途

本包用于在新对话中无损延续 R15。它把研究现状、问题判断、候选 Idea/Method、历史踩坑、当前进度和证据边界拆成独立文件，防止新对话把历史候选、审计意见和最终结论混在一起。

二、阅读顺序

1. 01_CURRENT_REQUIREMENTS_AND_PROGRESS.txt
   用户目标、已完成工作、当前冻结结论、尚未完成事项。

2. 02_FIELD_STATE_ANALYSIS_AND_FINAL_RESEARCH_PROBLEM.txt
   研究现状、strongest-prior 成对比较、真实领域痛点、最终研究问题与 venue-level 判断。

3. 03_IDEA_AND_METHOD_BLUEPRINT.txt
   从最终研究问题自然推导出的核心 Idea 和候选 Method 结构。该文件是下一阶段推导起点，不是已经冻结的 Final Method。

4. 04_PITFALLS_AND_NONREGRESSION_RULES.txt
   历史路线为何失效、已发现的数学/工程错误、禁止回退项和审稿人攻击点。

5. 05_EVIDENCE_AND_COMPARATOR_INDEX.txt
   核心论文、各自真实作用、venue calibration 与 strongest prior 的使用边界。

6. 06_E2E_HARNESS_REUSE_GUIDE.txt
   E2E 的可直接复用资产、需重新验证内容、禁止继承项、本地独有论文目录，以及 Harness“默认不用、仅可选参考”的明确边界。

三、当前一句话结论

连续机器人超声中的患者特异、空间变化软接触确实构成重要研究问题；应继续该方向，但问题必须从“在线辨识组织参数后选择阻抗”收敛为：在不进行独立逐点触诊的情况下，利用连续扫描中的安全交互同步更新与阻抗决策相关的接触信念，并生成满足任务、接触、执行器和能量约束的可执行法向—切向阻抗。

四、当前科学状态

- 已完成：领域与 strongest-prior 调研、顶会顶刊经验标尺核对、研究问题去留判断、候选科学增量收敛。
- 未完成：Final Method 的公式级完整推导、结构可辨识性证明、必要性命题、统一目标、实时近似、安全发布定理与完整工程接口。
- 因此，新对话不得声称“Final Method 已完成”或“已达到顶刊标准”；下一阶段目标正是完成并审计这些内容。

五、新对话建议首条指令

请完整读取本交接包 00—06。不要重复泛化领域综述，也不要继承任何历史 Method 为正确结论。按 06 检查本地 E2E 与 GitHub scientific-core，优先复用原论文、Evidence、strongest-prior 和失败历史，但不得继承旧 accepted/final 状态；Harness 默认不用，若 Academic Research Suite 已有更成熟能力则完全忽略，只可按需参考少量普通文件。以 02 中冻结的研究问题和 strongest-prior 边界为起点，严格按《顶会顶刊_Method_经验标尺_Reference_Distribution_v1_1_Audited》执行：先证明 Method necessity 与最小可辨识参数化，再从第一性原理推导单一、闭合、可执行的 Final Method；每个新增结构必须对应一个已确认 failure，禁止方法拼接。实验设计暂后置，但必须保留可验证的 claim 和 killer condition。完成 Candidate Method Seed、Validity/Novelty/Method Worthiness、公式、理论、求解与执行闭合并重新审计后再停止。

六、优先复用文件

- 顶会顶刊_Method_经验标尺_Reference_Distribution_v1_1_Audited
- R15_Strongest_Prior_Research_Package_2026-09-07.txt
- R15_STRONGEST_PRIOR_FINAL_AUDIT_2026-09-07.txt
- R15_v4_1_Field_Foundation_TopTier_Novelty_Reaudit_2026-09-06.txt
- GitHub：zyanzhi90-dot/e2e，scientific-core 分支
- GitHub：zyanzhi90-dot/harness，scientific-core 分支
- 本地：D:\桌面\科研Agent Harness设计\ARIS
- 本地：D:\桌面\科研Agent Harness设计\impedance-control-e2e
- 本包 01—06

若这些文件之间出现矛盾，优先级为：原论文全文与公式证据 > 审计者重建的 Scientific Delta > 作者自述 > 历史路线结论。
