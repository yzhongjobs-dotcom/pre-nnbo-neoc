# ARIS 正式文献补充：等待 source policy 审批

记录日期：2026-10-31。**正式 run 已启动；政策经 Controller 校验 PASS，按用户要求停在人工审批关口，尚未检索。本文是关口报告，研究现状补充尚未完成。** 项目仓库未修改，未执行 pull、提交或推送。

## A. Controller 最终状态及所需审批

| 项目 | 实际状态 |
|---|---|
| run_id | `recent-related-work-20261031` |
| RUN_ROOT | `D:/pip_env/temp/scan_impedance_literature_20261031/aris_run` |
| Controller-managed | `true` |
| 当前阶段 | `WAITING_FOR_HUMAN` |
| 停止原因 | `waiting_for = source_policy_approval` |
| 政策校验 | `PASS`，2026-10-31T16:03:06Z |
| 实时审批请求 ID | `ba52f967af6b4c879576416eba9e8d95` |
| 政策 SHA-256 | `b467b81f955d1690249d3e803f7c4fca75db9115f80f704d254a53e91cdc8665` |
| allowed-actions | `human_approve`、`request_source_policy_revision` |
| allowed-agents | 空；此阶段未派出 reader 或 coverage reviewer |

采用新 run：旧 `c1-related-work-20261026` 当前仍在 `QUERY_PLANNING`，累计 32 次查询，并绑定了尚未解决的 2018 tele-echography 原文/谱系缺口。本轮是新的 2022–2026 研究现状补充范围；独立绑定 source policy 和统计，可保留旧 run 未完成的工作，不混用旧查询与临时调查结果。

待审政策：[SOURCE_ADMISSION_POLICY.yaml](D:/pip_env/temp/scan_impedance_literature_20261031/aris_run/idea-stage/SOURCE_ADMISSION_POLICY.yaml)。该规范路径由 `submit-source-policy` 生成，没有直接写入。

**请审批的具体内容：**

- 用途是补研究现状、最近邻与领域进展；不能因方法采用图像、视觉、学习、专用末端或力矩接口，或难以做仿真对照而排除。
- 2022–2026 为新增范围；区分在线年、卷期年和预印本；历史锚点只作背景与引文扩展。明确覆盖 TUM/Navab/Jiang、非静止接触参数估计、力源运动估计以及用户指定的六条前向引文链。
- 阅读门槛采用用户列举的 venues，并提请将 **TIM、IEEE Sensors Journal、AIM** 一并列入，以涵盖已知直接超声/呼吸路线。其他来源按年龄校准引用门槛，或 Controller 记录的狭义最近邻/反证/诊断例外判定，不能仅凭相关或课题组声望放行。引用门槛严格大于：2000年前300；2000–2009年200；2010–2017年100；2018–2020年60；2021–2022年35；2023年20；2024年10；2025年3；2026年起1。
- 检索由 Controller 记录，优先 SerpApi Google Scholar；合法少结果/零结果继续细化与引文扩展。provider 故障与后备路由完全遵守 SKILL，不临时拼接其他发现源。
- 候选逐条真实摘要筛查；Field Map 与实时 Controller 请求下的独立原生 coverage reviewer 决定覆盖。全文自动获取仅限身份声明有 arXiv 的已准入论文，其余一次性人工批次，遇人工关口停止。

执行审批前可先检查状态。以下命令应在 **ARIS 目录**运行；审批命令由用户执行，或用户明确批准上述具体版本后再由代理记录：

```powershell
Set-Location -LiteralPath 'D:/桌面/科研Agent Harness设计/ARIS'
python -X utf8 -m arisctl --root 'D:/pip_env/temp/scan_impedance_literature_20261031/aris_run' status recent-related-work-20261031
python -X utf8 -m arisctl --root 'D:/pip_env/temp/scan_impedance_literature_20261031/aris_run' allowed-actions recent-related-work-20261031
python -X utf8 -m arisctl --root 'D:/pip_env/temp/scan_impedance_literature_20261031/aris_run' human-approve recent-related-work-20261031 source_policy_approval --decision approve
```

停止依据是本次用户要求及 [research-lit/SKILL.md](D:/桌面/科研Agent%20Harness设计/ARIS/skills/skills-codex/research-lit/SKILL.md) 的明确要求：`retrieval remains forbidden until the user approves that exact validated candidate.` 这不是权限故障或自动审批拒绝；代理未执行 `human-approve`，未伪造批准。

## B. 实际检索及覆盖

| 指标 | 本 run 实际值 |
|---|---:|
| provider 实际调用 | 无 |
| Controller query_count / query_events | 0 / 0 |
| 检索返回记录 | 0 |
| 去重候选 | 0 |
| admitted 新文献 | 0 |
| 全文获取 / 阅读 / Evidence | 0 / 0 / 0 |
| coverage reviewer | 未到阶段；无实时请求、无结论 |

仅检查到 `SERPAPI_KEY` 已配置，未联网探测 provider。旧报告的 9 条 R 记录与 35 条 L 记录已登记为 **44 条已知记录（含跨表重复），不是 44 篇去重论文或正式候选**；不计本轮新发现。指定 literature_index 条目中找到 16 行；当前索引没有 N01–N04、M01–M03 的独立行，其已知身份线索限于获准的草稿第一节和两份调查，后续须经 Controller 核实。

查询轴和六个引文锚点已写入政策。**尚未提交或执行 schema-v2 Query Plan**，因为当前阶段只允许政策审批/修订。正式 Field Map、综述与原始论文阅读、逐条筛查、六链前向扩展、定向细化、分页及年份/venue 覆盖、独立覆盖审查均未完成；不得宣称饱和。

## C. admitted 的新文献表

当前没有条目。审批前不得把上一轮候选复制为本轮 admitted，也不得填入未经本轮证据核实的机制或接口信息。后续表按以下列逐条补齐：

| 编号 | 完整题名 | 作者 | 年份 | venue | DOI | 主题类 | 刚度/运动信息来源 | 接口 | 研究现状中的作用（最近邻/路线代表/背景） | 全文状态（arXiv已下/待人工） |
|---|---|---|---|---|---|---|---|---|---|---|

## D. 对三条路线与缺口表述的影响

**本轮尚无新正式证据，不能据此修改三条路线或判断缺口成立。** 政策把以下机制问题保留为后续核查目标，不预设答案：

1. 激励辨识路线：在运动中辨识的究竟是材料刚度、接触几何还是相对形变；运动如何与参数区分；参数随后进入哪个控制通道。
2. 力反馈补偿路线：把运动合并进扰动与显式估计运动分别实现何种能力；是否在线更新参数；所需动力学及执行接口是什么。不能预设此路线均采用固定刚度。
3. 运动观测补偿路线：区分外部视觉、US图像、独立接触压头和受控探头自身力觉；分别核查连续反馈、重定位/重规划、预测前馈，不把它们混成同一种能力，也不预设所有预测均依赖周期。

针对本文的待核问题是：同一受控探头的力被在线刚度转换为等效运动量测时，参数误差如何耦合进预测前馈与反馈，现有工作已分析到哪些条件。三路线仅作为组织视角；分类边界、最近邻与剩余缺口须由正式 Evidence、Field Map 和 reviewer 后续修订，不能先写“首次”或“没有相关工作”。

## E. 人工下载清单

**当前没有下载条目；现在只需审批政策，不需下载。** 尚未发生 admission，Controller 的 `human_fulltext_request` 为 `null`，不能从旧临时报告移植一个正式批次。

预留目录已创建：`D:/pip_env/temp/scan_impedance_literature_20261031/aris_run/source-materials/`。后续对非 arXiv 的未读 admitted IDs 一次性执行 `defer-fulltext-batch`；完成 arXiv 读取/Evidence 后调用一次 `finish-reading`，在 `HUMAN_SEARCH_REQUIRED` 停下。届时才根据实时 Controller 批次，逐条提供 **完整题名、作者、年份、venue、DOI、目标目录和文件命名/manifest 要求**；一个 paper_id 对应一个文件，用户提交一份 manifest。当前不指定虚构文件名或要求提前下载。
