# 本轮待人工获取的原文

更新：2026-10-29。用途是补齐**动态、非均匀软组织扫描中的接触控制研究现状**，不是重新检索idea或判断C1已经达到顶刊标准。下列条目尚不能支持未核读的方法、条件或性能结论；拿到全文也须先核准相关性和版本，再进入相应现状位置。

本轮已完成7张原文选段证据卡。ARIS现为`HUMAN_SEARCH_REQUIRED / human_fulltext_batch`，正式请求为C01–C08共8篇；S01–S03是项目此前已知、仍影响最近近邻判断的3项补充缺口。一并列出，方便一次获取，不把后3篇伪装成Controller已选中的候选ID。优先级表示核查顺序，不是期刊质量或创新性评分。

可将PDF直接放入本项目的 `literature/manual_inbox/`（建议按下表命名）。该目录不自动提交Git；我会核对题名、DOI、版本、正文可读性及哈希，再复制到ARIS的`source-materials/`并生成一次批量登记。你不必手填JSON，也不必修改运行目录。PDF或完整作者稿均可；只有摘要、网页截图或补充材料不等于正文全文。若确实有拿不到的，集中注明编号和原因，不用逐篇往返。

## ARIS正式请求：8篇

| 编号／建议文件名 | 题名、来源与DOI | 优先级及需要原文回答的问题 | 将来的科学位置 |
|---|---|---|---|
| C01.pdf | **Machine Learning in Robotic Ultrasound Imaging: Challenges and Perspectives**；Annual Review of Control, Robotics, and Autonomous Systems，2024；[10.1146/annurev-control-091523-100042](https://doi.org/10.1146/annurev-control-091523-100042) | P2：以任务专门综述交叉检查学习接触控制分支和引用链，不能用综述作者的分类替代原始方法核读。 | 调研地图与线索，不能独自证明具体控制能力。 |
| C02.pdf | **Sensorless environment stiffness and interaction force estimation for impedance control tuning in robotized interaction tasks**；Autonomous Robots，2021；[10.1007/s10514-021-09970-z](https://doi.org/10.1007/s10514-021-09970-z) | P1：核清双EKF的状态、力/刚度关系、环境运动假设及估计如何改变阻抗。 | 在线接触辨识与控制更新的直接机制近邻。 |
| C03.pdf | **A Lightweight Robotic Ultrasound System With Hybrid Force Control and Modular Force-Range Switching**；T-RO，元数据2026；[10.1109/TRO.2026.3734315](https://doi.org/10.1109/TRO.2026.3734315) | P1：核清力范围切换、硬件/反馈作用、呼吸与变形条件；同时核准最早公开时间及版本。 | 超声接触控制任务近邻；硬件能力不改写成纯软件控制收益。 |
| C04.pdf | **Adaptive Robust Interaction Force Control of a Robotic Manipulator in Uncertain Environments**；TIE，2025；[10.1109/TIE.2024.3525103](https://doi.org/10.1109/TIE.2024.3525103) | P1：核清估计/前馈、鲁棒反馈与内层运动控制的分工，以及动态环境和未知刚度的实际假设。 | 不确定接触中的力控制机制，避免只看与C1相同的估计结构。 |
| C05.pdf | **Energy-Aware Adaptive Force-Admittance Control for Robotic Manipulators: Application in High-Precision Aircraft Assembly**；TIE，元数据2026；[10.1109/TIE.2025.3616403](https://doi.org/10.1109/TIE.2025.3616403) | P2：核清自适应导纳与能量罐如何协调、冻结/限制哪些更新、力跟踪成立条件及计算代价。 | 本域相关力控制与稳定执行机制；装配结果不当作组织扫描实证。 |
| C06.pdf | **Safe Learning by Constraint-Aware Policy Optimization for Robotic Ultrasound Imaging**；TASE，2025；[10.1109/TASE.2024.3378915](https://doi.org/10.1109/TASE.2024.3378915) | P1：核清策略的真实动作、力幅约束、训练与部署成本、脊柱模体实验。摘要不能证明它在线优化阻抗。 | 超声接触策略的直接任务先例；按动作和目标归类，不能因没有阻抗名称便排除。 |
| C07.pdf | **Concurrent nonparametric estimation of organ geometry and tissue stiffness using continuous adaptive palpation**；ICRA，2016；[10.1109/ICRA.2016.7487609](https://doi.org/10.1109/ICRA.2016.7487609) | P1：核清几何/刚度表示、触诊激励、连续更新和估计到动作的关系；区别静态触诊与移动组织恒力扫描。 | 联合接触估计先例及CARE引用链。 |
| C08.pdf | **Safety-Critical Passive Control for Robotic Soft Tissue Interaction via Control Barrier Functions**；TMECH，元数据2026；[10.1109/TMECH.2026.3725669](https://doi.org/10.1109/TMECH.2026.3725669) | P2：目前仅核准身份，未取得实际摘要。需原文判断CBF/无源性针对什么闭环、何种约束及组织实验，并核准公开时间。 | 待核的稳定执行基础；题名不能直接支撑安全保证或认定覆盖C1。 |

## 一并获取的既有关键缺口：3篇

| 编号／建议文件名 | 题名与DOI | 优先级、用途及现有边界 |
|---|---|---|
| S01.pdf | **Observer-based sliding mode admittance control for ultrasound robot force-tracking in complex interaction environments**（DOSMAC）；CEP，2026；[10.1016/j.conengprac.2026.106764](https://doi.org/10.1016/j.conengprac.2026.106764) | P1：本项目最近的动态接触替代路线之一。现有出版社摘要/贡献段支持把刚度和接触平衡位置变化作为总扰动补偿；必须取得完整控制式、假设和实验，不能据缺全文声称它不能处理共同变化。 |
| S02.pdf | **Observer-Based Variable Impedance Control Using Moving Horizon Estimation for Robot Machining Thin-Walled Workpieces**；TIE，2024；[10.1109/TIE.2023.3292856](https://doi.org/10.1109/TIE.2023.3292856) | P1：核清移动时域估计如何同时处理环境阻抗和机器人扰动，并驱动VIC；官方摘要/引言已有该关系，仍缺完整方程与条件。 |
| S03.pdf | **Autonomous robotic ultrasound scanning on complex surfaces: Automatic trajectory generation and disturbance observer-based model predictive admittance control**；[10.1016/j.conengprac.2026.107240](https://doi.org/10.1016/j.conengprac.2026.107240) | P1：直接超声DOB＋预测导纳线索。既有出版社索引为**2027年1月期次**，最早在线公开日期未核准。先核准日期/版本，再判断是否落在已批准的2026-10-26文献截止日内；不能仅凭DOI含2026便当作2026已公开能力。 |

正式请求的题名、作者、ID及DOI来自[Controller批次](aris_human_fulltext_request.json)。[登记模板](aris_fulltext_manifest_template.json)只包含C01–C08，路径是未来复制到ARIS后的预定位置，目前不代表文件已存在、已登记或已阅读；S01–S03取得后需另按正常身份、筛选和阅读流程进入。

## 接续方式与尚未完成的部分

收到原文后，继续当前文献阶段：核准身份/时间与来源，阅读方法和实验，完成证据卡，修订同一研究现状及ARIS地图，补必要的前后向引用链，并检查是否还出现新的直接近邻。不会直接重启idea、方法或仿真。

本轮受[ARIS research-lit/SKILL.md](<D:/桌面/科研Agent Harness设计/ARIS/skills/skills-codex/research-lit/SKILL.md>)中的明确规则约束："Automatic full-text acquisition is arXiv-only." 该文件同时要求将其他已选原文汇总为一次`defer-fulltext-batch`，由`finish-reading`产生人工下载批次。此次清单按该流程生成，并符合用户允许批量人工获取的指示；它不是再次请求检索授权。
