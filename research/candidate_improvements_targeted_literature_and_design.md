# 两候选的定向文献、设计与快验
<!-- project-update:start -->

**历史适用范围（2026-10-20审查）：** 本文保存早期定向分析与MHE预检的当时安排，第4节要求至少落实NNBO、暂停训练不是当前默认指令。用户后来已明确：先由研究问题/方法决定各候选代码和仿真，再选适合的公平比较；不强行使用不适合的基线。后续进度见[HANDOFF](../project_records/HANDOFF.md)及[双候选任务仿真](../experiments/dual_candidate_task_simulation/README.md)。

<!-- project-update:end -->

2026-10-19。用户已授权在原问题内检索、设计、快速验证并持续迭代。本文承接[候选定义](optimal_impedance_question_and_methods.md)，不修改长期要求或另设顶刊录用门槛。结论随证据修订。

## 1. 围绕具体设计检索和核读

本轮检索包括“MHE variable impedance control”“moving parameters optimization adaptive impedance”“Sobolev approximate nonlinear MPC”“impedance guidance learning”“model predictive impedance gradient neural”等组合，并查具体标题及作者原始来源。检索日期为2026-10-19；搜索排序/未检出不能证明没有相同工作。没有将原摘要筛查写成全文审计。

| 论文与原始来源 | 本轮核读范围 | 成熟能力、实际调整与本项目启发 |
|---|---|---|
| Roveda等，*Q-Learning-based model predictive variable impedance control for physical human-robot collaboration*，Artificial Intelligence312(2022)103771，[作者PDF](https://ipg.idsia.ch/preprints/Roveda2022a.pdf)，[DOI](https://doi.org/10.1016/j.artint.2022.103771) | PDF引言、控制/网络与更新关键段，尤其§2.6–2.8、Algorithm3、实验设置；未完整审计全部稳定性证明 | 预测控制与actor–critic已直接结合，动作含参考位置及阻尼；周期更新模型/critic/actor。可学习其问题—优化—反馈近似的组织方式，不能把网络替代阻抗优化单独称创新。它的任务/输出/更新时序不同于本项目，不能由摘要推断它忽略某一接触困难。 |
| Lu等，*You Only Plan Once: A Learning-Based One-Stage Planner With Guidance Learning*，RA-L2024，[本地原文](../literature/papers/method_transfer/YOPO.pdf)，[作者PDF](https://tju-air-lab.github.io/projects/YOPO/YOPO.pdf)，[作者代码](https://github.com/TJU-Aerial-Robotics/YOPO) | 原文训练目标/数值梯度/特权地图与部署观测；本轮另核作者仓库说明和loss文件入口。不是完整复现 | 继承轨迹优化及运动原语，将代价梯度用于学习，部署生成候选/评分。迁移重点是用任务后果训练，而非只拟合某个最优增益标签；需要把梯度改为接触和实际增益执行的后果。没有照搬其多候选飞行几何，也不称PPO。 |
| Vahidi-Moghaddam等，*Extended Neighboring Extremal Optimal Control with State and Preview Perturbations*，[arXiv2306.04830v1](https://arxiv.org/html/2306.04830v1) | §III-B/C、式11–36、算法1及约束奇异性说明；非全文证明审计，本文也不把预印本身份写成顶刊 | 成熟NE扩展到状态和预瞄扰动的两反馈矩阵，并处理非最优名义解/较大扰动。候选1可直接借鉴，不再声称“参数/预瞄敏感更新”本身独有；正定性、活跃约束变化和状态约束奇异性须检查。尚未将该算法实现为接触控制器。 |
| Zhang等，*Observer-Based Variable Impedance Control Using Moving Horizon Estimation for Robot Machining Thin-Walled Workpieces*，TIE71(6):5972–5982(2024)，[出版社](https://doi.org/10.1109/TIE.2023.3292856) | 出版身份、出版社可检索摘要；本轮未取得可核读全文 | MHE驱动变阻抗已有直接先例；不能将MHE＋变阻抗组合本身列为新贡献。尚不能断言其不具备本项目某项能力，不根据不可靠二手自动摘要推定细节。 |
| Wang等，*Moving Parameters Optimization-Based Model Predictive Adaptive Impedance Control*，[出版社](https://doi.org/10.1109/TIE.2025.3587119) | 身份/索引摘要；未核全文 | 有面向磨削的移动参数优化/预测自适应阻抗近邻；需继续取得方法细节，再界定与快速参数更新的区别。 |
| Lüken等，*Sobolev Training for Data-efficient Approximate Nonlinear MPC*，IFAC2023，[出版社](https://www.sciencedirect.com/science/article/pii/S2405896323009126) | 出版/摘要索引，直接页面获取受限；未核原文推导 | 用参数灵敏度提供导数监督是可查的外域线索。不能把未读推导或稳定性结论移用本项目；本轮未将其选为必须新增的训练模块。 |

## 2. 如何形成本轮方法，而不是只改论文故事

从已有两候选出发，第一步保留共同估计前端和实际阻抗执行。第二步用普通热启动预测优化建立可比较的模型路线，同时将候选2落实为真正的神经网络和任务代价引导训练。第三步让训练面对共同估计误差和实际增益执行；不能用理想材料标签、瞬时K/D跳变的好成绩代替闭环能力。联合参考/力补偿作为合法动作保留，并与纯K/D分开判断。

这里继承了MHE、接触模型、滚动优化、神经反馈和能量约束等成熟机制。具体适配关系是：**有限足迹带延迟接触信息 → 局部联合估计 → 执行状态与估计共同条件化决策 → 以真实过渡后的力/运动代价训练**。这是一项可检验的方法改造，还不是已建立的顶刊贡献；组件组合不是其价值证明，也不因组件成熟而否定其潜力。

实现、划分及边界见[实验协议](../experiments/mhe_guided_impedance_pilot/design_and_protocol.md)。本轮没有改写B0三核心或冻结快照；普通打靶不冒充HJB/NEOC，约化滚动窗估计也不冒充完整B0 MHE。

## 3. 当前迭代证据

两例训练分布工程预检已打通约化前端和执行链；预测当前Fn的RMSE约0.023–0.026N，单次MHE约3ms、模型打靶决策约31ms（本机、零发布等待；不是实时认证）。模型路线一例略改善、一例变差。对这两条轨迹抽取28个相关状态，模型与真值任务梯度有25%方向余弦为负，模型局部优化后35.7%的真值短时成本反而上升；只能支持本原型的局部模型决策诊断，不能推广为领域未解问题。

固定辅助调优使用训练12例，之后采集训练24例固定/解析各一条轨迹。网络训练已开始，但收到用户暂停纠正后终止，未生成完成actor；没有验证集或独立测试结果。上述结果不能评价完整候选2，更不能评价顶刊贡献。

## 4. 用户纠正：主对照必须落实到领域方法

用户明确指出当前对比基线层次错误：主要比较应至少包含NNBO原文及相关方法，不能以调优固定阻抗作为顶刊创新的主要依据。**这项纠正已接受并暂停当前训练。** 固定控制可以揭示工程实现是否有明显问题，或作为消融参照，但已有变阻抗/最优阻抗能力才是本项目贡献比较的主要对象。“变阻抗比固定好”不能代替“相对已有变阻抗具体改进了什么”。本轮缺失的是领域强方法的落实，不应先增加更多固定控制调参或继续网络训练。

主对照按具体能力选择，当前确定的落实次序如下；这是待实现工作，不是已完成比较：

| 主要对象 | 原始能力及核查入口 | 要落实的比较与避免的弱化 |
|---|---|---|
| **Kong等NNBO，TSMC2025**，[原文](../literature/papers/NNBO.pdf)、[已有最小机制复现](../experiments/nnbo_source_reproduction/README.md) | HC辨识＋single critic＋Newton Downhill位置参考；原式22权衡位置误差与作用力平方，Remark6的期望力为零 | 先补原辨识/机器人参考执行链核对，再接共同扫描任务。保留其位置参考输出，不强制改成四K/D。非零目标力/移动接触适配须明确推导、记录与原版区别，不能直接把目标失配造成的劣势算创新；已有最小复现未接位置内环，不冒充完整基线。 |
| **Jin等预测变阻抗，TMECH2023，E41**，[原文](../literature/papers/field_papers/E41.pdf)，DOI10.1109/TMECH.2022.3204350 | 将增益作为控制量；OS单步QP、MS多步预测及能量修正，直接处理精度—柔顺权衡 | 核读并落实其有代表性的可变增益优化接口/约束，比较闭环性能与全成本；不能将当前任意打靶脚本直接命名为E41复现。 |
| **Li等未知地形阻抗学习，T-RO2025，E27**，[原文](../literature/papers/field_papers/E27.pdf) | 随机地形仿真预学习＋变阻抗反馈、sim-to-real条件；本地已有全文 | 与候选2直接相关，须核清其网络输出、反馈、训练和迁移条件后落实；不能用旧短历史PPO代替这篇强方法，也不能只复用名称。 |

模型学习路线还以Zhao等TIE2022（E38，DOI10.1109/TIE.2021.3134082）及Roveda2022为设计/比较近邻，按任务匹配和实际可实现范围选择补充，不以方法数量替代最接近比较。

比较分两层核准：先在原论文设置检查关键算法与执行接口；再做清楚标记的共同扫描适配。共同任务、调参预算、训练/辨识准备、可得传感器、实际动作自由度、执行约束及计算/发布成本都须交代。既检查各方法原生信息配置，也用必要的信息/动作对齐区分贡献来自何处；不把原方法剥去其优势后称公平。没有完成这些实现和验证前，本轮原型只处于工程预检，不能写“已证明优于NNBO/强变阻抗”或“已满足顶刊标准”。
