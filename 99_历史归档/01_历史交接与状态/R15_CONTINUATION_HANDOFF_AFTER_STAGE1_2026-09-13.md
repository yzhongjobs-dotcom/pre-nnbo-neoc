# R15 后续算法续写 Handoff

## 1. 接续任务

在现有正文中继续完成 R15 的后续方法，使最终文件形成一篇从**总研究问题**出发、依次贯通以下主线的完整研究问题与技术路线：

$$
\text{柱面 Hertz 接触模型}
\rightarrow
\text{空间接触观测与 MHE}
\rightarrow
\text{未来环境参数和接触模型}
\rightarrow
\text{NNBO/single-critic HJB}
\rightarrow
\text{ANEOC 预测—校正}
\rightarrow
\text{实际最优阻抗}.
$$

第一阶段已经完成并获得用户认可。下一对话应从第一阶段输出如何进入最优控制开始，自然续写 NNBO、ANEOC 和实际阻抗生成，不重新改写已经完成的研究问题与接触建模主线。

## 2. 唯一续写基线

开始工作前必须完整阅读：

1. `R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_2026-09-12.md`：唯一正文基线；
2. `R15_RESEARCH_PROBLEM_SPATIAL_PREVIEW_CONTACT_MODEL_APPENDIX_2026-09-12.txt`：唯一补充推导基线。

正文目前完成至 **3.7 第一阶段形成的完整路线**，编号止于式（25）。后续正文应插入现有“参考文献”之前，原则上从第 4 节和式（26）继续。附录目前完成 A—G，后续证明从 H 开始续写。

以下文件只能用于查找后续思想、公式来源和历史问题，不能整段迁移，也不能覆盖当前基线：

- `R15_FINAL_METHOD_SPATIAL_PREVIEW_ANEOC_OPTIMAL_IMPEDANCE_2026-09-11.md`；
- `R15_METHOD_POST_STAGE1_PREVIEW_CONDITIONED_OPTIMAL_IMPEDANCE_2026-09-11.md`；
- `如何理解新method.md`；
- 其他旧版 audit、stage 或 handoff 文件。

旧文件中出现的有限元接触、接触 QP、ROM、全局区间外包、固定增益必要性门和多条备选路线，不属于当前已确认主线，不得因旧文件存在而重新引入。

## 3. 已锁定的研究含义

R15 的总问题是：机器人连续扫描空间变化的软组织时，怎样从接触过程获得新的环境信息，并在环境变化后快速得到相应的最优控制解和实际阻抗，而不是每次从头求解完整 HJB。

第一阶段已经闭合以下接口：

$$
({}^BT_{P,i},\mathbf F_{S,i},\mathbf M_{S,i})
\longrightarrow
(d_i,s_i,v_i,y_i)
\xrightarrow{\rm MHE}
(\hat\theta_f,\Delta\hat\theta,\hat b,\hat z,\Sigma_\zeta)
\longrightarrow
\widehat{\mathcal M}_{f}.
$$

后续必须保持以下对象含义：

- $\theta=(E^*,\mu,k_d)$ 是有效接触参数；$\theta_0,\theta_f$ 分别表示当前和前方区域。$\theta_0$ 由过渡前标定或上一窗口提供，并作为带协方差的到达先验进入 MHE；新窗口主要估计 $\theta_f$ 和边界 $b$；
- $\mathcal M(\theta)$ 或 $H_H$ 是接触模型；接触力/力矩是模型输出，不是 critic 权重或阻抗；
- $\eta$ 是后续最优控制问题的参数描述符，应由第一阶段输出和确有必要的任务量构成；不能不加说明地把 $\theta$、$\eta$ 和系统状态混为一体；
- $W$ 只表示 critic 权重；接触力旋量继续写作 $\mathcal W_c$，不得复用 $W$；
- $Z^*(\eta)$ 表示所定义有限维 HJB–KKT 问题的完整局部根，$W^*(\eta)=\Pi_W Z^*(\eta)$；
- critic 权重不是阻抗。它先确定价值函数和最优动作，再通过明确的控制/阻抗生成关系得到实际施加阻抗 $g_a$；
- 环境参数已估计并不意味着不再需要 HJB：$\theta$ 决定被控对象的接触动力学，HJB 仍用于在任务代价、机器人状态和约束下求最优反馈策略。

空间预瞄的本质是有限接触足迹在扫描过程中提前包含前方材料信息。方法不预设“必须在前缘到中心之间完成全部调整”这一人为截止标准；实际要求是环境发生接触变化后足够快地完成估计和阻抗调整。若中心到达前完成，只能作为更强性能结果，而不是方法成立条件。

## 4. 后续正文应完成的唯一主线

后续建议按以下三个连续阶段写入同一正文；小节可根据推导自然拆分，但不得增加平行方法路线。

### 第二阶段：由未来接触模型定义并求解最优阻抗

从 $\widehat{\mathcal M}_{f}$ 出发，明确最优控制的状态、动作、动力学、代价和约束，再写出 Hamiltonian、HJB、single critic 近似、策略恢复关系以及有限维根 $F(Z,\eta)=0$。必须把“接触模型怎样进入动力学”“为什么已知环境后仍要解 HJB”“critic 怎样决定控制动作和阻抗”逐步写清，不能用“使用 NNBO 求解”代替公式。

对原 NNBO 只能迁移其经核实且适合本问题的机制。原文中的 Hunt–Crossley 被控对象、球形接触关系或 force-as-action 语义不能自动照搬；本研究必须使用第一阶段已经建立的柱面 Hertz 接触模型，并最终明确一个一致的实际阻抗生成接口。

### 第三阶段：用 ANEOC 快速迁移同一个最优根

在同一 HJB 表达、critic 基、状态维度和约束布局下，由

$$
F(Z^*(\eta),\eta)=0
$$

推导根灵敏度线性方程

$$
F_Z C_\eta=-F_\eta,
\qquad
C_\eta=\frac{dZ^*}{d\eta},
$$

预测

$$
Z_{\rm pred}=Z^*(\eta_k)+C_\eta\Delta\eta,
$$

再用新描述符下的 HJB–KKT 残差进行 Newton/corrector。正文应说明：预测器消除局部参数变化的一阶影响，正则光滑分支内的理想截断误差为 $O(\lVert\Delta\eta\rVert^2)$；corrector 把预测点拉回新根附近。

必须简洁回答“为什么不重新求 HJB”：ANEOC 复用旧根及其局部灵敏度，使新问题从接近新根的位置开始；但“更快”只能由与**同一根问题**的 cold solve 和普通 warm-start corrector 的计算实验支持，不能仅由公式宣称。

### 第四阶段：从最优根到实际阻抗及完整在线算法

明确从 $Z^*$ 提取最优动作、形成阻抗命令并作用于机器人动力学的唯一链路。正文只需保留必要状态、约束和执行关系；完整在线过程最终用一段清晰伪代码汇总，不再逐句重复前文公式。

直接学习 $\eta\mapsto W^*(\eta)$ 只属于尚未选定的后续扩展，不得与 ANEOC 主线并列写入当前 Method，除非用户另行要求。

## 5. 每个模块的固定写法

每个模块按同一顺序写：

1. **先用一句直白的加粗句说明该模块做什么。** 例如“使用 MHE 从接触历史估计前方环境参数”“通过 HJB 定义给定环境下的最优阻抗”“使用 ANEOC 从旧最优根预测并校正新根”。
2. 紧接着定义本模块的研究子问题、输入和输出。
3. 按实际执行顺序给出完整核心公式流，并在公式后立即解释每个对象的物理或算法含义。核心关系必须显式写出，不能以“沿用某文献”“通过某方法求解”或抽象映射替代。
4. 在模块路线写完以后，最后用一个短段落写“**选择理由与适用条件。**”说明为什么采用该方法、它成立所需的关键条件及主要失效边界。
5. 对读者自然会提出的问题，在相关位置用一两句话回答；不要集中写成大段防守性说明，也不要让条件和例外淹没主方法。

正文应达到专业人员可沿公式一行一行读懂并复现逻辑的程度。核心 method、核心公式和关键接口放正文；较长推导、误差界、正则性证明、局部收敛证明及补充数值核对放附录。正文不能只给结论并把理解所需的公式全部推到附录。

## 6. 准确性、证据与写作边界

- 所有推导和符号必须逐项核对；没有理论依据的关系不得写成等式或定理。
- 对已有理论或继承机制，在必要位置引用已核实的原始论文；引用适量即可，不要句句引用，也不能完全没有支撑。
- 理论上成立、模型内部数值核对通过、真实组织实验有效和实际实时性达标必须分开表述。
- 当前跨区域柱面关系是正文已经声明的低维空间参数化近似。后续把它作为已定义的有效接触模型使用，并保留其验证边界；不得因此替换为更复杂的新接触模型。
- 不写“汇报版”“讲解版”等措辞。语言应像论文 Method，但保持直观、简洁，不堆叠抽象概念。
- 不为了极端情况增加支线、备用算法或大段安全机制。只有会直接破坏所选方法正确性或可执行性的条件，才以最小充分形式纳入。
- 不把尚未实验验证的速度、精度、稳定性或优越性写成既成事实。

## 7. 编辑纪律

严格遵守：

- **Surgical editing：**仅修改完成当前阶段所必需的内容，保留现有正确结构、文字、公式和编号关系；
- **Minimal sufficient change：**不增加推测性抽象、兼容层、重构、清理或相邻改进；
- **Target-state editing：**直接写最终正确状态，不保留旧方案、修改历史或“过去如何错误”的叙述；
- 如果局部补丁会留下公式、符号或结构冲突，只做解除该冲突所必需的重构；
- 除正文和对应附录外，不修改其他文件。

第一阶段原则上冻结。只有发现能够准确证明的数学、符号或接口错误时，才允许最小修正，并必须说明原因；不得借续写后续算法之机重写第一阶段。

## 8. 后续优先核实的原始材料

1. NNBO：
   `D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-mpvic/Neural-Network-Based_Optimal_Impedance_Control_for_Robots_in_Physical_Interaction_With_Soft_Environments.pdf`
2. ANEOC/NEOC：
   `D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-topvenue-calibration/Neighboring Extremal Optimal Control Theory for Parameter-Dependent Closed-loop Laws.pdf`
3. 隐式微分核对：
   `D:/桌面/科研Agent Harness设计/impedance-control-e2e/source-materials/user-batch-mpvic/Revisiting Implicit Differentiation for Learning Problems in Optimal Control.pdf`
4. 已提取文本可辅助定位：`tmp/pdfs/nnbo.txt`、`tmp/pdfs/neoc.txt`、`tmp/pdfs/idoc.txt`；正式引用和公式仍以原始论文为准。

使用旧 R15 文件时必须回到这些原始材料核实，不能把旧稿中的公式视为可靠来源。

## 9. 下一对话的完成标准

最终交付仍是当前同一份完整正文和同一份附录。完成后，读者应能不依赖旧文件回答：第一阶段得到了什么；未来接触模型怎样进入最优控制；single critic 怎样表示价值函数并产生阻抗动作；ANEOC 怎样利用参数变化更新同一个最优根；corrector 怎样恢复新根；最终阻抗怎样作用到机器人。

如某一环节存在必须推倒后续路线才能解决的硬伤，应停止写作并准确、简洁说明；除此之外，应在实验前把理论和接口做扎实后继续完成，不因一般不确定性或可实验验证项中断。
