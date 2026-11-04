以下为本次实际联网核查结果。未读取、修改项目文件，未提交；未下载大批 PDF。成功打开的全文/正式页面包括：Unser–Fageot–Ward 的 SIAM 页面、Aravkin 等广义 Kalman 平滑全文、Hartikainen–Särkkä 的 PMLR 页面、Lin–Nassar–Evans 全文、PhysioNet 数据页、Garg–Ross 正式机构记录等。Speckman 正式页发生重定向错误，但检索页显示了摘要；Forssell–Ljung 的 ScienceDirect PDF 被拒绝访问，只看到正式摘要页。

## 1. A 组表

| 编号 | 完整题名 | 作者 | 年份 | 期刊或会议 | DOI 或链接 | 打开了什么 | 与本问题的准确对应 | 成立条件 | 在线成本 |
|---|---|---|---:|---|---|---|---|---|---|
| A1 | *Splines Are Universal Solutions of Linear Inverse Problems with Generalized TV Regularization* | Michael Unser, Julien Fageot, John Paul Ward | 2017 | SIAM Review | [10.1137/16M1061199](https://doi.org/10.1137/16M1061199) | 正式摘要页，含定理摘要 | 若对 \(Ld\) 使用 \(L_1\)/TV，解具有有限个非零创新或样条结构；支持“稀疏创新”数学机制，但不是人体动力学证据 | 线性测量、凸一致性约束、适当的广义 TV 空间；结论是 representer/样条解，不等于真实过程必为 Lévy | 批量凸优化；不能直接据此声称实时 |
| A2 | *Sparse Representations*（书中章节） | Michael Unser, Pouya D. Tafti | 2014 | *An Introduction to Sparse Stochastic Processes* | [Cambridge DOI](https://doi.org/10.1017/CBO9781107415805.009) | 正式出版社章节页，摘要/章节说明 | 明确讨论 \(s=L^{-1}w\)、Lévy 过程、广义增量和离散算子；与 \(Ld=w\) 直接对应 | 需选定 \(L\)、边界条件及创新分布；采样后通常得到离散近似 | 离线分析；在线需递推状态空间或滑动窗口 |
| A3 | *Kernel Smoothing in Partial Linear Models* | Paul Speckman | 1988 | JRSS B | [10.1111/j.2517-6161.1988.tb01738.x](https://doi.org/10.1111/j.2517-6161.1988.tb01738.x) | 检索到正式题录和摘要；正式页面重定向失败 | 模型 \(y_i=\xi_i^\top\beta+f(t_i)+\epsilon_i\) 与 \(F_i=k y_i+d_i+e_i\) 形式高度对应；摘要明确报告 \(\beta\) 的渐近偏差和方差 | \(f\) 平滑、误差条件、核/样条平滑参数适当；不是稀疏跳变模型 | 局部残差化可在线近似；完整样条通常需批处理 |
| A4 | *Root-N-Consistent Semiparametric Regression* | Peter M. Robinson | 1988 | Econometrica | [10.2307/1912705](https://doi.org/10.2307/1912705) | 题录和摘要页，未打开全文 | 部分线性模型中，先消除非参数项再估计有限维 \(\beta\)；可作为估计 \(k\) 的理论基线 | 关键是 \(y\) 与 \(d\) 的残差化后仍有有效变化；若 \(d\) 与 \(y\) 共变则会偏差/不可辨识 | 在线可用递推均值/局部平滑，但理论保证需重新建立 |
| A5 | *Double/Debiased Machine Learning for Treatment and Structural Parameters* | Chernozhukov 等 | 2018 | *The Econometrics Journal* | [10.1111/ectj.12097](https://doi.org/10.1111/ectj.12097) | 预印本全文，含 Neyman 正交定义 | 可将 \(d\) 或其条件均值视为 nuisance，并构造对 nuisance 一阶误差不敏感的 \(k\) 估计；但本问题不是因果处理效应问题 | 需要独立/交叉拟合、正交 score、样本或时间依赖条件、残差化后的 \(y\) 有方差 | 原始 DML 通常批处理；在线版本需时间序列交叉拟合和延迟更新 |
| A6 | *Latent Force Models* | Mauricio Álvarez, David Luengo, Neil D. Lawrence | 2009 | AISTATS | [PMLR 正式页](https://proceedings.mlr.press/v5/alvarez09a.html) | 正式摘要页；另打开 PDF 摘要 | 物理微分方程由 GP latent force 驱动；可把人体运动/组织动力学作为 \(d\) 的动力学先验 | 主要是线性 ODE、GP 驱动、核函数和参数可识别；默认是高斯、连续型驱动，不是稀疏冲击 | GP 批推断通常 \(O(N^3)\)；状态空间核可降为递推 Kalman 规模 |
| A7 | *Sequential Inference for Latent Force Models* | Jouni Hartikainen, Simo Särkkä | 2011 | UAI | [PMLR 正式页](https://proceedings.mlr.press/r9/hartikainen11a.html) | 正式摘要和书目信息 | 将 GP-LFM 改写为线性状态空间模型，用 Kalman filter/smoother；还明确处理 switching LFM | GP 核需有状态空间表示；switching 需要 SLDS/离散状态推断 | 线性模型可递推；switching 需 IMM/粒子法或近似，成本随模式数增加 |
| A8 | *Generalized Kalman Smoothing: Modeling and Algorithms* | Aleksandr Aravkin, James Burke, Lennart Ljung, Aurélie Lozano, Gianluigi Pillonetto | 2017 | Automatica | [10.1016/j.automatica.2017.08.011](https://doi.org/10.1016/j.automatica.2017.08.011) | 全文 PDF | 覆盖非高斯、稀疏、鲁棒、约束 Kalman smoothing；与 \(L d=w\) 加 Huber/PLQ 惩罚最接近 | 线性状态模型、PLQ 密度、适当可行性和块三对角结构 | 全序列 IP 平滑复杂度为 \(O[N(n^3+m^3+\ell)]\)；过滤/固定滞后可在线近似 |
| A9 | *Smoothed State Estimates Under Abrupt Changes Using Sum-of-Norms Regularization* | Henrik Ohlsson, Fredrik Gustafsson, Lennart Ljung, Stephen Boyd | 2012 | Automatica | [10.1016/j.automatica.2011.08.063](https://doi.org/10.1016/j.automatica.2011.08.063) | 正式作者页面和摘要 | 对状态转移残差使用 sum-of-norms，直接针对突变/脉冲扰动；比“高斯+稀疏”更接近候选机制的优化版本 | 线性状态空间、凸正则化；突变被建模为稀疏状态跳变，不保证连续身体运动符合该模型 | 论文称算法规模随时间长度近似线性，但通常仍是平滑而非严格单步在线 |
| A10 | *Robust and Trend-Following Student’s t Kalman Smoothers* | Aleksandr Aravkin, James Burke, Gianluigi Pillonetto | 2014 | SIAM Journal on Control and Optimization | [10.1137/130918861](https://doi.org/10.1137/130918861) | 正式 DOI/题录页和摘要信息 | Student-t 重尾适合异常/大扰动；它表达“重尾”，不等于“连续高斯+少数稀疏脉冲” | Student-t 尺度混合、状态模型和数值迭代收敛条件 | 迭代加权 Kalman；每次线性规模，需多次迭代 |
| A11 | *l1 Trend Filtering* | Seung-Jean Kim, Kwangmoo Koh, Stephen Boyd, Dimitry Gorinevsky | 2009 | SIAM Review | [作者正式页面](https://www.web.stanford.edu/~boyd/papers/l1_trend_filter.html) | 正式作者页面检索到摘要；页面打开超时 | 对 \(d\) 的高阶差分施加 \(L_1\)，产生分段多项式；三次多项式投影的稀疏自适应替代方向 | \(d\) 的高阶差分应稀疏；不直接建模高斯连续创新 | 论文明确给出专用内点法近似线性复杂度；滑动窗口可在线 |
| A12 | *Adaptive Piecewise Polynomial Estimation via Trend Filtering* | Ryan J. Tibshirani | 2014 | Annals of Statistics | [10.1214/13-AOS1189](https://doi.org/10.1214/13-AOS1189) | arXiv 摘要和正式 DOI 题录；未打开期刊全文 | 给出趋势滤波的分段多项式结构和 minimax 理论；可解释为何高阶 \(L_1\) 比固定短窗多项式更自适应 | 真实函数的 \(k\) 阶导数具有有界变差等条件 | 批处理理论；在线需要递推/固定滞后版本 |
| A13 | *Impulsive Noise Mitigation in Powerline Communications Using Sparse Bayesian Learning* | Jing Lin, Marcel Nassar, Brian L. Evans | 2013 | IEEE JSAC | [10.1109/JSAC.2013.130702](https://doi.org/10.1109/JSAC.2013.130702) | 全文 PDF | 与“已知导频/已知激励 + 稀疏干扰联合估计”最相似；利用 null/pilot tones 估计脉冲噪声，并可联合检测 | 干扰在时间域稀疏；导频或其他观测对干扰提供额外约束 | SBL 需迭代；论文给出多种复杂度—性能折中，但非固定时间实时保证 |
| A14 | *Closed-loop Identification Revisited* | Urban Forssell, Lennart Ljung | 1999 | Automatica | [10.1016/S0005-1098(99)00022-9](https://doi.org/10.1016/S0005-1098(99)00022-9) | 正式摘要页；全文 PDF 被 ScienceDirect 拒绝 | 直接支持闭环下输入与噪声/扰动相关会影响开放环估计；本问题中 \(y\) 可能受力控回路影响，不能把 \(d\) 当外生误差 | 需明确闭环结构、噪声模型、参考输入和工具变量 | PEM/IV 可在线递推；普通最小二乘可能有系统偏差 |

### 部分相关但未作为已证实结论的文献

- Green–Silverman 的样条专著、Wahba 的样条理论：本次未核验到可直接打开的正式页面，因此不把具体定理写入结论。
- Särkkä–Hartikainen 的 LFM 结果已由 PMLR 正式页面核验；switching LFM 存在，但我未查到其用于人体接触刚度在线估计的证据。
- SBL 与 IRLS 的关系：本次查到的是 SBL 的稀疏脉冲估计和 PLQ/Huber 的迭代优化；未查到一篇同时对本问题的 \(k,d\) 联合模型给出严格等价定理的论文。

## 2. 最接近的已有应用及差异

1. **机器人超声呼吸补偿**  
   Yuan 与 Lee 的工作：[DOI 10.1109/AIM55361.2024.10637095](https://doi.org/10.1109/AIM55361.2024.10637095)。其对象是腹部机器人超声，使用自适应控制补偿伪周期呼吸扰动；摘要显示验证主要使用海绵腹部模型，不是人体 \(k,d\) 联合辨识，也没有 \(Ld=w\) 的稀疏创新模型。

2. **机器人/人工超声接触力控制**  
   Tele-operated robotic lung ultrasound 研究报告机器人正常呼吸下平均接触力 \(10.48\pm2.72\) N，人工操作为 \(9.52\pm1.02\) N，闭气时 \(10.13\pm0.38\) N：[PMC8442628](https://pmc.ncbi.nlm.nih.gov/articles/PMC8442628/)。其 \(d\) 实际上是呼吸造成的接触力变化，采用的是机构被动顺应/力控制，不是稀疏状态估计。

3. **手术机器人组织识别/刚度估计**  
   [Teleoperated Surgical Robot with Adaptive Interactive Control Architecture for Tissue Identification](https://pmc.ncbi.nlm.nih.gov/articles/PMC10603885/) 中存在接触力、位移和 stiffness estimation，但 \(d\) 未被建模为 Lévy、稀疏创新或 latent force；更接近交互控制中的局部组织参数估计。

4. **通信导频 + 稀疏脉冲干扰**  
   Lin–Nassar–Evans 是方法结构上最接近者：已知 null/pilot 信息与未知稀疏脉冲噪声联合估计。但通信信道模型、脉冲干扰和人体接触力的动力学含义不同，不能直接转移其统计假设。

**结论：未查到已有工作同时把 A1–A6 的结构用于“机器人超声接触刚度 \(k\) 在线估计 + 人体运动扰动 \(d\) 的稀疏动力学建模”。** 查到的是相邻的力控制、呼吸补偿、组织识别和通信稀疏干扰估计。

## 3. B 组：真实人体数据

| 来源 | 对象和条件 | 实测量和数值（原文数字） | 能否用于检验 |
|---|---|---|---|
| Dhyani 等，*A pilot study to precisely quantify forces applied by sonographers while scanning* | 10 名健康志愿者、12 名超声医师、47 次腹部检查 | 平均接触力 \(8.2\pm4.3\) N，范围 \(1.2\)–\(36.5\) N；高 BMI 组平均 \(10.5\) N，正常 BMI 组 \(7.9\) N；最大力均值分别 \(25.3\) 和 \(17.4\) N。 [正式 DOI](https://doi.org/10.3233/WOR-172611) | 可检验力幅值分布和个体差异；未给出原始时序、采样率、上升时间和频谱，不能直接检验 \(Ld=w\) |
| Tele-operated low-cost robotic lung ultrasound | 人体肺部超声，正常呼吸/闭气 | 机器人正常呼吸 \(10.48\pm2.72\) N，人工 \(9.52\pm1.02\) N，闭气 \(10.13\pm0.38\) N；文中说明力变化反映呼吸 | 可检验呼吸诱导力变化的量级；未给扰动时序、频谱、咳嗽/惊跳数据 |
| Yuan–Lee，机器人超声呼吸补偿 | 腹部干预场景；实验主要为海绵腹部模型 | 摘要没有给出力幅值、上升时间或频谱；只说明呼吸扰动和伪周期补偿 | 可证明该任务存在呼吸扰动；不能检验人体扰动统计模型 |
| Dhyani 等系统综述中汇总的超声压力研究 | 腹部、盆底、肌肉超声；多种人体研究 | 例：Flavell 等记录范围 \(0.88\)–\(5.26\) N，采样率 60 Hz；Schlosser 等机器人经腹超声 \(3.15\)–\(6.01\) N；Dall’Alba 等记录 \(3\)–\(14\) N。 [综述页面](https://www.mdpi.com/1424-8220/25/11/3415) | 可作为力幅值和采样率参考；多数不是连续原始数据，不能直接估计 \(d\) 的创新稀疏度 |
| 呼吸诱导腹部/器官位移 | 文献中的腹部运动建模 | 一个被引用的正常呼吸例子：S-I 位移 \(8.6\) mm、A-P 位移 \(1.1\) mm、呼吸频率 \(0.23\) Hz；页面明确这些是预定义正弦 profile，而非本研究原始测量。 [页面](https://www.tandfonline.com/doi/pdf/10.1080/02656736.2022.2091799) | 只能作为数量级和频率先验，不能作为人体接触力扰动的直接验证 |
| 胸壁/上腹壁咳嗽动力学 | 16 名健康志愿者；OEP 同时测胸壁和上腹壁运动、胸腔/腹腔压力、口腔气流 | 摘要确认测了这些变量，但没有给出可直接引用的幅值、上升时间或频谱。 [PubMed](https://pubmed.ncbi.nlm.nih.gov/22144580/) | 研究设计非常适合检验“连续呼吸 + 短时事件”，但本次未获得可下载原始数据 |
| Garg & Ross，人体垂直振动传递 | 12 人，8 男 4 女，1–50 Hz；足部输入位移 \(0.003\)–\(0.02\) inch | 测量人体频响并用 16 质量块模型拟合幅值和相位；未直接报告接触阻抗参数。 [Duke 正式记录](https://scholars.duke.edu/publication/683514) | 可支持“人体在 1–50 Hz 有多模态动态响应”，不能直接代替腹壁接触阻抗 |
| *In vivo intra-abdominal pressure of healthy participants* | 20 名健康参与者，呼吸、咳嗽、Valsalva | 公共数据页报告 36 个变量、60 个观测；数据为腹内压力，不是探头力或腹壁位移；许可证为自定义条款，页面显示可公开访问但部分文件需申请。 [DOI 10.57745/BLKGI0](https://doi.org/10.57745/BLKGI0) | 可检验呼吸/咳嗽事件是否连续、重尾或突发；不能直接验证 \(F=ky+d\) |
| PhysioNet PEEP respiratory dataset | 呼吸功能、胸部相关生理数据 | PhysioNet 页面说明 CC BY 4.0，公开文件包含 EIT、压力/流量相关目录；总解压大小 18.3 GB，未下载。采样率需读 README/数据文件，本次页面未给出。 [正式数据页](https://www.physionet.org/content/respiratory-dataset/1.0.0/) | 可作为呼吸动力学数据源，但不含探头接触力；不建议为本任务直接下载整包 |
| Structured Light Plethysmography dataset | 胸壁左右运动、胸腹呼吸贡献 | 数据集页面确认有胸壁和腹部呼吸贡献变量，但本次未核实采样率和许可证细节。 [数据页](https://pure.qub.ac.uk/en/datasets/structured-light-plethysmography-in-neuromuscular-disorders-proje/) | 可能用于检验 \(d\) 的低频连续部分；不能检验接触力耦合 |

### B 组明确缺口

- **机器人超声人体咳嗽、惊跳、突然移动的探头接触力：未查到同时公开原始力时序、采样率、上升时间和频谱的研究。**
- **人体腹部接触阻抗在 1–50 Hz 的直接实测质量、阻尼、刚度：未查到与探头—腹壁法向接触完全对应的公开数据。**
- 找到的 1–50 Hz 人体研究主要是全身振动传递或人体段落模型，不能直接当作腹壁接触阻抗。

## 4. C 组：对候选思路的挑错

### 4.1 最大问题不是“是否稀疏”，而是 \(k\) 与 \(d\) 的分解可辨识性

模型为

\[
F=ky+d+e.
\]

若 \(d\) 没有限制，则任意 \(k'\) 都可令

\[
d'=d+(k-k')y,
\]

得到完全相同的 \(F\)。因此必须依靠 \(L\) 惩罚、先验或外部激励区分二者。

若 \(L y\) 很小、很稀疏，或者 \(y\) 在局部窗口内可被低阶多项式很好解释，则改变 \(k\) 只会给 \(d\) 增加很小的惩罚，\(k\) 会发生系统性漂移。广义 TV 文献只说明给定正则化后的解具备样条/稀疏创新结构，并不自动保证 \(k\) 可辨识。[Unser–Fageot–Ward](https://doi.org/10.1137/16M1061199)

实际需要类似以下条件：

\[
\|L(y\Delta k)\|_{\text{稀疏代价}}
\]

必须明显大于允许的扰动代价，且残差化后的激励 \(y\) 仍有足够能量。20.8 Hz 正弦探测有利于此，但前提是：

- 探测信号真正进入接触位移；
- 接触刚度在该频率附近没有严重变化；
- 力传感器、控制器和人体运动不会把该频率成分同步吸收到 \(d\)；
- 闭环输入与扰动不相关，或使用正确的闭环辨识方法。

闭环识别文献明确指出，反馈会使输入与噪声/扰动相关，直接开放环估计可能有偏。[Forssell–Ljung](https://doi.org/10.1016/S0005-1098(99)00022-9)

### 4.2 “高斯 + 稀疏”未必是人体运动的正确生成模型

更可能存在三类成分：

1. **呼吸：** 低频、准周期、相位和幅值随时间变化，未必是平稳高斯过程。
2. **人体自主运动：** 可能持续数百毫秒到数秒，并非单点稀疏创新；其 \(Ld\) 可能是一段连续的大值，而非几个孤立脉冲。
3. **咳嗽、惊跳、探头滑移：** 可能是重尾、成簇、状态切换或有明显方向性的事件。

因此，单一 Huber \(\rho\) 可能只是在做鲁棒拟合，并不等于正确的“高斯背景 + 稀疏驱动”概率模型。Student-t 更接近重尾；sum-of-norms 更接近稀疏状态跳变；switching LFM 更接近不同运动状态切换，但三者的统计含义不同。

### 4.3 \(d\) 可能不是外生扰动

探头运动由力控制器产生，人体运动也可能因探头压力、疼痛或接触变化而改变。因此 \(d\) 可能与 \(y\)、控制输入和测量噪声相关。此时把 \(d\) 当作独立创新会产生偏差；必须记录控制输入、参考力、机器人位姿，并考虑工具变量、闭环 PEM 或正交化估计。

### 4.4 静态 \(k\) 可能过于简单

人体软组织具有：

- 位移相关的非线性；
- 速率相关的黏弹性；
- 接触面积和探头姿态变化；
- 压力导致的组织预应力变化；
- 迟滞和局部滑移。

所以 \(F=ky+d\) 的 \(k\) 可能只是短时局部斜率。如果窗口内 \(k\) 变化，算法会把刚度变化误认为 \(d\)，或反过来。

### 4.5 在线性限制

广义 Kalman、PLQ、趋势滤波都能做到随样本数近似线性或固定窗口计算，但：

- 完整 smoothing 使用未来数据，不是严格因果在线；
- Huber/TV/SBL 通常需要迭代；
- switching LFM 还要估计离散状态；
- 实际“在线”需要固定滞后、warm start、迭代上限和最坏时间预算。

因此现有论文支持“可构造低复杂度递推/固定滞后算法”，不支持直接宣称机器人实时认证。

## 5. 最值得尝试的机制（不超过 10 行）

1. **部分线性残差化 + 20.8 Hz 已知探测**：最直接估计 \(k\)，可减少 \(d\) 的平滑偏差。  
2. **二阶或三阶 \(L_1\)/Huber 趋势滤波**：比固定短窗三次投影更能适应局部运动和突变。  
3. **PLQ/广义 Kalman 固定滞后平滑**：把 \(d\) 写成低阶状态空间过程，Huber 或稀疏转移惩罚，在线成本可控。  
4. **暂不优先采用完整 GP-LFM 或 DML**：前者计算和模型选择较重，后者主要解决 nuisance 偏差而非本问题的物理可辨识性。  
5. **未核实内容**：人体探头接触力在咳嗽/惊跳下的原始频谱、上升时间；腹壁 1–50 Hz 直接接触阻抗；以及已有工作是否已完成 \(k,d\) 联合在线估计。