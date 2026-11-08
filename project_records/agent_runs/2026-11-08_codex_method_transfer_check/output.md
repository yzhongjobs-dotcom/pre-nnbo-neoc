**Q1．DOB 鲁棒稳定与预测 Q：机制成立，但必须检查完整闭环。**

**已核实条件。** Sun et al., IFAC PID 2018，§3、PDF p.3、式(11)–(12)，针对稳定乘性不确定性给出 \(\|\Delta Q\|_\infty\le1\)。设计保证应采用严格裕度：
\[
P=P_n(1+\Delta),\qquad \|\Delta Q\|_\infty<1 .
\]
这是 DOB **内环**的小增益充分条件，要求相关块稳定、互联适定且逆模型滤波后可实现。[原文](https://skoge.folk.ntnu.no/prost/proceedings/PID-2018/0031.PDF)

**完整闭环区别。** Shim 的作者综述 *Disturbance Observer*，PDF p.3、式(1)，分母为
\[
P_n(1+PC)+Q(P-P_n).
\]
代入上述乘性模型，令 \(L=P_nC\)，可直接推得
\[
(1+L)\left[1+\Delta\frac{L+Q}{1+L}\right].
\]
因此，在名义内部稳定等条件下，完整系统对应的小增益充分条件是
\[
\left\|\Delta\frac{L+Q}{1+L}\right\|_\infty<1 .
\]
后两式是**由原式推导**，不能把内环 \(\|\Delta Q\|<1\) 当作全部稳定性证明。该综述 pp.4–5、Theorem 1 重述 Shim–Jo 2009：名义闭环内部稳定、所有候选对象最小相位、固定相对阶及已知正高频增益区间，且
\[
p_f(s)=s^\nu+a_{\nu-1}s^{\nu-1}+\cdots+a_1s+(g/g^*)a_0
\]
对整个增益区间 Hurwitz，则充分小的 Q 时间常数保证稳定；这是另一类高带宽定理。[作者原文](https://arxiv.org/pdf/2101.02859)

**预测影响与时延先例。** 数学上，纯相位变化不改变内环乘积模值，但会改变 Nyquist 轨迹及完整闭环中的 \(L+Q\)；实际相位超前还可能增加高频增益。理想 \(e^{sT}\) 非因果，不能直接当作可实现 Q。Kato–Muis–Ohnishi，*Automatika* 2006，§3、pp.6–8、Fig.6，已将无时延逆模型 DOB 与串联 Smith predictor 组合；其名义式(15)给出 \(G_{ry}=CP_n e^{-sT_0}/(1+CP_n)\)，失配鲁棒性主要由实验支持。Natori–Ohnishi 2008 的具体定理本轮 **not verified**。[Smith–DOB 原文](https://hrcak.srce.hr/en/file/7208)

**适配：**可作为框架；须把 31 ms 伺服、IMM 预测及刚度误差回路共同纳入，冻结参数条件不能自动证明在线切换稳定。

**Q2．最小代价辨识：优化“为控制所需的信息”，有滚动时域扩展。**

Bombois et al., *Automatica* 2006，作者稿 §4、pp.4–5、Constraint 1／式(11)：固定数据长度 \(N\)，优化外加激励谱 \(\Phi_r\)，最小化激励引起的输入、输出扰动功率：
\[
\min_{\Phi_r\ge0}J_r
=\frac1{2\pi}\int_{-\pi}^{\pi}
\left(\alpha_y|G_0S_{id}|^2+\alpha_u|S_{id}|^2\right)\Phi_r\,d\omega ,
\]
约束为辨识不确定集内**每个对象**均满足新控制器的鲁棒性能要求：
\[
J(G,\hat C,W_l,W_r)\le1,\quad
\forall G\in D(\hat\theta_N,P_\theta).
\]
也讨论固定谱最小化 \(N\)，以及时间与扰动代价折中。条件：真实对象在可辨识模型集中、实验充分信息化、PEM 渐近协方差有效；设计计算使用未知参数的先验近似。置信集覆盖概率 \(\beta\) 不是确定性误差界。[作者稿 §§2–4、6](https://perso.uclouvain.be/michel.gevers/PublisMig/LCID_final.pdf)

在线扩展确实存在：Larsson et al., *European Journal of Control* 29 (2016)，*An application-oriented approach to dual control with excitation for closed-loop identification*，出版社摘要及引言明确为约束闭环生成低干扰激励，并使用滚动时域简化；§5介绍两种 output-error 系统的 MPC 控制器。**全文约束方程和在线保证 not verified**。Hjalmarsson 2005 作者页面确认文献，但原文具体优化式本轮 **not verified**。[Larsson 出版社摘要／节选](https://www.sciencedirect.com/science/article/pii/S0947358016300024)、[Hjalmarsson 作者页](https://hjalmars.netlify.app/publication/hjalmarsson-04-a/)

**适配：**以力波动作为激励成本，以允许的增益误差／控制性能作为精度约束；将固定 20.8 Hz 的幅值作为决策变量，是该思想的任务化应用，尚非现成保证。

**Q3．力限幅安全滤波：找到一篇直接相关工作，接口条件需保留。**

Wang et al., *Guarding Force*，UCL 记录为 **RA-L 2025 accepted/in press**；以下方程核自 **arXiv 2024 v1**，不假定终稿编号相同。§V、PDF pp.4–5，式(10)、(16)、(21)：
\[
h_e=K_e^{pri}(x-x_r^{pri})+\delta_e,\qquad b_{fc}=h_e^{max}-h_e,
\]
\[
\dot q^*=\arg\min_{\dot q}\tfrac12\|\dot q-\dot q_{nom}\|^2,\quad
-K_e^{pri}J\dot q-z_2-\sigma(t)+l b_{fc}\ge0 .
\]
实际刚度未知；需要已知先验刚度和静止位置，力传感器使残差 \(\delta_e\) 可测，微分器估计其导数。误差界依赖已知的残差二阶导数界和初始误差界。保证还要求初始安全、可行控制及相应正则性。主要实现为**关节速度接口**，UR3e 底层跟踪期望速度；§VII另给力矩接口扩展。这不是刚度区间端点式的一步位置 governor。[原文 §IV–VII](https://arxiv.org/pdf/2405.04859)、[RA-L 接收记录](https://discovery.ucl.ac.uk/id/eprint/10216915/)

**适配：**是直接可借鉴的力 CBF；位置接口与 31 ms 延迟需要额外预测裕度及采样间约束，不能沿用连续速度接口保证。

**Q4．Yuen 文件辨认：三个文件均不是该论文，只引用了它。**

本地核实：`E18.txt` 是 Cao et al. 2019（Yuen 为 PDF p.13、ref.[7]）；`G11.txt` 是 Xia et al. 2024（p.8、ref.[24]）；`C03.txt` 是 Wu et al. 2026（p.17、ref.[36]）。因此改查 Yuen 作者公开稿。

Yuen et al., **T-RO 2010**，作者稿 §IV、p.5、式(12)：
\[
f_a=f_d+K_f(f_d-f_e)+K_v(\hat{\dot x}_e-\dot x)
+b\hat{\dot x}_e+m\hat{\ddot x}_e .
\]
前馈来自 **3D 超声组织位置**；§V、pp.5–6、式(13)，以 Fourier 模型和 EKF 得到速度、加速度，接触后停止更新，避免受机器人变形影响的位置污染自由组织运动估计。它没有使用 \(F/\hat k-y\)，也没有分析该刚度误差诱发的前馈闭环；§VI-A、p.6只说明低增益考虑组织刚度不确定性。

带宽论据：§II、p.3、式(7)–(9)，临界阻尼设计下，拒绝约 15 Hz 组织运动需约 57 Hz 自然频率；§VII、p.7明确概括为约 **3.8 倍**。这是该模型与阻尼／扰动抑制指标的结果，不是普遍常数。[Yuen 作者稿](https://www.biorobotics.harvard.edu/pubs/2010/journal/yuen_ForceControl_TRo_2010.pdf)

**适配：**支持利用运动前馈绕开反馈带宽限制；其独立成像测量无法替代我们的刚度误差回路分析。

**Q5．快速新颖性扫描：未核实到完整满足 (a) 或 (b) 的论文；(c) 有必须比较的近邻。**

本地材料加 **10 条定向网络检索式**，结果如下：

- **Takakura et al., IEEE OJIES 2026，*Bilateral Teleoperation With Autoelimination of Probing Signal…***：作者接受稿 §IV、p.6、式(41)–(43)，
  \[
  f^{prob}=W^{prob}(e^{iden})\delta(t),\quad
  e^{iden}=\frac{g_{lpf}}{s+g_{lpf}}(\hat f_F^{RFOb}-\hat f_F)^2.
  \]
  以残差阈值开／关探测，环境变化后重新注入。条件为一维双边遥操作、RLS 接触模型及噪声相关阈值；不是由所需刚度置信精度连续优化幅值。**宽泛的“辨识好就减小激励”已存在。**已读公开作者稿，出版社终版 **not verified**。[作者稿全文](https://www.researchgate.net/publication/403145904_Bilateral_Teleoperation_With_Autoelimination_of_Probing_Signal_for_Environmental_Impedance_Identification_Considering_Operability)
- **Wang et al., TASE 2024，*Task Space Compliant Control and Six-Dimensional Force Regulation…***：本地 `F09.txt`，p.3655、§II-C、式(18)，前馈预测由给定期望位姿和 twist 推进；不是在线刚度反演的组织运动测量。
- **PG-DGPRL，Computers & Industrial Engineering 2026**：出版社摘要涉及不确定性感知前馈与延迟补偿；是否按在线参数区间调度前馈增益，**not verified**。[出版社摘要](https://www.sciencedirect.com/science/article/abs/pii/S0360835226000392)

**判定／适配：**本轮未找到并核实 (a)、(b) 或严格意义的 (c)；不能据此声明不存在。可继续考察的具体差异是：**刚度误差如何闭合预测前馈回路，以及如何由该回路的性能要求反推辨识精度与探测成本**。