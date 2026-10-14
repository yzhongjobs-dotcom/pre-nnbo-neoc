# 方法：预瞄条件下可执行最优阻抗的运输

> 论文草稿，2026 年 9 月 11 日。本版本将已完成的 Stage 1 模型与 Stage 1 之后的设计决策改写为论文式 Method。本文不声称已经获得实验验证。需要下一推导阶段才能确定的方程将被明确标出，而不会用历史 NNBO 方程填补。

## 工作论点

针对连续机器人超声扫查，本文构建一种预瞄条件下的阻抗控制方法：将有限接触足迹转化为对即将到来的组织转变的不确定性有界估计，将圆柱–组织接触响应与已施加控制器增益的动力学嵌入统一的最优控制模型，沿估计的未来接触方向运输参数相关的最优 critic，并在任务中心到达转变前渲染所得阻抗，同时满足模型有效性、时序与无源性条件。

## 3. 方法

### 3.1 问题表述

我们考虑一个在规定表面曲线持续移动、并保持压缩接触的机器人持超声探头。任务中心固定为扫查任务采用的控制点。令 \(s_c(t)\) 表示其弧长坐标，\(s_b(t)\) 表示局部组织转变的坐标，并定义

\[
r_b(t)=s_b(t)-s_c(t)
\tag{1}
\]

为有符号转变距离。探头沿扫查方向的实际接触半宽为 \(a(t)\)。因此，\(0<r_b<a\) 表示 spatial preview interval：接触足迹的前缘部分已经进入未来区域，而任务中心仍位于当前区域。该定义不意味着在前缘发生物理接触前即可感知。

在预瞄更新时间 \(t_p\)，Stage 1 估计器给出

\[
\mathcal I_p(t_p)=
\left\{
\hat r_b,\;
\Delta\hat q,\;
\hat T_{\mathrm{prev}},\;
\widehat{\mathcal M}_{c,1}
\right\},
\tag{2}
\]

其中，\(\Delta\hat q\) 是在当前预载和扫查运动下可观测的有效接触响应对比，

\[
\widehat{\mathcal M}_{c,1}
=
\left\{
\mathcal M_c(\eta_{c,1}):
\eta_{c,1}\in\widehat{\mathcal E}_1
\right\}
\tag{3}
\]

是未来区域的不确定性有界接触模型集合，\(\hat T_{\mathrm{prev}}\) 是到 \(r_b=0\) 前剩余时间的保守估计。物理接触参数

\[
\eta_c=
\left\{
\eta_e,\eta_v,\tau_v,\mu_k,c_t,R^*,L_{\mathrm{eff}}
\right\}
\tag{4}
\]

不假设可由单次转变事件完全辨识。

控制器必须将 \(\mathcal I_p\) 和当前机器人–接触状态映射为可执行阻抗。我们用

\[
g_a=
\begin{bmatrix}
k_n & k_t & d_n & d_t
\end{bmatrix}^{\!\top},
\tag{5}
\]

表示已施加增益向量，其中 \(k_n,d_n\) 分别为局部接触坐标系下的法向刚度和阻尼，\(k_t,d_t\) 为相应的切向量。目标不只是计算一个未来增益目标。所选增益必须在任务中心到达边界前已经成为底层控制器实际使用的增益：

\[
g_a(t_b^-)\in
\mathcal G_{\mathrm{safe}}
\!\left(\widehat{\mathcal E}_1\right),
\qquad
t_b=\inf\{t\ge t_p:r_b(t)=0\}.
\tag{6}
\]

集合 \(\mathcal G_{\mathrm{safe}}\) 纳入正刚度和正阻尼、执行器极限、接触维持、滑移极限以及时变阻抗实现的能量预算。

### 3.2 方法概览

该方法由四个耦合模块组成。首先，有限足迹估计器利用接触 wrench 的零阶和一阶空间矩估计即将到来的转变，并构建 \(\widehat{\mathcal M}_{c,1}\)。其次，边界渐进的圆柱–组织模型表征接触足迹跨越边界时当前与未来接触响应如何共存。第三，预瞄条件下的最优控制问题将已施加阻抗作为被控对象状态的一部分，并计算其可允许变化率。第四，自适应邻近极值最优控制（adaptive neighboring-extremal optimal control, ANEOC）更新将参数相关 critic 从当前接触模型运输至未来模型，并使所得增益轨迹通过时序、不确定性和无源性门。

执行顺序为

\[
\begin{aligned}
&\text{wrench 与机器人历史}
\longrightarrow
\{\hat r_b,\widehat{\mathcal E}_1,\hat T_{\mathrm{prev}}\}\\
&\longrightarrow
\text{边界渐进接触动力学}
\longrightarrow
\text{未来策略描述符}\\
&\longrightarrow
\text{critic 预测和校正}
\longrightarrow
\text{可允许的增益率策略}\\
&\longrightarrow
\text{在 }r_b=0\text{ 前施加阻抗}。
\end{aligned}
\tag{7}
\]

该顺序至关重要。接触 wrench 是物理输出，增益率命令是最优控制 action，而 \(g_a\) 是已施加控制器状态。这些对象不可互换。

### 3.3 边界渐进的圆柱–组织动力学

#### 3.3.1 转变进程

令 \(k(x;a,\chi)>0\) 为 spatial-preview 模型引入的已标定力灵敏度核，其中 \(\chi\) 包含局部预载和几何工作条件。我们将分配给未来区域的、力敏感接触足迹比例定义为

\[
\alpha(r_b,a,\chi)=
\begin{cases}
0, & r_b\ge a,\\[3pt]
\displaystyle
\frac{\int_{r_b}^{a}k(x;a,\chi)\,\mathrm dx}
{\int_{-a}^{a}k(x;a,\chi)\,\mathrm dx},
& -a<r_b<a,\\[12pt]
1, & r_b\le -a.
\end{cases}
\tag{8}
\]

因此，未来区域尚未进入足迹时 \(\alpha=0\)；足迹跨越边界时 \(0<\alpha<1\)；整个足迹进入未来区域后 \(\alpha=1\)。对于移动边界，

\[
\dot r_b=v_b-v_s.
\tag{9}
\]

式（8）采用生成预瞄观测的同一空间灵敏度，从而避免仅为控制便利而引入无关的转变调度。\(r_b\)、\(a\) 和 \(\chi\) 的不确定性诱导出 \(\alpha\) 的区间 \(\mathcal A(t)\)，并将其传播至下游接触集合和策略集合。

#### 3.3.2 降阶法向接触模型

对于均匀局部区域 \(j\in\{0,1\}\)，其中 \(j=0\) 表示当前区域、\(j=1\) 表示未来区域，Stage 1 给出

\[
w_j=h_{\infty,j}(\delta_n;\eta_{e,j})+q_{v,j},
\qquad
\dot q_{v,j}
=
k_{v,j}(\delta_n;\eta_{v,j})\dot\delta_n
-\frac{q_{v,j}}{\tau_{v,j}},
\tag{10}
\]

其中，\(w_j\) 是单位探头长度法向载荷，\(\delta_n\) 是压缩压入量，\(q_{v,j}\) 是最小黏弹内部状态。均匀模型在已标定工作域内局部单调：

\[
\frac{\partial h_{\infty,j}}{\partial\delta_n}
+k_{v,j}(\delta_n)>0.
\tag{11}
\]

当接触足迹跨越边界时，精确接触响应具有空间分布，通常不能用两组参数的算术平均代替。因此，我们定义降阶转变算子

\[
\mathcal M_{c,\mathrm{tr}}
=
\mathscr T_{\alpha}
\!\left(
\mathcal M_c(\eta_{c,0}),
\mathcal M_c(\eta_{c,1})
\right)
+d_{\mathrm{tr}},
\tag{12}
\]

其中，\(\mathscr T_{\alpha}\) 由两个子区域的压力核加权聚合获得，\(d_{\mathrm{tr}}\) 有界地表示被舍弃的空间模态。一个最小双分支实现为

\[
F_n
=L_{\mathrm{eff}}
\left[
(1-\alpha)h_{\infty,0}(\delta_n)
+\alpha h_{\infty,1}(\delta_n)
+Q_{v,0}+Q_{v,1}
\right]
+d_F,
\tag{13}
\]

\[
\dot Q_{v,j}
=
\omega_j(\alpha)
k_{v,j}(\delta_n)\dot\delta_n
-\frac{Q_{v,j}}{\tau_{v,j}}
+b_j(\delta_n,Q_{v,j},\eta_{c,j})\dot\omega_j
+d_{q,j},
\tag{14}
\]

其中 \(\omega_0=1-\alpha\)，\(\omega_1=\alpha\)。输运项 \(b_j\dot\omega_j\) 描述材料进入或离开降阶接触子域；必须由分布接触模型推导或在标定域内辨识，不能未经验证而设为零。若实验表明双分支降阶无法对转变响应给出有界表示，则式（12）必须由更高阶空间或有限元 surrogate 实现，并撤回低维 HJB 的主张。

法向运动学保持为

\[
\dot\delta_n
=-\mathbf n^\top(\mathbf v_p-\mathbf v_s)+d_\delta,
\tag{15}
\]

并满足

\[
F_n\ge0,\qquad g_n\ge0,\qquad F_ng_n=0.
\tag{16}
\]

主连续时间推导限制于持续接触分支 \(g_n=0\)。接触丢失作为混杂模式转变处理，并调用后备控制器。

#### 3.3.3 滑动过程中的切向接触

对于 \(\|\mathbf v_t\|\ge v_{\min}\)，接触足迹的未来与当前部分通过各自的法向载荷份额对切向 wrench 作出贡献。令 \(F_{n,j}\) 表示式（13）中的两个相应项，则滑动模型为

\[
\mathbf F_t
=-
\sum_{j=0}^{1}
\left[
\mu_{k,j}F_{n,j}
+\omega_j c_{t,j}\|\mathbf v_t\|
\right]
\frac{\mathbf v_t}
{\sqrt{\|\mathbf v_t\|^2+\varepsilon_v^2}}
+\mathbf d_t.
\tag{17}
\]

正则化 \(\varepsilon_v\) 仅在滑动分支内使用。低速粘着、扫查反向和 stick–slip 转换均作为独立模式，而不是式（17）的平滑延拓。

总接触力和机器人动力学为

\[
\mathbf F_c=F_n\mathbf n+\mathbf F_t,
\tag{18}
\]

\[
\mathbf M(q)\ddot q
+\mathbf C(q,\dot q)\dot q
+\mathbf g(q)
=
\tau+\mathbf J_c(q)^\top\mathbf W_c+\tau_d.
\tag{19}
\]

式（8）–（19）提供下游最优控制问题所需的物理状态和参数依赖关系。不保留 Hunt–Crossley 状态方程或 force-as-action 的解释。

### 3.4 将可执行阻抗作为增广被控对象状态

令 \(\mathbf e_x\) 和 \(\mathbf e_v\) 分别为在接触坐标系中表达的位姿误差和速度误差。平移阻抗矩阵为

\[
\mathbf K_c(g_a)=
\operatorname{diag}(k_t,k_t,k_n),
\qquad
\mathbf D_c(g_a)=
\operatorname{diag}(d_t,d_t,d_n).
\tag{20}
\]

在符号约定 \(\mathbf e_x=\mathbf x_p-\mathbf x_d\) 下，命令笛卡尔 wrench 为

\[
\mathbf F_{\mathrm{cmd}}
=
\mathbf F_{\mathrm{ff}}
-\mathbf K_c(g_a)\mathbf e_x
-\mathbf D_c(g_a)\mathbf e_v.
\tag{21}
\]

相应关节命令通过所规定的底层接口生成，例如

\[
\tau
=
\tau_{\mathrm{dyn}}
+\mathbf J^\top\mathbf R_c\mathbf F_{\mathrm{cmd}},
\tag{22}
\]

其中 \(\mathbf R_c\) 将接触坐标系下的力映射至机器人基坐标系。实际实现必须包含平台所用的转动阻抗和零空间控制器，但这些固定通道不是最小表述中的优化变量。

底层控制器不能瞬时改变 \(g_a\)。因此，我们采用可允许的增益率

\[
u_g=\dot g_a
\tag{23}
\]

作为最优控制 action，并满足

\[
g_a\in\mathcal G,\qquad
u_g\in\mathcal U_g,\qquad
\|u_g\|_{\infty}\le\bar u_g.
\tag{24}
\]

若硬件暴露的是命令增益 \(g_c\) 而非增益率，则以经辨识的 renderer

\[
\dot g_a=f_g(g_a,g_c),
\tag{25}
\]

替代式（23），并将 \(g_c\) 作为 action。无论采用何种实现，纳入 Markov state 的都是式（21）实际使用的增益，而不是发送给底层控制器的软件命令。

在持续接触和滑动模式内，一个完整增广状态为

\[
x=
\begin{bmatrix}
q^\top&
\dot q^\top&
\delta_n&
Q_{v,0}&
Q_{v,1}&
g_a^\top&
r_b
\end{bmatrix}^{\!\top}.
\tag{26}
\]

将式（8）–（25）结合，当增益率接口为 control-affine 时得到

\[
\dot x
=f_0(x;\eta)+G_g(x;\eta)u_g+d_x,
\tag{27}
\]

其中，\(\eta\) 是下文定义的策略描述符，\(d_x\) 汇总有界模型误差和观测误差。若式（25）中测得的 renderer 对命令并非 affine，则在低维 action 上对下游 Hamiltonian 进行数值最小化；不能仅为恢复闭式策略而强行假设 affine 模型。

### 3.5 从未来接触模型集合到策略描述符

物理参数向量 \(\eta_c\) 与策略描述符 \(\eta\) 的作用不同。前者参数化接触模型；后者必须仅包含在声明精度内决定最优反馈律所需的接触性质组合和预瞄变量。

在工作压入量 \(\bar\delta_n\) 处，名义法向响应量为

\[
\kappa_{e,j}
=
L_{\mathrm{eff}}
\frac{\partial h_{\infty,j}}
{\partial\delta_n}(\bar\delta_n),
\qquad
\kappa_{v,j}
=
L_{\mathrm{eff}}
k_{v,j}(\bar\delta_n),
\tag{28}
\]

以及 \(\tau_{v,j}\)、\(\mu_{k,j}\) 和 \(c_{t,j}\)。我们暂定

\[
\eta=
\left[
\kappa_{e,1},
\kappa_{v,1},
\tau_{v,1},
\mu_{k,1},
c_{t,1},
\hat r_b,
\widehat{\dot r}_b,
\varrho_{\mathcal E}
\right]^{\!\top},
\tag{29}
\]

其中，\(\varrho_{\mathcal E}\) 编码 \(\widehat{\mathcal E}_1\) 在会改变策略的方向上的半径或支撑函数值。式（29）是候选最小描述符，而非假设其中所有材料项都可点辨识。Stage 1 对比 \(\Delta\hat q\) 约束 \(\widehat{\mathcal E}_1\)，并不直接替代 \(\eta\)。

若对于任意两个被映射到同一描述符的可允许未来模型，在状态域 \(\mathcal X\) 上均满足

\[
\sup_{x\in\mathcal X}
\left\|
\pi^*(x;\eta_{c}^{(a)})
-\pi^*(x;\eta_{c}^{(b)})
\right\|
\le\varepsilon_\pi,
\tag{30}
\]

则该描述符在该域内对策略充分。该条件为式（29）提供降维检验。若在预瞄测量下不可区分的两个模型要求不相容的安全阻抗策略，控制器不发布唯一的未来策略，而是保留鲁棒当前增益、降低扫查速度或请求额外信息。

### 3.6 预瞄条件下的最优阻抗问题

期望扫查平衡具有非零法向力 \(F_d\)。因此，我们定义以平衡为中心的变量，而非惩罚绝对接触力。令 \(\tilde x=x-x_s(\eta)\)，其中 \(x_s(\eta)\) 对给定描述符满足接触平衡、期望扫查速度和 renderer 平衡。运行代价可写为

\[
\ell(\tilde x,u_g;\eta)
=
\tilde x^\top Q(\eta)\tilde x
+u_g^\top R_g u_g
+\ell_{\mathrm{slip}}
+\ell_{\mathrm{margin}},
\tag{31}
\]

其中 \(Q\succeq0\)，\(R_g\succ0\)，\(\ell_{\mathrm{slip}}\) 惩罚切向裕度损失，\(\ell_{\mathrm{margin}}\) 惩罚接近接触丢失和力安全边界。硬极限仍保留在可允许状态–action 集合中，而不被软惩罚替代。

对于式（27）中包含 \(r_b\) 和饱和转变变量 \(\alpha\) 的自治增广动力学，无限时域值函数为

\[
V^*(\tilde x;\eta)
=
\inf_{u_g(\cdot)}
\int_0^\infty
\ell(\tilde x(t),u_g(t);\eta)\,\mathrm dt,
\tag{32}
\]

并受式（24）、式（27）以及持续接触/滑动模式假设约束。相应的约束 Hamilton–Jacobi–Bellman 方程为

\[
0=
\min_{u_g\in\mathcal U_g(\tilde x,\eta)}
\left\{
\ell(\tilde x,u_g;\eta)
+\nabla_{\tilde x}V^*(\tilde x;\eta)^\top
\left[f_0(\tilde x;\eta)+G_g(\tilde x;\eta)u_g\right]
\right\}.
\tag{33}
\]

对于内部 action 和式（31）中的二次项，

\[
u_g^*(\tilde x;\eta)
=
-\frac{1}{2}
R_g^{-1}
G_g(\tilde x;\eta)^\top
\nabla_{\tilde x}V^*(\tilde x;\eta).
\tag{34}
\]

当 action 或安全边界激活时，以受约束 Hamiltonian 最小化器替代式（34）。式（33）–（34）的输出是增益率策略，其积分决定未来已施加阻抗。它不是先求最优接触力信号、再将其重新解释为阻抗的做法。

### 3.7 参数相关的单 critic

我们在已声明的紧状态域 \(\mathcal X\) 上用固定且无冗余的基函数逼近值函数：

\[
V_N(\tilde x;\eta)
=W(\eta)^\top\phi(\tilde x),
\tag{35}
\]

其中，基函数对所有接触描述符保持固定，并满足平衡条件 \(V_N(0;\eta)=0\)。将其代入式（33），并在固定测试点集上进行 collocation 或 Galerkin 投影，得到有限维残差方程

\[
F_{\mathrm{HJB}}(W,\eta)=0.
\tag{36}
\]

我们区分精确值函数 \(V^*\)、投影解 \(W_N^*(\eta)\) 和数值迭代 \(W\)。single critic 指式（35）中的单个值函数近似器；它并不意味着一般神经网络参数化拥有唯一的权重向量。为支持局部 continuation，所选参数化必须消除权重对称性，或以其他方式定义局部唯一的投影根。

critic 诱导的增益率策略为

\[
\hat u_g(\tilde x;\eta,W)
=
\arg\min_{u_g\in\mathcal U_g}
H\!\left(
\tilde x,u_g,
\nabla\phi(\tilde x)^\top W,
\eta
\right).
\tag{37}
\]

对于内部 control-affine 情形，式（37）退化为式（34），其中以 \(\nabla\phi^\top W\) 替代 \(\nabla V^*\)。

### 3.8 未来 critic 的 ANEOC 运输

设当前投影 critic 满足

\[
F_{\mathrm{HJB}}
\left(W_N^*(\eta_0),\eta_0\right)=0.
\tag{38}
\]

若 \(F_{\mathrm{HJB}}\) 连续可微，且

\[
F_W
=
\frac{\partial F_{\mathrm{HJB}}}{\partial W}
\tag{39}
\]

在局部解分支上非奇异，则隐式微分给出

\[
C_\eta
:=
\frac{\mathrm dW_N^*}{\mathrm d\eta}
=-F_W^{-1}F_\eta.
\tag{40}
\]

预瞄估计提供 \(\Delta\eta=\eta_1-\eta_0\)。ANEOC 首先预测

\[
W_{\mathrm{pred}}
=
W_N^*(\eta_0)+C_\eta\Delta\eta.
\tag{41}
\]

对于局部二次可微的解映射，

\[
\left\|
W_{\mathrm{pred}}-W_N^*(\eta_1)
\right\|
=O(\|\Delta\eta\|^2).
\tag{42}
\]

随后，带保护的校正器求解未来残差：

\[
W^{(m+1)}
=
W^{(m)}
-\gamma_m
F_W\!\left(W^{(m)},\eta_1\right)^{-1}
F_{\mathrm{HJB}}\!\left(W^{(m)},\eta_1\right),
\tag{43}
\]

直至残差、策略可允许性和条件数容限均得到满足。灵敏度因子在检测到转变前预先维护；若仅在检测后计算式（40），则相对于 Newton warm start 的预期 deadline 优势将消失。约束激活时，以完整局部 primal–dual root 替代 \(W\)，并且式（40）–（43）仅在 active set 不变时成立。active-set 变化、接触模式切换或 \(F_W\) 近奇异时，触发半光滑重求解或后备策略。

对于不确定性集合而非点描述符，式（41）运输策略相关描述符集合的支撑方向或选定顶点。仅当由此产生的策略离散度和 critic 近似误差均落入预设安全预算内时，才发布未来策略。

### 3.9 转变前渲染与可行性门

预测 critic 通过式（37）确定增益率策略。renderer 从测得的已施加增益 \(g_a(t_p)\) 出发，在式（24）约束下积分该策略。端到端时序条件为

\[
\tau_{\mathrm{filter}}
+\tau_{\mathrm{preview}}
+\tau_{\mathrm{transport}}
+\tau_{\mathrm{correct}}
+\tau_{\mathrm{comm}}
+\tau_{\mathrm{render}}
+\tau_{\mathrm{margin}}
\le
\underline T_{\mathrm{prev}}.
\tag{44}
\]

对于静止边界和保守距离估计 \(\underline r_b\)，

\[
v_s
\le
\frac{\underline r_b}
{\tau_{\mathrm{total}}}.
\tag{45}
\]

若名义扫查速度下不满足式（44），系统首先降低 \(v_s\)。若条件仍不可行，则保留已经验证的保守阻抗，或按临床安全状态机撤离。

时变刚度可注入能量。对于阻抗储能项

\[
S_K=\frac{1}{2}\mathbf e_x^\top\mathbf K_c(g_a)\mathbf e_x,
\tag{46}
\]

增益变化对 \(\dot S_K\) 的贡献为

\[
P_K
=
\frac{1}{2}
\mathbf e_x^\top
\dot{\mathbf K}_c
\mathbf e_x.
\tag{47}
\]

正的 \(P_K\) 必须由能量罐或等价的无源层提供。仅当

\[
\int_{t_p}^{t_b}
[P_K(t)]_+\,\mathrm dt
\le
E_{\mathrm{tank}}(t_p)-E_{\min}
+E_{\mathrm{harvest}},
\tag{48}
\]

并同时满足状态、力、滑移与执行器约束时，才施加候选增益轨迹。式（44）和式（48）使计算与渲染成为方法的一部分，而非事后安全检查。

### 3.10 假设与适用范围

连续推导假设：任务中心定义固定、扫查方向已知、持续压缩接触、滑动满足 \(\|\mathbf v_t\|\ge v_{\min}\)、圆柱–平面近似在局部有效，以及单个低维转变位于接触足迹内。平衡和黏弹映射必须在声明的压入量、速度、组织厚度和支撑范围内标定。未来 HJB root 必须位于局部正则分支上，critic 基函数必须定义局部唯一的投影解，且预瞄区间必须超过实测计算与渲染时延。

对于接触丢失、黏附、折皱、超出工作域的大压入、探头端部效应、多个未分辨转变、破坏 wrench-moment 可观测性的边界方向，或在 continuation 步骤中发生的 stick–slip 模式切换，本方法不经修改并不适用。若与预瞄相容的模型集合包含需要不相容安全增益的模型，也不发布唯一的运输策略。这些情形调用保守控制器或更高保真模型，而不是外推局部 HJB 分支。

## 拟议总览图图注

图 X | 预瞄条件下的最优阻抗运输与执行。有限圆柱接触足迹的前缘部分先于任务中心遭遇组织转变。力和力矩历史确定一个不确定性有界的转变距离和未来接触模型集合。边界渐进的圆柱–组织模型将该信息映射为策略描述符。ANEOC 将参数相关的 single-critic 解运输至未来描述符，随后受约束校正器恢复局部 HJB root。所得增益率策略由底层 renderer 积分，且仅在任务中心到达转变前满足不确定性、时序、接触和无源性门时施加。

## 面向作者的术语表

| 规范术语 | 首次使用时的定义 | 不可混用的术语 |
|---|---|---|
| 任务中心 | 扫查/控制任务使用的固定探头点 | 探头中心、声学中心，除非明确相同 |
| 空间预瞄 | 在 \(0<r_b<a\) 时获得的前缘接触信息 | 非接触预测 |
| 未来接触模型集合 \(\widehat{\mathcal M}_{c,1}\) | 与预瞄观测和先验标定一致的模型 | 点估计材料 |
| 物理接触参数 \(\eta_c\) | 圆柱–组织动力学参数 | 策略描述符 \(\eta\) |
| 策略描述符 \(\eta\) | 在 \(\varepsilon_\pi\) 内对最优策略充分的低维量 | \(\Delta q\) |
| 已施加增益 \(g_a\) | 底层阻抗律实际使用的增益 | 命令增益 |
| 增益率 action \(u_g\) | 改变 \(g_a\) 的有界速率 | 力 action |
| 投影最优 critic \(W_N^*(\eta)\) | 有限基 HJB 残差的根 | 精确 HJB 解、在线迭代 |
| ANEOC predictor–corrector | 对参数相关 critic root 的局部运输和校正 | 转变后重新训练 |
| 转变前执行 | 在 \(r_b=0\) 前准备好已施加增益 | 仅完成计算但尚未渲染 |

## 假设或缺失输入

1. 式（14）需要经过标定或推导的边界通量项 \(b_j\)；其精确形式是 Stage 1 后首个尚未完成的物理推导。
2. 式（29）的候选描述符必须通过式（30）的策略充分性检验。在对共同 HJB/NMPC 模型进行数值审查前，不能固定其最终维数。
3. 基函数 \(\phi\)、collocation 域、root 维数、残差容限和 critic-error certificate 尚待指定。
4. 相对于 warm-start Newton 和 preview-conditioned NMPC 的运行时间优势需要 P95/P99 时延测量；本草稿不作此断言。
5. 目标 venue、发布的代码/模型产物以及精确硬件 renderer 接口尚未指定。

## 主张–证据映射

| 方法主张 | 当前依据 | 状态 |
|---|---|---|
| 有限接触足迹在任务中心遭遇前提供信息 | Stage 1 几何与 wrench-moment 推导 | 在 Stage 1 假设下已获解析支持 |
| 未来模型应表示为集合 | Stage 1 可观测性边界 | 已支持 |
| 当增益率和渲染极限重要时，已施加增益必须是状态 | 闭环状态一致性 | 已推导 |
| 在新接触表述中，力不能继续作为 HJB action | 接触力由式（13）–（19）产生，是物理输出 | 已推导 |
| ANEOC 提供二阶局部 predictor | 正则性条件下的隐函数展开 | 有条件支持 |
| ANEOC 快于 warm-start Newton 或 NMPC | 尚无目标系统运行时间数据 | 需要证据 |
| 双分支转变模型可对真实组织转变给出有界表示 | 尚无标定数据 | 需要证据 |
| 所得方法可降低 first-hit force 和接触失败 | 实验尚未完成 | 需要证据 |
