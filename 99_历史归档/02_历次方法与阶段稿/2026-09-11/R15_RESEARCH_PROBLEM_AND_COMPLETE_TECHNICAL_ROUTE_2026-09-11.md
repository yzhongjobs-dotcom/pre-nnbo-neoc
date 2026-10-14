nnbo：
先进行环境辨识，然后把求阻抗转化为最优控制问题，再用single critic nn解hjb方程，得到最优阻抗


# R15 的研究问题与技术路线

## 1. 先把核心问题说清楚

我们重新梳理一下 R15 的技术路线。

机器人在人体软组织上连续扫描时，接触环境并不是固定的。不同位置的组织刚度、黏弹性、厚度和摩擦条件可能不同，肌肉状态、呼吸和身体微动也会改变局部接触状态。因此，可以先用一个参数向量

\[
\theta=(\beta_h,\beta_k,\beta_\tau,\beta_\mu,\beta_t)
\tag{1}
\]

描述局部组织的平衡载荷–压入关系、黏性加载、松弛时间、动摩擦和切向黏性牵引。

对于一个确定的环境参数 \(\theta\)，如果接触模型、优化目标、约束、时间域和 critic 表达都已经固定，就存在一个对应的局部最优价值函数，也存在一组对应的最优 critic 权重。直观上可以写成

\[
\theta
\longrightarrow
W^*(\theta)
\longrightarrow
\text{最优阻抗策略}.
\tag{2}
\]

所以最直观的理解是：**环境参数变了，最优 critic 权重应该跟着变，最终的最优阻抗也应该跟着变。**

不过，R15 中真正进入 HJB 的不只是 \(\theta\)。边界还有多久到达、扫描速度怎样变化、接触足迹多大、未来参数有多大不确定性，都会改变有限时域最优解。因此，后面更准确地使用固定维度描述符

\[
\eta=
\left(
\theta_0,ar\theta_f,r_b,eta_v,eta_a,eta_{\cal E},T_f
\right),
\tag{3}
\]

其中 \(\theta_0\) 是当前环境参数，\(\bar\theta_f\) 是未来环境的名义参数，\(r_b\) 是边界到任务中心的距离，\(\beta_v\) 和 \(\beta_a\) 描述扫描速度与接触足迹变化，\(\beta_{\cal E}\) 编码参数不确定性，\(T_f\) 是有限优化时域。

因此，R15 中更准确的关系是

\[
\eta
\longrightarrow
Z^*(\eta)=
\left(W^*(\eta),\{u_i^*,\lambda_i^*\}\right)
\longrightarrow
u_g^*(t,x)
\xrightarrow{\dot g_a=u_g^*}
g_a(t),
\tag{4}
\]

其中 \(Z^*\) 是完整的 HJB–KKT 最优根，\(W^*\) 是其中的 critic 权重，\(u_g^*\) 是最优阻抗变化率，\(g_a\) 才是机器人实际施加的阻抗。

这里要修正一个容易产生的误解：**critic 权重 \(W^*\) 并不直接等于刚度 \(K^*\)。** 它先决定价值函数及其梯度，再由 Hamiltonian 极小化得到最优阻抗变化率 \(u_g^*\)，最后通过 \(\dot g_a=u_g^*\) 得到实际阻抗。

所以，整个研究真正要解决的是：

> **当机器人连续扫描、未来接触环境发生变化时，如何提前得到新的环境描述符 \(\eta_f\)，并快速得到对应的最优根 \(Z^*(\eta_f)\) 和最优阻抗，而不是在任务中心进入新区域以后才重新从头求解 HJB。**

这条路线自然分成两个问题：

1. **未来环境怎么提前知道？**
2. **知道未来环境以后，新的最优阻抗怎么快速得到？**

空间预瞄解决第一个问题，NNBO/single critic 与 HJB 定义“什么是最优阻抗”，ANEOC 或直接学习参数到 critic 的映射解决第二个问题。

## 2. 第一个问题：未来环境参数怎么提前知道

### 2.1 为什么探头本身能够提供空间预瞄

超声探头不是一个点，它和组织之间存在有限接触足迹。令 \(s_c(t)\) 表示任务中心在扫描路径上的位置，\(s_b(t)\) 表示局部组织边界的位置，定义

\[
r_b(t)=s_b(t)-s_c(t),
\qquad
\dot r_b(t)=-v_{\rm rel}(t).
\tag{5}
\]

如果接触足迹在扫描方向上为 \(x\in[-a(t),a(t)]\)，那么边界在 \(r_b=a\) 时先进入探头前缘，在 \(r_b=0\) 时才到达任务中心。于是，从前缘接触新区域到任务中心到达新区域之间，存在剩余预瞄时间

\[
T_{\rm prev}(t)=
\inf\left\{T\ge0:
r_b(t)-\int_t^{t+T}v_{\rm rel}(\tau)d\tau=0
\right\}.
\tag{6}
\]

**空间预瞄不是非接触预测。** 新区域已经作用在探头前缘上，所以传感器已经能看到变化；只是这个变化还没有作用到我们定义的任务中心。这段时间正好可以用来估计未来环境、更新最优解并提前调整阻抗。

### 2.2 力和力矩为什么能够提供边界信息

六维 F/T 传感器测到的是整个接触斑上的合力和合力矩，而不是直接测得接触斑形状。传感器 wrench 需要先变换到任务中心处的接触坐标系：

\[
\mathbf F_c=\mathbf R_{cs}\mathbf F_s,
\qquad
\mathbf M_c=
\mathbf R_{cs}\mathbf M_s
-(\mathbf R_{cs}\mathbf r_{sc})\times\mathbf F_c.
\tag{7}
\]

完成重力、惯性、偏置和力臂补偿后，取法向力 \(F_n\) 和绕横向轴的俯仰力矩 \(M_y\)。它们分别是接触压力沿扫描方向的零阶和一阶空间矩：

\[
F_n=\int_{-a}^{a}\lambda(x)dx,
\qquad
M_y=\int_{-a}^{a}x\lambda(x)dx,
\tag{8}
\]

其中 \(\lambda(x)\) 是沿探头圆柱轴积分后的法向线载荷密度。

当未来区域从前缘进入时，相对于当前均匀区域的基线，力和力矩增量可以写成

\[
\widetilde{\Delta F}_n
=\Delta q\,A(r_b;a,\chi)+\varepsilon_F,
\qquad
\widetilde{\Delta M}_y
=\Delta q\,B(r_b;a,\chi)+\varepsilon_M,
\tag{9}
\]

其中

\[
A(r;a,\chi)=L\int_r^a\phi(x;a,\chi)dx,
\qquad
B(r;a,\chi)=L\int_r^a x\phi(x;a,\chi)dx.
\tag{10}
\]

\(\phi(x;a,\chi)>0\) 是根据探头几何、接触模型和均匀材料标定得到的足迹灵敏度核，\(\chi\) 包括压入量、姿态、预载和扫描速度等当前工作条件。接触半宽 \(a\) 和接触斑相关量不是 F/T 传感器直接给出的，而是由接触模型、位姿、载荷和历史信息共同估计得到。

在理想无噪声情况下，取力矩与力的比值

\[
\rho(r)=\frac{B(r)}{A(r)},
\qquad
\rho'(r)=
\frac{L\phi(r)}{A(r)}[\rho(r)-r]>0.
\tag{11}
\]

因此 \(\rho(r)\) 对边界距离严格单调，可以先求

\[
r_b=
ho^{-1}\left(
\frac{\widetilde{\Delta M}_y}
{\widetilde{\Delta F}_n}
\right),
\qquad
\Delta q=
\frac{\widetilde{\Delta F}_n}
{A(r_b;a,\chi)}.
\tag{12}
\]

这就是空间预瞄最核心的观测原理：**力告诉我们“变化了多少”，力矩进一步告诉我们“变化主要出现在接触足迹的什么位置”，两者联合把边界位置 \(r_b\) 和响应差异 \(\Delta q\) 分开。**

**为什么不能只用力？** 因为单独的力只给出 \(\Delta qA(r_b)\) 这个乘积。较大的响应差异和较短的进入长度，可能产生与较小响应差异和较长进入长度相同的力变化。力矩提供了额外的空间信息。

### 2.3 实际在线时怎样估计 \(r_b\) 和 \(\Delta q\)

式 (12) 解释了为什么可以估计，但不能直接作为实际求解器。边界刚进入前缘时，\(\widetilde{\Delta F}_n\) 很小，比值容易放大噪声；接触半宽、传感器偏置和组织运动也可能同时变化。因此，实际方法使用一段短时历史进行联合估计。

定义预瞄状态

\[
z_k=(r_{b,k},\Delta q_k,a_k,b_{F,k},b_{M,k}),
\tag{13}
\]

并用

\[
\begin{aligned}
r_{b,k+1}&=r_{b,k}-T_sv_{{\rm rel},k}+w_{r,k},\\
\Delta q_{k+1}&=\Delta q_k+w_{q,k},\\
a_{k+1}&=a_k+w_{a,k},\\
b_{F,k+1}&=b_{F,k}+w_{F,k},\\
b_{M,k+1}&=b_{M,k}+w_{M,k}
\end{aligned}
\tag{14}
\]

描述短时间内的状态变化。对应的观测模型是

\[
y_k=
\begin{bmatrix}
\widetilde{\Delta F}_{n,k}\\
\widetilde{\Delta M}_{y,k}
\end{bmatrix}
=
\begin{bmatrix}
\Delta q_kA(r_{b,k};a_k,\chi_k)+b_{F,k}\\
\Delta q_kB(r_{b,k};a_k,\chi_k)+b_{M,k}
\end{bmatrix}
+\varepsilon_{y,k}.
\tag{15}
\]

在窗口 \({\cal I}_k=\{k-N_e+1,\ldots,k\}\) 内，在线移动窗估计器求解

\[
\begin{aligned}
\hat z_{{\cal I}_k}
=\arg\min_{z_{{\cal I}_k}}
\quad&
\|z_{k-N_e+1}-\bar z_{k-N_e+1}\|_{P^{-1}}^2\\
&+\sum_{j\in{\cal I}_k}
\|y_j-h(z_j,\chi_j)\|_{R_j^{-1}}^2\\
&+\sum_{j=k-N_e+1}^{k-1}
\|z_{j+1}-f_z(z_j,v_{{\rm rel},j})\|_{Q_j^{-1}}^2,
\end{aligned}
\tag{16}
\]

同时满足

\[
0\le r_{b,j}\le a_j,
\quad
a_j\in[\underline a_j,\overline a_j],
\quad
\Delta q_j\in{\cal Q},
\quad
\varepsilon_{y,j}\in{\cal V}_j,
\quad
w_j\in{\cal W}_j.
\tag{17}
\]

这个 `argmin` 给出在线使用的名义估计 \(\hat r_b\)、\(\widehat{\Delta q}\) 和 \(\hat a\)。上一时刻的解可以作为 warm start，再用有界 Gauss–Newton 或 SQP 求解。

但名义估计不等于真实值已经被唯一确定。为了后续安全约束，还要保留所有与观测和误差界一致的状态集合

\[
{\cal Z}_k=
\left\{
z_{{\cal I}_k}:
(14),(15),(17)\ \text{在声明误差界内成立}
\right\},
\qquad
{\cal X}_{\rm pre}
=\Pi_{r_b,\Delta q,a,\chi}{\cal Z}_k.
\tag{18}
\]

根据这个集合，可以进一步得到保守的剩余预瞄时间

\[
\underline T_{\rm prev}=
\inf_{z_{{\cal I}_k}\in{\cal Z}_k}
\inf\left\{T\ge0:
r_{b,k}-\int_{t_k}^{t_k+T}v_{\rm rel}(\tau)d\tau=0
\right\}.
\tag{19}
\]

**为什么这里不用一个普通 EKF 就结束？** EKF 可以给式 (16) 提供初值，但本问题同时有边界距离、接触半宽、慢变偏置和明确的几何范围约束，而且后续安全性需要一个包含真实状态的不确定集合。因此，R15 使用移动窗优化得到名义值，再用集合成员外逼近保留仍未消除的不确定性。

### 2.4 怎样从 \(\Delta q\) 得到未来环境参数 \(\theta_f\)

这里还有一个必须说明的问题：\(\Delta q\) 只是当前工作条件下的有效响应差异，它不等于完整的未来材料参数 \(\theta_f\)。例如，仅靠一次法向自然扫描，可能能够看到平衡刚度发生变化，却不能同时唯一辨识松弛时间和摩擦参数。

离线均匀材料标定提供映射

\[
q_{\rm eff}(\theta;\chi),
\tag{20}
\]

它表示环境参数 \(\theta\) 在工作条件 \(\chi\) 下会产生怎样的有效接触响应。给定当前材料集合 \({\cal E}_0\)、未来参数先验 \({\cal E}_{\rm prior}\) 和移动窗估计得到的 \(\widehat{\Delta q}_i\)，未来环境的名义参数通过

\[
\begin{aligned}
(\hat\theta_0,\hat\theta_f)
=\arg\min_{\substack{
\theta_0\in{\cal E}_0\\
\theta_f\in{\cal E}_{\rm prior}}}
\quad&
\sum_{i\in{\cal I}_k^{\rm info}}
w_i\left[
q_{\rm eff}(\theta_f;\chi_i)
-q_{\rm eff}(\theta_0;\chi_i)
-\widehat{\Delta q}_i
\right]^2\\
&+\lambda_\theta
\|L_\theta(\theta_f-\theta_{f,{\rm prior}})\|^2
\end{aligned}
\tag{21}
\]

在线求得。这里 \({\cal I}_k^{\rm info}\) 只保留响应高于噪声阈值且具有信息量的样本，\(w_i\) 反映各样本的可信度，正则项只用于在不可辨识方向上选择一个稳定代表值。数值上使用上一窗口的参数解 warm start，再采用有界 SQP 或 trust-region least-squares 求解。

与此同时，所有与预瞄数据相容的未来参数构成

\[
\begin{aligned}
{\cal E}_f=
\bigg\{\theta_f\in{\cal E}_{\rm prior}:&\ \exists\,
z_{{\cal I}_k}\in{\cal Z}_k,
\ \exists\,\theta_0\in{\cal E}_0,\\
&\left|
q_{\rm eff}(\theta_f;\chi_i)
-q_{\rm eff}(\theta_0;\chi_i)
-\Delta q_i(z_{{\cal I}_k})
\right|
\le\epsilon_{q,i},\ \forall i
\bigg\}.
\end{aligned}
\tag{22}
\]

式 (22) 可以用区间分支定界或带认证误差的外逼近计算，保证没有把仍与观测相容的未来参数误删掉。

因此，空间预瞄最终不是只输出一个看起来很精确的 \(\hat\theta_f\)，而是同时输出

\[
\boxed{
\hat\theta_f\ \text{用于名义性能计算},
\qquad
{\cal E}_f\ \text{用于鲁棒安全约束}.
}
\tag{23}
\]

如果自然扫描不能辨识某些参数，这些参数就继续留在 \({\cal E}_f\) 中，而不是被强行估成唯一值。这也是为什么后续更准确地使用 \(\eta\) 和 \({\cal E}_\eta\)，而不是只使用一个点参数 \(\theta_f\)。

到这里，第一部分已经完成：

\[
\{F/T,\text{位姿},\text{接触模型},\text{历史状态}\}
\longrightarrow
\{\hat r_b,\widehat{\Delta q},\underline T_{\rm prev}\}
\longrightarrow
\{\hat\theta_f,{\cal E}_f\}.
\tag{24}
\]

## 3. 环境参数怎样进入后面的最优控制问题

得到 \(\theta_f\) 以后，还不能直接把它送进 critic。必须先说明：这个环境参数如何改变接触力、机器人运动和阻抗调整过程。

当边界位于接触足迹内部时，当前材料和未来材料同时存在。材料参数场写成

\[
\theta(x,t)=
\theta_0+H(x-r_b(t))(\theta_f-\theta_0).
\tag{25}
\]

为了描述组织的黏弹记忆，引入随材料在接触足迹中输运的内部状态 \(\bar q(x,t)\)：

\[
\partial_t\bar q-v_{\rm rel}\partial_x\bar q
=k_v(\delta_x;\beta_k)\dot\delta_x
-\tau_v(\delta_x;\beta_\tau)^{-1}\bar q.
\tag{26}
\]

法向载荷密度和接触 wrench 由

\[
\lambda_n(x,t)=
\phi(x;a,\chi)
\left[h_\infty(\delta_x;\theta(x,t))+\bar q(x,t)\right],
\tag{27}
\]

\[
F_n=L\int_{-a}^{a}\lambda_n(x,t)dx,
\qquad
M_y=L\int_{-a}^{a}x\lambda_n(x,t)dx
\tag{28}
\]

得到。

这里不能只把边界两侧的两个集总接触力按面积比例相加，因为旧材料正在离开、新材料正在进入，而黏弹状态仍携带历史信息。分布式模型保留了这个输运过程。

为了让 HJB 可以在线处理，用固定 Galerkin 基将 \(\bar q\) 降阶为

\[
\bar q_N(x,t)=\sum_{j=1}^{N_q}\xi_j(t)\psi_j(x),
\tag{29}
\]

从而得到固定维度模型

\[
M_\psi\dot\xi=
\left[C_\psi-T_\psi(\theta_0,\theta_f,r_b)\right]\xi
+B_\psi(\theta_0,\theta_f,r_b,e)\dot\delta_n
+b_{\rm in}+d_\xi.
\tag{30}
\]

降阶误差满足离线标定得到的界

\[
\|d_\xi\|\le\epsilon_{\xi,N},
\qquad
|d_{F,N}|\le\epsilon_{F,N},
\qquad
\|d_{t,N}\|\le\epsilon_{t,N}.
\tag{31}
\]

如果目标工作域内不能建立这些误差界，就必须提高降阶阶数或更换接触模型，不能继续用当前低维模型求 HJB。

## 4. NNBO 在这里到底做什么

### 4.1 先确定控制对象：调的是阻抗变化率，不是直接调力

R15 沿用 NNBO 中“用 single critic 近似价值函数，再由价值梯度恢复最优控制”的核心思路，但不能原样沿用其接触对象和动作定义。

在这里，接触力是机器人与组织相互作用产生的输出，而不是 HJB 可以直接指定的动作。机器人真正能够调节的是实际施加阻抗 \(g_a\)，并且阻抗不能瞬间跳变。因此定义

\[
\dot g_a=u_g,
\qquad
g_{\min}\le g_a\le g_{\max},
\qquad
|u_g|\le\bar r_g.
\tag{32}
\]

其中 \(u_g\) 是 HJB 的控制动作。完整状态为

\[
x=(e,v,\xi,g_a,E_T),
\tag{33}
\]

包括任务误差、机器人速度、接触内部状态、实际阻抗和能量罐状态。统一动力学写成

\[
\dot x=f_0(t,x;\eta)+G(t,x;\eta)u_g+d,
\qquad d\in{\cal D}_N.
\tag{34}
\]

能量罐用于限制时变刚度可能注入的能量：

\[
\dot E_T=
\sigma v_w^\top D_w(g_a)v_w
-\frac12e_w^\top\dot K_w(g_a,u_g,t)e_w
-\rho_E,
\qquad
E_{\min}\le E_T\le E_{\max}.
\tag{35}
\]

因此，后面求得的不是一个可以瞬间切换的 \(K^*\)，而是一条满足增益范围、变化率和能量条件的最优阻抗轨迹。

### 4.2 不是一开始就必须使用变量最优阻抗

如果存在一个固定阻抗，能够对未来参数集合中的所有环境保持安全，而且性能损失可以接受，就没有必要进入复杂的 HJB–ANEOC 路线。定义公共安全固定增益集合

\[
{\cal G}_{\rm robust}^{\rm fix}
=\bigcap_{\eta\in{\cal E}_\eta}
{\cal G}_{\rm safe}^{\rm fix}(\eta).
\tag{36}
\]

若存在 \(g_0\in{\cal G}_{\rm robust}^{\rm fix}\)，并且

\[
\sup_{\eta\in{\cal E}_\eta}
\left[
J_{\rm fix}(g_0,\eta)-J_{\rm var}^*(\eta)
\right]
\le\epsilon_J,
\tag{37}
\]

则保持这个固定增益。只有公共固定增益不存在，或者它相对变量最优阻抗的性能损失过大时，才继续求变量最优阻抗。

### 4.3 把“什么是最优阻抗”写成有限时域 HJB

由于阻抗必须在边界到达以前完成调整，优化问题是有限时域的。令 \(t_f=t_\times+T_{\rm post}\)，运行代价为

\[
\begin{aligned}
\ell(t,x,u_g;\eta)=
&\ e^\top Q_e e+v^\top Q_vv
+q_F(F_{n,N}-F_d)^2\\
&+q_s[\,m_{\rm slip}-\mu_{\rm lb}F_{n,N}
+\|F_{t,N}\|\,]_+^2\\
&+(g_a-g_{\rm soft})^\top Q_g(g_a-g_{\rm soft})
+u_g^\top R_gu_g.
\end{aligned}
\tag{38}
\]

它同时考虑轨迹误差、接触力、滑移风险、保持柔顺和阻抗调整速度。价值函数定义为

\[
\begin{aligned}
V^*(t,x;\eta)=
\min_{u_g(\cdot)}\quad&
\Phi_f(x(t_f);\eta)
+\int_t^{t_f}\ell(\tau,x,u_g;\eta)d\tau\\
\text{s.t.}\quad&
\dot x=f_0+Gu_g+d,\\
&u_g\in{\cal U}_{\rm safe}(t,x;{\cal E}_\eta),\\
&x(t_f)\in{\cal X}_f(\eta).
\end{aligned}
\tag{39}
\]

其中 \({\cal U}_{\rm safe}\) 已经包含力、接触、滑移、阻抗变化率、执行器和能量罐约束。对应的 HJB 为

\[
-\partial_tV^*=
\min_{u_g\in{\cal U}_{\rm safe}}
\left\{
\ell+
\nabla_xV^{*\top}(f_0+Gu_g)
\right\},
\qquad
V^*(t_f,x;\eta)=\Phi_f(x;\eta).
\tag{40}
\]

**为什么不是无限时域 HJB？** 因为边界还有多久到达本身会影响阻抗应该多快变化。有限时域问题可以同时描述越界前的阻抗调整和越界后的短时接触，稳态无限时域问题不能直接表达这个截止时间。

### 4.4 single critic 怎样把 HJB 变成有限维求解

直接求完整的 \(V^*(t,x;\eta)\) 很困难，所以选择固定基函数 \(\varphi(t,x)\)，用 single critic 近似

\[
V_N(t,x;W,\eta)=
\Phi_f(x;\eta)
+(t_f-t)W^\top\varphi(t,x).
\tag{41}
\]

这个形式自动满足终端条件。给定 \(W\) 后，在每个配置节点 \((t_i,x_i)\) 上，最优动作由

\[
u_i^*(W,\eta)=
\arg\min_{A_iu_i\le b_i}
\left\{
u_i^\top R_gu_i
+\nabla_xV_N^\top G_i u_i
+\ell_{0,i}
\right\}
\tag{42}
\]

得到。由于 \(R_g\succ0\)，该 Hamiltonian 子问题严格凸，所以给定价值梯度以后最优动作唯一，不需要再训练一个独立 actor。

把各节点的 Bellman 残差、动作最优性条件和互补约束叠在一起，可以写成一个有限维根：

\[
F(Z,\eta)=0,
\qquad
Z=(W,\{u_i,\lambda_i\}_{i=1}^{n_W}).
\tag{43}
\]

其中 \(\lambda_i\) 是约束乘子。更具体地，\(F\) 包含

\[
F(Z,\eta)=
\begin{bmatrix}
\{R_i(W,u_i;\eta)\}\\
\{2R_gu_i+G_i^\top\nabla_xV_N+A_i^\top\lambda_i\}\\
\{\phi_{\rm FB}(\lambda_i,b_i-A_iu_i)\}
\end{bmatrix},
\tag{44}
\]

其中

\[
R_i=\partial_tV_N+\ell_i+
\nabla_xV_N^\top(f_{0,i}+G_iu_i),
\qquad
\phi_{\rm FB}(\lambda,s)=
\sqrt{\lambda^2+s^2}-\lambda-s.
\tag{45}
\]

因此，NNBO/single critic 在 R15 中完成的是：

\[
\boxed{
\eta
\longrightarrow
\text{有限时域最优控制问题}
\longrightarrow
F(Z,\eta)=0
\longrightarrow
W^*,u_g^*
\longrightarrow
g_a^*.
}
\tag{46}
\]

## 5. 如果没有 ANEOC，会发生什么

当环境描述符从 \(\eta_k\) 变成 \(\eta_{k+1}\) 时，HJB 的动力学、代价和约束参数也随之变化。因此，需要从

\[
F(Z,\eta_k)=0
\tag{47}
\]

转到新的问题

\[
F(Z,\eta_{k+1})=0.
\tag{48}
\]

如果没有参数迁移机制，最直接的方法就是对式 (48) 重新求根。可以冷启动，也可以把旧根当作普通 warm start，但求解器仍需要在新问题上完成必要的 Jacobian 分解和 Newton 迭代。

严格来说，这里“重新求解 HJB”指的是重新求解离散后的投影 HJB–KKT 根 \(F(Z,\eta_{k+1})=0\)，而不是每次重新进行一遍符号上的 HJB 推导。

所以这里真正的问题变成：

> 已知旧环境下的最优根 \(Z^*(\eta_k)\)，空间预瞄又提前给出了新的 \(\eta_{k+1}\)，能不能利用两者之间的函数关系，直接预测新根在哪里？

这正是 ANEOC 要做的事情。

## 6. ANEOC 怎样从旧最优解得到新最优解

### 6.1 先求最优解对环境描述符的导数

在一个局部正则分支上，最优根是描述符的函数：

\[
Z^*=Z^*(\eta),
\qquad
F(Z^*(\eta),\eta)=0.
\tag{49}
\]

对上式关于 \(\eta\) 求导，根据隐函数定理得到

\[
F_Z\frac{dZ^*}{d\eta}+F_\eta=0,
\tag{50}
\]

因此

\[
C_\eta
=\frac{dZ^*}{d\eta}
=-F_Z^{-1}F_\eta.
\tag{51}
\]

如果只关注 critic 权重，这个矩阵中对应 \(W\) 的部分就是

\[
P_\eta
=\frac{dW^*}{d\eta}.
\tag{52}
\]

它告诉我们：当未来环境、边界时间或不确定性沿某个方向发生小变化时，最优 critic 权重和完整最优根应该向什么方向变化。

### 6.2 用一阶 Taylor 展开先预测新根

空间预瞄给出

\[
\Delta\eta=\eta_{k+1}-\eta_k.
\tag{53}
\]

于是 ANEOC 先计算

\[
Z_{\rm pred}
=Z^*(\eta_k)+C_\eta(\eta_k)\Delta\eta.
\tag{54}
\]

对应到 critic 权重，就是

\[
W_{\rm pred}
=W^*(\eta_k)+P_\eta(\eta_k)\Delta\eta.
\tag{55}
\]

从 Taylor 展开可以看到

\[
Z^*(\eta_{k+1})
=Z^*(\eta_k)
+C_\eta(\eta_k)\Delta\eta
+{\cal O}(\|\Delta\eta\|^2).
\tag{56}
\]

也就是说，一阶环境变化已经由 \(C_\eta\Delta\eta\) 处理，理想情况下剩下的主要是二阶余项。这里不是凭经验猜一个更新方向，而是利用 HJB–KKT 根 \(F(Z,\eta)=0\) 和隐函数定理直接得到最优解对环境参数的导数。

### 6.3 再用 corrector 消除剩余误差

Taylor 预测只能把解送到新最优点附近，不能保证已经严格满足新环境下的 HJB–KKT 条件。因此，以 \(Z_{\rm pred}\) 为初值，在 \(\eta_{k+1}\) 下进行 Newton 校正：

\[
Z^{j+1}=Z^j-
F_Z(Z^j,\eta_{k+1})^{-1}
F(Z^j,\eta_{k+1}),
\qquad
Z^0=Z_{\rm pred}.
\tag{57}
\]

当根残差、原始约束、互补条件和策略误差都低于阈值时，得到新的可接受根 \(Z^*(\eta_{k+1})\)。

所以，ANEOC 的作用可以概括为：

> **先利用最优解对环境描述符的导数，把旧最优根沿着最优解分支推进到新根附近；再用少量 Newton/corrector 把预测点拉回新环境下真正的 HJB–KKT 根。**

### 6.4 为什么不直接重新求解 HJB

重新求解式 (48) 当然可行，而且必须作为基线保留。使用 ANEOC 的原因是：空间预瞄给出的计算时间有限，而相邻环境下求解的是同一维度、同一 critic 基、同一配置节点和同一约束排列的参数化根。旧根及其导数已经包含新问题的一阶信息，把新问题完全当作无关问题重新求解会浪费这部分信息。

但这里有一个重要条件：\(F_Z\) 的分解和灵敏度 \(C_\eta\) 应该在预瞄事件之前或后台提前维护。这样，收到新的 \(\Delta\eta\) 后，关键路径主要是式 (54) 的矩阵–向量乘法和少量校正。如果等到检测到新环境以后才开始完整计算 \(-F_Z^{-1}F_\eta\)，那么它与从旧根开始做第一步 Newton 在局部一阶上是等价的，不能仅凭式 (56) 就声称更快。

因此，ANEOC 的优势不是“产生了一个 Newton 完全没有的方向”，而是：**把可提前完成的 Jacobian 分解和灵敏度计算移出事件后的关键路径。** 最终是否更快，需要与冷求解、普通 warm-start Newton 和直接有限时域 NMPC 在相同精度与硬件上比较。

### 6.5 ANEOC 在什么情况下成立

ANEOC 不是对任意大的环境突变都有效。它要求：

- 新旧问题使用相同的状态维数、基函数、配置节点和约束排列；
- 新根仍位于当前可以连续跟踪的局部正则分支；
- \(F_Z\) 非奇异且条件数没有恶化到不可接受；
- \(\Delta\eta\) 足够小，使预测点仍落在 Newton 可收敛区域内。

如果 \(\Delta\eta\) 较大，可以把它拆成多个同伦小步；如果活动约束改变，可以使用半光滑 Newton 和广义 Jacobian。但如果发生失去接触、冲击、粘滑模式切换或明显分支跳变，就不能继续沿当前曲线外推，而应切换安全模式或重新求解。

## 7. 空间预瞄与 ANEOC 怎样接起来

现在可以把前后两部分直接连起来：

\[
\boxed{
\begin{aligned}
&\text{探头前缘先进入新区域}\\
&\quad\longrightarrow
\text{F/T、位姿和历史状态发生可观测变化}\\
&\quad\longrightarrow
\text{MHE 求 }\hat r_b,\widehat{\Delta q},\underline T_{\rm prev}\\
&\quad\longrightarrow
\text{参数反演求 }\hat\theta_f\text{ 和 }{\cal E}_f\\
&\quad\longrightarrow
\text{形成未来描述符 }\eta_{k+1}\text{ 和 }\Delta\eta\\
&\quad\longrightarrow
Z_{\rm pred}=Z^*(\eta_k)+C_\eta\Delta\eta\\
&\quad\longrightarrow
\text{Newton/corrector 得到 }Z^*(\eta_{k+1})\\
&\quad\longrightarrow
u_g^*(t,x)\longrightarrow g_a(t)\\
&\quad\longrightarrow
\text{任务中心到达新区域以前完成阻抗调整}.
\end{aligned}}
\tag{58}
\]

因此，R15 解决的不只是“环境变化后实时调阻抗”，而是：**在未来环境真正作用到任务中心以前，提前估计未来接触、提前得到对应的最优解，并提前把阻抗调过去。**

## 8. 第二条路线：直接学习环境到最优 critic 的映射

前面的主路线仍然依赖局部延拓。如果环境变化很大，或者新旧描述符不够接近，ANEOC 的一阶预测可能不再可靠。由此可以得到第二种思路：离线直接学习

\[
\eta\longrightarrow W^*(\eta).
\tag{59}
\]

具体地，在离线参数域内选择样本 \(\{\eta_i\}\)，分别求得对应的最优根 \(Z_i^*\)、critic 权重 \(W_i^*\) 和灵敏度

\[
P_{\eta,i}=
\frac{dW^*}{d\eta}(\eta_i).
\tag{60}
\]

训练网络 \(\widehat W_\psi(\eta)\) 时，不只拟合离散的 critic 权重，还同时拟合其对环境描述符的导数：

\[
\begin{aligned}
{\cal L}(\psi)=
&\sum_i
\|\widehat W_\psi(\eta_i)-W_i^*\|^2\\
&+\lambda_J\sum_i
\left\|
\frac{\partial\widehat W_\psi}{\partial\eta}(\eta_i)
-P_{\eta,i}
\right\|_F^2.
\end{aligned}
\tag{61}
\]

这样，网络学习的不只是参数空间中的离散点，还学习了最优解分支在这些点附近的变化方向。在线得到未来 \(\eta_f\) 后，可以直接计算

\[
W_{\rm pred}=\widehat W_\psi(\eta_f).
\tag{62}
\]

随后根据式 (42) 恢复节点动作和乘子，组成预测根 \(Z_{\rm pred}\)，再用式 (57) 的 HJB–KKT corrector 检查和校正。

**为什么网络输出后仍然需要 corrector？** 因为网络拟合误差可能使 \(F(Z_{\rm pred},\eta_f)\neq0\)，也可能违反活动约束或安全裕量。网络适合提供更大范围内的快速初值，但不能仅凭训练损失替代 HJB 根残差和可行性检查。

所以两条路线的关系是：

1. **当前主路线：空间预瞄 + NNBO/single critic + ANEOC。** 它利用局部灵敏度从当前已接受根快速迁移到相邻未来根，并用 corrector 恢复精确性。
2. **第二路线：空间预瞄 + 导数增强的 \(\eta\rightarrow W^*\) 学习。** 它用离线全局近似直接给出 critic 初值，再接 HJB–KKT corrector 和安全检查。

第二条路线可能减少对新旧环境局部相邻性的依赖，但需要足够覆盖的离线样本、跨参数域的一致根分支和严格的域外检测。因此，它目前更适合作为替代预测器或后续扩展，而不是直接替代当前 ANEOC 主路线。

## 9. 算出最优阻抗以后，还要判断能不能及时执行

求出新的 \(Z^*\) 并不等于已经解决任务。估计、参数反演、预测校正、通信和真实阻抗渲染都需要时间。定义事件后的非渲染延迟

\[
\tau_{\rm nr}=
\tau_{\rm filter}
+\tau_{\rm est}
+\tau_{\rm pc}
+\tau_{\rm comm}
+\tau_{\rm margin}.
\tag{63}
\]

令 \({\cal R}_g(g_a(t_p),T;{\cal U}_{\rm safe})\) 表示在阻抗变化率、能量和执行器约束下，时间 \(T\) 内能够达到的实际阻抗集合；令

\[
{\cal G}_{\rm cross}({\cal E}_\eta)=
\bigcap_{\vartheta\in{\cal E}_\eta}
{\cal G}_{\rm cross}(\vartheta)
\tag{64}
\]

表示对全部未来环境都安全的越界阻抗集合。只有当

\[
{\cal R}_g\left(
g_a(t_p),
\underline T_{\rm prev}-\tau_{\rm nr};
{\cal U}_{\rm safe}
\right)
\cap
{\cal G}_{\rm cross}({\cal E}_\eta)
\neq\varnothing
\tag{65}
\]

时，新的阻抗策略才允许执行。

这一步回答了最后一个自然问题：即使未来最优阻抗已经算出来，如果剩余时间不够，或者真实执行器无法在速率与能量限制内到达它，那么这个“最优阻抗”仍然不能解决研究问题。

如果式 (65) 不成立，首先降低扫描速度以增加预瞄时间；如果仍然不可行，则保持保守鲁棒阻抗、停止或撤回。

## 10. 最终在线执行顺序

整个方法在线运行时，严格按照下面的顺序执行：

1. 读取六维 F/T、机器人位姿、速度和历史接触状态，并把传感器 wrench 转换到任务中心处的接触坐标系。
2. 根据当前均匀区域的基线，形成力和力矩增量 \(\widetilde{\Delta F}_n\)、\(\widetilde{\Delta M}_y\)。
3. 求解式 (16) 的移动窗 `argmin`，得到 \(\hat r_b\)、\(\widehat{\Delta q}\)、\(\hat a\) 和状态可行集合 \({\cal X}_{\rm pre}\)。
4. 根据式 (19) 计算保守预瞄时间 \(\underline T_{\rm prev}\)。
5. 求解式 (21) 的参数反演 `argmin`，得到名义未来参数 \(\hat\theta_f\)，同时计算未来参数集合 \({\cal E}_f\)。
6. 将未来参数、边界位置、剩余时间、扫描速度、接触足迹和不确定性编码成 \(\eta_{k+1}\) 与 \({\cal E}_\eta\)。
7. 先检查固定阻抗是否已经安全且性能足够；如果满足式 (36)–(37)，则不进入变量阻抗路线。
8. 否则，用提前维护的 \(C_\eta\) 和式 (54) 预测新的 HJB–KKT 根。
9. 用式 (57) 进行必要的 Newton/corrector，并检查根残差、约束、活动集和 Jacobian 条件数。
10. 从接受的 \(Z^*\) 中恢复最优阻抗变化率 \(u_g^*(t,x)\)，再检查能量条件和式 (65) 的越界前可达性。
11. 所有条件通过后，由 \(\dot g_a=u_g^*\) 提前调整实际阻抗；任何估计不一致、接触模式变化或求解失效都会取消当前自适应更新。

## 11. 最后把研究问题压缩成一句话

R15 的完整逻辑是：

\[
\boxed{
\text{空间预瞄提前得到未来环境}
\longrightarrow
\text{NNBO/HJB 定义未来最优阻抗}
\longrightarrow
\text{ANEOC 或学习映射快速得到新 critic}
\longrightarrow
\text{在任务中心越界前安全地建立新阻抗}.
}
\tag{66}
\]

**最核心的问题仍然是：对于不断变化的未来接触环境，如何快速得到对应的最优 critic 和最优阻抗。** 空间预瞄回答“未来环境怎么提前知道”；NNBO/single critic 回答“最优阻抗如何由 HJB 和 critic 表示”；ANEOC 或直接学习 \(\eta\rightarrow W^*(\eta)\) 回答“环境变化以后怎样快速得到新的最优解”；最后的可达性与能量检查回答“这个解能不能在边界到达以前真正执行”。

当前主路线是 **空间预瞄 + NNBO/single critic + ANEOC predictor–corrector**。直接学习 \(\eta\rightarrow W^*(\eta)\) 是逻辑上成立的第二路线，但只有在离线覆盖、根分支一致性、域外检测和 HJB corrector 均得到验证以后，才能替代 ANEOC 预测器。


更自然的做法：
样本积累、参数估计和阻抗更新可以并行滚动进行，不必等“参数完全估准”以后才开始调。