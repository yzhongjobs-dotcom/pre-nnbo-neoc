# R15：研究问题与技术路线

本文给出从接触过程中的空间观测、未来接触模型到 single-critic HJB、ANEOC 和实际阻抗生成的研究问题与方法路线。接触与控制推导见同目录 `technical_appendix.txt`；认证域、预测误差、闭环可行性、局部近最优性和计算优势的完整证明见 `theoretical_proofs.txt`（下称“认证证明”）。性能结论针对指定名义模型、比较域和有限时域；真实组织上的精度与实时性需通过实验验证。

## 1. 总研究问题

**机器人连续扫描空间变化的软组织时，环境参数的改变可能引起接触动力学和所需阻抗策略的相应改变。R15 要解决的是：怎样从接触过程中获得新的环境信息，并在环境变化后快速得到可认证的近最优控制策略与实际阻抗，而不是每次都从头求解完整的最优控制问题。**

不同位置的软组织可能具有不同的有效弹性和摩擦特性；呼吸或组织整体运动还会改变探头与参考表面的相对关系。因此，接触力和力矩会随扫描位置、压入量和相对运动变化。已有机器人超声研究已经将软组织黏弹性参数估计用于可变阻抗控制，也有工作将可变形接触模型用于机器人运动与接触力预测。[1,2] 这些研究分别支撑环境估计和接触感知控制的必要性，但不直接替代 R15 对整条技术路线的证明。

用 $\theta$ 表示环境参数，用 $\mathcal M(\theta)$ 表示相应接触模型；后续最优控制阶段将接触模型、任务状态和约束整理为描述符 $\eta$。在固定 HJB 表达、critic 基和约束布局，并处于同一局部正则解分支时，整体关系可写成

$$
\theta
\longrightarrow
\mathcal M(\theta)
\longrightarrow
\eta
\longrightarrow
Z_{\rm pred}(\eta)
\longrightarrow
\pi_N(t,x;W_{\rm pred},\eta)
\longrightarrow
g_a.
$$

其中 $Z^*(\eta)$ 表示所选有限维 HJB–KKT 问题的局部根；ANEOC 从认证锚点产生候选 $Z_{\rm pred}$，由 $W_{\rm pred}=\Pi_WZ_{\rm pred}$ 恢复反馈增益率。候选通过根误差、闭环可行性和局部性能证书后直接发布；需要时先校正或重新建立锚点。$g_a$ 是由增益率积分得到的实际阻抗参数。

因此，总研究问题可以分成三个连续子问题：

1. 怎样从接触信息估计未来环境参数并建立未来接触模型？
2. 怎样在该接触模型下定义最优阻抗控制问题，并求得其有限维 HJB–KKT 局部根？
3. 环境参数小幅改变后，怎样利用已有根直接生成可认证的近最优阻抗策略，并决定何时校正、重锚或重新求解？

第 3 节展开第一个子问题，第 4—6 节依次完成最优阻抗控制问题定义、有限根更新与实际执行。环境参数确定了接触动力学；最优控制还需根据任务代价和约束求解相应的阻抗动作。

## 2. 总体技术路线与当前范围

R15 的完整技术路线为

$$
\boxed{
\text{柱面 Hertz 与空间观测/MHE}
\rightarrow
\text{未来参数与接触模型}
\rightarrow
\text{NNBO/single-critic HJB}
\rightarrow
\text{ANEOC 预测与认证发布}
\rightarrow
\text{实际阻抗生成}
}.
$$

这条路线分为四个阶段：

1. **环境感知与建模：**基于柱面 Hertz 接触模型，用力、力矩、位姿和运动历史进行 MHE，得到未来参数及接触模型；
2. **最优阻抗控制问题定义：**把未来接触模型接入机器人动力学和约束最优控制，用 NNBO/single critic 表示有限维 HJB–KKT 局部根；
3. **邻近策略更新：**利用根灵敏度预测新参数下的候选，通过认证后直接发布；根精度不足时按局部收敛条件校正，锚点或分支失效时重新认证或求解；
4. **实际阻抗生成：**从已发布候选的 critic 恢复当前反馈，在增益、变化率、执行器和能量约束下形成实际阻抗。

第一阶段先定义接触模型，再用接触观测辨识其中的参数：

$$
\boxed{
\text{柱面 Hertz 接触模型}
\rightarrow
\text{空间观测}
\rightarrow
\text{MHE 在线求解}
\rightarrow
\text{未来参数反演}
\rightarrow
\text{未来接触模型建立}
}.
$$

## 3. 第一阶段：从接触信息到未来环境模型

### 3.1 第一阶段子问题与输入输出

**利用探头接触新区域过程中产生的力、力矩和运动信息，估计环境参数变化，为后续最优阻抗更新提供接触模型。**

采用与 T-RO 接触建模思路一致的低维参数：
$$
\theta_j=(E_j^*,\mu_j,k_{d,j}),\qquad j\in\{0,f\},
\qquad \Delta\theta=\theta_f-\theta_0.
\tag{1}
$$
$E^*$ 为约化有效弹性模量，$\mu$ 为滑动摩擦系数，$k_d$ 为探头当前接触长度下的等效切向黏性系数，单位为 N·s/m。下标 $0,f$ 分别表示当前区域和前方区域。探头半径、长度和位移参考基准由几何标定确定。

输入为同步的六维力/力矩、探头位姿和速度历史。第一阶段输出
$$
\bigl(\hat\theta_{f,k},\Delta\hat\theta_k,\hat b_k,
       \hat z_k,\Sigma_{\zeta,k}\bigr).
\tag{2}
$$
其中 $b$ 为材料分界位置，$z_k$ 为参考表面运动和传感器残余偏置状态，$\Sigma_{\zeta,k}$ 描述联合估计的局部不确定性。它们用于预测下一段运动中的接触响应。

“空间预瞄”来自有限尺寸接触斑的空间覆盖：探头在扫描中逐步接触前方材料，接触响应已包含该区域的信息。估计结果可用后便交给后续阻抗算法。

**选择理由与适用条件。** 本阶段采用准静态法向弹性与滑动摩擦模型，服务于接触变化后的快速估计。[2] 法向松弛、滞回和几何近似的残差需要在目标载荷、速度范围内足够小；辨识得到的 $E^*$ 是该工作范围内的有效参数。该范围的验证方法见附录 G。

### 3.2 接触模型：由球形 Hertz 转为柱面 Hertz

**沿用 T-RO“接触几何与压力分布 → 接触合力 → 力状态变化”的建模路线，将球形接触改为柱面线接触。**

设局部扫描方向为 $x$，圆柱轴向为 $y$，组织对探头的法向力沿 $+z$；三轴构成右手坐标系。圆柱半径为 $R$、有效接触长度为 $L$，接触半宽为 $a$。均匀材料上的约化模量为
$$
\frac{1}{E^*}=
\frac{1-\nu_p^2}{E_p}+
\frac{1-\nu_t^2}{E_t}.
\tag{3}
$$
刚性探头对应第一项可忽略。以 $P=F_n/L$ 表示单位长度法向载荷，柱面 Hertz 关系为
$$
a^2=\frac{4PR}{\pi E^*},
\qquad
p_H(x)=\frac{E^*}{2R}\sqrt{a^2-x^2},
\quad |x|\le a.
\tag{4}
$$
这里 $p_H$ 为压力，单位 Pa；沿圆柱轴向积分后得到线载荷 $q_n=Lp_H$。[3]

**柱面接触还需要明确压入量的参考基准。** 本文固定横向参考距离 $\ell_{\rm ref}>a$，以两侧参考位置的平均高度确定局部法向位移基准。均匀半空间的柱面解给出
$$
\delta=\Psi(a)=\frac{1}{2R}
\left[
\ell_{\rm ref}^2-\ell_{\rm ref}\sqrt{\ell_{\rm ref}^2-a^2}
+a^2\ln\frac{\ell_{\rm ref}+\sqrt{\ell_{\rm ref}^2-a^2}}{a}
\right].
\tag{5}
$$
$\ell_{\rm ref}$ 由几何与位移基准标定后固定；它不是待辨识的材料边界 $b$。式（5）在 $0<a<\ell_{\rm ref}$ 内严格递增，因而给定压入量便可唯一求得 $a=\Psi^{-1}(\delta)$。参考面的法向整体运动记为 $\gamma$，则 $\delta=d-\gamma$，$d$ 为探头相对标定表面的压入位置。推导与参考基准的含义见附录 A。[3]

由式（4）—（5）可得 T-RO 球形力—变形建模链在本文柱面几何下的对应关系。令
$g(x)=\sqrt{a^2-x^2}$、$I_0(a)=\pi a^2/2$，则均匀区域的法向线载荷和持续滑动时的切向线载荷为
$$
\begin{aligned}
q_n(x)&=\frac{LE^*}{2R}g(x),\\
q_t(x)&=-\operatorname{sgn}(v)\mu q_n(x)
        -v k_d\frac{g(x)}{I_0(a)},\qquad |x|\le a .
\end{aligned}
\tag{6a}
$$
沿接触宽度积分，得到本模型直接用于观测和动力学的法向合力、切向合力及绕 $y$ 轴力矩：
$$
\begin{aligned}
F_n&=\int_{-a}^{a}q_n(x)\,dx
    =\frac{\pi LE^*a^2}{4R},\\
F_t&=\int_{-a}^{a}q_t(x)\,dx
    =-\mu F_n\operatorname{sgn}(v)-k_dv,\\
M_y&=\int_{-a}^{a}\frac{x^2}{2R}q_t(x)\,dx
    =\frac{a^2}{8R}F_t .
\end{aligned}
\tag{6b}
$$
其中 $v$ 为探头相对组织的切向速度，负号表示组织对探头的阻力；$M_y$ 以圆柱截面最低点为参考，均匀区域的对称法向压力对该点的法向一阶矩为零。于是，给定 $(d-\gamma,v)$ 和参数 $(E^*,\mu,k_d)$，式（5）、（6a）和（6b）依次确定接触半宽、分布载荷以及完整接触输出。

为与后续最优控制使用的力状态一致，在材料参数和滑动方向不变的均匀区域内，
$$
\begin{aligned}
\dot a&=\frac{\dot d-\dot\gamma}{\Psi'(a)},\qquad
\Psi'(a)=\frac{a}{R}
\ln\frac{\ell_{\rm ref}+\sqrt{\ell_{\rm ref}^2-a^2}}{a},\\
\dot F_n&=k_n(\dot d-\dot\gamma),\qquad
k_n=\frac{\partial F_n}{\partial\delta}
=\frac{\pi LE^*a}{2R\Psi'(a)},\\
\dot F_t&=-\mu\operatorname{sgn}(v)\dot F_n-k_d\dot v,\\
\dot M_y&=\frac{a\dot a}{4R}F_t+\frac{a^2}{8R}\dot F_t .
\end{aligned}
\tag{6c}
$$
式（6a）—（6c）显式完成了“压力/牵引分布—接触合力—力状态变化”的柱面接触模型。它保留 T-RO 的滑动摩擦与速度相关阻力结构；T-RO 中由球形应力积分得到的球帽几何修正不适用于柱面，因而不直接沿用。下一节在同一模型中加入空间变化参数，形成跨区域接触关系。

**选择理由与适用条件。** 该模型用接触半宽、压力分布和少量参数描述法向与切向响应，便于在线反演。[2,3] 适用范围为小变形、局部平面应变、接触宽度相对圆柱半径和长度较小、保持接触且滑动方向在局部窗口内不反转的扫描过程。有限厚度、端部效应及非弹性残差按附录 G 验证。

### 3.3 空间观测：从接触分布得到力和力矩

**把前方材料逐步进入接触斑的过程写入同一个接触模型，用力与力矩的变化估计材料差异及其位置。**

令 $s_i$ 为探头中心在局部组织参考坐标中的扫描位置，局部 $+x$ 取扫描前进方向。在一个覆盖材料过渡的窗口内，取
$$
\theta(s;b)=
\begin{cases}
\theta_0,&s<b,\\
\theta_f,&s\ge b,
\end{cases}
\qquad r_i=b-s_i,\qquad
\mathcal C_i=[s_i-a_i,\ s_i+a_i].
\tag{7}
$$
$\mathcal C_i$ 是第 $i$ 个采样时刻在扫描方向上的接触足迹。探头沿 $+x$ 扫描时，$r_i\ge a_i$ 表示足迹尚未接触前方区域；当 $s_i+a_i=b$，即 $r_i=a_i$ 时，前缘首次到达材料边界；当 $-a_i<r_i<a_i$ 时，两种材料同时参与接触；当 $r_i\le-a_i$ 时，足迹完全进入前方区域。特别地，$0<r_i<a_i$ 时探头中心仍位于边界之前，但 $F_n,F_t,M_y$ 已含前方材料的贡献，这就是本文空间预瞄信息的来源。

沿用柱面 Hertz 的半椭圆形状，本文采用以下**空间参数化近似**：
$$
q_{n,i}(x)=c_HE^*(s_i+x;b)\,g_i(x),
\qquad
c_H=\frac{L}{2R},
\qquad
g_i(x)=\sqrt{a_i^2-x^2},\quad |x|\le a_i.
\tag{8}
$$
均匀材料时，式（8）恢复式（4）；跨区域时，它假定局部载荷强度随有效模量变化，接触半宽仍由式（5）确定。这是本文选定的低维空间接触近似，跨区域误差的含义和独立验证要求见附录 B、G。

定义总形状积分和前方区域的三个积分：
$$
\begin{aligned}
I_0(a)&=\int_{-a}^{a}g(x)\,dx=\frac{\pi a^2}{2},\\
I_f(a,r)&=\int_{-a}^{a}{\bf1}_{x\ge r}g(x)\,dx,\\
I_{1f}(a,r)&=\int_{-a}^{a}x{\bf1}_{x\ge r}g(x)\,dx,\\
I_2(a)&=\int_{-a}^{a}x^2g(x)\,dx=\frac{\pi a^4}{8},\\
I_{2f}(a,r)&=\int_{-a}^{a}x^2{\bf1}_{x\ge r}g(x)\,dx .
\end{aligned}
\tag{9}
$$
在部分覆盖区间 $-a<r<a$ 内，
$$
\begin{aligned}
I_f&=\frac12\left[a^2\arccos(r/a)-r\sqrt{a^2-r^2}\right],\\
I_{1f}&=\frac13(a^2-r^2)^{3/2},\\
I_{2f}&=\frac{a^4}{8}\arccos(r/a)
+\frac{r\sqrt{a^2-r^2}}{8}(a^2-2r^2).
\end{aligned}
\tag{10}
$$
足迹完全在当前区域时 $I_f=I_{1f}=I_{2f}=0$；完全在前方区域时 $I_f=I_0,\ I_{1f}=0,\ I_{2f}=I_2$。

于是法向合力和法向载荷一阶矩为
$$
\begin{aligned}
F_n&=c_H\left[E_0^*I_0+(E_f^*-E_0^*)I_f\right],\\
m_1&=c_H(E_f^*-E_0^*)I_{1f}.
\end{aligned}
\tag{11}
$$
**合力反映接触响应的大小，一阶矩反映载荷向前方或后方的偏移。** 对已知压入量和当前区域模量的理想情况，$\Delta F_n$ 与 $m_1$ 可以分别提供材料变化幅度与边界位置的信息；这一结论的解析证明及退化情形见附录 D。

切向载荷也按同一足迹积分。定义
$$
N_0=c_HE_0^*(I_0-I_f),\quad
N_f=c_HE_f^*I_f,\quad
\alpha=\frac{I_f}{I_0},
$$
则
$$
F_t=-\operatorname{sgn}(v)(\mu_0N_0+\mu_fN_f)
-v\left[(1-\alpha)k_{d,0}+\alpha k_{d,f}\right].
\tag{12}
$$
对应的线牵引为
$$
q_t(x)=-\operatorname{sgn}(v)\mu(s+x;b)q_n(x)
-v\,k_d(s+x;b)\frac{\sqrt{a^2-x^2}}{I_0(a)}.
\tag{13}
$$

力矩必须同时计入法向载荷不对称和切向力的力臂。以探头截面最低点为参考，柱面轮廓在小斜率下为 $z_c(x)=x^2/(2R)$。令 $\sigma=\operatorname{sgn}(v)$，则
$$
\begin{aligned}
M_y
&=\int_{-a}^{a}\left[z_c(x)q_t(x)-xq_n(x)\right]dx\\
&=-m_1-\frac{\sigma c_H}{2R}
\left[\mu_0E_0^*(I_2-I_{2f})+\mu_fE_f^*I_{2f}\right]\\
&\quad-\frac{v}{2RI_0}
\left[k_{d,0}(I_2-I_{2f})+k_{d,f}I_{2f}\right].
\end{aligned}
\tag{14}
$$
接下来把机器人位姿和六维力/力矩转换成上述接触模型的输入。设 $\,{}^BT_{P,i}$ 为测得的探头位姿，$\,{}^PT_C$ 为标定的探头坐标系到接触参考点变换，$\,{}^Bp_o$ 为局部参考面的标定原点，$\mathbf e_x,\mathbf e_z$ 为用基坐标表示的局部扫描和法向单位向量，采样周期为 $h$，则
$$
\begin{aligned}
{}^Bp_{c,i}&=\operatorname{trans}({}^BT_{P,i}{}^PT_C),\\
s_i&=\mathbf e_x^\top({}^Bp_{c,i}-{}^Bp_o),\\
d_i&=-\mathbf e_z^\top({}^Bp_{c,i}-{}^Bp_o),\\
v_i&=\frac{s_i-s_{i-1}}{h}.
\end{aligned}
\tag{14a}
$$
$d_i$ 的正方向取探头压入组织的方向；$\gamma_i$ 使用同一符号约定，因此实际压入量为 $\delta_i=d_i-\gamma_i$。

六维传感器数据先完成重力和惯性补偿。若 $\mathbf r_{SC}$ 为传感器原点指向接触参考点的向量并以传感器坐标表示，$R_{CS}$ 将传感器坐标旋转到局部接触坐标，则
$$
\mathbf F_{C,i}=R_{CS}\mathbf F_{S,i},\qquad
\mathbf M_{C,i}=R_{CS}
\left(\mathbf M_{S,i}-\mathbf r_{SC}\times\mathbf F_{S,i}\right).
\tag{14b}
$$
取局部扫描中的三个建模分量：
$$
y_i=
\begin{bmatrix}F_{C,z,i}&F_{C,x,i}&M_{C,y,i}\end{bmatrix}^{\!\top}
=
\begin{bmatrix}F_{n,i}^{\rm meas}&F_{t,i}^{\rm meas}&M_{y,i}^{\rm meas}\end{bmatrix}^{\!\top},
\qquad
H_H(d_i-\gamma_i,s_i,v_i;\vartheta)=
\begin{bmatrix}F_{n,i}&F_{t,i}&M_{y,i}\end{bmatrix}^{\!\top},
\quad
\vartheta=(\theta_0,\theta_f,b).
\tag{15}
$$
$H_H$ 由式（5）、（8）—（14）直接计算，全部积分都有闭式表达。接触斑与压力分布是模型输出；传感器直接测得的是合力和合力矩。

**选择理由与适用条件。** 保留力矩使材料变化的“幅度”与“空间位置”有机会被区分，采用同一压力模型保证观测与后续预测的一致性。局部扫描方向和材料参考坐标需明确；显著的未观测切向组织运动会影响边界定位。式（14）避免把全部测得力矩都解释为材料变化。

### 3.4 MHE 在线求解：由观测历史估计环境变化

**使用 MHE 在最近一段接触数据中联合估计两侧材料参数、边界位置和必要的表面运动状态。**

设采样周期为 $h$，窗口 $\mathcal I_k=\{j_k,\ldots,k\}$，$j_k=k-N+1$。运行时，当前区域参数 $\theta_0$ 由进入过渡区之前的标定或上一估计窗口提供；MHE 用该结果作为带协方差的到达先验，并在其不确定范围内修正 $\theta_0$，主要利用新进入窗口的数据估计 $\theta_f$ 和 $b$。以
$$
z_i=(\gamma_i,\dot\gamma_i,\beta_i^\top)^\top,\qquad \beta_i\in\mathbb R^3
$$
表示法向参考表面运动与三个观测分量的残余偏置。短窗口状态模型为
$$
z_{i+1}=Az_i+w_i,\qquad
A=\operatorname{blkdiag}\!\left(
\begin{bmatrix}1&h\\0&1\end{bmatrix},I_3\right),
\tag{16}
$$
观测模型为
$$
y_i=H_H(d_i-\gamma_i,s_i,v_i;\vartheta)+\beta_i+e_i.
\tag{17}
$$
几何和位姿标定误差计入观测误差预算；$R_y$ 为用于加权的有效观测协方差，$Q$ 为状态推进协方差。已知的系统性模型偏差不能仅靠增大 $R_y$ 消除，其影响在附录 F 单独分析。

联合变量 $\mathcal Z_k=(\vartheta,z_{j_k:k})$ 的 MHE 问题为
$$
\begin{aligned}
\mathcal Z_k^*&\in\arg\min_{\mathcal Z_k}J_k(\mathcal Z_k),\\
J_k={}&
\|z_{j_k}-\bar z_{j_k}\|_{P^{-1}}^2
+\|\vartheta-\bar\vartheta\|_{\Lambda}^2\\
&+\sum_{i=j_k}^{k}
\|y_i-H_H(d_i-\gamma_i,s_i,v_i;\vartheta)-\beta_i\|_{R_y^{-1}}^2\\
&+\sum_{i=j_k}^{k-1}\|z_{i+1}-Az_i\|_{Q^{-1}}^2,\\
\text{s.t. }&
E_{\min}\le E_0^*,E_f^*\le E_{\max},\\
&0\le\mu_0,\mu_f\le\mu_{\max},\quad
0\le k_{d,0},k_{d,f}\le k_{d,\max},\\
&b\in\mathcal B_k,\quad z_i\in\mathcal Z_{\rm adm},\quad
\delta_{\min}\le d_i-\gamma_i\le\delta_{\max}.
\end{aligned}
\tag{18}
$$
$P,Q,R_y$ 为正定矩阵，$\Lambda$ 为先验精度矩阵；其中 $\bar\vartheta$ 的 $\theta_0$ 分量取上一窗口均值 $\bar\theta_{0,k}$，$\Lambda$ 的对应块取其协方差逆 $\Sigma_{\theta_0,k}^{-1}$。$\theta_f$ 和 $b$ 使用较宽但有物理边界的初始先验，使新数据能够更新它们。$\mathcal B_k$ 是窗口覆盖和历史信息给出的有限边界位置区间，$\mathcal Z_{\rm adm}$ 是标定允许的紧状态域。参数边界、压入范围和初始表面状态来自标定及历史数据，并取 $\delta_{\min}=\Psi(a_{\min})$、$\delta_{\max}=\Psi(a_{\max})$，其中 $0<a_{\min}<a_{\max}<\ell_{\rm ref}$。窗口共有 $N$ 个观测与 $N-1$ 个状态推进残差。

**在线用窗口移位热启动的 SQP 求解式（18）。** 对每个候选参数和表面状态，先用带区间保护的 Newton 法解一维方程 $\Psi(a_i)=d_i-\gamma_i$，再用闭式积分得到力、力矩与导数。SQP 根据整段观测残差修正参数和状态，接受步长后重新计算观测模型。局部导数、初值和停止准则见附录 C、E。

得到满足数值准则的局部解 $\widehat{\mathcal Z}_k$ 后，直接提取
$$
\hat\theta_{f,k}=\Pi_f\widehat{\mathcal Z}_k,\qquad
\Delta\hat\theta_k=(\Pi_f-\Pi_0)\widehat{\mathcal Z}_k,\qquad
(\hat b_k,\hat z_k)=\Pi_{b,z_k}\widehat{\mathcal Z}_k.
\tag{19}
$$
这就是前方参数反演；不需要再解一个独立反演问题。前方区域在本窗口中已经产生接触信息，其估计值用于后续运动的模型预测。

**选择理由与适用条件。** MHE 将多步观测、表面运动和物理参数边界放入同一优化问题。[4] 估计窗口需包含有信息量的空间覆盖变化；同时估计 $\mu$ 与 $k_d$ 还要求载荷与速度的回归量不共线。所需样本数由这些信息、噪声与计算时间共同决定。相关秩条件和正则化的作用见附录 D、E。

### 3.5 估计结果的可信度与未来模型建立

**用同一次联合估计的参数和表面状态实例化接触模型，并给出与该局部解对应的不确定性。**

令 $\zeta=(\vartheta,z_k)$，以 $\Sigma_{\zeta,k}$ 表示由加权残差 Jacobian 得到的局部协方差近似。参数变化的局部协方差为
$$
\Sigma_{\Delta\theta,k}
=D_\Delta\Sigma_{\zeta,k}D_\Delta^\top,
\qquad
D_\Delta=\frac{\partial(\theta_f-\theta_0)}{\partial\zeta}.
\tag{20}
$$
它保留两侧参数与表面状态的相关性。该指标依赖当前局部解、噪声模型与辨识信息；计算方式及模型偏差的影响见附录 F。

给定后续候选运动 $(d,s,v)$，接触预测直接写成
$$
\widehat{\mathcal M}_{f,k}:
\qquad
\widehat{\mathcal W}_c
=H_H(d-\hat\gamma,s,v;
       \hat\theta_0,\hat\theta_f,\hat b),
\qquad
(\hat\gamma,\widehat{\dot\gamma})\text{按式（16）推进}.
\tag{21}
$$
这里 $\widehat{\mathcal W}_c$ 表示建模的三个力/矩分量，偏置 $\beta$ 不加入物理接触力。足迹完全进入前方区域后，式（21）自然变为参数 $\hat\theta_f$ 下的均匀柱面模型。

**选择理由与适用条件。** 名义参数用于计算，局部不确定性用于说明当前估计的精度；两者来自同一次 MHE。预测路径须位于本次两区域模型的空间及载荷范围内。局部协方差不等于全局包含保证，未辨识充分的参数不能被当作精确真值交接。

### 3.6 接触动力学与 NNBO/ANEOC 的接口

**将柱面接触模型写成后续最优控制可使用的力与力变化关系，连接环境辨识和最优阻抗更新。**

在固定材料场、固定滑动方向和同一光滑接触区间内，由 $F_n=F_n(\delta,s;\vartheta)$ 得
$$
\dot F_n=
\frac{\partial F_n}{\partial\delta}(\dot d-\dot\gamma)
+\frac{\partial F_n}{\partial s}v.
\tag{22}
$$
第一项表示压入变化，第二项表示扫描经过不同材料时的接触响应变化。对三个力/矩分量统一写成
$$
\dot{\mathcal W}_c=
H_{H,\delta}(\dot d-\dot\gamma)
+H_{H,s}v
+H_{H,v}\dot v.
\tag{23}
$$
这保留了 T-RO 的力状态动力学构造，并显式加入空间参数变化的影响。导数来自当前柱面模型，具体表达和分支条件见附录 C。

在同一局部接触参考点，将三个分量映射成六维力旋量
$$
\mathcal W_c^{(6)}=
(F_t,0,F_n,0,M_y,0)^\top,
$$
再经坐标变换与机器人 Jacobian 接入动力学：
$$
M_r(q_r)\ddot q_r+C_r(q_r,\dot q_r)\dot q_r+g_r(q_r)
=\tau_r+J_c^\top\mathcal W_c^{(6)}.
\tag{24}
$$
后续 NNBO 在该动力学及任务代价下定义有限维最优控制问题；ANEOC 使用环境描述符的变化和相应根灵敏度更新最优解。接触力若作为动态状态加入，初值需与式（21）一致。

**选择理由与适用条件。** 接触参数及其导数共同形成后续算法的模型输入，避免估计与优化使用不同接触关系。Hertz 空间积分在足迹跨过材料边界的起止点仅具有有限光滑性，后续灵敏度分析须在适用的局部分支内进行；具体限制见附录 C。此处完成接口定义，后续最优控制和根更新在下一阶段展开。

### 3.7 第一阶段形成的完整路线

$$
\boxed{
\begin{gathered}
\text{柱面 Hertz 参数化接触模型}\\
\downarrow\\
\text{接触过程中的力、力矩、位姿与速度}\\
\downarrow\\
\text{MHE 联合估计 }(\theta_0,\theta_f,b,z_k)\\
\downarrow\\
\text{提取 }\hat\theta_f,\Delta\hat\theta
\text{ 并实例化未来接触模型}\\
\downarrow\\
\text{提供给后续 NNBO 与 ANEOC}
\end{gathered}}
\tag{25}
$$

本阶段的核心结果是同一接触模型下的“观测—参数—预测”连接。其可用精度由空间辨识信息和模型残差决定，计算速度由实际 MHE 实现验证；附录 D—G 给出相应依据与验证标准。

## 4. 第二阶段：由未来接触模型求解最优阻抗

### 4.1 从第一阶段输出构造控制问题

**把未来环境参数和扫描任务写成固定维度描述符，使每次更新对应同一种最优控制问题。**

本阶段输入为式（19）—（21）的估计及未来接触模型、当前机器人状态和扫描参考，输出为最优阻抗反馈的有限维表示。每次模型更新以时间戳 $t_k$ 为起点，用 $t\in[0,T]$ 表示相对时间，$T$ 是固定控制时域。在这个时域内，采用第一阶段的常参数材料场与常速度表面预测：
$$
\begin{aligned}
\eta_k&=(\hat\theta_{0,k},\hat\theta_{f,k},\hat b_k,
                 \hat\gamma_k,\widehat{\dot\gamma}_k,\rho_k),\\
\gamma(t;\eta_k)&=\hat\gamma_k+t\widehat{\dot\gamma}_k,\\
\bigl(\xi_d(t),F_d(t)\bigr)&=\mathcal R(t;\rho_k),\\
\mathcal R(t;\rho_k)&=\sum_{j=0}^{p_r}\rho_{k,j}(t/T)^j .
\end{aligned}
\tag{26}
$$
$\xi_d\in\mathbb R^6$ 是局部位姿参考，$F_d>0$ 是目标法向接触力，$\rho_k$ 为二者在固定阶数 $p_r$ 下的系数。控制时域、代价权重、机器人参数和限制值在一次实验中固定。机器人状态是价值函数的自变量，不重复装入 $\eta$。偏置 $\beta$ 用于观测校正；$\Sigma_\zeta$ 用于附录 L 的局部误差传播，二者均不被当作真实接触载荷。

参数改变和时间窗口移位都通过式（26）更新。最近认证锚点的描述符记为 $\eta_a$，预测使用 $\eta_k-\eta_a$；锚点随其根误差和灵敏度证书一起保存。由 $\Delta\hat\theta$ 与 $\hat\theta_f$ 也可恢复 $\hat\theta_0$，计算中直接保留同次 MHE 的两个相关估计。

**选择理由与适用条件。** 固定有限时域直接表示扫描与表面运动，不需要把持续前进的位置当作无限时域平衡点。$T$ 应覆盖对当前动作有影响的前方接触变化，且处于式（21）的预测有效范围；它不构成中心越界前必须完成的截止时间。以下优化采用名义估计，局部协方差不提供鲁棒全局保证。

### 4.2 机器人—接触—实际阻抗状态空间

**将柱面接触力作为机器人动力学的输出反馈，并把实际阻抗及其变化率分别定义为状态和动作。**

本模块输入 $\eta$、关节位置 $q_r$、速度 $w_r=\dot q_r$ 和实际增益，输出闭合状态方程。令 $\xi(q_r)$ 为接触参考点的局部位姿坐标：前三项沿第 3 节固定的 $x,y,z$ 轴，后三项取无奇异的局部姿态坐标。定义
$$
\begin{aligned}
e&=\xi(q_r)-\xi_d(t),\\
J_\xi(q_r)&=\partial\xi/\partial q_r,
\qquad \dot e=J_\xi w_r-\dot\xi_d(t),\\
\delta&=d(q_r)-\gamma(t;\eta),\\
s&=s(q_r),\qquad v=J_s(q_r)w_r,
\quad J_s=\partial s/\partial q_r,\\
\mathcal W_c&=H_H(\delta,s,v;\hat\theta_0,\hat\theta_f,\hat b),\\
\mathcal W_c^{(6)}&=(F_t,0,F_n,0,M_y,0)^\top .
\end{aligned}
\tag{27}
$$
$s(q_r),d(q_r)$ 按式（14a）计算；因 $\xi_3$ 沿 $+z$ 而 $d$ 沿压入方向，二者符号相反。几何半宽每次仍由 $\Psi(a)=\delta$ 求得。接触力矩通过物理接触 Jacobian $J_c$ 进入关节动力学；$J_\xi$ 则与位姿坐标的广义作用力功率共轭，两种 Jacobian 不混用。

优化扫描切向和法向的刚度、阻尼，其他通道取固定正值：
$$
\begin{aligned}
g_a&=(\kappa_t,\kappa_n,d_t,d_n)^\top,\qquad \dot g_a=u_g,\\
K(g_a)&=\operatorname{diag}(\kappa_t,\kappa_y,\kappa_n,
                                  \kappa_{r1},\kappa_{r2},\kappa_{r3}),\\
D(g_a)&=\operatorname{diag}(d_t,d_y,d_n,d_{r1},d_{r2},d_{r3}),\\
f_d(t)&=(0,0,-F_d(t),0,0,0)^\top,\\
\tau_{\rm imp}&=g_r(q_r)-B_rw_r
      +J_\xi^\top\bigl[f_d-K(g_a)e-D(g_a)\dot e\bigr],\\
a_r(t,q_r,w_r,g_a;\eta)&=M_r^{-1}
 \bigl[\tau_{\rm imp}+J_c^\top\mathcal W_c^{(6)}-C_rw_r-g_r\bigr].
\end{aligned}
\tag{28}
$$
$B_r\succ0$ 为固定关节阻尼；$f_d$ 提供方向与组织反力相反的法向前馈。$\kappa_t,\kappa_n$ 是控制器刚度，与环境模量 $E^*$、接触斜率 $k_n$ 不同；$d_t,d_n$ 也不同于组织切向黏性系数 $k_d$。式（28）规定本文唯一的力矩执行接口，所优化的是该控制器族中的阻抗。

时变刚度可向机械系统注入功率。用一个能量状态 $E_T$ 记录该交换，令
$$
\begin{aligned}
p_g(e)&=\tfrac12(e_1^2,e_3^2,0,0)^\top,\qquad
P_K=p_g^\top u_g,\\
P_D&=\dot e^\top D(g_a)\dot e+w_r^\top B_rw_r,\\
\dot E_T&=\chi P_D-p_g^\top u_g,\qquad 0\le\chi<1 .
\end{aligned}
\tag{29}
$$
$E_T$ 为控制器能量账本，初值给定且不因环境更新而重置。刚度增加所需的 $P_K$ 从账本扣除，比例 $\chi$ 的阻尼耗散可回收；$E_T\ge E_{T,\min}>0$ 由下一节动作约束维持。该能量核算与变量阻抗的无源设计思想一致，[1] 针对式（28）的证明见附录 H。

于是 $n_r$ 关节机器人的控制状态与动作通道明确为
$$
\boxed{
\begin{aligned}
x&=(q_r^\top,w_r^\top,g_a^\top,E_T)^\top\in\mathbb R^{2n_r+5},\\
\dot x&=f_0(t,x;\eta)+G_g(t,x;\eta)u_g,\\
f_0&=\begin{bmatrix}w_r\\a_r\\0_4\\\chi P_D\end{bmatrix},
\qquad
G_g=\begin{bmatrix}0_{n_r\times4}\\0_{n_r\times4}\\I_4\\-p_g^\top\end{bmatrix}.
\end{aligned}}
\tag{30}
$$
给定 $(t,x,u_g,\eta)$，式（27）—（30）逐项确定所有状态导数。这里表面预测和参考已作为显式时间函数给出，无需再复制为状态。接触输出也不重复作为独立状态；其导数仍严格满足式（23），其中 $\dot v=J_s a_r+\dot J_s w_r$。增益率先改变增益，再改变加速度、运动和接触响应。

**选择理由与适用条件。** 该增广系统保持第一阶段接触关系，且 action 与实际增益有直接执行对应。它要求力矩接口可实现、局部姿态图有效、接触几何和滑动方向满足第 3 节条件。原 NNBO 的 Hunt–Crossley 模型能反解为“力输入—速度输出”，[5] 本文准静态法向模型则依靠机器人惯性闭合动力学；不添加未经辨识的法向黏性来仿写原文式（23）。

### 4.3 任务代价、约束与有限时域 HJB

**用接触任务的累计代价定义最优增益变化策略，并在同一个问题中约束增益、变化率和能量。**

本模块输入式（30）和任务限制，输出价值函数问题及 Hamiltonian。定义允许动作
$$
\begin{aligned}
\mathcal U(t,x;\eta)=\{u_g:\;&-\bar r_g\le u_g\le\bar r_g,\\
&-\alpha_g(g_a-g_{\min})\le u_g
                       \le\alpha_g(g_{\max}-g_a),\\
&p_g^\top u_g\le\chi P_D+\alpha_E(E_T-E_{T,\min})\}\\
&=\{u_g:A_u(t,x;\eta)u_g\le b_u(t,x;\eta)\},
\end{aligned}
\tag{31}
$$
不等式逐分量理解，$\alpha_g,\alpha_E>0$。增益下界严格为正；在初始增益和能量可行时，连续执行式（31）维持这些边界。这个集合包含 $u_g=0$，因此在该域内非空，但保持增益不变本身不保证接触任务可行。

接触与执行器限制是整条轨迹的限制。以 $h_j(t,x;\eta)\ge0$ 统一表示
$$
\begin{gathered}
q_{\min}\le q_r\le q_{\max},\qquad |w_r|\le\bar w_r,
\qquad |\tau_{\rm imp}|\le\bar\tau_r,\\
\delta_{\min}\le\delta\le\delta_{\max},\qquad
F_{\min}\le F_n\le F_{\max},\qquad |F_t|\le F_{t,\max},\\
0<v_{\min}\le v\le v_{\max},
\qquad x\text{ 位于标定的局部位姿与接触域内}.
\end{gathered}
\tag{32}
$$
持续扫描允许滑动；$F_{t,\max}$ 限制切向载荷，不把静摩擦锥条件当作滑动摩擦约束。式（32）随预测轨迹检查，不能由式（31）单独推出。

用固定尺度归一化各误差后，运行代价和终端代价取
$$
\begin{aligned}
\ell(t,x,u_g;\eta)&=\ell_0(t,x;\eta)+u_g^\top R_gu_g,\\
\ell_0&=e^\top Q_e e+\dot e^\top Q_v\dot e
       +q_F(F_n-F_d)^2
       +(g_a-g_{\rm ref})^\top Q_g(g_a-g_{\rm ref}),\\
\Phi(x;\eta)&=e(T)^\top Q_{eT}e(T)
             +q_{FT}\bigl[F_n(T,x;\eta)-F_d(T)\bigr]^2 .
\end{aligned}
\tag{33}
$$
$R_g,Q_g\succ0$，其余权重非负；$g_{\rm ref}$ 为预先给定的柔顺偏好，位于增益范围内。阻抗偏好与跟踪、力误差竞争，使代价明确表达“怎样的阻抗更合适”。$Q_g\succ0$ 不保证整个非线性控制问题具有唯一最优解。

价值函数为
$$
V^*(t,x;\eta)=\inf_{u_g(\cdot)}
\left\{\Phi(x(T);\eta)+\int_t^T\ell(\tau,x(\tau),u_g(\tau);\eta)\,d\tau\right\},
\tag{34}
$$
其中 $x(t)=x$，满足式（30）—（32）。在存在可行续接轨迹的域内、且价值函数可微的局部区域，定义
$$
\begin{aligned}
\mathcal H(t,x,u_g,p;\eta)&=\ell_0+u_g^\top R_gu_g
                                    +p^\top(f_0+G_gu_g),\\
0&=\partial_tV^*+
\min_{u_g\in\mathcal U(t,x;\eta)}\mathcal H(t,x,u_g,\nabla_xV^*;\eta),\\
V^*(T,x;\eta)&=\Phi(x;\eta).
\end{aligned}
\tag{35}
$$
本文的局部 HJB 近似与正则根分析针对式（32）保持严格裕量的轨迹邻域；增益、速率和能量约束可以激活。若物理状态约束成为最优轨迹的活动边界，则需要相应的状态约束边界解，不能用内部式（35）的有限配置点根声称已解决该情形。轨迹检查和这一边界的含义见附录 I、L。

环境参数已知后仍需解式（35），因为同一种组织在不同运动状态、目标力、增益和能量条件下可能要求不同动作。接触模型提供 $f_0$ 和代价中的 $F_n$；HJB 决定未来累计代价最小的反馈。

**选择理由与适用条件。** 二次动作代价与仿射动作通道使每个 Hamiltonian 极小化成为严格凸的四维问题。有限时域保证非零目标力和持续扫描的累计代价有明确含义；这里只声称有限时域名义最优性，不从终端代价直接推出滚动闭环渐近稳定或递归可行。

### 4.4 single critic、策略恢复与有限维根

**用一个固定基的 critic 表示价值函数，通过 Hamiltonian 恢复动作，并把待求量组成明确的 HJB–KKT 方程。**

本模块输入式（35），输出根 $Z^*(\eta)$ 及其 critic 分量。迁移 NNBO 的“价值近似—价值梯度—动作恢复”机制，[5] 在固定归一化坐标 $\bar x=S_x^{-1}(x-x_c)$、$\bar t=t/T$ 上取线性无关光滑基 $\phi(\bar t,\bar x)\in\mathbb R^{n_W}$，写成
$$
\begin{aligned}
\psi(t,x)&=(1-t/T)\phi(t/T,S_x^{-1}(x-x_c)),\\
V_N(t,x;W,\eta)&=\Phi(x;\eta)+\psi(t,x)^\top W,\\
\partial_tV_N&=\psi_t^\top W,\qquad
\nabla_xV_N=\nabla_x\Phi+J_\psi^\top W,
\quad J_\psi=\partial\psi/\partial x .
\end{aligned}
\tag{36}
$$
该形式精确满足终端条件。基、尺度、中心与节点在根更新期间固定，只有输出权重 $W$ 待求；不在线改变隐藏层参数或状态维度。附录 I 给出可直接构造的多项式基与节点规则。

critic 产生的反馈为
$$
\begin{aligned}
\pi_N(t,x;W,\eta)
&=\arg\min_{u_g\in\mathcal U(t,x;\eta)}
  \left[u_g^\top R_gu_g+
  (G_g^\top\nabla_xV_N)^\top u_g\right],\\
u_{\rm free}&=-\tfrac12R_g^{-1}G_g^\top\nabla_xV_N\quad
                       \text{（当该动作可行时）}.
\end{aligned}
\tag{37}
$$
动作唯一性来自 $R_g\succ0$，而非 critic 权重唯一性。由于式（30）含能量通道，$G_g^\top\nabla_xV_N=\nabla_{g_a}V_N-p_g\partial_{E_T}V_N$；忽略后一项会求解另一个问题。受约束时求解式（37）的同一极小化，不在求解后另行剪裁动作。

在 $n_W$ 个固定节点 $(t_i,x_i)$、$t_i<T$ 处，令 $u_i\in\mathbb R^4$、$\lambda_i\in\mathbb R^{m_u}$，其中式（31）共有 $m_u=17$ 行。记
$$
\begin{aligned}
p_i&=\nabla_xV_N(t_i,x_i),\\
\mathcal B_i&=\partial_tV_N(t_i,x_i)+\ell_{0,i}
                  +u_i^\top R_gu_i+p_i^\top(f_{0,i}+G_{g,i}u_i),\\
\mathcal S_i&=2R_gu_i+G_{g,i}^\top p_i+A_{u,i}^\top\lambda_i,\\
s_i^u&=b_{u,i}-A_{u,i}u_i .
\end{aligned}
\tag{38}
$$
$s_i^u$ 是动作不等式松弛量，与扫描位置 $s_i$ 不同。用互补关系构造固定长度根：
$$
\boxed{
\begin{aligned}
Z&=(W,\{u_i,\lambda_i\}_{i=1}^{n_W}),\\
F(Z,\eta)&=\operatorname{col}_{i=1}^{n_W}
 \bigl(\mathcal B_i,\mathcal S_i,\lambda_i\odot s_i^u\bigr)=0,\\
\lambda_i&\ge0,\qquad s_i^u\ge0,\qquad
W^*(\eta)=\Pi_WZ^*(\eta).
\end{aligned}}
\tag{39}
$$
这里方程数与未知数均为 $n_W(1+4+m_u)$；非负性必须同时检查，单有互补乘积为零不够。$Z^*$ 是这个有限维方程的接受局部根，节点动作是对应 critic 的 Hamiltonian 极小值；不将它称为精确无限维 HJB 解，也不声称 $W$ 是某个未定义优化目标的全局最优权重。

首次根用阻尼 Newton 求解：
$$
F_Z(Z^j,\eta_0)\Delta Z^j=-F(Z^j,\eta_0),\qquad
Z^{j+1}=Z^j+\alpha_j\Delta Z^j,
\quad 0<\alpha_j\le1 .
\tag{40}
$$
以 $\frac12\|F\|^2$ 回溯，求得数值锚点后认证其根误差、局部正则域及域内逼近质量。首次求解允许离线完成；初值、导数块及局部收敛条件见附录 I—K，局部域的可计算充分条件见认证证明第 1 项。后续环境更新使用该锚点及灵敏度产生候选。

**选择理由与适用条件。** single critic 省去单独的 actor，动作由明确的四维凸极小化恢复。有限配置点使未知量有限，但基阶数和节点数会随状态维度增长；局部根存在、$F_Z$ 非奇异和域内近似质量都需检查，不能由“神经网络可逼近”直接推出。正文选定一种配置点根，不并列引入另一种学习算法。

## 5. 第三阶段：ANEOC 预测与认证发布

### 5.1 从同一根方程计算参数灵敏度

**对已接受的 HJB–KKT 根隐式求导，计算环境描述符变化时根应当如何变化。**

本模块输入认证锚点 $(\widehat Z_a,\eta_a)$ 与式（39）的导数，输出数值灵敏度 $\widehat C_a$ 及其误差界。NEOC 用参数相关的价值方程推导邻近反馈变化；[6] 本文将这一机制用于固定的有限维 HJB–KKT 根，约束根微分采用隐式函数关系。[7] 精确根分支的灵敏度满足
$$
\begin{aligned}
F(Z^*(\eta),\eta)&=0,\\
F_Z\frac{dZ^*}{d\eta}+F_\eta&=0,\\
F_ZC_\eta&=-F_\eta,
\qquad C_\eta=\frac{dZ^*}{d\eta},\qquad
\frac{dW^*}{d\eta}=\Pi_WC_\eta .
\end{aligned}
\tag{41}
$$
在数值锚点处分解 $F_Z$，以多个右端项求得 $\widehat C_a$，用锚点残差、导数误差和线性求解残差认证 $\rho_a\ge\|\widehat Z_a-Z^*(\eta_a)\|$、$\varepsilon_{C,a}\ge\|\widehat C_a-D_\eta Z^*(\eta_a)\|$。$F_\eta$ 包括接触响应、参考与终端代价的参数导数；具体导数和界见附录 J 及认证证明第 1—2 项。

**选择理由与适用条件。** 灵敏度复用已求出的最优解。需要活动动作约束满足严格互补与独立性，且整个 $F_Z$ 非奇异；这些条件不同于单个 Hamiltonian 的严格凸性。二阶预测误差还要求所用根方程在参数邻域内二阶光滑。接触足迹边缘或活动集切换处不套用该结论；保留根的维度不等于自动保留其正则性。

### 5.2 候选预测、发布条件与按需校正

**从最近认证锚点预测候选，满足可行性和性能预算时直接发布。**

输入为新描述符 $\eta$ 和最近认证锚点，输出为候选 critic 及其发布证书。令
$$
\Delta\eta=\eta-\eta_a,\qquad
Z_{\rm pred}=\widehat Z_a+\widehat C_a\Delta\eta .
\tag{42}
$$
在认证参数域 $\mathcal P_{\rm cert}$ 内，若 $\|D^2Z^*\|\le M_Z$，计入预测运算误差 $\varepsilon_p$ 后有
$$
\|Z_{\rm pred}-Z^*(\eta)\|
\le\rho_a+\varepsilon_{C,a}\|\Delta\eta\|
 +\tfrac12M_Z\|\Delta\eta\|^2+\varepsilon_p
 =:\rho_Z .
\tag{43}
$$
由式（37）在当前查询状态恢复 $\pi_p=\pi_N(t,x;\Pi_WZ_{\rm pred},\eta)$。附录 L 将 $\rho_Z$ 传到策略和轨迹误差；候选 rollout 的积分缺陷及轨迹管给出约束下界，域内 Bellman 残差给出局部性能界。以下均在固定尺度下计算。

若候选尚未满足发布条件，而根分支有效且预测误差处于附录 K 的认证收敛域，可从 $Z^0=Z_{\rm pred}$ 进行
$$
\begin{aligned}
F_Z(Z^j,\eta)\Delta Z^j&=-F(Z^j,\eta),\\
Z^{j+1}&=Z^j+\alpha_j\Delta Z^j .
\end{aligned}
\tag{44}
$$
设 $\underline h_{j,l}^{\rm cand}$ 是候选在第 $l$ 个验证区间的连续约束下界，$\varepsilon_{\rm pred}$ 为候选与有限根策略的代价差界，$\varepsilon_{\rm conf}$ 为有限根策略与局部连续最优值的差界。发布准则为
$$
\begin{aligned}
\eta&\in\mathcal P_{\rm cert},\qquad \rho_Z\le\varepsilon_{Z,\rm allow},\\
\underline h_{j,l}^{\rm cand}&\ge m_j\quad\text{对所有 }j,l,\\
\varepsilon_J:=\varepsilon_{\rm pred}+\varepsilon_{\rm conf}
 &\le\varepsilon_{J,\rm allow}.
\end{aligned}
\tag{45}
$$
式（45）同时要求动作恢复可行、闭环适定并留在认证轨迹管内，且候选与有限根参考策略到终点 $T$ 的可行性及整个比较域的残差覆盖已认证。此时相对于明确的局部比较值 $V_{\mathcal D}^*$，有 $0\le J(\pi_p)-V_{\mathcal D}^*\le\varepsilon_J$。完整条件见附录 I、L 和认证证明第 3—6 项。预测候选满足这些条件即直接发布；校正后按同一条件重新检查。

发布记录与锚点记录分别保存。根误差过大且局部校正条件成立时使用 correct；灵敏度老化或原域不能覆盖后续查询时，建立通过认证的新锚点并重算灵敏度；配置或分支不再受原证书覆盖时重新求解并认证。校正只保证缩小根误差，性能残差或物理可行性不足须按其实际原因处理。触发条件与多步界见附录 K 和认证证明第 7 项。

**选择理由与适用条件。** 从认证锚点直接预测使小参数变化后的策略发布可复用已有灵敏度。该结论限于证书覆盖的名义模型和比较域；辨识、模型及执行误差按有依据的误差预算纳入。

### 5.3 计算收益的判定

**以相同根、相同精度的求解比较，检验预测器是否减少环境更新后的计算时间。**

输入为固定描述符序列与式（39），比较 ANEOC、复用旧 Jacobian 分解的 cached-J chord，以及新参数下的一次 Newton；完整 cold/warm solve 可作为补充。各方法使用同等可行性、性能精度和发布证书，均允许提前准备及结构化求解。串行执行时总延迟记为
$$
\tau_{\rm total}=\tau_{\rm data}+\tau_{\rm MHE}
 +\tau_{\eta}+\tau_{\rm pred}+\tau_{\rm corr}
 +\tau_{\rm check}+\tau_{\rm comm}+\tau_{\rm apply} .
\tag{46}
$$
预测直接发布时 $\tau_{\rm corr}=0$；并行实现按实际依赖图计算关键路径。另记录准备、维护、等待和失败恢复成本。设算法 $i$ 的额外准备成本为 $P_i$、每查询含认证和恢复的服务成本为 $S_{i,l}$，则相对基线 $B$ 的摊销优势条件为 $\sum_{l=1}^{Q}(S_{B,l}-S_{A,l})>P_A-P_B$。事件延迟还须计入参数到达时未完成的准备工作。证明与比较边界见认证证明第 8 项；附录 M 保留已有机制核对。

**选择理由与适用条件。** 三种候选均可能具有二阶根误差。优势取决于空间预瞄提供的实际准备时间、参数维数、有效查询次数和认证成本；具体 R15 实现须核验相应成本条件。

## 6. 第四阶段：实际阻抗生成与在线执行

### 6.1 从已发布候选恢复当前反馈

**在当前机器人状态重新计算 Hamiltonian 动作，把它积分为实际施加的刚度和阻尼。**

输入为最近已发布候选 $Z_{\rm pub}$、其描述符及参考时间原点 $t_k$，输出为每个控制周期的阻抗和关节力矩。当前时间为 $t_{\rm now}$，用实际相位 $t_a=t_{\rm now}-t_k<T$ 计算
$$
\begin{aligned}
\widehat W&=\Pi_W Z_{\rm pub},\\
u_g(t_{\rm now})&=\pi_N(t_a,x_{\rm meas};\widehat W,\eta_k),\\
g_a(t_{\rm now}+h_c)&=g_a(t_{\rm now})+
                      \int_{t_{\rm now}}^{t_{\rm now}+h_c}u_g(\tau)\,d\tau,\\
(K_a,D_a)&=(K(g_a),D(g_a)),\\
\tau_r&=g_r-B_rw_r+J_\xi^\top[f_d-K_a e-D_a\dot e] .
\end{aligned}
\tag{47}
$$
$h_c$ 是低层更新周期。节点动作用于配置方程；实际执行在当前状态求解式（37）。发布时保留实际 $g_a,E_T$，用最新状态和实际相位对齐证书。

式（47）按名义连续反馈解释。发布前，从最新状态 rollout 至 $T$，用候选轨迹管认证式（32），并按附录 I、L 完成域内性能检查。数字实现将动作求解、采样保持、增益积分和力矩误差计入相应预算。旧策略仅在其证书仍有效时继续执行；无法及时得到有效候选时，使用已有独立 supervisor 的经验证接口，其进入条件不满足则当前方法不授权继续执行。

对名义力矩实现，定义储能
$$
\begin{aligned}
\mathcal S&=\tfrac12 w_r^\top M_rw_r+\tfrac12e^\top K(g_a)e+E_T,\\
\dot{\mathcal S}&=(J_cw_r)^\top\mathcal W_c^{(6)}+P_{\rm ref}
                      -(1-\chi)P_D,\\
P_{\rm ref}&=(J_\xi w_r)^\top f_d
                -\dot\xi_d^\top(K e+D\dot e).
\end{aligned}
\tag{48}
$$
能量账本抵消变刚度功率，移动参考和前馈功率显式保留。该关系证明声明端口下的能量核算；它既不把运动参考当作无功输入，也不单独证明移动组织上的渐近稳定。力矩误差与参考切换的处理见附录 H、L。

**选择理由与适用条件。** 从 critic 到增益、力矩的链路与 HJB 内的动力学完全相同。实际平台须实现这套力矩与增益接口，且闭环轨迹留在验证域内；如果平台使用另一种增益响应动态，需先辨识并替换式（30），不能继续沿用当前根宣称同一最优性。

### 6.2 完整在线算法

**持续更新环境模型，预测并认证候选，按需恢复根精度，在低层周期施加已发布阻抗策略。**

算法输入为同步传感器历史、机器人模型、扫描参考和固定数值设置；输出为 $g_a(t)$ 与 $\tau_r(t)$。首次根、验证域、容差与能量初值在运行前建立；首次可用环境描述符可取均匀模型 $\theta_f=\theta_0$，保持全部参数槽位与根维度。

```text
算法 1  空间接触估计驱动的 ANEOC 最优阻抗
准备：固定 T、参考基、critic 基、配置节点、归一化尺度和约束行顺序；
      标定接触域、机器人接口、增益/速率/力矩限制；
      在首个描述符下求数值锚点，认证根域、误差及灵敏度；
      固定性能比较域与允许误差，分别建立锚点记录和发布记录。

每个估计窗口：
  1. 同步并补偿位姿与力/力矩，按式（14a）—（17）形成观测。
  2. 移位热启动 MHE，检查数值准则，提取参数、边界、表面状态及协方差。
  3. 建立新模型及描述符，计算相对最近认证锚点的变化。
  4. 原分支证书覆盖时按式（42）预测并计算式（43）的误差界。
  5. 从发布时的状态/相位建立 rollout，检查式（45）及其完整条件。
     通过则直接保存候选、参数、时间标签和发布证书，保留实际增益与能量。
  6. 未通过且局部校正条件成立时，用式（44）校正后重新检查；
     需要重锚时重新认证锚点和灵敏度；原配置/分支失效时重新求解并认证。
     未获得发布证书则不启用该候选；发布本身不重置锚点误差。

每个低层控制周期：
  7. 读取最新 q_r、w_r、g_a、E_T 以及已发布策略的实际时间相位。
  8. 求式（37）的四维受约束动作，按式（29）、（47）更新能量与增益并发送力矩。
  9. 检查发布证书的有效范围；失效时停止凭该证书执行，按既有经验证接口处理。
记录：估计误差、根残差、轨迹裕量、阻抗响应、全部计算及执行延迟。
```

**选择理由与适用条件。** 估计窗口更新参数化问题，低层周期执行同一模型中的反馈；认证锚点支持多次候选查询。附录 H—M 给出核心推导和机制核对，认证证明第 1—8 项给出完整条件。实际平台还需核验认证域、误差界和运行成本。

## 参考文献

[1] Beber, L. et al. *A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography*. ICRA (2024). [作者公开版本](https://arxiv.org/abs/2309.14893)。

[2] Wijayarathne, L., Zhou, Z., Zhao, Y. & Hammond III, F. L. *Real-Time Deformable-Contact-Aware Model Predictive Control for Force-Modulated Manipulation*. IEEE Transactions on Robotics 39(5), 3549–3566 (2023). III-A 为接触模型，VII-F 为材料参数辨识。[DOI](https://doi.org/10.1109/TRO.2023.3286070)。

[3] Lubarda, V. A. & Lubarda, M. V. *On the depth of cylindrical indentation of an elastic half-space for different types of displacement boundary conditions*. European Journal of Mechanics A/Solids 107, 105362 (2024). [DOI](https://doi.org/10.1016/j.euromechsol.2024.105362)；[作者公开版本](https://arxiv.org/abs/2309.02361)。

[4] Rao, C. V., Rawlings, J. B. & Mayne, D. Q. *Constrained state estimation for nonlinear discrete-time systems: Stability and moving horizon approximations*. IEEE Transactions on Automatic Control 48(2), 246–258 (2003). [DOI](https://doi.org/10.1109/TAC.2002.808470)。

[5] Kong, H., Peng, G., Li, G. & Yang, C. *Neural-Network-Based Optimal Impedance Control for Robots in Physical Interaction With Soft Environments*. IEEE Transactions on Systems, Man, and Cybernetics: Systems 55(10), 6463–6475 (2025). III-A—D 给出状态空间、HJB、single critic 与原文执行接口。[DOI](https://doi.org/10.1109/TSMC.2025.3579017)。

[6] Rai, A., Mou, S. & Anderson, B. D. O. *Neighboring Extremal Optimal Control Theory for Parameter-Dependent Closed-loop Laws*. arXiv:2312.04733v1 (2023). IV 与 VI 为本处核对的参数灵敏度和有限基方程来源；本文的受约束根 corrector 另行推导。[原文](https://arxiv.org/abs/2312.04733v1)。

[7] Xu, M., Molloy, T. & Gould, S. *Revisiting Implicit Differentiation for Learning Problems in Optimal Control*. Advances in Neural Information Processing Systems 36 (2023). 第 4 节给出微分 KKT 系统；其轨迹稀疏复杂度不直接移用于本文配置点根。[原文](https://papers.nips.cc/paper_files/paper/2023/file/bcfcf7232cb74e1ef82d751880ff835b-Paper-Conference.pdf)。
