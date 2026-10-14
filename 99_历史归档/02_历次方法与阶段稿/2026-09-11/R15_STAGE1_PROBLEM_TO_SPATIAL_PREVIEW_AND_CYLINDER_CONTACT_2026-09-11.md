# 研究问题与前半部分技术路线：有限接触空间预瞄与圆柱–组织接触动力学

日期：2026-09-11  
定位：面向 T-RO/IJRR 方法论文的前半部分 Method 草稿。  
范围：研究问题、方法概览、spatial preview、圆柱–组织法向/切向接触模型；不涉及 NNBO/HJB、ANEOC 和后续阻抗求解。

## 1. 研究问题

### 1.1 连续机器人超声中的接触适应问题

机器人超声需要在探头沿体表连续移动的同时维持稳定、适度且可重复的探头–组织接触。接触力不足会破坏声学耦合并导致图像中断；接触力过大则会增加患者不适并改变组织几何。探头压力及其空间分布还会影响图像质量、组织定位和弹性测量结果，因此接触控制并不是独立于成像任务的底层伺服问题 [1,2]。

连续扫查中的困难来自接触环境同时具有空间变化和动态性。组织刚度、厚度、表面曲率和支撑条件会沿扫描路径改变，患者呼吸和身体微动还会改变相对表面位置。现有机器人超声方法已经能够利用导纳控制、显式力/力矩控制、扰动观测、在线组织参数估计或模型预测控制改善当前接触处的力跟踪 [3,4]。另一类方法通过扫描前的点式触诊建立患者特异性刚度地图，再将该地图用于扫描或形变校正 [5]。这些方法分别强化了反应速度、当前环境估计或先验环境建模，但没有直接利用一个在连续扫查中始终存在的物理结构：超声探头与组织之间不是点接触，而是具有非零空间范围的接触域。

设探头上一个固定的任务参考点沿扫描路径运动。当组织性质发生空间变化时，探头接触域的前缘先进入新区域，任务参考点随后才到达该区域。前缘接触已经会改变合力和力矩，但新区域尚未作用于任务参考点。因此，有限接触域在连续运动中提供了一个相对于任务参考点的内生空间预瞄窗口。本文研究的问题是：

> 在不进行患者特异性离线逐点触诊的条件下，能否从标准六维力/力矩传感器测得的聚合接触 wrench 中，恢复即将到达任务参考点的组织边界距离和有效接触响应，并进一步建立对应的未来局部接触动力学？

这里的 spatial preview 严格指“边界到达任务参考点之前”的预瞄，而不是“任何接触扰动或图像变化发生之前”的预测。该定义既保留了有限接触域提供的真实时间优势，也避免把前缘已经发生的接触解释为非接触式预知。

### 1.2 方法的基本思路

方法前半部分包含两个连续的物理与计算环节。

第一，利用有限接触域中的分布压力与聚合 wrench 之间的空间矩关系，将探头前缘进入新区域后的合力和力矩变化转化为边界相对距离与有效接触响应的在线估计。该模块的输入为六维 F/T、探头位姿和速度、局部表面几何及当前接触状态，输出为边界距离、响应变化、剩余预瞄时间及相应的不确定性。

第二，根据实际超声探头的圆柱曲面几何，建立法向线接触、组织黏弹变形、切向滑动和机器人动力学之间的耦合模型。圆柱 Hertz 理论用于确定名义接触半宽与压力分布，而有限组织厚度、边界约束和速率效应通过局部标定的黏弹动力学描述。该模型将 spatial preview 的几何–力学估计转换为未来局部接触模型，为后续接触策略计算提供状态和参数接口。

整体计算关系为

\[
\{\mathbf F_s,\mathbf M_s,\mathbf T(q),\dot q,\widehat{\mathcal S}\}
\longrightarrow
\{\hat r_b,\Delta\hat q,\hat T_{\mathrm{prev}}\}
\longrightarrow
\widehat{\mathcal M}_{c,\mathrm{future}},
\]

其中 \(\widehat{\mathcal S}\) 表示局部表面几何，\(r_b\) 是边界相对任务参考点的距离，\(\Delta q\) 是未来区域相对于当前区域的有效接触响应变化，\(\widehat{\mathcal M}_{c,\mathrm{future}}\) 是与观测一致的未来局部接触模型。

## 2. 有限接触域产生的空间预瞄

### 2.1 接触坐标与预瞄区间

令 \(s_c(t)\) 表示任务参考点在组织表面扫描曲线上的弧长位置，扫描速度为

\[
v_s(t)=\dot s_c(t)>0.
\]

任务参考点必须在系统设计中固定为探头接触足迹的几何中心、控制作用点或指定声学参考线之一。以下称其为任务中心。建立探头局部接触坐标系 \(\{C\}\)：\(x\) 轴沿扫描方向，\(y\) 轴沿探头圆柱轴，\(z\) 轴沿组织外法向。设当前接触域在扫描方向上的投影为

\[
\mathcal P(t)=
\{s_c(t)+x\mid x\in[-a(t),a(t)]\},
\]

其中 \(a(t)\) 是实际接触半宽，而不是探头外壳的名义尺寸。扫描前缘对应 \(x=a(t)\)。

设局部组织转变边界位于 \(s_b(t)\)，定义它相对于任务中心的有符号距离

\[
r_b(t)=s_b(t)-s_c(t).
\]

当

\[
0<r_b(t)<a(t)
\]

时，接触域前缘已经越过边界，而任务中心仍位于当前区域。这一区间构成 spatial preview interval。若边界与组织表面一起运动，则其相对运动满足

\[
\dot r_b(t)=v_b(t)-v_s(t),
\]

其中 \(v_b\) 是边界沿扫描曲线的表观速度；对静止组织边界，\(\dot r_b=-v_s\)。

### 2.2 从传感器 wrench 到接触 wrench

六维 F/T 传感器通常位于探头后方，传感器力矩同时包含接触压力矩和力臂产生的力矩。令 \(O_s\) 为传感器原点，\(O_c\) 为任务中心处的接触参考点，\(\mathbf r_{sc}\) 为从 \(O_s\) 指向 \(O_c\) 的向量。将传感器 wrench 旋转到接触坐标系后，接触点处的力和力矩为

\[
\mathbf F_c=\mathbf R_{cs}\mathbf F_s,
\]

\[
\mathbf M_c=
\mathbf R_{cs}\mathbf M_s
-\left(\mathbf R_{cs}\mathbf r_{sc}\right)\times\mathbf F_c.
\]

重力、机器人运动惯性、传感器偏置和线缆载荷在该变换前后分别补偿。由此得到的法向力 \(F_n\) 与绕 \(y\) 轴的俯仰力矩 \(M_y\) 才能用于分析扫描方向上的压力分布。

令 \(p(x,y,t)\) 为法向接触压力，将其沿圆柱轴方向积分得到线载荷密度

\[
\lambda(x,t)=\int_{\mathcal Y(t)}p(x,y,t)\,\mathrm dy.
\]

在接触坐标原点位于任务中心、并将切向牵引造成的力矩单独补偿后，

\[
F_n(t)=\int_{-a(t)}^{a(t)}\lambda(x,t)\,\mathrm dx,
\]

\[
M_y(t)=\int_{-a(t)}^{a(t)}x\lambda(x,t)\,\mathrm dx.
\]

因此，\(F_n\) 是接触压力的零阶空间矩，\(M_y\) 是扫描方向的一阶空间矩。合力描述接触响应变化的总量，力矩同时携带该变化在接触域内的位置。这种由分布压力/牵引积分得到等效合力与力矩的表示与一般接触 wrench 建模一致 [8]。

### 2.3 边界进入前缘后的观测模型

在一个短空间窗口内，先考虑单一、局部近似横穿扫描方向的 step transition。令 \(\lambda_0(x,t)\) 是接触域完全位于当前区域时的基线线载荷密度。新区域进入前缘后，将增量压力在当前工作点附近线性化为

\[
\lambda(x,t)=
\lambda_0(x,t)
+\Delta q(t)\,k(x;a,\chi)\,
H\!\left(x-r_b(t)\right)
+e_\lambda(x,t).
\]

其中 \(H(\cdot)\) 是 Heaviside 函数，\(\Delta q\) 表示新旧区域在当前预载、压入速度和探头几何下的有效接触响应差异，\(k(x;a,\chi)>0\) 是经名义接触模型和均匀材料标定得到的压力灵敏度核，\(\chi\) 收集当前预载与几何条件，\(e_\lambda\) 表示局部线性化和跨边界弹性耦合造成的残差。这里的 \(\Delta q\) 不是完整 Young 模量或完整黏弹参数向量，而是当前自然扫查激励下可由聚合 wrench 分辨的有效响应对比。

从实测 wrench 中减去基线预测和已补偿项，得到增量观测

\[
\widetilde F_n
=\Delta q\,A(r_b;a,\chi)+\varepsilon_F,
\]

\[
\widetilde M_y
=\Delta q\,B(r_b;a,\chi)+\varepsilon_M,
\]

其中

\[
A(r;a,\chi)=\int_r^a k(x;a,\chi)\,\mathrm dx,
\qquad
B(r;a,\chi)=\int_r^a xk(x;a,\chi)\,\mathrm dx.
\]

\(\varepsilon_F\) 和 \(\varepsilon_M\) 汇总传感噪声、基线误差以及 \(e_\lambda\) 的积分效应。

### 2.4 边界距离与响应变化的局部可辨识性

忽略有界观测误差，并假设 \(\Delta q\neq0\)。由两个增量观测构造

\[
\rho(r_b)=
\frac{\widetilde M_y}{\widetilde F_n}
=\frac{B(r_b)}{A(r_b)}.
\]

因为

\[
A'(r)=-k(r),
\qquad
B'(r)=-rk(r),
\]

所以

\[
\begin{aligned}
\rho'(r)
&=
\frac{k(r)\left[B(r)-rA(r)\right]}{A(r)^2}\\
&=
\frac{
k(r)\int_r^a(x-r)k(x)\,\mathrm dx
}{
A(r)^2
}>0,
\qquad 0<r<a.
\end{aligned}
\]

因此，只要灵敏度核在接触区内为正，\(\rho(r)\) 对边界距离严格单调。对于已知 \(a\) 和 \(\chi\)，边界距离可由

\[
r_b=\rho^{-1}
\left(
\frac{\widetilde M_y}{\widetilde F_n}
\right)
\]

唯一确定；随后有效响应差异由

\[
\Delta q=
\frac{\widetilde F_n}
{A(r_b;a,\chi)}
\]

恢复。这一结果给出了有限接触域空间预瞄的核心数学依据：零阶矩确定响应变化的总量，一阶矩提供变化在接触域中的质心位置，两者联合分离边界位置与响应幅值。

以均匀灵敏度核 \(k(x)=c_\delta>0\) 为例，令新区域进入接触域的长度为

\[
\ell=a-r_b.
\]

则

\[
\widetilde F_n=c_\delta\Delta q\,\ell,
\]

\[
\widetilde M_y
=c_\delta\Delta q
\left(a\ell-\frac{\ell^2}{2}\right),
\]

进而

\[
\ell=
2\left(
a-\frac{\widetilde M_y}{\widetilde F_n}
\right),
\qquad
\Delta q=
\frac{\widetilde F_n}{c_\delta\ell}.
\]

该特例直观展示了单时刻法向力的乘积歧义：同一 \(\widetilde F_n\) 可以由较大的响应差异与较短的侵入长度产生，也可以由较小的响应差异与较长的侵入长度产生。加入一阶力矩后，两者在非零侵入和非零响应对比下能够分离。已知扫描运动下，较长的 force–pose history 也可能逐渐提供空间信息，因此本文使用力矩的目的不是宣称 force-only 在所有时间尺度上不可辨识，而是在有限预瞄窗口内提高最早阶段的空间分离能力。

### 2.5 基于短时历史的 spatial-preview estimator

上述闭式关系解释了可辨识性，但不能直接作为实际估计器。一方面，边界刚进入前缘时 \(\widetilde F_n\) 很小，比值对噪声敏感；另一方面，接触半宽、传感偏置、边界方向和组织运动均可能在扫描过程中变化。因此，实际算法在长度为 \(H\) 的短时窗口内联合拟合几何状态与接触响应。

在离散时刻 \(k\)，定义预瞄状态

\[
\mathbf z_k=
\begin{bmatrix}
r_{b,k} &
\Delta q_k &
a_k &
b_{F,k} &
b_{M,k}
\end{bmatrix}^{\!\top},
\]

其中 \(b_F,b_M\) 是剩余的慢变 wrench 偏置。状态演化采用

\[
r_{b,k+1}
=r_{b,k}
-\Delta s_k
+\Delta s_{b,k}
+w_{r,k},
\]

\[
\Delta q_{k+1}
=\Delta q_k+w_{q,k},
\]

\[
a_{k+1}
=a_k+w_{a,k},
\]

\[
b_{F,k+1}=b_{F,k}+w_{b_F,k},
\qquad
b_{M,k+1}=b_{M,k}+w_{b_M,k}.
\]

\(\Delta s_k\) 是任务中心在相邻采样间隔内的扫描位移，\(\Delta s_{b,k}\) 是边界随组织运动产生的位移。接触半宽在短窗口内按慢变状态传播，同时由第 3 节圆柱接触模型给出的名义关系

\[
a_k=g_a(F_{n,k},\delta_{n,k},\hat\eta_{n,k})+e_{a,k}
\]

作为软约束；\(e_{a,k}\) 吸收有限长度、组织边界和压力核失配。过程项 \(w_{\cdot,k}\) 允许小范围模型变化。

相对于当前区域基线的观测模型为

\[
\widetilde{\mathbf y}_k=
\begin{bmatrix}
\widetilde F_{n,k}\\
\widetilde M_{y,k}
\end{bmatrix}
=
\begin{bmatrix}
\Delta q_k A(r_{b,k};a_k,\chi_k)+b_{F,k}\\
\Delta q_k B(r_{b,k};a_k,\chi_k)+b_{M,k}
\end{bmatrix}
+\mathbf v_k.
\]

在每个时刻求解有界移动窗估计问题

\[
\begin{aligned}
\min_{\mathbf z_{k-H:k}}
\quad&
\|\mathbf z_{k-H}-\bar{\mathbf z}_{k-H}\|_{P^{-1}}^2\\
&+
\sum_{i=k-H}^{k}
\|
\widetilde{\mathbf y}_i-h(\mathbf z_i,\mathbf u_i)
\|_{R^{-1}}^2\\
&+
\sum_{i=k-H}^{k-1}
\|
\mathbf z_{i+1}-f_z(\mathbf z_i,\mathbf u_i)
\|_{Q^{-1}}^2
\end{aligned}
\]

并施加

\[
0\le r_{b,i}\le a_i,\qquad
a_i\in[\underline a_i,\overline a_i],\qquad
\Delta q_i\in\mathcal Q,
\]

以及由传感器标定给出的偏置和噪声界。估计器输出

\[
\widehat{\mathcal Z}_{\mathrm{future}}
=
\{
\hat r_b,\,
\Delta\hat q,\,
\hat a,\,
\hat m_c,\,
\mathcal U_p
\},
\]

其中 \(\hat m_c\) 是接触/转变模式，\(\mathcal U_p\) 是由估计协方差或 set-membership 可行域得到的不确定性描述。当观测不足以区分转变与偏置、曲率或摩擦变化时，估计器保持不确定模式，不输出确定的未来材料参数。

### 2.6 预瞄时间

对于局部静止边界和近似恒定的相对扫描速度，

\[
T_{\mathrm{prev}}=
\frac{r_b}{v_s}.
\]

若速度变化或组织边界运动不可忽略，则预瞄时间定义为边界与任务中心的首次相遇时间：

\[
T_{\mathrm{prev}}
=
\inf
\left\{
T\ge0:
\int_t^{t+T}
\left[v_s(\tau)-v_b(\tau)\right]
\,\mathrm d\tau
\ge r_b(t)
\right\}.
\]

估计器输出距离区间时，后续计算使用该区间与相对速度不确定性共同得到的保守预瞄时间下界。

## 3. 面向超声探头的圆柱–组织接触模型

### 3.1 几何选择与模型角色

原 T-RO deformable-contact-aware formulation 将球形刚性工具与软平面之间的 Hertz 法向接触、滑动摩擦和机器人动力学共同纳入预测模型 [6]。本研究沿用其“接触动力学与机器人运动联合建模”的结构，但将球形工具替换为更符合目标超声探头的圆柱曲面。球–平面 Hertz 的 \(F_n\propto\delta_n^{3/2}\) 关系和圆形接触斑上的切向积分不能直接用于圆柱线接触。

设探头接触面在扫描方向 \(x\) 上具有曲率半径 \(R_p\)，在 \(y\) 方向具有有效长度 \(L_{\mathrm{eff}}\)。组织表面在局部接触尺度内近似为曲率半径 \(R_t\) 的圆弧。圆柱轴近似沿 \(y\) 方向，且 \(L_{\mathrm{eff}}\) 明显大于扫描方向接触半宽，使接触域中部可采用二维线接触近似。若最终硬件为平面线阵面、双曲率面，或者扫描方向不沿圆柱弯曲主方向，则应替换相应几何接触模型；后续推导中的圆柱模型只适用于上述几何。

### 3.2 Hertz 线接触的接触半宽与压力分布

探头和组织的约化模量为

\[
\frac{1}{E^*}
=
\frac{1-\nu_p^2}{E_p}
+
\frac{1-\nu_t^2}{E_t},
\]

等效曲率半径为

\[
\frac{1}{R^*}
=
\frac{1}{R_p}
+
\frac{1}{R_t}.
\]

对于相对组织近似刚性的探头和局部平面组织，\(E_p\gg E_t\)、\(R_t\rightarrow\infty\)，因此

\[
E^*\simeq\frac{E_t}{1-\nu_t^2},
\qquad
R^*\simeq R_p.
\]

令总法向载荷为 \(F_n\)，单位圆柱轴长度上的载荷为

\[
w=\frac{F_n}{L_{\mathrm{eff}}}.
\]

在小接触宽度、线弹性、无粘附和近似无摩擦的法向接触条件下，Hertz 线接触给出接触半宽 [7]

\[
a_H=
\sqrt{
\frac{4wR^*}{\pi E^*}
},
\]

以及半椭圆压力分布

\[
p_H(x)
=
p_0
\sqrt{
1-\frac{x^2}{a_H^2}
},
\qquad
|x|\le a_H,
\]

\[
p_0=
\frac{2w}{\pi a_H}.
\]

该分布满足

\[
\int_{-a_H}^{a_H}p_H(x)\,\mathrm dx=w.
\]

因此，Hertz 线接触为 spatial-preview estimator 提供两个直接量：随载荷和局部响应变化的名义接触半宽 \(a_H\)，以及用于构造正灵敏度核 \(k(x)\) 的半椭圆压力形状。

### 3.3 接触区高度与总压痕深度

圆柱线接触中，圆柱帽接触区的几何高度为

\[
\delta_0
=
\frac{a_H^2}{2R^*}
=
\frac{2w}{\pi E^*}.
\]

\(\delta_0\) 与单位长度载荷 \(w\) 呈线性关系，但它不是相对于远处未变形组织表面的总压痕深度。二维线载荷引起的位移具有对数远场特性，因此总压痕深度依赖位移基准、组织厚度和侧向边界 [7]。

例如，将组织自由表面上 \(x=\pm B_0\) 处定义为零位移，令

\[
\beta=\frac{B_0}{a_H}>1,
\]

理想弹性半空间中的总压痕深度为

\[
\delta_n
=
\delta_0
\left[
\beta^2
-\beta\sqrt{\beta^2-1}
+\ln\!\left(
\beta+\sqrt{\beta^2-1}
\right)
\right].
\]

当 \(\beta\gg1\) 时，

\[
\delta_n
\simeq
\delta_0
\left(
\ln\frac{2B_0}{a_H}
+\frac{1}{2}
\right).
\]

这说明，圆柱接触不能通过简单地把球形 Hertz 的指数从 \(3/2\) 改成 \(1\) 得到真实的总力–压痕关系。对有限厚度软组织，法向接触映射必须包含组织厚度、支撑边界、有限探头长度和当前预载。

### 3.4 有限组织上的局部黏弹法向动力学

为了在控制时间尺度上描述真实探头–组织交互，本文使用 Hertz 线接触作为几何与压力先验，并在规定工作域内建立有限几何的局部黏弹模型。令

\[
w(t)
=
h_\infty
\left(
\delta_n(t);
\eta_e
\right)
+q_v(t)
+r_w(t),
\]

其中 \(h_\infty\) 是在实际 \(R_p,L_{\mathrm{eff}}\)、组织厚度、支撑条件和预载范围下标定的平衡载荷–压入映射，\(q_v\) 是描述速率相关响应和应力松弛的内部变量，\(r_w\) 是有界建模残差。

采用一个最小 generalized-Maxwell 分支，

\[
\dot q_v
=
k_v(\delta_n;\eta_v)\dot\delta_n
-\frac{q_v}{\tau_v},
\qquad
k_v\ge0,\quad \tau_v>0.
\]

于是

\[
F_n
=
L_{\mathrm{eff}}
\left[
h_\infty(\delta_n;\eta_e)
+q_v
\right]
+d_F,
\]

且在局部参数不发生瞬时变化时，

\[
\dot F_n
=
L_{\mathrm{eff}}
\left[
\left(
\frac{\partial h_\infty}{\partial\delta_n}
+k_v
\right)
\dot\delta_n
-\frac{q_v}{\tau_v}
\right]
+d_n.
\]

\(d_F=L_{\mathrm{eff}}r_w\)，\(d_n\) 汇总残差变化和未建模动态。为保证局部接触映射可逆并避免负切线刚度，工作域内要求

\[
\frac{\partial h_\infty}{\partial\delta_n}
+k_v(\delta_n)>0.
\]

该模型保留了 T-RO formulation 中“法向接触力随探头运动和组织内部状态演化”的动力学结构，同时将球形 Hertz 特定的指数替换为由圆柱几何、有限边界和黏弹数据共同确定的局部映射。若实验显示一个松弛时间常数不能覆盖目标速度范围，可增加 Maxwell 分支；在此之前不引入额外不可辨识状态。

### 3.5 法向接触运动学

令组织外法向 \(\mathbf n\) 从组织指向探头，探头与组织表面在接触点的速度分别为 \(\mathbf v_p\) 和 \(\mathbf v_s\)。定义向组织内部的压入为正，则

\[
\dot\delta_n
=
-\mathbf n^\top
\left(
\mathbf v_p-\mathbf v_s
\right)
+d_\delta.
\]

\(d_\delta\) 表示局部法向估计误差和未建模表面运动。令 \(g_n\) 表示刚性探头表面与未变形组织表面之间的非负法向间隙，则接触互补条件为

\[
F_n\ge0,\qquad
g_n\ge0,\qquad
F_ng_n=0.
\]

在 contact 模式下 \(g_n=0\)，\(\delta_n\ge0\) 表示组织内部的压缩变形；在 no-contact 模式下 \(F_n=0\)，黏弹接触模型停止产生法向输出。这样可以避免把组织压缩量 \(\delta_n\) 与几何间隙 \(g_n\) 混为同一变量，也不会在接触丢失后产生拉伸法向力。

### 3.6 切向滑动

探头在凝胶耦合条件下沿组织表面连续移动，切向作用同时包含界面摩擦和组织表层的速率相关拖曳。令

\[
\mathbf v_t
=
\left(
\mathbf I-\mathbf n\mathbf n^\top
\right)
\left(
\mathbf v_p-\mathbf v_s
\right)
\]

为相对切向速度。接触切向模型采用 stick–slip 分支：

\[
\mathrm{stick}:
\qquad
\mathbf v_t\approx\mathbf 0,
\quad
\|\mathbf F_t\|\le\mu_sF_n,
\]

\[
\mathrm{slip}:
\qquad
\mathbf F_t
=
-\left(
\mu_kF_n+c_t\|\mathbf v_t\|
\right)
\frac{\mathbf v_t}{\|\mathbf v_t\|}
+\mathbf d_t.
\]

\(\mu_s\) 和 \(\mu_k\) 分别是静、动摩擦系数，\(c_t\) 是局部切向黏性/形变拖曳系数，\(\mathbf d_t\) 包含凝胶膜、方向依赖和未建模剪切变形。本文的连续扫查主模型限定在 \(\|\mathbf v_t\|\ge v_{\min}\) 的 slip 分支；低速、反向和粘着转换通过接触模式单独处理。若后续求解需要光滑动力学，可以在 slip 分支内用

\[
\frac{\mathbf v_t}
{\sqrt{\|\mathbf v_t\|^2+\varepsilon_v^2}}
\]

近似单位方向，但该数值正则化不改变真实接触模式的定义。

### 3.7 机器人–接触耦合动力学

令

\[
\mathbf W_c=
\begin{bmatrix}
\mathbf F_c\\
\mathbf M_c
\end{bmatrix}
\]

表示环境作用在探头上的接触 wrench，其中

\[
\mathbf F_c=F_n\mathbf n+\mathbf F_t.
\]

机器人动力学为

\[
\mathbf M(q)\ddot q
+\mathbf C(q,\dot q)\dot q
+\mathbf g(q)
=
\tau
+\mathbf J_c(q)^\top\mathbf W_c
+\tau_d.
\]

这里使用“环境作用于机器人”的 wrench 符号约定；若实现采用机器人作用于组织的 wrench，接触项符号应整体反向。

在持续接触和滑动分支内，一个最小物理状态可写为

\[
\mathbf x_c=
\begin{bmatrix}
q^\top &
\dot q^\top &
\delta_n &
q_v
\end{bmatrix}^{\!\top}.
\]

\(\delta_n\) 由法向接触运动学演化，\(q_v\) 描述组织黏弹记忆，\(F_n\) 和 \(\mathbf F_t\) 由接触状态计算。若表面静止且几何完全已知，\(\delta_n\) 可以由机器人位姿代数计算并从状态中消去；若存在组织运动或表面估计误差，则保留 \(\delta_n\) 能形成更清楚的动态状态。

## 4. Spatial preview 与未来接触模型的接口

spatial-preview estimator 给出的 \(\Delta\hat q\) 是在当前预载和自然扫描速度下可观测的有效响应变化，它通常不足以唯一确定 \(E_t,\nu_t,\tau_v,\mu_k\) 等完整参数。因此，预瞄模块不直接输出一个无不确定性的未来材料参数向量，而是利用 \(\hat r_b\)、\(\Delta\hat q\)、当前圆柱接触模型和参数先验，构造与 wrench 历史一致的未来局部接触模型集合：

\[
\widehat{\mathcal M}_{c,\mathrm{future}}
=
\left\{
\mathcal M_c(\eta_c):
\eta_c\in
\widehat{\mathcal E}_{\mathrm{future}}
\right\}.
\]

其中

\[
\eta_c=
\{
\eta_e,\eta_v,\tau_v,\mu_k,c_t,R^*,L_{\mathrm{eff}}
\}
\]

描述圆柱–组织接触动力学，\(\widehat{\mathcal E}_{\mathrm{future}}\) 由当前模型、预瞄响应差异及其不确定性共同约束。对只提供法向自然激励的连续扫查，未来集合首先在局部有效法向切线响应方向上收紧；未被当前观测激励的黏弹时间常数或摩擦参数继续保留为区间，而不是由一次边界进入事件强行辨识。

至此，前半部分技术路线形成闭合关系：

\[
\boxed{
\begin{aligned}
&\text{finite probe contact}\\
&\quad\longrightarrow
\text{wrench spatial moments}\\
&\quad\longrightarrow
\{\hat r_b,\Delta\hat q,\hat T_{\mathrm{prev}}\}\\
&\quad\longrightarrow
\text{cylinder–tissue normal/tangential dynamics}\\
&\quad\longrightarrow
\widehat{\mathcal M}_{c,\mathrm{future}}.
\end{aligned}
}
\]

后续方法只需以 \(\widehat{\mathcal M}_{c,\mathrm{future}}\) 和 \(\hat T_{\mathrm{prev}}\) 为输入，计算并在任务中心到达边界前实现相应的接触策略。本稿不预设该后续策略的具体求解形式。

## 5. 本部分的假设与适用范围

上述推导针对连续压缩接触下的局部空间转变。解析可辨识性结果要求短窗口内存在单个低维 step transition、扫描方向已知、灵敏度核为正、接触半宽有界且响应对比非零。对于有限宽度 ramp、斜边界或变化曲率，同一分布压力模型仍可使用，但边界位置、方向和宽度需要在短时历史中联合估计，不能仅靠单时刻的 \(F_n/M_y\) 比值闭式恢复。

圆柱 Hertz 部分要求探头几何与扫描方向符合线接触近似，并且接触半宽相对曲率半径较小。真实软组织的有限厚度、近不可压缩性、非线性和黏弹性由 \(h_\infty,q_v\) 与有界残差描述；当压入、组织折皱、粘附或端部效应超出标定域时，应更换有限面超弹/黏弹接触模型，但有限接触域产生 spatial preview 的几何机制本身不依赖球形 Hertz 公式。

## 参考文献

[1] S. Suchoń, M. Burkacki, M. Chrzan, and M. Winder, “What Ranges of Probe Pressure Are Applied During Ultrasound Examinations? A Systematic Review,” Sensors, 2025. https://doi.org/10.3390/s25113415

[2] A. Y. Huang and B. W. Anthony, “An Instrumented Ultrasound Probe for Shear Wave Elastography With Uneven Force Distribution,” EMBC, 2019. https://doi.org/10.1109/EMBC.2019.8856782

[3] A. Piedra, R. B. Jeffrey, and O. Khatib, “Explicit force-based control strategies for compliant robot-assisted ultrasound imaging,” The International Journal of Robotics Research, 2026. https://doi.org/10.1177/02783649261441635

[4] L. Beber et al., “A Passive Variable Impedance Control Strategy with Viscoelastic Parameters Estimation of Soft Tissues for Safe Ultrasonography,” ICRA, 2024. https://doi.org/10.1109/ICRA57147.2024.10610167

[5] Z. Jiang et al., “Deformation-Aware Robotic 3D Ultrasound,” IEEE Robotics and Automation Letters, 2021. https://doi.org/10.1109/LRA.2021.3099080

[6] L. Wijayarathne, Z. Zhou, Y. Zhao, and F. L. Hammond III, “Real-Time Deformable-Contact-Aware Model Predictive Control for Force-Modulated Manipulation,” IEEE Transactions on Robotics, 2023. https://doi.org/10.1109/TRO.2023.3286070

[7] V. A. Lubarda and M. V. Lubarda, “On the depth of cylindrical indentation of an elastic half-space for different types of displacement boundary conditions,” European Journal of Mechanics - A/Solids, 2024. https://doi.org/10.1016/j.euromechsol.2024.105362

[8] R. Elandt, E. Drumwright, M. Sherman, and A. Ruina, “A pressure field model for fast, robust approximation of net contact force and moment between nominally rigid objects,” IROS, 2019. https://arxiv.org/abs/1904.11433

## 补充：spatial preview 的 prior、定量可实现性与失败边界

**Prior 边界。** 车辆 wheelbase preview 已证明“上游先遭遇的扰动可为下游受限速机构提供重构时间”，且性能取决于 preview time [9,10]；接触研究已证明合力/力矩是压力场的低阶矩 [8]，并可用于估计接触斑 [11]；超声中也已有压力阵列反馈 [12]。这些工作分别支持“预瞄控制”“wrench 空间信息”和“超声压力反馈”，但都没有给出本文所需的映射：由单个探头的标准聚合 wrench 在线恢复 *pre-task-center* 组织转变距离。因此本文不能声称首次 preview control 或首次 pressure feedback，核心仍是该目标域映射及其时序用途。

**量级与时延算例。** 以公开 4C-RS 凸阵探头的 \(R_p=60\,\mathrm{mm}\)、接触面 \(66.2\times18.3\,\mathrm{mm}\) [13]，公开超声机器人采用的 \(F_n=5\,\mathrm N\) [14]，以及人体臀部肌肉在 \(5\,\mathrm N\) 下报告的代表性 \(E_t=62.2\,\mathrm{kPa}\) [15] 为跨文献量级算例；取近不可压缩的 \(\nu_t=0.49\) 和 \(L_{\mathrm{eff}}=18.3\,\mathrm{mm}\)，则

\[
E^*=81.9\,\mathrm{kPa},\qquad
w=273.2\,\mathrm{N/m},\qquad
a_H=16.0\,\mathrm{mm},\qquad
\delta_0=2.13\,\mathrm{mm}.
\]

这说明毫米至厘米级空间窗口在量级上可能存在，但 \(a_H/R_p=0.27\)、\(L_{\mathrm{eff}}/a_H=1.15\)，并不满足严格的 \(a_H\ll R_p\) 与 \(L_{\mathrm{eff}}\gg a_H\)。因此 \(a_H\) 只能作 observer 的软先验；控制保证必须使用压痕实验得到的实际承载半宽及检测侵入量，令 \(r_{b,\mathrm{lower}}=a_{\mathrm{meas,lower}}-\ell_{\mathrm{det,upper}}\)。

公开机器人超声实验已覆盖 \(v_s=10\,\mathrm{mm/s}\) [16]；若实测保守剩余距离仅为 \(r_{b,\mathrm{lower}}=5\,\mathrm{mm}\)，则 \(T_{\mathrm{prev,lower}}=0.50\,\mathrm s\)。沿现有 RealMan Gen-4 实现，其官方 UDP 周期为 \(5\,\mathrm{ms}\) [17]。采用待实机核验的保守 P99 预算

\[
\tau_{\mathrm{total}}
=5\ \text{(采样)}+100\ \text{(滤波/MHE)}+20\ \text{(求解/通信)}+300\ \text{(实际响应重构)}
=425\,\mathrm{ms},
\]

得到 \(T_{\mathrm{margin}}=75\,\mathrm{ms}\)，即预瞄窗口的 \(15\%\)：名义可行，但尚不足以满足建议的 \(25\%\) 鲁棒余量。无余量时 \(v_s\le11.8\,\mathrm{mm/s}\)；要求 \(25\%\) 余量时应有 \(v_s\le8.8\,\mathrm{mm/s}\)。若 \(v_s=15\,\mathrm{mm/s}\)，则 \(T_{\mathrm{margin}}=-92\,\mathrm{ms}\)，不可行。上述 100/20/300 ms 是设计验收预算而非已测事实，最终只能用同步日志得到的 P95/P99/max 时延替换。传感幅值也必须过门槛：例如 \(0.5\,\mathrm N\) 增量作用于 \(5\,\mathrm{mm}\) 偏心仅产生 \(2.5\times10^{-3}\,\mathrm{N\,m}\)，对 \(0.001\,\mathrm{N\,m}\) 数值粒度仅约 2.5 个计数，实际 loaded noise、线缆和重力残差可能否决该路线。

**假设、适用范围与失败边界。** 本方法适用于任务中心定义固定、连续压缩接触、扫描方向已知、短窗口内单个低维转变、正灵敏度核、非零响应对比且 wrench/几何 nuisance 可补偿的扫查；圆柱模型还要求扫描方向与圆柱曲率方向一致，并在实测工作域内完成有限面、有限厚度和黏弹标定。出现以下任一结果即不能声称 spatial preview 可实现：

- 在固定误报率下，\(F+M\) 不能比最强的 force-history estimator 更早且更准地恢复边界；
- 实测 transition moment 不高于均匀材料上的 P99 residual，或滤波群时延耗尽提前量；
- \(T_{\mathrm{prev,lower}}\le\tau_{\mathrm{total,P99}}/(1-\gamma)\)（\(\gamma\) 为规定余量），即减速后仍无正时延裕度；
- 转变宽度远大于接触斑、多边界/斜边界不可分、接触丢失、rocking、gel pocket、呼吸动态或摩擦力矩与材料转变同构；
- Hertz 小参数不成立，且有限面标定模型仍不能把接触宽度与几何残差约束到可用范围。

补充参考：[9] Theunissen et al., *Annual Reviews in Control*, 2021, https://doi.org/10.1016/j.arcontrol.2021.03.010；[10] Li et al., *Journal of Sound and Vibration*, 2014, https://doi.org/10.1016/j.jsv.2014.06.017；[11] Ota et al., ICRA 2024, https://arxiv.org/abs/2309.14552；[12] Mousavi et al., *IJCARS*, 2013, https://doi.org/10.1007/s11548-013-0898-3；[13] GE HealthCare 4C-RS specification, https://www.gehealthcare.com/en-sg/-/jssmedia/global/products/files/ultrasound/versanaactive/versana-active-spec-sheet---uscan-jb714712us-oct-2019.pdf；[14] Bao et al., *IEEE TIM*, 2023, https://doi.org/10.1109/TIM.2023.3239925；[15] Isogai et al., 2022, https://doi.org/10.1589/jpts.34.777；[16] Min et al., 2026, https://doi.org/10.4271/2026-99-0751；[17] RealMan Gen-4 UDP documentation, https://develop.realman-robotics.com/en/robot4th/json/udpConfig/。
