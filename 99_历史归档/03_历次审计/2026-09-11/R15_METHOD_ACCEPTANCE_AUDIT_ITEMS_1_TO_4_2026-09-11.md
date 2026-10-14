# R15 Method 验收审计：最终顺序第 1–4 项

日期：2026-09-11  
审计范围：Problem fit、全链条数学推导、模块接口一致性、Method fit / necessity。  
本轮不审计 Scientific strength、novelty 或 venue positioning。  
验收规则：只输出 CLOSED 或 NOT CLOSED；不使用“条件通过”。

## 总结论

# NOT CLOSED

当前研究方向具有真实应用动机，Stage 1 的理想 step-transition 零阶/一阶 wrench 矩可辨识性推导成立，圆柱线接触的名义几何也基本正确。但是，当前完整 Method 仍有 12 个非实验性闭合缺口。其中 5 个直接阻断 HJB/ANEOC 主干：

1. 当前根 failure 把标准探头 F/T 近似成 task-center point sensing，忽略了聚合 wrench 在前缘进入新区域时已经发生变化；
2. Stage 1 的 \(\Delta q\) 尚未通过一个可计算的逆像映射生成 future contact-model set；
3. 跨边界黏弹接触的状态约化与 transport term 尚未推导；
4. preview phase、无限时域 equilibrium HJB 和 ANEOC 的自治/可控性假设彼此不一致；
5. \(F_{\mathrm{HJB}}(W,\eta)\)、critic root uniqueness、约束处理、ANEOC 误差传播与 applied-impedance 可达性定理均未闭合。

这些均属于方法本身的数学问题，不能留给实验回答。

---

## 1. Problem fit

### 1.1 已经有依据的领域痛点

机器人超声确实要求在扫描中同时维持适当接触力、声学耦合与患者安全。过大力会造成组织变形甚至安全风险，过小力会破坏耦合和图像质量。已有机器人超声研究也直接将 force overshoot、steady-state tracking error、变化组织刚度和患者运动列为控制难点。例如，Xie and Yan 的腹部扫描研究明确把瞬时接触力超调和稳态误差作为共同问题，并展示了自适应阻抗方案中的大幅超调；Beber 等在机器人超声中根据黏弹参数在线优化阻抗并加入 energy tank；近期 explicit force control 也以 heterogeneous tissue 和 dynamic environment 下的快速力调节为目标。对应原始资料：

- [Xie and Yan, 2022](https://doi.org/10.5194/ms-13-559-2022)
- [Beber et al., ICRA 2024](https://arxiv.org/abs/2309.14893)
- [Piedra et al., IJRR 2026](https://doi.org/10.1177/02783649261441635)

因此，“连续机器人超声中的组织变化会引起危险或影响成像的接触瞬态”是一个真实问题，不是由当前算法倒推出来的伪问题。

### 1.2 当前 causal failure 的结构性不准确

当前路线把传统控制器的时序写成

\[
\text{future tissue reaches task center}
\rightarrow
\text{wrench changes}
\rightarrow
\text{estimation and control update}.
\tag{1}
\]

这一表述与 Stage 1 的核心物理机制冲突。机器人腕部 F/T 传感器测量的是整个有限接触斑的聚合 wrench。只要未来区域进入前缘，即 \(0<r_b<a\)，Stage 1 模型已经给出

\[
\widetilde F_n=\Delta q A(r_b),\qquad
\widetilde M_y=\Delta q B(r_b).
\tag{2}
\]

因此，使用相同 F/T 信号的普通反馈控制器可能在 task center 到达边界之前就对 \(\widetilde F_n\) 作出反应。不能把所有既有方法都描述为严格的 current-point sensing，也不能声称它们必然到 \(r_b=0\) 后才开始调整。

真正由现有推导支持的 failure 更窄：

> 当未来区域只进入有限接触斑的一部分时，聚合力已经变化，但 force-only response 将 transition extent 与 response contrast 混合为乘积；若控制器没有恢复边界距离、未来模型和剩余执行时间，它可能只能对混合响应作反应，而不能确定 task-center encounter 前应达到的 future impedance。

均匀核特例清楚显示这种混合：

\[
\widetilde F_n=c_\delta\Delta q\ell.
\tag{3}
\]

同一个 \(\widetilde F_n\) 可以来自不同的 \((\Delta q,\ell)\)。这些组合可以具有不同的 future mechanics 和不同的剩余时间。加入

\[
\widetilde M_y
=c_\delta\Delta q
\left(a\ell-\frac{\ell^2}{2}\right)
\tag{4}
\]

后才可在理想条件下分离二者。

### 1.3 必须补上的问题级命题

需要建立一个与真实聚合 wrench 信息结构一致的 impossibility / value-of-preview 命题。其最小形式应为：

令两个转变假设 \(\theta_a=(r_a,\Delta q_a,\eta_a)\) 与
\(\theta_b=(r_b,\Delta q_b,\eta_b)\) 在 force-only 容差内满足

\[
\left\|
\widetilde F_n^{(a)}(0:t)
-\widetilde F_n^{(b)}(0:t)
\right\|
\le\varepsilon_F,
\tag{5}
\]

但它们要求的 task-center-crossing safe-gain sets 不相容：

\[
\mathcal G_{\mathrm{cross}}(\theta_a)
\cap
\mathcal G_{\mathrm{cross}}(\theta_b)
=\varnothing.
\tag{6}
\]

若剩余时间小于从共同初始 applied gain 分别到达两个集合所需的最大重构时间，则任何只使用该 force-only history 的因果策略至少在一个假设下不能保证 crossing constraint。随后再证明 \((F_n,M_y)\) 在声明的 transition family 内打破式 (5) 的歧义，或将歧义缩小到具有共同安全策略的集合。

当前文件没有该问题级命题，也没有把根 failure 从“看不见变化”修正为“看见混合变化但无法及时确定 future action”。

### 第 1 项结论

**NOT CLOSED。**

领域痛点真实，但当前 RCA 对标准聚合 wrench 的信息时序描述过强。必须先修正 failure，再给出式 (5)–(6) 类型的因果不可区分性与 action incompatibility 命题。该修正不会推翻 spatial preview，而是使它解决的 failure 与实际传感器一致。

---

## 2. 从 spatial preview 到 pre-transition execution 的全链条数学审计

### 2.1 Spatial preview 的理想解析部分

在 \(a,\chi,k\) 已知、\(k(x)>0\)、\(\Delta q\neq0\) 且 \(0<r_b<a\) 时，

\[
\rho(r)=\frac{B(r)}{A(r)}
\tag{7}
\]

满足

\[
\rho'(r)
=
\frac{k(r)\int_r^a(x-r)k(x)\,\mathrm dx}{A(r)^2}
>0.
\tag{8}
\]

因此 \(r_b\) 由 \(M_y/F_n\) 局部唯一确定，随后 \(\Delta q\) 可恢复。该推导正确。退化点也基本明确：\(\Delta q=0\)、侵入长度趋近零、核不为正、SNR 不足、边界方向退化或接触斑塌缩。

### 2.2 实际 MHE 并未继承上述可辨识性

实际 estimator 同时估计

\[
z_k=[r_{b,k},\Delta q_k,a_k,b_{F,k},b_{M,k}]^\top
\tag{9}
\]

而每个时刻只有两维观测。解析证明只覆盖已知 \(a,\chi\) 且无未知偏置的二维对象 \((r_b,\Delta q)\)。允许 \(a,b_F,b_M\) 均作 random walk 后，多个状态轨迹可以由偏置吸收而产生相同观测。加权最小二乘目标本身不构成 set-membership containment theorem。

必须补充：

1. augmented-state 的局部可观测秩或有限窗 injectivity 条件；
2. \(a\) 的独立几何/载荷约束强度；
3. 偏置变化率和噪声的硬界；
4. 从 residual threshold 到
   \(\Pr(z^\star\in\widehat{\mathcal Z})\ge1-\alpha_z\)
   或确定性 \(z^\star\in\widehat{\mathcal Z}\) 的构造；
5. \(\widehat r_b\) 到保守 \(\underline T_{\mathrm{prev}}\) 的误差传播。

### 2.3 \(\Delta q\rightarrow\widehat{\mathcal M}_{c,1}\) 尚无计算定义

Stage 1 声明

\[
\widehat{\mathcal M}_{c,1}
=
\{\mathcal M_c(\eta_{c,1}):
\eta_{c,1}\in\widehat{\mathcal E}_1\},
\tag{10}
\]

但没有定义 \(\widehat{\mathcal E}_1\) 如何由 \(\Delta\hat q\)、当前模型、自然扫描激励与先验得到。当前文字“共同约束”不能用于计算，也不能保证真实 future model 被集合包含。

需要给出明确逆像，例如

\[
\widehat{\mathcal E}_1
=
\left\{
\eta_{c,1}\in\mathcal E_{\mathrm{prior}}:
\exists z\in\widehat{\mathcal Z},\;
\|Y_H-\mathcal H_H(z,\eta_{c,1})\|_{\mathcal V}
\le1
\right\},
\tag{11}
\]

并证明在 Stage 1 的噪声、remainder 和 prior assumptions 下

\[
\eta_{c,1}^{\mathrm{true}}
\in\widehat{\mathcal E}_1.
\tag{12}
\]

如果自然扫描只收紧一个法向有效方向，则其余方向必须作为显式区间保留在式 (11) 中。

### 2.4 Homogeneous cylinder contact model 基本一致

以下部分在声明的局部工作域内数学上相容：

\[
a_H=\sqrt{\frac{4wR^*}{\pi E^*}},
\qquad
p_H(x)=\frac{2w}{\pi a_H}
\sqrt{1-\frac{x^2}{a_H^2}},
\tag{13}
\]

\[
w=h_\infty(\delta_n;\eta_e)+q_v+r_w,
\qquad
\dot q_v=k_v(\delta_n;\eta_v)\dot\delta_n-\frac{q_v}{\tau_v},
\tag{14}
\]

以及持续接触和 sliding branch 内的机器人–接触耦合。这里的 Hertz 部分只是 contact width 与 pressure-shape prior；真实 force–indentation relation 由有限几何黏弹标定给出。这一角色划分正确。

尚需在正式定理中固定 \(h_\infty,k_v\) 的光滑性、Lipschitz 界、单调域和 remainder 界，但这可以通过模型假设完成。

### 2.5 跨边界 contact dynamics 未闭合

当前 Method 使用

\[
\alpha(r_b,a,\chi)
=
\frac{\int_{r_b}^a k(x;a,\chi)\,\mathrm dx}
{\int_{-a}^a k(x;a,\chi)\,\mathrm dx}
\tag{15}
\]

作为 future-region 权重，并提出

\[
F_n=L_{\mathrm{eff}}
\big[
(1-\alpha)h_{\infty,0}
+\alpha h_{\infty,1}
+Q_{v,0}+Q_{v,1}
\big]+d_F,
\tag{16}
\]

\[
\dot Q_{v,j}
=
\omega_jk_{v,j}\dot\delta_n
-\frac{Q_{v,j}}{\tau_{v,j}}
+b_j\dot\omega_j+d_{q,j}.
\tag{17}
\]

式 (15) 中的 \(k\) 是增量压力的观测灵敏度核，不自动等于对完整平衡载荷、黏弹记忆和摩擦载荷都适用的 mixture weight。式 (17) 的 \(Q_{v,j}\) 是局部应力、单位宽度状态还是加权总载荷未定义；\(b_j\dot\omega_j\) 也未从分布 Maxwell state 的 moving-boundary balance 推出。其初值、边界流入状态和材料离开接触域时的 transport 均未定义。

必须从分布状态

\[
q_v=q_v(x,t),\qquad
\eta_c=\eta_c(s_c+x)
\tag{18}
\]

先写出空间积分与 moving-boundary transport，再通过 moment projection 或 quadrature 得到有限维状态，并给出 reduction remainder bound。否则式 (16)–(17) 只是未验证的模型猜测。

### 2.6 Preview phase 与 infinite-horizon HJB 不相容

当前 augmented state 同时包含 \(r_b\)，且

\[
\dot r_b=v_b-v_s.
\tag{19}
\]

在持续向前扫描时，\(r_b\) 不趋向平衡点，而会继续减小。Method 随后定义

\[
\tilde x=x-x_s(\eta)
\tag{20}
\]

和无限时域稳态 HJB，但包含 \(r_b\) 的 \(x_s(\eta)\) 不存在。即便 \(\alpha\) 在 1 处饱和，\(r_b\) 本身仍漂移。若把 \(r_b\) 纳入正定 running cost，积分通常发散；若不惩罚它，原 NEOC 所要求的原点正定 cost、自治平衡与可控性也不再直接成立。

此外，\(r_b\) 同时出现在 state 和 descriptor \(\eta\) 中，造成 state–parameter duplication。

必须在以下两种数学对象中选定一个，不能混写：

1. 有限时域 transition HJB/OCP，以 encounter time 为终点并连接 future homogeneous terminal value；
2. 以饱和 phase \(\zeta\in[0,1]\) 替代 \(r_b\) 的混杂自治系统，严格定义 \(\dot\zeta\)、终端平衡和不可控 phase 对 HJB/NEOC 定理的影响。

第二种对象不能直接引用现有 NEOC 的 complete-controllability 假设，因为 \(\zeta\) 是外生扫描相位，gain-rate action 无法控制它。

### 2.7 Gain-rate action 的语义正确，但最优问题未闭合

把 applied gain

\[
g_a=[k_n,k_t,d_n,d_t]^\top
\tag{21}
\]

作为状态，把

\[
u_g=\dot g_a
\tag{22}
\]

作为 action，解决了 commanded gain、applied gain 与 HJB action 的语义混淆。该修正正确。

但当前又同时保留“理想 gain-rate interface”和一般 renderer

\[
\dot g_a=f_g(g_a,g_c),
\tag{23}
\]

尚未选定最终 plant。两者对应不同 action、不同 \(G_g\) 和不同 HJB。

更重要的是，若 running cost 只惩罚 tracking state 与 \(u_g\)，那么在零 tracking error 时通常存在一族 \(u_g=0\)、任意 \(g_a\) 的平衡。此时 value function 对 augmented state 不正定，\(W^*(\eta)\) 也未必对应唯一的“optimal impedance”。当前 \(x_s(\eta)\) 中的 \(g_s\) 尚未由独立定义得到，不能先假设一个最优 gain 再用 HJB 求它。

必须明确：

- 是优化 transition gain trajectory，还是求 homogeneous environment 的 feedback-gain law；
- gain 本身通过何种 performance term、disturbance model 或 terminal condition获得唯一决策意义；
- baseline low-level stabilizer 与 gain-rate optimization 的 stabilizability/admissibility；
- non-zero desired force equilibrium 的存在与局部唯一性。

### 2.8 Constrained HJB 尚未定义为可解对象

当前写出

\[
0=\min_{u_g\in\mathcal U_g(x,\eta)}
\{\ell+\nabla V^\top(f_0+G_gu_g)\},
\tag{24}
\]

但 \(\mathcal U_g\) 同时承载 gain bounds、contact maintenance、force、slip、actuator 和 passivity constraints。对 state constraints，式 (24) 需要 viability domain 与边界条件；对 contact/slip switches，需要 hybrid HJB 或明确排除这些 mode；对 energy tank，需要把 tank state放入 plant，而不是只在求解后检查。

若先求 unconstrained HJB、再由 safety/passivity filter 修改 action，最终 applied policy 不再是式 (24) 的 optimal policy，随后 \(W^*(\eta)\) 和 ANEOC sensitivity 也不再对应真实执行闭环。

### 2.9 Single critic root 尚未闭合

当前只有抽象式

\[
V_N(x;\eta)=W(\eta)^\top\phi(x),
\qquad
F_{\mathrm{HJB}}(W,\eta)=0.
\tag{25}
\]

仍缺少：

1. 固定 state domain、basis 和 normalization；
2. Galerkin/collocation test functions 与积分/采样规则；
3. \(F_{\mathrm{HJB}}\) 的逐项展开；
4. projected root 的存在和局部唯一性；
5. \(F_W\) 非奇异条件；
6. value-gradient error 到 policy error 的界；
7. critic residual 到 closed-loop stability/constraint error 的界；
8. online trajectory learning还是 model-based collocation 的选择。

一般 neural-network weights 存在 permutation/scale symmetry，不能默认 \(F_W\) 可逆。若使用 fixed linear-in-weight basis，应明确 single critic 是有限基 value approximator，并避免在患者上依赖额外 exploration/PE。

### 2.10 ANEOC 局部公式正确，但还没有形成可执行定理

若确有光滑且局部唯一的 root，

\[
F(W^*(\eta),\eta)=0,
\tag{26}
\]

则

\[
\frac{\mathrm dW^*}{\mathrm d\eta}
=-F_W^{-1}F_\eta
\tag{27}
\]

和

\[
W_{\mathrm{pred}}
=W^*(\eta_0)
-F_W^{-1}F_\eta\Delta\eta
\tag{28}
\]

成立，理想 truncation error 为 \(O(\|\Delta\eta\|^2)\)。这部分微分本身正确，也与 parameter-dependent closed-loop NEOC 的基本机制一致。[Rai et al.](https://arxiv.org/abs/2312.04733)明确要求光滑动力学、正定 cost、admissible stabilizing law 和局部 regularity，并用 homotopy 处理较大参数变化。

但当前 route 还缺：

- observation/set error、current critic error、sensitivity error、Taylor remainder 与 corrector residual 的联合界；
- descriptor set 到 policy set 的可计算传播；
- trust-region/homotopy step rule；
- active-set、contact-mode 或 friction-mode 切换后的 root definition；
- future transition HJB 与 current homogeneous HJB 是否属于同维、同 basis、同 root family；
- sensitivity 何时计算和保存。

而且从旧 root 对 future problem 做一次 Newton 更新满足

\[
\begin{aligned}
W_{\mathrm{N},1}
&=
W_0
-F_W(W_0,\eta_1)^{-1}F(W_0,\eta_1)\\
&=
W_0-F_W^{-1}F_\eta\Delta\eta
+O(\|\Delta\eta\|^2),
\end{aligned}
\tag{29}
\]

与式 (28) 一阶相同。因此，若 \(F_W^{-1}F_\eta\) 在事件后才计算，ANEOC 并未形成独立的 deadline advantage。

### 2.11 Pre-transition execution theorem 尚未闭合

当前 deadline inequality

\[
\tau_{\mathrm{filter}}+\tau_{\mathrm{preview}}
+\tau_{\mathrm{transport}}+\tau_{\mathrm{correct}}
+\tau_{\mathrm{comm}}+\tau_{\mathrm{render}}
\le\underline T_{\mathrm{prev}}
\tag{30}
\]

只是必要的时间条件，不保证 gain 可达。完整条件至少需要

\[
\mathcal R_g
\left(
g_a(t_p),
\underline T_{\mathrm{prev}}-\tau_{\mathrm{compute}};
\mathcal U_g
\right)
\cap
\mathcal G_{\mathrm{cross}}
\left(\widehat{\mathcal E}_1\right)
\neq\varnothing,
\tag{31}
\]

其中 \(\mathcal R_g\) 是考虑 rate、actuator、communication 和 renderer error 后的 applied-gain reachable set。

当前 passivity inequality 只统计

\[
\frac12 e^\top\dot K e,
\tag{32}
\]

尚未定义完整 energy-tank dynamics、tank switching、damping harvest、feedforward power、旋转 contact frame 产生的 \(\dot R_c\) 项，以及 renderer tracking error。也没有从式 (31)、critic error 与 model-set containment 推出 crossing constraint satisfaction。

### 第 2 项结论

**NOT CLOSED。**

理想 step-transition wrench-moment 可辨识性和 homogeneous cylinder-contact backbone 已建立；从 practical estimator 开始，经 future model set、transition dynamics、HJB、critic、ANEOC 到 applied impedance 的其余主干尚未形成完整定理链。

---

## 3. 前后接口与数学对象一致性

| 接口 | 上游实际输出 | 下游要求 | 当前问题 |
|---|---|---|---|
| wrench \(\rightarrow\) preview | \(\hat r_b,\Delta\hat q,\hat a\) 与不确定性 | future contact-model set | 缺式 (11) 的逆像与 containment |
| preview \(\rightarrow\) transition model | boundary/model set | finite-dimensional Markov dynamics | mixture/transport 未推导 |
| transition model \(\rightarrow\eta\) | \(\eta_c\) 的集合 | smooth policy parameter | \(\eta\) 含 \(r_b\) 和 set radius，既是 state 又是 parameter，且可能 nonsmooth |
| \(\eta\rightarrow\) HJB | uncertain future model + preview phase | fixed-dimensional autonomous plant and equilibrium | phase drift、model set 与 point-HJB 不一致 |
| HJB \(\rightarrow\) critic | constrained/hybrid value problem | differentiable finite root \(F(W,\eta)=0\) | projection、boundary condition、uniqueness 未定义 |
| critic \(\rightarrow\) ANEOC | approximate root | regular solution branch | \(F_W\) nonsingularity 与 total error bound 缺失 |
| ANEOC \(\rightarrow\) impedance | future critic/policy | admissible gain trajectory | policy只给 gain rate；target/reachability 未闭合 |
| impedance \(\rightarrow\) robot | commanded/desired gains | measured applied gains | renderer 二选一、tracking error 未固定 |
| robot \(\rightarrow\) safety claim | contact trajectory | force/slip/contact/passivity guarantee | 缺 end-to-end robust theorem |

另外有三个明确的对象冲突：

1. \(\Delta q\) 是当前自然激励下的 observable contrast，不是完整 \(\eta_c\)，更不是直接的 HJB parameter；
2. \(r_b\) 不能同时作为 HJB state 和 \(\eta\) 的固定分量；
3. force 是 contact model 的 output，\(g_a\) 是 applied state，\(u_g\) 或 \(g_c\) 才是 action。当前 Method 已修正第一层语义，但 renderer action 尚未最终选择。

### 第 3 项结论

**NOT CLOSED。**

接口链方向自然，但类型还未闭合；尤其是 model set \(\rightarrow\) smooth descriptor \(\rightarrow\) fixed HJB root 这一段目前不可直接执行。

---

## 4. Method fit / necessity

本项只判断各机制是否由目标 failure 自然推出，以及是否存在承担同一任务或同一模块职责的更简单充分方案。用途或场景不同的方法不作为否决依据；找不到有效对比时不制造对比。

### 4.1 Spatial preview 与 wrench moment

在修正后的 failure 中，需要从已经变化的 mixed footprint wrench 中分离：

\[
\text{transition location/remaining time}
\quad\text{与}\quad
\text{future response contrast}.
\]

单时刻 \(F_n\) 只有乘积信息，加入 \(M_y\) 后在理想 step family 内可分离二者，因此 zero/first moment 的选择与 failure 直接对应。这里不需要把一般视觉 preview、车辆 preview 或压力阵列当作同任务竞争方法。

尚未闭合的是：使用较长 force–pose history、直接使用完整 6-D wrench 的高带宽反馈、或 end-to-end causal policy 是否能在相同实际 preview window 内达到相同 crossing-safe action。该问题只需在相同传感、transition、latency 和执行限制下比较。当前 Method 应把 moment estimator 的必要性限定为“提供显式位置—幅值分离与可验证时间证书”，而不是声称其他反馈完全没有前缘信息。

### 4.2 Cylinder–plane contact model

该模块必须同时服务 contact width、pressure moment、normal viscoelasticity、tangential sliding 和 robot coupling。原 NNBO 的一维 Hunt–Crossley 模型没有有限圆柱接触斑和空间矩，球–平面 Hertz 也不匹配目标探头形态。因此以 cylinder line contact 作为几何 prior、以局部有限组织黏弹映射作为实际 constitutive model，是由目标传感与控制对象共同推出的自然选择。

但跨边界模型尚未从分布接触约化，故该模块当前仍不能向 HJB 提供可信的 transition dynamics。

### 4.3 Variable optimal impedance

环境变化并不自动证明必须求 optimal variable impedance。最简单的同任务方案是一个 fixed robust gain \(g_0\)。必须定义

\[
J^*(\eta)=\inf_{\pi\in\Pi_{\mathrm{safe}}(\eta)}J(\pi;\eta)
\tag{33}
\]

和固定增益 regret

\[
\mathcal R(g_0,\eta)=J(g_0;\eta)-J^*(\eta).
\tag{34}
\]

若存在

\[
g_0\in\bigcap_{\eta\in\mathcal E}\mathcal G_{\mathrm{safe}}(\eta)
\quad\text{且}\quad
\sup_{\eta\in\mathcal E}\mathcal R(g_0,\eta)
\le\varepsilon_J,
\tag{35}
\]

则为该工作域引入 future optimal impedance 没有必要。当前 Method 尚未定义 \(\varepsilon_J\)、\(\mathcal E\) 和 incompatibility witness，也未证明式 (35) 失败。

这不是要求用不相干方法做泛化比较，而是检查同一任务中固定阻抗是否已经充分。

### 4.4 HJB 与 single critic

若任务确实要求一个在状态域上重复执行的 nonlinear feedback law，而不是单条 open-loop trajectory，那么 value-function/HJB 表示具有直接作用；single critic 在 Hamiltonian minimization 可由 value gradient恢复时也比额外 actor 更简洁。原 NNBO 支持“HJB value approximation \(\rightarrow\) single critic \(\rightarrow\) policy recovery”这一局部求解结构，但其 Hunt–Crossley force action 不能迁移。[Kong et al., 2025](https://doi.org/10.1109/TSMC.2025.3579017)

然而，本研究的真正求解对象尚未在“finite-horizon transition problem”和“autonomous phase-augmented infinite-horizon HJB”之间确定。因此无法判断 HJB 是否是最自然的数学对象。一个使用完全相同 transition dynamics、cost、gain-rate action 和 constraints 的 finite-horizon direct solver 承担相同求解职责，是有效的内部替代；它不需要来自同一应用论文才能作为求解器对照。

已有 variable-stiffness MPC 工作证明，有限时域内直接优化 stiffness/damping 并处理其执行稳定性是可实现的局部求解机制，但其应用场景不等同于本研究，故只用于说明 solver capability，而不用于否定整体路线。[Thelenberg and Ott, 2024](https://www.paperhost.org/proceedings/controls/ECC24/files/0838.pdf)

HJB 的选择必须由最终问题对象决定。目前该对象未闭合，因此 HJB necessity 未建立。

### 4.5 ANEOC predictor–corrector

当目标确实是更新同一个 parameterized HJB/Galerkin root 时，合理对比是：

- cold re-solve；
- old-root warm-start Newton；
- 预维护 Jacobian/sensitivity 的 predictor–corrector；
- 大变化下的 homotopy continuation。

这四者承担完全相同的 root-update 任务。式 (29) 表明 ANEOC predictor 与一次 warm-start Newton 在一阶上等价。只有当 \(F_W\) 的 factorization 或 \(C_\eta\) 在 preview event 前已经维护、能够复用，并减少 future-Jacobian factorization/corrector 次数时，ANEOC 才形成额外计算价值。Sensitivity warm-start 在受约束 NMPC 中也有成熟的同模块机制，说明 active constraints 不能只靠 fixed-active-set IFT 处理。[Liao-McPherson et al., 2020](https://doi.org/10.1109/TAC.2019.2954359)

当前没有给出预维护算法、operation count、corrector convergence radius 或与 warm-start Newton 的同 root 比较，因此 ANEOC necessity 未建立。

### 4.6 Rate、passivity 与 deadline gates

这些模块不是为了增加方法数量。只要 stiffness 随时间变化，

\[
\dot S_K
\supset
\frac12e^\top\dot K e
\tag{36}
\]

就可能产生正能量；只要低层 gain reconstruction 有限速率，完成计算也不等于 applied impedance 已经到位。因此 passivity state、gain reachable set 和 deadline gate 是 pre-transition execution claim 的必要组成。

当前问题不是这些模块是否需要，而是它们尚未被并入同一个 plant、HJB 和 end-to-end theorem。

### 第 4 项结论

**NOT CLOSED。**

spatial moments、cylinder-contact backbone、applied-gain state 和 execution gates 与修正后的 failure 具有直接因果对应。尚未证明 variable optimal impedance 相对 fixed robust gain 的必要性，尚未确定 HJB 相对同对象 finite-horizon solver 的方法选择，也尚未证明 ANEOC 相对同 root warm-start Newton 的净作用。这三个缺口必须在最终数学对象固定后解决。

---

## 非实验性闭合清单

以下 12 项必须在进入最终 Method 验收前完成：

1. 重写 causal failure：从“task center 到达前没有信息”改为“mixed footprint 信息不能直接确定 future model、remaining time 与 future action”。
2. 给出 force-only indistinguishability + incompatible crossing action 的问题级命题。
3. 证明 practical MHE augmented state 的有限窗可辨识/集合 containment，或收窄其输出声明。
4. 明确定义 \(\Delta q\rightarrow\widehat{\mathcal E}_1\) 的 inverse feasible set。
5. 从 distributed cylinder–tissue viscoelastic field 推导 transition reduced state、boundary transport 和 remainder bound。
6. 在 finite-horizon transition problem 与 autonomous phase HJB 之间选定唯一数学对象。
7. 从 \(\eta\) 中移除 state duplication，并证明 descriptor 的 policy sufficiency 或保留不可约 model set。
8. 选定真实 renderer action，定义 non-zero-force equilibrium、admissible baseline 和 gain decision uniqueness。
9. 完整定义 constrained/hybrid HJB 的 domain、boundary condition、safe set 和 passivity state。
10. 展开 \(F_{\mathrm{HJB}}(W,\eta)\)，证明 projected root regularity，并给出 critic/value-gradient/policy error bounds。
11. 给出 ANEOC total error、trust region、homotopy、active-set/mode-switch 与 warm-start complexity 定理。
12. 给出 applied-gain reachability + energy tank + uncertainty propagation + crossing constraint 的 end-to-end theorem。

## 仅由实验回答、不得伪装成数学缺口的项目

- 实际 \(F_n/M_y\) SNR 是否足以在所需距离检测 transition；
- 两分支或更高阶 transition model 的标定误差；
- P95/P99 estimation、solver、communication 和 renderer latency；
- fixed robust、matched reactive controller、warm-start solver 与完整方法的实际 first-hit force/contact/slip表现；
- 图像质量、患者舒适性和医生可接受性；
- 一个 Maxwell branch 是否覆盖目标组织与速度范围。

上述实验项可以暂不完成，但在实验前必须先完成前一节的 12 个数学闭合项。

## 最终结论

# NOT CLOSED

本轮不否定研究方向，也不改做其他问题。当前最先要修正的不是 HJB 细节，而是第 1–5 项：准确重写 causal failure，并把 practical preview observation 严格映射到可包含真实模型的 future contact set，再从分布接触推导跨边界有限维动力学。只有这一入口闭合后，\(\eta\)、HJB、critic、ANEOC 和 pre-transition execution 才有确定的数学对象。

