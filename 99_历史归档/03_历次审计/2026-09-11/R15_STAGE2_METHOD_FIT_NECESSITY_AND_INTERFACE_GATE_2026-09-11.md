# Stage 2：Method-fit / Necessity 与下游接口门

日期：2026-09-11  
状态：**REVISE-AND-CONTINUE；允许进入完整动力学推导，但 NNBO/HJB 与 ANEOC 尚未获得无条件 necessity PASS**  
上游基准：R15_STAGE1_PROBLEM_TO_SPATIAL_PREVIEW_AND_CYLINDER_CONTACT_2026-09-11.md  
目标：固定 Stage 1 已确认的研究问题和物理模型，判断后续 optimal-control / single-critic / ANEOC 路线是否自然、必要、可执行，并在正式推导前修正接口语义。

---

## 1. 固定参照系

### 1.1 Accepted failure

连续机器人超声跨越局部组织转变时，task center 已进入新区域后，current-point sensing、参数更新、优化求解、通信和阻抗重构才依次发生，因而即使稳态控制最终正确，仍可能先经历 first-hit force peak、contact loss、slip、path deviation 或 image interruption。该 failure 是 **causal adaptation lag**，不是一般的稳态 force-tracking error。

### 1.2 Required mechanism changes

本阶段固定四个 RMC：

1. **RMC-I：改变信息时序。** 将只在 task-center encounter 后可用的接触变化信息，提前到有限接触斑前缘已经接触、task center 尚未接触的时间窗。
2. **RMC-II：改变 future contact information 的控制角色。** future descriptor 不能只进入监测或预测，而必须改变实际执行的反馈力学。
3. **RMC-III：改变最优解更新时间。** 若未来接触会改变 optimal feedback law，必须在有限 preview time 内获得 future solution，而不是 encounter 后重新训练或完整重求。
4. **RMC-IV：闭合 applied action。** 未来策略必须在 gain-rate、actuator、communication、contact 与 energy/passivity 条件下，于 task-center encounter 前成为 applied impedance。

### 1.3 Target constraints

- 不使用患者特异的离线逐点触诊地图；
- 仅依赖目标系统中可提供的 6-axis F/T、robot pose/velocity、表面几何及短历史；
- 持续压缩接触、已知扫描方向、主定理先限定于 slip branch；
- cylinder–plane nominal geometry，有限厚度/有限长度/黏弹性通过已标定局部模型和有界 remainder 表示；
- 不把一次或一段弱激励的 wrench 解释成完整材料参数；
- 输出必须通过真实 low-level controller 成为 applied normal/tangential impedance；
- end-to-end deadline 和 passivity/energy 不能由事后 safety narration 代替。

---

## 2. 已有机制边界

以下原子能力已被 prior 覆盖，不能单独承担 novelty：

- online/direct stiffness–damping optimization；
- finite-horizon variable compliance planning；
- uncertainty-aware predictive impedance；
- online environment estimation + VIC/MPVIC；
- contact-aware robot dynamics and MPC；
- HJB/critic-based nonlinear optimal interaction；
- parameter sensitivity / neighboring-extremal continuation；
- energy tank、passivity gate、constraint handling 和 solver acceleration。

与本阶段直接相关的原始证据包括：

- Kong et al., NNBO：Hunt–Crossley 模型被改写为 control-affine dynamics，并明确将 interaction force \(f\) 作为 HJB control；critic 逼近 value function，最后利用 HC-specific algebraic structure 把 force policy 解释为 nonlinear impedance，再通过 Newton-downhill 求 position reference。DOI: https://doi.org/10.1109/TSMC.2025.3579017
- Rai, Mou, and Anderson：NEOC 对已知参数扰动下的 closed-loop optimal law 求一阶变化；Galerkin 数值版本通过线性方程求 coefficient sensitivity，并用 homotopy 处理较大扰动。DOI: https://doi.org/10.1016/j.automatica.2025.112693；公开版本：https://arxiv.org/abs/2312.04733
- Thelenberg and Ott：已能把 time-varying stiffness/damping 置于 finite-horizon MPC 中，实际执行第一个参数解，并讨论 \(\frac12\tilde x^\top\dot K\tilde x\) 的非被动风险。DOI: https://doi.org/10.23919/ECC64448.2024.10591005
- Pollayil et al.：已能以较集中、解析/数值结合的优化直接规划可执行 stiffness/damping；因此 optimal impedance 本身不要求 HJB。DOI: https://doi.org/10.1109/TRO.2022.3216078
- Liao-McPherson et al.：对 constrained NMPC 的参数化 OCP 已有 sensitivity-based warm start，而且使用 semiderivative 处理 polyhedral constraints，不要求活动约束状态完全不变。DOI: https://doi.org/10.1109/TAC.2019.2954359

结论：当前论文若成立，Scientific Delta 只能位于 **finite-patch future-contact observability → future-contact-conditioned optimal feedback mechanics → preview-time solution transport → pre-transition physical execution** 这一有向关系；HJB、critic、ANEOC 和 passivity 均只能作为实现该关系的必要机制或支撑机制。

---

## 3. 决定性 Method choice A：是否需要 variable/optimal impedance

令 \(\eta\in\mathcal E\) 表示未来局部接触描述，\(g\) 表示可执行阻抗参数，\(\mathcal G_{\rm safe}(\eta)\) 表示在接触、力、滑移、执行器和能量约束下的安全阻抗集合，\(J(g;\eta)\) 表示同一任务目标下的性能。定义

\[
J^*(\eta)=\min_{g\in\mathcal G_{\rm safe}(\eta)}J(g;\eta),
\qquad
\mathcal R(g;\eta)=J(g;\eta)-J^*(\eta).
\]

固定鲁棒阻抗足够的反事实为：存在 \(g_0\)，使

\[
g_0\in\bigcap_{\eta\in\mathcal E}\mathcal G_{\rm safe}(\eta),
\qquad
\sup_{\eta\in\mathcal E}\mathcal R(g_0;\eta)\le\varepsilon_J.
\]

若该条件在目标 operating envelope 内成立，future optimal impedance 没有 necessity。因而必须在实验前定义、实验中验证至少一个判别条件：

\[
\min_{g\in\cap_{\eta\in\mathcal E}\mathcal G_{\rm safe}(\eta)}
\sup_{\eta\in\mathcal E}\mathcal R(g;\eta)>\varepsilon_J,
\]

或者不同 future contacts 的安全集合不存在共同可接受元素。它对应“两种未来接触条件需要不可兼容或显著不同的 feedback mechanics”。

**判定 A：CONDITIONAL PASS。** Variable/optimal impedance 由 failure 自然推出，但只有上述 non-common-near-optimum witness 成立时才不可被 fixed robust impedance 支配。

---

## 4. 决定性 Method choice B：为什么不能直接使用 reactive force control / VIC

对任意只依赖 task-center 当前历史 \(\mathcal H_t^c\) 的 causal controller，若两个 future transitions \(\eta_1,\eta_2\) 在 encounter 前产生相同的 \(\mathcal H_t^c\)，则控制器在该时刻必须发布相同动作。若两种未来条件要求不同的 pre-transition feasible gain sets，controller 至少在一种条件下出现 transition loss。

Stage 1 的有限接触斑 moment observer 提供额外信息 \(\mathcal H_t^{FM}\)。其必要性应通过 matched comparison 建立：

\[
\mathcal H_t^c(\eta_1)\approx\mathcal H_t^c(\eta_2),
\qquad
\mathcal H_t^{FM}(\eta_1)\not\approx\mathcal H_t^{FM}(\eta_2),
\]

且两个未来条件的 target gain regions 不同。如果 matched high-bandwidth explicit force control 在 encounter 后的反应仍可消除全部 first-hit loss，preview 的任务必要性失败。

**判定 B：CONDITIONAL PASS。** Previewed adaptation 与 reactive control 的 intervention position 不同，但 superiority 不能由“更早”自动推出，必须由 first-hit transition outcome 和实际时延共同证明。

---

## 5. 决定性 Method choice C：是否必须使用 HJB/NNBO

### 5.1 主要竞争机制

| Mechanism class | 能实现的共同 RMC | 主要边界 |
|---|---|---|
| 固定鲁棒阻抗 | 安全覆盖一组环境 | 若同时近最优，则否决 variable-optimal route |
| preview-conditioned lookup / static law | future descriptor 直接映射到 gains | 适合低维、准静态、约束可解析情形 |
| preview-conditioned current-step QP | 直接选择物理 gains 并处理局部约束 | 不能自动处理黏弹记忆和跨时段状态后果 |
| finite-horizon VIC/NMPC | 处理 gain dynamics、hard constraints、time-to-boundary | 每次需解 OCP，但与预过渡任务最直接 |
| infinite/long-horizon HJB + critic | 得到覆盖一片状态域的 nonlinear feedback law，在线 policy evaluation 便宜 | HJB/critic 构造困难；hard constraints、hybrid modes 和 time-to-boundary 会增加维度 |

### 5.2 HJB 的必要条件

HJB 只有同时满足以下条件才比更简单机制更自然：

1. optimal gain decision 对 contact memory、robot state、force/path error 和 normal–tangential coupling 显著依赖，不能由一个静态 map 或小 QP 近似；
2. 任务需要一个可在状态扰动下持续反馈的 policy，而不是只需要 boundary 时刻的单个 gain target；
3. 同一 parameterized problem 会被反复查询，使 critic 的预计算/延续成本能被摊销；
4. finite-horizon direct solver 在相同模型、约束和硬件下不能满足 preview deadline，或其解/梯度质量明显低于 critic policy；
5. critic approximation 的维度、误差和 regularity 可被验证，而不是把 HJB 求解困难转移成不可审计的 NN 误差。

### 5.3 当前判定

**HJB/NNBO：INDETERMINATE，保留为 primary candidate 进入下一阶段推导。**

理由不是 HJB 更先进，而是新的 cylinder-contact system 含黏弹内部状态、robot/contact coupling 和状态依赖约束，可能需要 nonlinear feedback law；但 finite-horizon gain-planning MPC 在 time-to-boundary 与 hard constraints 上更直接，是必须保留的 strongest simpler rival。在建立最小 Markov dynamics 后，如果一个低维 QP/NMPC 已经同样实时并达到相同 transition capability，则 HJB 必须删除。

---

## 6. 结构性修正 1：HJB action 不能与 applied impedance 混同

令

\[
g_a=[k_n,k_t,d_n,d_t]^\top
\]

表示实际施加的阻抗参数。真实系统存在 command filtering、rate limit、communication delay 和 inner-loop bandwidth，因此一般只能写成

\[
\dot g_a=f_g(g_a,g_c),
\]

或在低层直接支持 slew command 时写成

\[
\dot g_a=r_g,
\qquad r_g\in\mathcal R_g(g_a,E_{\rm tank}).
\]

接触坐标中的命令 wrench 为

\[
w_{c,\rm cmd}=w_{ff}-K(g_a)e_x-D(g_a)e_v,
\]

再通过 Jacobian 与 torque loop 作用到 robot/contact dynamics。

因此：

- \(g_a\) 是 physical state / applied feedback mechanics；
- \(g_c\) 或 \(r_g\) 才是上层 optimizer 能发布的 action；
- 若使用 \(r_g\)，增广系统可自然形成对 action affine 的结构；
- 只有当实测 renderer bandwidth 远高于相关 contact dynamics，才允许把 \(g_c\approx g_a\) 作为受控近似，并必须给误差界。

原候选“直接把 \(K,D\) 当作瞬时 HJB action”与本研究自己的 gain-rate failure 冲突。后续默认采用 **applied gains as state，gain rate/command as action**，除非实机辨识证明准静态近似足够。

---

## 7. 结构性修正 2：future HJB optimum 与 pre-transition trajectory 不能混同

仅在均匀 future environment \(\eta_f\) 上求一个 stationary optimal law，并在边界前直接启用，未必对当前区域安全或最优。必须选择下列两个语义之一：

### Route H1：future target + feasibility renderer

HJB 只给出 future-environment target policy / target gain region；pre-transition 阶段只求一个满足 current-side safety、gain-rate 和 energy 条件的可达轨迹。此时论文能声称“未来最优阻抗在 encounter 前完成预配置”，不能声称整个过渡轨迹全局最优。

### Route H2：preview-augmented optimal control

为避免与 Stage 1 中的力矩–合力比 \(\rho(r_b)=\widetilde M_y/\widetilde F_n\) 混淆，将边界 phase 记为 \(\sigma_b=r_b\)，并纳入状态：

\[
\dot\sigma_b=v_b-v_s,
\]

并用 Stage 1 的有限接触模型描述 \(\eta_c(\sigma_b;\eta_0,\eta_f)\)。HJB state 同时包含 robot/contact state、applied gain state 和 preview phase；future descriptor \(\eta_f\) 是 parameter。这样得到的 policy 才直接优化 pre-transition gain motion 和 post-transition interaction。

**当前优先推导 H2，H1 作为最小 fallback。** H2 能使 spatial preview、gain-rate action 和 future-contact policy 成为同一个 optimal-control problem；但若维度或 hybrid complexity 使 critic/HJB 不可验证，则回退 H1，而不是伪称 H2 已闭合。

---

## 8. 结构性修正 3：Stage 1 的 \(\Delta q\) 不是 HJB descriptor \(\eta\)

Stage 1 的 \(\Delta q\) 是自然扫描工作点下可观测的 effective response contrast。改变 impedance 会改变 indentation、velocity 和状态轨迹，因此同一 \(\Delta q\) 未必唯一决定 off-nominal contact dynamics。

设完整局部接触参数为 \(\theta_c\)，observer 输出一致集合 \(\Theta_f\)。下游 descriptor \(\eta=T(\theta_c)\) 只有在下列意义下才是 policy-sufficient：对于所有满足 \(T(\theta_1)=T(\theta_2)\) 的模型，

\[
\sup_{x\in\Omega}
\|\pi^*(x;\theta_1)-\pi^*(x;\theta_2)\|
\le \varepsilon_\pi,
\]

并且两者的 hard-constraint margin 差异不超过已分配预算。否则该降维会把需要不同动作的模型错误合并。

下一阶段应从 cylinder–viscoelastic transfer relation 和 reachable action directions 反推 \(\eta\)。当前可预期但尚未冻结的候选是：

- 法向平衡切线响应；
- 与控制带宽相关的黏弹增量和 relaxation time；
- slip branch 的 \(\mu_k,c_t\) 或其保守区间；
- geometry/robot inertia 等已知量作为 state/scheduling variables，而不是混入待辨识 material descriptor。

未被自然扫描激励辨识的方向必须留在 nuisance uncertainty set 中。不能为了让 ANEOC 接受点参数而伪造 full-ID。

---

## 9. 决定性 Method choice D：single critic 是否成立

若采用

\[
V_N(x;\eta)=W(\eta)^\top\phi(x)
\]

则 single critic 只是 HJB 数值表示，不是独立科学贡献。它可保留的必要条件为：

1. basis \(\phi\) 固定且对目标状态域足够；
2. projected/Galerkin root \(F_N(W,\eta)=0\) 与 exact HJB 明确区分；
3. basis 无冗余或采用明确 gauge/regularization，使 coefficient root 局部唯一；
4. \(F_W\) 条件良好，critic error 可传播到 policy/constraint margin；
5. 不依赖患者扫描过程中的 unsafe exploration noise 或 PE 才能获得 future critic。

一般多层 NN 权重存在 permutation/scaling non-uniqueness，可能导致 \(F_W\) 奇异；因此后续 ANEOC 主定理不能对任意“神经网络权重”直接使用隐函数定理。优先使用 **fixed features + linear-in-weight critic / nonredundant Galerkin coefficients**。原 NNBO 的 online PE learning 不能原样迁移到患者扫描。

**判定 D：CONDITIONAL PASS。** 若 HJB 保留且有限基表示可通过 residual、conditioning 和 policy-error 检查，single critic 是合理的最小数值结构；否则使用直接 collocation/OCP，不为保留 NNBO 名称强行使用 critic。

---

## 10. 决定性 Method choice E：ANEOC 是否必要

设 future critic root 为

\[
F_N(W^*(\eta),\eta)=0,
\qquad
C_\eta=-F_W^{-1}F_\eta.
\]

ANEOC predictor 为

\[
W_{\rm pred}=W^*(\eta)+C_\eta\Delta\eta.
\]

但是从旧根对新参数问题做一次 Newton step：

\[
W_1=W^*(\eta)-F_W(W^*(\eta),\eta+\Delta\eta)^{-1}
F_N(W^*(\eta),\eta+\Delta\eta).
\]

Taylor 展开给出

\[
F_N(W^*(\eta),\eta+\Delta\eta)
=F_\eta\Delta\eta+O(\|\Delta\eta\|^2),
\]

从而

\[
W_1=W^*(\eta)-F_W^{-1}F_\eta\Delta\eta
+O(\|\Delta\eta\|^2)
=W_{\rm pred}+O(\|\Delta\eta\|^2).
\]

因此 **ANEOC predictor 与 warm-start future solve 的第一 Newton 步在局部一阶上等价**。若每次 preview 后才在线计算 \(F_W^{-1}F_\eta\)，再运行同样的 Newton corrector，则不能仅凭 \(O(\|\Delta\eta\|^2)\) predictor error 宣称更快。

ANEOC 的非支配作用只能来自：

1. 在当前 regular branch 上后台预计算或持续维护 \(C_\eta\)，使 preview arrival 后 predictor 只需 matrix–vector product；
2. 一次 sensitivity 可廉价传播多个 future descriptor vertices/candidates；
3. predictor 显著减少 residual evaluations、factorizations、line searches 或 corrector iterations；
4. 在相同精度、admissibility 与硬件上，P95/P99 runtime 明显优于 stale、warm-start Newton 和 full solve；
5. branch conditioning、trust region 和 re-anchor semantics 可在线检查。

若存在 active constraints，应对完整 residual/KKT variable \(z=[W,\lambda,\ldots]\) 求

\[
R_z z_\eta+R_\eta=0,
\]

而不是只对 \(W\) 求普通导数。固定 active set 内可使用 classical IFT；active-set switch、contact-mode switch 或 nonsmooth saturation 处需要 semiderivative/generalized equation 或直接 re-solve。

**判定 E：INDETERMINATE，条件性保留。** ANEOC 是合理 source mechanism，但其 necessity 必须由“预维护 sensitivity 后的真实 deadline advantage”证明。若 warm-start Newton/full solve 已同样快，ANEOC 删除；若 constrained NMPC sensitivity 更自然，则 HJB-ANEOC 路线也必须重新比较。

---

## 11. 最小自然架构（当前 provisional selection）

在尚未完成下一阶段动力学推导前，最小且不自相矛盾的候选架构为：

1. **Finite-patch preview observer**：输出 \(r_b,T_{\rm prev}\) 和 contact-model consistency set，不输出虚假 full-ID。
2. **Preview-augmented cylinder-contact dynamics**：状态含 robot/contact memory、applied gains 和 boundary phase；future descriptor 仅参数化被未来区域改变的动力学方向。
3. **Executable-gain optimal control**：applied gains 为状态，gain command/rate 为 action；单一 cost 同时表达 force、path、slip/contact 与 gain-motion trade-off。
4. **Conditional single critic**：只对固定特征、局部唯一的 projected HJB root 使用 \(W^*(\eta)\)。
5. **Conditional ANEOC**：当前 branch 上预维护 sensitivity，preview 到来后 predictor–corrector；不把在线求 sensitivity 的成本隐藏掉。
6. **Execution gate**：统一使用 preview lower bound、actual gain dynamics、actuator、energy/passivity、model/critic/corrector error；失败时减速、保持 verified incumbent 或退出接触。

该架构不是模块拼接。其一条因果链是：

\[
\text{future transition becomes observable}
\rightarrow
\text{preview phase and future dynamics enter one control problem}
\rightarrow
\text{rate-limited feedback mechanics are optimized as physical states}
\rightarrow
\text{the parameterized solution is transported before encounter}
\rightarrow
\text{only verified applied gains are released}.
\]

---

## 12. 需要在后续阶段执行的 kill tests

### K1 — Variable-impedance necessity

若 fixed robust impedance 对全部目标 transitions 同时满足 hard constraints，且 worst-case regret 不超过预定 \(\varepsilon_J\)，删除 optimal variable-impedance layer。

### K2 — Preview necessity

若 matched explicit force controller 或 force-history estimator 在相同 sensor、bandwidth 和 delay 下消除全部 first-hit loss，收窄或否决 preview-control claim。

### K3 — HJB necessity

若 preview-conditioned static map、QP 或 finite-horizon NMPC 以更弱假设、更低复杂度达到相同状态反馈、约束和 deadline capability，删除 HJB/critic。

### K4 — Critic viability

若 projected root 非唯一、\(F_W\) 近奇异、basis error 无法压入 safety margin，或 critic state dimension 不可计算，删除 single critic。

### K5 — ANEOC necessity

若 ANEOC 与 warm-start/full solve 的 P95/P99 runtime、residual 和成功率无显著差异，删除 ANEOC；不能用平均 iteration 数替代 deadline evidence。

### K6 — Descriptor sufficiency

若 observer-equivalent future models 会产生不相容 optimal actions 或 safety margins，则当前 \(\eta\) 不充分；必须扩充可辨识 descriptor、保留 robust set，或缩小 claim/operating envelope。

### K7 — Pre-transition feasibility

若在临床有意义的最低扫描速度下仍有

\[
T_{\rm prev,lower}<
\tau_{\rm estimate}+\tau_{\rm solve}+\tau_{\rm communicate}+\tau_{\rm render}
\]

或 gain/energy reachable set 不含 future target，路线不具备 pre-transition capability。

---

## 13. Devil's-advocate checkpoint

### Verdict：REVISE-AND-CONTINUE

没有发现必须推翻 Stage 1 的 critical flaw，但发现三个若不修正会使整条 Method 失效的 major issues：

1. **Action–execution equivocation**：把 commanded gains、applied gains 和 HJB action 当成同一对象。已通过“applied gains as state，gain rate/command as action”修正。
2. **Future optimum–transition optimum equivocation**：future homogeneous HJB 不自动优化边界前的 gain trajectory。优先采用 preview-phase-augmented problem；无法闭合时降级为 future target + feasibility renderer，并收窄 claim。
3. **ANEOC speed circularity**：若 sensitivity 在事件后在线求出，其核心线性求解与 warm-start Newton 第一迭代同阶、局部等价。必须预维护 sensitivity 并以真实 runtime 证明差异。

最强反对意见是：**一个带 preview state 和 gain dynamics 的有限时域 NMPC，可能比 infinite-horizon HJB + critic + ANEOC 更直接地处理同一个问题。** 后续只有在最小动力学建立后，通过 state-feedback capability、求解时延、约束处理和 approximation error 的同条件比较，才能解除这一反对意见。

---

## 14. Stage 2 结论与下一 Gate

### 已通过

- Stage 1 spatial preview 与 cylinder-contact backbone 具有明确 RMC；
- future information 必须改变 applied feedback mechanics，而不是只改变预测；
- 原 NNBO 的 force-as-action 与 HC-specific impedance mapping 不可迁移；
- parameter-dependent optimal law + local continuation 在概念上可形成连续链。

### 尚未通过

- HJB 相对 preview-conditioned QP/NMPC 的 necessity；
- single critic 的可计算维度、root uniqueness 和 error bound；
- ANEOC 相对 warm-start Newton / parametric NMPC sensitivity 的 deadline advantage；
- Stage 1 observer 输出到最小 policy-sufficient \(\eta\) 的映射。

### 下一阶段

从 Stage 1 圆柱–黏弹接触模型出发，完成：

1. contact-frame reduced robot dynamics；
2. applied-gain dynamics；
3. preview phase \(\sigma_b\) 与 transition mixing；
4. equilibrium-centered minimal Markov state；
5. candidate descriptor \(\eta\) 及其 observable/unobservable decomposition；
6. HJB 与 finite-horizon NMPC 两种 formulation 的同对象对照。

只有完成这些对象后，才能决定 HJB/critic/ANEOC 是最小充分解，还是应改成更简单的 parametric NMPC sensitivity route。
