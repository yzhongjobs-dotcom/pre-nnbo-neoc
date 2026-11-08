1. **第2节公式：仅在条件下正确；数值正确。**  
   \(\tau=\sqrt{N/q}\)、噪声方差 \(N/(2\tau)\)、斜坡滞后 \(\rho\tau\)，适用于连续时间或小步长近似，且 \(q\) 是单位时间过程噪声强度。这里的方差仅计测量噪声贡献；若真实对象也按随机游走变化，总跟踪误差还包含过程噪声，不能仍写成 \(N/(2\tau)\)。

   对确定性斜坡，
   \[
   \mathrm{MSE}=\frac{N}{2\tau}+\rho^2\tau^2
   \]
   确实给出所写的 \(\tau^*\)。取精确值 \(N=0.00768\)，三组最优 RMSE 为 **0.07931、0.14249、0.19384**，与笔记一致。\(2\mathrm{sd}\le0.09\) 给出 **\(\tau\ge1.8963\) s、\(q\le0.002136\)**，约写为 1.9 s、0.002 正确。这是该滤波器族的近似结果，不是实际闭环或所有估计器的可达下界。

2. **Huang 闭式解：仅在条件下正确；“32%／0.33r”必须区分含义。**  
   按要求用 pdftotext 核对了[原论文](E:/github_project/scan-impedance-update/literature/manual_inbox/Novel%20sparseness-inducing%20dual%20Kalman%20filter%20and%20its%20application%20to%20tracking%20time-varying%20spatially-sparse%20structural%20stiffness%20changes%20and%20inputs.pdf)第6页。固定精度、标量观测系数为1、关闭参考层时，若 \(\alpha\) 表示**实际变化方差**、\(S=P_0+r\)，则
   \[
   \alpha^*=(\nu^2-S)_+
   \]
   正确。原文式14使用的是缩放方差 \(\lambda\)，对应 \(\lambda^*=(\beta\nu^2-\Gamma)_+\)。增益 \(1-r/\nu^2\) **只适用于触发分支**；未触发时为 \(P_0/(P_0+r)\)。

   校准的高斯新息满足 \(P(|\nu|>\sqrt S)=\mathbf{31.731\%}\)。进一步取 \(P_0\approx0、S\approx r\)，更新量的方差为：

   - **每一步无条件：\(0.333262r\)**；
   - **每次触发条件下：\(1.050271r\)**；
   - 若指注入的过程方差 \(\alpha\)，其每步期望则是 **\(0.483941r\)**。

   因此笔记的 \(0.33r\) 仅能解释为**每步更新量方差**。32%也不能直接当作当前重叠膨胀代码的实际触发率。另外，“原文§4.2关闭参考层”的描述错误：第15页关闭该层的是对照试验，论文报告其出现漂移。

3. **代码实现：固定精度 HuangARD 基本正确；精度学习错误；Forgetting 描述需限定时序。**  
   [HuangARD](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/memory_filters.py:135)固定精度分支与上述闭式解一致，零方差以 \(10^{-12}\) 近似。学习精度分支并非原文式15／17的正确标量实现：

   - 原文所有相关协方差统一按 \(1/\beta\) 缩放；代码只把观测方差改为 \(R/\rho\)。
   - 以 \(E\) 表示未除精度的协方差，正确标量关系是  
     \[
     b_i=b_{i-1}+\frac{\nu^2}{2(R+E_i^-)},\qquad
     P_i=\frac{E_i^-R}{(E_i^-+R)\rho_i}.
     \]
     代码把实际方差 \(P_0+\alpha\) 与基准 \(R\) 混加，尺度不一致。
   - 新精度未用于当拍重新更新，也未与变化方差联合迭代；初始化实际为 \(b_0=1\)，不是原文的 \(10^{-10}\)。

   因而“学习精度几乎无记忆”的结果不能归因为原文算法。[Forgetting](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/memory_filters.py:83)确实依据 \(\nu^2/S\) 调节遗忘，但顺序是**先更新均值，再膨胀方差**；当前异常使下一次观测增益增大。返回的 \(P\) 已膨胀，并非当前测量更新后的纯后验方差。

4. **attach：重叠系数有条件正确；“门控不变”错误。**  
   \(nw/hop=24/6=\mathbf4\) 算术正确。对等权移动平均白噪声，相关系数为 \(0.75、0.5、0.25\)，长期方差膨胀恰为4；对实际 Huber-D⁴ 窗口，这只是信息率近似，不能保证协方差精确校准。

   [实际路径](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/memory_filters.py:209)返回原始 χ，因此 χ 判据保留；但外层 `precise`、`meas_bound` 等判据使用的是**滤后 k 和 σ**。跳变检查通常恢复原始窗口分子，分母 `k_good` 仍由滤后结果产生；若始终没有符合原始精度规则的窗口，`state['raw']` 为空，滤后 `k_check` 可以直接进入跳变检查。

   `usable` 的 χ 筛选合理，但精度条件用本窗口 k 自归一化，会偏向接纳高估值，与宿主已有的先验参考门控不同；它还没有排除 \(k=200\) 的截断值。

   **额外实质错误：实验中的“现有窗口”也经过 attach。** `Window` 输出被每6拍保持一次，且 σ 约变为原始值的 **2倍**。因此结果表中的 Huber-D⁴ 并非未经修改的现有逐拍窗口基线。

5. **聚类与指定两个区间：正确；默认 report 不完全复现笔记。**  
   [boot](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/memory_closed_loop.py:117)确实整簇抽样，重复抽到的簇也正确重复计入。筛选有探头质量、按物理场景 `(scene,i)` 聚类，对[原始结果](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/memory_closed_loop_test.json)以5000次、seed=0复算：

   | 配对差 | 回合／簇数 | 均值差及95%区间 |
   |---|---:|---|
   | Huber-D⁴ − Proposed，nominal全部场景 | 84／21 | **+0.086153 [0.015051, 0.169994]** |
   | Rossi-ID − Proposed，slow_thin | 36／9 | **−0.127515 [−0.175959, −0.080907]** |

   均与笔记舍入值一致。不过 `report` 当前调用2000次，默认参照也是窗口，而非 Proposed；2000次所得区间分别为 **[0.015229, 0.170626]、[−0.175319, −0.080658]**。

6. **第1节表格正确；“f=1.0正好在悬崖边”仅获部分支持。**  
   [gain_need_sweep.json](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/gain_need_sweep.json)中，有探头质量条件下，表内全部 L 均值均与所写精度一致。低于1时明显恶化、A的28°条件在1.2时明显改善，支持低侧陡变；但离散扫描不足以定位“正好”的边界，B中1.0反而优于1.2。

   “A最大力35→11 N、均值约减半”的参照实际是**现有估计器→f=1.2**：35.47→10.68 N，L为0.15245→0.07655。若指 **f=1.0→1.2**，正确值应为 **79.06→10.68 N，L为0.43694→0.07655**。

全程只读，未读取 HANDOFF、未修改仓库文件、未运行闭环批次，未提交或推送。