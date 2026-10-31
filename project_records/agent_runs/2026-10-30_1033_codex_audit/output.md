**总体判断：当前仿真不能原样支撑第五节的全部结论。** 完整方法的已保存性能结果具有参考价值，但发现了一个足以改变核心消融结论的确定代码 bug，以及场景混杂、对照适配和统计报告问题。不能据此否定完整方法，也不能继续用现有解释宣布机制已被充分验证。

全程未修改文件；诊断仅在内存中执行。以下按严重程度排序。

1. **【阻断】“去掉保守刚度”消融实际同时破坏了前馈和激励。**

   **位置：** [proposed.py:78](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:78)、[proposed.py:165](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:165)。

   `conservative=False` 提前返回，跳过 `k_good` 更新。于是 `k_good` 始终为 NaN，`k_seek` 最终降到 **200 N/m**，前馈与激励都被错误地按这个尺度换算。该消融并非只把反馈中的保守刚度换成估计值。

   我在内存中仅恢复“精确窗口更新 `k_good`”这条逻辑，结果如下，单位均为 N：

   | 场景 | 原消融 RMSE | 修正后消融 RMSE | 完整方法 RMSE |
   |---|---:|---:|---:|
   | T2，9601 | 391.357 | 2.089 | 1.784 |
   | T2，9602 | 65.043 | 0.793 | 0.892 |
   | T3，9700 | 86.426 | 0.265 | 0.265 |

   这不是完整重跑，但已足以否定“原消融的巨大恶化证明保守刚度必要”的解释。

   **方向：偏向所提方法。必须修正：是。最小修正：** 将 `k_good` 的精度判断与更新移到提前返回之前，保持前馈和激励逻辑一致，重跑这一项消融。第五节关于“去掉保守刚度后发散”“各部件必要”的结论应暂撤。

2. **【阻断】T1 实际有呼吸，不能隔离刚度归一化贡献。**

   **位置：** [scenarios.py:39](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/scenarios.py:39)、[run_experiments.py:202](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:202)、[论文:235](E:/github_project/scan-impedance-update/research/paper_draft_zh.md:235)。

   T1 默认呼吸幅值为 2–5 mm，实验入口没有清零。主测试12个 T1 场景实际幅值为 **2.1–4.9 mm**。因此，所提方法拥有运动前馈，而 Rossi、Wij 型没有；目前的 T1 优势同时包含呼吸补偿收益。

   **方向：偏向所提方法的机制归因。必须修正：是。最小修正：** 如果保留“静止组织上的刚度跨度检验”，设置 `A=0` 重跑；现有结果可以保留，但应称为“呼吸下的非均匀组织扫描”，不能据此单独证明归一化使响应一致。

3. **【重要】Rossi 型存在确定的信号处理失配，其低估不能全部归因于组织运动。**

   **位置：** [comparators.py:48](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/comparators.py:48)、[plant.py:145](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:145)、[tune_equal.py:28](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/tune_equal.py:28)。

   Rossi 回归使用“已经过30 Hz低通的力”和“未经同样处理的编码器位移”。之后两个通道再经过相同带通，**不能消除之前的幅相失配**。在固定20.8 Hz激励下，这一误差不可忽略。

   静止、均匀、无呼吸的内存抽查中，Rossi 型估计/真值中位数仍只有 **0.809、0.616**；这已经排除了呼吸作为唯一解释。

   此外，原文实验采用0.5 Hz力参考激励及不同带通实现，本代码的20.8 Hz位移激励属于任务适配，不能笼统称为“按原文方程忠实复现”。原文见 p.1144。:codex-file-citation{path="E:/github_project/scan-impedance-update/literature/manual_inbox/Implicit_force_control_for_an_industrial_robot_based_on_stiffness_estimation_and_compensation_during_motion.pdf" purpose="source"}

   **方向：偏向所提方法。必须修正：是，至少修正比较与归因。最小修正：** 为对照提供测力链路一致的位移预处理，并允许选择合理激励频率；分别报告静止时的处理偏差和运动带来的新增偏差。

4. **【重要】DO-MPAC 的调参自由度没有覆盖其主要性能旋钮，现有结果不足以代表“扰动补偿路线”。**

   **位置：** [tune_equal.py:33](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/tune_equal.py:33)、[comparators.py:215](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/comparators.py:215)。

   `Np` 固定15，`Q、R` 独立搜索，但无约束标量预测控制主要取决于 **Q/R 和预测时域**。独立搜索两个尺度浪费预算，也没有检验更长时域。

   已选 DO-MPAC 的 `Q/R≈70.55`，对应名义双积分器反馈增益约为 `[0.5073, 0.04257]`，连续等效极点约为：

   \[
   -0.0203\pm0.7120j .
   \]

   这是衰减很慢的名义闭环，与所提方法经过稳定性预筛的控制器并非相近的动态能力。相同61次评估不能自动排除搜索空间削弱对照。

   **方向：偏向所提方法。必须修正：若保留路线级优越性结论，是。最小修正：** 用 `Q/R` 参数化，开放 `Np`，先验证静止接触下的合理响应范围；再做同预算选择。结论应暂限定为“优于当前 DO-MPAC 适配实例”。

5. **【重要】A₁ 设计脚本与保存的设计、候选集不一致，无法按说明复现。**

   **位置：** [design_A1.py:19](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/design_A1.py:19)、[design_A1.py:120](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/design_A1.py:120)、[A1_design.json:1](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/A1_design.json:1)。

   脚本使用 `KI=30`，质量搜索范围为 \(10^{-6}\)–\(10^{-4}\)；保存设计却为 `ki=0、M1=0.01`，最终调优实例为 `ki=0、M1=0.001`。当前脚本也没有生成所用535个可行设计的完整流程。

   对最终参数直接检查：`KI=30` 时冻结稳定性检查失败、圆检验失败；设为实际使用的 `KI=0` 后，两项才通过。

   **方向：中性，影响可复现性与理论对应。必须修正：是。最小修正：** 统一设计脚本、参数范围和候选生成流程，并明确最终实际使用的设计。

6. **【重要】圆检验是简化模型上的数值筛选，不能称为实际 plant 的稳定性认证。**

   **位置：** [design_A1.py:31](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/design_A1.py:31)、[design_A1.py:90](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/design_A1.py:90)。

   当前模型保留了控制周期、延迟和滤波的大体结构，但没有 SLS、速率项、饱和和接触恢复；伺服离散化也不同于 plant。`circle_ok` 检查1500个频点的圆外距离，配合有限个冻结增益检查，是合理的工程筛选，**不是连续频域和变参数执行系统的严格证书**。

   **方向：偏向所提方法的稳定性表述。必须修正：表述必须；控制器未必需要改。最小修正：** 改称“简化模型数值稳定性筛选”，补充实际模型的扇区越界与饱和统计。

7. **【重要】列8的 \(k_{\mathrm{eff}}\) 与窗口回归刚度只有近似对应。**

   **位置：** [plant.py:108](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:108)、[plant.py:133](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:133)、[common.py:77](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/common.py:77)。

   列8定义为弹性足迹切线刚度的1.4倍，即 SLS 的高频极限刚度。局部线性化、静止工作点下，实际接触更接近：

   \[
   G(s)=k_0\!\left(1+0.4\frac{s}{s+1/0.3}\right)+0.3F_0s .
   \]

   因此它还包含频率相关的损耗项和速率阻尼。20.8 Hz下，SLS储能刚度确实接近 \(1.4k_0\)，所以列8不是随意定义；但它不是包含 `lam` 的完整动态刚度，窗口标量回归也不保证恰好估计它。

   **方向：主要偏向所提方法的“无偏”和覆盖率解释。必须修正：是，定义与结论必须限定。最小修正：** 明确列8是“瞬时弹性切线刚度参考”；把 \(\hat k/k\)、κ和覆盖率称为相对此参考的经验统计，不据此直接证明实际动态回路满足稳定扇区。

8. **【重要】direct 模式没有发现未来泄露，但在线观测与离线 S1 数据并非完全一致。**

   **位置：** [plant.py:229](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:229)、[plant.py:142](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:142)、[common.py:65](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/common.py:65)。

   策略读取临时行时，力值还是上一次1 kHz采样后的滤波值；随后 `advance()` 在当前边界采样并覆盖该行。在线窗口末端读取的值，与实验结束后离线回归看到的同一行不同。

   因而末端力最多落后编码器约1 ms，且 `_build` 对末端区间假定的8次滤波更新与在线实际可用更新不完全一致。抽查场景中，在线/离线估计比值中位数约 **0.9968**，说明影响通常较小，但“S1完全检验在线估计器”不准确。

   **方向：对实际控制共同影响；可能美化所提方法的离线估计证据。必须修正：是。最小修正：** 固定采样—控制—记录顺序，保存策略实际消费的观测，S1使用该数据。

9. **【重要】S1 的推拉检验包含截断事件，并被大量非运动时段稀释。**

   **位置：** [run_experiments.py:98](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:98)、[run_experiments.py:104](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:104)。

   先按8 s的 T2生成事件，再强制只运行6 s。有些事件在6 s后才开始，有些回撤被截断。一个中等组织“推拉”场景中，真正推拉斜坡只占统计区间的 **3.2%**。

   全时段中位数和10%–90%分位数容易遗漏最困难的运动阶段，不能支持“推拉时估计仍准确”的充分结论。

   **方向：偏向所提方法。必须修正：是。最小修正：** 按实际时长生成完整事件，并单独报告推入、保持、回撤及恢复阶段的估计误差和覆盖率。

10. **【重要】硬组织快速推拉属于强压力测试，不能直接当作代表性超声扫描工况。**

    **位置：** [scenarios.py:45](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/scenarios.py:45)、[scenarios.py:47](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/scenarios.py:47)、[plant.py:166](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/plant.py:166)。

    平滑阶跃的峰值速度为 \(1.5A/t_r\)，给定范围对应 **9–80 mm/s**；主测试实际达到约54–60 mm/s。4 N平衡点下，0.5–1.5 MPa均匀介质的压入量仅约 **0.455–0.175 mm**，瞬时刚度约14.3–36.4 kN/m。几十毫米每秒的表面运动叠加伺服滞后，容易产生远大于4 N目标力的瞬态。

    此外，plant假定机器人在大接触力下仍按一阶位置伺服运动，没有负载相关柔顺或加速度限制。它允许评价算法压力响应，但不足以支撑临床代表性。

    主测试中，所提方法硬组织误差平方总量约 **80%来自两个场景**。

    呼吸幅值2–5 mm本身不能判错，但这里是**单边幅值，峰峰值4–10 mm**；需要区分扫描部位和体位。人体胸壁研究确实观察到毫米级呼吸位移，但不能直接标定本项目的脊柱工况。[原始研究](https://journals.physiology.org/doi/full/10.1152/jappl.1997.83.5.1531)

    **方向：混合；有利于预测前馈，也会严重惩罚其预测误差。必须修正：代表性主张必须。最小修正：** 将该组合明确列为压力测试，分阶段报告误差贡献，避免少数极端事件主导整体解释。

11. **【重要】“缩放即最优更新”超过了 scaling_check 的检验能力。**

    **位置：** [scaling_check.py:18](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/scaling_check.py:18)、[scaling_check.py:34](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/scaling_check.py:34)、[论文:252](E:/github_project/scan-impedance-update/research/paper_draft_zh.md:252)。

    四档模量都选中 `D1×0.5`，恰好是搜索下边界；`M1`固定，只搜索两个参数的粗网格，每档只有4个场景。因此结果只能说明“四档条件的最佳已测网格点相同”，不能说明连续最优参数相同，更不能证明最优阻抗只差比例因子。

    **方向：偏向所提方法的理论验证表述。必须修正：是。最小修正：** 当前结论降为“与缩放关系相容”；若继续声称最优，至少扩展边界并检验网格分辨率。

12. **【重要】Wijayarathne 型是特权刚度 PI 代理，不是原方法复现。**

    **位置：** [comparators.py:114](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/comparators.py:114)、[comparators.py:157](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/comparators.py:157)。

    原文辨识弹簧、阻尼两个参数，并采用协方差映射；这里替换成带噪真值的标量RLS，还把 `tanh(innovation)` 用作常规分母遗忘因子。其适应速度含义与原映射不等价。原文的20 Hz是**探头激励频率**，不能直接推出刚度标签应每50 ms才更新一次。[原文§III-B、IV-B–C](https://arxiv.org/html/2004.09729v1)

    真值信息有利于对照，但人为采样与遗忘机制也可能削弱它，不能简单认为这些改动全部“更有利”。

    **方向：混合。必须修正：比较身份必须。最小修正：** 明确称为“带噪特权刚度PI代理”；检验更新率与遗忘机制敏感性，避免把结果推广到原论文整体方法。

13. **【重要】窗口有效性保护没有覆盖前馈与激励通道。**

    **位置：** [proposed.py:83](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:83)、[proposed.py:96](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:96)、[proposed.py:130](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:130)。

    `precise` 成立但 `chi2` 不通过时，仍更新 `k_good`；`drift_pred` 只要求估计有限，就允许使用漂移外推。也就是说，反馈刚度虽然保持保守，失配窗口仍可能通过前馈产生大动作。

    恢复接触时又只重置部分状态，激励相位差状态没有同步重置。这些并非已证明导致全部大峰值，但说明“无效窗口受到保护”只对部分通道成立。

    **方向：可能低估所提方法性能，也削弱安全性解释。必须修正：逻辑说明必须；控制修正需验证。最小修正：** 分别定义反馈、前馈和激励的有效条件，并明确恢复时各状态如何衔接。

14. **【重要】失败场景被 NaN 聚合删除，配对统计变成了成功者比较。**

    **位置：** [run_experiments.py:161](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:161)、[run_experiments.py:175](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:175)、[S3_main.json:2440](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/results/S3_main.json:2440)。

    主测试 Rossi＋DO-MPAC 的9710场景失败，重复测试也有一个失败。RMSE使用 `nanmean`，bootstrap删除非有限差值，丢失次数使用 `nansum`。因此表中的均值、次数和区间均没有完整包含这些失败。

    **方向：偏向发生失败的对照。必须修正：是。最小修正：** 单列失败率与有效配对数，明确RMSE是条件于成功的统计；失败不能在安全指标中贡献零次事件。

15. **【重要】接触丢失次数不能代表接触保持能力，峰值也只是125 Hz采样峰值。**

    **位置：** [scenarios.py:75](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/scenarios.py:75)、[common.py:149](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/common.py:149)、[run_experiments.py:74](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:74)。

    短暂反复过阈值会产生很多次丢失，长时间失联可能只有一次。主测试硬组织中，所提方法虽然只有28次、Wij型162次，但累计失联时间约为 **8.59 s与7.70 s**，次数并不能支持前者接触保持更好。

    plant内部步长为0.5 ms，但峰值和接触事件只在125 Hz日志上统计，可能漏掉周期内尖峰与短暂失联。`warm=1 s`则排除了启动辨识代价，不能用于证明启动更新快。

    **方向：当前接触保持叙述偏向所提方法；峰值漏采共同影响所有方法。必须修正：是。最小修正：** 并列报告丢失比例、最长失联、次数和失败率；明确采样峰值及预热范围。

16. **【中等】激励必须计入主 RMSE；现有“>10 Hz振荡”和纹波分解不够严格。**

    **位置：** [common.py:130](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/common.py:130)、[scenarios.py:81](E:/github_project/scan-impedance-update/experiments/cobot_scan_benchmark/scenarios.py:81)。

    激励产生的真实力变化会影响接触和成像，**计入主 RMSE是合理且必要的**。可以并列报告去激励谱线后的误差，用于区分调节误差与辨识代价，但不能替代原始RMSE。

    当前用“减去约0.1 s移动平均”近似高通，并不等价于严格的>10 Hz频带；随后只减去全程固定幅相的20.8 Hz正弦，不能完整捕捉幅度变化、停启和边带。

    **方向：去线后的指标可能美化有激励的方法。必须修正：指标命名必须。最小修正：** 保留原始RMSE，增加清楚定义的频带或同步谱线分解，并同时报告原始振荡量。

17. **【中等】σ目前是设计噪声尺度下的精度指标，不是已校准的统计标准差。**

    **位置：** [common.py:111](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/common.py:111)、[proposed.py:76](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/proposed.py:76)。

    理想白噪声条件下，当前滤波创新的标准差约 **0.0239 N**，代码使用0.04 N，偏保守；但模型失配、闭环相关性和第8项的末端时序差异没有进入σ。

    因此经验覆盖率可以报告，却不能把 `k_hat+cσ`直接解释为具有标称概率的统计置信上界。

    **方向：偏向所提方法的不确定性保证表述。必须修正：是。最小修正：** 将σ称为“模型条件下的精度尺度”，覆盖率称为经验覆盖率，区分噪声与模型失配。

18. **【中等】配对 bootstrap 的单位选对了，但样本不足以支持“相当”，重复验证也有范围限制。**

    **位置：** [run_experiments.py:171](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:171)、[run_experiments.py:94](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:94)、[论文:273](E:/github_project/scan-impedance-update/research/paper_draft_zh.md:273)。

    以场景为配对重采样单位是正确选择。问题是T2只有5个软、7个硬场景，又存在明显重尾，百分位bootstrap区间不宜承担强等效性判断。

    **区间包含0只表示没有检出差异，不证明性能相当。** 重复测试中，软组织相对Rossi的差异也不再显著。

    另外，S1的随机种子固定9800，不受 `REPLICATE`影响，两份S1结果逐字节相同；S1没有得到独立重复验证。代码中未发现直接把测试种子用于调参的证据。

    **方向：偏向所提方法的结论强度。必须修正：是。最小修正：** 改为“未检出差异”，报告实际样本数和逐场景差值；明确哪些实验真正重复。

19. **【中等】计算时间的统计口径与论文描述不符。**

    **位置：** [run_experiments.py:63](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:63)、[run_experiments.py:167](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/run_experiments.py:167)、[make_tables.py:85](E:/github_project/scan-impedance-update/experiments/motion_aware_scan/make_tables.py:85)。

    实验在 `Pool(16)` 中运行；报告的是场景P99的中位数，再取任务最大值，不是单进程总体P99或最坏场景P99。

    原始记录中，所提方法某场景P99为 **16.63 ms**，重复测试另有 **9.36 ms**，均超过8 ms。因此“均远低于8 ms”不受当前数据支持。它也可能受并行调度影响，不能反过来据此判定控制器必然超时。

    **方向：偏向所提方法的计算时限表述。必须修正：是。最小修正：** 在独立串行环境重新计时，报告总体P99、最坏场景P99和超周期比例。

另外，以下几项**没有查到对应的确定错误**：

- direct模式的指令队列确实实现了一拍延迟，没有额外叠加一拍；第8项是采样边界不一致，不是未来信息。
- 非特权方法代码未读取真值列。将所提方法可见的真值列全部替换成NaN后，抽查轨迹完全一致。Wij与Oracle的真值读取是显式特权。
- 当前Rossi采用固定位置幅值激励，没有发现持续锁在 `K_MAX` 的现象。
- DO-MPAC的坐标符号转换、凝聚预测增益和ESO矩阵指数离散化没有发现明显代数错误，不能指责为显式欧拉失稳。
- 对6个抽查场景，在4 N平衡点附近将足迹积分从12点提高到192点，力与刚度误差很小，未发现足迹求积导致主要结论失真的证据。

**修正后预计会改变的结论是：**

- “保守刚度是避免发散的必要部件”很可能明显削弱，现有巨大消融差距主要受bug污染。
- T1的现有排名可以描述为该混合工况下的结果，但“优势来自刚度归一化”的归因必须重新验证。
- Rossi估计失效的程度可能减轻；相对Rossi和DO-MPAC的优势幅度、显著性需要重新判断，不能预先保证保持。
- “跨模量最优参数只需缩放”“与对照相当”“接触保持更好”“均满足8 ms”都应降级或改写。
- 漂移投影和前馈具有任务收益的证据仍值得保留，但应限定为当前合成模型与被测实例，不能用现有第五节宣布全部机制与公平优势已经成立。