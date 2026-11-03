# 第一个模块：带宽、接触信息与局部非平稳估计

科学设计、文献来源、数学假设、正负结果及贡献边界见 [研究正文](../../research/first_module_information.md)。用户已授权联合研究主动激励/接触表示/估计；未改原 V1＋Proposed、对象、IMM/导纳参数或冻结快照。当前候选没有正式采纳。

核心文件：local_contact.py 为纯测量接口的滤波对齐回归和数据无关探测优化；joint_policy.py 将合格的当前增益接入原控制器，并维护独立原增益流；run_joint.py 负责真值屏蔽、固定场景与评分。只有评分读取参考刚度和真实力，且在回合结束后添加。

当前较强实例：原20.8 Hz、小幅探测；32拍一阶时变增益、四次漂移、无速率项；原反馈/IMM/恢复不变。保留原探测只因尚未得到经公平检验的替代，不代表原探测最优。

继续接触模型的最新证据见[正文 §8–10](../../research/first_module_information.md)：已实现端点/移动参考曲率、正参数局部 SLS、相同指数漂移消融与12个定向闭环实例。曲率筛选只在硬组织带来很小差值；SLS改善该实例的条件预测，标量增益接入却使闭环变差。尚未采纳，下一项仍是动态接触观测的物理/信息闭合，不能用这些结果裁决整条接触路线。

## 有限且有目的的证据

**共享运动与动态响应候选（最新）：** 具体设计/推导/正负结果见[研究正文](../../research/shared_contact_response.md)。`structured_memory.py`用同一局部运动作用弹性和记忆，`dynamic_compliance.py`实现时变逆、作用时域信息及历史/噪声传播，`response_policy.py`接入原控制器并模拟32 ms识别发布延迟。首版、历史重构、观测方差三轮共10次闭环，均未形成可采纳的通用改进；严重历史版负结果保留。没有改变对象、原20.8 Hz或反馈/IMM参数。

```powershell
python -X utf8 experiments/contact_information/verify_response.py
python -X utf8 experiments/contact_information/check_response_results.py
# 重做闭环时用新标签；保存七通道输入以供无真值前缀重放。
python -X utf8 experiments/contact_information/run_response.py --label response_reproduced
python -X utf8 experiments/contact_information/run_response.py --models shared_history --label history_reproduced
python -X utf8 experiments/contact_information/run_response.py --models shared_history_uncertain --label uncertainty_reproduced
# audit_response.py当前默认主标签为含同自由度对照的v2，拒绝覆盖。
```

主数据：`shared_response_v1.json`、`shared_response_history_v2.json`、`shared_response_uncertainty_v3.json`及对应测量NPZ；`shared_response_same_input_matched_v2.json`为同输入主对照（包含同四个干扰参数的自由表示），未加该对照的旧回放也保留。每个模型信息、目标状态、因果发布和计时边界见正文。源码新增兼容分支之后，前期源码哈希不等同当前文件；`check_response_results.py`对v1与v3代表回合仅给测量前缀，指令逐点重现，未以当前哈希冒称原文件逐字节未变。协方差是局部条件参考，不是初态/模型误差或真实闭环的认证。

**已移除的激励实验：** 降频、等幅对照、有界偏移、过零更新、自适应调幅、指令预算波形和曲率探测设计已移出工作区（Git `dd2884a` 可恢复）。2026-11-03复核：波形比较只抽192个候选、没有真正优化，主实验无相位自由度、最高31.25 Hz，且峰值约束不利于多频；换激励时反馈刚度仍用为20.8 Hz配的16拍窗口；候选一律沿用基线参数、未重新调参；目标只算名义方差；每族仅1个开发实例，没有快速推拉。这些结果不能说明正弦更好，也不能否定激励优化路线。保留 `response_band_information.py` 这一固定设计计算，尚未复核。

**原激励的任务价值：** `audit_probe_information.py`在冻结测量中核对探测支持的信息，恢复实际探测指令并逐点验证原控制指令；`probe_value_intervention.py`保持精简估计器，只有启动后半幅/关闭两个干预。科学推理见[主动激励§7](../../research/active_excitation.md)。这些只是单实例诊断。

```powershell
python -X utf8 experiments/contact_information/verify_probe_value.py
# 重做时使用新标签，保留原始输出；重叠窗口不是独立重复。
python -X utf8 experiments/contact_information/audit_probe_information.py --recorded-probe --label probe_information_reproduced
python -X utf8 experiments/contact_information/probe_value_intervention.py --label probe_value_reproduced
```

主结果是 `probe_information_actual_command_v2`、`probe_value_intervention_v2`，汇总/核验为`probe_value_checkpoint`。关激励只在共同启动后，不代表全程无辨识。

本次接触细化的代码：`contact_jet.py`（局部曲率及端点功能量），`curved_policy.py`（原精简候选回退），`memory_contact.py`（解析SLS卷积、未知表面/初值剖面及正参数约束），`memory_policy.py`（只验证探测实部增益的有限接口，尚无动态逆）。这些估计源均无对象导入；运行器使用原有测量屏蔽器。

```powershell
python -X utf8 experiments/contact_information/verify_contact_jet.py
python -X utf8 experiments/contact_information/verify_memory_contact.py
python -X utf8 experiments/contact_information/audit_refinement.py
# 已有 lean_input 文件会保留；freeze 只生成缺少的精简候选输入。
python -X utf8 experiments/contact_information/run_refinement.py freeze
python -X utf8 experiments/contact_information/run_refinement.py replay --label contact_jet_reproduced
python -X utf8 experiments/contact_information/run_refinement.py replay --label contact_memory_reproduced --models lean drift_040 memoryc_040
python -X utf8 experiments/contact_information/run_refinement.py closed --label contact_curvature_closed_reproduced
python -X utf8 experiments/contact_information/run_refinement.py memory_closed --label contact_memory_closed_reproduced
```

本次原始结果按迭代分开保留：`contact_jet_replay_v1`/`contact_jet_replay_centered_v2`、`contact_curvature_closed_v1`；`contact_memory_replay_v1`/`contact_memory_constrained_v2`/`contact_memory_drift_ablation_v3`、`contact_memory_gain_closed_v1`。预测评分的未来编码器是给定条件，未来力隐藏；不是在线运动预测。重叠窗口不作为独立重复。`contact_refinement_audit`汇总完整性、数值核验和记忆闭环 v1 的继承描述字段纠正。

以下为上一轮联合优化的已完成证据，仍保留：

- joint_response_pilot_v2 的“仅新探测”“联合”两列不作为新探测无效的证据，原因见上。
- joint_response_pilot_v2：物理响应预算及长常数窗/去速率/探测单改/联合对照；两个原场景。
- joint_abrupt_motion_v3：只增加一个硬组织骤变运动实例，检查长窗边界。
- identical_replay_v1：同一冻结测量，分离估计表示与闭环数据变化；日志瞬时参考只评分，不称最优控制量。
- joint_value_attribution_v4：当前值与接受选择规则的诊断对照；实际闭环时刻不相同，不作严格可加归因。
- joint_noise_pair_v5：同一物理 S2，换一次噪声实现，只检验偶然噪声替代解释。
- joint_lag_pair_v6：同一物理 S2，已有滞后×1.5条件；双方从名义标定出发，V1 从测量辨识。

调用各脚本时使用新 label 保留主结果；下面以复现标签示例，不重跑全部矩阵：

```powershell
python -X utf8 experiments/contact_information/verify_local.py
python -X utf8 experiments/contact_information/replay_identical.py --label identical_replay_reproduced
python -X utf8 experiments/contact_information/run_joint.py --probe-budget response --include-original --models baseline joint_lean --label joint_lean_reproduced
python -X utf8 experiments/contact_information/run_joint.py --probe-budget response --models baseline long_constant varying_no_rate fit_only probe_only joint --label joint_response_reproduced
python -X utf8 experiments/contact_information/run_joint.py --probe-budget response --include-original --families T2m-stiff --models baseline joint_lean --label joint_abrupt_reproduced
python -X utf8 experiments/contact_information/run_joint.py --probe-budget response --include-original --noise-offset 100000 --families S2 --models baseline joint_lean --label joint_noise_reproduced
python -X utf8 experiments/contact_information/run_joint.py --probe-budget response --include-original --servo-lag-scale 1.5 --families S2 --models baseline joint_lean --label joint_lag_reproduced
python -X utf8 experiments/contact_information/summarize_joint.py
```

汇总只读取完整原批次，核验最终选中的原探测＋精简估计逐值复现前一轮相同实例，不挑各版本最好结果。数值核验包括独立影响行的精确协方差与高斯序列；不作为真实组织上的覆盖率保证。所有调用计时是 Python 原型，未模拟计算延迟，未上最终确认集；未来重复数量按科学主张决定。

v4–v6 原始元数据的 probe.selected_kind=original_sine、wave=null 指实际保留原探测；这些批次的 rms/peak/fourier 仍描述未采用的周期候选，不可据其宣称最终探测功率下降。后续代码将未选方案明确放入 periodic_candidate。实际力的探测相干功率始终由回合后的测量评分，主结果不受该显示字段影响。
