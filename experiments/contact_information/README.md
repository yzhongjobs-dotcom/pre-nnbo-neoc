# 第一个模块：带宽、接触信息与局部非平稳估计

科学设计、文献来源、数学假设、正负结果及贡献边界见 [研究正文](../../research/first_module_information.md)。用户已授权联合研究主动激励/接触表示/估计；未改原 V1＋Proposed、对象、IMM/导纳参数或冻结快照。当前候选没有正式采纳。

核心文件：local_contact.py 为纯测量接口的滤波对齐回归和数据无关探测优化；joint_policy.py 将合格的当前增益接入原控制器，并维护独立原增益流；run_joint.py 负责真值屏蔽、固定场景与评分。只有评分读取参考刚度和真实力，且在回合结束后添加。

当前较强实例：原20.8 Hz、小幅探测；32拍一阶时变增益、四次漂移、无速率项；原反馈/IMM/恢复不变。原探测纳入可行集后，其名义当前增益信息优于所搜周期 Fourier 候选。保留原探测是比较结果，不是预设必须替换。

## 有限且有目的的证据

- joint_pilot_v1：指令预算/增大幅值的首轮正负结果；联合方案反而变差。
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
