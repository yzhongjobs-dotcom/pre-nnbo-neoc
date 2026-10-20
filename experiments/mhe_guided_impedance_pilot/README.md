# 约化MHE—阻抗执行工程预检（已暂停）
<!-- project-update:start -->

**历史状态说明（2026-10-20审查）：** “已暂停”、未训练完成和下文“恢复前先补基线”只指本目录当时原型。后续两候选已有独立训练、仿真与比较，见[最新任务仿真](../dual_candidate_task_simulation/README.md)。当前不恢复强制NNBO或先完成全部对照的关卡；执行顺序见[长期要求](../../project_records/long_term_requirements.md)，当前任务见[HANDOFF](../../project_records/HANDOFF.md)。

<!-- project-update:end -->

2026-10-19。用户指出本轮将固定阻抗调优置于主要比较位置，未落实NNBO及领域强方法；已立即停止NN训练。保留本目录以复用估计/执行工具和复核实际正负结果，**不作为候选有效性、领域优势或顶刊贡献验证**。恢复工作优先落实[领域主对照](../../research/candidate_improvements_targeted_literature_and_design.md#4-用户纠正主对照必须落实到领域方法)，不先续训。

[协议及边界](design_and_protocol.md)：局部七参数约束滚动窗联合估计，非完整B0 MHE；普通两段打靶，非HJB/NEOC；拟任务梯度引导NN，非PPO。真值接触和增益/能量执行只读复用此前修正模型，旧结果未覆盖。

已完成两例训练分布pilot、固定辅助在12例调优、24例固定/解析轨迹各一条及3600个相关训练上下文。训练进程中断，无actor权重、验证集或独立测试结果。28个相关状态的模型梯度诊断不等于28例独立闭环试验。零发布等待成绩不构成部署实时性。

核心文件为`core.py`、`frontend.py`、`run.py`；`tune.py`仅复现固定辅助调优；`diagnose.py`复现两例局部梯度诊断。`first_pilot*`、`tuned_fixed.json`、`train_contexts*`和`contexts.npz`保留原始数据；不要按运行新默认覆盖它们。首pilot轨迹29列，后采集轨迹34列（另加anchor/sensed）；索引在`run.py`末尾、诊断使用首pilot格式。依赖NumPy/SciPy/Numba/Torch；运行环境见`environment.json`。

数值检查（不等于完整控制验证）：

```bash
python experiments/mhe_guided_impedance_pilot/checks.py
```

检查均匀场接触模型一致、任务梯度的独立方向差分、延迟未到达样本不能影响估计。结果见`checks.json`。NN脚本的部署路径不读取材料真值；仿真训练代价允许使用特权信息。未来继续前仍需补全领域主基线、真实计时发布、积分/足迹收敛和强任务测试。
