# 接触模型研究

历史方案和纠偏见 [专项研究](../../research/contact_model_study.md)、[方法审计](../../project_records/handoff/03_方法演进审计.md)。`p3_*`、`vp_map.py` 等是历史候选；`finite_patch_*` 与 `line_patch_model.py` 是 2026-11-02 的接触专项原型，未并入 Proposed Method。

新原型用同一接触模型的逆压入一致性与压力矩，比较局部线性接触、分布弹性层有限足迹及 E3 相关柱面关系。它是离线机制研究，不是可按 8 ms 部署的估计器。Ft 目前只作为力矩观测的力臂补偿输入，不是模型预测的输出；不要称为完整六维接触模型。

## 当前主要证据

- `results/finite_patch_input_*.npz`：五族各一个既定诊断场景的输入。`measured` 严格为 t/x/y/Fn/cmd/Ft/My；`scoring_truth_k` 仅评分。原始轨迹由 `ServoID(**CE.FINAL)` 产生。
- `finite_patch_ablation_v4.json`：S2/T3m 连续、不重叠 16 拍窗口；同记忆简单模型、去力矩和完整干扰模型消融。
- `finite_patch_audit.json`：物理推导、前缀因果性、偏置/双力臂、条件力验证、求解预算与局部几何/标定敏感性。
- `finite_patch_line_v5.json`：E3 柱面近似在同接口的稀疏初筛，与 v4 的采样协议不同，不直接排名。
- v0/v1/v2/v3 保留了被修正的实例，科学解释见专项研究 §4，不能挑其中最优结果作为最终方法。

所有计算使用原对象，不修改参数来适配候选；未使用 final 确认集。窗口不是独立场景。以上为初轮静态研究；后续动态候选闭环见下。

## 动态一致候选及闭环（未采纳）

`dynamic_contact.py` 分离弹性关系和局部松弛假设，直接拟合前向力；`dynamic_contact_policy.py` 比较量测增益和目标力压入修正接口，保留原方法独立回退。四个不重叠短窗共享接触参数，每窗独立表面运动；退出窗口只计入一次到达信息。`validate_next` 只在后续新窗口检验旧接触参数，残差小不代表弱激励下模型必然正确。原对象、反馈、IMM、激励和原调参不变。

```powershell
python -X utf8 experiments/contact_model/run_dynamic_contact.py --verify-only
python -X utf8 experiments/contact_model/audit_dynamic_contact.py
python -X utf8 experiments/contact_model/run_dynamic_closed_loop.py --scenes 1 --workers 2 --families S2 T2m-stiff --models memory_linear_checked memory_patch_checked_gain memory_patch_checked_target --label dynamic_contact_checked_smoke_reproduced
```

`dynamic_contact_replay_v1/v2` 保留错误目标/逆根饱和的失败实例，`v3` 修为前向力后才是当前方程。`dynamic_contact_closed_loop_smoke`、`smoke_v2` 是已修的回退/先验置信度实例；`smoke_v3` 六条和 `smoke_v3_completion` 两条合成完整两场景记录；`dynamic_contact_checked_smoke_v4` 是下一短窗检验后的六条。旧实例精确源码可从 Git 检查点恢复，不用新代码冒充逐位复现旧轨迹。v4 起记录源码哈希。所有同步求解结果均未把墙钟耗时计入对象延迟，不能声称实时或实际部署成功。

## 锚点论文痛点在本项目是否存在

后续按用户明确的研究顺序补做需求对应及对象通道验证，见专项研究 §6。`check_paper_painpoints.py` 用原控制器和同噪声做摩擦参数敏感性验证；其 `--normal-state-diagnostic` 模式只在评分路径分解弹性力与动态参考增益，隐状态不提供给候选。

```powershell
python -X utf8 experiments/contact_model/check_paper_painpoints.py
python -X utf8 experiments/contact_model/check_paper_painpoints.py --normal-state-diagnostic
```

结果是 `results/paper_painpoint_sensitivity.json` 与 `results/normal_contact_state_diagnostic.json`。前者检验当前仿真可表达的机制，后者定位静态逆模型的接触量混用；都不是新模型的闭环性能结果。当前不修改被控对象。

## 复现（仓库根目录，Windows）

```powershell
python -X utf8 experiments/contact_model/run_finite_patch_pilot.py --verify-only
python -X utf8 experiments/contact_model/audit_finite_patch_pilot.py
python -X utf8 experiments/contact_model/run_finite_patch_pilot.py --reuse --scenes 1 --window 16 --stride 16 --max-windows 0 --families S2 T3m --models uniform_p1 uniform_arrival patch_normal_arrival patch_wrench_arrival --label finite_patch_ablation_v4_reproduced
python -X utf8 experiments/contact_model/run_finite_patch_pilot.py --reuse --window 16 --stride 32 --max-windows 8 --families S2 T3m --models line_uniform_arrival line_patch_normal_arrival line_patch_wrench_arrival --label finite_patch_line_v5_reproduced
```

移除 `--reuse` 从头生成相同 V1 基线轨迹；若输入存在，`--reuse` 只使用已保存的测量。`--label` 防止覆盖主结果。结果 JSON 中的源码哈希记录生成时的代码；当前源码保留各旧实例的选项，旧源码的精确状态未单独归档，科学数值可按下列配置复现，耗时和 JSON 源码哈希不比较。

| 实例 | 公共原五模型以外的参数 |
|---|---|
| inverse_v0 | `--window 64 --stride 32 --max-windows 8 --objective inverse_unweighted --legacy-arrival --wrench-nuisance bias_tangent` |
| force64_v1 | `--window 64 --stride 32 --max-windows 8 --legacy-arrival --wrench-nuisance bias_tangent` |
| force16_v2 | `--window 16 --stride 32 --max-windows 8 --wrench-nuisance bias_tangent` |
| continuous_v3 | `--window 16 --stride 16 --max-windows 0 --families S2 T3m --models uniform_p1 patch_wrench patch_wrench_arrival --wrench-nuisance bias_tangent` |

“原五模型”需显式指定 `--models uniform_p1 uniform_p1.5 patch_normal patch_wrench patch_wrench_arrival`，并使用新的 `--label`。
