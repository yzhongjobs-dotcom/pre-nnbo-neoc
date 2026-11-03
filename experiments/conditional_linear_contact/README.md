# conditional_linear_contact

**结论：不采纳，详见 [研究文档](../../research/conditional_linear_contact.md)。** `results/` 中 Codex 的 V0 结果（`closed_loop_*`、`mismatch`、`open_loop`、`tuning_*`、`final_report.md`）有参数串扰：调参时 `p_stay`、`r_scale`、`noise` 同时改了运动 IMM，违反“运动 IMM 不变”的规格，不能直接引用。

验收后的诊断和迭代由以下脚本完成：`diagnose.py`（串扰、隔离、同输入回放）、`failure_trace.py`（失控前状态）、`bank_policy.py`、`bank_eval.py`、`bank_trace.py`（V2 滤波器组）。主要结果：`results/diagnose_*.json`、`results/bank_v2_first.json`、`results/bank_v2_probe.json`。

```powershell
python -X utf8 diagnose.py tracker
python -X utf8 diagnose.py isolated
python -X utf8 diagnose.py replay --scenes 12
python -X utf8 bank_eval.py --configs v1 lean gain_only noprobe slow.3 slow.3_fb --label reproduced
```


本目录是“以当前编码器位置为支点”的条件线性接触响应 V0 实验。它只新增实验代码和结果，不修改 `plant.py`、`common.py`、`ServoID`、已有实验结果或 `00_基线快照_2026-10-14/`。本轮不访问确认集；所有闭环诊断场景使用 `110000+`，调参使用 `190000+`，不会使用 `130000+`。

## 方法

状态为 `x=[k,D,r,a]`。对 8 个 1 ms 样本，`w_m` 与 `common._build` 相同，`tau_m=(m-1) ms` 是区间起点相对时间：

```text
z_j = F[j+1] - a_f^8 F[j]
z_j = [Y_j, W0, W1, W2] x_j + eps_j
W0 = sum(w_m), W1 = sum(w_m*tau_m), W2 = sum(w_m*tau_m^2/2)
```

`Y_j` 直接调用 `common._build(..., enc=True)` 的端点插值分支，并在每次 V1 伺服辨识后使用当前 `tau`。时间转移为

```text
k+ = k + w_k
D+ = D + k*Delta_y + r*dt + a*dt^2/2 + w_D
r+ = r + a*dt + w_r
a+ = a + w_a
```

`[D,r,a]` 使用连续白急动模型的精确离散协方差；其谱密度为 `q_F=k_hat^2*q`。`k` 的方差为 `(k_hat*q_lnk)^2*dt`，即对数刚度随机游走的相对尺度。量测噪声基础方差是 `0.04^2`，`r_scale` 只作为未建模相关性/近似误差的尺度项。

四个 IMM 模式是 `{q_smooth,q_maneuver} × {q_lnk0,q_lnk1}`。每一拍执行混合、四个并行线性 Joseph 更新、似然归一化、模式合并和时间传播；记录后验模式概率、NIS 和短时一致性平均。

V0 继承 `ServoID` 和 `MotionAware` 的控制接口，只替换 `_stiffness_est` 的反馈刚度估计来源。反馈刚度仍遵守：精确且一致时 `max(k_hat+c*sigma,(1-rho_d)*k_bar_prev)`，精确但不一致时只能上调，不精确时保持；一致性由 NIS 短时平均判断。`k_good` 只有后验相对标准差和一致性同时通过才更新，默认相对带宽为 `ln(1.2)`。前 1 s 仍逐值调用原 `MotionAware._stiffness_est`，新滤波器只并行记录；切换有一拍保持，避免跳到不可信值。`F<=.5` 时冻结 `k`、重置 `D/r/a` 并放大 `k` 协方差，重新接触后恢复。

控制器通过本目录中的 `MeasuredOnly` 包装，只暴露时间、编码器端点、滤波力、上一条命令和测量通道；真值列只用于离线评分。

## 明确排除的错误实现

此前一轮实现已丢弃，结果无效，不能作为负结果。本目录不复用它。必须避免：

1. `Y` 的参考点每拍变化时不变换扰动状态；本实现把参考点固定为当前编码器端点，并将 `D/r/a` 作为当前端点的状态，观测方程与状态转移一致。
2. 把 `z` 中 `D` 的系数写成 1；本实现使用滤波权重和 `W0`（由当前 `common._build` 的 `afc` 实际计算，约为 0.73–0.78，当前 `fc=30 Hz` 时为 0.77864），并计算 `W1,W2`。
3. 反馈刚度直接使用 `k_hat`；本实现保留保守裕量 `c*sigma` 和下降限速 `rho_d`。
4. 用 `|Y|` 阈值切模式或漏用 `q_k1`；本实现是四模式标准 IMM，刚度两模式均进入过程噪声和转移概率。
5. 用错误过程噪声离散化；本实现使用白急动连续模型给定的精确五次/四次/三次矩阵。
6. 只跑两个场景或不调参；脚本包含每族 2 场景的默认试通、61 次调参、五族诊断集、开环回放、失配和最多三轮迭代记录。

## 复现

在仓库根目录运行（Windows）：

```powershell
python -X utf8 experiments/conditional_linear_contact/smoke.py
python -X utf8 experiments/conditional_linear_contact/theory_checks.py
python -X utf8 experiments/conditional_linear_contact/tune.py --workers 14
python -X utf8 experiments/conditional_linear_contact/closed_loop.py --params results/tuned_params.json --workers 14
python -X utf8 experiments/conditional_linear_contact/open_loop.py --params results/tuned_params.json --workers 14
python -X utf8 experiments/conditional_linear_contact/mismatch.py --params results/tuned_params.json --workers 14
python -X utf8 experiments/conditional_linear_contact/report.py
```

脚本会逐行/逐场景写 JSON，已存在的同名结果默认续跑而不覆盖。调参候选为默认 1 次、随机 40 次、局部 20 次，共 61 次；每个候选在五族各 6 个 `190000+` 场景上计算 `mean(log(1+min(J,1000)))`。可调范围记录在 `tune.py` 和 `results/tuning_protocol.json`，包括 `q_lnk0/q_lnk1`、两个 IMM 转移概率、噪声尺度、`k_good` 接受阈值、一致性阈值、`c` 和 `rho_d`。

## 结果解释边界

诊断集只用于方法诊断与迭代，不是确认集。计时是单进程 Python 控制调用时间，不含通信，不构成实时认证。噪声为零孪生只用于分解结构误差/噪声误差，不是额外的真实性保证。若某轮不理想，`results/iteration_*.json` 和 `results/iteration_notes.md` 必须区分原理、设计、实现、参数和接口原因，不把本轮的错误实现算作负结果。
