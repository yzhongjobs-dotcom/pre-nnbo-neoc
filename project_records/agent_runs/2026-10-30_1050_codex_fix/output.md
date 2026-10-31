已仅修改指定目录中的 [comparators.py](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/87afbb68-b141-417b-8521-15bc71118679/scratchpad/codex_out/comparators.py)，并新增 [space_patch.py](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/87afbb68-b141-417b-8521-15bc71118679/scratchpad/codex_out/space_patch.py)、[smoke.py](D:/pip_env/temp/claude/E--github-project-scan-impedance-update/87afbb68-b141-417b-8521-15bc71118679/scratchpad/codex_out/smoke.py)。未修改仓库文件。

修改及理由：

- **Rossi 两类**：用编码器和延迟指令在 1 ms 子步重建位移，施加与测力相同的低通，再进入相同带通和 EWRLS，修正处理域失配。保留固定 `A_E`，开放 `f_exc=2–30 Hz`；docstring 已说明与原文实验的差异。
- **DO-MPAC 两类**：改用 `q_over_r`、固定 `R=1`，开放 `Np∈{8,15,25,40}`，保留其余物理参数。搜索起始范围为 `q_over_r∈[10³,10⁸]`，逐组合检查名义极点，并仅保留衰减时间常数为 **0.04–2 s** 的组合。
- **WijType**：明确标为“带噪特权刚度 PI 代理”；开放 `probe_hz=20–125 Hz`，支持 `forget_mode='tanh'/'constant'`，常数模式的 `lam∈[0.8,0.999]`。
- **smoke**：优先导入副本，检查调参空间并运行默认控制器、静止估计前后对照；已用 `python -X utf8` 运行通过。

静止均匀硬组织保留噪声和黏弹性，设置 `p[2]=p[15]=p[5]=0`。在 `t≥2 s`，**median(k_hat/k_eff)：0.609163 → 0.967996**；修正前数值另用仓库原版核验一致。

冒烟结果如下；T2 使用允许的软组织版本。指标采用 `t≥1 s` 的真实接触力，RMSE、peak 单位为 N。

| 场景 | 方法 | RMSE | peak | n_loss |
|---|---|---:|---:|---:|
| T2，seed 11 | RossiType | 0.172444 | 4.859784 | 0 |
| | WijType | 0.336812 | 4.862816 | 0 |
| | DOMPAC | 4.008505 | 13.955754 | 1 |
| | RossiDOMPAC | 12.363125 | 24.036958 | 0 |
| T3，seed 12 | RossiType | 0.995905 | 8.917453 | 3 |
| | WijType | 1.669037 | 13.568034 | 2 |
| | DOMPAC | 7.597682 | 42.782130 | 3 |
| | RossiDOMPAC | 4.000000 | 0.000000 | 1 |

仍存在三点局限：

- 默认 `q_over_r=20` 保留旧增益，其名义时间常数约 **174 s**，已被新搜索空间排除。上述默认结果不能代替公平重调参；T3 的组合方法在整个计分段失去接触。
- 名义极点筛选不保证含 ESO、伺服、延迟和饱和的实际闭环稳定；诊断中 `q_over_r=10⁵` 也出现过发散。
- 原 `tune_equal.py` 的第二阶段 `perturb()` 尚未遵守新增边界，需要同步约束；已在 `space_patch.py` 注明。