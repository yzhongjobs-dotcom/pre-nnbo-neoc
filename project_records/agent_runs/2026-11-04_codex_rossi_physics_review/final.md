审查范围：5 个目标文件逐行检查；只读核对旧 `cobot_scan_benchmark/plant.py`、`motion_aware_scan/proposed.py`，并运行了短小数值测试。没有修改、创建或删除仓库文件，也没有提交。

## 一、会改变结论的错误

### 1. 第 8 列不是真实的“编码器到力”的探测频率增益

位置：`physics_plant.py:111–119、133–136`

第 8 列实际计算的是：

```python
ke × Prony 实部频响 × kr/(ke×Prony实部频响 + kr)
```

但实际闭环中的编码器到力传递还包含：

- 身体支撑动力学：`physics_plant.py:77`
- 身体受力反作用；
- 接触阻尼项：`physics_plant.py:90`
- 身体动态柔度与相位。

因此它不是实际 `F/y_enc` 的频率增益。短测试结果：

| 频率 | 实际 `F/y` 同相增益 | 第 8 列 |
|---|---:|---:|
| 1 Hz | 约 3166 N/m | 约 4412 N/m |
| 5 Hz | 约 3695 N/m | 约 4412 N/m |
| 20.8 Hz | 约 4480 N/m | 约 4412 N/m |
| 40 Hz | 约 4711 N/m | 约 4412 N/m |

20.8 Hz 尚有约 1.5% 偏差，但 1 Hz 已约 40% 偏差。现有 `failure_replay.py:120–124` 用第 8 列作为 `k` 的评分真值，而 Rossi/投影估计器实际观测的是测得力与编码器位移，因此评分目标错位，会直接改变“估计失败原因”和方法排名。

另外，Prony 频响只取实部；若“频率刚度”意指幅值或复频响，则还遗漏相位。若只想评分接触材料切线，应明确改成“接触刚度真值”，不能把它当作编码器—力闭环增益。

### 2. `failure_replay.py` 没有真正构造“有探测/无探测”的双胞胎

位置：`failure_replay.py:65–89`

代码只做了一次：

```python
log = PP.episode(...原闭环...)
clean = PP.replay(sc, log[:, 11])
```

随后 `no_slow` 和 `ideal` 都是对同一条 clean 信号离线滤波或重构，并没有：

- 去掉探测后的命令序列；
- 第二次无探测物理 replay；
- `F_probe = F_with_probe - F_without_probe`；
- `y_probe = y_with_probe - y_without_probe`。

`cum` 仅作为 `ref` 传给估计器，未参与物理 replay。因此文件头部关于“with and without probe”“probe own response”的描述不成立。由此得到的四类误差不是物理因果分解，只是同一轨迹上的信号代数分解。

### 3. “tracking”反事实把真值刚度直接送入估计器输入

位置：`failure_replay.py:85–89`

```python
kc = log[:, 8]
Fideal = kc * (...)
```

这里 `kc` 是物理真值第 8 列，随后构造的 `Fideal` 被送入：

```python
est = {k: ES.run_named(v, NAMES, ref) for k, v in sets.items()}
```

因此估计器看到的 `ideal` 力通道已经包含真值刚度。文件头声称“Truth is read only for scoring”不正确。`tracking=id_[ok]` 也不只是 tracking，而是“在真值刚度重构力、零相位滤波、非物理命令通道下的残差”。

### 4. 离线零相位滤波污染了“因果失败归因”

位置：`failure_replay.py:82–88`

```python
Fh = sg.sosfiltfilt(...)
yh = sg.sosfiltfilt(...)
ch = sg.sosfiltfilt(...)
```

`sosfiltfilt` 使用未来数据。虽然它确实只出现在离线分析中，但结果随后又被送入因果估计器，尤其：

- `no_slow` 的输入已含未来信息；
- `ch` 被作为命令重新进入 `matched_displacement`；
- 该命令不再是实际执行过的命令，且与 `yh` 的伺服状态不一致。

所以不能据此声称“因果带通/多项式方法的误差分别来自 slow leak、in-band、tracking、noise”。

### 5. 场景实际上没有测试 20.8 Hz 附近的真实未知运动

位置：`scenes.py:123–132`，调用位置：`closed_loop.py:27`、`failure_replay.py:140–143`

`tremor` 频率为 8–12 Hz：

```python
f0 = rng.uniform(8., 12.)
```

而 Rossi 探测频率是 20.8 Hz。它属于明显的带外扰动，带通方法天然受益。

同时，工具 ripple 只有在传入 `overrides['ripple_amp']` 时才生成：

```python
if o.get('ripple_amp'):
```

当前 `closed_loop.py` 和 `failure_replay.py` 都没有传入该参数，因此 `r(t)` 全部为零。也就是说，场景没有实际提供最直接的同频未知干扰，却用它判断 Rossi 频率分离在真实扫描中是否失效，结论不能外推。

## 二、可能影响结论的问题

### 6. 物理积分是显式近似，不是严格串联耦合求解

位置：`physics_plant.py:77、84–98`

身体更新使用上一子步的接触力：

```python
... + fx - s[12]
```

而当前子步的 `F` 到 `90–102` 行才求出。这样身体反作用有一个 0.5 ms 数值滞后。静态符号和平衡关系是正确的；短测试中 `-F/k` 身体下沉也复现正确。但在快速 flinch/cough 或强阻尼瞬态中，这会改变高频相位。

牛顿斜率：

```python
slope = 1. + (kt_guess + lam/dt*F)/kr
```

没有完整包含

```text
dF/dd = d(fe+q)/dd · (1+λ d') + (fe+q)λ/dt
```

而且 `kt_guess` 还是上一次日志更新的近似值。三次迭代在当前样例中残差很小，但它不是严格的串联柔度牛顿迭代。

### 7. Prony 递推是合理的中点近似，但不是精确递推

位置：`physics_plant.py:62、88、96`

```python
q_new = q_old*exp(-dt/tau) + alpha*(fe-fe_prev)*exp(-dt/(2*tau))
```

符号和松弛方向正确，恒定 `fe` 时会正确衰减；但对于线性变化的 `fe`，精确积分系数应为：

```text
alpha * (tau/dt) * (1-exp(-dt/tau))
```

当前中点系数只是近似。它不一定推翻当前结论，但应在文档中称为离散近似，而不是未说明的 Prony 真值实现。

### 8. `matched_displacement` 只在默认平台上近似匹配新物理 plant

位置：`estimators.py:27–35`

估计器使用：

```python
tau=.019, ndel=1
yy = tgt + (yy-tgt)*exp(-.001/tau)
```

而新 plant 使用 0.5 ms 显式 Euler：

```python
physics_plant.py:123–124
s[2] += dt/tau_a*(buf[0]-s[2])
```

默认 PC125 下误差较小，但它不是严格同一离散系统；若平台滞后或延迟改变，`estimators.py` 的固定参数会失配。`closed_loop.py` 中的在线 `_MatchedDisplacement` 能从平台读取参数，离线 `estimators.py` 则不能。

另外，`run_grid()` 中：

```python
rf = matched_displacement(ref, ref)
```

把参考探测命令同时当成“编码器轨迹”和“命令”处理，并非实际 probe encoder 响应。它只影响 IV 变体，当前主 `NAMES` 未使用 IV。

### 9. BandPassProposed 并非完全等价于原 `MotionAware._stiffness_est`

位置：`closed_loop.py:65–92`；对照 `motion_aware_scan/proposed.py:79–140`

当前默认参数下，主要接受规则大致相同；但实现省略或硬编码了：

- `kgood_rule`；
- `conservative` 分支；
- `rel_meas`；
- `meas_better`；
- `kg_rel/f_set` 的完整更新；
- `dof = nw-degree-2`，改成了硬编码 `11`。

因此它只能说是“当前 FINAL 参数下的近似替代”，不能说两者规则完全对应。

### 10. Band-pass 的噪声模型低估了新 plant 的传感器噪声

位置：`physics_plant.py:170`、`closed_loop.py:72`、`proposed.py:55`

新 plant 固定使用：

```python
sig = .08
```

而 `MotionAware` 默认 `noise=.04`。BandPassProposed 使用：

```python
noise2 = p['noise']**2 * self.est.gain
```

因此不确定度约被低估为实际标准差的一半，可能使 `precise`、`consistent` 和 `k_good` 过早通过。V1/Recording 也继承同一 `.04` 假设。该问题会差异化影响带通和投影方法的闭环行为。

### 11. Oracle 不是严格的“任何估计器上限”

位置：`closed_loop.py:130–131、261–263`；`proposed.py:74–76`；`variants.py:119–120`

该 oracle 使用第 8 列真值刚度作为反馈和运动观测参考，但：

- 第 8 列本身不是完整 `F/y_enc` 动态增益；
- 没有使用真实未来表面；
- `oracle_base` 直接采用 V1 的调参结果，而不是 oracle 自己调参。

因此它是“给 Proposed 结构注入真值刚度的参考”，不是严格上界。用于估计器诊断可以，但不能写成任何估计器可达到的上限。

### 12. 调参次数相同，但物理激励预算并不公平

位置：`closed_loop.py:139–144、200–219`

各方法确实都是 61 次评估，训练/测试种子也分离；但参数空间不是同一物理预算：

- Rossi：`A_E=5e-6–2e-4 m`；
- Proposed：`a_force=.05–.5 N`，再除以估计刚度；
- Rossi 允许 `f_exc=12–30 Hz`；
- Proposed 带通固定 20.8 Hz；
- `rossi_bp` 还预置了一个 `20.8 Hz、half=4、order=2` 的人工起点：`closed_loop.py:206–207`。

所以“相同 61 次”成立，但“相同探测力/位移代价”不成立，结果可能受调参范围影响。

### 13. 非标称 lag 评估时，测量包装器仍按标称平台 reset

位置：`closed_loop.py:170–171、189–194`

```python
wrap = pol if name == 'oracle' else _Measured(pol, plat)
```

但：

```python
_Measured.reset:
    self.inner.reset(dict(PP.PLATFORM))
```

忽略了传入的 `plat`。因此 `plat_scale != 1` 时，plant 使用改变后的滞后，控制器仍使用标称滞后。当前主评估 `scale=1`，所以不影响主结果；但会使 lag 鲁棒性测试无效。

### 14. 误差组不是互斥组

位置：`failure_replay.py:94–102`

`change` 可以同时属于 `event`，而且 `.05 ≤ dln < .15` 的样本不进入 `steady` 或 `change`。这适合描述性分组，但不适合把各组统计解释成完整、互斥的误差来源。

另外，`noise = noisy-clean` 只表示固定原始命令轨迹下的输入噪声影响；没有在无噪声条件下重新闭环，因此不包含噪声通过控制器改变命令后的二次影响。

## 三、次要问题

- `physics_plant.py:117` 的 `E_under` 是中心位置的 `exp(ln_e)`，不是有限足迹下的平均或等效模量；只是诊断列时问题较小。
- `physics_plant.py:58` 的 `direct_log` 参数未使用。
- `scenes.py:27` 的 `a0` 未使用。
- `combo` 只随机抽取 1–2 个事件，并不保证 slow/fast/flinch/cough 都出现：`scenes.py:90–96`。
- 场景随机数对每个 cell 基本独立，且方法间正确共享同一场景和同一传感器噪声；但 `base+i` 使不同 cell 间的传感器噪声按索引复用，若将所有 cell 当完全独立样本计算置信区间，会略显乐观。
- `estimators.py` 的 FIR 设计没有发现问题：零直流、20.8 Hz 归一化和 `lfilter` 因果性均通过短测试；Rossi 折扣和式确实等价于未归一化 EWRLS；`projection`、`local_linear` 没有发现未来数据泄漏。
- RossiLoop 的方向和离散单位基本正确：`closed_loop.py:113` 是 `K_I/k_hat * (FD-F) * Tc` 加上正弦增量，符合原式 (11) 的积分力误差加探测运动思想。主要问题在于其输入估计和真值/场景设计，而不是这一行的符号。

## 总结

最需要阻止当前结论发布的是：

1. 第 8 列真值没有包含身体动态反作用，却被用于估计误差评分和 oracle；
2. 场景没有真正产生 20.8 Hz 附近的未知扰动，带通方法受到结构性有利条件；
3. `failure_replay.py` 没有有/无探测物理双胞胎，并把真值刚度用于构造估计器输入；
4. 零相位处理和固定命令 replay 不能支持闭环因果归因；
5. 新 plant 的 `.08 N` 噪声与控制器使用的 `.04 N` 不一致。

因此，当前代码可以支持“若干信号处理实现和默认 replay 时序可运行”的判断，但不能可靠支持“Rossi 频率分离在真实超声滑动扫描中已经被证明失效、且失效原因已由四类误差分解确定”的结论。