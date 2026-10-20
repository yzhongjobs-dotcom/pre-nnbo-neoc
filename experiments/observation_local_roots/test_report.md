> 既有局部验证记录。过时选题与续做安排已移除，原配置和数值结论保留。后续安排见[HANDOFF](../../project_records/HANDOFF.md)。

# R15 小型可行性测试

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run + reproducibility check
- Origin Date: 2026-10-12
- Verification Status: VERIFIED，仅指固定种子的数值结果独立重跑完全一致，计时不要求一致
- Version Label: quick_probe_v1
- Type: simulation / algebraic feedback check
- Status: completed
- Python: 3.10.11；NumPy: 1.26.4；SciPy: 1.13.1
- 运行命令：`python -X utf8 run_quick_checks.py`，在本目录运行。报告为保留结果说明，旧汇总脚本未留在工作区。
- 最后一次脚本用时：5.65 秒，退出码 0。

## 目的和结论范围

用户授权用仿真或理论做一个简单可行性测试。本次测试两个局部机制，不重写论文，也不以新玩具模型代替 R15 的柱面机器人主线。

**结果支持空间力矩信息及小变化下直接接受候选的可行基础，但没有证明整套 R15 端到端加速、真实组织效果或具体平台 PublishOK。**

两个测试彼此独立，估计结果没有接入后面的控制根。接触测试采用正文模型内部的理想法向一阶矩，控制测试采用附录 M.3 的单状态例子。它们没有完整机器人、真实传感器、表面/偏置状态、能量账本或 MHE—控制闭环。

## 测试 A：力矩是否帮助提前分辨前方材料

来源为当前正文式（9）—（11）和附录 D.1。用跨区域半椭圆载荷的闭式积分合成法向力 F 和理想法向载荷一阶矩 m1，取 R=30 mm、L=20 mm、a=3 mm、E0*=30 kPa、Ef*=45 kPa、边界 b=0.5 mm。12 个扫描位置从 -2 mm 到 -0.2 mm，均在中心越界之前。

对同一历史和相同法向噪声，用加权有界非线性最小二乘分别拟合 F-history 和 F+m1-history。前者不是单点弱阈值基线。估计 Ef 与 b；其余几何、当前模量、表面和偏置已知。每种场景 100 次配对试验，随机种子固定。

噪声假设：nominal 的 σF=0.003 N、σm=2×10^-6 N·m；higher 的 σF=0.01 N、σm=8×10^-6 N·m。宽度失配场景的真值仍为 a=3 mm，拟合模型使用 1.05a。噪声和参数是本次人工设置，不声称符合真实传感器。

| 场景 | 方法 | 边界 RMSE（mm） | 前方模量 RMSE（%） |
|---|---|---:|---:|
| nominal_noise | F_history | 0.2105 | 5.801 |
| nominal_noise | F+M_history | 0.0484 | 1.068 |
| higher_noise | F_history | 0.9456 | 25.804 |
| higher_noise | F+M_history | 0.1823 | 4.644 |
| width_5percent_mismatch | F_history | 1.6006 | 41.073 |
| width_5percent_mismatch | F+M_history | 0.3230 | 2.790 |

所有 600 次拟合数值结束成功，不代表所有估计准确。标定噪声场景中，F+m1 的边界 RMSE 由约 0.210 mm 降为 0.048 mm。宽度仅错 5% 时，F+m1 的边界误差仍升至约 0.323 mm，主要是系统偏差。

**解释：力矩信息在这个理想空间模型中确有用，但模型/几何失配可能比随机噪声更重要。** 标准腕部俯仰力矩还含切向牵引力臂和坐标变换，本次没有验证能否在未知摩擦条件下可靠分离 m1，也没有验证所有参数的联合可辨识性。生成与拟合共用同一接触近似，因此匹配场景仅是模型内部检查。

## 测试 B：预测候选是否必须进一步 Newton

复用附录 M.3 的归一化例子：z_dot=u，u≥-0.7，T=1，运行代价 η²z²+u²，critic V=(1-t)Wz²。固定配置点 (t,z)=(0,1)，有限根为

```text
F(Z,η)=(-W+η²+u²+2Wu, 2u+2W-λ, λ(u+0.7))=0
Z=(W,u,λ)
```

锚点 ηa=0.8，非活动分支 W*=(-1+sqrt(1+4η²))/2、u*=-W*、λ*=0。ANEOC 预先计算根灵敏度；chord 缓存锚点 Jacobian 的 LU 分解；一次 Newton 使用新参数残差。

本例 F_Z 在旧根处不直接依赖新 η，因此 chord 与一次 Newton 数学相同。这不是所有问题中的等价结论。另加入按该结构消元的 structured_chord 和同样简化的 structured_ANEOC，避免靠小矩阵库调用开销制造加速。

### 共用解析筛选

仅在非活动参数分支、0≤W≤0.7、u=-W、λ=0 时接受。实际反馈恢复为 u(t)=-(1-t)Wz(t)，z(t)=exp[-W(t-t²/2)]，从 z0=1 出发始终满足动作下界。

标量残差 B=η²-W-W² 给出

```text
|W-W*| = |B|/(1+W+W*) ≤ |B|/(1+W)
|J(W)-J(W*)| ≤ L_J |B|/(1+W)
L_J=2U/3+2η²/3+2U²/15，U=0.7
```

最后的界来自对成本关于 W 的导数积分作显式上界。所有候选采用同一检查，预算 0.001。精确根只用于事后核对；不进入筛选算法。实际根误差与代价差均未超过这些界。

**这只是单状态问题相对于有限根反馈的解析误差与动作可行筛选，不是完整的连续 HJB 近最优认证。**

### 小变化和变大后的表现

每个半径使用同一组 2,000 个均匀随机描述符查询。

| 相对锚点变化范围 | 方法 | 直接通过 | 通过率 |
|---|---|---:|---:|
| ±0.02 | ANEOC | 2000/2000 | 100.00% |
| ±0.02 | cached_chord | 2000/2000 | 100.00% |
| ±0.02 | one_Newton | 2000/2000 | 100.00% |
| ±0.08 | ANEOC | 1842/2000 | 92.10% |
| ±0.08 | cached_chord | 1137/2000 | 56.85% |
| ±0.08 | one_Newton | 1137/2000 | 56.85% |
| ±0.20 | ANEOC | 687/2000 | 34.35% |
| ±0.20 | cached_chord | 431/2000 | 21.55% |
| ±0.20 | one_Newton | 431/2000 | 21.55% |

Δη=0.02 时，ANEOC 根误差约 8.27×10^-5，相对有限根反馈的实际代价差约 4.47×10^-6，解析代价上界约 7.54×10^-5，不用校正就通过 0.001 的预算。Δη=0.08 时解析上界约 0.001247，拒绝直接接受，即便实际代价差仍较小。这说明保守界会影响可用区间，拒绝不等于实际策略一定不好。

Δη=0.4 跨越该例的活动分支，当前筛选拒绝，不能继续沿用原分支保证。没有据此实现恢复或保证恢复成功。

### 同等筛选下的耗时

每次执行候选生成加同一个解析检查，3,000 个相同小变化查询，7 批交错顺序测量，表中为批均值的中位数（μs/查询）。同机两次运行如下。

| 实现 | 第一次 | 复现运行 |
|---|---:|---:|
| ANEOC | 6.886 | 7.493 |
| cached_chord | 19.335 | 20.912 |
| one_Newton | 22.082 | 23.550 |
| structured_ANEOC | 6.704 | 6.843 |
| structured_chord | 7.240 | 7.558 |

通用矩阵实现中预测较快；结构化 chord 消除了大部分差距。剩余微秒差受 Python、函数调用和桌面调度影响，未做显著性分析或实时最坏界，不能据此声称强基线下已证明加速。灵敏度增量准备成本另存在 results.json；公共初始根求解未计入，完整认证、状态恢复、MHE、通信与执行也未计入。

### 一个比预测误差更重要的发现

当 η=0.82 时，连续标量最优值可独立计算为 η tanh(η)，且对应最优动作满足限制。有限 critic 根反馈相对这个最优值的实际代价差约 **0.005133**，远大于本次 ANEOC 的额外预测误差 **0.00000447**。

因此，通过“预测误差预算 0.001”并不意味着总近最优误差也小于 0.001。即使重新精确求根，也仍有这份 critic/configuration 差距。参数较大时无约束连续解会先违反动作约束，脚本不再把它当作约束最优值比较。

## 实现与复现检查

- 解析灵敏度与闭式根中心差分核对通过。
- 96 点 Gauss 积分与独立自适应积分的成本核对通过。
- 候选的实际误差没有超过解析上界。
- 两个结构化实现与对应矩阵实现输出一致。
- 固定种子重跑的根误差、通过计数、感知汇总和设置完全一致。计时自然变化，保留两次记录。
- 初次开发检查发现“非活动配置根”不保证无约束连续最优动作可行，已修正比较范围；未把该失败输出作为结果。

## 适用边界

理想观测、小根预测与有限critic质量是不同层次；完整MHE、实际多维阻抗和整条控制路线尚未验证。本报告不设置后续研究优先级。

## 文件

- [运行脚本](run_quick_checks.py)、[汇总/复现检查脚本](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/summarize_results.py)。
- [最后一次结果](results.json)、[第一次完整运行结果](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/results_run1.json)、[运行日志](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/run_log.txt)。
- [根步长结果](root_steps.csv)、[筛选通过统计](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/acceptance.csv)。
- [全部感知试验](sensing_trials.csv)、[感知汇总](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/sensing_summary.csv)。
- [复现记录](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/%E5%BF%AB%E9%80%9F%E9%AA%8C%E8%AF%81_2026-10-12/reproducibility.json)。
