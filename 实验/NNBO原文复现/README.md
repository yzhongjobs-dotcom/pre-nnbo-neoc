> 既有局部验证记录。过时选题与续做安排已移除，原配置和数值结论保留。后续安排见[HANDOFF](../../项目记录/HANDOFF.md)。

# NNBO 原文最小机制复现

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent；PDF 阅读；nature-figure（Python 绘图脚本）
- Origin Mode: run + source-equation checks + reproducibility
- Origin Date: 2026-10-12
- Verification Status: REPRODUCED，四组原始轨迹与指标重跑完全一致；不表示精确复现所有原图或证明 HJB 最优
- Version Label: nnbo_source_minimal_v1
- Status: PASS_MINIMAL_MECHANISM_NOT_EXACT_FIGURE_REPLICATION
- Python 3.10.11，NumPy 1.26.4，SciPy 1.13.1
- 最近主运行用时 16.54 秒，退出码 0。

## 1. 已跑通什么

已实现并执行原文的一维 Hunt–Crossley 力输入增广系统、三项 single critic、在线权重更新和原文探索信号。两个代价权重都完成 0—7 s 数值交互，未丢失接触，状态与权重未出现非有限值；更高位置误差权重产生更靠近目标位置、同时更大的接触力。

另外，对记录的力与 critic 做了原文式（52）Newton Downhill 参考位置反解，逐点用独立括区间求根核对。反解是独立执行接口检查，**没有将它接成“测量力→参考位置→理想位置机器人”的数字闭环**。

因此本次完成的是 **原文 §III-A—C/§IV 连续增广系统的最小机制仿真及 §III-D 反解检查**。没有复现参数辨识、实机、文献 [44]/[45] 对照或全部原图，不能称为整篇论文精确复现。它也不是当前 R15 五状态接触 HJB–KKT 求根的复现。

## 2. 原文对应与参数

原始论文为 Kong et al., *Neural-Network-Based Optimal Impedance Control for Robots in Physical Interaction With Soft Environments*, IEEE TSMC: Systems 55(10), 6463—6475 (2025)，DOI 10.1109/TSMC.2025.3579017。

本项目原文：[NNBO PDF](../../文献/原文/NNBO.pdf)。

PDF及原运行代码SHA256见results.json；提取文本位于文献目录。重复公式截图已清理，核对时直接读PDF。

| 原文位置 | 本次实现 |
|---|---|
| §II-D，式（21）与 §IV | r_dot=-r+0.03，r0=0.08，xd=r |
| §III-A，式（23）—（25） | x_dot=-Ke/Be+f/(Be(x-xe)^n)+ε；X=(e,r)，e=x-r |
| §III-C，式（35）—（37） | Jhat=W1 e²+W2 er+W3 r²；fhat=-(2W1 e+W2 r)/(2RBe(x-xe)^n) |
| 式（38）—（41） | Hhat=W^Tν+Qe²+Rfhat²+ζ²；W_dot=-κνHhat/(1+ν^Tν)² |
| §IV（PDF p7—8） | Ke=1149，Be=567，n=1.2，xe=0.01 m，x0=0.011 m；W0=(20,20,20)，κ=5；Q=500000/1000000，R=1 |
| §IV 探索与扰动 | 按原文完整三角探索信号加入力输入；ε=1000 e sin(r) |
| §III-D，式（52）/Algorithm 1 | 独立执行参考位置 Newton Downhill，固定测量力、r、W，并与 brentq 核对 |

ν 使用原文辅助模型 F+Gfhat，物理增广系统使用探索后的 fhat+N 和声明扰动；不把真实扰动并入名义 critic 动力学。

## 3. 未报告设置与复现假设

1. **ζ 的具体形式原文未明确给出。** 主运行采用已声明仿真扰动的紧上界 |1000e sin(r)|；另运行保守上界 1000|e|，两者均满足本例 |ε|≤ζ。两组结果全部保留，不按图形相似程度选优。
2. **探索终止时刻原文未明确。** 本次整个 7 s 持续加入原文信号，没有自行截断或缩小幅值。末端仍有探索导致的波动，不称为严格稳态。
3. 7 s 取自 Fig.2—3 的显示范围，输出间隔按原文 0.001 s。采用连续学习律的自适应 RK45 内部积分，最大内部步长 0.001 s；这不同于将权重每 1 ms 显式 Euler 更新并保持不变。
4. 初始浅接触导致强刚性。开发时固定 1 ms RK4 出现内部接触域越界；自适应积分解决数值步长问题，没有裁剪力、权重、位置或修改纸面学习律。独立 Radau 核对通过。
5. 直接积分原文 force-as-input 增广系统是 §III-A—C 的闭环数值实例，不宣称已经实现原文的机器人内位置环。该区别在后续 R15 迁移时不能省略。

## 4. 实测结果

表中为 t=7 s 的实际输出；applied force 包含探索。该时刻 xd≈30.046 mm，组织表面为 10 mm。

| Q | 扰动上界假设 | 最终 x（mm） | 最终施加力（N） | 峰值力（N） | 最终 Hhat |
|---:|---|---:|---:|---:|---:|
| 500000 | tight_known_disturbance | 22.468 | 6.636 | 14.330 | 62.938 |
| 1000000 | tight_known_disturbance | 24.708 | 7.849 | 18.126 | 79.754 |
| 500000 | uniform_sine_bound | 26.063 | 8.562 | 20.619 | 87.214 |
| 1000000 | uniform_sine_bound | 27.047 | 9.068 | 22.529 | 91.160 |

主假设下 Q 增大时 x 从约 22.47 mm 增到 24.71 mm，施加力从约 6.64 N 增到 7.85 N，与原文所描述的位置精度—作用力权衡方向一致。它只支持定性机制一致性，不代表数值曲线完全吻合原图。

**重要限制：Hamiltonian 残差并未趋近零，权重在最后一秒仍变化。** 运行没有崩溃、交互轨迹合理，不等于 critic 已找到精确 HJB 解或满足 R15 的高精度根门槛。上界选择明显改变策略，说明未报告设置会影响精确复现。

## 5. 数值和接口检查

- 标量实现的策略/学习律与独立矩阵表达逐项核对通过。
- HC 名义力加扰动项与施加力的最大恒等误差约 7.11×10^-15 N。
- 四组仿真均完成 7001 个输出点，压入量始终为正；未添加人工饱和来通过检查。
- 4×701 次参考位置反解均与独立 brentq 一致，方程残差小于 10^-10 N，最多 5 次 Newton 迭代。
- 独立 Radau 相对 RK45 的最大位置差为 2.58e-10 m，最大权重差为 2.51e-07。
- 四组原始 t/y/data 数组及指标重跑逐项完全一致，耗时不要求一致。
- 记录了有限 7 s 的激励 Gram 矩阵，不将这一个区间的正特征值冒充对任意窗口成立的 PE 证明。

## 6. 成本、最优性及与 R15 的关系

CSV 分别保存运行代价 Qe²+Rfhat²、0—t 的有限累计代价、critic 价值估计及 Hhat。原文 Fig.2(d)—3(d) 展示下降的“cost”曲线，但非负运行代价的前向累计积分不应下降；原文缺少足够实现细节使这些曲线可以逐值重建，本次不将任何一种量改名以制造吻合。

由于末端仍有非零跟踪误差/作用力，本次不声称其无限时域累计成本有限或已经求得无限时域最优解。原文关于权重有界的论证与 PE 假设也未由本次有限仿真独立证明。

这次成功建立了可运行的原始 single-critic 学习机制参照。此前 R15 的失败发生在另一套固定配置点根的初始求解，因此 **原文仿真跑通不会自动证明 R15 根可解**。原学习机制与B0固定根的状态、动作和表达不同，不能直接迁移其三维权重。

## 7. 文件与重跑

- [运行脚本](run_reproduction.py)、[冻结配置](config.json)、[完整结果](results.json)。
- 当前保留四个 `Q*_*.npz` 的完整轨迹及results.json；旧CSV副本和反解表可从Git回查或由driver重新生成。
- [复现检查](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/reproducibility.json)、[第一次运行结果](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/run1/results.json)、[运行日志](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/run_log.txt)。
- [绘图脚本](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/plot_diagnostics.py)：当前 Python 未安装 Matplotlib，安装下载未完成，未生成或宣称已检查诊断图。数值运行和复现检查不依赖它；完整轨迹已保存为 CSV/NPZ。
- [汇总与重跑比对脚本](https://github.com/zyanzhi90-dot/pre-nnbo-neoc/blob/b3046ce04a99d1896beaf6b7d57d4c10cb85dc5c/03_%E5%AE%9E%E9%AA%8C%E4%B8%8E%E5%90%8E%E7%BB%AD%E5%B7%A5%E4%BD%9C/NNBO%E5%8E%9F%E6%96%87%E6%9C%80%E5%B0%8F%E5%A4%8D%E7%8E%B0_2026-10-12/summarize_reproduction.py)。

在本目录执行：

```powershell
python -X utf8 run_reproduction.py
```

首轮重复副本已移出工作区，历史结果仍在Git。数值运行需要NumPy/SciPy；重跑请使用独立副本或新输出位置。
