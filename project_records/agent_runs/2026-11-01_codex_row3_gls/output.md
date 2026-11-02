已完成，**不建议直接替换现有投影**。主目录未改动，未提交或推送，未使用 final 确认集。

交付：[补丁](<D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/row3_gls_delivery/row3_gls.patch>) · [代码与结果包](<D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/row3_gls_delivery/row3_gls.zip>) · [完整报告](<D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/row3_gls_delivery/experiments/method_redesign/results/row3_gls.md>) · [JSON](<D:/pip_env/temp/claude/E--github-project-scan-impedance-update/03dd8b7f-853e-4b3c-8a7b-88d5bd165f74/scratchpad/row3_gls_delivery/experiments/method_redesign/results/row3_gls.json>)。共6个新增文件；Windows应用补丁请使用 `git -c core.autocrlf=false apply`，以保持源码哈希一致。

开环总体结果，每族20场景；误差为 \(|\ln(\hat k/k)|\)。

| 方法 | 中位/P95 | 危险偏低 | 危险偏高 | 2σ覆盖 |
|---|---:|---:|---:|---:|
| E0 | 0.104/0.418 | 10.1% | 6.1% | 93.3% |
| E2 | 0.102/0.408 | 9.7% | 5.8% | 96.1% |
| E4 | 0.108/0.393 | 9.8% | 6.8% | 94.6% |

E4改善硬组织推拉的尾部：中速/快速族P95从 **0.418/0.517降至0.283/0.356**，但总体中位误差略升。

闭环诊断，每族50场景；RMSE为中位/均值，单位N；计数为Proposed/E4。

| 族 | Proposed | E4 | 失控 | 长失联 | 配对均值差［95%区间］ |
|---|---:|---:|---:|---:|---:|
| S2 | 0.124/0.281 | 0.369/0.800 | 0/1 | 2/3 | +0.519［+0.240,+0.951］ |
| 混合扫描 | 0.151/0.157 | 0.316/0.380 | 0/0 | 0/1 | +0.223［+0.159,+0.300］ |
| 软组织推拉 | 0.077/0.078 | 0.084/0.087 | 0/0 | 0/0 | +0.010［+0.008,+0.012］ |
| 硬组织中速 | 0.283/2.305 | 0.275/1.393 | 2/1 | 1/1 | −0.912［−2.794,+0.048］ |
| 硬组织快速 | 0.597/1.085 | 0.630/1.179 | 0/2 | 16/15 | +0.094［−0.094,+0.386］ |

选定超参：**ℓ=120 ms，标准差比例a=3，\(v=(3F_w)^2\)**；使用12次调参评估，接受阈值固定为 \(2\sigma/\hat k\le\ln1.2\)。

主要限制是**后验门不够可靠**：开环接受率仅6.0%，接受子集的2σ覆盖只有80.2%。GLS缓解部分漂移泄漏，同时改变更新时机与保守反馈；本轮未分离这些作用。sweep仅最慢0.6 s推拉档明确小幅改善，刚度比3、10、20倍均变差。代表场景串行E4完整控制P99为0.52–0.66 ms，属于软件计时。

代数退化、真值/未来泄漏检查、E0/E2复现、配对统计重算及代表场景串行复跑均通过。