已完成实现并仅冒烟验证，未提交、未推送、未运行全量矩阵。

文件：

- [physics_plant2.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/physics_plant2.py)
- [scenes2.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/scenes2.py)
- [contact_candidates.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/contact_candidates.py)
- [check_plant2.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/check_plant2.py)
- [test_candidates.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/test_candidates.py)
- [contact_select.py](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/contact_select.py)

结果写入：

- [plant2_checks.json](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/plant2_checks.json)
- [candidate_checks.json](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/candidate_checks.json)
- [contact_select.json](E:/github_project/scan-impedance-update/experiments/rossi_failure_study/results/contact_select.json)

实际偏差：

- 新场景均使用 `xcmd`；为兼容旧场景，`physics_plant2.py` 保留了缺失 `xcmd` 时由 `sc['v']` 合成轨迹的回退，不读取 `v_scan`。
- `--smoke` 使用每工况两个场景，共 24 场景；全量 `matrix` 未运行。
- 低刚度 stick-slip 检查用近似停滞比例辅助判定；指定的平滑 `tanh` 摩擦式下未出现速度符号跳变。
- 其余规格未发现实现偏差。

`check_plant2.py` 完整输出：

```text
zero-new-physics replay max diff: col3 1.77635684e-15, col5 8.8817842e-16, col8 7.55790097e-10
REGRESSION PASS: all requested differences < 1e-6
uniform sliding |Ft|/Fn median 0.139663; set mu_k 0.139954; n=374
reversal Fx truth jumps: median 1.24823, expected 2 mu F median 1.38175, n=1
low-Kt stick-slip: vc range [0.000284627, 0.0115635], |vc| P10/P90 0.000396255/0.0100305, near-stick fraction 0.412, sign changes 0, mu_s/mu_k 1.086
F_axis exact decomposition residual: median 2.77556e-17, max abs 4.44089e-16; small-angle residual median -0.0400131
v=.03 pressure-asymmetry F(v)/F(0): observed 1.01422801, eq(14) 1.01347218, ratio 1.00074578
2 kHz finite-state check: True rows 500
```

`test_candidates.py` 完整输出：

```text
L vs contact_compare elastic command max gap: 3.49436192e-14
L regression PASS
W phi=0 -> L max numerical gap: 0
S phi=0 -> L max numerical gap: 0
H phi=0 -> L max numerical gap: 0
phi=0 degeneration PASS
NaN truth-column mask command max gap: 0
truth masking PASS
L   phi mem 0.12216/0.07746 rate 0.06215/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23358 | call 2.1596s
W   phi mem 0.13386/0.07746 rate 0.04478/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23358 | call 2.3190s
S   phi mem 0.10262/0.07746 rate 0.02993/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23358 | call 2.2324s
H   phi mem 0.13733/0.07746 rate 0.04997/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23358 | call 2.3301s
Hz  phi mem 0.12379/0.07746 rate 0.06019/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23358 | call 2.1996s
AS  phi mem 0.12211/0.07746 rate 0.06224/0.07746 | mu 0.00000/0.28359 kd 1.00000/12.17384 | c 0.48167 target 0.23358 | call 2.2434s
TR  phi mem 0.12379/0.07746 rate 0.06019/0.07746 | mu 0.00000/0.28359 kd 0.99359/12.17384 | c 0.48169 target 0.23358 | call 2.1197s
LF  phi mem 0.12216/0.07746 rate 0.06215/0.07746 | mu 0.00000/0.28359 kd 1.00000/12.17384 | c 0.48167 target 0.23358 | call 2.1368s
WF  phi mem 0.12930/0.07746 rate 0.07274/0.07746 | mu 0.00000/0.28359 kd 1.00000/12.17384 | c 0.48167 target 0.23359 | call 2.1025s
WA  phi mem 0.12913/0.07746 rate 0.07283/0.07746 | mu 0.00000/0.28359 kd 1.00000/12.17384 | c 0.48167 target 0.23359 | call 2.2005s
WT  phi mem 0.12755/0.07746 rate 0.07349/0.07746 | mu 0.00000/0.28359 kd 0.99451/12.17384 | c 0.48169 target 0.23358 | call 2.1391s
P5  phi mem 0.10100/0.07746 rate 0.03777/0.07746 | mu 0.25000/0.28359 kd 1.00000/12.17384 | c 0.00000 target 0.23100 | call 2.0503s
```

单次控制调用耗时中位数 / P99（μs）：

```text
L 518 / 923
W 611 / 1050
S 585 / 975
H 569 / 1008
Hz 534 / 937
AS 547 / 1077
TR 552 / 1080
LF 517 / 1109
WF 551 / 1079
WA 580 / 1062
WT 574 / 1132
P5 473 / 666
```

`contact_select --smoke`：24 场景×12 候选全部返回日志，失控均为 0；`stop_go` 工况的 S 候选有 1 次长失联，其余为 0。全量比较留给 Claude 验收后运行。