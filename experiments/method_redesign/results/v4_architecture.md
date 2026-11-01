# V4 架构验证

本结果来自指定诊断集与剂量集；未调用 final 确认集。调参只用 190000+。

## 设计与实现

- **stiffness**: C1 repaired projection unchanged; a separate FINAL projection fit supplies exactly the existing k_good acceptance and motion measurement
- **drift**: nominal rollout only: delta_f(h)=weight*k_good*(v*h+a*(0.5*Tc*h+0.5*h*h)); 1 ms table, full horizon, no legacy tv clipping; no extra tau*sd because servo explicit
- **recovery**: C1 unchanged, IMM cleared and disabled while in C1 recovery/F<0.3; no additive external FF; C1 last-motion recovery rule retained
- **solver**: Original source defaults: H=16, iters=1 (not the 3 claimed in prose); kept identical for original C1 and V4 per task preserve-settings instruction
- **variables**: log D, KPF/(1e-4), uff_0/(1e-3), uff_slope/(1e-2); original warm start, ramp-then-hold, probes, residuals, bounds
- **ablation**: V4-noFF retains C1 polynomial drift, adds V1 only, and shares the selected V4 ru to avoid confounding the motion ablation with a different control-rate penalty; no motion observer work
- **timing**: Separate single-process pass, no physical real-time claim
- **information**: only t, encoder y, filtered F, command; no Ft/My or contact/scene truth; existing gain and direct interfaces retained

## 调参

预算 61/61；最优参数 `{'motion_weight': 0.3075269106650179, 'ru': 0.0004615657809345654}`；目标 2.416445。
默认+40随机+20局部；每族10；五族等权 mean log(1+min(J,1000))。所有试验分数见 v4_tuning.json。

原C1与Proposed参数按任务冻结，因此本轮是给定基线的架构试验，不能声称所有方法都在新调参集上做过同预算重新优化。

## 指标（RMSE N，中位/均值；r=失控，L=长失联）

|条件|V4|V4-noFF|C1|Proposed|Proposed+V1|
|---|---|---|---|---|---|
|diag/S2|0.148/0.159 r0 L0 d0|0.083/0.089 r0 L0 d0|0.083/0.089 r0 L0 d0|0.124/0.281 r0 L2 d0|0.126/0.271 r0 L1 d0|
|diag/T2f-stiff|5.735/8.873 r22 L45 d0|4.725/8.222 r19 L44 d0|4.176/7.872 r18 L43 d0|0.597/1.085 r0 L16 d0|0.621/1.092 r0 L17 d0|
|diag/T2m-soft|0.149/0.152 r0 L0 d0|0.084/0.087 r0 L0 d0|0.085/0.087 r0 L0 d0|0.077/0.078 r0 L0 d0|0.077/0.077 r0 L0 d0|
|diag/T2m-stiff|4.018/5.951 r8 L45 d0|2.849/5.014 r8 L31 d0|2.988/5.128 r9 L26 d0|0.283/2.305 r2 L1 d0|0.287/3.257 r2 L1 d0|
|diag/T3m|0.587/1.051 r0 L15 d0|0.185/0.227 r0 L0 d0|0.184/0.332 r0 L2 d0|0.151/0.157 r0 L0 d0|0.151/0.157 r0 L0 d0|
|sweep/ramp:0.1|6.758/9.301 r12 L16 d0|7.259/7.888 r12 L15 d0|6.152/6.744 r10 L17 d0|2.497/6.120 r5 L14 d0|2.511/6.309 r5 L15 d0|
|sweep/ramp:0.15|8.251/10.542 r15 L18 d0|6.204/8.980 r11 L18 d0|8.454/10.655 r13 L18 d0|2.421/14.497 r5 L17 d0|2.483/8.404 r6 L17 d0|
|sweep/ramp:0.2|5.973/10.369 r10 L16 d0|4.008/9.671 r9 L16 d0|6.231/9.165 r10 L17 d0|0.996/1.479 r1 L8 d0|1.011/1.481 r1 L8 d0|
|sweep/ramp:0.3|4.796/7.339 r6 L11 d0|4.949/7.059 r7 L11 d0|4.397/6.683 r6 L11 d0|0.302/0.474 r0 L0 d0|0.308/0.478 r0 L0 d0|
|sweep/ramp:0.45|3.357/5.430 r2 L13 d0|2.015/4.124 r2 L6 d0|2.114/4.168 r2 L6 d0|0.232/0.286 r0 L1 d0|0.237/0.290 r0 L1 d0|
|sweep/ramp:0.6|4.072/5.688 r3 L16 d0|2.735/4.697 r3 L12 d0|2.723/4.684 r3 L12 d0|0.243/0.336 r0 L1 d0|0.245/0.345 r0 L1 d0|
|sweep/ratio:10|0.098/0.111 r0 L0 d0|0.077/0.078 r0 L0 d0|0.077/0.079 r0 L0 d0|0.105/0.112 r0 L0 d0|0.106/0.113 r0 L0 d0|
|sweep/ratio:20|0.144/0.160 r0 L0 d0|0.086/0.091 r0 L0 d0|0.086/0.092 r0 L0 d0|0.136/0.275 r0 L0 d0|0.128/0.196 r0 L1 d0|
|sweep/ratio:3|0.063/0.065 r0 L0 d0|0.075/0.075 r0 L0 d0|0.075/0.076 r0 L0 d0|0.080/0.082 r0 L0 d0|0.081/0.082 r0 L0 d0|
|sweep/ratio:30|0.159/0.196 r0 L0 d0|0.085/0.098 r0 L0 d0|0.085/0.099 r0 L0 d0|0.135/1.122 r2 L1 d0|0.136/0.495 r0 L2 d0|
|sweep/ratio:45|0.178/0.212 r0 L0 d0|0.086/0.101 r0 L0 d0|0.086/0.100 r0 L0 d0|0.158/1.299 r2 L0 d0|0.148/0.602 r0 L1 d0|
|sweep/ratio:60|0.286/0.613 r0 L2 d0|0.108/0.132 r0 L0 d0|0.109/0.131 r0 L0 d0|0.254/4.468 r2 L6 d0|0.482/3.945 r3 L6 d0|

## 配对bootstrap：V4 − 对照，均值RMSE差（N）

|条件|Proposed：差 [95% CI]，有效配对数|C1：差 [95% CI]，有效配对数|
|---|---|---|
|diag/S2|-0.123 [-0.252, -0.020] n=50|0.069 [0.056, 0.083] n=50|
|diag/T2f-stiff|7.789 [6.031, 9.795] n=50|1.001 [-0.626, 2.877] n=50|
|diag/T2m-soft|0.074 [0.063, 0.086] n=50|0.065 [0.055, 0.075] n=50|
|diag/T2m-stiff|3.647 [0.341, 6.366] n=50|0.824 [0.260, 1.307] n=50|
|diag/T3m|0.894 [0.663, 1.125] n=50|0.719 [0.534, 0.920] n=50|
|sweep/ramp:0.1|3.181 [-2.407, 7.973] n=20|2.557 [0.546, 4.770] n=20|
|sweep/ramp:0.15|-3.955 [-28.974, 10.057] n=20|-0.113 [-1.415, 1.050] n=20|
|sweep/ramp:0.2|8.890 [5.919, 12.337] n=20|1.204 [-0.520, 3.353] n=20|
|sweep/ramp:0.3|6.865 [4.187, 10.092] n=20|0.656 [-0.327, 1.697] n=20|
|sweep/ramp:0.45|5.144 [3.156, 7.750] n=20|1.262 [0.639, 1.975] n=20|
|sweep/ramp:0.6|5.351 [3.708, 7.206] n=20|1.003 [0.465, 1.593] n=20|
|sweep/ratio:10|-0.001 [-0.013, 0.014] n=20|0.032 [0.021, 0.045] n=20|
|sweep/ratio:20|-0.115 [-0.373, 0.023] n=20|0.068 [0.050, 0.087] n=20|
|sweep/ratio:3|-0.017 [-0.022, -0.011] n=20|-0.011 [-0.015, -0.005] n=20|
|sweep/ratio:30|-0.927 [-2.196, -0.045] n=20|0.096 [0.070, 0.126] n=20|
|sweep/ratio:45|-1.087 [-2.320, 0.015] n=20|0.112 [0.084, 0.145] n=20|
|sweep/ratio:60|-3.855 [-10.832, 0.009] n=20|0.482 [0.172, 0.901] n=20|

## 单拍计算时间（单进程，ms）

|方法|P50|P99|
|---|---:|---:|
|V4|0.588|0.917|
|V4-noFF|0.142|0.361|
|C1|0.136|0.284|
|Proposed|0.385|0.653|
|Proposed+V1|0.368|0.598|

## 结论与下一环

未同时保住C1在S2的优势与Proposed在运动族的优势；不采纳当前V4拼接实例。

- S2：V4均值0.159 N，较C1 0.089 N差0.069 N，配对95% CI [0.056,0.083]。V4较Proposed均值0.281 N更好，但中位数0.148 N差于Proposed 0.124 N；不能只报均值称兼得。
- T3m、软推拉、硬推拉：V4均值分别较Proposed差0.894、0.074、3.647 N，配对区间均高于0。快速硬推拉差7.789 N，95% CI [6.031,9.795]；失控22/50、长失联45/50。
- V4-noFF的S2为0.083/0.089 N，与C1基本一致；V4变差出现在接入运动预测后的闭环链条，而非仅增加V1辨识。
- 剂量扫描：V4在ratio=3时优于C1；ratio>=10时没有保住C1精度。运动扫描各档中位数均大幅差于Proposed；ramp=0.15时Proposed均值14.497 N反而高于V4的10.542 N，说明其尾部也严重，不能宣称逐档均值全胜。

### 失败证据与机理边界

- 恢复/启动阻塞（直接观测）：硬推拉8/50、快速推拉5/50从未启用运动预测；29/50和40/50超过半数拍处于恢复。场景113026在0.064 s进入恢复，991/999记录拍处于恢复，运动预测从未启用。失联关闭FF的规则与原C1恢复组合后，运动信息无法帮助这些场景。
- 长时域预测失真（直接观测）：T3m、硬推拉、快速推拉的名义力漂移预测绝对值最大6481、5032、3486 N，分别12/50、37/50、36/50场景出现>50 N预测漂移。它们是模型预测数值，不是真实接触峰值。
- 离线已实现轨迹诊断：场景111015的32/128 ms等效力漂移误差中位数2.465/9.595 N；114030为2.984/12.338 N。该指标采用冻结真实切线刚度，包含非线性/刚度变化，不能当成纯组织运动误差。
- 估计不是唯一解释：114030的有效量测增益kg/ktrue中位1.031，越过历史筛选区间的比例仅2.2%，但128 ms预测误差仍大；S2两个审计场景存在9.2%/26.1%的增益越界，现有接受准则也不保证该闭环分布下始终有效。
- 机制推断：FINAL原本只用于执行延迟附近的短前馈预测；将同一CA后验外推到完整128 ms、且把极端力漂移作为确定模型交给GN，会驱动不合适的参考整形，并与长恢复相互强化。当前证据未通过oracle预测或恢复替换独立量化每一环的因果份额。

### 下一步（未在本轮实施）

- 先改恢复与启动接口：使没有k_good的早期失联也能可靠重建接触，再在可信接触后重启量测；不能靠继续调运动权重解决从未启用的场景。
- 随后改预测：为长时域建立可信范围/不确定性衰减，避免把短时CA外推原样视为128 ms确定漂移。用固定增益/已知运动的分环核验区分增益失配与时域失配。
- 再检验优化如何使用不确定预测与恢复状态；本轮没有证据把增加GN迭代或换优化器列为第一优先项。保留估计增益偏差问题，但不将所有失败归给估计。

### 执行与预算说明

- Windows multiprocessing Pipe creation denied; persistent process workers use temporary file IPC.
- At 39 saved trials an IPC PermissionError interrupted the next candidate; that same candidate was repeated after adding file-access retries. Exactly 61 distinct candidate configurations received completed objective scores; one interrupted attempt incurred extra simulation work, with no additional parameter selection opportunity.
- All metric computation uses the existing common.task_metrics warmup=1 s and true-force scoring convention.

## 边界与疑点

- Prespecified historical baseline parameters frozen, not independently retuned 61 times on new seeds. V4 alone gets <=61 new two-parameter trials as requested.
- Diagnostic and sweep seeds are historically used: no fresh-confirmation claim.
- Divergent RMSE undefined; finite-only moments and complete-finite-pair bootstrap explicitly report exclusions; divergence separately counted.
- Gain policy observes the prior logged row; direct policies consume the provisional current row as in their original implementations. Same sensor types, but existing timing/interface differences remain.
- Nominal plant evaluation only; V1 mismatch validation is inherited, not newly established here.

发散时RMSE与最长失联无定义，不补零、不剔除失败后宣称稳健；均值仅在有限结果上，发散数另列。Bootstrap每条件4000次，无多重比较校正。

复现：`python -X utf8 experiments/method_redesign/v4_evaluate.py --stage all --workers 12`。删除本任务新增的结果JSON后可从头运行；已有检查点默认续跑。
