from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent


def read(path):return json.loads((HERE/path).read_text(encoding='utf-8'))


def main():
    nonlinear=read('01_复用动作计算灵敏度/results.json')
    smooth=read('02_平滑二次critic对照/results.json')
    dense=read('03_小邻域与密集更新_探索性/results.json')
    verified=read('reproducibility.json')
    assert nonlinear['status']=='COMPLETED_DIAGNOSTIC_NOT_CERTIFIED'
    assert verified['status']=='VERIFIED_STATED_SUBSETS'
    dn=next(row for row in dense['scores'] if row['method']=='local_NEOC')
    dr=next(row for row in dense['scores'] if row['method']=='full_solve')
    rows=[]
    for case in nonlinear['scenarios']:
        local=next(row for row in case['scores'] if row['method']=='local_NEOC')
        cont=next(row for row in case['scores'] if row['method']=='continuous_NEOC')
        rows.append(f"| {case['path']}，{case['duration_s']:g}s / {case['period_s']:g}s | {100*local['cost_extra_fraction']:.3f}% | {100*cont['cost_extra_fraction']:.3f}% | {local['anchor_solves']} | {'通过' if local['screen_pass'] else '未通过'} |")
    smooth_rows=[]
    for case in smooth['scenarios']:
        local=next(row for row in case['scores'] if row['method']=='local_NEOC')
        nominal=next(row for row in case['scores'] if row['method']=='fixed_gain')
        robust=next(row for row in case['scores'] if row['method']=='fixed_midrange')
        smooth_rows.append(f"| {case['path']}，{case['duration_s']:g}s / {case['period_s']:g}s | {100*local['cost_extra_fraction']:.4f}% | {100*nominal['cost_extra_fraction']:.2f}% | {100*robust['cost_extra_fraction']:.2f}% |")
    report=f'''# 持续缓变接触闭环：第二步报告

2026-10-13。结论：**持续缓变环境中直接 NEOC 的机制得到理论小模型与接触仿真的支持，不能沿用此前“大跳变失败所以 idea 不成立”的判断。** 在更小预测邻域、更多查询的一组非线性接触测试中，NEOC 同时保持精度并减少总计算量；不过简单固定控制也足够好，目前尚未证明整个 R15 的任务必要性或顶刊增量。

## 1. 本轮对准了什么

环境在物理仿真中持续变化，各方法以相同周期获得准确的**当前**参数；未来参数不提供给控制器。普通 NEOC 更新不查询真实解、不计算候选验收残差、不做候选 rollout 检查。只根据预设参数邻域维护锚点，或在上一预测点维护灵敏度。后者确有计算成本，没有视为免费。

先使用有限时域、单价值 critic 的非线性柱面 Hertz 法向接触投影。状态为压入量、速度、实际刚度；动作是刚度变化率，实际积分刚度并记录能量。阻尼固定 8、质量 0.2、目标力约 0.1414 N；材料 E(t) 从 30 缓慢增至 36 kPa，另测平滑往复。控制器预测同一个 1.5s 冻结当前参数问题，物理环境在采样间持续变化。没有完整 MHE、空间预瞄、切向/图像目标或压缩 NN critic。

对照局部复用 NEOC、沿途维护 NEOC、cached-J chord、一次 Newton、完整结构化求解和固定刚度。完整求解采用一次向后 Bellman 扫描；此结构不能靠热启动减少扫描次数，不能用低效通用 Newton 冒充最强完整求解。全部使用同一个问题、状态网格、终端代价和动作约束。

## 2. 定位并修复数值载体

初始 Euler 预测的网格加密不能消除控制代价偏差；已完整保存在 [初始记录](00_初始Euler载体记录/pilot_results.json)。将控制模型的接触运动预测及阶段代价改为 RK4 后，三种网格相对独立非线性轨迹优化的代价差为 0.333%、0.294%、0.168%，见 [网格核对](resolution_audit.json)。独立参照为可行局部解，不是全局最优证明；其解析梯度方向核对误差约 2.4e−8。

另外复用旧极小动作计算 F_η，避免灵敏度维护时反复求新极小动作。两锚点上的导数与原差分实现相对差约 4.4e−10 / 6.0e−10，与完整扰动根差分相对差约 0.0122% / 0.000705%；准备耗时由约 3.3–3.8s 降至约 0.5s。见 [导数核对](01_复用动作计算灵敏度/derivative_audit.json)。这属于实现优化，未给 NEOC 添加新控制模块。

该节点插值载体仍是分段表示，不能宣称满足原正文整个邻域的 C² 根条件；以下参数半径是预设诊断配置，尚不是原 R15 的严格部署证书。

## 3. 四组预设测试与一组探索测试

原配置半径 0.025、更新周期 0.25/0.5s；所有方法均通过事后采样接触、增益/速率和能量检查。代价差以相同信息完整求解为参照。联合筛查仍用原定精度门槛（代价增加≤1%、力 RMSE 增加≤5%）和总计算减少≥30%。

| 路径，总时长 / 参数周期 | 局部复用代价差 | 沿途维护代价差 | 局部复用锚点数（含初始） | 局部复用联合门槛 |
|---|---:|---:|---:|---|
{chr(10).join(rows)}

沿途维护的控制表现接近参照，但每次维护仍需全局计算，其成本没有优于本载体的一次结构化完整求解。局部复用在单调路径上约 1.6% 的代价差，不能仅以“参数每步小”就默认其范围足够；这里还涉及有限邻域近似和策略对价值梯度的敏感性。

据此，**明确另开探索性测试**：保持材料路径、物理模型、代价和原门槛，把参数半径减半至 0.0125，更新周期改为 0.1s；不是把原失败结果删除，也不是改宽判据。见 [配置与结果](03_小邻域与密集更新_探索性/results.json)。

| 方法 | 累计代价 | 力 RMSE（N） | 总控制计算量（s） |
|---|---:|---:|---:|
{chr(10).join(f"| {r['method']} | {r['cost']:.8f} | {r['force_rmse_N']:.8f} | {r['total_controller_seconds']:.4f} |" for r in dense['scores'])}

此组 NEOC 代价增加 **{100*dn['cost_extra_fraction']:.3f}%**、力 RMSE 增加 **{100*dn['force_rmse_extra_fraction']:.3f}%**；60 个查询点仅初始加 14 次锚点重建，45 次直接预测。包含初始求解、灵敏度、维护、预测和低层动作恢复的总计算量，从完整求解 **{dr['total_controller_seconds']:.2f}s** 降到 **{dn['total_controller_seconds']:.2f}s**，减少 **{100*dn['total_time_reduction']:.1f}%**。该探索组复测一次 Newton/完整求解/固定刚度，未重新运行 chord；不能冒充对所有可能优化器的全面速度证明。

这是虚拟时间 6s 的质量仿真，求解延迟没有注入物理轨迹。NEOC 14 次维护超过 0.1s 查询周期，完整求解 59 次超期；总计算量仍大于虚拟任务时长。因此这里验证了相对计算摊销，**没有验证实时机器人部署**。上述计时都是本机 Python 原型测量，不是压缩 NN critic 的耗时。

## 4. 用平滑根单独验证“先证明范围，再直接 NEOC”

新增独立局部线性接触对照：用三系数二次价值 critic 的精确 Riccati 根恢复法向 PD 阻抗，实际增益按周期因果插值，材料有效刚度从 300 连续变化到 450 N/m；各控制器还受到相同的已知当前表面位移/速度。其物理模型和代价与上面的增益率优化不同，不能把两套结果拼成完整 R15 认证。

[解析条件](平滑对照的解析条件.md) 先给出光滑根、隐式灵敏度及增益曲率界；预设刚度误差≤0.1 N/m、阻尼误差≤0.002 Ns/m，得到参数半径 **{smooth['budget']['radius_N_per_m']:.3f} N/m**。还给出包含材料缓变、采样和增益插值偏差的齐次稳定条件；四组衰减裕量均为正。范围内直接预测，不依据候选结果决定是否发布。浮点实现做独立核对，未提供区间算术数值认证。

| 路径，总时长 / 参数周期 | 局部 NEOC 代价差 | 初始参数固定阻抗 | 区间中点固定阻抗 |
|---|---:|---:|---:|
{chr(10).join(smooth_rows)}

参照为同信息、冻结参数的完整解析控制解。局部 NEOC 四组差异均小于 0.03%；最大刚度误差约 0.0902 N/m，在事先界内，采样约束均通过。沿途维护也保持接近参照，但往复路径的累积误差较大；其结果见 [完整数据](02_平滑二次critic对照/results.json)。

中点固定阻抗比这个冻结模型参照还好，说明后者并非整条时变、有外部表面运动任务的全局最优解。不能把负代价差称作“NEOC 超越最优”，也不能只用初始固定基线的劣势宣称自适应必要。本小模型有极便宜的闭式解，不能靠与通用 CARE 求解器比较来证明 NEOC 的普遍加速。

## 5. 当前判断及后续

**机制支持继续推进：** 在先限定范围、正确维护灵敏度的条件下，缓变环境可直接 NEOC；少量维护配合较密查询能够获得计算摊销收益。不是“NEOC 后每次 Newton”，也不是“必须在线重新验收预测”。

**尚未建立顶刊级必要性：** 非线性简化任务中固定刚度与优化参照相当；平滑对照中中点固定方案甚至更好。下一步先建立与真实扫描任务对齐、能让同信息完整优化明确改善强固定/反应式基线的有限时域问题，并落实质量合格的平滑压缩 critic；再检验 NEOC 是否低成本保留那份收益。不要以此推翻 NEOC 原理，也不要现在扩展硬件或堆叠算法。预瞄控制与噪声估计尚未执行。

复现：16 条非线性基线轨迹跨原/优化实现字节一致；32 条平滑对照轨迹两次运行字节一致；密集组 NEOC 独立复跑轨迹字节一致。原/优化 NEOC 的四组差异均低于 1e−6 相对代价。计时未要求一致。见 [复现记录](reproducibility.json)；不声称所有非线性 NEOC 轨迹都重复验证。

主要入口：[共同载体](continuous_contact.py)、[四组运行](run_contact_validation.py)、[结构优化](run_structured_check.py)、[密集探索](run_dense_exploration.py)、[平滑对照](smooth_contact_check.py)、[复现脚本](verify_stage2.py)。原结果、失败和三份核心理论原件均保留。
'''
    (HERE/'README.md').write_text(report,encoding='utf-8')
    # Recover the intermediate smooth source corresponding to saved run1 if
    # only the later half-step stability budget edit separates the versions.
    first=read('02_平滑二次critic对照/results_run1.json')
    source=(HERE/'smooth_contact_check.py').read_text(encoding='utf-8')
    older=source.replace('    # Midpoint gain samples are held over DT; include their half-step offset.\n    lag=period+DT/2\n    delta_k=.1+lk*rate*lag;delta_d=.002+ld*rate*lag',
                         '    delta_k=.1+lk*rate*period;delta_d=.002+ld*rate*period')
    if hashlib.sha256((HERE/'smooth_contact_check.py').read_bytes()).hexdigest()==first['script_sha256']:
        (HERE/'02_平滑二次critic对照/results_run1_source.py').write_bytes((HERE/'smooth_contact_check.py').read_bytes())
    elif hashlib.sha256(older.encode('utf-8')).hexdigest()==first['script_sha256']:
        (HERE/'02_平滑二次critic对照/results_run1_source.py').write_text(older,encoding='utf-8',newline='\n')
    manifest={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='artifact_manifest.json'}
    parent=HERE.parents[1]
    inherited=[parent/'快速验证修复与闭环_2026-10-12/bellman_carrier.py',
               parent/'快速验证修复与闭环_2026-10-12/independent_reference.py',
               parent/'最小闭环验证_2026-10-12/run_experiment.py']
    record=dict(files=manifest,inherited_source_sha256={str(p.relative_to(parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inherited},
        history_note='Archived source copies retain original bytes; relative imports assume original folder layout when restoring a snapshot.')
    (HERE/'artifact_manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Wrote stage 2 report and artifact hashes.',flush=True)


if __name__=='__main__':main()
