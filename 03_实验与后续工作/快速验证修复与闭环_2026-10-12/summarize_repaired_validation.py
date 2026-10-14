"""Check reproducibility and write the final bounded scientific decision."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
a=json.loads((HERE/'results_run1.json').read_text(encoding='utf-8'))
b=json.loads((HERE/'results.json').read_text(encoding='utf-8'))

def numerical_only(x):
    if isinstance(x,dict):
        return {k:numerical_only(v) for k,v in x.items()
                if not k.endswith('_seconds') and k not in ['latency_reduction_vs_fused_exact']}
    if isinstance(x,list):return [numerical_only(v) for v in x]
    return x

assert numerical_only(a)==numerical_only(b)
old_hashes=json.loads((HERE/'csv_hashes_run1.json').read_text(encoding='utf-8'))
for name,h in old_hashes.items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,name
repro=dict(numerical_metrics_and_candidate_hashes_exact=True,
           raw_CSV_files_exact=True,raw_CSV_count=len(old_hashes),
           wallclock_excluded=True,source_hashes_equal=a['source_hashes']==b['source_hashes'])
(HERE/'reproducibility.json').write_text(json.dumps(repro,indent=2),encoding='utf-8')
sr=b['stages']['shared_root'];up=b['stages']['parameter_update']
first=a['stages']['parameter_update'];local=b['local_change_check'];est=b['common_estimate']
reference=b['independent_reference'];fixed=up['fixed_gain_metrics']
table='\n'.join(f"| {r['method']} | {r['metrics']['cost']:.6f} | {100*r['relative_cost_extra_to_full']:.2f}% | {r['metrics']['force_rmse_N']:.6f} | {'通过' if r['metrics']['sampled_constraints_pass'] else '未通过'} |" for r in up['results'])
timing=[]
for r1,r2 in zip(first['results'],up['results']):
    timing.append(f"| {r2['method']} | {r1['generation_seconds']:.3f} / {r2['generation_seconds']:.3f} | {r1['with_standalone_checks_seconds']:.3f} / {r2['with_standalone_checks_seconds']:.3f} |")
resolutions='\n'.join(f"| {r['nd']}×{r['nv']} / {r['steps']}步 | {'集中于实际工作域' if r['focused'] else '均匀'} | {r['cost']:.6f} | {100*r['cost_gap_to_shooting']:.3f}% |" for r in b['resolution_checks'])
report=f'''# R15 快速验证修复与实际闭环结果

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: carrier repair + numerical validation + reproducibility
- Origin Date: 2026-10-13（工作开始于 2026-10-12）
- Verification Status: REPRODUCED，两次结果、候选哈希及 {len(old_hashes)} 个原始 CSV 完全一致；耗时单独保留
- Status: {b['status']}
- 本轮不改正文、附录、八项证明或既有失败记录。
- 主运行约 {a['duration_seconds']:.1f} / {b['duration_seconds']:.1f} s；两次进程正常结束，科研判据未通过。

## 1. 结论与真正停点

用户要求自主修复快速验证中的实现、设置和载体问题；只有证据反推出核心设计或必要性不足，才讨论调整 idea。按此要求，本轮**没有停在初始根失败**，完成了有限时域载体修复、独立策略质量核对和 NEOC 的实际接触闭环。

现在暴露的是三个科学问题：

1. 模量增加 50% 的主要任务中，一阶 NEOC 直接更新保持约束，但代价比可靠完整求解高约 29%；一次 Newton 的差距只有约 0.8%。
2. 增加 5% 时直接更新很准确，但固定刚度已达到完整求解的相同表现。
3. 主要任务中预先固定的刚度 150 仅比独立优化参照高约 0.15%。精确优化本身尚未带来有意义的必要收益；全局误差检查也吃掉大部分候选生成的加速。

因此当前最小任务没有提供足够理由投入完整顶刊级开发。**这是该任务中“直接更新的适用范围及必须在线优化阻抗的必要性”未成立的结果，不是程序无法运行，也不是全体 NEOC/R15 已被证明无效。** 尚未比较空间预瞄控制收益，不能将早期感知成功当作整条方法的成功。

## 2. 修复了哪些根因

| 旧问题 | 本轮处理 | 限制 |
|---|---|---|
| 原 NNBO 无限时域目标在持续接触中的含义未闭合 | 回到 R15 已有的 1.5 s 有限时域任务和终端代价 | 不声称修复或否定原论文全部理论 |
| 三项二次 critic 精确插值节点，却不能表达域内价值 | 使用固定网格的单价值函数表示和半拉格朗日 Bellman 方程 | 是数值验证载体，不是已冻结正文的精确配置点实现 |
| 未获得合格共同根 | 因时间方向三角结构，用逆时间 Bellman 求出唯一离散根，再求根灵敏度 | 离散根唯一不等于连续 HJB 精确最优 |
| 初期均匀网格策略受数值扩散影响 | 对慢法向速度/实际压入区集中分辨率，并细化时间 | 保留粗网格负结果；没有按方法分别调权重 |
| 只用同一个求解器评价自身策略 | 加入独立、带解析梯度的约束直接轨迹优化 | 该优化给出可行局部解，不冒充全局最优证明 |

最终离散共同根残差为 {sr['root_residual']:.1f}，方向灵敏度与独立微扰根的相对差约 {100*sr['directional_derivative_relative_error']:.3f}%。NEOC 全程零在线 Newton；终端代价的直接参数影响保持精确，不把末端条件也线性预测掉。

## 3. 物理场景与简化边界

仍是原定 30→45 kPa 的柱面软→硬扫描，R=30 mm、L=20 mm、接触初始半宽 3 mm、速度 8 mm/s、路径 -6→6 mm、总时长 1.5 s。目标法向力约 0.141372 N（见冻结配置/源码），是合成模型尺度，不是临床建议。

保留法向惯性、实际刚度状态及其变化率、原柱面空间接触核、增益/速率约束和能量账本。阻尼固定为 8 Ns/m，本轮不优化切向或姿态，也未完整验证原四维阻抗动作。

根状态为 (压入量、法向速度、实际刚度)。只在能量约束不激活的当前包络内消去能量的价值维度：即便忽略所有阻尼回收，能量下界仍为 {b['energy_inactive_envelope_bound']['lower_bound']:.6f} J，高于 0.001 J；实际轨迹另外核对能量。没有重置能量或将增益瞬时跳变。

价值函数 V=Φ(η)+W 的固定节点表示，终端 W=0；参数变化时 Φ 的直接变化与 critic 灵敏度同时处理。半拉格朗日下一状态价值为分段线性插值；每段二次动作代价的极小值枚举后取全局最小，连续动作是刚度变化率，不是手写增益轨迹。

网格边缘会截断**价值查询**，不截断物理运动。实际 RK4 轨迹始终在声明域内；下述细化与独立参照用于核对截断/离散影响，但不构成所有状态下的连续证明。该网格方程分段光滑，本次只测软→硬方向灵敏度，不冒充正文的全域 C²/PublishOK 认证。

## 4. 共同根策略质量

独立直接优化用 30 个真实刚度节点、显式增益/变化率约束和真实物理 RK4 积分；解析梯度的独立方向差分相对误差为 {reference['gradient_direction_relative_error']:.3g}，求解成功，最大变化率约 {reference['maximum_rate']:.3f}，低于 500 的限制。

| 网格/时间 | 布局 | 物理闭环代价 | 相对直接优化参照 |
|---|---|---:|---:|
{resolutions}

直接优化代价为 {reference['cost']:.9f}，最终 Bellman 策略约 0.104691，二者只差约 0.0021%。微小负差并不证明 Bellman 更优，参照的动作时间表示与运行代价积分规则稍有不同。采用它们的一致性作为快速质量核对，而非连续问题全局最优证书。

## 5. 同一因果接触历史与前方估计

两种方法共用旧策略运行的真实模型轨迹。法向力相对旧均匀模型超过 0.5%Fd 后，收集 12 个间隔 5 ms 的原始法向力/俯仰力矩样本。完整力矩含已知切向牵引力臂，拟合器没有读取生成数据的真值。

在 t={est['event_time']:.3f} s 得到前方模量 {est['estimated_front_E_Pa']:.1f} Pa、边界约 0 m；中心到边界约还有 {est['nominal_lead_time_s']:.3f} s。这里无噪声且几何、表面与切向参数已知，是有限历史拟合，不是完整 MHE 或真实组织验证。

## 6. 真正执行的参数更新与闭环对照

共同锚点对应均匀材料；新描述符来自相同拟合输出。所有方法更新同一个离散根，采用相同代价、动作恢复、物理积分、增益与能量检查。

| 方法 | 闭环代价 | 相对完整求解增加 | 力 RMSE N | 采样约束 |
|---|---:|---:|---:|---|
{table}
| 固定刚度 150 | {fixed['cost']:.6f} | {100*up['fixed_cost_extra_to_full']:.2f}% | {fixed['force_rmse_N']:.6f} | 通过 |

NEOC 保持了接触与约束，但代价增加 {100*up['results'][0]['relative_cost_extra_to_full']:.2f}%，未满足原任务 5% 的策略损失预算。这里并非估计误差造成：拟合误差接近浮点精度。方向导数核对也支持局部切线实现；50% 变化超出了一阶近似在本任务里的有效性能范围。没有用 Newton 修复 NEOC 后仍记作成功。

所有闭环策略质量对照采用**理想零更新延迟**，为了先隔离策略误差。下面真实 Python 计算时延另列，没有悄悄插入人为延迟让基线变差，也没有把离线轨迹宣称为实时部署结果。

## 7. 成本比较不能只看预测乘法

| 方法 | 两次候选生成 s | 两次模型准备+生成+独立全局检查+采样 rollout s |
|---|---:|---:|
{chr(10).join(timing)}

NEOC 候选生成确实很快，但全局残差检查占主导。完整逆时间求解可以在构造每层根时复用已计算的极小值，不应强迫它另做同样的完整残差扫描来弱化基线。允许这一结构融合后，完整求解为 {first['full_single_event_total_seconds']:.3f}/{up['full_single_event_total_seconds']:.3f} s，NEOC 端到端只减少约 {100*first['latency_reduction_vs_fused_exact']:.1f}%/{100*up['latency_reduction_vs_fused_exact']:.1f}%，达不到本轮 30% 的投入判据。

NEOC 额外灵敏度准备耗时约 {sr['sensitivity_seconds']:.3f} s，也单独计入。本轮只有一个新描述符事件，没有证据把准备成本摊到虚构的无限次查询；若后续主张多次复用，基线也应允许缓存及同等预计算。

这个原型有约 778 万个 critic 网格系数，计算量不代表正文的特定低维配置实例。所有带全局检查的更新耗时都超过实际扫描剩余时长，故**未验证实时性**。若改成离线误差包络或局部证书，可改变成本，但尚未实现，不能把它的假想收益计入当前结果。

非扩张的离散 Bellman 算子给出根误差界 E≤Σ_k||F_k||∞、离散 MDP 贪心策略差界≤2E；NEOC 的该全局界约 {up['results'][0]['discrete_MDP_policy_gap_bound']:.3f}，远不能认证本轮小代价预算。这个界针对投影网格模型，不是物理 ODE 的近最优证书，不能把“采样约束通过”写成 PublishOK 通过。

## 8. 小变化是否有效

模量只增加 5% 时，NEOC 代价为 {local['NEOC']['cost']:.9f}，完整求解为 {local['full']['cost']:.9f}，额外损失约 {100*local['relative_NEOC_cost_extra']:.4f}%，约束通过。这支持局部直接更新机制，而非所有幅度变化。

但此时固定刚度的代价也为 {local['fixed']['cost']:.9f}，与完整求解相同。主要 50% 任务中固定刚度也只比独立参照高 {100*up['fixed_cost_extra_to_independent_reference']:.3f}%。因此，本任务还不能证明“必须实时求最优阻抗”具有足够实际收益。

## 9. 接下来的科学决策

载体和实现问题已继续修复；现在应讨论的是核心问题—机制匹配，而非再次因求根失败中止：

- 若主张大幅跨区变化后一次 NEOC、零校正即可近最优，本轮结果不支持。
- 若主张小变化内维护策略，应同时找到确实需要策略调整、且强简单方法不能等价解决的代表性真实 failure；本轮小变化例子没有做到。
- 若主张可认证端到端加速，应实例化比全局扫描便宜的有效证书，并计入真实准备与复用；不能只报告候选生成速度。

不为了获得正结果而放宽预算、提高基线难度或随意更换权重。未来更复杂法向—切向耦合、阻尼优化或动态图像任务可能不同，本轮没有验证它们，所以不否定整个 R15；但也不能把“也许复杂任务有用”当作已通过顶刊 idea 快速验证。

按约定未进入预瞄控制收益比较：已有早期感知结果，参数更新的必要收益尚未成立。本轮的停止理由现在是设计主张/必要性需讨论，而不是不愿继续修复验证载体。

## 10. 文件

- [主运行脚本](run_repaired_validation.py)、[Bellman 载体](bellman_carrier.py)、[独立优化参照](independent_reference.py)。
- [配置](config.json)、[完整结果](results.json)、[首次结果](results_run1.json)、[复现检查](reproducibility.json)。
- [共同传感历史](common_wrench_history.csv)、[独立优化轨迹](independent_shooting_trace.csv)、[刚度节点](independent_shooting_knots.csv)。
- `trace_NEOC.csv`、`trace_cached_chord.csv`、`trace_one_Newton.csv`、`trace_full_backward_solve.csv`、`trace_fixed_gain.csv` 为完整物理轨迹。
- [首次日志](run_log.txt)、[复现日志](run_log_reproduction.txt)、[整理脚本](summarize_repaired_validation.py)。

本目录执行 `python -X utf8 run_repaired_validation.py`。依赖 NumPy/SciPy，实验规模和成本已在配置中声明。旧失败目录和当前三份理论文件保持原样。
'''
(HERE/'README.md').write_text(report,encoding='utf-8')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir()
          if p.is_file() and p.name!='artifact_manifest.json'}
(HERE/'artifact_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('Repaired physical results and raw CSVs reproduced exactly; report written.')
