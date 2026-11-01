"""Render scientific interpretation from frozen summary results, without raw trajectories."""
import json
from step3b_core import *

def main():
    op=json.loads((RES/'step3b_predictors.json').read_text(encoding='utf8'))
    cl=json.loads((RES/'step3b_closed_loop.json').read_text(encoding='utf8'))
    for name,v in cl['linear_tolerance'].items():
        r=[r for r in v['records'] if r['kappa']<=1]
        v['conservative_feedback_subset']=dict(kappa_grid=[.5,1.],common_interval=[max(x['interval'][0] for x in r),min(x['interval'][1] for x in r)],all_accepted_stable=all(x['max_radius_accepted']<1-1e-9 for x in r),max_radius_accepted=max(x['max_radius_accepted'] for x in r))
    old=json.loads((HERE.parent/'motion_aware_scan/results_v2/ff_loop.json').read_text(encoding='utf8'))
    cl['baseline_existing_coarse_grid_interval']=[max(v['interval'][0] for v in old.values()),min(v['interval'][1] for v in old.values())]
    cl['recommendation']='Do not replace FINAL IMM with either P1 or P3 on this evidence. P1 has no consistent paired gain, worsens a fast-push tail and narrows frozen-mode gain tolerance. P3 reduces pooled absolute surface-target error but worsens closed-loop force tracking and contact retention. P2 was not selected and its61-evaluation covariance fit did not converge; do not infer a GP-family upper bound.'
    (RES/'step3b_closed_loop.json').write_text(json.dumps(cl,indent=2,allow_nan=False),encoding='utf8')
    extra='\n## Feedback-sector dependence\n\n| Predictor | κ=0.5/1 common interval | Accepted band covered? | κ=0.5/1/2 common interval |\n|---|---|---|---|\n'
    for n,v in cl['linear_tolerance'].items():
        r=v['conservative_feedback_subset'];extra+=f"| {n} | {r['common_interval']} | {r['all_accepted_stable']} | {v['common_interval']} |\n"
    extra+='\nThe legacy approximate[0.75,1.2] tolerance is consistent with conservative feedback κ≤1, not the full grid includingκ=2. The existing ff_loop.json already has a coarse all-grid common interval[0.95,1.1]; our finer grid gives[0.925,1.125]. This is an assumption boundary, not a newly discovered baseline failure. P1 fails full b≤ln1.2 coverage even in the conservative subset; P3 has no analytic certificate and suffers contact losses in gain stress. The tested gain multipliers act on the existing estimated gain, so they are not the actual κ_g of every tick.\n'
    path=RES/'step3b_closed_loop.md';s=path.read_text(encoding='utf8');s=s.split('\n## Feedback-sector dependence')[0]+extra;path.write_text(s,encoding='utf8')
    lines=['# 第3步续：运动预测模块对比结论','', '**结论：本轮不建议用 P1 或 P3 替换 FINAL IMM。** 训练、选择和评价已完成；确认集未使用。', '',
           '## 数据与公平性','',
           '训练每族100场景（200000+，共500）；验证严格每族10（190000+，共50）；原生诊断每族20（110000+，共100）；分布外每族每剖面20（共300）。训练剖面含smoothstep、余弦、线性、最小急动及呼吸叠加；梯形速度、双段动作、呼吸频率突变仅用于分布外评价。闭环诊断每族50；推拉斜坡和脊比扫描各6档、每档20；增益压力四档、每档每族50。全部最终闭环结果共5450场景运行，逐拍轨迹不交付。', '',
           '所有候选复用同一 FINAL 量测；不向新预测器提供在线真值。P1和重调双CA敏感性对照P0r各61个固定候选，GP61次训练似然评估，小网络61个epoch、每epoch16384窗口（3170参数）；每类验证查看≤61次。相同的是可用场景、信息上限与选型预算，非等价的计算量或监督方式；GP使用500个训练窗口拟合，小网络使用随机重采样窗口和离线真值标签。计时均在批次结束后单独串行复测。', '',
           '曾发现并纠正物理标签符号及额外增广验证集两个实现问题。最终仅用规定的50个验证场景重放原先固定候选/训练过程，未依据初步评价增删架构、搜索范围或候选。重复计算的成本在协议中披露。本轮为探索性诊断，非预注册确认实验。', '',
           '## 开环：物理目标绝对误差','',
           '单位mm；ID=原生诊断，OOD=独立未见剖面。绝对误差包含累计量测/接触偏差，不能直接等同于短时预测能力。', '',
           '| 预测器 | ID RMSE / P95 | OOD RMSE / P95 | ID 1σ / 2σ覆盖 | 预测 P50/P99 ms |','|---|---:|---:|---:|---:|']
    for n in ('P0','P1','P2','P3','P0r'):
        a=op['evaluation']['id'][n]['groups']['all|all'];b=op['evaluation']['ood'][n]['groups']['all|all'];t=op['timing'][n]
        lines.append(f"| {n} | {a['rmse_um']/1000:.3f}/{a['p95_um']/1000:.3f} | {b['rmse_um']/1000:.3f}/{b['p95_um']/1000:.3f} | {a['coverage_1']:.3f}/{a['coverage_2']:.3f} | {t['p50_ms']:.3f}/{t['p99_ms']:.3f} |")
    lines+=['', 'P3总体绝对RMSE下降约11%，但ID硬组织推拉从0.545增至1.144mm、快速推拉从0.209增至1.071mm；这两族的呼吸、起始前50ms、斜坡、保持、回撤五阶段均变差。P1与P0的ID误差基本相同，回撤阶段仅有极小改善；OOD快速族绝对误差从1.091增至1.186mm。P2未入选闭环；其OOD快速族总体误差较低（0.754mm），但不是全阶段改善，训练似然优化也未收敛。逐族、逐阶段、逐OOD剖面以及消除常数后的局部时域/前馈增量诊断见step3b_predictors.json/md。', '',
           '## 闭环','', '单位N；表内为中位/均值，所有有限高峰值轨迹均保留。', '', '| 族 | Proposed | P1 | P3 | TrueK | Oracle+exc |','|---|---:|---:|---:|---:|---:|']
    for fam,g in cl['diagnostic'].items():
        names=('Proposed','P1','P3','TrueK','Oracle+exc')
        lines.append('| '+fam+' | '+' | '.join(f"{g[n]['median']:.3f}/{g[n]['mean']:.3f}" for n in names)+' |')
    lines+=['', '| 新版本－Proposed | S2 | T3m | 软推拉 | 硬推拉 | 快速推拉 |','|---|---|---|---|---|---|']
    for n in ('P1','P3'):
        cells=[]
        for f in FAMS:
            r=cl['diagnostic'][f][n]['paired'];cells.append(f"{r['mean_delta']:+.3f} [{r['ci95'][0]:+.3f},{r['ci95'][1]:+.3f}]")
        lines.append('| '+n+' | '+' | '.join(cells)+' |')
    lines+=['', '配对bootstrap以场景为单位、10000次、百分位95%区间，无多重比较校正，n=50/族；没有发散配对被排除。P1快速推拉中位数下降约10%，但均值上升、区间跨0；不能以中位数宣称改进。P3在S2、T3m、软推拉和快速推拉的均值均显著变差，硬推拉区间跨0。', '', '| 版本 | 失控>50N / 250 | 长失联>0.5s / 250 |','|---|---:|---:|']
    for n in ('Proposed','P1','P3','TrueK','Oracle+exc'):
        g=[v[n] for v in cl['diagnostic'].values()];lines.append(f"| {n} | {sum(x['runaway'] for x in g)} | {sum(x['long_loss_05'] for x in g)} |")
    lines+=['', '剂量扫描并未给出一致替换收益。P1在0.1s斜坡上的中位数较低（2.061对2.497N），但0.45s档出现1/20失控、均值8.087对0.286N；脊比45×有局部均值/尾部收益，30×及60×未形成一致优势。P3在高脊比的均值或失控计数有局部改善，同时中位数通常严重变差；较低失控不能掩盖失联或力偏差。所有剂量和配对区间保留于closed_loop.json/md。', '',
           '## κ_g容许性','', 'b≤ln1.2对应κ_g∈[0.8333,1.2]。冻结模式/稳态增益分析给出网格区间，不是切换系统证明。', '',
           '| 预测器 | 保守反馈 κ=0.5/1 | 加入 κ=2 的共同区间 | 覆盖接受区间？ |','|---|---|---|---|']
    for n,v in cl['linear_tolerance'].items():
        r=v['conservative_feedback_subset'];lines.append(f"| {n} | [{r['common_interval'][0]:.3f},{r['common_interval'][1]:.3f}] | [{v['common_interval'][0]:.3f},{v['common_interval'][1]:.3f}] | {r['all_accepted_stable']} |")
    lines+=['', '| 量测增益倍率 | Proposed 失控/250 | P1 失控/250 | P3 失控/250 | P3长失联/250 |','|---|---:|---:|---:|---:|']
    for gain,g in cl['gain_stress'].items():lines.append(f"| {gain} | {g['Proposed']['runaway']} | {g['P1']['runaway']} | {g['P3']['runaway']} | {g['P3']['long_loss_05']} |")
    lines+=['', 'P1/P2不能直接沿用原接受门限。P3没有解析容许区间，增益压力中低失控率往往伴随61–81/250长失联，不构成安全收益。原IMM近似[0.75,1.2]依赖保守反馈范围；κ=2时原结果也更窄（既有ff_loop.json粗网格共同区间[0.95,1.1]），这不是本轮才发现的新基线缺陷。倍率乘在原有估计增益上，不能当作每一拍的真实κ_g。', '',
           '## 代价、局限与下一步','', '| 控制器整体（含估计、预测、控制） | 串行P50/P99 ms |','|---|---:|']
    for n,r in cl['pooled_serial_timing'].items():lines.append(f"| {n} | {r['p50_ms']:.3f}/{r['p99_ms']:.3f} |")
    lines+=['', '新预测器仅替换前馈均值，保留旧IMM提供FINAL恢复速度，因此整体计时包含两者。P1/P0r用Numba、P0用NumPy；单独预测速度差是当前实现差异，不是复杂度定理。控制器计时包括实际失联/恢复阶段；失联时跳过预测也会降低用时，不能单凭P50判优。计时不是实时认证。', '',
           '物理目标存在适用范围：若F=k(t)(y+h)，则dF/k−dy=dh+(y+h)dlnk；黏弹、量测刚度误差及滤波进一步改变等效运动。因此跨刚度时η并非仅为物理表面加常数。原IMM中的累计“偏差”可能部分承担有用的等效扰动补偿；压低物理表面标签误差会损害闭环。这是两种目标的边界，不能据此否定物理运动补偿在硬组织推拉中的潜力。Oracle+exc也不是跨刚度族的普适性能下界。', '',
           '其余限制：仅一次小网络初始化和一个网络结构；GP只拟合500个训练窗口，61次似然评估未收敛；GP/网络未显式输入量测异方差，网络0.5s历史不能恢复任意长期偏置；所有模型沿用原量测时标，区间加权η与点值真值的时标差值得单独检查；预测数据来自旧策略，闭环状态分布变化是本轮已观察到的局限；参数型IMM方差未纳入刚度/黏弹偏差；OOD按独立场景而非与ID配对，不能用其总体RMSE较低宣称泛化更好；本轮诊断尾部区间较宽，且修正实现后重放了固定候选，不能视作预注册确认结论。', '',
           '建议保留FINAL IMM。后续优先检查量测有效时间与预测目标：分别预测物理运动、等效前馈增量，或以可辨识的短时增量构造监督损失，再用相同闭环协议验证；这是待检验方向，不预设收益。没有改写现有研究问题或主方法。', '',
           '## 复现与交付','', '依赖Python、NumPy、SciPy、Numba，无PyTorch依赖。所有新代码与汇总均在experiments/method_redesign/。主目录无修改、未提交或推送；初始pull因.git写权限被拒。补丁基于6a60559df7ec4d200e46c4807812a7e7377b2eae。缓存放在交付目录之外；Windows沙箱不允许multiprocessing命名管道，故使用可复现的文件传输工作进程。', '',
           '```powershell', 'python -X utf8 experiments/method_redesign/step3b_predictors.py --cache <data-cache> --workers 8', 'python -X utf8 experiments/method_redesign/step3b_analyse.py --cache <data-cache>', 'python -X utf8 experiments/method_redesign/step3b_closed_loop.py --cache <closed-cache> --workers 8', 'python -X utf8 experiments/method_redesign/step3b_finalize.py --cache <data-cache>', 'python -X utf8 experiments/method_redesign/step3b_validate.py', 'python -X utf8 experiments/method_redesign/step3b_report.py', '```', '', '新增校验覆盖被动记录不改闭环轨迹、原IMM重放、替换接口、因果前缀、Oracle标签、八种运动剖面导数与ff_loop谱半径；详见step3b_validation.json。']
    (RES/'step3b_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf8');print('summary written')

if __name__=='__main__':main()
