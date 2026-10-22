"""Build the final Chinese experiment section from the completed raw results."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import complete_task_evidence as E

ROOT = E.ROOT.parents[1]
DOC = ROOT / 'research/paper_experiments_draft.md'
LABELS = {
    'Best_fixed': '调优固定控制 B1', 'EtR-A': '估计＋规则 EtR-A',
    'EtR-A-Haninger': 'Haninger 阻尼规则适配', 'EtR-R': 'RLS＋规则',
    'Duan': 'Duan 适配', 'Song': 'Song 适配', 'Jin_OS': 'Jin 单步', 'Jin_MS': 'Jin 多步',
    'Xue': 'Xue：初始表面先验', 'Xue_nominal': 'Xue：名义表面',
    'Learned_VIC': '学习式 VIC（ES）', 'Candidate1_v2': 'C1 v2：模型优化',
    'Candidate2': 'C2：学习'}


def rows(name):
    return list(csv.DictReader((E.OUT / name).open(encoding='utf-8', newline='')))


def table(lines, heads, values):
    lines += ['| ' + ' | '.join(heads) + ' |', '| ' + ' | '.join(['---'] * len(heads)) + ' |']
    lines += ['| ' + ' | '.join(map(str, r)) + ' |' for r in values]
    lines += ['']


def main():
    E.verify()
    summary = {(r['task'], r['method']): r for r in rows('summary.csv')}
    raw = rows('per_scene.csv')
    ablations = rows('ablations.csv')
    assert len(raw) == 2040 and len(ablations) == 1040
    assert all(int(r['failed']) == 0 for r in raw + ablations)
    lines = ['# V. 仿真实验：任务设置与最终结果', '',
        '围绕“利用新接触信息及时更新阻抗、改善接触控制”的问题，设置刚度变化、独立表面运动和组合扫描三类任务。'
        '任务结构参照 Xue 的撤离／推压、Song 与 Li 的刚度／表面变化实验；各任务的作用见下表，原文映射保留于附录 A。'
        '本节使用新完成的独立场景，不混用旧版本结果。', '',
        '## A. 实验设置与任务目的', '',
        '采用位置控制机器人模型 PC125：外环周期 8 ms，传输延迟 8 ms，位置伺服时间常数 19 ms；力测量 1 kHz、低通 30 Hz、噪声标准差 0.08 N。'
        '目标接触力 4 N，积分步长 0.5 ms。接触包含有限足迹、变模量、黏弹性和摩擦；这些是合成工程模型参数。', '']
    table(lines, ['任务', '设置', '目的'], [
        ['T1：独立刚度变化', '扫描 8 s；软组织 10–40 kPa、软硬模量比 30；关闭呼吸与推拉', '隔离材料变化，检验环境变化后的控制调节'],
        ['T2：独立表面推拉', '均匀软组织 10–40 kPa／硬组织 0.5–1.5 MPa；关闭呼吸；2 s 开始 ±5 mm 推拉，0.5 s 过渡、保持 2 s 后返回；共 8 s', '隔离表面运动，检验接触保持、冲击及恢复；软／硬×推／拉四格各 10 场景'],
        ['T3：组合扫描', '扫描 10 s；骨样软硬区域模量比 10–30，呼吸幅值 1 mm，随机推拉 3–8 mm、过渡 0.5 s', '检验空间异质性与组织运动同时存在时的综合表现']])
    lines += ['测试种子 7200，每任务 40 个相同场景及噪声实现，协议先冻结、随后运行。'
        '控制器沿用 train 1000／dev 2000 选定的参数和权重：T1 使用原 T1 配置，T2／T3 使用原 T3 配置。'
        'T1 去呼吸及 T2 为冻结参数的结构隔离／迁移检查，不将它们称为每种新任务重新调参后的算法能力上限。'
        'T3 保持原训练任务分布，提供新的主比较。神经方法各 3 个训练种子，先在场景内平均，独立场景数仍为 40。', '',
        '对照沿用已修正适配和扩展调参后的版本，包括 Xue 刚度／投影修正、Jin 单步动作可见端点修正、学习式对照的参考输出和训练校准。'
        'Haninger 只适配其阻尼规则；学习式 VIC 是 E27 启发的 ES 训练对照。当前 C1 是估计＋直接优化＋参考补偿＋监督层，C2 是历史信息／锁相估计驱动的蒸馏网络；二者不冒称完整 B0 实现。', '',
        '指标使用去掉前 1 s 后、8 ms 日志中的真实接触力。报告力 RMSE、峰值、零力时长；峰值是采样峰值。'
        'T2 恢复时间从每段运动斜坡结束计算：首次进入 3.5–4.5 N 并保持 25 个连续采样。'
        '每段观察 2 s，未恢复者保留并按 2 s 右删失；“截断平均恢复时间”包含这些事件，不是只对成功事件求均值。', '',
        '## B. 全方法主结果', '',
        '**表 1：力 RMSE，单位 N，场景均值 ± 场景标准差。**', '']
    table(lines, ['方法', 'T1：刚度变化', 'T2：均匀组织推拉', 'T3：组合扫描'],
          [[label] + [f"{float(summary[t, name]['rmse_N_mean']):.4f} ± {float(summary[t, name]['rmse_N_std']):.4f}"
                      for t in E.TASKS] for name, label in LABELS.items()])
    lines += ['主比较 2,040 次运行和关键消融 1,040 次运行均未出现积分失败。'
        'C1／C2 在 T1 与 T3 的平均力 RMSE 均低于全部对照；在 T2 则不是最优，Haninger 规则适配的平均 RMSE 最低。'
        '有限均值排名不等于所有场景或全部物理指标占优。', '',
        '**表 2：T2 的动态与接触表现。峰值为场景峰值的平均，零力时长为单次完整评估的平均；恢复指标合并去程与返回，保留删失。**', '']
    selected = ['Best_fixed', 'EtR-A', 'EtR-A-Haninger', 'Jin_MS', 'Candidate1_v2', 'Candidate2']
    vals = []
    for name in selected:
        r = summary['T2_motion', name]
        recovery = (float(r['out_recovered_mean']) + float(r['return_recovered_mean'])) / 2
        restricted = (float(r['out_restricted_recovery_s_mean']) + float(r['return_restricted_recovery_s_mean'])) / 2
        vals.append([LABELS[name], f"{float(r['peak_force_N_mean']):.3f}", f"{float(r['no_contact_s_mean']):.4f}",
                     f'{100 * recovery:.2f}', f'{restricted:.3f}'])
    table(lines, ['方法', '平均峰值／N', '零力时长／s', '2 s 内恢复率／%', '截断平均恢复时间／s'], vals)
    lines += ['T2 的材料分层说明了失败条件：C1 在均匀软组织上的 RMSE 为 0.0572 N，在硬组织上增至 1.3717 N；'
        'Haninger 适配分别为 0.1778／0.4330 N，C2 为 0.1428／0.4926 N。'
        '因此 C1 的综合优势没有延伸到持续硬接触推拉，不能只展示软组织结果。'
        '这组结果评价冻结配置的迁移，不证明所有方法在独立重新调参后的排名。', '',
        '## C. 实际运动响应', '',
        '图 1 展示 T2 的首个软／硬组织推压场景（编号 0／2；神经模型 0），同时给出接触力与实际法向位移。'
        '选择由材料类别和场景顺序确定，不按误差挑选。四条代表曲线用于解释响应，全部 13 种方法仍在表 1 与逐场景数据中。'
        '无平滑、无轨迹截除、无超范围隐藏。', '',
        '![独立表面推拉下的接触力与实际法向运动](../experiments/cobot_scan_benchmark/completed_tasks/figures/motion_response.png)', '',
        '法向速度由 8 ms 实际位置差分获得。T3 中 C1／C2／B1 的平均速度 RMS 分别为 5.818／5.880／5.919 mm/s，'
        '数值接近，不能据此宣称显著的运动平顺性优势；T2 对应为 6.805／4.738／13.333 mm/s。'
        '速度 RMS 本身同时包含必要跟随与不必要运动，不作为“越小越好”的单独排序指标。'
        '当前横向通道由共同位置伺服独立控制，因此现有模型不能证明横向轨迹跟踪改善。', '',
        '## D. 最终方法的机制消融', '',
        '以下在相同新场景上对最终冻结实现做运行时干预，不重新训练或调参；它们检验该实现的依赖关系，不代表每种删减结构的充分优化性能。'
        '其中 C2 的力特征遮蔽会改变训练输入分布，故不把其大幅退化当作重新训练的无力信息方法性能。', '',
        '**表 3：关键消融的力 RMSE，单位 N。**', '']
    ab = {}
    for r in ablations:
        ab.setdefault((r['task'], r['method']), []).append(float(r['rmse_N']))
    configs = [('C1 完整', 'Candidate1_v2'), ('C1 去表面运动估计', 'C1_no_surface_estimate'),
               ('C1 只调增益', 'C1_gains_only'), ('C1 刚度信息延迟 100 ms', 'C1_information_delay100ms'),
               ('C1 每 40 ms 更新', 'C1_update40ms'), ('C2 完整', 'Candidate2'),
               ('C2 只调增益', 'C2_gains_only'), ('C2 每 96 ms 决策', 'C2_update96ms')]
    values = []
    for label, name in configs:
        values.append([label] + [f"{float(summary[t, name]['rmse_N_mean']) if name in LABELS else np.mean(ab[t, name]):.4f}"
                                for t in ('T1_static', 'T3_combined')])
    table(lines, ['配置', 'T1', 'T3'], values)
    lines += ['C1 在 T3 中只调增益时从 0.1471 N 增至 0.4888 N，说明参考补偿对完整方法的优势具有关键作用。'
        '因此不能把全部收益归因于快速更新刚度／阻尼。去表面运动估计在无表面运动的 T1 中也明显退化，'
        '意味着该状态还可能补偿接触／执行模型误差，不能把这一消融解释成已准确辨识真实表面运动。'
        '放慢 C1／C2 决策均使误差上升；额外 100 ms 刚度信息延迟的影响则在 T1 明显、T3 较小。', '',
        '## E. 在线更新时间与局部数值质量', '',
        '计时在同一机器、单个仿真进程、预热后进行，每任务每模型 3 个完整场景；从观测对象组装至估计、决策和输出限幅。'
        '不包含仿真动力学积分，也不伪称已测量实物传感器通信与机器人发布。神经方法合并三个模型的原始调用时间后计算分位数。', '',
        '**表 4：在线更新时间 P99，单位 ms。**', '']
    times = np.load(E.OUT / 'timing_raw.npz')
    values = []
    for name in ('EtR-A', 'EtR-A-Haninger', 'Jin_OS', 'Jin_MS', 'Xue', 'Learned_VIC', 'Candidate1_v2', 'Candidate2'):
        cells = []
        for t in ('T1_static', 'T3_combined'):
            keys = [f'{t}_{name}_{seed}' for seed in (-1, 0, 1, 2)]
            a = np.concatenate([times[k] for k in keys if k in times.files])
            cells.append(f'{np.percentile(a, 99):.4f}')
        values.append([LABELS[name], *cells])
    table(lines, ['方法', 'T1', 'T3'], values)
    tt = np.concatenate([times[k] for k in times.files])
    lines += [f'全部 {len(tt):,} 次已记录软件更新中，超过 8 ms 的次数为 {int((tt > 8).sum())}。'
        '固定控制和 Duan 为预计算增益，只有初始化接口计时，不将三个初始化样本与逐周期更新混作等样本统计。'
        '当前所有在线对照均明显快于控制周期，C2 快于 C1，但模型优化耗时尚未构成该任务的主要瓶颈。', '',
        '最终 C1 在两任务各取前 6 个场景，以固定间隔抽取 234 个实际优化快照。'
        '保持同一状态、估计、目标、约束和动作类，从当前解、前次解及规则先验三个初值，分别运行 25／50／100 次 GN，取较低代价作为数值参考。'
        '参考在 50→100 次迭代时的最大相对代价变化低于 3×10⁻⁸%，表明该比较中的数值代价已经稳定；这不是全局最优证明。', '']
    refs = rows('reference.csv')
    values = []
    for t in ('T1_static', 'T3_combined'):
        a = np.array([float(r['gap_percent']) for r in refs if r['task'] == t])
        values.append([t, len(a), f'{a.mean():.3f}', f'{np.median(a):.3f}', f'{np.percentile(a, 95):.3f}', f'{a.max():.3f}'])
    table(lines, ['任务', '快照数', '平均代价差／%', '中位数／%', 'P95／%', '最大／%'], values)
    lines += ['大多数快照接近这一局部参考，但个别差值显著，不能声称每次更新都保证近优。'
        '代价差包含该 OCP 已定义的正则项，不冒称物理误差百分比。C2 没有相同动作类与目标下的最优性保证，本次不以力 RMSE 排名替代该保证。', '',
        '## F. 统计与结论', '',
        '场景是统计单位。全部 66 个“候选－对照－任务”的 RMSE 差异使用场景配对 bootstrap（100,000 次），'
        '并以 Bonferroni 调整得到约 95% 家族覆盖的区间；它们是条件于已有训练模型的区间，不额外推断训练算法的总体分布。'
        '完整区间见 completed_tasks/paired.csv，不能把均值排名自动称为每项差异均有统计证据。', '',
        '调整后的区间中，C1 在 T1／T3 相对全部 11 个对照的 RMSE 差均低于零；C2 在 T1 也如此，'
        '但 T3 相对 EtR-A 和学习式 VIC 的区间跨零。因此 C2 的 T3 均值领先尚不足以支持其对全部对照的差异均成立。', '',
        '1. **控制效果：**C1／C2 在独立刚度变化和原分布组合扫描中降低平均力误差；均匀硬组织推拉存在明确反例，故不支持所有任务／指标全面领先。',
        '2. **收益来源：**当前 C1 的优势来自估计、参考补偿与及时更新的共同作用，证据不足以单独归因于阻抗更新算法。',
        '3. **计算与近优：**软件更新均满足本机测得的 8 ms 预算；现有平台下尚无主要优化计算瓶颈，且局部参考存在显著个别差值，不支持统一的近优保证。',
        '4. **证据范围：**本轮完成三类任务、全对照、关键机制、实际法向运动和软件时效比较。未增加跨平台／能耗／成像指标；这些不构成本轮接触控制仿真的必做清单。真实机器人应用与完整 B0 认证仍不由本仿真替代。', '',
        '## 附录 A：逐篇“研究问题—实验任务—评价指标”核读映射', '']
    # Keep the focused original-paper reading, rather than regenerating new claims.
    old = DOC.read_text(encoding='utf-8')
    marker = '## 附录 A：逐篇“研究问题—实验任务—评价指标”核读映射'
    appendix = old.split(marker, 1)[1].split('## 附录 B：', 1)[0] if marker in old else ''
    lines += [appendix.strip(), '', '## 附录 B：复算与原始证据', '',
        '- 新主结果、消融、恢复时间、计时原始样本及数值参考：[completed_tasks](../experiments/cobot_scan_benchmark/completed_tasks/README.md)。',
        '- 旧任务与适配修正证据：[refined_baselines](../experiments/cobot_scan_benchmark/refined_baselines/results.md)。',
        '- 生成正文：`python -X utf8 experiments/cobot_scan_benchmark/write_completed_report.py`。', '']
    DOC.write_text('\n'.join(lines), encoding='utf-8', newline='\n')
    sources = ['per_scene.csv', 'ablations.csv', 'timing.csv', 'timing_raw.npz', 'reference.csv', 'summary.csv', 'paired.csv', 'protocol.json', 'checks.json']
    E.save('report_manifest.json', dict(primary_runs=len(raw), ablation_runs=len(ablations), failed_runs=0,
           original_log_exact_checks=51, numerical_reference_snapshots=len(refs),
           timing_samples=len(tt), updates_over8ms=int((tt > 8).sum()),
           sources={s: hashlib.sha256((E.OUT / s).read_bytes()).hexdigest() for s in sources},
           report_sha256=hashlib.sha256(DOC.read_bytes()).hexdigest(),
           generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           analysis='scene-paired; neural seeds averaged within scenes; all primary/ablation runs retained'))
    print(DOC)


if __name__ == '__main__':
    main()
