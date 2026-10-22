"""Build the Chinese experiment section from primary results and paper evidence.

No controller, tuning record, original paper or frozen snapshot is modified.
Missing measurements remain explicit; historical cohorts/versions stay separate.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
B = ROOT / 'experiments/cobot_scan_benchmark'
OUT = B / 'paper_tables'
DRAFT = ROOT / 'research/paper_experiments_draft.md'

# Reading scope is research-question + experimental/evaluation passages, not
# a claim to have re-read all proofs or all 36 papers in the earlier corpus.
PAPERS = [
    ('Xue，RAS 2025', 'xue.txt', [(1, 45), (475, 630)],
     '未知软环境和突发推拉下兼顾接触保持、柔顺性与约束',
     '非均匀硅胶的撤离／推压、位置约束和手臂交互',
     '分扰动阶段的力 RMSE、最大绝对误差、推压后的力增量、每步计算时间；位置／力／阻抗曲线',
     '迁移分阶段接触与冲击指标、更新计时；不据力误差声称成像效果'),
    ('Jin，TMECH 2023', 'field_papers/E41.txt', [(1, 45), (629, 814)],
     '在线调整运动精度与柔顺性的权衡，同时处理约束',
     '轨迹遇障、状态约束、不同速度的自由空间轨迹／插入',
     '最大／RMS 位置误差与最大／RMS 接触力，速度和增益轨迹',
     '力与运动必须分开验证；只有存在独立运动目标和可受方法影响的运动通道，位置误差才有比较意义'),
    ('Song，RA-L 2025', 'field_papers/E04.txt', [(1, 38), (382, 477)],
     '环境位置和刚度未知时实现鲁棒力跟踪',
     '刚度突变平面、斜面、正弦面；仿真与实际接触',
     '力误差、超调和调节时间，位置／力响应曲线',
     '补环境变化后的恢复时间，不能只用整段平均误差掩盖短时迟滞'),
    ('Duan，RAS 2018', 'field_supplements/F02.txt', [(1, 44), (416, 527)],
     '未知环境和变化期望力下的自适应接触力跟踪',
     '平面刚度变化、斜面、正弦／复合面、动态期望力',
     '位置与力的时域响应、稳态偏差和适应过渡；不虚构原文不存在的统一统计总表',
     '采用能区分刚度变化与表面运动的任务；位置曲线不等于位置跟踪误差'),
    ('Li，T-RO 2025（E27）', 'field_papers/E27.txt', [(1, 38), (578, 718), (784, 874)],
     '通过仿真阻抗学习提高未知地形上的力跟踪与泛化',
     '平面／斜面／复杂面、训练迭代、未见地形及加工／剥皮／人机任务',
     '力误差、调节时间、最大超调与振荡；训练阶段误差；应用任务另报剥皮厚度或人用力',
     '训练曲线和未见环境对学习路线有意义；剥皮／人用力不适用于当前扫描模型'),
    ('Haninger，RAS 2023', 'haninger.txt', [(1, 44), (643, 748)],
     '用未知接触模型和任务代价调节交互阻抗',
     '未知刚度接触、导轨装配、人工引导抛光与振动扰动',
     '高／低频实际速度 RMS（表 3 以 15 Hz 分带）、计算规模；同时检查任务运动',
     '接触运动的振荡代价比套用无人机平顺性更直接；频带和激励应统一公开，不能单独按速度越小越好排序'),
    ('Kong，TSMC: Systems 2025（NNBO）', 'nnbo.txt', [(1, 39), (610, 779)],
     '改善软环境模型，并通过最优阻抗关系优化交互性能',
     'Franka 与硅胶／海绵交互、HC／KV 模型和不同目标权重',
     '力模型 ME／AE／MaxE，位置／力轨迹，收敛代价和累计代价',
     '最优控制论文可以报告已定义的任务代价；不能把本项目调参综合分数冒称物理量或最优性证明'),
    ('Matschek，TCST 2023', 'field_supplements/F04.txt', [(1, 43), (578, 655)],
     '模型不确定时同时控制运动和接触力并满足安全约束',
     '机器人接触路径、GP 混合模型、输入扰动及实物验证',
     '力 RMS 误差、位置输出与约束满足；总计算时间及传感器通信／OCP 分解',
     '完整链路计时与约束验证应独立于纯策略调用时间'),
    ('Jin，TMECH 2025', 'field_supplements/F06.txt', [(1, 44), (816, 868)],
     '稳定可变阻抗下保持精度、柔顺性和平顺性',
     '机器人抛光、能量罐对比、近似等跟踪精度的控制比较',
     '最大／平均控制力、圆轨迹误差、控制力变化平顺性指标',
     '评价需对应任务；圆轨迹误差不搬到无圆轨迹目标的扫描任务，也不混淆控制力与接触力'),
    ('YOPO，RA-L 2024', 'method_transfer/YOPO.txt', [(1, 45), (476, 545)],
     '把在线规划负担转移至离线学习，兼顾规划质量与低延迟',
     '同目标的优化专家比较、不同障碍密度／速度的森林飞行',
     '候选平均／最佳轨迹代价，分模块／总延迟，障碍距离，积分平方 jerk、路径长度和成功率',
     '迁移质量—时间联合评价和离线／在线分工；不强行搬入障碍距离和飞行成功率'),
    ('Wang，TASE 2024', 'field_supplements/F09.txt', [(1, 40), (475, 535)],
     '位置控制机器人上的六维柔顺运动和力调节',
     '独立自由空间轨迹、外力扰动及超声接触',
     '独立轨迹的平均／最大位置误差、相对力误差和各维响应',
     '运动评价来自可比较的运动任务，而不是给当前统一横向伺服轨迹补一列无区分度误差'),
    ('AuSoScan，TMECH 2025', 'field_supplements/F08.txt', [(1, 43), (515, 541), (568, 621)],
     '软接触建模、图像感知和脊柱超声自动扫描／重建',
     '阶跃／斜坡／波形力与路径跟踪、模体和人体脊柱扫描',
     '力／位置跟踪误差，图像定位置信和脊柱重建结果',
     '力控制收益与成像／重建质量分层；本项目无声学或图像模型，不能把力合格率当成成像质量'),
]

LABELS = {
    'Best_fixed': '调优固定控制 B1', 'EtR-A': '估计＋刚度规则 EtR-A',
    'EtR-A-Haninger': 'Haninger 阻尼规则适配', 'EtR-R': 'RLS＋规则',
    'Duan2018': 'Duan 适配', 'Song2025': 'Song 适配',
    'Jin_OS': 'Jin 单步', 'Jin_MS': 'Jin 多步',
    'Xue_A1': 'Xue：初始表面先验', 'Xue_nominal_ys0': 'Xue：名义表面',
    'Learned_VIC': '学习式 VIC（ES，3 种子）',
    'Candidate1_v2': 'C1 v2：模型优化路线', 'Candidate2': 'C2：学习路线（3 种子）',
}


def read(name):
    return json.loads((B / name).read_text(encoding='utf-8'))


def csvread(name):
    rows = list(csv.DictReader((B / name).open(encoding='utf-8', newline='')))
    for r in rows:
        for k, v in r.items():
            if v == '': r[k] = None
            elif k not in ('method', 'task'): r[k] = float(v)
    return rows


def write_csv(name, rows):
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)


def f(x, n=4):
    return '—' if x is None else f'{x:.{n}f}'


def table(lines, head, rows):
    lines += ['| ' + ' | '.join(head) + ' |', '|' + '|'.join(['---'] * len(head)) + '|']
    lines += ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows]
    lines += ['']


def bad_scenes(entry):
    return sum(v is None or not math.isfinite(v) for v in entry['per_scene']['rmse'])


def complete_mean(entry):
    return entry['summary']['rmse'] if not bad_scenes(entry) else None


def shown_mean(entry):
    n = bad_scenes(entry)
    return f'NA（{n}/40 场景含非有限值）' if n else f(entry['summary']['rmse'])


def main():
    OUT.mkdir(exist_ok=True)
    summary = csvread('refined_baselines/summary.csv')
    lookup = {(int(r['cohort']), r['task'], r['method']): r for r in summary}
    main_rows = [r for r in summary if r['method'] in LABELS]
    assert len(main_rows) == 52
    write_csv('main_outcomes.csv', main_rows)
    c1, c2, gaps = read('candidate1_results.json'), read('candidate2_results.json'), read('gaps_results.json')
    paired = list(csv.DictReader((B / 'refined_baselines/paired_differences.csv').open(encoding='utf-8')))
    primary_pairs = [r for r in paired if r['cohort'] == '4000' and r['metric'] == 'rmse_N'
                     and r['candidate'] in ('Candidate1_v2', 'Candidate2')
                     and r['reference'] == ('Jin_OS' if r['task'] == 'T1' else 'EtR-A')]
    assert len(primary_pairs) == 4
    write_csv('paired_primary_comparisons.csv', primary_pairs)

    solvers = []
    for task, methods in c1['solver_study'].items():
        for name, v in methods.items():
            solvers.append(dict(task=task, solver=name, reference_gap_percent=100*v['gap_rel_mean'],
                                mean_policy_ms=v['t_update_mean_us']/1000,
                                mean_estimator_ms=v['t_est_mean_us']/1000,
                                ocp_p99_ms=v['t_ocp_p99_us']/1000,
                                max_policy_ms=v['t_update_max_us']/1000))
    write_csv('solver_quality_cost.csv', solvers)
    sx = {(r['task'], r['solver']): r for r in solvers}

    abls = []
    c1_names = {'C1': 'C1 v1 完整', 'abl_a_no_surface_motion': 'C1 去表面运动估计',
                'abl_b_gains_only': 'C1 只调增益', 'abl_d1_klat_0.1s': 'C1 环境信息额外延迟 100 ms',
                'abl_d2_update_40ms': 'C1 每 40 ms 决策'}
    for task in ('T1', 'T3'):
        for name, label in c1_names.items():
            abls.append(dict(task=task, method=label, rmse_N=c1['test'][task]['methods'][name]['test']['rmse'], cohort=3000))
    c2_names = {'cand2_full': 'C2 完整', 'cand2_a': 'C2 无教师', 'cand2_b': 'C2 仅教师训练',
                'cand2_c': 'C2 只调增益', 'cand2_d': 'C2 移除力派生信息', 'cand2_e': 'C2 每 96 ms 决策'}
    for cell in c2['cells'].values():
        if cell['dose'].get('ramp', .5) != .5: continue
        for name, label in c2_names.items():
            abls.append(dict(task=cell['task'], method=label, rmse_N=cell['methods'][name]['mean']['rmse'], cohort=3000))
    write_csv('information_ablations.csv', abls)
    ax = {(r['task'], r['method']): r['rmse_N'] for r in abls}

    robustness = []
    for cell, data in gaps['c2_seed5000'].items():
        for name, val in data['methods'].items():
            robustness.append(dict(family='task', condition=cell, task=data['task'], cohort=5000,
                                   method=name, rmse_N=complete_mean(val), nonfinite_scene_count=bad_scenes(val),
                                   raw_summary_rmse_N=val['summary']['rmse']))
    incomplete = []
    for platform, cells in gaps['platform_seed6000'].items():
        for cell, data in cells.items():
            if 'methods' not in data:
                incomplete.append(platform + '/' + cell); continue
            for name, val in data['methods'].items():
                robustness.append(dict(family='platform', condition=platform, task=cell, cohort=6000,
                                       method=name, rmse_N=complete_mean(val), nonfinite_scene_count=bad_scenes(val),
                                       raw_summary_rmse_N=val['summary']['rmse']))
    write_csv('robustness_outcomes.csv', robustness)

    lines = ['# V. 仿真实验：实验设置、已有结果与待补证据', '',
             '本文档按用户要求，将原英文实验草稿更新为中文，并依据原论文实验部分和当前结果重新组织。'
             '本次仅核读、提取和整理已存结果，未运行新的控制仿真。结果有数值依据才填写；待补项不包装成已验证结论。', '',
             '## A. 研究问题与评价依据', '',
             '本研究考察：软组织扫描过程中，如何利用新接触信息及时更新阻抗，在有限在线计算条件下改善实际接触与运动表现。'
             '“快速”需同时看信息时效和在线计算；“近优”需针对明确的控制目标、约束及可实现策略类判断；控制结果还需分别展示精度、冲击和运动响应。', '',
             '本次回查原有 36 篇实验分析，并重点核读 12 篇原文的研究问题与实验／评价部分（附录 A）；'
             '不声称重新完整精读了 36 篇理论与实验全文。接触控制、运动柔顺性、最优控制、快速学习和超声应用论文的指标分别服务其主张，不能整套互相照搬。', '',
             '与直接相关原文相比，正文应形成下面的证据链。该选择是针对本研究的判断，不是任何期刊规定的统一指标清单。', '']
    metric_rows = [
        ('接触精度／保持', '力 RMSE（N）、实际峰值力（N）、零力时长（s）', '主结果表；已有', 'Xue、Song、E27'),
        ('环境变化后的适应速度', '分事件恢复时间（s）及事件对齐响应', '正文动态结果；缺统一测量', 'Song、E27、Xue'),
        ('实际运动行为', '法向实际速度的高／低频 RMS（mm/s），配合轨迹／频谱；存在独立目标时才报位置 RMSE（mm）', '运动结果；缺统一测量', 'Haninger；Jin、Wang 用独立运动任务'),
        ('在线时效', '估计＋决策＋发布的更新时间中位数／P99（ms）、超过控制期限的比例（%）', '正文时效表；目前只有调用计时和局部拆分', 'Matschek、Xue、YOPO'),
        ('决策质量与近优程度', '同一问题下的参考代价差（%）与计算时间；参考解收敛诊断', '独立质量—时间图／表；只有旧 C1 局部研究', 'NNBO 的代价评价、YOPO 的专家代价比较'),
        ('信息为何有效', '估计偏差／检测延迟；去环境运动项、加信息延迟、减慢决策的消融', '机制图／表；已有旧版本证据', 'NNBO、Matschek、E27；与本研究信息链对应'),
        ('适用范围', '各扰动／平台下继续报告同一组主指标，保留失败', '正文边界；已有部分结果', 'E27、Haninger、Wang、AuSoScan'),
    ]
    table(lines, ['研究主张', '指标与单位', '正文位置及现状', '原文依据'], metric_rows)
    lines += ['这里保留代价评价，但不把它冒称物理量。现有调参综合分数 J 混合归一化误差、低接触比例和超力项；'
              '它可以用于已声明的参数选择，却不自动等于完整 B0 的最优控制目标。数值参考、真实刚度规则和全局最优控制器是三个不同概念。', '',
              '运动指标也需与任务匹配：当前切向运动由所有方法共用的独立位置伺服推进，不能用该通道的误差证明阻抗方法改善运动跟踪。'
              '当前可补测的是**实际法向运动与振荡**。Haninger 表 3 的 15 Hz 分带可作为预先固定的参照，仍应结合扫描／呼吸频带和已知激励公开解释；'
              '不能通过事后换滤波器隐藏振荡，也不能把运动越慢直接认作越好。YOPO 的积分平方 jerk 可作补充平顺性量，但并非本接触任务必须堆入主表的指标。', '',
              '## B. 实验设置', '',
              '采用 PC125 位置控制机器人模型，外环频率为 125 Hz（周期 8 ms），传输延迟为一个控制周期，位置伺服时间常数为 19 ms。'
              '力传感器以 1 kHz 采样，经过 30 Hz 一阶低通，目标接触力为 4 N；积分步长为 0.5 ms。接触模型包含柱面有限足迹、'
              '空间变化模量、标准线性固体松弛、Hunt–Crossley 阻尼和摩擦。它是合成工程模型，尚未进行实物组织和机器人标定。', '',
              '所有控制方法采用同一被控对象、可用传感信息和执行接口。共同接口包含虚拟质量、阻尼、弹簧／积分／比例增益及可选参考偏移；'
              '各方法实际调整的子集不同，并保留自身机制和经过说明的任务适配。部署时 Xue 的初始表面先验单独标注，其余真实环境信息仅供数值参照使用。', '',
              '当前 C1 是局部刚度／表面运动估计、短时域直接优化及安全监督切换；C2 是接触历史／锁相估计输入的教师蒸馏网络。'
              'C2 的离线教师可以访问真实刚度，部署时不能访问；离线信息和准备成本并非各方法完全相同。'
              '这些是当前实际实现，不冒称完整 B0 固定节点 HJB／ANEOC、空间预瞄或真实机器人验证。', '']
    table(lines, ['任务', '设置', '目的'], [
        ('T1 主任务', '模量比 30 的刚度阶梯；扫描 8 s；呼吸幅值 1 mm', '快速经过软硬变化时的信息更新与接触调节'),
        ('T3 主任务', '骨样凸起；随机推拉幅值 3–8 mm、过渡 0.5 s；扫描 10 s；呼吸 1 mm', '同时应对空间异质性和组织运动'),
        ('T3 压力测试', '推拉过渡缩短至 0.25 s；另测呼吸幅值 4 mm', '快速扰动及较大生理运动下的失败边界'),
        ('N1／T2 对照', '均匀软／硬组织；另定义均匀组织推拉 T2', '区分适应收益与不必要的适应；N1 有小样本结果，T2 结果未完成'),
    ])
    lines += ['主对比使用训练种子 1000、开发种子 2000 选参，随后冻结控制器；测试种子 3000 和 4000 各含每任务 40 个相同场景及噪声实现。'
              '这两批均已在历史工作中查看。神经方法有三个训练种子，先在同一场景内平均，独立场景数仍为 40。'
              '后续 seed 5000 使用 C1 v3，seed 6000 用于部分平台测试，分别陈述，不能与修正后所有 baseline 的主比较混成一个试验。', '',
              '力指标使用真实仿真力、去除最初 1 s，按 8 ms 日志采样。零力时长计算 F≤0 的时间，是卸载接触代理量；'
              '峰值是采样峰值，不能冒称积分步长下的连续最大值。7 N 是原选参工程阈值，不是临床安全界限。', '',
              '## C. 修正适配并补调参后的主结果', '',
              '修正 Xue 的刚度支路与参数投影，补查伺服衔接；修正 Jin 单步预测在 8 ms 更新时的动作不可见端点；'
              '扩展经典控制对照的参数搜索，并补回学习式对照的参考输出、按任务训练及输出校准。全部修正测试共 1,440 次，未出现积分失败。'
              'C1 v2 与 C2 保持原有参数和权重。以下给出 seed 4000 主表；另一批场景、最大误差、配对区间和全部逐场景数据见原始结果入口。', '']
    for task in ('T1', 'T3'):
        lines += [f'**表 {"1" if task=="T1" else "2"}：{task}，40 个配对场景。**', '']
        rows = []
        for name, label in LABELS.items():
            r = lookup[4000, task, name]
            rows.append([label, f(r['rmse_N_mean']), f(r['peak_force_N_mean']),
                         f(r['no_contact_s_mean']), f(r['worst_peak_force_N'])])
        table(lines, ['方法', '力 RMSE（N）', '场景平均峰值（N）', '平均零力时间（s）', '最坏采样峰值（N）'], rows)
    lines += ['C1 v2 和 C2 的平均力 RMSE 在两批主任务场景中均低于全部竞争性对照。seed 4000 的 T1 最低对照为 Jin 单步 0.1161 N，'
              'C1 v2／C2 为 0.0673／0.0896 N；T3 最低对照为 EtR-A 0.2726 N，C1 v2／C2 为 0.1726／0.2308 N。'
              '但 T3 的 Jin 多步最坏采样峰值为 6.4988 N，低于 C1 v2 的 7.1640 N，因此平均精度收益并非所有物理量全面占优。'
              '神经方法的最坏峰值覆盖 120 个场景／种子运行，不能与单个模型的 40 次极值作等样本数统计推断。', '',
              '以下给出候选相对本批最低平均 RMSE 对照的配对差；参照是描述性选择，区间是条件于已训练模型、以场景为单位的 '
              '10,000 次 bootstrap 边际 95% 区间，未作多重比较校正。负值表示候选误差较低。', '']
    table(lines, ['任务', '候选', '参照', 'RMSE 差（N）', '边际 95% 区间（N）'], [
        [r['task'], LABELS[r['candidate']], LABELS[r['reference']], f(float(r['difference_mean'])),
         f"[{float(r['ci95_low']):.4f}, {float(r['ci95_high']):.4f}]"] for r in primary_pairs])
    lines += [
              '上述比较支持的是当前任务中的接触精度收益。学习式 VIC 是 ES 训练的 E27 启发对照，Haninger 仅移植阻尼规则；'
              '不称其整篇算法完整复现。C1 v2 的主任务配置与 C1 v3 的对应 1 mm 配置一致；不替换为另一个共享配置。'
              'C1 v1 的历史 T3/4000 运行曾出现 46.9164 N 峰值，该负结果仍保留在完整物理报告中。', '',
              '## D. 更新时效与数值决策质量', '',
              '**表 3：已有策略调用计时。数值为各场景调用 P95 的平均值（ms），不是合并全体调用得到的 P95，也不是完整在线链路 P99。**', '']
    table(lines, ['方法', 'T1：调用 P95 均值（ms）', 'T3：调用 P95 均值（ms）'], [
        [label, f(lookup[4000, 'T1', name]['policy_p95_ms_mean'], 5),
         f(lookup[4000, 'T3', name]['policy_p95_ms_mean'], 5)] for name, label in LABELS.items()])
    lines += ['固定增益直接在被控对象中执行，没有单独策略调用记录，故填“—”，不是零耗时。'
              '原候选计时与新增对照计时使用不同并发负载，不能据本表宣称严格的硬件加速比。'
              '全链路的 P99、期限错失率、通信／发布延迟及计算超时引起的真实闭环退化仍缺统一测量。', '',
              '**表 4：旧 C1 的求解器机制试验，开发集每任务 6 场景；独立单进程计时使用其中 3 场景。**', '']
    choices = ['ORD(1 iter)', 'ORD(2 iter)', 'ORD(5 iter)', 'NU(ncorr1,nref2)']
    rows = []
    for nm in choices:
        t1, t3 = sx['T1', nm], sx['T3', nm]
        rows.append([nm, f(t1['reference_gap_percent']), f(t1['mean_policy_ms'], 5),
                     f(t3['reference_gap_percent']), f(t3['mean_policy_ms'], 5)])
    table(lines, ['求解方式', 'T1 参考代价差（%）', 'T1 平均更新（ms）', 'T3 参考代价差（%）', 'T3 平均更新（ms）'], rows)
    lines += ['参考为相同状态、估计、目标和约束下，从上一解和当前解分别执行 25 次 GN 后取较低代价；'
              '它是双初始化多迭代数值参考，不是经过全局最优认证的解。表中间隙也不是力误差百分比。'
              '该试验隔离了旧 C1 内部的数值更新算法，不能直接当作当前 C1 v2／v3、C2 和所有 baseline 的统一近优性比较。', '',
              '普通单次 GN 更新约需 0.11–0.13 ms，在当前载体下没有证据表明求解耗时已成为主要瓶颈。'
              '邻近更新 NU 在 T3 并未比普通单次 GN 更快或更接近该参考；因而当前数据不支持 NU 的独特优势。'
              '若正文继续声称“快速近优”，需补最终实现的质量—时间曲线与参考解收敛诊断；C2 尚无同口径参考代价差。', '',
              '## E. 在线环境信息和决策频率为何有效', '',
              '**表 5：已有信息／频率消融，seed 3000；列中均为力 RMSE（N）。C1 是旧 v1 实例，C2 为原有三个训练种子。**', '']
    table(lines, ['机制配置', 'T1', 'T3'], [[label, f(ax['T1', label]), f(ax['T3', label])]
          for label in list(c1_names.values()) + list(c2_names.values())])
    lines += ['移除表面运动估计使 C1 的 T1／T3 RMSE 从 0.0672／0.1633 N 增至 0.7563／1.2417 N，'
              '显示组织运动与材料变化的分离对该实例具有实际控制作用。减慢 C1、C2 的决策频率均使两个主任务误差上升。'
              '但是 C1 的环境信息额外延迟 100 ms 在 T3 未增加平均 RMSE（0.1633→0.1610 N），T1 则变差（0.0672→0.1152 N）；'
              '信息时效的作用随任务而异，不能概括成任意信息越新越好。', '',
              '旧 C1 的 10 场景开发集估计诊断中，T1 首次骨样区域检测延迟中位数为 48 ms，移除表面运动项后为 368 ms；'
              'T3 分别为 0／144 ms。这些值已存于 JSON，但本次未在可执行来源中定位该诊断的生成代码和完整检测判据，'
              '因此仅作历史诊断记录，不把它们直接作为通用适应时间的正文主证据。'
              '完整消融和估计偏差均已存在，但最终 C1 v2／v3 的关键消融、统一事件恢复时间和运动振荡指标仍需补测。', '',
              '## F. 已有压力测试与适用边界', '',
              '**表 6：seed 5000，40 场景／任务。单位 N，均为力 RMSE；此处是 C1 v3，与前述修正后所有论文对照的主表分开。**', '']
    cell_labels = {'T1': 'T1，呼吸 1 mm', 'T1_4': 'T1，呼吸 4 mm', 'T3': 'T3，0.5 s／1 mm',
                   'T3_4': 'T3，0.5 s／4 mm', 'T3fast': 'T3，0.25 s／1 mm', 'T3fast4': 'T3，0.25 s／4 mm'}
    rows = []
    for cell in ('T1', 'T1_4', 'T3', 'T3_4', 'T3fast', 'T3fast4'):
        m = gaps['c2_seed5000'][cell]['methods']
        rows.append([cell_labels[cell]] + [shown_mean(m[n]) for n in ('B1', 'EtR-A', 'C1v3', 'C2')])
    table(lines, ['任务条件', 'B1', 'EtR-A', 'C1 v3', 'C2'], rows)
    lines += ['C1 v3 在这些场景中保持较低平均误差；C2 在快速推拉下明显退化，T3fast 为 1.0483 N，高于 B1 的 0.6822 N。'
              '因此不能把两个主任务的结果扩展为所有扰动条件下的普遍优势。该表未包含已修正并重新调参的全部论文 baseline；'
              '它支持候选的压力边界，不支持“全部对照上的泛化领先”。', '',
              '**表 7：已完成的平台压力测试，seed 6000。B1 为平台上重新调优；C1 在非标称平台也重新调优；C2 为原 PC125 模型直接迁移。单位 N，均为 RMSE。**', '']
    rows = []
    for platform in ('PC125', 'PC125_lag0.5', 'PC125dead', 'PC500'):
        for task in ('T1', 'T3'):
            m = gaps['platform_seed6000'][platform][task]['methods']
            cn = 'C1v3_tr' if platform == 'PC125' else 'C1v3_rt'
            rows.append([platform, task, shown_mean(m['B1']), shown_mean(m[cn]), shown_mean(m['C2_tr'])])
    table(lines, ['平台', '任务', 'B1：平台调优', 'C1 v3：标称／平台调优', 'C2：直接迁移'], rows)
    lines += ['PC500 的平台调优 B1 已达到 T1／T3 的 0.0414／0.0576 N，低于平台调优 C1 的 0.0482／0.0701 N。'
              'C2 的直接迁移结果也很差，不能隐藏。但跨采样周期后，固定历史窗口的物理时长、激励和归一化语义也改变；'
              '这些结果是未重新训练／校准的迁移压力结果，不能代替公平平台适配后的学习算法排名。'
              'PC125dead 的 C2 在 T3／T1 分别有 1／2 个场景的种子平均值含积分失败，整组 RMSE 因而记 NA；'
              '原始汇总曾忽略非有限值而给出有限均值，这里不把它当成完整 40 场景结果。'
              '因此现有证据将收益限定在特定信息、执行平台与更新条件下。', '']
    n1 = read('candidate1_v3_results.json')['known_behaviour']['C1v3_T3']['soft_15kPa']
    n1_c = sum(r['rmse'] for r in n1)/len(n1)
    n1_b = sum(r['rmse_B1'] for r in n1)/len(n1)
    lines += [f'均匀软组织的已有负对照也未显示普遍收益：T3 配置迁移到 15 kPa 均匀组织的 {len(n1)} 个场景时，'
              f'C1 v3 平均力 RMSE 为 {n1_c:.4f} N，比较的 B1 为 {n1_b:.4f} N。'
              '这只是小样本配置迁移检查，不能替代完成 T2 推拉任务及其独立调参比较。', '',
              '## G. 需要补哪些实验', '',
              '已有结果不需要重新发明一套总分；应补齐下列明确证据。前三项直接对应当前完整研究主张，后两项取决于希望声称的泛化范围。', '']
    needs = [
        ('1．最终方法的实际运动与全链路计时', '同一 PC125 场景重放所有最终控制器，保存高频法向位置／速度和逐调用／发布原始时间戳；统一硬件、单进程计时与预热。报告高／低频速度 RMS、轨迹／频谱、全链路中位数／P99、8 ms 期限错失率。', '目前只有力主表与调用计时；不需要先改横向模型就能补法向运动。若要声称横向跟踪收益，另需可受控制器影响的耦合运动任务。'),
        ('2．孤立环境变化的恢复与信息时效', '分别施加刚度变化、表面移动和组合事件，使用冻结的最终 C1／C2 及对照；测事件对齐响应、恢复时间，扫描信息延迟与决策周期。原有 ±0.5 N 容差带可沿用，持续保持时长须在看结果前固定。', '已有旧版单点延迟／频率消融，缺最终版的系统曲线和统一恢复时间；不得把未恢复事件丢掉。'),
        ('3．最终方法的数值决策质量—计算时间', '固定相同状态、可用估计、目标和可实现动作集合，构建多初始化、多迭代且有收敛检查的数值参考；给预算—参考代价差及对应闭环物理结果。C2 若动作类不同，须统一动作类或分别声明参考类，不能直接混排差值。', '旧 C1 有局部参考试验，C2 没有同口径结果；普通一步求解已很快，需直接检验效率收益而非预设瓶颈。'),
        ('4．统一最终对照的独立场景与边界', '在新的未查看场景中统一最终候选和已修正的论文 baseline；必要时覆盖快推拉／大呼吸。完成平台 lag2、额外两周期延迟和 T2 均匀组织推拉后再声称相应泛化。', 'seed 5000／6000 有候选和部分基线，但不是全部修正对照；lag2 仅部分选参，两周期延迟测试和 T2 没有结果。'),
        ('5．应用与理论主张需要独立证据', '只有声称图像／重建收益时才补超声或声学模型／实物指标；只有声称被动性或能量保证时才补适用执行模型下的理论与能量核查。', '当前模型无图像、关节力矩和执行器功率，不能编造成像质量、真实能耗或实物安全结论；这些不是当前接触控制表必须硬加的列。'),
    ]
    table(lines, ['待补试验', '具体操作和产出', '现有基础与缺口'], needs)
    lines += ['## 附录 A：逐篇“研究问题—实验任务—评价指标”核读映射', '',
              '下列核读记录是本次选择指标的原始依据。正文仅保留与自身主张相匹配的指标；不把 12 篇所有指标的并集当作投稿清单。', '']
    for name, path, ranges, question, task, metrics, transfer in PAPERS:
        line = ranges[-1][0]
        lines += [f'**{name}**（[原文实验段](../literature/extracted_text/{path}#L{line})；核读段落行号 {ranges}）', '',
                  f'- 研究问题：{question}。', f'- 实验任务：{task}。',
                  f'- 原文评价：{metrics}。', f'- 本研究如何使用：{transfer}。', '']
    lines += ['## 附录 B：来源与复算', '',
              '- [修正后完整物理结果](../experiments/cobot_scan_benchmark/refined_baselines/results.md)包含两个主测试批次、配对边际区间及原始逐场景数据。',
              '- [正文表格数据](../experiments/cobot_scan_benchmark/paper_tables/)由原始 JSON／CSV 直接提取；不混用不同版本、平台和场景的数值。',
              '- [candidate1_results.json](../experiments/cobot_scan_benchmark/candidate1_results.json)：旧 C1 消融、求解器质量—耗时和估计诊断。',
              '- [candidate2_results.json](../experiments/cobot_scan_benchmark/candidate2_results.json)：原 C2 消融和训练记录。',
              '- [gaps_results.json](../experiments/cobot_scan_benchmark/gaps_results.json)：新同步的 seed 5000／6000 结果；部分平台未完成。',
              '- 复算：`python -X utf8 experiments/cobot_scan_benchmark/build_paper_experiments.py`。', '',
              '主表的配对置信区间是以场景为单位、条件于已训练模型的边际区间，未作多重比较校正。'
              '旧消融和求解器表主要给描述量，不另制造显著性星号。全部论文原件、抽取文本、已有控制代码／权重和 B0 冻结快照均未改动。', '']

    # Explicit provenance makes re-extraction verifiable without copying source data.
    sources = ['refined_baselines/summary.csv', 'refined_baselines/paired_differences.csv',
               'candidate1_results.json', 'candidate1_v3_results.json', 'candidate2_results.json', 'gaps_results.json',
               'candidate1_model.py', 'candidate2_learned.py', 'gaps_eval.py', 'plant.py']
    hashes = {str((B / n).relative_to(ROOT)): hashlib.sha256((B / n).read_bytes()).hexdigest() for n in sources}
    literature = []
    for name, path, ranges, *_ in PAPERS:
        src = ROOT / 'literature/extracted_text' / path
        assert max(end for _, end in ranges) <= len(src.read_text(encoding='utf-8').splitlines())
        literature.append(dict(paper=name, source=str(src.relative_to(ROOT)), read_ranges=ranges,
                               sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
    audit = dict(main_groups=len(main_rows), solver_groups=len(solvers), ablation_groups=len(abls),
                 robustness_groups=len(robustness), papers_with_experiment_passages_read=len(PAPERS),
                 new_control_simulations_run=False, incomplete_platform_cells=incomplete,
                 missing_declared_platforms=['PC125_dead2'],
                 nonfinite_scene_groups=[{k: r[k] for k in ('condition', 'task', 'method', 'nonfinite_scene_count')}
                                         for r in robustness if r['nonfinite_scene_count']],
                 t2_results_present='T2_seed6000' in gaps, sources=hashes, literature=literature,
                 historical_diagnostic_generation_unlocated=['candidate1_results.json:estimator_bias_dev'],
                 generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'source_manifest.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8', newline='\n')
    # The evidence audit precedes drafting; the reader sees setup first.
    ia, ib, ic = (lines.index(h) for h in ('## A. 研究问题与评价依据', '## B. 实验设置', '## C. 修正适配并补调参后的主结果'))
    lines = lines[:ia] + ['## A. 实验设置'] + lines[ib+1:ic] + ['## B. 研究问题与评价依据'] + lines[ia+1:ib] + lines[ic:]
    DRAFT.write_text('\n'.join(lines), encoding='utf-8', newline='\n')
    print(json.dumps({k: v for k, v in audit.items() if k not in ('sources', 'literature')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
