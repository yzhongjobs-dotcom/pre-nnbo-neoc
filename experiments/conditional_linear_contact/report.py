"""把已完成的 JSON 结果整理成可复核的中文报告。"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
from closed_loop import summarize  # noqa: E402

RES = HERE / 'results'


def load(name):
    return json.loads((RES / name).read_text(encoding='utf8'))


def f(x, n=3):
    return 'NA' if x is None else f'{float(x):.{n}f}'


def main():
    tuned = load('closed_loop_tuned_lag1.json')
    default = load('closed_loop_default_lag1.json')
    tuned_summary = summarize(tuned['rows'], 1.)
    default_summary = summarize(default['rows'], 1.)
    theory = load('theory_checks.json')
    open_loop = load('open_loop.json')['rows']
    mismatch = load('mismatch.json')
    params = load('tuned_params.json')

    lines = [
        '# 条件线性接触 V0 实验报告', '',
        '本报告只使用诊断集 `110000+` 和调参集 `190000+`；没有访问 `130000+` 确认集。',
        '被控对象、已有实验代码和已有结果未修改。软件计时不含通信，不是实时认证。', '',
        '## 轮次结论', '',
        '| 轮次 | 改动/目的 | 结果 | 判断 |',
        '|---|---|---|---|',
        '| 1 | 默认参数试通，每族 2 场景；核对启动、IMM、失联接口 | 启动段 10/10 场景逐值相同，最大命令差 0；部分 S2/T3m 明显失控 | 接口可运行，但默认后验进入反馈后有尾部风险 |',
        '| 2 | 固定设计下完成 61 次调参；五族各 6 个 190000+ 场景 | 最优目标 0.762；调参参数在 `q_lnk1/gain_rel_max/c` 上界 | 中位收益可能来自宽参数，需诊断集检验，不能直接采纳 |',
        '| 3 | 在诊断集五族各 50 场景上复核调参点，并跑默认参数作为设计对照 | 调参点均值和峰值尾部恶化；默认点更差，接受率下降 | 当前 V0 尚未证明比两个对照更 work；问题定位到后验—反馈接口/模型失配，未再增加调参预算 |',
        '', '## 名义闭环：调参点，诊断集 50 场景/族', '',
        '| 族 | 方法 | RMSE 中位/均值 N | 峰值>50 N | 长失联>0.5 s | 接受率 | 调用 P50/P99 ms |',
        '|---|---|---:|---:|---:|---:|---:|',
    ]
    for x in tuned_summary:
        if 'method' not in x:
            continue
        lines.append(f"| {x['family']} | {x['method']} | {f(x['rmse_median'])}/{f(x['rmse_mean'])} | {x['peak_gt_50']} | {x['long_loss_gt_05']} | {f(x['acceptance_rate'])} | {f(x['call_ms_P50'],2)}/{f(x['call_ms_P99'],2)} |")
    lines += ['', '配对 bootstrap 95% 区间（候选−V1+Proposed，RMSE N）：', '',
              '| 族 | joint_lean | conditional_v0 |', '|---|---:|---:|']
    for fam in sorted({x['family'] for x in tuned_summary if 'comparison' in x}):
        a = next(x for x in tuned_summary if x.get('family') == fam and x.get('comparison') == 'joint_lean-v1_proposed')
        b = next(x for x in tuned_summary if x.get('family') == fam and x.get('comparison') == 'conditional_v0-v1_proposed')
        lines.append(f"| {fam} | {f(a['delta_rmse_mean'])} [{f(a['delta_rmse_bootstrap95'][0])},{f(a['delta_rmse_bootstrap95'][1])}] | {f(b['delta_rmse_mean'])} [{f(b['delta_rmse_bootstrap95'][0])},{f(b['delta_rmse_bootstrap95'][1])}] |")

    lines += ['', '## 默认参数第 3 轮诊断摘要', '',
              '| 族 | V1 中位/均值 | joint_lean 中位/均值 | conditional_v0 中位/均值 | 候选峰值>50 |',
              '|---|---:|---:|---:|---:|']
    for fam in sorted({x['family'] for x in default_summary if 'method' in x}):
        a = {x['method']: x for x in default_summary if x.get('family') == fam and 'method' in x}
        lines.append(f"| {fam} | {f(a['v1_proposed']['rmse_median'])}/{f(a['v1_proposed']['rmse_mean'])} | {f(a['joint_lean']['rmse_median'])}/{f(a['joint_lean']['rmse_mean'])} | {f(a['conditional_v0']['rmse_median'])}/{f(a['conditional_v0']['rmse_mean'])} | {a['conditional_v0']['peak_gt_50']} |")

    lines += ['', '## 失配：伺服滞后 ×0.6、×1.5，各族 20 场景', '',
              '| 滞后 | 族 | 方法 | RMSE 中位/均值 | 峰值>50 N |', '|---:|---|---|---:|---:|']
    for lag in (.6, 1.5):
        sm = summarize(mismatch['rows'], lag)
        for x in sm:
            if 'method' in x:
                lines.append(f"| {lag:g} | {x['family']} | {x['method']} | {f(x['rmse_median'])}/{f(x['rmse_mean'])} | {x['peak_gt_50']} |")

    lines += ['', '## 开环回放与零噪声孪生', '',
              '| 分组 | 方法 | 结构误差 RMSE | 噪声误差 RMSE | 零噪声 95% 覆盖 | 接受率 |', '|---|---|---:|---:|---:|---:|']
    for group in ('calm', 'stiffness_change', 'fast_motion'):
        for method in ('projection16', 'lean32', 'conditional_v0'):
            a = [r[method] for r in open_loop if r.get('groups') == group and method in r]
            if not a:
                continue
            lines.append(f"| {group} | {method} | {f(np.nanmean([x['structure_rmse'] for x in a]))} | {f(np.nanmean([x['noise_rmse'] for x in a]))} | {f(np.nanmean([x['twin_coverage95'] for x in a]))} | {f(np.nanmean([x['acceptance_rate'] for x in a]))} |")

    s = theory['synthetic']['results']
    lines += ['', '## 理论与合成检验', '',
              f"- 合成闭环（40 次、450 拍、含小激励）：无反馈归一化误差方差 {f(s['no_feedback']['normalized_error_variance'])}、95% 覆盖 {f(s['no_feedback']['coverage95'])}；有反馈分别为 {f(s['feedback']['normalized_error_variance'])}、{f(s['feedback']['coverage95'])}。该检验支持“在一致增量模型下条件后验可校准”，不支持真实接触对象上的普适保证。",
              '- 频率响应、三阶 Butterworth 比较、`k` 通道的创新积分等数值表见 `theory_checks.json`；它们是离散滤波器的数值比较，不是闭式等价证明。',
              '- 无限扰动先验的有限窗口极限数值趋于多项式投影；成立条件已写入 JSON：多项式满秩、Y 不在扰动子空间、加权最小二乘模型且先验方差趋于无穷。',
              '- q_s/q_m 与 k=1e3/1e4/1e5 的转折频率表见 `theory_checks.json`。',
              '', '## 参数与未解决问题', '',
              f'- 调参最优参数：`{json.dumps(params, ensure_ascii=False)}`。调参严格使用 61 次，不因迭代增加。',
              '- 当前 V0 尚未被采纳：它在部分族的中位数改善，但诊断集均值、峰值尾部、长失联和失配下不稳定，完整任务收益没有超过 V1+Proposed 和 joint_lean。',
              '- 尚未解决：真实接触模型失配时，四状态后验的一致性门与保守反馈规则之间的接口；低信息/刚度突变时如何避免把合法“只能上调”的反馈路径变成失稳尾部；以及更严格接受率与任务收益的权衡。',
              '- 没有做：确认集 130000+、V1（D/r/a 前馈）接口、真机验证；这些不应从本轮结果推断。',
              '']
    report = '\n'.join(lines)
    (HERE / 'results' / 'final_report.md').write_text(report, encoding='utf8')
    notes = ('# 迭代记录\n\n'
             '1. 默认试通：启动段逐值一致，但默认滤波后验进入反馈后出现尾部失稳；不把旧错误实现算作负结果。\n'
             '2. 61 次调参：最优点降低调参目标，但参数触及宽模式/接受带上界；诊断集暴露尾部风险。\n'
             '3. 诊断复核：调参点和默认点均未达到可采纳的完整闭环收益。后续应先修正后验—反馈接口和模型失配问题，不能追加本轮调参预算或访问确认集。\n')
    (HERE / 'results' / 'iteration_notes.md').write_text(notes, encoding='utf8')
    print(HERE / 'results' / 'final_report.md')


if __name__ == '__main__':
    main()
