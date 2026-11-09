"""Huber-D4 against D4 under one configuration: same window (24 ticks), fc, gate and sigma scales.

Three arms each (gated, ungated, ungated + IMM), the three tables' scoring functions, nothing retuned.
Replay observations come from the gated Huber-D4 closed loop, as in the main tables.
Predictions fixed before running, from huber_prior_check.py: at equal window, plateau raw errors identical, event differences
only in cough and combo, B identical. Tick data goes to --ticks outside Git.
Writes results/huber_unified_compare.json and huber_unified_tables.md.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import identification_completion as IC
import overlay_compare as O

HERE = Path(__file__).resolve().parent
ARMS = ('Huber-D4', 'Huber-D4_nogate', 'Huber-D4_imm', 'D4', 'D4_nogate', 'D4_imm')
WINDOWS = (24,)
TICKS = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE/'results/huber_unified_ticks'
ORIGINAL = O.chosen


def configure(nw):
    configs = ORIGINAL()
    same = dict(nw=nw, fc=14., gate=.551, sigma_scale=.796)
    configs['Huber-D4'] = dict(same); configs['D4'] = dict(same)
    O.chosen = lambda: configs


def closed(task):
    nw, method, case = task
    configure(nw)
    row = CC.job((method, case))
    return dict(nw=nw, method=method, index=case['index'], group=case['group'], **{k: row[k] for k in ('L', 'peak', 'force_rmse', 'loss_time', 'long_loss')})


def replay(task):
    nw, source, case = task
    configure(nw)
    IC.TICKS = TICKS/str(nw)
    return IC.case_job((source, case, ARMS))


def score(path):
    return IC.score_case(path)


def table1(rows):
    head = ['窗口', '方法', 'L 均值', 'L 中位', 'A', 'A慢薄', 'B', 'A+B', '咳嗽', '综合', 'A_large', 'RMSE 中位', '最大力', '>15 N', '>36.5 N', '丢失时间', '长失联']
    out = [IC.md_row(head), '|'+'---|'*len(head)]
    for nw in WINDOWS:
        for m in ARMS:
            rr = [r for r in rows if r['nw'] == nw and r['method'] == m]
            L = lambda g: float(np.mean([r['L'] for r in rr if r['group'] == g]))
            out.append(IC.md_row([nw, IC.LABELS[m], IC.f(float(np.mean([r['L'] for r in rr]))), IC.f(float(np.median([r['L'] for r in rr])))]+[IC.f(L(g)) for g in CC.SCENES]+[
                IC.f(float(np.median([r['force_rmse'] for r in rr]))), IC.f(max(r['peak'] for r in rr), 1), IC.pct(np.mean([r['peak'] > 15. for r in rr])), IC.pct(np.mean([r['peak'] > 36.5 for r in rr])),
                IC.pct(float(np.mean([r['loss_time'] for r in rr]))/9., 2), IC.pct(np.mean([r['long_loss'] for r in rr]))]))
    return out


def table2(seg):
    head = ['窗口', '方法', '通道', '|误差| 中位', 'P90', '偏差', '组内波动', '低于 0.85', '接受率', '错误更新', '覆盖率']
    out = []
    for s, title in (('plateau', '平稳段'), ('transition', '刚度变化段'), ('event', '突发扰动段')):
        out += ['', f'**{title}**', '', IC.md_row(head), '|'+'---|'*len(head)]
        for nw in WINDOWS:
            for m in ARMS:
                x = next(i for i in seg if i['nw'] == nw and i['method'] == m and i['segment'] == s and i['group'] == 'all')
                for c, label in (('raw', '原始'), ('k_good', '控制用')):
                    if c == 'k_good' and m not in ('Huber-D4', 'D4'):
                        continue
                    p = x[c]
                    extra = [IC.pct(x['accepted']) if m in ('Huber-D4', 'D4') else '无门控', IC.pct(x['false_update'], 1), IC.pct(x['coverage'])] if c == 'raw' else ['—']*3
                    out.append(IC.md_row([nw, IC.LABELS[m], label, IC.f(p.get('abs_median')), IC.f(p.get('abs_p90')), IC.f(p.get('bias')), IC.f(p.get('within_case_sd')), IC.pct(p.get('below_085'))]+extra))
    for s, title, groups in (('transition', '刚度变化段', ('A', 'A_slow_thin', 'A+B', 'cough', 'combo', 'A_large')), ('event', '突发扰动段', ('B', 'A+B', 'cough', 'combo'))):
        h = ['窗口', '方法']+list(groups)
        out += ['', f'**{title}，原始输出 |误差| 中位 / P90**', '', IC.md_row(h), '|'+'---|'*len(h)]
        for nw in WINDOWS:
            for m in ARMS:
                cells = []
                for g in groups:
                    x = next(i for i in seg if i['nw'] == nw and i['method'] == m and i['segment'] == s and i['group'] == g)
                    cells.append(f"{IC.f(x['raw'].get('abs_median'))} / {IC.f(x['raw'].get('abs_p90'))}")
                out.append(IC.md_row([nw, IC.LABELS[m]]+cells))
    return out


def table3(track):
    head = ['窗口', '方法', '方向', '事件数', '变化段误差中位', '最大低估中位', '最大低估最差', '恢复中位 (s)', '已恢复', '未恢复', '下次变化截断', '回合结束截断']
    out = []
    for family, title in (('no_event', '既有场景，|Δln k|≥0.3'), ('no_event_small', '既有场景，0.15≤|Δln k|<0.3'), ('A_large_no_event', 'A_large，|Δln k|≥0.3')):
        for c, label in (('raw', '原始输出'), ('k_good', '控制用 k_good')):
            out += ['', f'**{title}，{label}**', '', IC.md_row(head), '|'+'---|'*len(head)]
            for nw in WINDOWS:
                for m in ARMS:
                    if c == 'k_good' and m not in ('Huber-D4', 'D4'):
                        continue
                    for d, dl in (('rising', '上升'), ('falling', '下降')):
                        x = next(i for i in track if i['nw'] == nw and i['method'] == m and i['direction'] == d and i['family'] == family)
                        p = x[c]
                        rec = '—' if x['n'] == 0 else '未达到' if p['recovery_median'] is None else IC.f(p['recovery_median'], 2)
                        out.append(IC.md_row([nw, IC.LABELS[m], dl, x['n'], IC.f(p['change_mae']['median']), IC.f(p['under_raw']['median']), IC.f(p['under_raw']['max']), rec, p['recovered'], p['unrecovered'], p['truncated_next_change'], p['truncated']]))
    return out


def main():
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    cases = CC.cases()
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(closed, [(nw, m, c) for nw in WINDOWS for m in ARMS for c in cases], chunksize=1))
    worst = max(abs(r['L']-stored[(r['method'], r['index'])]['L']) for r in rows if r['nw'] == 24 and r['method'].startswith('Huber'))
    assert worst == 0., worst
    print('closed done; Huber at 24 equals stored', flush=True)
    for nw in WINDOWS:
        (TICKS/str(nw)).mkdir(parents=True, exist_ok=True)
    report = dict(config='fc 14, gate .551, sigma_scale .796; nw as listed', closed=rows, replay={})
    seg_all, track_all = [], []
    for source in ('Huber-D4',):
        for nw in WINDOWS:
            for p in (TICKS/str(nw)).glob('*.npz'):
                p.unlink()
        with ProcessPoolExecutor(8) as ex:
            list(ex.map(replay, [(nw, source, c) for nw in WINDOWS for c in cases], chunksize=1))
        print('replay done', source, flush=True)
        for nw in WINDOWS:
            with ProcessPoolExecutor(8) as ex:
                results = list(ex.map(score, sorted((TICKS/str(nw)).glob('*.npz')), chunksize=1))
            seg, track, _ = IC.summarize([r for r, _ in results], [s for _, s in results], source, ARMS)
            for item in seg+track:
                item['nw'] = nw
            report['replay'][f'{source}_{nw}'] = dict(segments=seg, tracking=track)
            if source == 'Huber-D4':
                seg_all += seg; track_all += track
    md = ['# Huber-D⁴ 与 D⁴ 同配置三张表', '', '两者相同：fc 14 Hz，门控尺度 .551，σ 尺度 .796，窗口如表所列。未重新调参。共同观测来自同窗口的 Huber-D⁴ 门控闭环。42 例、每组 6 例，只作描述。', '', '## 表 1：闭环', '']+table1(rows)+['', '## 表 2：辨识精度']+table2(seg_all)+['', '## 表 3：变化跟踪']+table3(track_all)
    (HERE/'huber_unified_tables.md').write_text('\n'.join(md)+'\n', encoding='utf-8')
    (HERE/'results/huber_unified_compare.json').write_text(json.dumps(report, default=float), encoding='utf-8')
    print('tables written', flush=True)


if __name__ == '__main__':
    main()
