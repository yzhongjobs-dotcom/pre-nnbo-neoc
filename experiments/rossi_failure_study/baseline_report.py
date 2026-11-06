"""Recompute all report tables from saved primary JSON/NPZ, no simulations."""
import json
import math
from pathlib import Path

import numpy as np

from baseline_compare import OUT,HERE,METHODS,ID_METHODS,key,load

REPORT=Path('D:/pip_env/temp/claude/8fbcb2674095e6be/841f52cd-14a8-4db4-b56f-36be6ffcc337/scratchpad/baseline_report.md')


def finite(v):return np.array([x for x in v if x is not None and np.isfinite(x)],float)
def median(v):
    a=finite(v);return float(np.median(a)) if len(a) else None
def mean(v):
    a=finite(v);return float(np.mean(a)) if len(a) else None
def fmt(x,d=3):
    if x is None or not np.isfinite(x):return '—'
    if d==3 and 0<abs(x)<.0005:return f'{x:.2g}'
    return f'{x:.{d}f}'
def stat(rows):
    return dict(n=len(rows),L=mean([math.log1p(r['J']) for r in rows]),rmse_med=median([r['rmse'] for r in rows]),
                rmse_mean=mean([r['rmse'] for r in rows]),event_med=median([r['event_rmse'] for r in rows]),
                event_mean=mean([r['event_rmse'] for r in rows]),over15=sum((r['peak'] or 0)>15 for r in rows),
                over36=sum((r['peak'] or 0)>36.5 for r in rows),long=sum(r['long_loss'] for r in rows),
                fail=sum(r['failed'] for r in rows),nonfinite=sum(r.get('nonfinite_commands',0) for r in rows),
                peak=max(finite([r['peak'] for r in rows]),default=np.nan),
                calls50=median([r['call_p50'] for r in rows]),calls99=median([r['call_p99'] for r in rows]),
                sat=sum(r.get('saturated_commands',0) for r in rows))


def paired(rows,base):
    d=np.array([math.log1p(r['J'])-math.log1p(base[key(r)]['J']) for r in rows])
    if not len(d):return dict(n=0,mean=None,ci=None,stratified_ci=None)
    rng=np.random.default_rng(574501)
    boot=d[rng.integers(0,len(d),(5000,len(d)))].mean(1)
    # Retain the protocol's condition/dynamics/phase cell sizes in each draw.
    strata={}
    for j,r in enumerate(rows):strata.setdefault((r['cond'],r['tissue'],r['dynamics'],r['phase']),[]).append(j)
    st=np.zeros(5000)
    for ix in strata.values():
        values=d[ix];st+=values[rng.integers(0,len(ix),(5000,len(ix)))].sum(1)/len(d)
    clusters={}
    for j,r in enumerate(rows):clusters.setdefault((r['cond'],r['tissue'],r['i']),[]).append(j)
    # Shared physical scene draws across dynamics/phase are not independent.
    sums=np.array([d[ix].sum() for ix in clusters.values()]);sizes=np.array([len(ix) for ix in clusters.values()])
    draws=rng.integers(0,len(clusters),(5000,len(clusters)))
    clustered=sums[draws].sum(1)/sizes[draws].sum(1)
    return dict(n=len(d),clusters=len(clusters),mean=float(d.mean()),ci=np.percentile(clustered,[2.5,97.5]).tolist(),
                ordinary_ci=np.percentile(boot,[2.5,97.5]).tolist(),stratified_ci=np.percentile(st,[2.5,97.5]).tolist())


def difftext(a):
    signed=lambda x:f'{x:+.2g}' if 0<abs(x)<.0005 else f'{x:+.3f}'
    return '—' if a['mean'] is None else f"{signed(a['mean'])} [{signed(a['ci'][0])}, {signed(a['ci'][1])}]"


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,r))+' |' for r in rows])+'\n'


def module_tables(module):
    rows=module['rows'];summary={};text=[]
    identity=[r['source_command_max_error'] for r in rows if r.get('source_command_max_error') is not None]
    assert len(identity)==336 and max(identity)<1e-10
    summary['source_command_reproduction']=dict(n=len(identity),max_error=max(identity))
    text.append(f'336 份源策略自身回放均通过命令一致性检查，最大命令差 {max(identity):.3g} m；使用回调观测而非把最终日志最新行直接当作实时观测。模块中的 Huber-D4 是本轮调后配置，σ/χ 尺度见 B，不能把其覆盖率直接当作当前默认配置的覆盖率。')
    text+=['模块误差是 |log(khat/k_true)|；有限率、门控可用率分母包含所有 warm 后且真实 k>200 的采样时刻。“接受”指 accept_trace 中精度和 χ 基础门控通过；实际更新前馈 k_good 还有公共 rel_meas/bound 等规则，不等于每个通过窗口都会替换 k_good。静止列按协议称 quiet，实际仍有呼吸/滑动，不能解读成完全静止实验。覆盖率为被接受窗口中满足题定 `|log(khat/k)|≤2σ/khat` 的比例；不是默认 95% 的保证。']
    for source in ('Huber-D4','D4'):
        text+=['\n### 源日志：'+source+'\n']
        tab=[]
        for method in ID_METHODS:
            rr=[r for r in rows if r['source']==source and r['method']==method]
            for kind in ('quiet','event','changing'):
                restricted=[dict(r,groups={kind:r['groups'][kind]}) for r in rr]
                # Inline aggregation without assuming the removed groups exist.
                s=aggregate_group(restricted);summary[f'{source}:{method}:{kind}']=s
                tab.append([method,kind,s['n'],fmt(s['valid']),fmt(s['available']),fmt(s['error_med'])+'/'+fmt(s['error_p90']),fmt(s['accepted_med']),fmt(s['coverage']),fmt(s['t50']*1e6 if s['t50'] else None,0)+'/'+fmt(s['t99']*1e6 if s['t99'] else None,0)])
        text.append(table(['估计器','片段','采样数','有限率','可用率','误差中位/P90','接受误差','覆盖率','估计调用 p50/p99 μs'],tab))
    text+=['\n### 按工况（两个日志来源合并，完整分来源统计见 JSON）\n']
    tab=[]
    for method in ID_METHODS:
        for cond in ('breath','shift_fast','flinch','cough','rib','combo'):
            rr=[r for r in rows if r['method']==method and r['cond']==cond]
            s=aggregate_group(rr);summary[f'all:{method}:{cond}']=s
            tab.append([method,cond,fmt(s['available']),fmt(s['accepted_med']),fmt(s['coverage']),str(s['recovered'])+'/'+str(s['recovery_n']),fmt(s['recovery_med'])])
    text.append(table(['估计器','工况','可用率','接受误差','覆盖率','入界/过渡数','已入界延迟中位 s'],tab))
    text+=['未入界的过渡已单独计入分母，上表的延迟中位只针对已入界者，会乐观；原始 JSON 保留每次右删失时长。真实 k 本身不是阶跃，且相邻事件/刚度过渡可重叠；自然过渡统计没有要求过渡前已经准确，尤其 Xue 的事后入界不一定代表主动跟踪恢复。因此这只描述自然过渡响应，带前置准确性资格的人工阶跃见附加检查。']
    return '\n'.join(text),summary


def aggregate_group(rows):
    n=sum(g['n'] for r in rows for g in r['groups'].values())
    accepted=sum(g['accepted'] for r in rows for g in r['groups'].values());valid=sum(g['valid'] for r in rows for g in r['groups'].values())
    covered=sum(g['covered'] for r in rows for g in r['groups'].values());errors=[];accerr=[];times=[]
    for r in rows:
        with np.load(OUT/r['samples_file']) as samples:
            for kind in r['groups']:
                errors.extend(samples[kind+'_error_samples']);accerr.extend(samples[kind+'_accepted_errors']);times.extend(samples[kind+'_times'])
    rec=[v for r in rows for v in r.get('recovery',[])];delays=[v['delay'] for v in rec if v['delay'] is not None]
    return dict(n=n,valid=valid/n if n else None,available=accepted/n if n else None,error_med=median(errors),
                error_p90=float(np.percentile(errors,90)) if len(errors) else None,accepted_med=median(accerr),
                coverage=covered/accepted if accepted else None,t50=median(times),t99=float(np.percentile(times,99)) if len(times) else None,
                recovery_n=len(rec),recovered=len(delays),recovery_med=median(delays))


def validate_primary(dev,test,module,diagnostic):
    """Check budgets, frozen choices, paired cells and retained evidence."""
    canonical=lambda cfg:json.dumps(cfg,sort_keys=True)
    for method in METHODS:
        score=dev['scores'][method];configs=score['configs']
        assert len(configs)==len({canonical(c) for c in configs})==24
        episodes=[r for r in dev['rows'] if r['method']==method];values=[]
        for cfg in configs:
            rr=[r for r in episodes if r['cfg']==cfg]
            assert len(rr)==len({key(r) for r in rr})==8
            values.append(float(np.mean([math.log1p(r['J']) for r in rr])))
        assert np.allclose(values,score['all'],atol=1e-12,rtol=0)
        assert int(np.argmin(values))==score['index']
        assert dev['chosen'][method]==configs[score['index']]
    cell_sets=[]
    for method,cfg in test['chosen'].items():
        rr=[r for r in test['rows'] if r['method']==method]
        assert all(r['cfg']==cfg for r in rr)
        cell_sets.append({key(r) for r in rr})
    assert all(s==cell_sets[0] for s in cell_sets)
    identities=[r['source_command_max_error'] for r in module['rows'] if r['source_command_max_error'] is not None]
    assert len(identities)==336 and max(identities)<1e-10
    assert len({(r['source'],r['method'],key(r)) for r in module['rows']})==2688
    assert len({r['samples_file'] for r in module['rows']})==2688
    for r in module['rows']:
        assert (OUT/r['samples_file']).exists()
        for g in r['groups'].values():assert 0<=g['covered']<=g['accepted']<=g['valid']<=g['n']
    capture=load('baseline_capture.json');calibration=load('baseline_calibration.json')
    assert len(capture['rows'])==336 and all(r['physical_log_bitwise_equal'] for r in capture['rows'])
    assert len(calibration['positions'])==176
    assert all(abs(r['J_reproduction_error'])<1e-9 for r in diagnostic['cases'])
    return dict(passed=True,development_episodes=len(dev['rows']),configs_per_method=24,
                test_episodes=len(test['rows']),paired_cells_per_method=len(cell_sets[0]),
                source_logs=336,calibrations=176,module_replays=len(module['rows']),
                source_command_max_error=max(identities),physical_source_logs_bitwise_equal=True,
                synthetic_cases=len(diagnostic['rows']),frozen_config_case_checks=len(diagnostic['cases']),
                frozen_config_max_J_error=max(abs(r['J_reproduction_error']) for r in diagnostic['cases']))


def generate():
    dev=load('baseline_dev.json');test=load('baseline_test.json');checks=load('baseline_checks.json')
    rows=test['rows'];methods=list(test['chosen'])
    assert set(dev['chosen'])==set(METHODS) and len(dev['rows'])==len(METHODS)*24*8
    assert all(len([r for r in rows if r['method']==m])==168 for m in methods)
    assert all(len({key(r) for r in rows if r['method']==m})==168 for m in methods)
    assert checks.get('source_replay_command_error') is not None
    assert (OUT/'baseline_module.json').exists() and (OUT/'baseline_diagnostics.json').exists()
    assert len(load('baseline_module.json')['rows'])==len(ID_METHODS)*2*168
    validation=validate_primary(dev,test,load('baseline_module.json'),load('baseline_diagnostics.json'))
    refname='Huber-current' if 'Huber-current' in methods else 'Huber-D4'
    base={key(r):r for r in rows if r['method']==refname};tuned={key(r):r for r in rows if r['method']=='Huber-D4'}
    lines=['# 可靠接触辨识：有竞争力 baseline 的公平仿真比较','',
           '2026-11-06。状态：实现和运行结果待 Claude 验收；未提交、未推送。固定对象/主场景不变，所有测试使用本轮新种子。',
           '',(HERE/'baseline_protocol.md').read_text(encoding='utf-8').split('## A.',1)[1].join(['## A.','']),
           '\n### 单元检查实际结果\n',table(['检查','实际结果'],[[k,str(v)] for k,v in checks.items()]),
           '\n## B. 调参\n',
           '每方法严格 24 配置×8 开发回合。默认是本轮网格第 0 项，未把适配参数称为原文所有参数的默认值；完整 Method 第 0 项对应当前默认。公共 tracker 不变。以下选择仅依开发场景，测试之后不改参数。']
    tab=[]
    for method,s in dev['scores'].items():tab.append([method,fmt(s['default']),fmt(s['best']),s['index'],json.dumps(dev['chosen'][method],ensure_ascii=False)])
    lines.append(table(['方法','默认 L','选中 L','配置序号','选中配置'],tab))
    lines+=['Xue 的 Δ 敏感性（每个 Δ 内其他参数的开发最小 L；不可当测试结果）：\n']
    tab=[]
    for method in ('Xue-const','Xue-chi'):
        s=dev['scores'][method]
        for delta in (.005,.05,.5):
            vv=[v for cfg,v in zip(s['configs'],s['all']) if cfg['delta']==delta]
            tab.append([method,delta,fmt(min(vv)),fmt(np.mean(vv))])
    lines.append(table(['方法','Δ','最小开发 L','所有该 Δ 配置平均 L'],tab))
    lines+=['\n## C. 模块级结果\n']
    module_summary={}
    if (OUT/'baseline_module.json').exists():
        text,module_summary=module_tables(load('baseline_module.json'));lines.append(text)
    else:lines.append('模块回放尚未完成。')
    if (OUT/'baseline_diagnostics.json').exists():
        diagnostic=load('baseline_diagnostics.json')
        lines+=['\n### 独立人工刚度阶跃\n',
                '独立构造测力方程：3 s 时 2000→8000 N/m，探头平均位置距初始表面 4 mm，加 20.8 Hz、0.1 mm 探测；另给 0.3 Hz、1 mm 表面运动。固定 H/宿主相位和门控，仅回放估计器；不用该检查重新选参。保持 5 拍进入真实刚度的 1/1.2–1.2 倍范围后记恢复；— 表示截至 6 s 未恢复。误差仍为绝对对数误差；它只检验这个线性方程，不代替原对象闭环。Xue 在此获得精确初始接触点，避免把静止检查的失败归给标定误差。']
        tab=[]
        for r in diagnostic['rows']:
            tab.append([r['method'],'运动' if r['moving'] else '静止',r['pre_step_accurate'],fmt(r['before_error']),fmt(r['after_error']),
                        fmt(r['available_before']),fmt(r['available_after']),fmt(r['delay']['finite']),fmt(r['delay']['accepted']),
                        fmt(r['serial_estimator_p50_us'],0)+'/'+fmt(r['serial_estimator_p99_us'],0)])
        lines.append(table(['估计器','表面','阶跃前已入界','阶跃前误差','阶跃后末秒误差','前可用率','后可用率','有限入界 s','接受入界 s','单进程估计 p50/p99 μs'],tab))
        lines.append('入界时间保留原始统计；只有阶跃前已连续 5 拍入界的样例才能解释成恢复。Xue-const 静止样例的 0 s 是阶跃前约 3.4 倍的高估恰好接近阶跃后真值，不能称为即时跟踪；表中阶跃前误差和资格标志用于暴露这一点。')
        lines+=['\n### 冻结配置的代表性轨迹核查\n',
                '事后选取 Rossi 相对默认 Method 最占优的跨肋骨样例、默认 Method 最大 J 样例及两者配对，并取调后 Huber、DO-MPAC、Xue-const 的最大力样例。只说明尾部机制，不新增统计胜负或选参。每个重跑须复现原 J；原始物理日志和实际补偿后观测/内部状态保存在 `baseline_cases/*.npz`。k_good/真实增益比是在峰前 0.25 s 的描述统计，不等于已证明前馈误差造成该峰。',
                table(['方法','场景','J','峰值 N','峰值时刻 s','峰前 k_good/k','重跑 J 差'],
                      [[r['method'],str(r['cell']),fmt(r['J']),fmt(r['peak'],1),fmt(r['peak_time']),fmt(r['preceding_k_good_to_true_gain']),fmt(r['J_reproduction_error'],10)] for r in diagnostic.get('cases',[])])]
    lines+=['\n## D. 闭环结果与配对差\n',f'主配对参考为 **{refname}**（当前默认完整 Method）。另列相对本轮等预算调后的 Huber-D4，避免把我方重调收益归给新 baseline。每方法 n=168；RMSE 单位 N。数值失败的 L 仍按 J=1000 计入，不删除失败回合；RMSE/峰值只能对有日志回合计算，失败单列，不能把未定义 RMSE 的消失当作优势。事件 RMSE 的聚合是逐回合 RMSE 的中位/均值，不是把所有事件采样点混为一个独立样本。']
    tab=[];stats={};diffs={};groups={}
    for method in methods:
        rr=[r for r in rows if r['method']==method];s=stat(rr);d=paired(rr,base);dt=paired(rr,tuned)
        stats[method]=s;diffs[method]=dict(current=d,tuned=dt)
        tab.append([method,fmt(s['L']),fmt(s['rmse_med'])+'/'+fmt(s['rmse_mean']),fmt(s['event_med'])+'/'+fmt(s['event_mean']),s['over15'],s['over36'],s['long'],s['fail'],fmt(s['peak'],1),difftext(d),difftext(dt)])
    lines.append(table(['方法','平均 L','RMSE 中位/均值','事件 RMSE 中位/均值','>15','>36.5','长失联','数值失败','最大力','ΔL 对当前 [95% CI]','ΔL 对调后 [95% CI]'],tab))
    for label,predicate in [('六工况',lambda r:r['cond']),('动态模型',lambda r:r['dynamics']),('相位',lambda r:str(r['phase']))]:
        lines+=['\n### 按'+label+'\n'];tab=[]
        names=sorted({predicate(r) for r in rows})
        for method in methods:
            for name in names:
                rr=[r for r in rows if r['method']==method and predicate(r)==name];s=stat(rr);d=paired(rr,base)
                groups[f'{method}:{label}:{name}']=dict(stats=s,difference=d)
                tab.append([method,name,len(rr),fmt(s['L']),fmt(s['rmse_med']),fmt(s['event_med']),s['over15'],s['over36'],s['long'],s['fail'],difftext(d)])
        lines.append(table(['方法','分组','n','L','RMSE 中位','事件中位','>15','>36.5','长失联','失败','ΔL 对当前 [95% CI]'],tab))
    lines+=['\n### 六工况等权与调用时间\n'];tab=[]
    for method in methods:
        rr=[r for r in rows if r['method']==method];s=stats[method]
        eq=np.mean([stat([r for r in rr if r['cond']==c])['L'] for c in ('breath','shift_fast','flinch','cough','rib','combo')])
        stats[method]['equal_condition_L']=float(eq)
        tab.append([method,fmt(eq),fmt(s['calls50'],0)+'/'+fmt(s['calls99'],0),s['sat'],s['nonfinite']])
    lines.append(table(['方法','六工况等权 L','每回合调用 p50/p99 的中位 μs','饱和命令总数','非有限命令总数'],tab))
    lines+=['以上调用时间来自 15 进程并行竞争，p99 一列是回合内 p99 的跨回合中位。完整时延分布可从原始 rows 复算。不能和单核实时截止期直接比较。主表显示按 (工况,组织,种子 i) 聚类的配对 bootstrap 区间：总体 42 个物理场景簇，每簇四个动态/相位版本；分工况为 3 或 9 簇。普通 168 行配对与保持结构格样本数的分层 bootstrap 亦保存到 `baseline_summary.json` 供核查，不能忽略共同组织参数带来的相关性。']
    if (OUT/'baseline_diagnostics.json').exists():
        lines+=['\n### 单进程暖机全策略调用\n',
                '在同一份 Huber-D4 呼吸日志上，逐方法完整回放含 ForceSensing 的策略，舍弃前 125 次调用。以下是固定源日志计时，不是该控制器自己的闭环性能，也不是实时认证。',
                table(['方法','p50 μs','p99 μs'],[[r['method'],fmt(r['p50_us'],0),fmt(r['p99_us'],0)] for r in diagnostic['serial_full_calls']])]
    stage_rows=[]
    for filename,label in [('baseline_dev.json','开发'),('baseline_test.json','新种子闭环'),('baseline_capture.json','固定命令观测补录/逐位验证'),('baseline_module.json','双源模块回放')]:
        if (OUT/filename).exists():
            saved=load(filename);stage_rows.append([label,len(saved['rows']),fmt(saved['seconds']/60.,2)])
    lines+=['\n### 实际运行量与耗时\n',table(['阶段','回合/回放数','墙钟 min'],stage_rows),
            '各主阶段使用 15 个工作进程、BLAS 单线程；独立标定另 176 回合，串行代表性核查至多 7 个冻结配置回合。长跑前做了 11 方法小样计时。实际耗时受本机可用 CPU、进程启动与任务竞争影响，未把论文中毫秒级单步代价当作整组实验耗时。']
    lines.append('数据完整性校验已核对每方法 24 个唯一配置、每配置 8 个开发回合、所选配置确为原始开发 J 的最小平均 L、测试配置冻结、场景逐一配对、2688 个回放文件无重名、336 份物理日志逐位一致、全部自身源策略命令复现及 7 个尾部重跑 J 一致。机器可读结果见 `baseline_validation.json`。')
    lines+=['\n## E. 结论与机制解释\n']
    interpretation=HERE/'baseline_interpretation.md'
    lines.append(interpretation.read_text(encoding='utf-8') if interpretation.exists() else '等待全部结果后的机制分析。')
    lines+=['\n## F. 不确定或可能不公平之处\n',
            '1. 这些是明确适配的论文机制，不能宣称逐篇完整复现实验。尤其 Haninger 的动力学 EKF 和 Zhang 的已知表面 MHE 被适配到相同力/编码器信息；原论文的额外信息不能保留。',
            '2. Xue 有测量标定得到的初始接触位置，其他方法未使用该先验；标定时间另计。切线、探测增益及初始深度的差别影响 Xue，失败不能简单归因于死区本身。固定 σ 的覆盖率没有概率意义。',
            '3. 24 组是相同次数预算，不是相同搜索难度。Kalman、MHE 和 HC 的参数维度更高；边界最优不证明方法能力到此为止。Kalman 选到 q_k 下边界、q_d 上边界，Xue 选到 Δ 上边界；Rossi 控制器选到 KI 上边界，DO-MPAC 选到 ESO 带宽上边界和 Q/R 下边界。这些结果只裁决本次配置集；MHE 没有将物理不等式加入轨迹优化，输出才用公共 k 范围裁剪，不能声称覆盖约束 MHE 的全部能力。不同方法在结构/门控维度上的配置分配不同，内部基线搜索了更多门控尺度。不能在看测试后给任何一方补参。',
            '4. 同一 H/宿主隔离了估计器效应，也限定了结论；联合设计的完整架构可有不同优势。DO-MPAC 与 Xing 没有被塞入 Method。',
            '5. 实测 σ 覆盖率包含噪声、模型偏差及已有真值定义误差。P3 col21 的标称探测增益未包含探头惯性经机械臂柔顺产生的全部路径，既有研究指出真实编码器到力的增益因此比 col21 高约 4–9%，估计相对该列可表现为正偏；另有周期内插值误差。它是所有回放共同的评分边界，不能把小偏差都归为估计器不准。',
            '6. 初始场景与接触/肌肉活化均为仿真构造，未验证人体；bootstrap 是此场景分布下的经验区间。事件段相互重叠，窗口非独立，模块表作描述统计，不对数十万窗口计算虚假的极窄置信区间。',
            '7. 168 格子对快速变化工况重复更多，pooled 分数并非六工况等权；两者都报告。多分组区间未作多重校正；单个显著区间不是通用胜负。',
            '8. 固定增益 MPC 在该伺服上的实际闭环极点不等于名义 A−BK；论文内部的离散/符号歧义已按实际采用公式声明。时延没有被计时数字认证。',
            '\n## G. 现有文件改动\n',
            '现有仿真、Method、相位检测、估计器及场景文件均未修改；新增局部适配器用 phase-only D⁴ 矩阵接入已有 update_phase。新增源码见 `id_baselines.py`、`ctrl_baselines.py`、`baseline_compare.py`、`baseline_checks.py`、`baseline_diagnostics.py`、`baseline_report.py`；原文卡片与协议见 `baseline_protocol.md`，结果解释见 `baseline_interpretation.md`。`project_records/decision_log.md` 和 `project_records/agent_runs/2026-11-06_codex_baselines/` 仅补充项目约定的审阅记录；未运行同步脚本，未提交或推送。冻结快照、原文及既有提取文本保持不变。',
            '\n### 复现命令\n','```powershell',
            'python -X utf8 experiments/rossi_failure_study/baseline_checks.py',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage calibrate',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage dev',
            'python -X utf8 experiments/rossi_failure_study/baseline_checks.py --selected',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage test',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage capture',
            'python -X utf8 experiments/rossi_failure_study/baseline_checks.py --replay',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage replay',
            'python -X utf8 experiments/rossi_failure_study/baseline_diagnostics.py',
            'python -X utf8 experiments/rossi_failure_study/baseline_compare.py --stage report','```']
    # Stage JSON is already strict; convert numpy scalars for the analysis JSON.
    summary=dict(stats=stats,differences=diffs,groups=groups,module=module_summary,reference=refname)
    from baseline_compare import dump
    dump('baseline_summary.json',summary)
    dump('baseline_validation.json',validation)
    REPORT.parent.mkdir(parents=True,exist_ok=True);REPORT.write_text('\n\n'.join(lines),encoding='utf-8')
    print(REPORT)


if __name__=='__main__':generate()
