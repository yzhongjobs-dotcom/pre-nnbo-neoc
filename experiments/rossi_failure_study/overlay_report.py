"""Recompute and report only the three additive scene comparisons."""
import csv
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np
import baseline_compare as B
import baseline_report as F
import overlay_compare as O

DISCUSSION = '''## 结果分析与剩余问题

**总体：当前Huber-D⁴平均控制损失最低，但不是处处最好。** 84个实例的平均L为0.181，三场景等权为0.171。相对D⁴、Rossi-ID、Kalman、MHE的总体配对区间均跨零，因此本批不能证明对这些主要对手有稳定优势。对P5、两种Xue和三个完整控制器的总体ΔL区间不跨零，仍须连同小样本和未校正的比较口径解释。

**呼吸＋跨肋骨：** Xing的平均L=0.114，低于Huber-D⁴的0.152，其评分段最大力也仅9.9 N，对比Huber的35.5 N。这是实质保留的反例，不能写成Huber在跨肋骨一定最好。两者ΔL=-0.039，95%区间[-0.141,0.024]，尚不支持稳定排名。Huber的RMSE中位/均值为0.117/0.179 N，均值受一个28°实例的尖峰影响。

**呼吸＋突然推拉：** D⁴=0.119、Rossi-ID=0.132、Huber=0.137、Kalman/MHE≈0.225。Huber与Rossi-ID基本接近，内部D⁴略好（ΔL=-0.018，区间[-0.032,-0.001]，仅3个物理簇）。Huber平均RMSE=0.187 N，fast事件平均RMSE=0.545 N、最大力10.3 N。完整Rossi、DO-MPAC、Xing在fast事件中平均RMSE约1.83、5.07、3.65 N，明显高于此处的估计器替换组；这是本任务中已有完整实现的表现，不是原论文所有任务的裁决。

**呼吸＋跨肋骨＋突然推拉：** Huber的L=0.223为最低点估计，Kalman=0.245、D⁴=0.252、MHE=0.259、Rossi-ID=0.284；这些主要对照的配对区间都跨零。Huber平均RMSE=0.259 N，事件平均RMSE=0.390 N，但最大力仍30.1 N。Rossi-ID出现85.9 N峰值；完整Rossi、DO-MPAC、Xing的平均L较高，但本批最大力均低于Huber，平均精度与尾部不能混为一谈。

**刚度精度并非Huber全面领先。** 双源回放的有限误差中位数，Huber为0.179/0.204/0.185，Kalman/MHE为0.136/0.287/0.127（依次为三个场景）。含肋骨时Kalman/MHE中位误差更低；但Huber的P90为0.382/0.494/0.400，低于Kalman/MHE的0.471/0.581/0.538，可用率也高于它们。Rossi-ID中位误差约0.183/0.195/0.193，可用率约97–99%，覆盖率却只有41–54%。D⁴的点误差也略低于Huber，但可用率仅13–29%。这些差别说明需要共同评价误差分布、门控和控制接口，不能只挑接受窗口误差或把控制收益全部归给刚度点估计。

**Xue仍不适配本批持续运动。** 原接法的有限误差中位数约0.601/0.824/0.768；加χ后约0.427/0.536/0.514，但可用率仅0.06–0.13%，其较低控制损失不表示在线辨识已恢复正常。未知运动造成的参考压入偏差、HC与当前对象失配、初始化误差及切线/探测频率增益目标不同都仍是原协议的解释边界，不能只归因于死区，也不能把常数σ的极低覆盖率说成原文概率保证失效。

**仍未解决：** Huber本批0次长失联、0次>36.5 N，但有2回合>15 N，最高35.5 N；两个大峰均在28°层，尚无因果分解。其整体覆盖率约82–84%，三因素event段约77%，呼吸＋跨肋骨changing段约68%。平均控制较好不等于瞬态或不确定性已可靠，刚度的典型绝对对数误差约0.18–0.20也不等于精确真值。本批只补证据，不据此改参数、改方法或另开实验。
'''


def validate(test,module):
    plan = O.ORIGINAL_LOAD('overlay_plan.json')
    assert test['chosen'] == plan['chosen'] == O.chosen()
    assert all(O.sha(B.HERE/p) == h for p,h in plan['protected_hashes'].items())
    assert O.sha(B.OUT/'baseline_test_r2.json') == plan['config_source_sha256']
    expected = set(O.cells()); errors = []; crossing = {}
    for m,cfg in test['chosen'].items():
        rows = [r for r in test['rows'] if r['method'] == m]
        assert len(rows) == len(expected) == len({B.key(r) for r in rows})
        assert {B.key(r) for r in rows} == expected
        for r in rows:
            assert r['cfg'] == cfg and r['nonfinite_commands'] == 0
            assert not r['failed'], (m,B.key(r))
            with np.load(B.OUT/r['log_file']) as z:
                cols = z['cols'].tolist(); log = z['log']
                get = lambda c: log[:,cols.index(c)]
                t, force = get(0),get(18)
                assert len(t) == 1250 and np.allclose(np.diff(t),.008,rtol=0,atol=1e-12)
                warm = t >= 1.; f = force[warm]; err = f-8.
                J = min(1000.,float(np.mean((err/.5)**2)+20*np.mean(f<.5)+np.mean(np.maximum(f-11.,0.)**2)))
                diffs = [abs(J-r['J']),abs(np.sqrt(np.mean(err**2))-r['rmse']),abs(f.max()-r['peak'])]
                normal=get(3)[warm]-8.
                diffs += [abs(np.sqrt(np.mean(normal**2))-r['normal_rmse']),
                          abs(np.mean(err)-r['mean_err']),abs(np.mean(normal)-r['mean_normal_err']),
                          abs(np.mean(f>15.)*(t[-1]-1.)-r['over_15']),
                          abs(np.mean(f>11.)*(t[-1]-1.)-r['over_band']),
                          abs(np.mean(f<.5)*(t[-1]-1.)-r['loss_time'])]
                assert r['over_36']==bool(f.max()>36.5)
                longest=cur=0
                for loss in f<.5:
                    cur = cur+1 if loss else 0;longest=max(longest,cur)
                assert r['long_loss'] == (longest*.008>.5)
                sc = B.scene(B.key(r),O.BASE,0.)
                event = np.zeros(len(t),bool)
                for _,a,b in sc['events']:
                    event |= (t>=a)&(t<=b)
                if np.any(event & warm):
                    diffs.append(abs(np.sqrt(np.mean((force[event & warm]-8.)**2))-r['event_rmse']))
                else:
                    assert r['event_rmse'] is None
                errors.extend(diffs)
                # Measured scan-axis encoder path, not merely xcmd. The plant
                # additionally retains its original tangential contact state.
                x=get(1)
                crossed=sum(x.min()<=center<=x.max() for center,*_ in sc['ribs'])
                assert bool(crossed) == (r['cond']!='breath_fast')
                crossing[m+':'+':'.join(map(str,B.key(r)))]=int(crossed)
    assert max(errors)<1e-10
    sources = O.ORIGINAL_LOAD('overlay_capture.json')['rows']
    assert len(sources)==168 and all(r['physical_log_bitwise_equal'] for r in sources)
    assert len(module['rows'])==1344
    assert len({(r['source'],r['method'],B.key(r)) for r in module['rows']})==1344
    assert len({r['samples_file'] for r in module['rows']})==1344
    identities=[];sample_errors=[]
    for r in module['rows']:
        for g in r['groups'].values():
            assert 0<=g['covered']<=g['accepted']<=g['valid']<=g['n']
        if r['source_command_max_error'] is not None:
            identities.append(r['source_command_max_error'])
        with np.load(B.OUT/r['samples_file']) as z:
            assert len(z['t'])==len(z['estimate'])==len(z['truth'])==len(z['accepted'])==1250
            sc=B.scene(B.key(r),O.BASE,0.);t=z['t'];k=z['truth'];kh=z['estimate'];ok=z['accepted'];sig=z['sigma']
            event=np.zeros(len(t),bool)
            for e in sc['event_info']:
                event|=(t>=e['start'])&(t<=e['end']+.3)
            changing=np.convolve(abs(np.gradient(np.log(np.maximum(k,1.)),t))>2.,np.ones(30),'same')>0
            for kind in ('quiet','event','changing'):
                mask=(t>=1.)&(k>200.)
                mask &= event if kind=='event' else (~event & (changing if kind=='changing' else ~changing))
                finite=mask&np.isfinite(kh)&(kh>0)&np.isfinite(sig)
                g=r['groups'][kind]
                assert g['n']==int(mask.sum())
                assert g['valid']==int(finite.sum())
                ae=z[kind+'_accepted_errors']
                assert len(ae)==g['accepted'] and len(z[kind+'_error_samples'])==g['valid']
                assert int((mask&ok).sum())==g['accepted']
                expected_errors=np.abs(np.log(kh[mask&ok]/k[mask&ok]))
                assert np.allclose(ae,expected_errors,rtol=0,atol=1e-12)
                assert np.allclose(z[kind+'_error_samples'],abs(np.log(kh[finite]/k[finite])),rtol=0,atol=1e-12)
                assert g['covered']==int(np.sum(expected_errors<=2*sig[mask&ok]/kh[mask&ok]))
                sample_errors.extend(abs(ae-expected_errors).tolist())
    assert len(identities)==168 and max(identities)<1e-10
    return dict(passed=True,episodes=len(test['rows']),paired_cells=84,replays=1344,
                physical_replays_bitwise_equal=168,source_policy_replays=168,
                max_command_error=max(identities),max_metric_error=max(errors),
                max_accepted_sample_error=max(sample_errors,default=0.),actual_rib_crossings=crossing,
                frozen_core_unchanged=True,old_scene_count_bitwise_equal=plan['old_scenes_bitwise_equal'],
                python=platform.python_version(),numpy=np.__version__)


def generate():
    O.setup()
    test=O.ORIGINAL_LOAD('overlay_test.json');module=O.ORIGINAL_LOAD('overlay_module.json')
    qa=validate(test,module);B.dump('overlay_checks.json',qa)
    rows=test['rows'];methods=list(test['chosen']);base={B.key(r):r for r in rows if r['method']=='Huber-D4'}
    summary=dict(overall={},by_scene={},module={})
    lines=['# 三个叠加工况：原协议增量比较（2026-11-07）','',
           '## 范围与复现','',
           '仅补三个场景；不改物理对象、方法或参数，不做开发选参。`overlay_compare.py` 是原 `baseline_compare.py` 的场景/种子适配入口，仍调用原 execute、标定、capture、replay_job、M.metric 和 bootstrap。论文机制的适配边界见 [原协议](baseline_protocol.md)，原比较见 [原分析](baseline_interpretation.md)。已核对 `structured_tvk_audit.py` 的同轨迹、真值仅评分原则；本次维持 baseline 协议的8 N/10 s，不照搬该专项的4 N/6 s，不引入其候选估计模型。',
           '',
           '- `breath_rib`：原 breath 的连续呼吸＋原肋骨组织场；原 rib 组织层，无额外 mixed 夹杂。',
           '- `breath_fast`：原 shift_fast 的直接别名，物理参数与波形逐项一致。原 shift_fast 本来就含呼吸。',
           '- `breath_rib_fast`：原 shift_fast＋原肋骨场，只有一个 fast 推拉事件，未混入 slow/flinch/cough，也无肌肉活化。',
           '- 三者保留原有心搏小扰动、表面起伏、非线性接触、记忆、摩擦、机器人/人体柔顺、探头惯性、伺服、测量噪声和失配。没有用 combo 替代。',
           '- 原 CONDITIONS 九位置平衡和旧种子不变；新别名借用 breath/shift_fast 的随机化层。两类肋骨场景每格9种子，快推拉每格3种子，均为 neural_white/viscous × 14°/28°，共36/12/36=84配对场景。物体base=674500，测量base=957000；不据这批结果调参。',
           '- FD=8 N、10 s、1 kHz运动驱动、0.5 ms积分、8 ms控制周期及全部评分设置沿用原协议；11种方法共924回合。Xue另用原静止接触标定（同一场景副本、2 s、0.5 N阈值），并非真值接触位置。',
           '- Huber-D4 表示已采纳的当前默认（原 Huber-current），没有采用第一轮重调后被拒绝的 Huber-D4配置；其余全部复制 baseline_test_r2.json 的选中配置。完整参数及输入哈希在 results/overlay_plan.json。',
           '- 运行：`python -X utf8 experiments/rossi_failure_study/overlay_compare.py --stage STAGE`，依次 check、pilot、test、capture、replay、report；结果为 overlay_*.json。逐回合压缩日志在原 results/baseline_logs/，模块样本在原 results/baseline_replay/（均按仓库约定忽略、可重建），本次所有场景/方法均保存日志。',
           '', '## 恒力控制','',
           '评分t≥1 s，L=ln(1+J)，J、RMSE、事件RMSE、超限和长失联完全沿用原M.metric。越低越好；均值L与RMSE不是同一指标。配对ΔL=对手−当前Huber-D⁴，95%分位区间按(工况,组织,i)物理场景簇bootstrap 5000次；总21簇，分场景9/3/9簇。四个动态/相位版本不当作独立物理重复；小簇区间不稳定，多个比较未校正，不能据此作人体或尾部安全推断。']
    alltab=[];csvrows=[]
    for c in (None,*O.SCENES):
        tab=[]
        for m in methods:
            rr=[r for r in rows if r['method']==m and (c is None or r['cond']==c)]
            s=F.stat(rr);d=F.paired(rr,base)
            item=dict(stats=s,delta=d)
            if c is None:summary['overall'][m]=item
            else:summary['by_scene'].setdefault(c,{})[m]=item
            tab.append([m,F.fmt(s['L']),F.fmt(s['rmse_med'])+'/'+F.fmt(s['rmse_mean']),
                        F.fmt(s['event_med'])+'/'+F.fmt(s['event_mean']),F.fmt(s['peak'],1),
                        s['over15'],s['over36'],s['long'],s['fail'],F.difftext(d)])
            csvrows.append(dict(scene=c or 'all',method=m,**s,delta_L=d['mean'],ci_low=d['ci'][0],ci_high=d['ci'][1]))
        lines+=['','### '+('全部（36/12/36加权）' if c is None else O.LABELS[c]),'',
                F.table(['方法','L','RMSE中位/均值 N','事件RMSE中位/均值 N','最大力 N','>15 N回合','>36.5 N回合','长失联','数值失败','ΔL [95% CI]'],tab)]
    equal={m:float(np.mean([summary['by_scene'][c][m]['stats']['L'] for c in O.SCENES])) for m in methods}
    summary['scene_equal_L']=equal
    lines+=['','三场景等权L：'+ '；'.join(m+' '+F.fmt(equal[m]) for m in methods)+'。',
            '', '呼吸＋跨肋骨没有sc.events时间事件，故原协议事件RMSE为缺失，不能填0；它的快速刚度过渡由模块changing统计覆盖。另两类事件RMSE仅在原fast事件起止内，模块event则延至end+0.3 s。峰值为评分段采样峰值，并非积分子步最大值。',
            '', '## 刚度估计（相同观测双源回放）','',
            '使用当前Huber-D⁴和D⁴各84份固定源日志，8个原辨识器共1344次回放。每个源经相同命令物理重放捕获真实回调观测，并验证完成日志逐位相同。估计器只见测力、编码器和已发命令；真实k只在回放结束后评分。点误差=|ln(khat/|K20.8Hz|)|；接受误差统计只含接受窗口，可能存在筛选偏差，需连同可用率读。这里评价的是探测频率有效刚度模，不是静态切线刚度。',
            '', '可用率分母含全部t≥1 s且真值k>200的采样时刻，包括未调用/无效估计。覆盖率沿用原|ln(khat/k)|≤2σ/khat，是近似校准诊断。Xue固定σ=1 N/m没有不确定性模型，其覆盖率不宜解释成原论文概率保证的失败。两个源合并的样本不是独立重复，以下为描述统计。']
    moduletab=[]
    for c in O.SCENES:
        for m in B.ID_METHODS:
            rr=[r for r in module['rows'] if r['method']==m and r['cond']==c]
            for kind in ('all','quiet','event','changing'):
                restricted=rr if kind=='all' else [dict(r,groups={kind:r['groups'][kind]}) for r in rr]
                s=F.aggregate_group(restricted)
                summary['module'][f'{c}:{m}:{kind}']=s
                # Retain source-specific summaries so feedback distribution
                # cannot be silently confused with estimator superiority.
                for source in ('Huber-D4','D4'):
                    summary['module'][f'{source}:{c}:{m}:{kind}']=F.aggregate_group([r for r in restricted if r['source']==source])
                moduletab.append(dict(scene=c,method=m,group=kind,**s))
        tab=[]
        for m in B.ID_METHODS:
            s=summary['module'][f'{c}:{m}:all']
            tab.append([m,F.fmt(s['valid']),F.fmt(s['available']),F.fmt(s['error_med'])+'/'+F.fmt(s['error_p90']),
                        F.fmt(s['accepted_med']),F.fmt(s['coverage']),str(s['recovered'])+'/'+str(s['recovery_n']),F.fmt(s['recovery_med'])])
        lines+=['','### '+O.LABELS[c],'',F.table(['估计器','有限率','可用率','有限误差中位/P90','接受误差中位','2σ覆盖率','自然过渡入界/总数','入界延迟中位 s'],tab)]
    segments=[]
    for c in O.SCENES:
        for m in ('Huber-D4','Rossi-ID','Kalman','MHE'):
            ev=summary['module'][f'{c}:{m}:event'];ch=summary['module'][f'{c}:{m}:changing']
            segments.append([O.LABELS[c],m,F.fmt(ev['available']),F.fmt(ev['coverage']),F.fmt(ch['available']),F.fmt(ch['coverage'])])
    lines+=['','### 事件与刚度变化段的可靠性','',
            F.table(['场景','估计器','event可用率','event覆盖率','changing可用率','changing覆盖率'],segments),
            '', '完整quiet/event/changing、分源统计见overlay_summary.json与overlay_identification.csv。自然过渡入界沿用连续5拍、对数20%界、1.5 s等待和右删失；延迟只对已入界者取中位，未要求过渡前准确，因此不等同严格阶跃恢复。CSV的恢复统计属于整条轨迹，分段行不应再次累加恢复数。Rossi-control的刚度属于自身未作H校正的回归目标，DO-MPAC为固定名义刚度，Xing不估刚度，均不混入上述估计器表。',
            '', '## 图表与读法','',
            '- [配对恒力控制差值](figures/overlay/paired_control.png)：全部对手的ΔL及原协议物理场景簇95%区间；正值表示当前Huber-D⁴更低。',
            '- [同观测刚度误差](figures/overlay/identification.png)：全部8种辨识器的有限误差中位/P90与接受误差中位，和表中的可用率同时读。',
            '- [固定实例恒力轨迹](figures/overlay/example_force.png)、[固定实例刚度轨迹](figures/overlay/example_stiffness.png)：预定neural_white、14°、i=0，展示当前Huber-D⁴、Rossi-ID、MHE和完整Rossi；完整Rossi的刚度目标不同，刚度图仅显示前三者。阴影为fast运动区间；刚度图全部使用同一Huber-D⁴源观测。固定实例不能代替全体统计，主要表和误差图保留全部方法。',
            '- 同目录另存可编辑SVG/PDF、最终轴几何JSON；表格源数据为results/overlay_control.csv和overlay_identification.csv，逐回合分数及参数在overlay_test.json。',
            '', '三个场景之间并未额外做因子配对或选择事件与肋骨边缘同时到达的时机，因此不从场景均值差推断三因素交互的因果量。范围是原随机化与物理假设中的增量仿真，不是人体实验。',
            '', DISCUSSION, '', '## 完成前核验','',
            f"全部{qa['episodes']}回合独立重算J、RMSE、事件RMSE和峰值，最大差{qa['max_metric_error']:.3g}；长失联和采样轴一致。按各方法扫描轴编码器位置逐回合确认轨迹覆盖肋骨中心；接触点还受原对象切向柔顺影响，编码器不冒充接触点真值。168个源物理重放逐位相同、168个源策略自身回放命令一致，最大差{qa['max_command_error']:.3g} m。1344份模块结果的键、时序、采样数、接受误差和分组均核对通过。36个旧场景参数/波形逐项一致；物理对象、方法与原选参文件哈希未变。",'',
            '计算时间见overlay_test.json和逐回合seconds/call_p50/call_p99；15进程、BLAS单线程，原型并发计时不作实时认证。']
    B.dump('overlay_summary.json',summary)
    for name,data in [('overlay_control.csv',csvrows),('overlay_identification.csv',moduletab)]:
        with (B.OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    (B.HERE/'overlay_interpretation.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Validated metrics, scene crossings, causal replay and frozen implementations.',flush=True)
    from overlay_plots import plot
    plot(summary,test,module)


if __name__=='__main__':
    generate()
