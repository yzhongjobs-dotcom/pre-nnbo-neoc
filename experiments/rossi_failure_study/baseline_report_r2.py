"""Round 2/3 report from retained primary results; no selection or simulation."""
import hashlib
import json
import math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import baseline_compare as B
import baseline_round2 as R
import baseline_report as F


def validate(dev,test,module,mechanism):
    original=B.load('baseline_dev.json');counts={}
    for m,s in dev['scores'].items():
        cfgs=s['configs'];assert len(cfgs)==len({R.canon(m,c) for c in cfgs})
        vals=[]
        for cfg in cfgs:
            rr=[r for r in dev['rows'] if r['method']==m and r['cfg']==cfg]
            assert len(rr)==len({B.key(r) for r in rr})==8,(m,cfg,len(rr))
            vals.append(B.L(rr))
        assert np.allclose(vals,s['all'],atol=1e-12,rtol=0)
        idx=int(np.argmin(vals));assert idx==s['index'] and dev['chosen'][m]==cfgs[idx]
        counts[m]=len(cfgs)
        if m not in R.METHODS:assert s==original['scores'][m]
    for marker,cfgs in dev['added_configs'].items():assert len(cfgs)<=24
    old=B.load('baseline_test.json');cells=set(B.cells())
    for m,cfg in test['chosen'].items():
        rr=[r for r in test['rows'] if r['method']==m]
        assert len(rr)==len({B.key(r) for r in rr})==168 and {B.key(r) for r in rr}==cells
        assert all(r['cfg']==cfg for r in rr)
        if m not in R.changed():assert rr==[r for r in old['rows'] if r['method']==m]
    assert len(module['rows'])==2688
    assert len({(r['source'],r['method'],B.key(r)) for r in module['rows']})==2688
    assert len({r['samples_file'] for r in module['rows']})==2688
    for r in module['rows']:
        assert (B.OUT/r['samples_file']).exists()
        for g in r['groups'].values():assert 0<=g['covered']<=g['accepted']<=g['valid']<=g['n']
    oldmanifest=B.load('baseline_manifest.json');original_results={}
    for p,h in oldmanifest['result_hashes'].items():
        original_results[p]=hashlib.sha256((B.OUT/p).read_bytes()).hexdigest()==h
    assert all(original_results.values()),original_results
    manifest=B.load('baseline_manifest_r2.json')
    assert hashlib.sha256((B.OUT/'baseline_plan_r2.json').read_bytes()).hexdigest()==manifest['rules_sha256']
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in manifest['protected_core_hashes'].items())
    untouched=['project_records/decision_log.md']+[p for p in manifest['initial_hashes'] if p.endswith('/ctrl_baselines.py')]
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==manifest['initial_hashes'][p] for p in untouched)
    assert mechanism['kalman_nis_max_J_error']<1e-9
    # Verify original retained binary primary evidence, not just JSON indexes.
    trees={}
    for folder,expected in oldmanifest['data_groups'].items():
        # Round 1 used Windows/PowerShell's case-insensitive filename sort.
        digest=hashlib.sha256();paths=sorted((B.OUT/folder).glob('*'),key=lambda p:p.name.casefold())
        paths=[p for p in paths if p.is_file()]
        for p in paths:
            digest.update(p.name.encode('utf-8')+b'\0'+hashlib.sha256(p.read_bytes()).digest()+b'\n')
        trees[folder]=dict(files=len(paths),sha256_tree=digest.hexdigest())
        assert len(paths)==expected['files'] and digest.hexdigest()==expected['sha256_tree'],(folder,trees[folder])
    plan=B.load('baseline_plan_r2.json');assert plan['rules']==R.RULES
    return dict(passed=True,configs_per_method=counts,development_episodes=len(dev['rows']),
                added_per_round={k:len(v) for k,v in dev['added_configs'].items()},
                paired_cells_per_method=168,test_episodes=len(test['rows']),module_replays=len(module['rows']),
                rerun_methods=R.changed(),original_primary_results_unchanged=original_results,
                original_binary_evidence=trees,protected_core_unchanged=True,rules_unchanged=True,
                decision_log_unchanged=True,mechanism_episodes=len(mechanism['rows']),
                kalman_nis_max_J_error=mechanism['kalman_nis_max_J_error'])


def generate():
    dev=B.load('baseline_dev_r2.json');test=B.load('baseline_test_r2.json');module=B.load('baseline_module_r2.json')
    diag=B.load('baseline_diagnostics_r2.json');mechanism=B.load('baseline_mechanism_r2.json')
    validation=validate(dev,test,module,mechanism);rows=test['rows'];methods=list(test['chosen']);old=B.load('baseline_dev.json')
    base={B.key(r):r for r in rows if r['method']=='Huber-current'}
    tuned={B.key(r):r for r in rows if r['method']=='Huber-D4'}
    lines=['# 可靠接触辨识：baseline 第二轮边界扩展与更新比较',
           '2026-11-06。实现、运行和报告由 Codex 完成，待 Claude 验收；未提交、未推送、未改 decision_log。第一轮原始结果保留；新文件以 `_r2` 或独立轮次计划命名。',R.RULES]
    tab=[]
    for m,s in dev['scores'].items():
        tab.append([m,24,len(dev['added_configs'].get(m+':r2',[])),len(dev['added_configs'].get(m+':r3',[])),F.fmt(old['scores'][m]['best'],6),F.fmt(s['best'],6),json.dumps(dev['chosen'][m],ensure_ascii=False)])
    lines+=['### 新的选中配置与开发分数\n',F.table(['方法','第一轮数','第二轮新增','第三轮新增','第一轮最优 L','全部开发最优 L','最终选中配置'],tab)]
    lines+=['### 第三轮后残留边界\n']
    remaining={m:R.boundary_axes(m,dev['chosen'][m],dev['scores'][m]['configs']) for m in R.METHODS}
    lines.append('下表只描述已搜索正值数轴的边缘，不等于下一轮运行计划；本轮停止扩展。门控与结构菜单的边界另属固定离散选择。')
    lines.append(F.table(['方法','仍在边缘的连续轴'],[[m,', '.join(a) or '无'] for m,a in remaining.items()]))
    previous_summary=B.load('baseline_summary.json');tab=[]
    for m in R.METHODS:
        rr=[r for r in rows if r['method']==m];d=F.paired(rr,base)
        verdict='均值更低' if d['mean']<0 else '均值更高'
        verdict+=('，区间跨零' if d['ci'][0]<=0<=d['ci'][1] else '，区间不跨零')
        tab.append([m,F.fmt(previous_summary['stats'][m]['L']),F.fmt(B.L(rr)),F.difftext(d),verdict])
    lines+=['### 哪些测试结论改变\n',
            F.table(['外部方法','第一轮测试L','更新测试L','更新ΔL对当前[95%CI]','更新后判断'],tab),
            '具体工况与机制结论见E。第一轮“MHE开发选择q=0”的描述仅保留为旧网格事实，其机制解释撤回；是否仍选q=0直接按最终配置表判断。内部Method不增加预算，测试结果完全沿用第一轮。']
    lines+=['### 每个方法实际新增配置\n','以下完整清单按运行前冻结顺序列出；每份原始计划保留 frozen_utc，可与开发 JSON 逐条对应。']
    for m in R.METHODS:
        for round_number in (2,3):
            p=B.OUT/f'baseline_plan_{m}_r{round_number}.json'
            if p.exists():
                entry=json.loads(p.read_text(encoding='utf-8'))
                lines+=['#### '+m+' / 第'+str(round_number)+'轮\n',
                        '冻结时间：'+entry['frozen_utc']+'；新增 '+str(len(entry['configs']))+' 项。',
                        '```json\n'+json.dumps(entry['configs'],ensure_ascii=False,indent=1)+'\n```']
    original=R.REPORT.with_name('baseline_report.md').read_text(encoding='utf-8')
    section_a='## A.'+original.split('## A.',1)[1].split('## B.',1)[0]
    lines+=['## A. 第一轮原文与适配核对的沿用范围\n',
            '以下引用第一轮原文卡片和单元检查；其中“24组”“一阶MHE”“N=16/24”等描述的是第一轮。第二轮以开头冻结规则和下一段的明确增量为准，不将旧预算描述当作当前预算。',
            '二阶MHE新增状态 [k/1000,d,d_dot]，过程 A[1,2]=Tc，Q_d=q_d²[[Tc³/3,Tc²/2],[Tc²/2,Tc]]，与Kalman共用精确离散模型。窗口 Hessian 用块带状Cholesky解，左端预测边际只吸收窗前数据；χ为平滑轨迹的测量残差平方和/[噪声方差·max(N−状态维数,1)]，不含到达/过程代价。对一阶仍保留原max(N−2,1)，二阶为max(N−3,1)，这是统计接口近似，不能声称严格χ²分布。',
            section_a.replace('## A. 原文卡片与单元检查','### 第一轮卡片（历史设置）',1)]
    checks=B.load('baseline_checks_r2.json');cr=checks['second_order_finite_q']
    lines+=['### 第二轮单元检查\n',
            f"27组有限正q、N=16/24/48的二阶MHE/Kalman端点检查全部通过。最大状态绝对误差 {max(r['endpoint_abs_error'] for r in cr):.3g}；k方差最大相对误差 {max(r['variance_relative_error'] for r in cr):.3g}。同时保留原一阶有限q、q=0常数/RLS、原控制公式、真值遮蔽检查。检查脚本不写第一轮JSON。"]
    lines+=['## B. 选参及复用\n',
            '开发场景与第一轮完全一致，所有最终选择都由开发原始J重算最小平均L。测试不参加选择。每方法实际新增数见开头表；Huber-current、调后Huber-D4、D4、P5全部不重调。',
            '重新闭环的方法：'+', '.join(test['rerun_methods'])+'。其余方法逐条沿用第一轮结果，验证脚本核对完全一致。',
            '所有条件包含第一轮测试格子的同样168个新对象/传感种子；“新”相对于开发及更早研究，不是第二轮另换测试集。第二轮为对已知调参不足的协议修正，不能再把重复使用的测试集称为完全未查看的最终盲测集。']
    lines+=['## C. 模块级结果\n']
    moduletext,modulesummary=F.module_tables(module);lines.append(moduletext)
    lines+=['### 独立人工刚度阶跃\n',
            '与第一轮同一线性方程和随机种子，3s时2000→8000 N/m；静止/0.3Hz、1mm表面运动。当前表仅更新配置改变的估计器，其余沿用第一轮。入界需连续5拍满足|log(khat/k)|≤log(1.2)，未入界右删失。只有阶跃前已入界者才能称为恢复；不能把静止Xue的初始偏差碰巧接近新真值当作即时恢复。',
            F.table(['估计器','表面','阶跃前已入界','前误差','后末秒误差','前可用率','后可用率','有限入界s','接受入界s','估计p50/p99 μs'],[
                [r['method'],'运动' if r['moving'] else '静止',r['pre_step_accurate'],F.fmt(r['before_error']),F.fmt(r['after_error']),F.fmt(r['available_before']),F.fmt(r['available_after']),F.fmt(r['delay']['finite']),F.fmt(r['delay']['accepted']),F.fmt(r['serial_estimator_p50_us'],0)+'/'+F.fmt(r['serial_estimator_p99_us'],0)] for r in diag['rows']])]
    lines+=['## D. 闭环结果与配对差\n',
            '主参考为当前默认Huber-current，另列第一轮调后Huber-D4。每方法168回合；L=ln(1+J)，RMSE单位N。按第一轮warm≥1s评分，失败按J=1000保留；峰值/失联口径不变。事件RMSE是逐回合中位/均值。配对区间均为5000次物理场景簇bootstrap，总体42簇，每簇四个动态/相位版本；分工况3或9簇。']
    stats={};diffs={};groups={};tab=[]
    for m in methods:
        rr=[r for r in rows if r['method']==m];s=F.stat(rr);d=F.paired(rr,base);dt=F.paired(rr,tuned)
        stats[m]=s;diffs[m]=dict(current=d,tuned=dt)
        tab.append([m,F.fmt(s['L']),F.fmt(s['rmse_med'])+'/'+F.fmt(s['rmse_mean']),F.fmt(s['event_med'])+'/'+F.fmt(s['event_mean']),s['over15'],s['over36'],s['long'],s['fail'],F.fmt(s['peak'],1),F.difftext(d),F.difftext(dt)])
    lines.append(F.table(['方法','平均L','RMSE中位/均值','事件RMSE中位/均值','>15','>36.5','长失联','数值失败','最大力','ΔL对当前[95%CI]','ΔL对调后[95%CI]'],tab))
    for label,predicate in [('六工况',lambda r:r['cond']),('动态模型',lambda r:r['dynamics']),('相位',lambda r:str(r['phase']))]:
        lines+=['### 按'+label+'\n'];tab=[]
        for m in methods:
            for name in sorted({predicate(r) for r in rows}):
                rr=[r for r in rows if r['method']==m and predicate(r)==name];s=F.stat(rr);d=F.paired(rr,base)
                groups[f'{m}:{label}:{name}']=dict(stats=s,difference=d)
                tab.append([m,name,len(rr),F.fmt(s['L']),F.fmt(s['rmse_med']),F.fmt(s['event_med']),s['over15'],s['over36'],s['long'],s['fail'],F.difftext(d)])
        lines.append(F.table(['方法','分组','n','L','RMSE中位','事件中位','>15','>36.5','长失联','失败','ΔL对当前[95%CI]'],tab))
    tab=[]
    for m in methods:
        rr=[r for r in rows if r['method']==m];s=stats[m]
        eq=float(np.mean([F.stat([r for r in rr if r['cond']==c])['L'] for c in ('breath','shift_fast','flinch','cough','rib','combo')]))
        stats[m]['equal_condition_L']=eq
        tab.append([m,F.fmt(eq),F.fmt(s['calls50'],0)+'/'+F.fmt(s['calls99'],0),s['sat'],s['nonfinite']])
    lines+=['### 六工况等权与调用时间\n',F.table(['方法','等权L','回合调用p50/p99中位μs','饱和命令总数','非有限命令'],tab),
            '并行调用时间来自15进程竞争，p99为每回合p99的跨回合中位。未变方法的并行时间沿用第一轮，不能把两轮计时差解释为代码优化。下表全部按同一源日志重新做单进程暖机完整策略调用，包含ForceSensing，舍弃前125次；用于比较原型计算代价，不作实时认证。',
            '### 单进程暖机完整调用\n',F.table(['方法','p50 μs','p99 μs'],[[r['method'],F.fmt(r['p50_us'],0),F.fmt(r['p99_us'],0)] for r in diag['serial_full_calls']])]
    lines+=['### 第二轮实际运行量与耗时\n']
    stages=[['小样',len(B.load('baseline_pilot_r2.json')['rows']),B.load('baseline_pilot_r2.json')['seconds']/60],
            ['新增开发',len(dev['rows'])-len(old['rows']),sum(dev['round_seconds'].values())/60],
            ['配置变化闭环',168*len(test['rerun_methods']),sum(test['round_seconds'].values())/60],
            ['变化估计器双源回放',336*len(module['rerun_methods']),sum(module['round_seconds'].values())/60],
            ['同模型χ机制核查',len(mechanism['rows']),sum(mechanism['round_seconds'].values())/60]]
    lines.append(F.table(['阶段','回合/回放数','墙钟min'],[[a,b,F.fmt(c,2)] for a,b,c in stages]))
    lines+=['### MHE / Kalman 同模型与χ机制核查\n',
            '三个事后冻结的参数交换只用于解释，不进入选参。Kalman@MHE复制MHE最终过程、gate、σ；MHE@Kalman复制Kalman最终过程、gate、σ，N固定为MHE选中N；MHE-NIS@Kalman再将χ换为同样24拍的新息均值，并对齐宿主残差bound所用自由度22及启动。其余宿主相同。每项168个原配对格。对同一过程模型，无约束线性高斯MHE端点与Kalman相等；历史平滑轨迹只进入窗口残差门控（χ及其自由度），不作为额外端点状态输入控制器。',
            F.table(['固定比较','L','ΔL对当前[95%CI]','ΔL对原Kalman[95%CI]'],[
                [m,F.fmt(B.L([r for r in mechanism['rows'] if r['method']==m])),F.difftext(F.paired([r for r in mechanism['rows'] if r['method']==m],base)),F.difftext(F.paired([r for r in mechanism['rows'] if r['method']==m],{B.key(r):r for r in rows if r['method']=='Kalman'}))] for m in mechanism['completed']]),
            'MHE-NIS@Kalman对原Kalman的最大逐回合J差：'+str(mechanism['kalman_nis_max_J_error'])+'。这是统计量/启动对齐后的数值等价核查，不是另一个候选方法的胜负。']
    lines+=['## E. 结论与机制解释\n',(B.HERE/'baseline_interpretation.md').read_text(encoding='utf-8')]
    lines+=['## F. 不确定或可能不公平之处\n',
            '1. 已修正第一轮边界调参不足、MHE缺少二阶和Kalman优选邻域、外部辨识器门控菜单较窄的问题。实际新增配置及残留边界完整公开，不能再沿用“MHE机制导致选常数q=0”等第一轮归因。对手比内部Method拥有更多开发预算，预算不等对我方保守；24项上限仍限制高维轴的联合搜索，不是全局最优证书。',
            '2. 第二轮协议是在Claude检查第一轮开发边界之后制定；测试沿用相同168格，因而属于对已知调参缺口的复核，不是全新的最终盲测。严格只使用开发分数选择并冻结配置，仍不能消除第一轮测试已查看带来的研究者自由度。第三轮后残留边界按表披露，不继续延伸。',
            '3. gate、sigma_scale选择集一致，但受上限约束并非每个结构都搜索完整4×2。尤其MHE需分配预算给N和过程模型；不能声称所有方法搜索难度等同。平滑残差的自由度只是接口近似，强过程噪声可能以轨迹自由度吸收测量，较低χ不自动说明不确定度可信。',
            '4. MHE无物理不等式约束、到达信息精确线性高斯、只输出窗口端点；这是保留历史及到达代价的任务适配，不代表原Zhang控制器或约束MHE全部能力。同模型下所谓平滑收益通过χ进入宿主；不能把端点等价的两方法包装成独立估计能力。Haninger原动力学EKF也被适配为力/编码器联合KF，不是完整机器人复现。',
            '5. Xue额外使用同对象静止副本的2s接触位置测量标定，其准备成本和先验优势保留。真实参考深度未知、HC与黏弹对象不匹配、切线与20.8Hz探测增益不同；死区大、很少更新可能成为开发所选策略，不能因此推广“Xue路线失败”。固定σ没有概率覆盖保证。',
            '6. 模块真值col21没有计入探头惯性经机械臂柔顺的全部增益路径，已有约4–9%的共同评分偏差；H相位、相关噪声、模型误差和Huber冻结IRLS权重也没有完全进入σ。基础门控接受不等于必然更新前馈k_good。窗口高度相关，模块统计作描述，不能按数十万窗口当独立样本。',
            '7. 各分组CI未作多重比较校正；工况3或9簇较少。pooled快速工况权重更大，六工况等权另列。这里没有人体或机器人实物验证，不能将经验bootstrap区间推广到人体可靠性。warm前1秒不计力峰是共同评价边界；J截断不删除峰值、失联和失败。',
            '8. 完整Rossi、DO-MPAC、Xing采用各自控制架构，不能与Method中的同名辨识器互换归因。DO-MPAC名义极点、ESO极点检查不认证真实有延迟的闭环稳定性。15进程及串行Python计时均不是实时认证。',
            '## G. 文件与复现\n',
            '仅修改/新增baseline源码、协议和结果解释；未改仿真、场景、当前Method、冻结快照、原文。所有第一轮原始JSON哈希与原manifest一致，原日志和原模块NPZ保留。decision_log没有改动；未提交、未推送。第一轮尾部案例只保留在原报告，未将旧配置案例冒充新配置诊断。',
            '```powershell\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage checks\n# plan仅首次运行；已冻结时禁止重写\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage plan\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage pilot\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage dev\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage test\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage replay\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage diagnostics\npython -X utf8 experiments/rossi_failure_study/baseline_round2.py --stage mechanism\npython -X utf8 experiments/rossi_failure_study/baseline_report_r2.py\n```']
    mbase={B.key(r):r for r in rows if r['method']=='Kalman'}
    summary=dict(round=2,stats=stats,differences=diffs,groups=groups,module=modulesummary,reference='Huber-current',
                 mhe_vs_kalman=F.paired([r for r in rows if r['method']=='MHE'],mbase),
                 mechanism={m:dict(stats=F.stat([r for r in mechanism['rows'] if r['method']==m]),
                                  versus_kalman=F.paired([r for r in mechanism['rows'] if r['method']==m],mbase)) for m in mechanism['completed']})
    B.dump('baseline_summary_r2.json',summary);B.dump('baseline_validation_r2.json',validation)
    R.REPORT.write_text('\n\n'.join(lines),encoding='utf-8')
    manifest=B.load('baseline_manifest_r2.json');manifest['completed_utc']=datetime.now(timezone.utc).isoformat()
    manifest.update(workers=15,blas_threads=1,validation=validation,
                    final_code_hashes={str(p.relative_to(B.HERE)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in B.HERE.glob('*baseline*.py')},
                    result_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in B.OUT.glob('*r2*.json') if p.name!='baseline_manifest_r2.json'},
                    report_sha256=hashlib.sha256(R.REPORT.read_bytes()).hexdigest())
    B.dump('baseline_manifest_r2.json',manifest);print(R.REPORT,flush=True)


if __name__=='__main__':generate()
