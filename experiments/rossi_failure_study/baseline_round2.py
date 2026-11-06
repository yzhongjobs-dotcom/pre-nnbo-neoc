"""Predeclared boundary extension, restartable round 2/3; no Git operations.

Only development scores select models. Test, replay and diagnostic outputs are
separate from round 1. All worker pools have 15 processes, BLAS one thread.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import hashlib
import json
import math
import time
from datetime import datetime,timezone
from itertools import product
from pathlib import Path
import numpy as np
import baseline_compare as B
import id_baselines as I

REPORT=Path('D:/pip_env/temp/claude/8fbcb2674095e6be/841f52cd-14a8-4db4-b56f-36be6ffcc337/scratchpad/baseline_report_r2.md')
METHODS=('Kalman','MHE','Rossi-ID','Rossi-control','DO-MPAC','Xue-const','Xue-chi','Xing-DAHIC')
GATES=(.2,.551,1.,2.)
SIGMAS=(.796,1.5)
# Positive scale parameters use the original neighbouring multiplicative step;
# KI and n use the original outer arithmetic interval. Discrete structural
# choices (order,N,Np) and prescribed gate/sigma menus are not extended.
AXES={
 'Kalman':{'q_k':('log',10.),'q_d':('log',10.)},
 'MHE':{'q_k':('log',10.),'q_d':('log',10.)},
 'Rossi-ID':{'half':('log',2.),'memory':('log',4.)},
 'Rossi-control':{'K_I':('linear',30.),'memory':('log',4.)},
 'DO-MPAC':{'eso_bandwidth':('log',2.5),'q_over_r':('log',100.),'k_h0':('log',10.)},
 'Xue-const':{'delta':('log',10.),'gain':('log',1/.3),'initial_k':('log',10.),'initial_n':('linear',.26)},
 'Xue-chi':{'delta':('log',10.),'gain':('log',1/.3),'initial_k':('log',10.),'initial_n':('linear',.26)},
 'Xing-DAHIC':{'m':('log',3.5),'b':('log',3.),'k1':('log',10.),'k2':('log',100.)}}

RULES="""## 第二轮变化

### 运行前冻结的扩展规则

本节及 `results/baseline_plan_r2.json` 在任何第二轮开发运行之前写入。第二轮规则冻结后不据结果修改；第三轮按下述确定规则生成清单，并在第三轮开发之前落盘。

1. 每个外部方法第一轮选中连续参数若在该轴边界，沿边界方向延伸两步，与其他轴最优的一个值交叉，另加入各轴同一步数的联合延伸。每方法每轮最多 24 个**新增唯一配置**；第一轮配置继续参加选择。优先确保每条边界轴的两步，其次完整门控集合，再补与第二优结构的交叉。不追加内部 Huber-D4、D4、P5 的预算。当前 Method 固定 Huber-current，第一轮调后 Huber-D4 仅作参照；因此外部对手的预算更多，对我方保守。
2. 正值尺度轴按原相邻值倍率：Kalman/MHE q_k、q_d 为 ×10；Rossi-ID half ×2、memory ×4；Rossi-control K_I 按最外段等差 +30，遗忘时间 memory ×4（运行时转回 λ=exp(−0.008/memory)）；DO-MPAC 带宽 ×2.5、Q/R ×100、k_h0 ×10；Xue Δ ×10、gain ×(1/0.3)、initial_k ×10、initial_n 等差 0.26；Xing m ×3.5、b 按最外相邻近似对数间隔 ×3、k1 ×10、k2 ×100。只保留物理定义域内的正值。第一轮非等距轴采用上述最邻近外段步长，避免临时挑选。alpha=0.99 不在边界。order、N、Np 是离散结构，不以连续轴规则延伸；MHE 的 N 专门固定为 {16,24,48}。
3. Rossi-ID、Kalman、MHE、Xue-chi 均可选 gate∈{0.2,0.551,1,2}、sigma_scale∈{0.796,1.5}。每轮先给每个边界结构一个原门控，再在最优旧结构及最外联合结构上补完整 4×2；预算不足优先保留轴延伸。余量按固定结构顺序补门控。该例外意味着门控集合一致、结构交叉不必穷举；报告真实配置清单，不声称网格完全相同。Xue-const 无门控，不补此维度。
4. Kalman 先完成第二轮及有触发的第三轮，再生成 MHE 第二轮清单。MHE 必须覆盖 Kalman 第一轮及扩展后最终最优过程模型；若两者相同，另取开发集次优的不同过程模型。每个过程模型与 N={16,24,48} 交叉，先保留原两个 gate，再在扩展后最佳过程的 N=24/48 上补完整 4×2，合计不超过24。这是依赖开发集的预声明规则，不读取测试。MHE 第三轮若过程参数仍在其正值网格边界，则按两步扩展与三个 N 交叉，再在所选结构补门控到24。
5. 第二轮选中值若仍在新增数值边界（不含固定门控菜单和离散结构），允许第三轮再沿该轴延伸两步，仍最多24个；第三轮后不扩展，如实报告残留边界。只在8个原开发场景上，按全部已运行配置的平均 L=mean(log1p(J)) 选；不因测试分数回头改参数。相同分数按第一轮、第二轮、第三轮及配置原顺序稳定取最早者。
6. 选中配置改变的方法，在与第一轮相同的168个新种子测试格上重跑；未变者逐条复用第一轮测试结果。改变的辨识器在原 Huber-D4、D4 两个固定日志源各168份上回放，并重做同一人工阶跃。第一轮原始 JSON/NPZ 不覆盖。统计沿用5000次、按(工况,组织,i)物理场景簇聚类的配对bootstrap，其他评分/事件/删失口径不变。
7. 在有限正 q 下，新增二阶 MHE 采用与二阶 Kalman 完全相同的积分随机游走和精确离散非对角 Q。检查端点均值及方差一致性、窗口平移到达代价。MHE 通过平滑轨迹形成窗口测量残差 χ，控制器只接端点 k、σ、χ，不消费历史平滑轨迹。新增固定参数机制核查：同一过程模型、门控尺度下比较 Kalman 与 MHE；将 MHE χ 替换为同样的24拍新息统计应复现 Kalman 的端点和闭环。此核查不参与选参，以区分窗口 χ 与端点估计收益。
8. 先小样计时；正式计算用15进程、BLAS单线程。串行暖机调用另计；本机原型计时不作实时认证。不提交、不推送、不改 decision_log，交由 Claude 验收。
"""

def canon(method,cfg):
    c=dict(cfg)
    if method in ('Kalman','MHE','Rossi-ID','Xue-chi'):
        c.setdefault('gate',.551);c.setdefault('sigma_scale',.796)
    if method=='MHE':c.setdefault('d_order',1)
    return json.dumps(c,sort_keys=True)

def expanded_cfg(method,cfg):
    c=dict(cfg)
    if method=='Rossi-control':c['memory']=-.008/math.log(c.pop('lam'))
    return c

def executable_cfg(method,cfg):
    c=dict(cfg)
    if method=='Rossi-control':c['lam']=math.exp(-.008/c.pop('memory'))
    return c

def external(method):return method in ('Kalman','MHE','Rossi-ID','Xue-chi')

def boundary_axes(method,chosen,configs):
    selected=expanded_cfg(method,chosen);cfgs=[expanded_cfg(method,c) for c in configs];axes={}
    for axis,(kind,step) in AXES[method].items():
        if axis not in selected:continue
        vals=[c[axis] for c in cfgs if axis in c and c[axis]>0]
        v=selected[axis]
        if not vals or v<=0:continue
        lo,hi=min(vals),max(vals)
        direction=-1 if math.isclose(v,lo,rel_tol=1e-8) else (1 if math.isclose(v,hi,rel_tol=1e-8) else 0)
        if direction:
            vv=[v*step**(direction*j) if kind=='log' else v+direction*j*step for j in (1,2)]
            axes[axis]=[float(x) for x in vv if x>0]
    return axes

def general_configs(method,chosen,configs,scores,round_number):
    axes=boundary_axes(method,chosen,configs)
    if round_number==3 and not axes:return [],axes
    root=expanded_cfg(method,chosen);structures=[root]
    for axis,values in axes.items():
        structures.extend(dict(root,**{axis:x}) for x in values)
    for j in (0,1):
        structures.append(dict(root,**{a:v[j] for a,v in axes.items() if len(v)>j}))
    structures=[executable_cfg(method,c) for c in structures]
    seen={canon(method,c) for c in configs};out=[]
    def add(c):
        s=canon(method,c)
        if s not in seen and len(out)<24:seen.add(s);out.append(c)
    for c in structures:add(c)
    if external(method):
        for c in (structures[0],structures[-1]):
            for g,s in product(GATES,SIGMAS):add(dict(c,gate=g,sigma_scale=s))
        for c in structures:
            for g,s in product(GATES,SIGMAS):add(dict(c,gate=g,sigma_scale=s))
    # Cross each boundary axis with the second best existing configuration,
    # then all previously declared boundary values in stable product order.
    order=np.argsort(scores,kind='stable')
    for ix in order[:2]:
        c=expanded_cfg(method,configs[int(ix)])
        for axis,values in axes.items():
            for v in values:add(executable_cfg(method,dict(c,**{axis:v})))
    choices=[(a,[root[a]]+v) for a,v in axes.items()]
    for values in product(*(v for _,v in choices)):
        add(executable_cfg(method,dict(root,**dict(zip((a for a,_ in choices),values)))))
    return out,axes

def mhe_configs(dev):
    first=B.load('baseline_dev.json')['chosen']['Kalman'];best=dev['chosen']['Kalman']
    fields=('q_k','q_d','d_order');models=[{k:c[k] for k in fields} for c in (first,best)]
    if models[0]==models[1]:
        ss=dev['scores']['Kalman']
        models[0]=next({k:ss['configs'][j][k] for k in fields} for j in np.argsort(ss['all'],kind='stable') if any(ss['configs'][j][k]!=best[k] for k in fields))
    out=[];seen=set()
    def add(c):
        s=canon('MHE',c)
        if s not in seen and len(out)<24:seen.add(s);out.append(c)
    for model,n,g in product(models,(16,24,48),(.551,1.)):
        add(dict(model,N=n,gate=g,sigma_scale=.796))
    for n,g,s in product((24,48),GATES,SIGMAS):
        add(dict(models[-1],N=n,gate=g,sigma_scale=s))
    return out,dict(kalman_first_or_neighbour=models[0],kalman_expanded_best=models[-1])

def plan():
    assert not (B.OUT/'baseline_plan_r2.json').exists(),'plan already frozen'
    dev=B.load('baseline_dev.json');entries={}
    for m in METHODS:
        if m=='MHE':continue
        s=dev['scores'][m]
        cfgs,axes=general_configs(m,dev['chosen'][m],s['configs'],s['all'],2)
        entries[m]=dict(configs=cfgs,axes=axes)
    obj=dict(round=2,frozen_utc=datetime.now(timezone.utc).isoformat(),rules=RULES,methods=entries,
             mhe='deterministic development-only dependent plan; materialized after Kalman rounds 2/3')
    B.dump('baseline_plan_r2.json',obj)
    old=REPORT.with_name('baseline_report.md').read_text(encoding='utf-8')
    listing='\n### 已冻结的第二轮新增配置（MHE 按上述依赖规则补入）\n\n```json\n'+json.dumps(entries,ensure_ascii=False,indent=1)+'\n```\n'
    REPORT.write_text(old.split('## A.')[0]+RULES+listing+'\n## A.'+old.split('## A.',1)[1],encoding='utf-8')
    protocol=B.HERE/'baseline_protocol.md'
    protocol.write_text(protocol.read_text(encoding='utf-8')+'\n\n'+RULES+listing,encoding='utf-8')
    print({m:len(v['configs']) for m,v in entries.items()},flush=True)

def all_dev():
    path=B.OUT/'baseline_dev_r2.json'
    return B.load(path.name) if path.exists() else B.load('baseline_dev.json')

def run_round(method,number):
    dev=all_dev();s=dev['scores'][method]
    marker=f'{method}:r{number}'
    if marker in dev.get('completed_rounds',[]):return
    if number==2:
        if method=='MHE':configs,axes=mhe_configs(dev)
        else:
            entry=B.load('baseline_plan_r2.json')['methods'][method];configs=entry['configs'];axes=entry['axes']
    elif method=='MHE':
        axes=boundary_axes(method,dev['chosen'][method],s['configs']);configs=[]
        if axes:
            root=dev['chosen'][method]
            for a,vs in axes.items():
                for v,n in product(vs,(16,24,48)):configs.append(dict(root,**{a:v,'N':n}))
            # Fill the menu on selected/boundary structures, preserving all axis trials.
            for c in [root]+configs.copy():
                for g,scale in product(GATES,SIGMAS):
                    candidate=dict(c,gate=g,sigma_scale=scale)
                    if len(configs)<24 and canon(method,candidate) not in {canon(method,x) for x in s['configs']+configs}:configs.append(candidate)
    else:configs,axes=general_configs(method,dev['chosen'][method],s['configs'],s['all'],number)
    old={canon(method,c) for c in s['configs']};configs=[c for c in configs if canon(method,c) not in old]
    assert len(configs)<=24 and len({canon(method,c) for c in configs})==len(configs)
    plan_name=f'baseline_plan_{method}_r{number}.json'
    payload=dict(round=number,method=method,configs=configs,axes=axes,frozen_utc=datetime.now(timezone.utc).isoformat())
    if (B.OUT/plan_name).exists():
        previous=B.load(plan_name);assert previous['configs']==configs
    else:
        B.dump(plan_name,payload)
        with REPORT.open('a',encoding='utf-8') as f:f.write('\n\n<!-- 运行前冻结 '+plan_name+' -->\n```json\n'+json.dumps(payload,ensure_ascii=False,indent=1)+'\n```\n')
    if configs:
        cells=[(c,t,'neural_white',21.,i) for c,t in B.D.DEV for i in range(2)]
        started=time.perf_counter()
        rows=B.run([(method,c,cell,B.D.BASE,0.,False,True) for c in configs for cell in cells])
        for r in rows:r['round']=number
        vals=[B.L(rows[j*8:(j+1)*8]) for j in range(len(configs))]
        s['configs']+=configs;s['all']+=vals;index=int(np.argmin(s['all']))
        s.update(index=index,best=s['all'][index]);dev['chosen'][method]=s['configs'][index]
        dev['rows']+=rows
        dev.setdefault('round_seconds',{})[marker]=time.perf_counter()-started
    dev.setdefault('completed_rounds',[]).append(marker);dev['round']=2
    dev.setdefault('added_configs',{})[marker]=configs
    B.dump('baseline_dev_r2.json',dev)
    print(marker,'ADDED',len(configs),'BEST',s['best'],dev['chosen'][method],flush=True)

def develop():
    # Explicit ordering is scientific: MHE covers the final expanded Kalman model.
    for m in METHODS:
        run_round(m,2);run_round(m,3)

def changed():
    d=all_dev()['chosen'];old=B.load('baseline_dev.json')['chosen']
    return [m for m in METHODS if canon(m,d[m])!=canon(m,old[m])]

def pilot():
    d=B.load('baseline_plan_r2.json');jobs=[]
    for m in METHODS:
        c=(dict(N=48,q_k=300.,q_d=30.,d_order=2,gate=.551) if m=='MHE' else d['methods'][m]['configs'][0])
        jobs.append((m,c,('breath','relaxed','neural_white',21.,0),B.D.BASE,0.,False,True))
    tic=time.perf_counter();rows=B.run(jobs)
    B.dump('baseline_pilot_r2.json',dict(rows=rows,seconds=time.perf_counter()-tic,workers=15,blas_threads=1))
    for r in rows:print(r['method'],r['seconds'],'seconds',flush=True)

def test():
    original=B.load('baseline_test.json');dev=all_dev();changes=changed()
    saved=B.load('baseline_test_r2.json') if (B.OUT/'baseline_test_r2.json').exists() else dict(original,rows=[r for r in original['rows'] if r['method'] not in changes],round=2,rerun_methods=[],reused_methods=[m for m in original['chosen'] if m not in changes],round_seconds={})
    for m in changes:
        if m in saved['rerun_methods']:continue
        tic=time.perf_counter();rr=B.run([(m,dev['chosen'][m],cell,B.TEST_BASE,0.,False,False) for cell in B.cells()])
        for r in rr:r['round']=2
        saved['rows']+=rr;saved['chosen'][m]=dev['chosen'][m];saved['rerun_methods'].append(m)
        saved['round_seconds'][m]=time.perf_counter()-tic
        saved['seconds']=sum(saved['round_seconds'].values());B.dump('baseline_test_r2.json',saved)
        print(m,'TEST L',B.L(rr),flush=True)

def replay():
    changes=[m for m in changed() if m in B.ID_METHODS];original=B.load('baseline_module.json')
    saved=B.load('baseline_module_r2.json') if (B.OUT/'baseline_module_r2.json').exists() else dict(round=2,rows=[r for r in original['rows'] if r['method'] not in changes],rerun_methods=[],round_seconds={})
    sources=[r for r in B.load('baseline_test.json')['rows'] if r['method'] in ('Huber-D4','D4') and 'log_file' in r]
    chosen=all_dev()['chosen'];folder=B.OUT/'baseline_replay_r2';folder.mkdir(exist_ok=True)
    for m in changes:
        if m in saved['rerun_methods']:continue
        tic=time.perf_counter();rr=B.run([(s,m,chosen[m],0.) for s in sources],worker=B.replay_job)
        for r in rr:
            data={}
            for kind,g in r['groups'].items():
                for field in ('error_samples','accepted_errors','times'):data[kind+'_'+field]=np.array(g.pop(field))
            path=folder/f"{m}_{r['source']}_{r['cond']}_{r['dynamics']}_{int(r['phase'])}_{r['i']}.npz"
            np.savez_compressed(path,**data);r['samples_file']=str(path.relative_to(B.OUT)).replace('\\','/');r['round']=2
        saved['rows']+=rr;saved['rerun_methods'].append(m);saved['round_seconds'][m]=time.perf_counter()-tic
        saved['seconds']=sum(saved['round_seconds'].values());B.dump('baseline_module_r2.json',saved)

def checks():
    from baseline_checks import checks as oldchecks
    rng=np.random.default_rng(724);Y=rng.normal(0,.0002,64);z=2300*Y+1.2+rng.normal(0,.04,64)
    records=[]
    for qk,qd,n in product((3.,30.,300.),(3.,30.,300.),(16,24,48)):
        x=np.array([30.,0.,0.]);P=np.diag([90000.,10000.,10000.]);A,Q=I.process_model(qk,qd,.008,2)
        priors=[];post=[]
        for y,zz in zip(Y,z):
            priors.append((A@x,A@P@A.T+Q))
            x,P,_=I.kf_step(x,P,np.array([1000*y,1.,0.]),zz,.04**2,A,Q);post.append((x.copy(),P.copy()))
        a=len(Y)-n;tr,var,_=I.mhe_solve(Y[a:],z[a:],.04**2,*priors[a],qk,qd,.008,2)
        err=float(np.max(abs(tr[-1]-x)));ve=abs(var-P[0,0]*1e6)
        assert np.allclose(tr[-1],x,atol=2e-5,rtol=2e-6),(qk,qd,n,err)
        assert np.isclose(var,P[0,0]*1e6,atol=.01,rtol=2e-5),(qk,qd,n,ve)
        records.append(dict(q_k=qk,q_d=qd,N=n,endpoint_abs_error=err,variance_abs_error=ve,variance_relative_error=ve/(P[0,0]*1e6)))
    old=json.loads(json.dumps(oldchecks(),default=str))
    result=dict(first_order_regression=old,second_order_finite_q=records)
    B.dump('baseline_checks_r2.json',result);print('PASSED',len(records),'second order checks',flush=True)

def diagnostics():
    import baseline_diagnostics as D
    original=B.load('baseline_diagnostics.json');changes=changed();chosen=all_dev()['chosen']
    rows=[r for r in original['rows'] if r['method'] not in changes]
    D.diagnostic('D4',chosen['D4'],False)
    for m in changes:
        if m in B.ID_METHODS:
            for moving in (False,True):rows.append(D.diagnostic(m,chosen[m],moving))
    saved=B.load('baseline_test_r2.json');source=next(r for r in saved['rows'] if r['method']=='Huber-D4' and r['cond']=='breath')
    with np.load(B.OUT/source['log_file']) as z:
        log=np.full((len(z['log']),26),np.nan);log[:,z['cols']]=z['log'];observed=z['observed'].copy()
    sc=B.scene(B.key(source),B.TEST_BASE,0.)
    sc['initial_contact_observed']=B.load('baseline_calibration.json')['positions'][B.calkey(B.key(source),False)]
    serial=[]
    for m,c in saved['chosen'].items():
        pol=B.build(m,c,sc);pol.reset(dict(B.P3.PLATFORM,tau_a=sc['tau_a_hat']))
        for n in range(1,len(log)+1):pol(B.source_observation(log,observed,n))
        calls=np.array(pol.t_call)[125:]*1e6
        serial.append(dict(method=m,p50_us=float(np.median(calls)),p99_us=float(np.percentile(calls,99))))
    B.dump('baseline_diagnostics_r2.json',dict(rows=rows,serial_full_calls=serial,cases=[],note='Unchanged synthetic rows reused; all full calls freshly timed serially; round-1 cases remain only in round-1 report.'))


def mechanism_job(job):
    label,method,cfg,cell=job
    r=B.execute((method,cfg,cell,B.TEST_BASE,0.,False,False));r['method']=label
    return r


def mechanism():
    """Frozen post-selection process/chi swaps on the paired test grid.

    These are explanations of selected comparisons, never new selections.
    """
    chosen=all_dev()['chosen'];k=chosen['Kalman'];m=chosen['MHE']
    km={a:v for a,v in m.items() if a in ('q_k','q_d','d_order','gate','sigma_scale')}
    mk=dict(k,N=m['N'])
    variants=[('Kalman@MHE','Kalman',km),('MHE@Kalman','MHE',mk),('MHE-NIS@Kalman','MHE',dict(mk,chi_mode='nis'))]
    path=B.OUT/'baseline_mechanism_r2.json'
    saved=B.load(path.name) if path.exists() else dict(round=2,chosen=chosen,variants=variants,rows=[],completed=[],round_seconds={})
    for label,method,cfg in variants:
        if label in saved['completed']:continue
        tic=time.perf_counter();rr=B.run([(label,method,cfg,c) for c in B.cells()],worker=mechanism_job)
        saved['rows']+=rr;saved['completed'].append(label);saved['round_seconds'][label]=time.perf_counter()-tic
        if label=='MHE-NIS@Kalman':
            base={B.key(r):r for r in B.load('baseline_test_r2.json')['rows'] if r['method']=='Kalman'}
            errors=[abs(r['J']-base[B.key(r)]['J']) for r in rr]
            saved['kalman_nis_max_J_error']=max(errors);assert max(errors)<1e-9,max(errors)
        B.dump(path.name,saved);print(label,'L',B.L(rr),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True,choices=['plan','checks','pilot','dev','test','replay','diagnostics','mechanism']);a=p.parse_args()
    globals()['develop' if a.stage=='dev' else a.stage]()
