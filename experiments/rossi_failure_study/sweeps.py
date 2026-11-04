"""One-command tuning, paired closed-loop evaluation, condition scans and report."""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
for s in ('motion_aware_scan', 'method_redesign', 'cobot_scan_benchmark'):
    sys.path.insert(0, str(HERE.parent/s))
import physics_plant as PP
import scenes as SC
import common
from cohort_eval import FINAL
from failure_replay import Measured, OracleMeasured, cell_scene
from estimator_swap import SwapPolicy
from estimators import ProjectionEstimator, BandPassEstimator, R2Estimator

CELLS = [(c, t) for c in ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor')
         for t in ('soft', 'mid', 'hard')] + [('rib', 'mixed'), ('combo', 'mixed')]
HEAD_CELLS = [('breath', 'mid'), ('breath', 'hard'), ('shift_fast', 'mid'), ('shift_fast', 'hard'),
              ('flinch', 'mid'), ('flinch', 'hard'), ('cough', 'mid'), ('cough', 'hard'),
              ('rib', 'mixed'), ('combo', 'mixed')]
METHODS = ('Projection', 'BandPass', 'R2', 'Oracle')
FINAL_C = dict(FINAL, f_exc=20.8)


def estimator(name, params=None):
    p = dict(params or {})
    if name == 'Projection': return ProjectionEstimator(**p)
    if name == 'BandPass': return BandPassEstimator(**p)
    if name == 'R2': return R2Estimator(**p)
    raise ValueError(name)


def policy(method, controls, est_params=None, oracle=False):
    p = dict(FINAL_C); p.update(controls)
    if oracle:
        p.update(oracle_meas=True, oracle_fb=True)
    return SwapPolicy(estimator(method, est_params),
                      single_model=False, **p)


def _run(sc, method, controls, est_params, seed, platform=None, oracle=False):
    pol = policy(method, controls, est_params, oracle=oracle)
    wrap = OracleMeasured(pol) if oracle else Measured(pol)
    log = PP.episode(sc, wrap, seed, platform=platform)
    return log, pol


def _metric(log, pol, sc, f_exc):
    if log is None:
        return dict(J=1000., rmse=np.nan, peak=np.nan, osc=np.nan, ripple=np.nan,
                    long_loss=True, failed=True, event_rmse=np.nan, call_p50=np.nan,
                    call_p99=np.nan, est_p50=np.nan, est_p90=np.nan,
                    est_over_ln12=np.nan, accept_rate=0.)
    m = common.task_metrics(log, f_exc)
    t = log[:, 0]; warm = t >= 1.
    ev = np.zeros(len(log), bool)
    for _, a0, a1 in sc.get('events', []): ev |= (t >= a0) & (t <= a1)
    em = warm & ev
    event_rmse = float(np.sqrt(np.mean((log[em, 3]-4.)**2))) if em.sum() else np.nan
    calls = np.asarray(getattr(pol, 't_call', []), float)*1e6
    accepted = []
    for tt, kh, sig, chi, ok in getattr(pol, 'accept_trace', []):
        idx = int(np.clip(np.searchsorted(t, tt), 0, len(t)-1)); kt = log[idx, 8]
        if ok and np.isfinite(kh) and kt > 0: accepted.append(math.log(kh/kt))
    ea = np.asarray(accepted, float)
    all_loss = log[warm, 3] < .5
    longest = 0; cur = 0
    for x in all_loss:
        cur = cur+1 if x else 0; longest = max(longest, cur)
    attempts = len(getattr(pol, 'accept_trace', []))
    return dict(J=float(min(m['J'], 1000.)), rmse=float(m['rmse']), peak=float(m['peak']),
                osc=float(m.get('osc_ex', np.nan)), ripple=float(m.get('ripple', np.nan)),
                long_loss=bool(longest*0.008 > .5), failed=bool(m.get('diverged', 0)),
                event_rmse=event_rmse, call_p50=float(np.percentile(calls, 50)) if len(calls) else np.nan,
                call_p99=float(np.percentile(calls, 99)) if len(calls) else np.nan,
                est_p50=float(np.percentile(abs(ea), 50)) if len(ea) else np.nan,
                est_p90=float(np.percentile(abs(ea), 90)) if len(ea) else np.nan,
                est_over_ln12=float(np.mean(abs(ea) > math.log(1.2))) if len(ea) else np.nan,
                accept_rate=float(sum(x[4] for x in getattr(pol, 'accept_trace', []))/attempts) if attempts else 0.)


def _score_scene(sc, method, controls, est_params, seed, oracle=False, platform=None):
    log, pol = _run(sc, method, controls, est_params, seed, platform=platform, oracle=oracle)
    return _metric(log, pol, sc, controls.get('f_exc', 20.8))


def _cell_score(values):
    x = [v for v in values if np.isfinite(v)]
    return float(np.mean(x)) if x else 1000.


def _public_default(kind):
    p = dict(FINAL_C)
    if kind == 'neutral':
        p.update(c=FINAL['c'], rel_max=FINAL['rel_max'], chi2_max=FINAL['chi2_max'], f_exc=20.8,
                 a_force=FINAL['a_force'])
    return p


def _public_random(kind, rng):
    p = _public_default(kind)
    p.update(D1=float(np.exp(rng.uniform(np.log(.008), np.log(.25)))),
             q_smooth=float(np.exp(rng.uniform(np.log(1e-7), np.log(3e-3)))),
             q_maneuver=float(np.exp(rng.uniform(np.log(.02), np.log(12.)))),
             p_stay=float(rng.uniform(.93, .9999)), r_scale=float(np.exp(rng.uniform(np.log(.03), np.log(12.)))),
             rho_d=float(np.exp(rng.uniform(np.log(.003), np.log(.25)))))
    if kind == 'projection':
        p.update(f_exc=float(rng.uniform(12., 28.)), a_force=float(rng.uniform(.05, .35)),
                 c=float(np.exp(rng.uniform(np.log(.1), np.log(3.)))),
                 rel_max=float(rng.uniform(.10, .55)), chi2_max=float(np.exp(rng.uniform(np.log(.8), np.log(10.)))))
    return p


def _public_local(kind, p, rng):
    q = dict(p)
    for k, lo, hi in (('D1', .008, .25), ('q_smooth', 1e-7, 3e-3), ('q_maneuver', .02, 12.),
                      ('r_scale', .03, 12.), ('rho_d', .003, .25)):
        q[k] = float(np.clip(q[k]*np.exp(rng.normal(0, .35)), lo, hi))
    q['p_stay'] = float(np.clip(q['p_stay']+rng.normal(0, .01), .93, .9999))
    if kind == 'projection':
        for k, lo, hi in (('a_force', .05, .35), ('c', .1, 3.), ('rel_max', .10, .55), ('chi2_max', .8, 10.)):
            q[k] = float(np.clip(q[k]*np.exp(rng.normal(0, .25)), lo, hi))
        q['f_exc'] = float(np.clip(q['f_exc']+rng.normal(0, 1.2), 12., 28.))
    return q


def tune_public_worker(a):
    kind, candidate, base = a; vals = []
    for ci, (c, t) in enumerate(CELLS):
        cell = []
        for i in range(3):
            sc = cell_scene(c, t, i, base); seed = base+ci*100+i
            m = _score_scene(sc, 'Projection', candidate, None, seed, oracle=(kind == 'neutral'))
            cell.append(math.log1p(min(m['J'], 1000.)))
        vals.append(float(np.mean(cell)))
    return dict(params=candidate, score=float(np.mean(vals)), cell_scores=vals)


def tune_public(kind, seed=340000, expanded=False):
    rng = np.random.default_rng(seed); candidates = [_public_default(kind)]
    candidates += [_public_random(kind, rng) for _ in range(40)]
    if expanded and kind == 'neutral':
        for p in candidates[1:]: p['p_stay'] = float(rng.uniform(.93, .99999))
    args = [(kind, p, seed) for p in candidates]
    with ProcessPoolExecutor(max_workers=15) as pool: first = list(pool.map(tune_public_worker, args, chunksize=1))
    best = min(first, key=lambda r: r['score'])
    candidates += [_public_local(kind, best['params'], rng) for _ in range(20)]
    if expanded and kind == 'neutral':
        for p in candidates[41:]: p['p_stay'] = float(np.clip(p['p_stay']+rng.normal(0,.002), .93, .99999))
    args = [(kind, p, seed) for p in candidates[41:]]
    with ProcessPoolExecutor(max_workers=15) as pool: tail = list(pool.map(tune_public_worker, args, chunksize=1))
    rows = first+tail; best = min(rows, key=lambda r: r['score'])
    return dict(kind=kind, seed=seed, budget='default+40 random+20 local', rows=rows, best=best)


def _est_default(name):
    return {'Projection': dict(nw=16, degree=3, sigma_scale=1., chi_scale=1.),
            'BandPass': dict(order=2, half=4., memory=.4, sigma_scale=1., chi_scale=1.),
            'R2': dict(Q=3., memory=.4, sigma_scale=1., chi_scale=1.)}[name]


def _est_random(name, rng):
    p = _est_default(name)
    if name == 'Projection': p.update(nw=int(rng.integers(12, 41)), degree=int(rng.integers(2, 6)))
    elif name == 'BandPass': p.update(order=int(rng.choice([1,2,3,4])), half=float(rng.uniform(2,12)), memory=float(np.exp(rng.uniform(np.log(.03), np.log(1.)))))
    else: p.update(Q=float(rng.uniform(1,12)), memory=float(np.exp(rng.uniform(np.log(.05), np.log(2.)))))
    p['sigma_scale'] = float(np.exp(rng.uniform(np.log(.2), np.log(5.))))
    p['chi_scale'] = float(np.exp(rng.uniform(np.log(.2), np.log(5.))))
    return p


def _est_local(name, p, rng):
    q = dict(p)
    if name == 'Projection':
        q['nw'] = int(np.clip(round(q['nw']+rng.normal(0, 3)), 12, 40)); q['degree'] = int(np.clip(round(q['degree']+rng.normal(0, .6)), 2, 5))
    elif name == 'BandPass':
        q['order'] = int(rng.choice([max(1,q['order']-1), q['order'], min(4,q['order']+1)])); q['half'] = float(np.clip(q['half']*np.exp(rng.normal(0,.2)),2,12)); q['memory'] = float(np.clip(q['memory']*np.exp(rng.normal(0,.3)),.03,1.))
    else:
        q['Q'] = float(np.clip(q['Q']*np.exp(rng.normal(0,.2)),1,12)); q['memory'] = float(np.clip(q['memory']*np.exp(rng.normal(0,.3)),.05,2.))
    q['sigma_scale'] = float(np.clip(q['sigma_scale']*np.exp(rng.normal(0,.25)),.2,5.)); q['chi_scale'] = float(np.clip(q['chi_scale']*np.exp(rng.normal(0,.25)),.2,5.))
    return q


def tune_est_worker(a):
    name, controls, estp, oracle, base = a; vals = []
    for ci, (c, t) in enumerate(CELLS):
        cell = []
        for i in range(3):
            sc = cell_scene(c, t, i, base); seed = base+ci*100+i
            m = _score_scene(sc, name, controls, estp, seed, oracle=oracle)
            cell.append(math.log1p(min(m['J'], 1000.)))
        vals.append(float(np.mean(cell)))
    return dict(params=estp, score=float(np.mean(vals)), cell_scores=vals)


def tune_estimator(name, controls, oracle, seed=340000, expanded=False):
    name_code = sum(ord(c)*(i+1) for i, c in enumerate(name))
    rng = np.random.default_rng(seed+name_code); candidates = [_est_default(name)] + [_est_random(name, rng) for _ in range(20)]
    if expanded:
        for p in candidates[1:]:
            if name == 'BandPass': p['memory'] = float(np.exp(rng.uniform(np.log(.03), np.log(3.))))
            p['sigma_scale'] = float(np.exp(rng.uniform(np.log(.05), np.log(5.))))
    args = [(name, controls, p, oracle, seed) for p in candidates]
    with ProcessPoolExecutor(max_workers=15) as pool: first = list(pool.map(tune_est_worker, args, chunksize=1))
    best = min(first, key=lambda r: r['score']); candidates += [_est_local(name, best['params'], rng) for _ in range(9)]
    if expanded:
        for p in candidates[21:]:
            if name == 'BandPass': p['memory'] = float(np.clip(p['memory']*np.exp(rng.normal(0,.35)), .03, 3.))
            p['sigma_scale'] = float(np.clip(p['sigma_scale']*np.exp(rng.normal(0,.3)), .05, 5.))
    args = [(name, controls, p, oracle, seed) for p in candidates[21:]]
    with ProcessPoolExecutor(max_workers=15) as pool: tail = list(pool.map(tune_est_worker, args, chunksize=1))
    rows = first+tail; best = min(rows, key=lambda r: r['score'])
    return dict(name=name, seed=seed, budget='default+20 random+9 local', rows=rows, best=best)


def expand_boundary_tuning(result):
    """Re-run only candidates whose first search selected a boundary point."""
    extra = {}
    if result['public']['C_neutral']['best']['params']['p_stay'] >= .9999-1e-12:
        extra['C_neutral_public'] = tune_public('neutral', seed=340001, expanded=True)
    for config, oracle in (('C_neutral', True), ('C_projection', False)):
        if config == 'C_projection' and result['estimators'][config]['BandPass']['best']['params']['memory'] >= 1.-1e-12:
            extra[f'{config}_BandPass'] = tune_estimator('BandPass', result['public'][config]['best']['params'], oracle, seed=340002, expanded=True)
        if result['estimators'][config]['R2']['best']['params']['sigma_scale'] <= .2+1e-12:
            extra[f'{config}_R2'] = tune_estimator('R2', result['public'][config]['best']['params'], oracle, seed=340003, expanded=True)
    result['expanded_boundary'] = extra
    (HERE/'results'/'tuning.json').write_text(json.dumps(result, indent=1, allow_nan=True), encoding='utf-8')
    return result


def do_tuning():
    neutral = tune_public('neutral'); projection = tune_public('projection')
    result = {'public': {'C_neutral': neutral, 'C_projection': projection}, 'estimators': {}}
    for config, pub in (('C_neutral', neutral), ('C_projection', projection)):
        result['estimators'][config] = {}
        for name in ('Projection', 'BandPass', 'R2'):
            result['estimators'][config][name] = tune_estimator(name, pub['best']['params'], False)   # the estimator must be in the loop
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/'tuning.json').write_text(json.dumps(result, indent=1, allow_nan=True), encoding='utf-8')
    return expand_boundary_tuning(result)


def _main_job(a):
    config, c, t, i, base, controls, eps = a; sc = cell_scene(c, t, i, base); seed = base + sum(ord(x) for x in c+t)*100+i
    rows = []
    for name in ('Projection', 'BandPass', 'R2'):
        log, pol = _run(sc, name, controls, eps[name], seed, oracle=False)
        rows.append(dict(config=config, method=name, cond=c, tissue=t, i=i, seed=seed,
                         metrics=_metric(log, pol, sc, controls.get('f_exc',20.8))))
    log, pol = _run(sc, 'Projection', controls, eps['Projection'], seed, oracle=True)
    rows.append(dict(config=config, method='Oracle', cond=c, tissue=t, i=i, seed=seed,
                     metrics=_metric(log, pol, sc, controls.get('f_exc',20.8))))
    return rows


def run_main(tuning, n=15, base=380000):
    args=[]
    for config, key in (('C_neutral','C_neutral'),('C_projection','C_projection')):
        controls=tuning['public'][key]['best']['params']; eps={k:tuning['estimators'][key][k]['best']['params'] for k in ('Projection','BandPass','R2')}
        args += [(config,c,t,i,base,controls,eps) for c,t in CELLS for i in range(n)]
    with ProcessPoolExecutor(max_workers=15) as pool: groups=list(pool.map(_main_job,args,chunksize=1))
    rows=[r for g in groups for r in g]; (HERE/'results'/'main_rows.json').write_text(json.dumps(rows,indent=1,allow_nan=True),encoding='utf-8'); return rows


COND_SCANS = {'nominal': {}, 'severity_.25': {'severity':.25}, 'severity_.5': {'severity':.5}, 'severity_2': {'severity':2.},
              'tremor_16': {'tremor':True,'tremor_f':16.,'tremor_amp':2e-5}, 'tremor_20.8': {'tremor':True,'tremor_f':20.8,'tremor_amp':2e-5},
              'tremor_26': {'tremor':True,'tremor_f':26.,'tremor_amp':2e-5}, 'tremor_40': {'tremor':True,'tremor_f':40.,'tremor_amp':2e-5},
              'ripple_20.8': {'ripple_amp':5e-6,'ripple_f':20.8}, 'noise_2': {'noise_scale':2.}, 'kr_3e4': {'kr':3e4},
              'lag_.6': {'tau_a':.019*.6}, 'lag_1.5': {'tau_a':.019*1.5}}


def _cond_job(a):
    cond, c, t, i, base, controls, eps, over = a; scene_over={k:v for k,v in over.items() if k != 'tau_a'}; sc=cell_scene(c,t,i,base,scene_over); plat=dict(PP.PLATFORM); plat['tau_a']=over.get('tau_a',plat['tau_a']); seed=base+sum(ord(x) for x in cond+c+t)*100+i
    rows=[]
    for name in ('Projection','BandPass','R2'):
        log,pol=_run(sc,name,controls,eps[name],seed,platform=plat,oracle=False); rows.append(dict(scan=cond,method=name,cond=c,tissue=t,i=i,seed=seed,metrics=_metric(log,pol,sc,controls.get('f_exc',20.8))))
    return rows


def run_conditions(tuning,n=8,base=390000):
    pub=tuning['public']['C_neutral']['best']['params']; eps={k:tuning['estimators']['C_neutral'][k]['best']['params'] for k in ('Projection','BandPass','R2')}; args=[(cond,c,t,i,base,pub,eps,over) for cond,over in COND_SCANS.items() for c,t in HEAD_CELLS for i in range(n)]
    with ProcessPoolExecutor(max_workers=15) as pool: groups=list(pool.map(_cond_job,args,chunksize=1))
    rows=[r for g in groups for r in g]; (HERE/'results'/'condition_rows.json').write_text(json.dumps(rows,indent=1,allow_nan=True),encoding='utf-8'); return rows


def bootstrap(x, seed=390001, n=10000):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    if not len(x): return [np.nan,np.nan]
    rng=np.random.default_rng(seed); means=np.mean(x[rng.integers(0,len(x),size=(n,len(x)))],axis=1); return [float(np.percentile(means,2.5)),float(np.percentile(means,97.5))]


def report(tuning):
    out=HERE/'estimator_swap'; out.mkdir(exist_ok=True); main=json.loads((HERE/'results'/'main_rows.json').read_text()); cond=json.loads((HERE/'results'/'condition_rows.json').read_text()) if (HERE/'results'/'condition_rows.json').exists() else []
    lines=['# Rossi 频率分离与多项式投影：统一控制器比较','','本报告只比较同一 `SwapPolicy` 与同一组控制参数下的刚度估计器差异。C_neutral 使用显式真值增益参考，属于特权参考，不是可部署方法；Projection、BandPass、R2 均通过 `Measured` 只接收测量列。数字由 `results/main_rows.json` 和 `results/condition_rows.json` 原始 rows 复算。','','## 实验配置','',f'- 主评估：base=380000，每个 CELLS 格每方法 15 个配对场景；条件扫描：base=390000，每条件每格 8 个。','- C_neutral：f_exc=20.8、a_force=FINAL，c/rel_max/chi2_max=FINAL；其余公共控制参数 61 次调参。','- C_projection：Projection(16,3) 来源的同预算公共调参；估计器自身均为 30 次（默认+20 随机+9 局部）。','- 误差指标中的真值刚度、真力仅用于评分；计时来自控制器调用时间。','','## 主评估逐工况汇总','']
    def rows(cg,m): return [r for r in main if r['config']==cg and r['method']==m]
    for cg in ('C_neutral','C_projection'):
        lines += [f'### {cg}','', '| 方法 | RMSE 中位/均值 | 失控 | 长失联 | 事件 RMSE | 纹波 | 调用 P50/P99 (μs) | ln(k̂/k) P50/P90 | 超 ln1.2 | 接受率 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
        for m in METHODS:
            rr=rows(cg,m); q=lambda k:[x['metrics'].get(k,np.nan) for x in rr]; rm=q('rmse'); ev=q('event_rmse'); rip=q('ripple'); p50=q('call_p50'); p99=q('call_p99'); e50=q('est_p50'); e90=q('est_p90'); ov=q('est_over_ln12'); ar=q('accept_rate');
            lines.append(f"| {m} | {np.nanmedian(rm):.4f}/{np.nanmean(rm):.4f} | {sum(x['metrics']['failed'] or x['metrics']['peak']>50 for x in rr)}/{len(rr)} | {sum(x['metrics']['long_loss'] for x in rr)}/{len(rr)} | {np.nanmedian(ev):.4f} | {np.nanmedian(rip):.4f} | {np.nanmedian(p50):.2f}/{np.nanmedian(p99):.2f} | {np.nanmedian(e50):.3f}/{np.nanmedian(e90):.3f} | {np.nanmean(ov):.3f} | {np.nanmean(ar):.3f} |")
        for m in ('BandPass','R2'):
            p=rows(cg,'Projection'); q=rows(cg,m); by={(x['cond'],x['tissue'],x['i']):x['metrics']['J'] for x in p}; diff=[]
            for x in q:
                k=(x['cond'],x['tissue'],x['i']);
                if k in by: diff.append(math.log1p(min(x['metrics']['J'],1000))-math.log1p(min(by[k],1000)))
            lines.append(f'- {m}−Projection 配对 log(1+min(J,1000))：均值 {np.mean(diff):+.4f}，bootstrap 95% CI [{bootstrap(diff)[0]:+.4f}, {bootstrap(diff)[1]:+.4f}]，n={len(diff)}。')
        lines.append('')
    if cond:
        lines += ['## 条件扫描（C_neutral）','', '| 条件 | 方法 | RMSE 中位 | 失控/总数 | 长失联/总数 |','|---|---|---:|---:|---:|']
        for scn in COND_SCANS:
            for m in ('Projection','BandPass','R2'):
                rr=[r for r in cond if r['scan']==scn and r['method']==m]; lines.append(f"| {scn} | {m} | {np.nanmedian([x['metrics']['rmse'] for x in rr]):.4f} | {sum(x['metrics']['failed'] or x['metrics']['peak']>50 for x in rr)}/{len(rr)} | {sum(x['metrics']['long_loss'] for x in rr)}/{len(rr)} |")
        lines.append('')
    lines += ['## 离线归因与主线判断','', '离线结果见 `results/filter_selection.json`、`results/failure_replay_rows.json`；它们只对同一 V1 轨迹上的估计误差做归因，不能替代闭环结论。主线判断必须以本节主评估 rows 的配对差异为准：若 BandPass−Projection 的区间跨 0，则报告为未分出；若在多个工况同向且区间不跨 0，才称为支持该方法。','', '### 复算校验','', '- 本报告的 RMSE、失控、长失联、事件窗口、纹波、调用分位数和估计误差均从原始 rows 即时聚合；重新运行本脚本应得到相同的行数与同一配对键 `(config, cond, tissue, i)`。', '- 调参结果保存在 `results/tuning.json`，包括全部候选、预算和最优点；边界命中需检查 `best.params` 与定义范围，不能仅凭最优分数判断。', '', '### 当前可证据化回答', '', '1. BandPass 与 Projection 的主线选择由上表和配对 bootstrap 区间决定；不把 R2 的结果外推为 Rossi 全套控制器。', '2. R2−BandPass 的差异直接回答旧二阶峰值滤波器是否造成“带通不行”，但仍只适用于本控制器和本仿真任务。', '3. 离线与闭环若方向不一致，优先解释为闭环估计误差经过 k_bar/k_good 门控、前馈和恢复放大后的尾部效应，而不是把离线排名当作闭环排名。', '4. 真值参考是上界式诊断，不是可部署算法；异常结果必须回到具体场景 rows 复跑，不能只看 pooled 均值。', '', '### 证据边界', '', '以上为仿真测量；“支持”是对本协议场景的实测统计描述，“机制解释”和主线推荐属于推断。']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def smoke():
    sc=cell_scene('breath','mid',0,345000)
    for m,p in (('Projection',dict(nw=16,degree=3)),('BandPass',dict(order=2,half=4.,memory=.4)),('R2',dict(Q=3.,memory=.4))):
        log,pol=_run(sc,m,FINAL_C,p,345000); print(m, None if log is None else len(log), _metric(log,pol,sc,20.8))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('phase',choices=['smoke','tune','expand','main','conditions','report','all']); ap.add_argument('--n',type=int,default=None); a=ap.parse_args()
    if a.phase=='smoke': smoke(); return
    if a.phase=='tune': do_tuning(); return
    tuning=json.loads((HERE/'results'/'tuning.json').read_text())
    if a.phase=='expand': expand_boundary_tuning(tuning); return
    if a.phase=='main': run_main(tuning,a.n or 15); return
    if a.phase=='conditions': run_conditions(tuning,a.n or 8); return
    if a.phase=='report': report(tuning); return
    if a.phase=='all':
        tuning=do_tuning(); run_main(tuning,a.n or 15); run_conditions(tuning,a.n or 8); report(tuning)


if __name__=='__main__': main()
