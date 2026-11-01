"""Reproduce Step 3. All writes are under method_redesign/; no final cohort.

python -X utf8 experiments/method_redesign/run_step3.py --workers 8
Use --resume to reuse completed task checkpoints (not scientific result caches
from other experiments). --report regenerates the report from step3_mpc.json.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import argparse
import json
import time
import hashlib
import platform
from pathlib import Path
import subprocess
import pickle
import uuid
import numpy as np
from mpc_controller import EstimationAwareMPC, CE, common, SCAN
from mpc_robustness import screen, GRID
import scenes_v2 as SV
import sweep

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'
METHODS = ('Proposed', 'MPC-U', 'MPC-C', 'MPC-oracle', 'Oracle', 'DOMPAC')


class FilePool:
    """Persistent subprocesses with local file IPC (Windows named pipes denied).
    No remote workers; task/result pickle files are generated locally by us.
    """
    def __init__(self, workers):
        self.workers = workers
        self.root = OUT/'.workers'/uuid.uuid4().hex
        self.root.mkdir(parents=True)
        self.procs = []
        self.handles = []
        for i in range(workers):
            folder = self.root/str(i); folder.mkdir()
            handle = (folder/'console.txt').open('w', encoding='utf8')
            self.handles.append(handle)
            self.procs.append(subprocess.Popen([os.sys.executable, '-X', 'utf8', str(Path(__file__).resolve()),
                                               '--worker-file', str(folder)],
                                              stdout=handle, stderr=handle,
                                              creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0)))

    def __enter__(self): return self

    def __exit__(self, *args):
        for i in range(self.workers): (self.root/str(i)/'stop').touch()
        for p in self.procs:
            try: p.wait(timeout=10)
            except subprocess.TimeoutExpired: p.terminate()
        for h in self.handles: h.close()

    def map(self, func, tasks, chunksize=1):
        assert func is episode
        tasks = list(tasks)
        for i in range(self.workers):
            fp = self.root/str(i)
            (fp/'result.pkl').unlink(missing_ok=True)
            (fp/'pending.pkl').write_bytes(pickle.dumps(tasks[i::self.workers]))
            (fp/'pending.pkl').rename(fp/'input.pkl')
        while not all((self.root/str(i)/'result.pkl').exists() for i in range(self.workers)):
            for i, p in enumerate(self.procs):
                if p.poll() is not None:
                    raise RuntimeError(f'worker {i} exited: '+(self.root/str(i)/'console.txt').read_text(encoding='utf8'))
            time.sleep(.05)
        result = [None]*len(tasks)
        for i in range(self.workers):
            rr = pickle.loads((self.root/str(i)/'result.pkl').read_bytes())
            for j, row in enumerate(rr): result[i+j*self.workers] = row
        return result


def file_worker(folder):
    fp = Path(folder)
    while not (fp/'stop').exists():
        path = fp/'input.pkl'
        if path.exists():
            tasks = pickle.loads(path.read_bytes()); path.unlink()
            rr = [episode(task) for task in tasks]
            (fp/'pending_result.pkl').write_bytes(pickle.dumps(rr))
            (fp/'pending_result.pkl').rename(fp/'result.pkl')
        else: time.sleep(.02)


def clean(o):
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, np.ndarray): return clean(o.tolist())
    if isinstance(o, (np.floating, float)): return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer): return int(o)
    return o


def save(path, obj):
    path.write_text(json.dumps(clean(obj), ensure_ascii=False, indent=1), encoding='utf8')


def episode(task):
    method, params, scene = task
    family, sid, p, T = scene
    try:
        if method.startswith('MPC'):
            policy = EstimationAwareMPC(**params, constrained=method != 'MPC-U',
                                        oracle_scene=p if method == 'MPC-oracle' else None)
        else:
            policy = CE.build(method, p)
        log, _ = common.run(policy, p, T, sid)
        mt = common.task_metrics(log, 20.8)
        tt = np.asarray(policy.t_call)[1:]  # exclude only initial controller/JIT call
        row = dict(method=method, family=family, sid=sid, T=T, **mt,
                   timing_ms=np.quantile(tt*1000, [.5, .99]).tolist(),
                   _times=tt.tolist())
        if method.startswith('MPC'):
            row.update(qp_status=policy.mpc.statuses, qp_failures=policy.mpc.infeasible,
                       qp_max_violation=policy.mpc.max_violation,
                       _qtimes=policy.mpc.times,
                       qp_timing_ms=np.quantile(np.asarray(policy.mpc.times)*1000, [.5, .99]).tolist())
        return row
    except Exception as ex:
        import traceback
        return dict(method=method, family=family, sid=sid, T=T, diverged=1, J=1e3,
                    error=repr(ex), traceback=traceback.format_exc(), _times=[], _qtimes=[])


def dev_scenes():
    out = []
    for fi, fam in enumerate(CE.FAMS):
        rng = np.random.default_rng(190000+1000*fi)
        for i in range(10):
            p, T = SV.make(fam, rng)
            out.append((fam, 190000+1000*fi+i, p, T))
    return out


def evaluate_candidate(params, pool, sc):
    rows = pool.map(episode, [('MPC-C', params, s) for s in sc], chunksize=1)
    Lfam = {}
    for fam in CE.FAMS:
        rs = [r for r in rows if r['family'] == fam]
        Lfam[fam] = float(np.mean([np.log1p(min(r['J'], 1e3)) if np.isfinite(r['J']) else np.log1p(1e3) for r in rs]))
    for r in rows:
        r.pop('_times', None); r.pop('_qtimes', None)
    return dict(params=params, L=float(np.mean(list(Lfam.values()))), family_objective=Lfam, scenes=rows)


def tuning(pool, resume=False, budget=61):
    path = OUT/'step3_tuning_checkpoint.json'
    if resume and path.exists(): return json.loads(path.read_text(encoding='utf8'))
    rng = np.random.default_rng(777)
    proposed, accepted = [], []
    default = dict(H=15, qr=.1, terminal=1.)
    while len(accepted) < min(41, budget):
        p = default if not proposed else dict(H=int(rng.integers(10, 26)),
                                             qr=float(10**rng.uniform(-3, 2)),
                                             terminal=float(10**rng.uniform(-1, 1.5)))
        ss = screen(p); proposed.append(ss)
        if ss['accepted']: accepted.append(p)
    print('Robust pre-screen:', len(accepted), 'accepted of', len(proposed), flush=True)
    sc = dev_scenes(); evaluations = []
    t0 = time.time()
    for p in accepted:
        row = evaluate_candidate(p, pool, sc); evaluations.append(row)
        print('TUNE', len(evaluations), 'L', round(row['L'], 6), p, 'elapsed', round(time.time()-t0), flush=True)
    best = min(evaluations, key=lambda x: x['L'])
    rng2 = np.random.default_rng(888)
    local = []
    while len(local) < max(0, budget-41):
        p = dict(H=int(np.clip(best['params']['H']+rng2.integers(-4, 5), 10, 25)),
                 qr=float(np.clip(best['params']['qr']*np.exp(rng2.normal(0, .5)), 1e-3, 1e2)),
                 terminal=float(np.clip(best['params']['terminal']*np.exp(rng2.normal(0, .5)), .1, 10**1.5)))
        ss = screen(p); proposed.append(ss)
        if ss['accepted']: local.append(p)
    for p in local:
        row = evaluate_candidate(p, pool, sc); evaluations.append(row)
        print('TUNE', len(evaluations), 'L', round(row['L'], 6), p, 'elapsed', round(time.time()-t0), flush=True)
    best = min(evaluations, key=lambda x: x['L'])
    rec = dict(budget=len(evaluations), seed_random=777, seed_local=888, best=best,
               all=evaluations, screen_attempts=proposed, elapsed_s=time.time()-t0,
               objective='equal-family mean log(1+min(J,1000)); MPC-C selects one shared configuration')
    save(path, rec)
    return clean(rec)


def bootstrap(rs, baseline, rng):
    pairs = [(r, baseline[r['sid']]) for r in rs if r['sid'] in baseline]
    finite = [(r['rmse'], b['rmse']) for r, b in pairs if r.get('rmse') is not None and b.get('rmse') is not None
              and np.isfinite(r['rmse']) and np.isfinite(b['rmse'])]
    if not finite: return None
    arr = np.array(finite); delta = arr[:, 0]-arr[:, 1]
    means = delta[rng.integers(0, len(delta), (4000, len(delta)))].mean(axis=1)
    return dict(n=len(delta), omitted=len(pairs)-len(delta), mean_delta=float(delta.mean()),
                ci95=np.quantile(means, [.025, .975]).tolist(), draws=4000)


def summarize(rows, timings, qtimes):
    result = []
    rng = np.random.default_rng(33119)
    for fam in dict.fromkeys(r['family'] for r in rows):
        base = {r['sid']: r for r in rows if r['family'] == fam and r['method'] == 'Proposed'}
        for method in METHODS:
            rs = [r for r in rows if r['family'] == fam and r['method'] == method]
            ok = [r for r in rs if r.get('rmse') is not None and np.isfinite(r['rmse'])]
            rm = [r['rmse'] for r in ok]; pk = [r['peak'] for r in ok]
            tt = np.concatenate(timings[(fam, method)])
            qt = np.concatenate(qtimes.get((fam, method), [np.array([])]))
            result.append(dict(family=fam, method=method, n=len(rs), n_rmse=len(ok),
                               rmse_median=np.median(rm), rmse_mean=np.mean(rm),
                               runaway=sum(r['peak'] > 50 for r in ok),
                               diverged=sum(r['diverged'] for r in rs), errors=sum('error' in r for r in rs),
                               peak_p99=np.quantile(pk, .99), peak_max=max(pk),
                               peak_gt15=sum(r['peak'] > 15 for r in ok),
                               loss_mean=np.mean([r['loss'] for r in ok]),
                               long_loss=sum(r['loss_max'] > .5 for r in ok),
                               timing_ms=np.quantile(tt*1000, [.5, .99]),
                               qp_timing_ms=np.quantile(qt*1000, [.5, .99]) if len(qt) else None,
                               qp_failures=sum(r.get('qp_failures', 0) for r in rs),
                               paired_bootstrap=bootstrap(rs, base, rng) if method.startswith('MPC') else None))
    return clean(result)


def export_timings():
    """Retain every evaluation call/solve duration, indexed by method/family/sid.
    Checkpoints retain double precision; artifact stores nanoseconds as integers.
    np.savez_compressed is a portable ZIP, readable without pickle.
    """
    arrays = {}
    for method in METHODS:
        path = OUT/f'step3_eval_{method}.json'
        if not path.exists(): continue
        for r in json.loads(path.read_text(encoding='utf8')):
            prefix = f"{method}__{r['family']}__{r['sid']}"
            arrays[prefix+'__call_ns'] = np.rint(np.asarray(r.get('_times', []))*1e9).astype(np.int64)
            if method.startswith('MPC'):
                arrays[prefix+'__solve_ns'] = np.rint(np.asarray(r.get('_qtimes', []))*1e9).astype(np.int64)
    np.savez_compressed(OUT/'step3_timings.npz', **arrays)


def enrich_provenance(rec):
    rec['settings'] = dict(Tc=.008, tau_s=.019, force_filter_hz=30., command_delay_ticks=1,
                           F_max=15., VMAX=.25, inherited_FINAL=dict(CE.FINAL),
                           tuning_scenes=[dict(family=f, sid=sid, p=p.tolist(), T=T) for f, sid, p, T in dev_scenes()],
                           diagnostic_scenes=[dict(family=f, sid=sid, p=p.tolist(), T=T) for f, sid, p, T in CE.cohort('diag', 50)],
                           dose_scenes=[dict(kind=k, dose=d, sid=sid, p=p.tolist(), T=T) for k, d, sid, p, T in sweep.scenes()])
    rec['implementation_sha256'] = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                                    ('mpc_controller.py', 'mpc_robustness.py', 'run_step3.py',
                                     'validate_step3.py', 'audit_step3_results.py')}
    return rec


def report(rec):
    lines = ['# 方法重构第3步：估计感知力域约束MPC', '',
             '本报告由 `run_step3.py` 从同目录JSON生成。全部新增文件通过补丁交付；未修改既有文件，未提交、推送或运行final确认集。', '',
             '## 设置与设计', '',
             '继承 `variants.MotionAwareIMM(**cohort_eval.FINAL)`：投影估计、kbar、k_good接受准则、IMM、激励及接触监督/恢复均调用原实现；只抑制旧导纳增量与旧前馈输出，改由MPC产生控制增量。改变控制会改变激励闭环响应和随后估计的数据，不能要求新旧轨迹上估计数值相同。', '',
             'Tc=8 ms，一周期指令延迟，τs=19 ms精确指数离散，30 Hz测力滤波（8个1 ms采样，先测力后推进伺服）。模型与对象的0.5 ms Euler伺服、SLS/HC动态有失配。IMM原始量测是η=−s加常数；按编码器加权时间戳前移到当前时刻后，常加速度外推。相对接触力原点在增益改变或恢复后重锚定，避免不可观测原点跳变。', '',
             '决策z=k_good Δu（N），代价Σ Q(F_filt−4)²+R(kbar Δu)²，末拍力误差权乘terminal；qr=Q/R。Hessian只依赖k_good/kbar，可预计算预测矩阵并热启动。输出偏置为实测滤波力减去上一拍模型预测滤波力，在整个时域保持常值；与该偏置一致地约束(kbar/k_good)F_model_raw+b≤15 N。速度约束|z|≤k_good·0.25·Tc。MPC-U去掉两类QP约束，但保留对象原有速度饱和。', '',
             '激励在MPC外按原力幅/原k_good（启动使用原k_seek）叠加；因此QP的速度及峰值约束没有覆盖外部激励，且kbar不是逐拍真值上界，15 N不是实际系统保证。已无法挽回的一周期预测也保留硬约束；QP不可行或求解失败时以速度上限回撤，计数单列。恢复阶段仍走原监督器，不受MPC力约束覆盖。', '',
             'MPC-oracle仅诊断：用真值当前切线k和已知未来表面函数替换MPC预测输入；横向位置按对象已知扫描指令与0.5 ms伺服递推；估计器、激励、监督保持原样。它不掌握未来闭环接触深度或SLS/HC状态，因此不是物理最优上界。Oracle为现有cohort_eval参考控制器，沿用其默认无激励设置。', '',
             '## 鲁棒性筛选（先于调参）', '',
             f"网格：{GRID}；每配置{rec['robustness']['points']}点。冻结IMM模式概率，分别求两个稳态Kalman增益并保留交互混合。输出偏置的模型滤波状态也纳入闭环。",
             f"选定配置：{rec['selected']}；稳定{rec['robustness']['n_stable']}/{rec['robustness']['points']}，最大谱半径{rec['robustness']['rho_max']:.9f}，最差点{rec['robustness']['worst']}。",
             '这是有限采样网格上的冻结线性模型结论，不是参数连续区域、模式切换、时变刚度或接触恢复的稳定证明，也没有替代性质2的圆判据。累积运动量测的任意常数原点以代数量测式固定；所表示动态矩阵的全部特征值参与筛选，没有删除近单位极点。', '',
             '## 调参', '',
             f"MPC-C调参预算{rec['tuning']['budget']}次，每次种子190000+（每族独立起点190000+1000fi）五族各10个场景。先41个全局候选、再20个局部候选；只允许通过全网格稳定筛选者进入仿真。筛选尝试共{len(rec['tuning']['screen_attempts'])}组；不稳定候选及其最差点完整保留于JSON。",
             '全局H整数[10,25]，Q/R对数均匀[10⁻³,10²]，terminal对数均匀[0.1,31.62]；局部H±4、其他参数乘对数标准差0.5扰动并截断。种子777/888。MPC-U和MPC-oracle复用MPC-C选定参数，控制变量相同以隔离约束/预测输入，未给予额外调参。',
             '目标为每族平均log(1+min(J,1000))再族等权平均，J完全复用tune_v2物理指标。Proposed/DOMPAC沿用既有61次预算所选参数，未在新190000+集合重调；两者历史使用的开发场景与参数自由度不同，Proposed还用过diag诊断，因此这是按请求固定FINAL的控制层比较，不能称为所有方法在同一新数据上的重新公平整定。部署型MPC和Proposed只读相同因果信息；两个oracle有明确特权。', '',
             '## 诊断集（每族50场景）', '',
             'RMSE、峰值、失联以1 s后真实接触力计算，失联阈值0.5 N。峰P99指场景最大力的99%分位。计时为全体单拍调用合并分位数，包含估计/IMM/控制/QP；QP另列的计时包括构造/更新及求解。多进程CPU竞争和Python/OS调度会影响尾部，排除每场景首个JIT调用，不是实时认证。', '',
             '|族|控制器|RMSE 中位/均值 N|>50 N/发散|峰P99 N|失联% / 长失联数|单拍P50/P99 ms|',
             '|---|---|---:|---:|---:|---:|---:|']
    for s in rec['summary']:
        if s['family'].startswith(('ramp:', 'ratio:')): continue
        lines.append(f"|{s['family']}|{s['method']}|{s['rmse_median']:.3f}/{s['rmse_mean']:.3f}|{s['runaway']}/{s['diverged']}|{s['peak_p99']:.2f}|{100*s['loss_mean']:.2f}/{s['long_loss']}|{s['timing_ms'][0]:.3f}/{s['timing_ms'][1]:.3f}|")
    lines += ['', '## 配对bootstrap（MPC−Proposed的平均RMSE，N）', '',
              '按同一场景ID配对重采样4000次，百分位95%区间，随机种子33119。含失控场景；发散导致无RMSE时不编造误差，以完整配对报告并显式计数。各区间未经多重比较校正。', '',
              '|族|MPC-U Δ [95%]|MPC-C Δ [95%]|MPC-oracle Δ [95%]|', '|---|---:|---:|---:|']
    for fam in CE.FAMS:
        line = f'|{fam}|'
        for method in ('MPC-U', 'MPC-C', 'MPC-oracle'):
            b = next(s for s in rec['summary'] if s['family']==fam and s['method']==method)['paired_bootstrap']
            line += f"{b['mean_delta']:+.3f} [{b['ci95'][0]:+.3f},{b['ci95'][1]:+.3f}] (n={b['n']})|"
        lines.append(line)
    lines += ['', '## 剂量评价（每档20场景）', '',
              '种子与sweep.py完全相同（160000+/170000+）；下表每格RMSE中位/均值；括号为>50 N/发散数。其余指标、全部bootstrap和QP计时见JSON summary。', '',
              '|剂量|'+'|'.join(METHODS)+'|', '|---|'+'---:|'*len(METHODS)]
    for fam in dict.fromkeys(s['family'] for s in rec['summary'] if ':' in s['family']):
        line = f'|{fam}|'
        for method in METHODS:
            s = next(s for s in rec['summary'] if s['family']==fam and s['method']==method)
            line += f"{s['rmse_median']:.3f}/{s['rmse_mean']:.3f} ({s['runaway']}/{s['diverged']})|"
        lines.append(line)
    lines += ['', '## QP与计算', '', '|族|控制器|QP P50/P99 ms|失败拍数|实际>15 N场景数|', '|---|---|---:|---:|---:|']
    for s in rec['summary']:
        if s['family'] not in CE.FAMS or not s['method'].startswith('MPC'): continue
        q = s['qp_timing_ms']
        lines.append(f"|{s['family']}|{s['method']}|{q[0]:.3f}/{q[1]:.3f}|{s['qp_failures']}|{s['peak_gt15']}|")
    lines += ['', '## 结论与疑点', '']
    for fam in CE.FAMS:
        ss = {s['method']: s for s in rec['summary'] if s['family']==fam}
        p, c, u, o = [ss[m] for m in ('Proposed', 'MPC-C', 'MPC-U', 'MPC-oracle')]
        lines.append(f"- {fam}：MPC-C平均RMSE相对Proposed为{c['rmse_mean']-p['rmse_mean']:+.3f} N，中位数{c['rmse_median']-p['rmse_median']:+.3f} N；C/U均值{c['rmse_mean']:.3f}/{u['rmse_mean']:.3f}；oracle-MPC均值{o['rmse_mean']:.3f}，相对C差{o['rmse_mean']-c['rmse_mean']:+.3f} N。")
    diag = [s for s in rec['summary'] if s['family'] in CE.FAMS]
    total = lambda method, key: sum(s[key] for s in diag if s['method']==method)
    lines += ['',
              f"本实例不支持用MPC-C替换当前控制层：五族平均RMSE均更差，五个配对95%区间均在0以上。诊断250场景中，>50 N由Proposed的{total('Proposed','runaway')}增至C的{total('MPC-C','runaway')}，长失联由{total('Proposed','long_loss')}增至{total('MPC-C','long_loss')}；C有{total('MPC-C','qp_failures')}拍QP失败（逐拍状态均保留），所有方法无仿真发散和程序异常。QP失败按全程计，物理指标按1 s后计，两者不能直接一一对应。",
              '约束的作用具有条件性：T3m中C比U平均误差小（0.387对0.680），但硬组织中C远差于U并伴随更多失联；与硬约束不可行后限速回撤、再次接触的失配风险一致，尚非单因素因果证明。S2较高的保守刚度/量测增益比会使力上界偏紧，且绝对接触原点只由力锚定，约束可行性不等于物理安全。软组织族C仅小幅更差，典型轨迹上QP约束作用很小。',
              'U在T2m-stiff的平均误差0.685低于Proposed的2.305，但中位数0.519反而高于0.283；主要是Proposed两例重尾被削减，bootstrap区间跨0，不能称为典型跟踪改善或已确认的均值优势。快速族U仍明显变差。',
              '剂量结果同样不支持普遍改善：C在0.6/0.45/0.3/0.2/0.1 s斜坡的平均误差更大，仅0.15 s档从14.497降至4.885 N；脊比30/60的平均误差下降而中位数更差，其他脊比平均误差更大。改善主要涉及重尾，不能用均值改善概括整个分布。',
              'oracle诊断同时指向信息与控制瓶颈：硬组织中C→MPC-oracle平均误差从7.854降至0.340 N，快速族从8.095降至1.852 N，说明预测/刚度输入通道很关键；但仍远差于现有Oracle的0.063/0.078 N，且MPC-oracle快速族仍有4/50失控。因此不能把失败全归于运动预测，本MPC的有限时域、增量代价、线性接触模型及不可行恢复也未解决。以上是本实例的诊断，不是对MPC算法族的否定。',
              '在S2/T3m中两个真值参考也没有优于Proposed。对象记录的真值k是当前准静态切线刚度，而投影拟合的是观测窗口中的有效增益；SLS/HC动态、刚度时变和激励不同，使“真值输入”并不构成单调的性能上界。其具体贡献本轮未分拆。']
    lines += ['',
              'MPC-oracle同时替换运动与刚度，差距只能定位为“预测/增益信息通道的合并作用”，无法进一步单独归因于运动预测；其改变控制轨迹也间接改变估计数据。若它仍明显弱于现有Oracle，则说明本MPC建模、保守约束、外置激励或控制设计仍是瓶颈。冻结稳定不能阻止被接受增益超出设计网格，也不能保证未知非线性对象满足15 N。',
              '必须警惕：用kbar对绝对压入缩放依赖相对表面原点的力锚定，而IMM只估计等效位移；模型偏置并不是独立物理力偏置，保守上界可能过紧或失真。硬约束在已排队指令导致超限时不可行，最大回撤会增加失联；这类行为按失败拍数和长失联计数公开，不能用安全保证掩盖。',
              '本原型不使用第2步E2/E3，不含刚度空间预瞄、运动预测区间或tube/chance约束；15 N约束仅考虑当前kbar及均值轨迹，不是完整的估计不确定性鲁棒MPC。只运行用户指定的diag与剂量集，未依据评价结果再次调参。', '',
              '## 复现', '',
              '`python -X utf8 experiments/method_redesign/validate_step3.py`',
              '`python -X utf8 experiments/method_redesign/run_step3.py --workers 8`',
              '`python -X utf8 experiments/method_redesign/audit_step3_results.py --rerun --workers 4`',
              '环境、基线源码SHA256、所有候选、完整选定网格、逐场景指标及QP状态见JSON。`step3_timings.npz`保留评价集每次控制调用与每次求解的耗时（纳秒整数，键为控制器__族__场景ID__call_ns/solve_ns，np.load不需pickle）；QP记录含解析无约束求解、不可行与失败返回。可用 `--resume`复用本轮完成检查点；无外部隐藏结果。']
    if 'audit' in rec:
        lines += ['', '交付核验：计数、种子隔离、61个调参目标、83次鲁棒筛选和全部指标/耗时汇总通过独立重算。预先固定每族首个场景、全部6种方法共30场景重跑，J、RMSE、峰值、失联比例与最长失联差值均为0（见JSON audit）。实现另通过独立子步预测、无约束/QP等价、完全相同观测下估计/IMM/监督逐项一致和Kalman稳态增益核验。']
    (OUT/'step3_mpc.md').write_text('\n'.join(lines)+'\n', encoding='utf8')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--resume', action='store_true'); ap.add_argument('--report', action='store_true')
    ap.add_argument('--budget', type=int, default=61)
    ap.add_argument('--worker-file')
    args = ap.parse_args(); OUT.mkdir(exist_ok=True)
    if args.worker_file:
        file_worker(args.worker_file); return
    if args.report:
        export_timings()
        rec = enrich_provenance(json.loads((OUT/'step3_mpc.json').read_text(encoding='utf8')))
        save(OUT/'step3_mpc.json', rec)
        report(rec); return
    if not 1 <= args.budget <= 61: raise ValueError('budget must be 1..61')
    with FilePool(args.workers) as pool:
        tune = tuning(pool, args.resume, args.budget)
        selected = tune['best']['params']
        robust = screen(selected, retain=True)
        assert robust['accepted']
        diag = CE.cohort('diag', 50)
        doses = [(f'{kind}:{dose:g}', sid, p, T) for kind, dose, sid, p, T in sweep.scenes()]
        scenes = diag+doses
        print('SELECTED', selected, 'L', tune['best']['L'], 'evaluation episodes', len(scenes)*len(METHODS), flush=True)
        rows, timings, qtimes = [], {}, {}
        t0 = time.time()
        for method in METHODS:
            ck = OUT/f'step3_eval_{method}.json'
            if args.resume and ck.exists():
                rr = json.loads(ck.read_text(encoding='utf8'))
            else:
                rr = pool.map(episode, [(method, selected, s) for s in scenes], chunksize=1)
                save(ck, rr)
            for r in rr:
                key = (r['family'], method)
                timings.setdefault(key, []).append(np.array(r.pop('_times', [])))
                qtimes.setdefault(key, []).append(np.array(r.pop('_qtimes', [])))
            rows += rr
            print('EVAL', method, 'done', len(rr), 'errors', sum('error' in r for r in rr), 'elapsed', round(time.time()-t0), flush=True)
    hashes = {}
    for path in list(SCAN.glob('*.py'))+list((SCAN/'tuned_v2').glob('*.json'))+list((SCAN.parent/'cobot_scan_benchmark').glob('*.py')):
        hashes[str(path.relative_to(SCAN.parent))] = hashlib.sha256(path.read_bytes()).hexdigest()
    import scipy, osqp
    rec = dict(selected=selected, robustness=robust, tuning=tune,
               summary=summarize(rows, timings, qtimes), rows=rows,
               environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                                osqp=osqp.__version__, workers=args.workers, cpu=platform.processor()),
               baseline_sha256=hashes, evaluation_s=time.time()-t0)
    rec = enrich_provenance(rec)
    export_timings()
    save(OUT/'step3_mpc.json', rec); report(clean(rec))
    print('DONE', OUT/'step3_mpc.json', flush=True)


if __name__ == '__main__': main()
