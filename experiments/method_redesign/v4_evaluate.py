"""Reproduce V4 without writing outside method_redesign. No final cohort access.

python -X utf8 experiments/method_redesign/v4_evaluate.py --stage all --workers 12
Stages: check, tune, evaluate, timing, report, all. No per-tick files are saved.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, hashlib, json, math, platform, time, sys, tempfile, subprocess, queue, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from multiprocessing import get_context
import numpy as np
import v4_architecture as V
from v4_architecture import P, CE, common, variants, ContactUpdate, ServoID
import scenes_v2 as SV
import sweep

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'
METHODS = ('V4', 'V4-noFF', 'C1', 'Proposed', 'Proposed+V1')


def fs_retry(fn):
    for i in range(500):
        try:
            return fn()
        except PermissionError:
            if i == 499:
                raise
            time.sleep(.02)


class FilePool:
    """Persistent Python processes using file IPC (sandbox denies Windows pipes).

    Independent episodes run in processes; threads only wait for file replies.
    IPC and error files live in a temporary directory and are not artifacts.
    """
    def __init__(self, workers, initializer=None):
        self.workers = workers

    def __enter__(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='v4_ipc_')
        self.slots = queue.Queue(); self.processes = []; self.errors = []
        self.threads = ThreadPoolExecutor(self.workers)
        for i in range(self.workers):
            d = Path(self.tmp.name)/str(i); d.mkdir()
            err = (d/'stderr.txt').open('w',encoding='utf-8'); self.errors.append(err)
            proc = subprocess.Popen([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--worker-dir',str(d)],
                                    stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=err,
                                    creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            self.processes.append(proc); self.slots.put((d,proc))
        return self

    def call(self, job):
        d,proc = self.slots.get()
        try:
            tmp=d/'request.tmp'
            tmp.write_text(json.dumps(job,default=lambda x:x.tolist()),encoding='utf-8')
            fs_retry(lambda:tmp.replace(d/'request.json'))
            while not (d/'response.json').exists():
                if proc.poll() is not None:
                    raise RuntimeError((d/'stderr.txt').read_text(encoding='utf-8'))
                time.sleep(.03)
            reply=json.loads(fs_retry(lambda:(d/'response.json').read_text(encoding='utf-8')))
            fs_retry(lambda:(d/'response.json').unlink())
            if 'error' in reply:
                raise RuntimeError(reply['error'])
            return reply['row']
        finally:
            self.slots.put((d,proc))

    def map(self, fn, jobs, chunksize=1):
        return list(self.threads.map(self.call,jobs))

    def imap_unordered(self, fn, jobs, chunksize=1):
        futures=[self.threads.submit(self.call,j) for j in jobs]
        for f in as_completed(futures):
            yield f.result()

    def __exit__(self, *exc):
        self.threads.shutdown(wait=True)
        for d in Path(self.tmp.name).iterdir():
            (d/'stop').touch()
        for p in self.processes:
            p.wait()
        for f in self.errors:
            f.close()
        self.tmp.cleanup()


def worker_loop(d):
    d=Path(d)
    warmup()
    while not (d/'stop').exists():
        req=d/'request.json'
        if not req.exists():
            time.sleep(.03); continue
        try:
            job=json.loads(fs_retry(lambda:req.read_text(encoding='utf-8')))
            fs_retry(lambda:req.unlink())
            reply=dict(row=episode(job))
        except Exception:
            reply=dict(error=traceback.format_exc())
        tmp=d/'response.tmp'
        tmp.write_text(json.dumps(reply,allow_nan=False),encoding='utf-8')
        fs_retry(lambda:tmp.replace(d/'response.json'))


def save(name, data):
    OUT.mkdir(exist_ok=True)
    target = OUT/name
    tmp = target.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=1, allow_nan=False), encoding='utf-8')
    tmp.replace(target)


def build(name, params):
    if name == 'V4':
        return V.V4(**params)
    if name == 'V4-noFF':
        return V.V4(motion=False, ru=params.get('ru', 0.))
    if name == 'C1':
        return ContactUpdate()
    if name == 'Proposed':
        return variants.build('Proposed-IMM', dict(CE.FINAL))
    if name == 'Proposed+V1':
        return ServoID(single_model=False, **CE.FINAL)
    raise ValueError(name)


def train_scenes():
    out = []
    for fi, fam in enumerate(CE.FAMS):
        rng = np.random.default_rng(190000+1000*fi)
        for i in range(10):
            p, T = SV.make(fam, rng)
            out.append((fam, 190000+1000*fi+i, p, T))
    return out


def eval_scenes():
    out = [('diag', fam, sid, p, T) for fam, sid, p, T in CE.cohort('diag', 50)]
    out += [('sweep', f'{kind}:{dose:g}', sid, p, T) for kind, dose, sid, p, T in sweep.scenes()]
    assert all(not 130000 <= s[2] < 140000 for s in out)
    return out


def episode(job, keep_times=False):
    idx, name, params, scene = job
    split, fam, sid, p, T = scene
    pol = build(name, params)
    measured_times = []
    if keep_times:
        inner = pol
        class TimedPolicy:
            direct = getattr(inner, 'direct', False)
            def reset(self, platform):
                inner.reset(platform)
            def __call__(self, o):
                started = time.perf_counter()
                result = inner(o)
                measured_times.append(time.perf_counter()-started)
                return result
            def __getattr__(self, name):
                return getattr(inner, name)
        pol = TimedPolicy()
    log = P.episode(p, pol, T=T, seed=sid)
    m = common.task_metrics(log, 20.8)
    finite = not m['diverged'] and np.isfinite(m['rmse'])
    row = dict(index=idx, method=name, split=split, family=fam, seed=sid,
               rmse=float(m['rmse']) if finite else None,
               peak=float(m['peak']) if finite else None,
               runaway=bool(not finite or m['peak'] > 50),
               diverged=not finite,
               long_loss=bool(m['loss_max'] > .5) if finite else None,
               loss_max=float(m['loss_max']) if finite else None,
               L=math.log1p(min(float(m['J']), 1e3)) if finite else math.log1p(1e3))
    if hasattr(pol, 'counts'):
        row['diagnostic'] = dict(pol.counts, force_drift_max_N=pol.force_drift_max,
                                 tau_hat=pol.tau_hat, ndel_hat=pol.ndel_hat)
    if keep_times:
        return row, measured_times
    return row


def warmup():
    # JIT is warmed before timing; kernels cache only in the writable copy.
    sc = train_scenes()[0]
    for name in METHODS:
        episode((0, name, {}, ('warm', sc[0], sc[1], sc[2], 1.3)))


def check():
    import candidate1_model as C
    # The modified kernel's 16-state path must be numerically identical to C1.
    c = ContactUpdate(); c.reset(P.PLATFORMS['PC125'])
    x = np.array([math.log(200.), 1., .2, -.1])
    st = np.array([.001, .001, .2, .001, .001, 4., .001, 4., 7000., .002, .03, .001, 0., 0., 0., 0.])
    q = np.array([.001]); dith = np.zeros(c.H)
    prev = x.copy(); prior = x.copy()
    r0 = np.zeros(c.nr); r1 = r0.copy()
    C._resid(x, c.H, c.per, st, q, dith, c.pr, c.scen, prev, prior, c.rw, r0)
    V._resid(x, c.H, c.per, st, q, dith, c.pr, c.scen, prev, prior, c.rw, r1)
    np.testing.assert_allclose(r0, r1, rtol=0, atol=0)
    # Tabulated quadratic motion reproduces the polynomial rollout exactly
    # with k_good==k_model, ga=1 and tv beyond the full horizon.
    pr = c.pr.copy(); pr[8] = 1.; pr[9] = 10.
    h = (np.arange((c.H+1)*c.ns)+1)*.001
    tab = st[8]*(st[9]*h+.5*st[10]*h*h)
    # wv=0 in actual C1 defaults, so the nominal terminal residual is also identical.
    C._resid(x, c.H, c.per, st, q, dith, pr, c.scen, prev, prior, c.rw, r0)
    V._resid(x, c.H, c.per, np.r_[st, tab], q, dith, pr, c.scen, prev, prior, c.rw, r1)
    np.testing.assert_allclose(r0, r1, rtol=1e-11, atol=1e-11)
    # An identity actuator-ID hook makes noFF identical to original C1 closed loop.
    class NoID(V.V4):
        def _apply(self, tau, ndel):
            pass
    fam, sid, p, T = train_scenes()[0]
    a = P.episode(p, ContactUpdate(), T=1.5, seed=sid)
    b = P.episode(p, NoID(motion=False), T=1.5, seed=sid)
    np.testing.assert_allclose(a, b, rtol=0, atol=0)
    # Measured gain and recovery gates, no scene truth.
    c = V.V4(); c.reset(P.PLATFORMS['PC125'])
    c.observer.k_good = 8000.; c.observer.tracker.reset(0.)
    c.observer.tracker.x[:] = [0., .01, .2]; c.motion_ready = True
    c.kuse = 4000.
    out = c._motion_state(st.copy())
    expected = 8000.*(.01*h+.2*(.5*c.Tc*h+.5*h*h))
    np.testing.assert_allclose(out[16:], expected)
    c.motion_ready = False
    assert np.all(c._motion_state(st.copy())[16:] == 0)
    save('v4_checks.json', dict(passed=True, legacy_residual_exact=True,
         quadratic_motion_equivalent=True, noID_noFF_closed_loop_exact=True,
         measurement_gain_force_conversion=True, disabled_motion_zero=True,
         actual_C1_iters=c.iters, H=c.H, FINAL=CE.FINAL))
    print('CHECKS PASSED', flush=True)


def tune(workers):
    path = OUT/'v4_tuning.json'
    rec = json.loads(path.read_text(encoding='utf-8')) if path.exists() else dict(configs=[], trials=[])
    rng = np.random.default_rng(777)
    if not rec['configs']:
        rec['configs'] = [dict(motion_weight=1., ru=0.)]
        rec['configs'] += [dict(motion_weight=float(rng.uniform(.25, 1.25)),
                                ru=float(10**rng.uniform(-4, 0))) for _ in range(40)]
    sc = [('tune', f, sid, p, T) for f, sid, p, T in train_scenes()]
    with FilePool(workers, initializer=warmup) as pool:
        for stage in range(2):
            if stage == 1 and len(rec['configs']) == 41:
                best = min(rec['trials'], key=lambda r: r['objective'])['params']
                rng = np.random.default_rng(888)
                rec['configs'] += [dict(motion_weight=float(np.clip(best['motion_weight']*math.exp(rng.normal(0,.2)), .25, 1.25)),
                                       ru=float(np.clip(max(best['ru'],1e-4)*math.exp(rng.normal(0,.7)), 1e-4, 1.)))
                                   for _ in range(20)]
            target = 41 if stage == 0 else 61
            # One candidate at a time permits exact checkpoint/resume without
            # selecting on any evaluation seed. Parallelism is over scenes.
            while len(rec['trials']) < target:
                i = len(rec['trials']); params = rec['configs'][i]
                t0 = time.perf_counter()
                rows = pool.map(episode, [(i, 'V4', params, s) for s in sc], chunksize=1)
                family_L = {f: float(np.mean([r['L'] for r in rows if r['family']==f])) for f in CE.FAMS}
                rec['trials'].append(dict(index=i, params=params, objective=float(np.mean(list(family_L.values()))),
                                         family_L=family_L, runaways=sum(r['runaway'] for r in rows)))
                rec['best'] = min(rec['trials'], key=lambda r:r['objective'])
                rec['protocol'] = dict(budget=61, train_seed_base=190000, n_per_family=10,
                     objective='equal-family mean log(1+min(J,1e3))', random_seed=777, local_seed=888,
                     bounds=dict(motion_weight=[.25,1.25],ru=[1e-4,1.]),
                     baseline_note='C1 defaults and FINAL frozen by task; no baseline retuning')
                save('v4_tuning.json', rec)
                print(f'TUNE {i+1}/61 L={rec["trials"][-1]["objective"]:.5f} best={rec["best"]["objective"]:.5f} seconds={time.perf_counter()-t0:.1f}', flush=True)


def evaluate(workers):
    params = json.loads((OUT/'v4_tuning.json').read_text(encoding='utf-8'))['best']['params']
    path = OUT/'v4_episode_summary.json'
    rows = json.loads(path.read_text(encoding='utf-8')) if path.exists() else []
    done = {(r['method'],r['split'],r['family'],r['seed']) for r in rows}
    sc = eval_scenes()
    jobs = [(i,n,params,s) for i,n in enumerate(METHODS) for s in sc if (n,s[0],s[1],s[2]) not in done]
    with FilePool(workers, initializer=warmup) as pool:
        for i, row in enumerate(pool.imap_unordered(episode, jobs, chunksize=1)):
            rows.append(row)
            if (i+1)%50 == 0:
                save('v4_episode_summary.json', rows)
                print(f'EVAL {len(rows)}/{len(METHODS)*len(sc)}',flush=True)
    save('v4_episode_summary.json', sorted(rows,key=lambda r:(r['split'],r['family'],r['seed'],r['method'])))


def timing():
    params = json.loads((OUT/'v4_tuning.json').read_text(encoding='utf-8'))['best']['params']
    warmup()
    out = {}
    # Same five first diagnostic scenes, single process, no pool contention.
    sc = [s for s in eval_scenes() if s[0]=='diag' and s[2]%1000 == 0]
    for name in METHODS:
        times = []
        for s in sc:
            _, t = episode((0,name,params,s), keep_times=True)
            times.extend(t[1:])
        out[name] = dict(p50_ms=float(np.percentile(times,50)*1000),
                         p99_ms=float(np.percentile(times,99)*1000), ticks=len(times))
        print('TIMING',name,out[name],flush=True)
    save('v4_timing.json',dict(methods=out, protocol='Single process, 5 fixed diagnostic episodes, all ticks except first; includes estimator+IMM+ID+GN+recording; excludes JIT, plant, communication.',
         cpu=platform.processor(), python=platform.python_version(), logical_cpus=os.cpu_count()))


def summary(rows):
    finite = [r['rmse'] for r in rows if r['rmse'] is not None]
    return dict(n=len(rows), rmse_median=float(np.median(finite)) if finite else None,
                rmse_mean=float(np.mean(finite)) if finite else None,
                runaway=sum(r['runaway'] for r in rows), diverged=sum(r['diverged'] for r in rows),
                long_loss=sum(r['long_loss'] is True for r in rows), loss_unknown=sum(r['long_loss'] is None for r in rows))


def bootstrap(rows, ref):
    a = {r['seed']:r for r in rows if r['method']=='V4'}
    b = {r['seed']:r for r in rows if r['method']==ref}
    d = [a[s]['rmse']-b[s]['rmse'] for s in sorted(a.keys()&b.keys())
         if a[s]['rmse'] is not None and b[s]['rmse'] is not None]
    if not d:
        return dict(n_pairs=0, excluded_pairs=len(a), mean_difference=None, ci95=None)
    rng = np.random.default_rng(20261101)
    d = np.asarray(d)
    means = np.mean(d[rng.integers(0,len(d),size=(4000,len(d)))],axis=1)
    return dict(n_pairs=len(d), excluded_pairs=len(a)-len(d), mean_difference=float(d.mean()),
                ci95=[float(v) for v in np.percentile(means,[2.5,97.5])], resamples=4000)


def report():
    rows = json.loads((OUT/'v4_episode_summary.json').read_text(encoding='utf-8'))
    groups = sorted(set((r['split'],r['family']) for r in rows))
    tables = {}
    for split,fam in groups:
        rr = [r for r in rows if r['split']==split and r['family']==fam]
        tables[f'{split}/{fam}'] = dict(metrics={n:summary([r for r in rr if r['method']==n]) for n in METHODS},
                                      paired={n:bootstrap(rr,n) for n in ('Proposed','C1')})
    architecture_diagnostics = {}
    for fam in CE.FAMS:
        rr=[r for r in rows if r['split']=='diag' and r['family']==fam and r['method']=='V4']
        ds=[r['diagnostic'] for r in rr]
        architecture_diagnostics[fam]=dict(
            zero_motion_episodes=sum(d['motion_ticks']==0 for d in ds),
            recovery_over_half_episodes=sum(d['recovery_ticks']>d['ticks']/2 for d in ds),
            median_recovery_tick_fraction=float(np.median([d['recovery_ticks']/max(d['ticks'],1) for d in ds])),
            max_force_drift_N=max(d['force_drift_max_N'] for d in ds),
            drift_above_50N_episodes=sum(d['force_drift_max_N']>50 for d in ds))
    source_files = [HERE/'v4_architecture.py',HERE/'v4_evaluate.py', HERE.parent/'cobot_scan_benchmark/contact_update.py',
                    HERE.parent/'cobot_scan_benchmark/candidate1_model.py', HERE/'v1_servo_id.py',
                    HERE.parent/'motion_aware_scan/variants.py',HERE.parent/'motion_aware_scan/proposed.py',
                    HERE.parent/'motion_aware_scan/imm.py',HERE/'v4_failure_analysis.py',
                    HERE.parent/'motion_aware_scan/common.py',HERE.parent/'motion_aware_scan/cohort_eval.py',
                    HERE.parent/'motion_aware_scan/scenes_v2.py',HERE.parent/'motion_aware_scan/sweep.py',
                    HERE.parent/'motion_aware_scan/tuned_v2/Proposed-IMMb.json',
                    HERE.parent/'cobot_scan_benchmark/plant.py',HERE.parent/'cobot_scan_benchmark/scenarios.py']
    result = dict(design='C1 GN + FINAL force-equivalent motion/IMM + V1 online servo ID',
        seed_policy='tune 190000+, diag 110000+, sweep 160000+/170000+; no final 130000+ calls',
        implementation=dict(stiffness='C1 repaired projection unchanged; a separate FINAL projection fit supplies exactly the existing k_good acceptance and motion measurement',
             drift='nominal rollout only: delta_f(h)=weight*k_good*(v*h+a*(0.5*Tc*h+0.5*h*h)); 1 ms table, full horizon, no legacy tv clipping; no extra tau*sd because servo explicit',
             recovery='C1 unchanged, IMM cleared and disabled while in C1 recovery/F<0.3; no additive external FF; C1 last-motion recovery rule retained',
             solver='Original source defaults: H=16, iters=1 (not the 3 claimed in prose); kept identical for original C1 and V4 per task preserve-settings instruction',
             variables='log D, KPF/(1e-4), uff_0/(1e-3), uff_slope/(1e-2); original warm start, ramp-then-hold, probes, residuals, bounds',
             ablation='V4-noFF retains C1 polynomial drift, adds V1 only, and shares the selected V4 ru to avoid confounding the motion ablation with a different control-rate penalty; no motion observer work',
             timing='Separate single-process pass, no physical real-time claim',
             information='only t, encoder y, filtered F, command; no Ft/My or contact/scene truth; existing gain and direct interfaces retained'),
        tuning=json.loads((OUT/'v4_tuning.json').read_text(encoding='utf-8')),
        timing=json.loads((OUT/'v4_timing.json').read_text(encoding='utf-8')) if (OUT/'v4_timing.json').exists() else None,
        checks=json.loads((OUT/'v4_checks.json').read_text(encoding='utf-8')),
        failure_analysis=json.loads((OUT/'v4_failure_analysis.json').read_text(encoding='utf-8')) if (OUT/'v4_failure_analysis.json').exists() else None,
        architecture_diagnostics=architecture_diagnostics,
        baseline_fact_recheck_first20={f:{n:summary([r for r in rows if r['split']=='diag' and r['family']==f and r['method']==n and r['seed']%1000<20]) for n in ('C1','Proposed')} for f in CE.FAMS},
        execution_notes=['Windows multiprocessing Pipe creation denied; persistent process workers use temporary file IPC.',
                         'At 39 saved trials an IPC PermissionError interrupted the next candidate; that same candidate was repeated after adding file-access retries. Exactly 61 distinct candidate configurations received completed objective scores; one interrupted attempt incurred extra simulation work, with no additional parameter selection opportunity.',
                         'All metric computation uses the existing common.task_metrics warmup=1 s and true-force scoring convention.'],
        tables=tables, source_sha256={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        limitations=['Prespecified historical baseline parameters frozen, not independently retuned 61 times on new seeds. V4 alone gets <=61 new two-parameter trials as requested.',
                     'Diagnostic and sweep seeds are historically used: no fresh-confirmation claim.',
                     'Divergent RMSE undefined; finite-only moments and complete-finite-pair bootstrap explicitly report exclusions; divergence separately counted.',
                     'Gain policy observes the prior logged row; direct policies consume the provisional current row as in their original implementations. Same sensor types, but existing timing/interface differences remain.',
                     'Nominal plant evaluation only; V1 mismatch validation is inherited, not newly established here.'])
    result['conclusion'] = dict(
        adopted=False,
        verdict='未同时保住C1在S2的优势与Proposed在运动族的优势；不采纳当前V4拼接实例。',
        supported_findings=[
            'S2：V4均值0.159 N，较C1 0.089 N差0.069 N，配对95% CI [0.056,0.083]。V4较Proposed均值0.281 N更好，但中位数0.148 N差于Proposed 0.124 N；不能只报均值称兼得。',
            'T3m、软推拉、硬推拉：V4均值分别较Proposed差0.894、0.074、3.647 N，配对区间均高于0。快速硬推拉差7.789 N，95% CI [6.031,9.795]；失控22/50、长失联45/50。',
            'V4-noFF的S2为0.083/0.089 N，与C1基本一致；V4变差出现在接入运动预测后的闭环链条，而非仅增加V1辨识。',
            '剂量扫描：V4在ratio=3时优于C1；ratio>=10时没有保住C1精度。运动扫描各档中位数均大幅差于Proposed；ramp=0.15时Proposed均值14.497 N反而高于V4的10.542 N，说明其尾部也严重，不能宣称逐档均值全胜。'],
        failure_evidence=[
            '恢复/启动阻塞（直接观测）：硬推拉8/50、快速推拉5/50从未启用运动预测；29/50和40/50超过半数拍处于恢复。场景113026在0.064 s进入恢复，991/999记录拍处于恢复，运动预测从未启用。失联关闭FF的规则与原C1恢复组合后，运动信息无法帮助这些场景。',
            '长时域预测失真（直接观测）：T3m、硬推拉、快速推拉的名义力漂移预测绝对值最大6481、5032、3486 N，分别12/50、37/50、36/50场景出现>50 N预测漂移。它们是模型预测数值，不是真实接触峰值。',
            '离线已实现轨迹诊断：场景111015的32/128 ms等效力漂移误差中位数2.465/9.595 N；114030为2.984/12.338 N。该指标采用冻结真实切线刚度，包含非线性/刚度变化，不能当成纯组织运动误差。',
            '估计不是唯一解释：114030的有效量测增益kg/ktrue中位1.031，越过历史筛选区间的比例仅2.2%，但128 ms预测误差仍大；S2两个审计场景存在9.2%/26.1%的增益越界，现有接受准则也不保证该闭环分布下始终有效。',
            '机制推断：FINAL原本只用于执行延迟附近的短前馈预测；将同一CA后验外推到完整128 ms、且把极端力漂移作为确定模型交给GN，会驱动不合适的参考整形，并与长恢复相互强化。当前证据未通过oracle预测或恢复替换独立量化每一环的因果份额。'],
        next_step=[
            '先改恢复与启动接口：使没有k_good的早期失联也能可靠重建接触，再在可信接触后重启量测；不能靠继续调运动权重解决从未启用的场景。',
            '随后改预测：为长时域建立可信范围/不确定性衰减，避免把短时CA外推原样视为128 ms确定漂移。用固定增益/已知运动的分环核验区分增益失配与时域失配。',
            '再检验优化如何使用不确定预测与恢复状态；本轮没有证据把增加GN迭代或换优化器列为第一优先项。保留估计增益偏差问题，但不将所有失败归给估计。'])
    save('v4_architecture.json', result)
    lines = ['# V4 架构验证','', '本结果来自指定诊断集与剂量集；未调用 final 确认集。调参只用 190000+。', '',
             '## 设计与实现','', *[f'- **{k}**: {v}' for k,v in result['implementation'].items()], '',
             '## 调参','',f"预算 {len(result['tuning']['trials'])}/61；最优参数 `{result['tuning']['best']['params']}`；目标 {result['tuning']['best']['objective']:.6f}。",
             '默认+40随机+20局部；每族10；五族等权 mean log(1+min(J,1000))。所有试验分数见 v4_tuning.json。', '',
             '原C1与Proposed参数按任务冻结，因此本轮是给定基线的架构试验，不能声称所有方法都在新调参集上做过同预算重新优化。', '',
             '## 指标（RMSE N，中位/均值；r=失控，L=长失联）','',
             '|条件|'+ '|'.join(METHODS)+'|', '|---|'+ '|'.join(['---']*len(METHODS))+'|']
    for group,t in tables.items():
        cells=[]
        for n in METHODS:
            m=t['metrics'][n]
            cells.append(f"{m['rmse_median']:.3f}/{m['rmse_mean']:.3f} r{m['runaway']} L{m['long_loss']} d{m['diverged']}")
        lines.append('|'+group+'|'+'|'.join(cells)+'|')
    lines += ['', '## 配对bootstrap：V4 − 对照，均值RMSE差（N）','',
              '|条件|Proposed：差 [95% CI]，有效配对数|C1：差 [95% CI]，有效配对数|','|---|---|---|']
    for group,t in tables.items():
        cells=[]
        for n in ('Proposed','C1'):
            b=t['paired'][n]
            cells.append(f"{b['mean_difference']:.3f} [{b['ci95'][0]:.3f}, {b['ci95'][1]:.3f}] n={b['n_pairs']}" if b['ci95'] else 'undefined')
        lines.append('|'+group+'|'+'|'.join(cells)+'|')
    if result['timing']:
        lines += ['', '## 单拍计算时间（单进程，ms）','', '|方法|P50|P99|','|---|---:|---:|']
        for n,m in result['timing']['methods'].items():
            lines.append(f"|{n}|{m['p50_ms']:.3f}|{m['p99_ms']:.3f}|")
    lines += ['', '## 结论与下一环','', result['conclusion']['verdict'], '',
              *['- '+s for s in result['conclusion']['supported_findings']], '',
              '### 失败证据与机理边界','', *['- '+s for s in result['conclusion']['failure_evidence']], '',
              '### 下一步（未在本轮实施）','', *['- '+s for s in result['conclusion']['next_step']], '',
              '### 执行与预算说明','', *['- '+s for s in result['execution_notes']]]
    lines += ['', '## 边界与疑点','', *['- '+s for s in result['limitations']], '',
              '发散时RMSE与最长失联无定义，不补零、不剔除失败后宣称稳健；均值仅在有限结果上，发散数另列。Bootstrap每条件4000次，无多重比较校正。', '',
              '复现：`python -X utf8 experiments/method_redesign/v4_evaluate.py --stage all --workers 12`。删除本任务新增的结果JSON后可从头运行；已有检查点默认续跑。']
    (OUT/'v4_architecture.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('REPORT WRITTEN',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--stage', choices=('check','tune','evaluate','timing','audit','report','all'),default='all')
    ap.add_argument('--workers',type=int,default=12); ap.add_argument('--worker-dir'); args=ap.parse_args()
    if args.worker_dir:
        worker_loop(args.worker_dir); sys.exit(0)
    t0=time.perf_counter()
    def audit():
        subprocess.run([sys.executable,'-X','utf8',str(HERE/'v4_failure_analysis.py')], check=True)
    for name,fn in [('check',check),('tune',lambda:tune(args.workers)),('evaluate',lambda:evaluate(args.workers)),('timing',timing),('audit',audit),('report',report)]:
        if args.stage in (name,'all'):
            fn()
    print(f'FINISHED {args.stage} {time.perf_counter()-t0:.1f}s',flush=True)
