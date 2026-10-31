"""YuenType smoke/time-alignment audit and isolated confirm evaluation.

Run: python -X utf8 yuen_smoke.py [--confirm] [--workers 16]
Default smoke uses two development scenes; --confirm additionally reads the
frozen confirm/main.json Proposed rows, writes ONLY confirm_yuen.json, and
prints per-family paired bootstrap intervals (YuenType minus Proposed).
No Proposed controller is rerun: its current source may differ from that file.
"""
import argparse
import hashlib
import json
import os
from multiprocessing import Pool
from pathlib import Path
from unittest.mock import patch

import numpy as np
import common
import comparators as C
import scenes_v2 as SV


HERE = Path(common.HERE)
# An output override supports running a reviewable overlay when the repo is read-only.
OUTPUT = Path(os.environ.get('YUEN_OUTPUT', str(HERE)))
FAIRNESS = [
    '特权传感器直接量测未变形局部表面；接触中持续更新，比 Yuen 原文更强。',
    '空间配准和延迟已知，无视觉遮挡、接触形变或重建误差；噪声仅为 0.1 mm 高斯噪声。',
    '位置接口 PI 和位置前馈是任务适配，不是原文力矩控制器的完整复现。',
    '相邻周期预测位置之差包含 EKF/重拟合校正跳变；相较平滑速度/加速度前馈可能不利。',
    '无预热数据；2 s 后仅用已到达量测初始化，原文在接触前训练约 5 s。',
    '周期模型面对非周期推拉、静态扫描凸起存在失配；固定刚度不能适应空间刚度跨度。',
    '传感器 seed=0 固定，使所有参数共享噪声；确认配对来自既有 Proposed 文件。',
]


def effective_params(params):
    return dict(C.YuenType.params, **params)


def smoke(params=None):
    """Two full development episodes plus a guarded prefix replay per scene."""
    params = params or {}
    scenes = SV.scenes('dev')
    out = []
    for family in ('T3m', 'T2m-stiff'):
        fam, sid, p, T = next(s for s in scenes if s[0] == family)
        pol = C.YuenType(p, **params)
        log, _ = common.run(pol, p, T, sid)
        assert log is not None and len(log) > 100
        m = common.task_metrics(log, 20.8)
        assert all(np.isfinite(m[k]) for k in ('rmse', 'peak', 'J'))
        errors = []
        for r in pol.trace:
            if r['prediction'] is not None and 1. <= r['t'] and r['target_t'] <= log[-1, 0]:
                xx = float(np.interp(r['target_t'], log[:, 0], log[:, 1]))
                truth = common.P.surface(xx, r['target_t'], p)[0]  # OFFLINE evaluation only
                errors.append(r['prediction']-truth)
        st = pol.sensor_trace
        assert np.allclose([r['sample_t'] for r in st], np.arange(len(st))/pol.f_s, atol=1e-12)
        for r in st:
            assert r['used_t']+1e-12 >= r['available_t']
            assert r['used_t']-r['available_t'] < pol.Tc+1e-10
            assert r['used_t']-r['sample_t'] >= pol.tau_m-1e-12
        assert pol.trace[0]['prediction'] is None  # no sensor before 30 ms
        for r in pol.trace:
            assert abs(r['target_t']-r['t']-((pol.ndel+.5)*pol.Tc+pol.tau)) < 1e-10
            if r['sample_t'] is not None:
                assert r['sample_t']+pol.tau_m <= r['t']+1e-12

        # The guard sees every surface() call made by the replayed policy.
        # All future rows and unobservable truth columns are poisoned with NaNs.
        # Also change a FUTURE event in p: commands before it must be identical.
        cutoff = min(1.5, float(p[16])-.1) if p[15] else 1.5
        assert cutoff > .1
        changed = np.array(p, copy=True)
        changed[15] = 1; changed[16:20] = [cutoff+.5, .05, .2, .5]
        a, b = C.YuenType(p, **params), C.YuenType(changed, **params)
        platform = common.P.PLATFORMS['PC125']; a.reset(platform); b.reset(platform)
        real_surface = common.P.surface
        now = [0.]; queries = []

        def guarded_surface(x, ts, pp):
            assert ts+pol.tau_m <= now[0]+1e-12, 'future/unarrived surface queried'
            queries.append(float(ts))
            return real_surface(x, ts, pp)

        n_replay = int(np.searchsorted(log[:, 0], cutoff))
        poisoned = np.full_like(log, np.nan)
        with patch.object(common.P, 'surface', guarded_surface):
            for n in range(1, n_replay+1):
                now[0] = float(log[n-1, 0])
                poisoned[n-1, [0, 1, 2, 5, 11]] = log[n-1, [0, 1, 2, 5, 11]]
                clean = a(common.P.Obs(log[:n].copy(), n))
                dirty = b(common.P.Obs(poisoned, n))
                assert np.isfinite(dirty) and clean == dirty, 'future/truth dependence'
        assert len(queries) == len(a.sensor_trace)+len(b.sensor_trace) > 0
        row = dict(family=fam, sid=sid, prediction_rmse_mm=float(np.sqrt(np.mean(np.square(errors)))*1000),
                   prediction_peak_mm=float(np.max(np.abs(errors))*1000), force_rmse=m['rmse'],
                   force_peak=m['peak'], loss_max=m['loss_max'], refits=pol.refits,
                   sensor_samples=len(st), min_latency_ms=min(r['used_t']-r['sample_t'] for r in st)*1000,
                   max_latency_ms=max(r['used_t']-r['sample_t'] for r in st)*1000,
                   prefix_replay_ticks=n_replay, alignment_and_no_future_leak=True)
        out.append(row)
        print('smoke:', json.dumps(row, ensure_ascii=False), flush=True)
    return out


def _episode(args):
    params, (fam, sid, p, T) = args
    log, _ = common.run(C.YuenType(p, **params), p, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(method='YuenType', family=fam, sid=sid, ok=bool(not m['diverged']),
                rmse=m['rmse'], peak=m['peak'], loss=m['loss'], loss_frac=m['loss'],
                loss_max=m['loss_max'], J=m['J'])


def paired(yuen, proposed, key, n_boot=4000):
    """Paired mean difference; percentile CI, scene resampling, fixed seed 0."""
    a = {r['sid']: r for r in yuen}; b = {r['sid']: r for r in proposed}
    assert set(a) == set(b), 'scene IDs differ from frozen Proposed confirmation'
    ids = sorted(k for k in a if a[k]['ok'] and b[k]['ok'])
    d = np.array([a[k][key]-b[k]['loss_frac' if key == 'loss' else key] for k in ids])
    if not len(d):
        return dict(n=0, mean=None, ci=[None, None], excluded=len(a))
    assert np.all(np.isfinite(d))
    rng = np.random.default_rng(0)
    means = np.mean(rng.choice(d, (n_boot, len(d))), axis=1)
    return dict(n=len(d), mean=float(np.mean(d)), ci=np.percentile(means, [2.5, 97.5]).tolist(),
                excluded=len(a)-len(d), per_scene=d.tolist(), sid=ids)


def confirm(workers=16):
    source = HERE/'results_v2/confirm/main.json'
    tuned = OUTPUT/'tuned_v2/YuenType.json'
    dest = OUTPUT/'results_v2/confirm_yuen.json'
    if dest.exists():
        raise FileExistsError(f'Refusing to overwrite completed confirmation: {dest}')
    source_bytes = source.read_bytes()
    proposed = [r for r in json.loads(source_bytes) if r['method'] == 'Proposed']
    sc = SV.scenes('confirm')
    assert len(proposed) == len(sc) == len({(r['family'], r['sid']) for r in proposed})
    assert {(r['family'], r['sid']) for r in proposed} == {(f, sid) for f, sid, _, _ in sc}
    record = json.loads(tuned.read_text(encoding='utf8')); params = record['best']['params']
    assert record['budget'] == len(record['all']) == 61
    jobs = [(params, s) for s in sc]
    if workers == 1:
        rows = list(map(_episode, jobs))
    else:
        with Pool(workers) as pool:
            rows = pool.map(_episode, jobs, chunksize=1)
    summary = {}
    for fam in SV.FAMILIES:
        yy = [r for r in rows if r['family'] == fam]; good = [r for r in yy if r['ok']]
        pp = [r for r in proposed if r['family'] == fam]
        stats = dict(n=len(yy), n_fail=len(yy)-len(good),
                     rmse_mean=float(np.mean([r['rmse'] for r in good])) if good else None,
                     peak_max=max((r['peak'] for r in good), default=None),
                     loss_max=max((r['loss_max'] for r in good), default=None),
                     paired={k: paired(yy, pp, k) for k in ('rmse', 'peak', 'loss', 'loss_max', 'J')})
        summary[fam] = stats
    payload = dict(method='YuenType', split='confirm', params=effective_params(params),
                   sensor=dict(f_s=C.YuenType.f_s, tau_m=C.YuenType.tau_m, sigma_m=C.YuenType.sigma_m, seed=0),
                   tuning=dict(budget=61, split='dev', n_scenes=30, objective='mean log(1+min(J,1e3))',
                               random_seed=777, local_seed=888, best=record['best']),
                   reference=str(source.relative_to(HERE)), reference_sha256=hashlib.sha256(source_bytes).hexdigest(),
                   tuned_sha256=hashlib.sha256(tuned.read_bytes()).hexdigest(),
                   scene_definition_sha256=hashlib.sha256(Path(SV.__file__).read_bytes()).hexdigest(),
                   bootstrap=dict(difference='YuenType-Proposed', n=4000, seed=0, confidence=.95,
                                  statistic='paired mean difference', failures='reported; successful pairs only'),
                   rows=rows, summary=summary, fairness=FAIRNESS)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(payload, ensure_ascii=False, indent=1, allow_nan=False), encoding='utf8')
    assert source.read_bytes() == source_bytes
    print_report(payload)
    return payload


def print_report(result):
    print('\n实现：延迟特权表面传感器 + 在线傅里叶 EKF/残差重拟合 + 固定刚度增量 PI/预测位置前馈。')
    print('最优参数：', json.dumps(result['params'], ensure_ascii=False))
    print('开发集：61 次，30 场景/次；L=', result['tuning']['best']['L'])
    print('确认集：RMSE 为逐场景均值；peak/loss_max 为该族最坏值；差值/区间是逐场景差的均值及95%配对 bootstrap CI。')
    def fmt(value): return 'NA' if value is None else f'{value:.6f}'
    for fam, s in result['summary'].items():
        print(f"{fam}: n={s['n']} fail={s['n_fail']} RMSE={fmt(s['rmse_mean'])} N peak={fmt(s['peak_max'])} N loss_max={fmt(s['loss_max'])} s")
        for key, p in s['paired'].items():
            print(f"  {key}: YuenType-Proposed={fmt(p['mean'])}, 95% CI={p['ci']}, n={p['n']}")
    print('公平性/适用边界：')
    for item in FAIRNESS: print(' -', item)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--confirm', action='store_true')
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()
    tp = OUTPUT/'tuned_v2/YuenType.json'
    params = json.loads(tp.read_text(encoding='utf8'))['best']['params'] if tp.exists() else {}
    smoke(params)
    if args.confirm: confirm(args.workers)
