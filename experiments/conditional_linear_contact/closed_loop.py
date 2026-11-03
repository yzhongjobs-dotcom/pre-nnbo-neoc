"""诊断集闭环：五族各 50 场景，三种方法配对。"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent / 'contact_information'))
from cohort_eval import FINAL, cohort  # noqa: E402
import common  # noqa: E402
from joint_policy import JointPolicy  # noqa: E402
from conditional_filter import ConditionalLinearContact, MeasuredOnly  # noqa: E402
from v1_servo_id import ServoID  # noqa: E402


def load_params(path):
    if not path:
        return {}
    return json.loads(Path(path).read_text(encoding='utf8'))


def build(name, params):
    p = dict(FINAL); q = dict(params)
    if 'noise_scale' in q:
        q['noise'] = .04 * q.pop('noise_scale')
    p.update(q)
    if name == 'v1_proposed':
        return ServoID(single_model=False, **dict(FINAL))
    if name == 'joint_lean':
        return JointPolicy(wave=None, new_fit=True, rate=False, amplitude_scale=1., **dict(FINAL))
    if name == 'conditional_v0':
        return ConditionalLinearContact(single_model=False, **p)
    raise ValueError(name)


def one(job):
    name, scene, params, lag_scale = job
    fam, sid, plant_p, T = scene
    platform = dict(common.P.PLATFORMS['PC125'])
    platform['tau_a'] *= lag_scale
    pol = build(name, params)
    inner = pol
    tic = time.perf_counter()
    log, _ = common.run(MeasuredOnly(pol), plant_p, T, sid, platform=platform)
    elapsed = time.perf_counter() - tic
    m = common.task_metrics(log, 20.8)
    calls = np.asarray(inner.t_call[126:], dtype=float) * 1000.
    accepted = getattr(inner, 'accepted_count', 0)
    total = len(getattr(inner, 'accept_history', []))
    return dict(method=name, family=fam, sid=sid, lag_scale=lag_scale,
                rmse=None if not np.isfinite(m['rmse']) else float(m['rmse']),
                J=float(m['J']) if np.isfinite(m['J']) else 1000.,
                peak=None if not np.isfinite(m['peak']) else float(m['peak']),
                loss=float(m['loss']) if np.isfinite(m['loss']) else 1.,
                loss_max=float(m['loss_max']) if np.isfinite(m['loss_max']) else 1e9,
                diverged=bool(m['diverged']),
                elapsed_s=elapsed,
                call_ms_P50=None if len(calls) == 0 else float(np.percentile(calls, 50)),
                call_ms_P99=None if len(calls) == 0 else float(np.percentile(calls, 99)),
                acceptance_rate=None if total == 0 else float(accepted / total),
                accepted=int(accepted), accept_trials=int(total),
                mode_probability_mean=(np.mean([x['mu'] for x in getattr(inner, 'mode_history', [])], axis=0).tolist()
                                       if getattr(inner, 'mode_history', []) else None))


def bootstrap(values, B=2000, seed=7):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return [None, None]
    rng = np.random.default_rng(seed)
    means = np.empty(B)
    for i in range(B):
        means[i] = np.mean(values[rng.integers(0, len(values), len(values))])
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def summarize(rows, lag_scale=1.):
    out = []
    methods = ('v1_proposed', 'joint_lean', 'conditional_v0')
    for fam in sorted({r['family'] for r in rows if r['lag_scale'] == lag_scale}):
        for method in methods:
            R = [r for r in rows if r['family'] == fam and r['method'] == method and r['lag_scale'] == lag_scale]
            if not R:
                continue
            finite = [r for r in R if r['rmse'] is not None]
            out.append(dict(family=fam, method=method, n=len(R),
                            rmse_median=float(np.median([r['rmse'] for r in finite])) if finite else None,
                            rmse_mean=float(np.mean([r['rmse'] for r in finite])) if finite else None,
                            peak_median=float(np.median([r['peak'] for r in finite])) if finite else None,
                            peak_mean=float(np.mean([r['peak'] for r in finite])) if finite else None,
                            peak_gt_50=sum(r['peak'] is None or r['peak'] > 50 for r in R),
                            diverged=sum(r['diverged'] for r in R),
                            long_loss_gt_05=sum(r['loss_max'] > .5 for r in R),
                            acceptance_rate=float(np.mean([r['acceptance_rate'] for r in R if r['acceptance_rate'] is not None]))
                            if any(r['acceptance_rate'] is not None for r in R) else None,
                            call_ms_P50=float(np.median([r['call_ms_P50'] for r in R if r['call_ms_P50'] is not None]))
                            if any(r['call_ms_P50'] is not None for r in R) else None,
                            call_ms_P99=float(np.percentile([r['call_ms_P99'] for r in R if r['call_ms_P99'] is not None], 99))
                            if any(r['call_ms_P99'] is not None for r in R) else None))
        b = {(r['family'], r['sid']): r for r in rows
             if r['family'] == fam and r['method'] == 'v1_proposed' and r['lag_scale'] == lag_scale}
        for method in ('joint_lean', 'conditional_v0'):
            c = {(r['family'], r['sid']): r for r in rows
                 if r['family'] == fam and r['method'] == method and r['lag_scale'] == lag_scale}
            common_keys = sorted(set(b) & set(c))
            d_rmse = [c[k]['rmse'] - b[k]['rmse'] for k in common_keys
                      if c[k]['rmse'] is not None and b[k]['rmse'] is not None]
            d_J = [c[k]['J'] - b[k]['J'] for k in common_keys]
            out.append(dict(family=fam, comparison=f'{method}-v1_proposed', n=len(d_rmse),
                            delta_rmse_mean=float(np.mean(d_rmse)) if d_rmse else None,
                            delta_rmse_bootstrap95=bootstrap(d_rmse),
                            delta_J_mean=float(np.mean(d_J)) if d_J else None,
                            delta_J_bootstrap95=bootstrap(d_J, seed=11)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--params', default=None)
    ap.add_argument('--workers', type=int, default=14)
    ap.add_argument('--n', type=int, default=50)
    ap.add_argument('--families', default=','.join(('S2', 'T3m', 'T2m-soft', 'T2m-stiff', 'T2f-stiff')))
    ap.add_argument('--lag-scale', type=float, default=1.)
    ap.add_argument('--label', default='tuned')
    args = ap.parse_args()
    params = load_params(args.params)
    fams = args.families.split(',')
    scenes = [s for s in cohort('diag', args.n) if s[0] in fams]
    jobs = [(m, s, params, args.lag_scale) for m in ('v1_proposed', 'joint_lean', 'conditional_v0') for s in scenes]
    path = HERE / 'results' / f'closed_loop_{args.label}_lag{args.lag_scale:g}.json'
    rows = []
    with ProcessPoolExecutor(max_workers=min(args.workers, 14)) as pool:
        fs = [pool.submit(one, j) for j in jobs]
        for i, f in enumerate(as_completed(fs), 1):
            rows.append(f.result())
            if i % 10 == 0:
                path.write_text(json.dumps(dict(protocol=dict(seed_split='diagnostic 110000+',
                                                              confirmation_accessed=False, n=args.n,
                                                              lag_scale=args.lag_scale, params=params), rows=rows),
                                                     indent=2, default=float), encoding='utf8')
                print(i, '/', len(jobs), flush=True)
    payload = dict(protocol=dict(seed_split='diagnostic 110000+', confirmation_accessed=False,
                                 n=args.n, lag_scale=args.lag_scale, params=params), rows=rows,
                   summary=summarize(rows, args.lag_scale))
    path.write_text(json.dumps(payload, indent=2, default=float), encoding='utf8')
    print(path)


if __name__ == '__main__':
    main()
