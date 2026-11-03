"""条件线性 V0 的 1+40+20 等预算调参。

调参场景只用五族各 6 个 190000+ 场景；不读取也不生成 130000+ 确认集。
"""
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
from cohort_eval import FINAL  # noqa: E402
import common  # noqa: E402
import scenes_v2 as SV  # noqa: E402
from conditional_filter import ConditionalLinearContact, MeasuredOnly  # noqa: E402
from v1_servo_id import ServoID  # noqa: E402

OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)

DEFAULT = dict(q_lnk0=1e-4, q_lnk1=1e-2, p_stay_stiff=.98,
               gain_rel_max=math.log(1.2), consistency_max=float(FINAL['chi2_max']),
               r_scale=float(FINAL['r_scale']), p_stay=float(FINAL['p_stay']),
               noise_scale=1., c=float(FINAL['c']), rho_d=float(FINAL['rho_d']))

BOUNDS = {
    'q_lnk0': (1e-7, 3e-3), 'q_lnk1': (3e-5, .3),
    'p_stay_stiff': (.90, .9999), 'p_stay': (.95, .9999),
    'r_scale': (.03, 3.), 'noise_scale': (.5, 2.),
    'gain_rel_max': (math.log(1.10), math.log(1.40)),
    'consistency_max': (1.2, 5.), 'c': (.2, 2.), 'rho_d': (.005, .2),
}


def scenes():
    out = []
    for fi, fam in enumerate(SV.FAMILIES):
        rng = np.random.default_rng(190000 + 1000 * fi)
        for i in range(6):
            p, T = SV.make(fam, rng)
            out.append((fam, 190000 + 1000 * fi + i, p, T))
    return out


def _build_params(overrides):
    p = dict(FINAL)
    q = dict(overrides)
    if 'noise_scale' in q:
        q['noise'] = .04 * q.pop('noise_scale')
    p.update(q)
    return p


def one(job):
    idx, params, scene, method = job
    fam, sid, plant_p, T = scene
    try:
        if method == 'conditional_v0':
            pol = ConditionalLinearContact(single_model=False, **_build_params(params))
        elif method == 'v1_proposed':
            pol = ServoID(single_model=False, **dict(FINAL))
        else:
            raise ValueError(method)
        log, _ = common.run(MeasuredOnly(pol), plant_p, T, sid)
        m = common.task_metrics(log, 20.8)
        J = float(m['J']) if np.isfinite(m['J']) else 1e3
        return dict(idx=idx, family=fam, sid=sid, L=float(math.log1p(min(J, 1000.))),
                    J=min(J, 1000.), rmse=float(m['rmse']) if np.isfinite(m['rmse']) else None,
                    peak=float(m['peak']) if np.isfinite(m['peak']) else None,
                    diverged=bool(m['diverged']))
    except Exception as ex:
        return dict(idx=idx, family=fam, sid=sid, L=math.log1p(1000.), J=1000.,
                    rmse=None, peak=None, diverged=True, error=repr(ex))


def evaluate(params, scenes_, workers, method='conditional_v0'):
    jobs = [(i, params, s, method) for i, s in enumerate(scenes_)]
    rows = []
    with ProcessPoolExecutor(max_workers=min(int(workers), 14)) as pool:
        fs = [pool.submit(one, j) for j in jobs]
        for f in as_completed(fs):
            rows.append(f.result())
    rows.sort(key=lambda x: x['idx'])
    return dict(L=float(np.mean([r['L'] for r in rows])),
                rmse=float(np.nanmean([r['rmse'] if r['rmse'] is not None else np.nan for r in rows])),
                rows=rows)


def random_candidate(rng):
    q = {}
    for k, (lo, hi) in BOUNDS.items():
        if k.startswith('p_stay') or k in ('consistency_max', 'c', 'rho_d', 'noise_scale'):
            q[k] = float(rng.uniform(lo, hi))
        else:
            q[k] = float(math.exp(rng.uniform(math.log(lo), math.log(hi))))
    return q


def local_candidate(base, rng):
    q = dict(base)
    for k, (lo, hi) in BOUNDS.items():
        v = float(q[k])
        if k.startswith('p_stay') or k in ('consistency_max', 'c', 'rho_d', 'noise_scale'):
            v += rng.normal(0., .12 * (hi - lo))
        else:
            v *= math.exp(rng.normal(0., .35))
        q[k] = float(np.clip(v, lo, hi))
    return q


def score_record(i, kind, params, result):
    return dict(index=i, kind=kind, params=params, L=result['L'], rmse=result['rmse'],
                rows=result['rows'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=14)
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    sc = scenes()
    protocol = dict(method='conditional_v0', random_evaluations=40, local_evaluations=20,
                    total_evaluations=61, random_seed=190000, local_seed=190001,
                    scenes_per_family=6, families=list(SV.FAMILIES), confirmation_accessed=False,
                    objective='mean over scenes of log(1+min(J,1000))', bounds=BOUNDS)
    (OUT / 'tuning_protocol.json').write_text(json.dumps(protocol, indent=2, default=float), encoding='utf8')
    path = OUT / 'tuning_candidates.json'
    records = json.loads(path.read_text(encoding='utf8')) if args.resume and path.exists() else []
    done = {int(r['index']) for r in records}
    candidates = [(0, 'default', DEFAULT)]
    rng = np.random.default_rng(190000)
    candidates.extend((i, 'random', random_candidate(rng)) for i in range(1, 41))

    def run_candidates(batch):
        nonlocal records, done
        for idx, kind, params in batch:
            if idx in done:
                continue
            tic = time.perf_counter()
            result = evaluate(params, sc, args.workers)
            records.append(score_record(idx, kind, params, result))
            records.sort(key=lambda x: x['index'])
            done.add(idx)
            path.write_text(json.dumps(records, indent=2, default=float), encoding='utf8')
            print(idx, kind, result['L'], result['rmse'], 'elapsed', round(time.perf_counter() - tic, 1), flush=True)

    run_candidates(candidates)
    base = min(records, key=lambda r: r['L'])['params']
    rng2 = np.random.default_rng(190001)
    local = [(41 + i, 'local', local_candidate(base, rng2)) for i in range(20)]
    run_candidates(local)

    best = min(records, key=lambda r: r['L'])
    params = dict(best['params'])
    (OUT / 'tuned_params.json').write_text(json.dumps(params, indent=2, default=float), encoding='utf8')
    (OUT / 'tuning_summary.json').write_text(json.dumps(dict(protocol=protocol, best=best,
                                                              evaluated=len(records)), indent=2, default=float), encoding='utf8')
    print('best', best['index'], best['L'], params)


if __name__ == '__main__':
    main()
