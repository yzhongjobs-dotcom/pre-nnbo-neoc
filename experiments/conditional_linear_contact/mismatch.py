"""伺服滞后失配：tau_a ×0.6 与 ×1.5，各族 20 场景、三方法配对。"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
from closed_loop import build, load_params, one, summarize  # noqa: E402
import cohort_eval as CE  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--params', default=None); ap.add_argument('--workers', type=int, default=14)
    args = ap.parse_args()
    params = load_params(args.params)
    scenes = CE.cohort('diag', 20)
    rows = []
    jobs = [(m, s, params, lag) for lag in (.6, 1.5)
            for m in ('v1_proposed', 'joint_lean', 'conditional_v0') for s in scenes]
    with ProcessPoolExecutor(max_workers=min(args.workers, 14)) as pool:
        fs = [pool.submit(one, j) for j in jobs]
        for i, f in enumerate(as_completed(fs), 1):
            rows.append(f.result())
            if i % 10 == 0: print(i, '/', len(jobs), flush=True)
    payload = dict(protocol=dict(seed_split='diagnostic 110000+', scenarios_per_family=20,
                                 lag_scales=[.6, 1.5], confirmation_accessed=False, params=params),
                   rows=rows,
                   summary=[x for lag in (.6, 1.5) for x in summarize(rows, lag)])
    path = HERE / 'results' / 'mismatch.json'
    path.write_text(json.dumps(payload, indent=2, default=float), encoding='utf8')
    print(path)


if __name__ == '__main__':
    main()
