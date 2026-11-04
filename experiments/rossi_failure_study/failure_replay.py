"""Offline attribution replay for the four retained estimator families."""
import argparse
import json
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
import estimators as ES
from cohort_eval import FINAL
from v1_servo_id import ServoID


class Measured:
    """Expose only causal measurement columns and reset on the nominal platform."""
    direct = True
    def __init__(self, inner): self.inner = inner
    def reset(self, p): self.inner.reset(dict(PP.PLATFORM))
    def __call__(self, o):
        a = np.full((o.n, 16), np.nan)
        a[:, [0, 1, 2, 5, 11, 12, 13]] = o.log[:o.n, [0, 1, 2, 5, 11, 12, 13]]
        return self.inner(PP.OLD.Obs(a, o.n))


class OracleMeasured(Measured):
    """Privileged reference wrapper used only by C_neutral/oracle rows."""
    def __call__(self, o):
        a = np.full((o.n, 16), np.nan)
        a[:, [0, 1, 2, 5, 8, 11, 12, 13]] = o.log[:o.n, [0, 1, 2, 5, 8, 11, 12, 13]]
        return self.inner(PP.OLD.Obs(a, o.n))


class Recording(ServoID):
    """V1 trajectory generator for offline attribution; estimator is not swapped."""
    def reset(self, p):
        super().reset(p); self.probe = []
    def step(self, o):
        prev = self.x_prev; cmd = super().step(o)
        k = self.k_good if np.isfinite(self.k_good) else self.k_seek
        self.probe.append((self.x_prev-prev)/k)
        return cmd


def cell_scene(cond, tissue, i, base, overrides=None, v=None):
    code = sum(ord(c)*(j+1) for j, c in enumerate(cond+tissue)) % 9973
    return SC.make(cond, tissue, np.random.default_rng([base, code, i]), T=10., v=v,
                   overrides=overrides)


def _summary(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x): return dict(n=0, p50=np.nan, mean=np.nan, p90=np.nan, rmse=np.nan)
    return dict(n=int(len(x)), p50=float(np.percentile(abs(x), 50)),
                mean=float(np.mean(x)), p90=float(np.percentile(abs(x), 90)),
                rmse=float(np.sqrt(np.mean(x*x))))


def job(a):
    cond, tissue, i, base = a
    sc = cell_scene(cond, tissue, i, base)
    pol = Recording(single_model=False, **FINAL)
    log = PP.episode(sc, Measured(pol), base+i)
    if log is None:
        return dict(cond=cond, tissue=tissue, i=i, failed=True)
    t = log[:, 0]; kt = log[:, 8]
    mask = (t >= 1.5) & (kt > 0) & (log[:, 5] > .5)
    out = {}
    for name in ES.OFFLINE_GRID:
        try:
            v = ES.offline_named(log, name)
            ok = mask & np.isfinite(v)
            out[name] = _summary(np.log(v[ok]/kt[ok]))
        except (ValueError, FloatingPointError) as e:
            out[name] = dict(n=0, error=str(e))
    return dict(cond=cond, tissue=tissue, i=i, failed=False, methods=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=12); ap.add_argument('--base', type=int, default=320000)
    ap.add_argument('--label', default='failure_replay')
    a = ap.parse_args()
    cells = [(c, t) for c in ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor')
             for t in ('soft', 'mid', 'hard')] + [('rib', 'mixed'), ('combo', 'mixed')]
    jobs = [(c, t, i, a.base) for c, t in cells for i in range(a.n)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, jobs, chunksize=1))
    (HERE/'results').mkdir(exist_ok=True)
    out = HERE/'results'/f'{a.label}_rows.json'
    out.write_text(json.dumps(rows, indent=1, allow_nan=True), encoding='utf-8')
    # A compact family summary is derived from rows, never from a separately
    # tuned table, so the report can reproduce it exactly.
    score = {}
    for name in ES.OFFLINE_GRID:
        cellscores = []
        for c, t in cells:
            vals = [r.get('methods', {}).get(name, {}).get('p90', np.nan)
                    for r in rows if r['cond'] == c and r['tissue'] == t]
            vals = [v for v in vals if np.isfinite(v)]
            if vals: cellscores.append(float(np.median(vals)))
        score[name] = dict(mean_cell_p90=float(np.mean(cellscores)) if cellscores else np.nan,
                           worst_cell_p90=float(np.max(cellscores)) if cellscores else np.nan)
    (HERE/'results'/f'{a.label}_summary.json').write_text(json.dumps(score, indent=1, allow_nan=True), encoding='utf-8')
    print('offline rows', len(rows), 'failed', sum(r.get('failed', False) for r in rows))
    for k in sorted(score, key=lambda z: score[z]['mean_cell_p90'])[:12]: print(k, score[k])


if __name__ == '__main__':
    main()
