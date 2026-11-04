"""Select offline estimator representatives without using closed-loop results."""
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

from failure_replay import HERE, ES, cell_scene, Measured, Recording, FINAL, PP

BASE = 310000
CELLS = [(c, t) for c in ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor')
         for t in ('soft', 'mid', 'hard')] + [('rib', 'mixed'), ('combo', 'mixed')]


def main():
    jobs = [(c, t, i, BASE) for c, t in CELLS for i in range(4)]
    from failure_replay import job
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, jobs, chunksize=1))
    score = {}
    for name in ES.OFFLINE_GRID:
        per = []
        for c, t in CELLS:
            x = [r.get('methods', {}).get(name, {}).get('p90', np.nan)
                 for r in rows if r['cond'] == c and r['tissue'] == t]
            x = [z for z in x if np.isfinite(z)]
            if x: per.append(float(np.median(x)))
        score[name] = dict(mean_cell_p90=float(np.mean(per)), worst_cell_p90=float(np.max(per)))
    def best(pred): return min((k for k in score if pred(k)), key=lambda k: score[k]['mean_cell_p90'])
    picks = {
        'R2': best(lambda k: k.startswith('R2|')),
        'Butterworth': best(lambda k: k.startswith('B')),
        'Kaiser': best(lambda k: k.startswith('K')),
        'Projection': best(lambda k: k.startswith('P')),
    }
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/'filter_selection_rows.json').write_text(json.dumps(rows, indent=1, allow_nan=True), encoding='utf-8')
    (HERE/'results'/'filter_selection.json').write_text(json.dumps(dict(base=BASE, score=score, picks=picks), indent=1, allow_nan=True), encoding='utf-8')
    for k, v in picks.items(): print(k, v, score[v])


if __name__ == '__main__': main()
