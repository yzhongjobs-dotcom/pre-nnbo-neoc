"""S2 supplement: does "update = scaling" lose much against re-optimising per stiffness?

On homogeneous breathing tissue at four stiffness levels, the controller uses the true
stiffness (to isolate the scaling claim from estimation) and the drift feed-forward off.
For each level we compare the composite cost of the scaled unit design A1 with the best
of a local grid of (D1, KPF1) re-optimised for that level on the actual plant.
Run: python -X utf8 scaling_check.py -> results/S2_scaling.json
"""
import json, os
import numpy as np
from multiprocessing import Pool
import common
from common import S
import proposed as Pm
import run_experiments as RX

LEVELS = (1.5e4, 6e4, 2.5e5, 1.0e6)       # tissue modulus [Pa]
GRID_D = (0.5, 0.7, 1.0, 1.4, 2.0)          # multipliers on D1
GRID_K = (0.0, 0.5, 1.0, 1.5)


def scenes(E):
    r = np.random.default_rng(int(E) % 9973)
    out = []
    for i in range(4):
        p, T = S.make('N1', r, dict(E=E*float(r.uniform(.9, 1.1))))
        out.append((p, T, 9900+i))
    return out


def _cost(args):
    E, dmul, kpf = args
    base = RX.tuned('Proposed')
    kw = dict(M1=base['M1'], D1=base['D1']*dmul, KPF1=kpf if kpf is not None else base['KPF1'], ki=0.0,
              feedforward=False)
    Js = []
    for p, T, sid in scenes(E):
        log, _ = common.run(Pm.Oracle(**kw), p, T, sid)
        Js.append(min(common.task_metrics(log, None)['J'], 1e3))
    return dict(E=E, dmul=dmul, kpf=kpf, J=float(np.mean(Js)))


if __name__ == '__main__':
    jobs = [(E, 1.0, None) for E in LEVELS]+[(E, d, k) for E in LEVELS for d in GRID_D for k in GRID_K]
    with Pool(16) as pool:
        res = pool.map(_cost, jobs, chunksize=1)
    out = []
    for E in LEVELS:
        a1 = [r for r in res if r['E'] == E and r['kpf'] is None][0]['J']
        grid = [r for r in res if r['E'] == E and r['kpf'] is not None]
        best = min(grid, key=lambda r: r['J'])
        out.append(dict(E=E, J_scaled_A1=a1, J_best_reoptimised=best['J'], best=best,
                        gap=(a1-best['J'])/max(best['J'], 1e-9)))
    json.dump(out, open(os.path.join(RX.RES, 'S2_scaling.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))
