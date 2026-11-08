"""Which gain does the closed loop need? Feed window estimate = factor * logged |K(20.8 Hz)|, gates unchanged, and sweep factor.

Run on scenes A and B of the overlay protocol, with and without the probe mass beyond the sensor. Also records how far the estimator's own output sits
from the logged truth in each scene (median of k_hat/|K| over the episode). Diagnostic only; truth is only used through the factor.
Usage: python -X utf8 gain_need_sweep.py
"""
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import overlay_compare as O
import method3 as M
import physics_plant3 as P3

B = O.B
HERE = Path(__file__).resolve().parent
FACTORS = (.7, .85, 1., 1.2, 1.5, 2., 3.)


def job(a):
    cell, mass, factor = a
    O.setup()
    cfg = O.chosen()['Huber-D4']
    sc = B.scene(cell, O.BASE, 0.)
    if not mass:
        sc['tool_mass'] = 0.; sc['tool_mass_hat'] = 0.
    pol = B.build('Huber-D4', cfg, sc)
    ratio = []
    est = pol.estimator; original = est.estimate
    def estimate(o, n, policy, original=original):
        k, sig, chi = original(o, n, policy)
        truth = getattr(o, 'truth_k', 0.)
        if np.isfinite(k) and truth > 200.:
            ratio.append(k/truth)
            if factor is not None:
                sig = sig*factor*truth/k if k > 0 else sig
                k = factor*truth
        return k, sig, chi
    est.estimate = estimate
    log = P3.episode(sc, pol, O.sensor_seed(cell), diagnostic=True)
    m = M.metric(log, pol, sc)
    c, t, dy, ph, i = cell
    w = log[:, 0] >= 1.
    return dict(scene=c, dynamics=dy, phase=ph, i=i, mass=mass, factor=factor, L=float(np.log1p(m['J'])), rmse=float(m['rmse']), peak=float(m['peak']),
                k_med=float(np.median(log[w, 21])), own_ratio=float(np.median(ratio)) if ratio else None, tool=float(sc['tool_mass']), kr=float(sc['kr']))


def main():
    start = time.perf_counter()
    cells = [c for c in O.cells() if c[0] in ('breath_rib', 'breath_fast')]
    jobs = [(c, mass, f) for mass in (True, False) for f in (None,)+FACTORS for c in cells]
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(job, jobs, chunksize=2))
    (HERE/'results'/'gain_need_sweep.json').write_text(json.dumps(rows, indent=1))
    for mass in (True, False):
        for s in ('breath_rib', 'breath_fast'):
            for ph in (14., 28.):
                q = [r for r in rows if r['mass'] == mass and r['scene'] == s and r['phase'] == ph]
                line = []
                for f in (None,)+FACTORS:
                    v = [r for r in q if r['factor'] == f]
                    line.append(f"{'own' if f is None else f}: {np.mean([r['L'] for r in v]):.3f}/{np.median([r['L'] for r in v]):.3f}/{max(r['peak'] for r in v):.0f}")
                own = [r['own_ratio'] for r in q if r['factor'] is None and r['own_ratio']]
                print(f"mass={mass} {s} {ph:.0f} deg own k/|K| {np.median(own):.2f} | "+' | '.join(line), flush=True)
    print('seconds', round(time.perf_counter()-start), flush=True)


if __name__ == '__main__':
    main()
