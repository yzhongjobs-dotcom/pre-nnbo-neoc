"""Perfect-estimator test in scenes A and B, with the logged truth and with the gain the sampled loop actually sees.

Open-loop replays show that the y->F gain at the 8 ms ticks exceeds the logged |K(20.8 Hz)| by about (1 + beta*tool_mass/kr), beta ~ 1.17e5 1/s^2:
the 125 Hz command hold puts a 104 Hz image into the probe's acceleration, the arm compliance turns it into depth, and tick sampling aliases it to 20.8 Hz.
Modes, all inside the frozen overlay protocol and the adopted Huber-D4 Method, gates unchanged:
  own       the estimator as is
  logged    window estimate = logged |K(20.8 Hz)|
  sampled   window estimate = logged |K(20.8 Hz)|*(1 + beta*tool_mass/kr)
With --no-tool-mass the probe mass beyond the sensor is set to zero, which removes that path and leaves the logged truth correct.
Diagnostic only. Usage: python -X utf8 gain_truth_ceiling.py [--no-tool-mass]
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import overlay_compare as O
import method3 as M
import physics_plant3 as P3

B = O.B
HERE = Path(__file__).resolve().parent
BETA = 1.17e5
NOMASS = '--no-tool-mass' in sys.argv
MODES = ('own', 'logged') if NOMASS else ('own', 'logged', 'sampled')


def job(a):
    cell, mode = a
    O.setup()
    cfg = O.chosen()['Huber-D4']
    sc = B.scene(cell, O.BASE, 0.)
    if NOMASS:
        sc['tool_mass'] = 0.; sc['tool_mass_hat'] = 0.
    pol = B.build('Huber-D4', cfg, sc)
    factor = 1.+BETA*sc['tool_mass']/sc['kr']
    if mode != 'own':
        est = pol.estimator; original = est.estimate
        def estimate(o, n, policy, original=original):
            k, sig, chi = original(o, n, policy)
            truth = getattr(o, 'truth_k', 0.)*(factor if mode == 'sampled' else 1.)
            if np.isfinite(k) and truth > 200.:
                sig = sig*truth/k if k > 0 else sig
                k = truth
            return k, sig, chi
        est.estimate = estimate
    log = P3.episode(sc, pol, O.sensor_seed(cell), diagnostic=mode != 'own')
    m = M.metric(log, pol, sc)
    c, t, dy, ph, i = cell
    return dict(scene=c, dynamics=dy, phase=ph, i=i, mode=mode, L=float(np.log1p(m['J'])), J=float(m['J']), rmse=float(m['rmse']), peak=float(m['peak']), factor=float(factor))


def boot(d, rng, clusters):
    keys = sorted(set(clusters))
    idx = {k: np.where(np.array(clusters) == k)[0] for k in keys}
    means = []
    for _ in range(5000):
        pick = rng.choice(len(keys), len(keys))
        means.append(np.mean(np.concatenate([d[idx[keys[p]]] for p in pick])))
    return np.percentile(means, [2.5, 97.5])


def main():
    start = time.perf_counter()
    cells = [c for c in O.cells() if c[0] in ('breath_rib', 'breath_fast')]
    jobs = [(c, mode) for mode in MODES for c in cells]
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(job, jobs, chunksize=1))
    rng = np.random.default_rng(0)
    key = lambda r: (r['scene'], r['dynamics'], r['phase'], r['i'])
    own = {key(r): r for r in rows if r['mode'] == 'own'}
    summary = {}
    for scene in ('breath_rib', 'breath_fast'):
        for mode in MODES:
            q = [r for r in rows if r['scene'] == scene and r['mode'] == mode]
            d = np.array([r['L']-own[key(r)]['L'] for r in q])
            ci = boot(d, rng, [r['i'] for r in q]) if mode != 'own' else [0., 0.]
            s = dict(n=len(q), L_mean=float(np.mean([r['L'] for r in q])), L_median=float(np.median([r['L'] for r in q])), rmse_median=float(np.median([r['rmse'] for r in q])), peak_max=float(np.max([r['peak'] for r in q])), over15=int(sum(r['peak'] > 15. for r in q)), dL=float(d.mean()), ci=[float(ci[0]), float(ci[1])])
            summary[f'{scene} {mode}'] = s
            print(f'{scene} {mode}', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()}, flush=True)
    print('factor median', np.median([r['factor'] for r in rows]), 'seconds', round(time.perf_counter()-start), flush=True)
    (HERE/'results'/('gain_truth_ceiling_no_mass.json' if NOMASS else 'gain_truth_ceiling.json')).write_text(json.dumps(dict(summary=summary, rows=rows, beta=BETA), indent=1))


if __name__ == '__main__':
    main()
