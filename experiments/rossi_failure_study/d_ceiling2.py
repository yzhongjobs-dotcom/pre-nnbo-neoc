"""A perfect estimator inside the unchanged Method: only the window estimate is replaced by the truth, all gating, conservatism and guards run as usual.

own      the D⁴ host as tuned
probe    window estimate = true gain at the probe frequency, |K(20.8 Hz)| with the series path
low      window estimate = true gain at 1 Hz with the series path, the band the force loop and the feed-forward work in
available  window estimate = true probe gain, also precise and consistent inside abrupt changes: the value of a d treatment that never loses a window
The estimator's own sigma and chi are kept, so acceptance is unchanged. Diagnostic only.
"""
import json
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import d_compare as D
import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
MODES = tuple(sys.argv[2].split(',')) if len(sys.argv) > 2 else ('own', 'probe', 'low', 'available')


def _job(a):
    cfg, cond, tissue, dyn, phase, i, mode = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=phase))
    pol = D.make('D4', cfg, sc)
    if mode != 'own':
        est = pol.estimator; original = est.estimate
        def estimate(o, n, policy, original=original):
            k, sig, chi = original(o, n, policy)
            truth = getattr(o, 'truth_k1' if mode == 'low' else 'truth_k', 0.)
            if np.isfinite(k) and truth > 200.:
                if mode == 'available':
                    return truth, .05*truth, 0.          # also usable inside abrupt changes: precise and consistent
                sig = sig*truth/k if k > 0 else sig       # same relative precision
                k = truth
            return k, sig, chi
        est.estimate = estimate
    log = P3.episode(sc, pol, 900+i, diagnostic=mode != 'own')
    m = M.metric(log, pol, sc); m.update(cond=cond, tissue=tissue, dynamics=dyn, phase=phase, i=i, mode=mode)
    if log is not None:
        warm = log[:, 0] > 1.
        m['k1_over_probe'] = float(np.median(log[warm, 9]/np.maximum(log[warm, 21], 1.)))
    return m


def main():
    chosen = pickle.loads(Path(sys.argv[1]).read_bytes())['chosen']['D4']
    jobs = [(chosen, c, t, dyn, ph, i, mode) for mode in MODES for c, t in D.TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    base = {key(r): r for r in rows if r['mode'] == 'own'}
    print('true 1 Hz gain over true probe-frequency gain, median per scene:', np.round(np.percentile([r['k1_over_probe'] for r in rows if r['mode'] == 'own'], [10, 50, 90]), 3))
    for mode in MODES:
        q = [r for r in rows if r['mode'] == mode]
        d = np.array([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in q])
        ci = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
        print(f"  {mode:5s} {D.line(q)}  vs own {d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]")
        for c, t in D.TEST:
            qq = [r for r in q if r['cond'] == c]
            print(f"       {c:10s} rmse {np.median([r['rmse'] for r in qq]):.3f} event {np.nanmedian([r['event_rmse'] for r in qq]):.3f} peak max {max(r['peak'] for r in qq):.1f} L diff {np.mean([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in qq]):+.3f}")
        for ph in (14., 28.):
            qq = [r for r in q if r['phase'] == ph]
            print(f"       phase {ph:.0f}: L diff {np.mean([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in qq]):+.3f}")
    (HERE/'results'/('d_ceiling2.json' if len(MODES) == 4 else 'd_ceiling2_'+'_'.join(MODES)+'.json')).write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
