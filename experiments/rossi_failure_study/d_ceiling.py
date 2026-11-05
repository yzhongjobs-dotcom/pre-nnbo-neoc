"""How much is left for the handling of d? The tuned D⁴ host against the same Method given the true probe-frequency gain (diagnostic only), on the d_compare test scenes."""
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


def _job(a):
    cfg, cond, tissue, dyn, phase, i, mode = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=phase))
    pol = D.make('D4', cfg, sc)
    if mode != 'own':
        original = pol.stiffness
        def stiffness(o_, n, original=original, pol=pol):
            kb, ok = original(o_, n)
            k = getattr(o_, 'truth_k', 0.)
            if k > 200.:
                pol.k_good = k
                if mode == 'both':
                    return k, True
            return kb, ok
        pol.stiffness = stiffness
    log = P3.episode(sc, pol, 900+i, diagnostic=mode != 'own')
    m = M.metric(log, pol, sc); m.update(cond=cond, tissue=tissue, dynamics=dyn, phase=phase, i=i, mode=mode)
    return m


def main():
    chosen = pickle.loads(Path(sys.argv[1]).read_bytes())['chosen']['D4']
    jobs = [(chosen, c, t, dyn, ph, i, mode) for mode in ('meas', 'both') for c, t in D.TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    own = json.loads((HERE/'results'/'d_compare.json').read_text(encoding='utf-8'))['test']
    own = [dict(r, mode='own') for r in own if r['host'] == 'D4']
    rows += own
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    print('mode: own = tuned D⁴; meas = true gain in the force-to-motion conversion; both = true gain also in the feedback')
    for mode in ('own', 'meas', 'both'):
        print(' ', mode, D.line([r for r in rows if r['mode'] == mode]))
    base = {key(r): r for r in rows if r['mode'] == 'own'}
    for mode in ('meas', 'both'):
        for c, t in D.TEST:
            q = [r for r in rows if r['mode'] == mode and r['cond'] == c]
            d = [np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in q]
            print(f"  {mode} {c:10s} rmse own {np.median([base[key(r)]['rmse'] for r in q]):.3f} -> {np.median([r['rmse'] for r in q]):.3f}  L diff {np.mean(d):+.3f}  peak own {max(base[key(r)]['peak'] for r in q):.1f} -> {max(r['peak'] for r in q):.1f}")
    (HERE/'results'/'d_ceiling.json').write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
