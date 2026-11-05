"""Can tuning the motion tracker rescue the linear relation at high phase, and what does each tuned form cost? Same budget for both, chosen on development seeds, tested on new ones."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
BASE = 480000
CONFIGS = [dict(r_scale=r, q_maneuver=q) for r in (.2264, 1., 4., 16.) for q in (3.335, .8)]


def _job(a):
    name, c, cond, tissue, dyn, phase, i, base = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=base, overrides=dict(slide='long', dynamics=dyn, phase=phase))
    pol = M.build(name, sc, CONFIGS[c])
    log = P3.episode(sc, pol, 800+i)
    m = M.metric(log, pol, sc); m.update(name=name, c=c, cond=cond, tissue=tissue, dyn=dyn, phase=phase, i=i)
    return m


def main():
    names = ('L', 'H')
    dev = [(n, c, cond, t, dyn, 23., i, BASE) for n in names for c in range(len(CONFIGS)) for cond in ('breath', 'shift_fast') for t in ('relaxed', 'active') for dyn in ('hc', 'power') for i in range(2)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, dev, chunksize=1))
    best = {}
    score = lambda q: float(np.mean([np.log1p(r['J']) for r in q]))
    print('development, high phase: mean L, axial RMSE median, peak max')
    for n in names:
        for c, cfg in enumerate(CONFIGS):
            q = [r for r in rows if r['name'] == n and r['c'] == c]
            print(f"  {n} {cfg} L {score(q):.3f} rmse {np.median([r['rmse'] for r in q]):.3f} peak {max(r['peak'] for r in q):.1f}")
            if n not in best or score(q) < best[n][0]:
                best[n] = (score(q), c)
    print('chosen', {n: CONFIGS[best[n][1]] for n in names})
    test = list({(n, c, cond, t, dyn, ph, i, BASE+500): None for n in names for c in (0, best[n][1]) for cond in ('breath', 'shift_fast', 'flinch', 'cough') for t in ('relaxed', 'active', 'rib', 'mixed') for dyn, ph in (('power', 11.), ('kb', 11.), ('hc', 23.), ('power', 23.)) for i in range(3)})
    print(len(test), 'test runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows2 = list(ex.map(_job, test, chunksize=1))
    print('\nnew seeds: mean L, axial RMSE median/mean, event RMSE median, peak>15, peak>36.5, long loss')
    for n in names:
        for c in sorted({0, best[n][1]}):
            for ph in (11., 23.):
                q = [r for r in rows2 if r['name'] == n and r['c'] == c and r['phase'] == ph]
                print(f"  {n} {'default' if c == 0 else 'tuned  '} phase {ph:.0f}: L {score(q):.3f} rmse {np.median([r['rmse'] for r in q]):.3f}/{np.mean([r['rmse'] for r in q]):.3f} event {np.nanmedian([r['event_rmse'] for r in q]):.3f} >15 {sum(r['peak']>15 for r in q)} >36.5 {sum(r['peak']>36.5 for r in q)} long {sum(r['long_loss'] for r in q)}")
    k = lambda r: (r['cond'], r['tissue'], r['dyn'], r['phase'], r['i'])
    for ph in (11., 23.):
        A = {k(r): r for r in rows2 if r['name'] == 'H' and r['c'] == best['H'][1] and r['phase'] == ph}
        B = {k(r): r for r in rows2 if r['name'] == 'L' and r['c'] == best['L'][1] and r['phase'] == ph}
        d = np.array([np.log1p(A[x]['J'])-np.log1p(B[x]['J']) for x in A])
        c = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
        print(f'tuned H - tuned L, phase {ph:.0f}: {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] n={len(d)}')
    (HERE/'results'/'select3_tune.json').write_text(json.dumps(dict(dev=rows, test=rows2, chosen={n: CONFIGS[best[n][1]] for n in names}), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
