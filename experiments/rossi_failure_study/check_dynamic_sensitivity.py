"""How much do the unsettled parts of the dynamic truth matter to a closed-loop run? One change at a time, same scenes, L only; not a comparison of forms."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
CASES = [(c, t, 'nominal', i) for c, t in (('breath', 'relaxed'), ('rest', 'relaxed'), ('breath', 'active'), ('shift_fast', 'relaxed'), ('breath', 'mixed'), ('rib', 'relaxed')) for i in range(4)]
CHANGES = {
    'white, 14 deg': dict(dynamics='neural_white', phase=14.),
    'white, 14 deg, no slow memory': dict(dynamics='neural_white', phase=14., slow_scale=0.),
    'white, 14 deg, full slow memory': dict(dynamics='neural_white', phase=14., slow_scale=1.2),
    'white, 14 deg, long settled': dict(dynamics='neural_white', phase=14., load_wait=60.),
    'white, 14 deg, just arrived': dict(dynamics='neural_white', phase=14., load_wait=2.),
    'grey, 14 deg': dict(dynamics='neural_grey', phase=14.),
    'viscous, 14 deg': dict(dynamics='viscous', phase=14.),
    'white, 8 deg': dict(dynamics='neural_white', phase=8.),
    'white, 28 deg': dict(dynamics='neural_white', phase=28.),
    'viscous, 28 deg': dict(dynamics='viscous', phase=28.),
    'white, 38 deg (stress: wave radiation at probe scale, not measured)': dict(dynamics='neural_white', phase=38.),
    'viscous, 38 deg (stress: wave radiation at probe scale, not measured)': dict(dynamics='viscous', phase=38.),
    'earlier hc, 11 deg': dict(dynamics='hc', phase=11.),
    'earlier kb, 11 deg': dict(dynamics='kb', phase=11.),
    'earlier power, 11 deg': dict(dynamics='power', phase=11.),
}


def job(a):
    case, change = a
    c, t, g, i = case
    sc = S3.cell_scene3(c, t, g, i, base=494000, overrides=dict(force=8., **CHANGES[change]))
    pol = M.build('L', sc)
    log = P3.episode(sc, pol, 960+i)
    m = M.metric(log, pol, sc); m.update(case=case, change=change, equilibrium=sc.get('equilibrium', 1.))
    return m


def main():
    jobs = [(c, n) for c in CASES for n in CHANGES]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(job, jobs, chunksize=1))
    base = {tuple(r['case']): r for r in rows if r['change'] == 'white, 14 deg'}
    print('change | axial RMSE median/mean | peak>15 | peak>36.5 | numeric | accept | gain error | equilibrium | depth mm | L against white 14')
    for n in CHANGES:
        q = [r for r in rows if r['change'] == n]; ok = [r for r in q if not r['failed']]
        d = np.mean([np.log1p(r['J'])-np.log1p(base[tuple(r['case'])]['J']) for r in q])
        print(f"{n:70s} {np.median([r['rmse'] for r in ok]):.3f}/{np.mean([r['rmse'] for r in ok]):.3f}  {sum(r['peak']>15 for r in ok):2d} {sum(r['peak']>36.5 for r in ok):2d} {sum(r['failed'] for r in q):2d}  {np.mean([r['accept'] for r in ok]):.2f}  {np.nanmedian([abs(r['est_abs']) for r in ok]):.3f}  {np.mean([r['equilibrium'] for r in q]):.2f} {1e3*np.median([r['depth'] for r in ok]):5.1f}  {d:+.3f}")
    (HERE/'results'/'dynamic_sensitivity.json').write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
