"""What makes a probe crossing thin-cover ribs lose force control? One change at a time on scenes that lose it: the setting, the gain the Method uses, the actuator, the feed-forward."""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
NAMES = ('L', 'H')
BASE = dict(cover=.007, slide='short')
CASES = [(c, t, g, i) for g in ('nominal', 'stop_go', 'reversal') for c, t in (('rest', 'rib'), ('breath', 'rib'), ('shift_slow', 'rib'), ('rib', 'relaxed'), ('rib', 'mixed'), ('rib', 'active')) for i in range(2)]
CHANGES = {
    'baseline': dict(),
    'true gain, measurement': dict(oracle='meas'),
    'true gain, measurement and feedback': dict(oracle='both'),
    'no feed-forward': dict(ff=False),
    'no feed-forward, true feedback gain': dict(ff=False, oracle='fb'),
    'actuator speed x4': dict(vmax=1.),
    'scan speed /3': dict(scene=dict(speed_div=3.)),
    'cover 18 mm': dict(scene=dict(cover=.018)),
    'ribs 3x wider edge': dict(edge=3.),
    'open feed-forward at gain jump 1.5': dict(jump=1.5),
    'open feed-forward at gain jump 2': dict(jump=2.),
}


def _job(a):
    case, name, change = a
    c, t, g, i = case
    spec = CHANGES[change]
    o = dict(BASE, **{k: v for k, v in spec.get('scene', {}).items() if k != 'speed_div'})
    sc = S3.cell_scene3(c, t, g, i, base=452000, overrides=o)
    if 'speed_div' in spec.get('scene', {}):
        sc = S3.cell_scene3(c, t, g, i, base=452000, overrides=dict(o, speed=sc['v']/spec['scene']['speed_div']))
    if 'edge' in spec and len(sc['ribs']):
        sc['ribs'] = sc['ribs'].copy(); sc['ribs'][:, 2] *= spec['edge']
    pol = M.build(name, sc)
    if 'jump' in spec:
        pol.params['ff_jump'] = spec['jump']
    if 'ff' in spec:
        pol.params['feedforward'] = spec['ff']
    if 'oracle' in spec:
        original = pol.stiffness; mode = spec['oracle']
        def stiffness(o_, n, original=original, pol=pol, mode=mode):
            kb, ok = original(o_, n)
            k = getattr(o_, 'truth_k', 0.)
            if k > 200.:
                if mode in ('meas', 'both'):
                    pol.k_good = k
                if mode in ('both', 'fb'):
                    return k, True
            return kb, ok
        pol.stiffness = stiffness
    platform = dict(P3.PLATFORM)
    if 'vmax' in spec:
        platform['VMAX'] = spec['vmax']
    log = P3.episode(sc, pol, 200+i, platform=platform)
    m = M.metric(log, pol, sc)
    tr = np.array(pol.accept_trace) if pol.accept_trace else np.zeros((0, 5))
    lag = np.nan
    if log is not None and len(tr):
        # how far the gain falls behind the true stiffness in the half second before the worst moment
        j = int(np.argmax(log[:, 18])); t0 = log[j, 0]
        ok = (tr[:, 4] > 0) & (tr[:, 0] < t0)
        if ok.any():
            lag = float(log[j, 21]/tr[ok][-1, 1])
    cmd = np.abs(np.diff(log[:, 11]))/.008 if log is not None else np.zeros(1)
    m.update(case=case, name=name, change=change, lag=lag, sat=float(np.mean(cmd > .98*platform['VMAX'])), k_peak=float(log[:, 21].max()) if log is not None else np.nan)
    return m


def main():
    jobs = [(c, n, ch) for c in CASES for n in NAMES for ch in CHANGES]
    print(len(jobs), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    (HERE/'results'/'failure_causes.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    base = {(tuple(r['case']), r['name']): r for r in rows if r['change'] == 'baseline'}
    bad = {k for k, r in base.items() if r['failed'] or r['peak'] > 36.5}
    print('scenes', len(CASES), 'baseline force above 36.5 N or numerically failed:', len(bad), 'of', len(base))
    print('\nchange | all scenes: numeric fail, peak>36.5, peak>15, median peak, rmse, time>15 N, saturated ticks | scenes bad at baseline: still bad, median peak')
    for ch in CHANGES:
        q = [r for r in rows if r['change'] == ch]
        qb = [r for r in q if (tuple(r['case']), r['name']) in bad]
        ok = [r for r in q if not r['failed']]
        print(f"{ch:38s} | {sum(r['failed'] for r in q):2d} {sum(r['failed'] or r['peak']>36.5 for r in q):2d} {sum(r['failed'] or r['peak']>15 for r in q):2d} {np.median([r['peak'] for r in ok]):6.1f} {np.median([r['rmse'] for r in ok]):.2f} {np.mean([r['over_15'] for r in ok]):.3f} {np.mean([r['sat'] for r in ok]):.3f} | {sum(r['failed'] or r['peak']>36.5 for r in qb):2d}/{len(qb)} {np.nanmedian([r['peak'] for r in qb]) if qb else float('nan'):6.1f}")
    print('\nbad baseline scenes: stiffness at the worst moment over last accepted gain, median', np.nanmedian([base[k]['lag'] for k in bad]) if bad else None, ' peak true stiffness N/mm', np.nanmedian([base[k]['k_peak'] for k in bad])/1e3 if bad else None)
    print('by name, baseline', {n: sum(1 for k in bad if k[1] == n) for n in NAMES})


if __name__ == '__main__':
    main()
