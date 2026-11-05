"""Close-out of the selection: can tuning rescue the linear relation, does the other host agree, and does the choice hold in the full Method across tissues and events."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
BASE = 470000
Q = dict(slide='long', pose_error=0., rho=0., probe='linear', contact_law='halfspace', force=8.)
CONFIGS = [dict(r_scale=r, D1=d) for r in (.226, .6, 1.5, 4.) for d in (.0365, .07)]


def _tune(a):
    name, c, cond, dyn, phase, i, base = a
    sc = S3.cell_scene3(cond, 'relaxed', 'nominal', i, base=base, overrides=dict(Q, dynamics=dyn, phase=phase))
    pol = M.build(name, sc); pol.params.update(CONFIGS[c])
    log = P3.episode(sc, pol, 600+i)
    m = M.metric(log, pol, sc); m.update(name=name, c=c, cond=cond, dyn=dyn, phase=phase, i=i)
    return m


def _full(a):
    name, cond, tissue, phase, dyn, i = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=BASE+1000, overrides=dict(dynamics=dyn, phase=phase))
    log, pol, _ = M.run(name, sc, seed=700+i)
    m = M.metric(log, pol, sc); m.update(name=name, cond=cond, tissue=tissue, phase=phase, dyn=dyn, i=i, force=sc['FD'], probe=sc['probe'])
    return m


def main():
    names = ('L', 'H')
    dev = [(n, c, cond, dyn, 23., i, BASE) for n in names for c in range(len(CONFIGS)) for cond in ('rest', 'breath') for dyn in ('hc', 'power') for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_tune, dev, chunksize=1))
    print('tuning at the high phase, development seeds: config, mean axial RMSE, max peak')
    best = {}
    for n in names:
        for c, cfg in enumerate(CONFIGS):
            q = [r for r in rows if r['name'] == n and r['c'] == c]
            score = float(np.mean([r['rmse'] for r in q]))
            print(f"  {n} {cfg} rmse {score:.3f} peak {max(r['peak'] for r in q):.1f}")
            if n not in best or score < best[n][0]:
                best[n] = (score, c)
    print('best', {n: CONFIGS[best[n][1]] for n in names})
    test = [(n, c, cond, dyn, ph, i, BASE+500) for n in names for c in (0, best[n][1]) for cond in ('rest', 'breath', 'shift_fast', 'cough') for dyn in ('hc', 'power', 'kb') for ph in (11., 23.) for i in range(3)]
    test = list({t: None for t in test})
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows2 = list(ex.map(_tune, test, chunksize=1))
    print('\nnew seeds: axial RMSE mean, event-free median, peak max, time over FD+3 mean, by phase')
    for n in names:
        for c in sorted({0, best[n][1]}):
            for ph in (11., 23.):
                q = [r for r in rows2 if r['name'] == n and r['c'] == c and r['phase'] == ph]
                print(f"  {n} {'default' if c == 0 else 'tuned  '} phase {ph:.0f}: mean {np.mean([r['rmse'] for r in q]):.3f} median {np.median([r['rmse'] for r in q]):.3f} peak {max(r['peak'] for r in q):.1f} over {np.mean([r['over_band'] for r in q]):.3f} long {sum(r['long_loss'] for r in q)}")
    # the full Method, both hosts, default parameters
    full = [(n, cond, t, ph, dyn, i) for n in ('L', 'H', 'P5:L', 'P5:H') for t in ('relaxed', 'active', 'rib', 'mixed') for cond in ('breath', 'shift_fast', 'flinch', 'cough') for ph, dyn in ((11., 'power'), (11., 'kb'), (23., 'power'), (23., 'hc')) for i in range(3)]
    print('\nfull Method', len(full), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows3 = list(ex.map(_full, full, chunksize=1))
    print('name phase: axial RMSE median/mean, normal median, event median, peak>15 scenes, peak>36.5 scenes, numeric, long loss, accept, gain error, call p50/p99')
    for ph in (11., 23.):
        for n in ('L', 'H', 'P5:L', 'P5:H'):
            q = [r for r in rows3 if r['name'] == n and r['phase'] == ph]; ok = [r for r in q if not r['failed']]
            print(f"  {n:5s} {ph:.0f}: {np.median([r['rmse'] for r in ok]):.3f}/{np.mean([r['rmse'] for r in ok]):.3f} {np.median([r['normal_rmse'] for r in ok]):.3f} {np.nanmedian([r['event_rmse'] for r in ok]):.3f} {sum(r['peak']>15 for r in ok)} {sum(r['peak']>36.5 for r in ok)} {sum(r['failed'] for r in q)} {sum(r['long_loss'] for r in q)} {np.mean([r['accept'] for r in ok]):.2f} {np.nanmedian([abs(r['est_abs']) for r in ok]):.3f} {np.median([r['call_p50'] for r in ok]):.0f}/{np.median([r['call_p99'] for r in ok]):.0f}")
    def paired(a, b, sel):
        A = {(r['cond'], r['tissue'], r['phase'], r['dyn'], r['i']): r for r in rows3 if r['name'] == a and sel(r)}
        B = {(r['cond'], r['tissue'], r['phase'], r['dyn'], r['i']): r for r in rows3 if r['name'] == b and sel(r)}
        d = np.array([np.log1p(A[k]['J'])-np.log1p(B[k]['J']) for k in A])
        c = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
        return f'{d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] n={len(d)}'
    print('\npaired L = log(1+J) differences')
    for a, b in (('H', 'L'), ('P5:H', 'P5:L')):
        print(' ', a, '-', b, 'phase 11', paired(a, b, lambda r: r['phase'] == 11.), '| phase 23', paired(a, b, lambda r: r['phase'] == 23.))
    print('by tissue, H - L', {t: paired('H', 'L', lambda r: r['tissue'] == t) for t in ('relaxed', 'active', 'rib', 'mixed')})
    (HERE/'results'/'select3_check.json').write_text(json.dumps(dict(tuning=rows, test=rows2, full=rows3), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
