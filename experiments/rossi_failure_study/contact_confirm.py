"""Confirmation of the contact-model selection on unused seeds, and where the linear relation stops holding.

confirm  the top two of contact_select plus references, six regimes plus 'power' (the phase of 'rate' from the weak power law's own branches,
         so the high-phase case is not tested only with the Hunt-Crossley law), base 430000.
sweep    the Hunt-Crossley coefficient lam of the plant over a grid on the same scenes, with the plant's phase measured per scene.
report   tables from the raw rows of both.
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import contact_candidates as C
import scenes2 as SC2
import sweeps as S
from check_candidates import plant_phase
from contact_select import _force_metrics

BASE = 430000
REGIMES = ('nominal', 'stop_go', 'reversal', 'curved', 'rate', 'memory', 'power')
NAMES = ('L', 'W', 'H', 'HA', 'P5:L', 'P5:W', 'P5:H', 'P5:HA')
SWEEP = ('L', 'H', 'P5:L', 'P5:H')
LAMS = (.2, .4, .6, .8, 1., 1.5, 2., 3.)


def _seed(base, regime, cond, tissue, i):
    return base+SC2.REGIME_CODE[regime]*100000+sum(ord(x) for x in cond+tissue)*100+i


def _run(name, sc, seed):
    log, pol, _ = C.run_one(name, sc, seed=seed, sensor=True)
    m = S._metric(log, pol, sc, 20.8)
    m.update(_force_metrics(log)); m['phi'] = float(getattr(pol, 'phi', 0.)); m['phi_count'] = int(getattr(pol, 'phi_count', 0))
    return m


def _confirm(arg):
    regime, (cond, tissue), i = arg
    sc = SC2.cell_scene2(cond, tissue, regime, i, base=BASE)
    seed = _seed(BASE, regime, cond, tissue, i)
    return [dict(regime=regime, cond=cond, tissue=tissue, i=i, candidate=n, metrics=_run(n, sc, seed)) for n in NAMES]


def _sweep(arg):
    lam, (cond, tissue), i = arg
    sc = SC2.cell_scene2(cond, tissue, 'rate_sweep', i, base=BASE)
    sc['lam'] = lam
    seed = _seed(BASE, 'rate_sweep', cond, tissue, i)
    true, gain = plant_phase(sc)
    return [dict(lam=lam, cond=cond, tissue=tissue, i=i, candidate=n, phase=true, gain=gain, metrics=_run(n, sc, seed)) for n in SWEEP]


L = lambda m: math.log1p(min(m['J'], 1000))
fail = lambda m: bool(m['failed']) or m['peak'] > 50


def boot(d):
    d = np.asarray(d, float)
    return np.percentile(d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1), [2.5, 97.5])


def report():
    rows = json.loads((HERE/'results'/'contact_confirm.json').read_text(encoding='utf-8'))
    T = {}
    for r in rows:
        T.setdefault(r['candidate'], {})[(r['regime'], r['cond'], r['tissue'], r['i'])] = r['metrics']
    print('rows', len(rows))
    print('\nmean L per regime; mean of the six selection regimes')
    print('name     '+' '.join(f'{g:>9s}' for g in REGIMES)+'      six')
    for n in NAMES:
        v = [float(np.mean([L(m) for k, m in T[n].items() if k[0] == g])) for g in REGIMES]
        print(f'{n:8s} '+' '.join(f'{x:9.3f}' for x in v)+f' {np.mean(v[:6]):8.3f}')
    print('\npaired differences, six regimes [95%]; * per regime excludes 0; power listed last')
    for a, b in (('H', 'L'), ('HA', 'L'), ('HA', 'H'), ('W', 'L'), ('H', 'W'), ('P5:H', 'P5:L'), ('P5:HA', 'P5:L'), ('P5:HA', 'P5:H'), ('P5:W', 'P5:L'), ('P5:H', 'P5:W'), ('P5:H', 'H'), ('P5:L', 'L')):
        cells = []
        for g in REGIMES:
            d = np.array([L(m)-L(T[b][k]) for k, m in T[a].items() if k[0] == g]); c = boot(d)
            cells.append(f"{d.mean():+.3f}{'*' if c[0] > 0 or c[1] < 0 else ' '}")
        d = np.array([L(m)-L(T[b][k]) for k, m in T[a].items() if k[0] != 'power']); c = boot(d)
        e = np.array([L(m)-L(T[b][k]) for k, m in T[a].items() if k[0] not in ('power', 'rate')]); f = boot(e)
        print(f'{a:6s}-{b:6s}'+' '.join(f'{x:>8s}' for x in cells)+f'  six {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]  no rate {e.mean():+.3f} [{f[0]:+.3f},{f[1]:+.3f}]')
    print('\nsix regimes pooled, then power: Fn RMSE median/P90, event median, runaways, long losses')
    for n in NAMES:
        for lab, keep in (('six', lambda k: k[0] != 'power'), ('power', lambda k: k[0] == 'power')):
            rr = [m for k, m in T[n].items() if keep(k)]
            print(f"{n:8s}{lab:6s} Fn {np.nanmedian([m['fn_rmse'] for m in rr]):.3f}/{np.nanpercentile([m['fn_rmse'] for m in rr], 90):.3f} event {np.nanmedian([m['event_rmse'] for m in rr]):.3f} fail {sum(fail(m) for m in rr)} long {sum(m['long_loss'] for m in rr)} of {len(rr)}")
    print('\nonline phase, degrees, median by regime')
    for n in ('H', 'P5:H', 'W', 'P5:W'):
        print(f'{n:6s}', ' '.join(f"{g}:{math.degrees(np.median([m['phi'] for k, m in T[n].items() if k[0] == g])):.1f}" for g in REGIMES))
    path = HERE/'results'/'contact_lambda_sweep.json'
    if path.exists():
        rows = json.loads(path.read_text(encoding='utf-8'))
        print('\nlam sweep: plant phase deg median [P10,P90], then mean L, runaways of each, share of scenes with J>1')
        print('lam   phase              '+' '.join(f'{n:>16s}' for n in SWEEP))
        for lam in LAMS:
            q = [r for r in rows if r['lam'] == lam]
            ph = np.degrees([r['phase'] for r in q if r['candidate'] == 'L'])
            cells = []
            for n in SWEEP:
                mm = [r['metrics'] for r in q if r['candidate'] == n]
                cells.append(f"{np.mean([L(m) for m in mm]):.2f}/{sum(fail(m) for m in mm)}/{np.mean([m['J'] > 1 for m in mm]):.2f}")
            print(f"{lam:4.1f}  {np.median(ph):5.1f} [{np.percentile(ph, 10):5.1f},{np.percentile(ph, 90):5.1f}]  "+' '.join(f'{x:>16s}' for x in cells))
        print('\nby measured phase band, pooled over lam: n, mean L and share J>1')
        edges = [0, 8, 12, 16, 20, 25, 30, 90]
        for a, b in zip(edges[:-1], edges[1:]):
            q = [r for r in rows if a <= math.degrees(r['phase']) < b]
            if q:
                print(f'{a:3d}-{b:3d} deg n {sum(r["candidate"] == "L" for r in q):4d} '+' '.join(f"{n}:{np.mean([L(r['metrics']) for r in q if r['candidate'] == n]):.2f}/{np.mean([r['metrics']['J'] > 1 for r in q if r['candidate'] == n]):.2f}" for n in SWEEP))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'report'
    if mode == 'report':
        return report()
    smoke = len(sys.argv) > 2 and sys.argv[2] == 'smoke'
    cells = [('rest', 'soft'), ('breath', 'mid')] if smoke else SC2.CELLS2
    reps = 1 if smoke else 3
    if mode == 'confirm':
        jobs = [(g, c, i) for g in REGIMES for c in cells for i in range(reps)]; fn = _confirm; out = 'contact_confirm.json'
    else:
        jobs = [(lam, c, 0) for lam in LAMS for c in cells]; fn = _sweep; out = 'contact_lambda_sweep.json'
    print(mode, len(jobs), 'scenes', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = [r for block in ex.map(fn, jobs, chunksize=1) for r in block]
    if not smoke:
        (HERE/'results'/out).write_text(json.dumps(rows, indent=1), encoding='utf-8')
    else:
        for r in rows:
            print(r.get('regime', r.get('lam')), r['cond'], r['tissue'], r['candidate'], round(r['metrics']['J'], 3), round(math.degrees(r['metrics']['phi']), 1), round(math.degrees(r.get('phase', 0.)), 1))


if __name__ == '__main__':
    main()
