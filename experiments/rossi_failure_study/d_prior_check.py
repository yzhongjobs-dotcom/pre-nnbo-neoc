"""A d treatment against D⁴ (break: one structural break; huber c: heavy-tailed prior), in the full Method: module-level windows by kind of moment, then the closed loop on the d_compare test scenes.

Pass: no clear loss in breathing and slow motion, a real gain in the fast-change scenes, flinch not worse. Truth only for scoring.
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import contact_candidates as C
import d_break as B
import d_robust as R
import d_compare as D
import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
LN12 = np.log(1.2)


def _job(a):
    method, cond, tissue, dyn, ph, i = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=ph))
    if method == 'D4':
        pol = D.make('D4', (0, 0, 0), sc)
    elif method.startswith('huber'):
        pol = R.build(sc, c=float(method.split()[-1]))
    else:
        pol = B.build(sc, penalty=float(method.split()[-1]))
    log = P3.episode(sc, pol, 900+i)
    m = M.metric(log, pol, sc); m.update(method=method, cond=cond, tissue=tissue, dynamics=dyn, phase=ph, i=i)
    if log is None:
        return m
    t = log[:, 0]; k = log[:, 21]
    ev = np.zeros(len(t), bool)
    for e in sc['event_info']:
        ev |= (t >= e['start']) & (t <= e['end']+.3)
    changing = np.convolve(np.abs(np.gradient(np.log(np.maximum(k, 1.)), t)) > 2., np.ones(30), 'same') > 0
    w = []
    for tk, kh, sig, chi, ok in pol.accept_trace:
        j = int(np.clip(np.searchsorted(t, tk), 0, len(t)-1))
        if t[j] < 1. or k[j] <= 200. or not np.isfinite(kh) or kh <= 0.:
            continue
        w.append(('event' if ev[j] else ('changing' if changing[j] else 'quiet'), float(np.log(kh/k[j])), bool(ok)))
    m['windows'] = w
    return m


def main():
    methods = ['D4']+(sys.argv[1].split(',') if len(sys.argv) > 1 else ['break 1'])
    more = len(sys.argv) > 2 and sys.argv[2] == 'more'      # six new seeds on the fast-change cells, to judge tails
    cells = [ct for ct in D.TEST if ct[0] in ('flinch', 'cough', 'rib', 'combo')] if more else D.TEST
    seeds = range(3, 9) if more else range(3)
    jobs = [(mt, c, t, dyn, ph, i) for mt in methods for c, t in cells for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in seeds]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    print(len(rows), 'runs', round(time.time()-t0), 's')
    print('module: kind, accepted share, |ln k/k| median of accepted, share of accepted beyond ln 1.2')
    for mt in methods:
        for kind in ('quiet', 'event', 'changing'):
            w = [x for r in rows if r['method'] == mt for x in r.get('windows', []) if x[0] == kind]
            ok = np.array([x[2] for x in w]); e = np.abs([x[1] for x in w])[ok]
            print(f"  {mt:10s} {kind:9s} accepted {ok.mean():.2f}  error {np.median(e):.3f}  beyond 1.2 {np.mean(e > LN12):.3f}")
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    base = {key(r): r for r in rows if r['method'] == 'D4'}
    print('closed loop, against D⁴, mean difference in L [95% CI], negative is better')
    for mt in methods:
        q = [r for r in rows if r['method'] == mt]
        print(f"  {mt:10s} {D.line(q)}")
        if mt == 'D4':
            continue
        def diff(sel):
            d = np.array([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in q if sel(r)])
            ci = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
            return f"{d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]"
        print(f"     all {diff(lambda r: True)}")
        for c, t in cells:
            qq = [r for r in q if r['cond'] == c]
            print(f"     {c:10s} {diff(lambda r, c=c: r['cond'] == c)}  rmse {np.median([base[key(r)]['rmse'] for r in qq]):.3f} -> {np.median([r['rmse'] for r in qq]):.3f}  event {np.nanmedian([base[key(r)]['event_rmse'] for r in qq]) if c not in ('breath', 'rib') else float('nan'):.3f} -> {np.nanmedian([r['event_rmse'] for r in qq]) if c not in ('breath', 'rib') else float('nan'):.3f}  peak {max(base[key(r)]['peak'] for r in qq):.1f} -> {max(r['peak'] for r in qq):.1f}")
    (HERE/'results'/(('d_break_check' if any(m.startswith('break') for m in methods) else 'd_prior_check')+('_more' if more else '')+'.json')).write_text(json.dumps([{k: v for k, v in r.items() if k != 'windows'} for r in rows], indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
