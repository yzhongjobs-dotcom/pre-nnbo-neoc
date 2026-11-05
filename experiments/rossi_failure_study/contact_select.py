"""Smoke/full matrix runner for the contact-model selection study."""
import argparse
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import contact_candidates as C
import physics_plant2 as PP2
import scenes2 as SC2
import sweeps as S

BASE = 420000
REGIMES = ('nominal', 'stop_go', 'reversal', 'curved', 'rate', 'memory')
ONLY = ('W', 'S', 'WF', 'WA', 'WT', 'P5:W', 'P5:S')    # the forms whose response changed; the other rows of the matrix stay


def _force_metrics(log):
    if log is None:
        return dict(fn_rmse=np.nan, faxis_rmse=np.nan, fn_p90=np.nan, faxis_p90=np.nan)
    warm = log[:, 0] >= 1.
    fn = log[warm, 3] - 4.; fa = log[warm, 18] - 4.
    return dict(fn_rmse=float(np.sqrt(np.mean(fn * fn))), faxis_rmse=float(np.sqrt(np.mean(fa * fa))),
                fn_p90=float(np.percentile(abs(fn), 90)), faxis_p90=float(np.percentile(abs(fa), 90)))


def _job(arg):
    regime, cell, i, base = arg
    cond, tissue = cell
    sc = SC2.cell_scene2(cond, tissue, regime, i, base=base)
    seed = base + SC2.REGIME_CODE[regime] * 100000 + sum(ord(x) for x in cond + tissue) * 100 + i
    rows = []
    for name in ONLY or C.NAMES:
        log, pol, elapsed = C.run_one(name, sc, seed=seed, sensor=True)
        metric = S._metric(log, pol, sc, 20.8)
        metric.update(_force_metrics(log))
        metric['seconds'] = float(elapsed)
        rows.append(dict(regime=regime, cond=cond, tissue=tissue, i=i, seed=seed, candidate=name, metrics=metric))
    return rows


def _boot(x, seed=420991, n=2000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x): return [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    means = x[rng.integers(0, len(x), (n, len(x)))].mean(axis=1)
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def _summary(rows):
    def med(x, q=50):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        return float(np.percentile(x, q)) if len(x) else np.nan

    grouped = {}
    for regime in REGIMES:
        grouped[regime] = {}
        rr = [r for r in rows if r['regime'] == regime]
        ref = {(r['cond'], r['tissue'], r['i']): r['metrics'] for r in rr if r['candidate'] == 'L'}
        for name in C.NAMES:
            q = [r for r in rr if r['candidate'] == name]
            J = np.array([r['metrics']['J'] for r in q], float)
            rm = np.array([r['metrics']['rmse'] for r in q], float)
            ev = np.array([r['metrics']['event_rmse'] for r in q], float)
            fn = np.array([r['metrics']['fn_rmse'] for r in q], float)
            fa = np.array([r['metrics']['faxis_rmse'] for r in q], float)
            dif = []
            for r in q:
                key = (r['cond'], r['tissue'], r['i'])
                if name != 'L' and key in ref:
                    dif.append(math.log1p(min(r['metrics']['J'], 1000.)) - math.log1p(min(ref[key]['J'], 1000.)))
            grouped[regime][name] = dict(n=len(q), L=med(J),
                                         rmse_median=med(rm),
                                         rmse_p90=med(rm, 90),
                                         event_rmse_median=med(ev),
                                         fail=int(sum(bool(r['metrics']['failed']) for r in q)),
                                         long_loss=int(sum(bool(r['metrics']['long_loss']) for r in q)),
                                         fn_rmse_median=med(fn),
                                         faxis_rmse_median=med(fa),
                                         call_us_median=med([r['metrics'].get('call_p50', np.nan) for r in q]),
                                         vs_L=[float(np.mean(dif)) if dif else np.nan, *_boot(dif)])
    return grouped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', nargs='?', default='matrix', choices=('matrix',))
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--workers', type=int, default=15)
    args = ap.parse_args()
    cells = [('rest', 'soft'), ('breath', 'mid')] if args.smoke else SC2.CELLS2
    reps = 2 if args.smoke else 3
    jobs = [(regime, cell, i, BASE) for regime in REGIMES for cell in cells for i in range(reps)]
    print('contact_select %s: %d scenes x %d candidates, workers=%d' % ('smoke' if args.smoke else 'matrix', len(jobs), len(C.NAMES), args.workers), flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = [r for block in ex.map(_job, jobs, chunksize=1) for r in block]
    path = HERE/'results'/'contact_select.json'
    rows += [r for r in json.loads(path.read_text(encoding='utf8'))['rows'] if r['candidate'] not in ONLY]
    summary = _summary(rows)
    for regime in REGIMES:
        print('[%s]' % regime)
        for name in C.NAMES:
            x = summary[regime][name]
            print('  %-3s L %.4f RMSE %.4f/P90 %.4f event %.4f fail %d long %d Fn %.4f Faxis %.4f | vs L %.4f [%.4f, %.4f] call %.0fus' %
                  (name, x['L'], x['rmse_median'], x['rmse_p90'], x['event_rmse_median'], x['fail'], x['long_loss'], x['fn_rmse_median'], x['faxis_rmse_median'], x['vs_L'][0], x['vs_L'][1], x['vs_L'][2], x['call_us_median']))
    out = dict(mode='smoke' if args.smoke else 'matrix', base=BASE, cells=cells, repetitions=reps, candidates=C.NAMES,
               rows=rows, summary=summary)
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results' / 'contact_select.json').write_text(json.dumps(out, indent=1, allow_nan=True), encoding='utf8')


if __name__ == '__main__':
    main()
