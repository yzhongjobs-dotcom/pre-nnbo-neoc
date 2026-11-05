"""Module-level reading of the two d treatments inside their own closed loops: every estimation window against the true probe-frequency gain.

Per window: estimate, whether the Method accepted it, and the truth at its end. Split by what the scene is doing: quiet, a motion event, the gain changing
(sliding over ribs or inclusions). Truth is read only here, for scoring.
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
LN12 = np.log(1.2)


def _job(a):
    host, cfg, cond, tissue, dyn, ph, i = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=ph))
    pol = D.make(host, cfg, sc)
    log = P3.episode(sc, pol, 900+i)
    if log is None:
        return None
    t = log[:, 0]; k = log[:, 21]
    ev = np.zeros(len(t), bool)
    for e in sc['event_info']:
        ev |= (t >= e['start']) & (t <= e['end']+.3)
    dk = np.abs(np.gradient(np.log(np.maximum(k, 1.)), t))
    changing = np.convolve(dk > 2., np.ones(30), 'same') > 0       # gain changing faster than 2/s, within about a quarter second
    rows = []
    for tk, kh, sig, chi, ok in pol.accept_trace:
        j = int(np.clip(np.searchsorted(t, tk), 0, len(t)-1))
        if t[j] < 1. or k[j] <= 200. or not np.isfinite(kh) or kh <= 0.:
            continue
        kind = 'event' if ev[j] else ('changing' if changing[j] else 'quiet')
        rows.append((kind, float(np.log(kh/k[j])), bool(ok), float(sig/kh) if np.isfinite(sig) else np.inf, float(chi)))
    return dict(host=host, cond=cond, tissue=tissue, dynamics=dyn, phase=ph, i=i, windows=rows)


def main():
    chosen = pickle.loads(Path(sys.argv[1]).read_bytes())['chosen']
    jobs = [(h, chosen[h], c, t, dyn, ph, i) for h in ('D4', 'P') for c, t in D.TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        out = [r for r in ex.map(_job, jobs, chunksize=1) if r]
    print('per host and kind of moment: windows, accepted share, |ln k/k| median and P90 of all windows, of accepted windows, accepted share beyond ln 1.2, median rel sigma, median chi')
    summary = {}
    for h in ('D4', 'P'):
        for kind in ('quiet', 'event', 'changing'):
            w = [x for r in out if r['host'] == h for x in r['windows'] if x[0] == kind]
            if not w:
                continue
            e = np.abs([x[1] for x in w]); ok = np.array([x[2] for x in w]); ea = e[ok]
            s = dict(n=len(w), accepted=float(ok.mean()), all_med=float(np.median(e)), all_p90=float(np.percentile(e, 90)),
                     acc_med=float(np.median(ea)) if len(ea) else np.nan, acc_p90=float(np.percentile(ea, 90)) if len(ea) else np.nan,
                     acc_bad=float(np.mean(ea > LN12)) if len(ea) else np.nan, rel=float(np.median([x[3] for x in w])), chi=float(np.median([x[4] for x in w])), signed_acc=float(np.median([x[1] for x in w if x[2]])) if len(ea) else np.nan)
            summary[f'{h} {kind}'] = s
            print(f"  {h:3s} {kind:9s} n {s['n']:6d} acc {s['accepted']:.2f} | all {s['all_med']:.3f}/{s['all_p90']:.3f} | accepted {s['acc_med']:.3f}/{s['acc_p90']:.3f} beyond 1.2 {s['acc_bad']:.3f} | rel {s['rel']:.3f} chi {s['chi']:.2f} | signed median of accepted {s['signed_acc']:+.3f}")
    print('by condition, accepted windows: median |ln k/k| and share beyond ln 1.2, D4 then P')
    for c, t in D.TEST:
        line = []
        for h in ('D4', 'P'):
            w = [x for r in out if r['host'] == h and r['cond'] == c for x in r['windows'] if x[2]]
            e = np.abs([x[1] for x in w])
            line.append(f"{np.median(e):.3f} {np.mean(e > LN12):.3f} acc {len(w)/max(sum(len(r['windows']) for r in out if r['host'] == h and r['cond'] == c), 1):.2f}")
        print(f"  {c:10s} D4 {line[0]} | P {line[1]}")
    (HERE/'results'/'d_module.json').write_text(json.dumps(summary, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
