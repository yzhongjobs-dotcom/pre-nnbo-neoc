"""Final confirmation on unused seeds of the configuration chosen after the ablation.

Simple    relaxed prior (D4, 24 ticks, 14 Hz) with weak power-law memory in the regressor and k_good centred in the stable band.
          The force-to-motion conversion and the first-gain rule stay as in the existing method.
Full      the same with memory in the conversion and the first-gain rule, for comparison.
Baseline  P5_16 tuned, the existing form with the same control parameters. All on six plant laws.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP
from contact_compare import ModelPolicy, RIDGE, CTL, P5, L, key, boot, scene, PLANTS, GROUPS
from dynamic_gain_value import Truth, stats
from ridge_estimator import RidgeProjection, POP_ALPHA, POP_TAU

HERE = Path(__file__).resolve().parent
BASE = 415000
METHODS = ('P5_16', 'Simple', 'Full')


def job(a):
    name, plant, c, t, i = a
    sc = scene(c, t, i, BASE, plant); seed = BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    if name == 'P5_16':
        pol = ModelPolicy(S.estimator('Projection', P5), single_model=False, **p); pol.centre = False; pol.acquire = False
    else:
        pol = ModelPolicy(RidgeProjection(memory=True, **RIDGE), single_model=False, **p)
        if name == 'Full':
            pol.zalpha = POP_ALPHA; pol.ztau = POP_TAU
        else:
            pol.acquire = False
    log = PP.episode(sc, Truth(pol), seed)
    return dict(name=name, plant=plant, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, 20.8))


def main():
    args = [(n, 'footprint', c, t, i) for n in METHODS for c, t in S.CELLS for i in range(10)]
    args += [(n, pl, c, t, i) for n in METHODS for pl in PLANTS if pl != 'footprint' for c, t in S.CELLS for i in range(5)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(job, args, chunksize=1))
    out = {}
    for pl in PLANTS:
        ref = {key(r): r['metrics'] for r in rows if r['plant'] == pl and r['name'] == 'P5_16'}
        out[pl] = {}
        print(f'[{pl}]')
        for n in METHODS:
            rr = [r for r in rows if r['plant'] == pl and r['name'] == n]; s = stats(rr)
            d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(d)
            g = {}
            for gn, f in GROUPS:
                q = [r for r in rr if f(r)]; dg = np.array([L(r['metrics'])-L(ref[key(r)]) for r in q]); cg = boot(dg)
                g[gn] = [float(dg.mean()), float(cg[0]), float(cg[1])]
            out[pl][n] = dict(s, vs_p5=[float(d.mean()), float(c[0]), float(c[1])], groups=g)
            print(f"  {n:7s} n{s['n']} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.2f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs P5 {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]")
            if pl == 'footprint':
                print('      '+' '.join(f"{gn}:{v[0]:+.2f}[{v[1]:+.2f},{v[2]:+.2f}]" for gn, v in g.items()))
    calls = [r['metrics']['call_p50'] for r in rows if r['name'] == 'Simple']; p5c = [r['metrics']['call_p50'] for r in rows if r['name'] == 'P5_16']
    print('call P50 us: Simple', round(float(np.nanmedian(calls))), 'P5_16', round(float(np.nanmedian(p5c))))
    print('mean L over plants', {n: round(float(np.mean([out[pl][n]['L'] for pl in PLANTS])), 3) for n in METHODS})
    (HERE/'results'/'final_config.json').write_text(json.dumps(dict(summary=out, rows=rows), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
