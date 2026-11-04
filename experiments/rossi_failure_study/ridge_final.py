"""Final acceptance on seeds not used for any choice: tuned projection against the full candidate, with and without centring k_good.

Centring uses the stable band of the repository's filter rule, k_good = k/sqrt(.75*1.2), an offset of +.053 in ln k, not a tuned value.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP
from deployable_contact import Deployed
from dynamic_gain_value import stats
from failure_replay import Measured, cell_scene
from ridge_estimator import RidgeProjection

HERE = Path(__file__).resolve().parent
BASE = 401000
N_EVAL = 10
CENTRE = 1./math.sqrt(.75*1.2)


class Candidate(Deployed):
    centre = False

    def _stiffness_est(self, o, n):
        before = self.k_good
        out = super()._stiffness_est(o, n)
        if self.centre and np.isfinite(self.k_good) and self.k_good != before:
            self.k_good *= CENTRE
        return out


def one(a):
    spec, controls, ci, c, t, i = a
    sc = cell_scene(c, t, i, 340000); seed = 340000+ci*100+i
    p = dict(S.FINAL_C); p.update(controls)
    pol = Candidate(RidgeProjection(**spec), single_model=False, **p); pol.sc = sc; pol.centre = True
    return math.log1p(min(S._metric(PP.episode(sc, Measured(pol), seed), pol, sc, 20.8)['J'], 1000.))


def job(a):
    name, spec, c, t, i, controls = a
    sc = cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(controls)
    if name == 'P5_16 tuned':
        pol = S.SwapPolicy(S.estimator('Projection', spec), single_model=False, **p)
    else:
        pol = Candidate(RidgeProjection(**spec), single_model=False, **p); pol.sc = sc; pol.centre = name.endswith('centred')
    return dict(name=name, cond=c, tissue=t, i=i, metrics=S._metric(PP.episode(sc, Measured(pol), seed), pol, sc, 20.8))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']['Projection']['best']['params']
    cands = [dict(nw=24, order=4, fc=14., memory=True, sigma_scale=a, chi_scale=b) for a in (.8, 1.5) for b in (.5, 1., 2.)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        scored = []
        for spec in cands:
            args = [(spec, controls, ci, c, ti, i) for ci, (c, ti) in enumerate(S.CELLS) for i in range(3)]
            scored.append((float(np.mean(list(ex.map(one, args, chunksize=1)))), spec))
        v, best = min(scored, key=lambda x: x[0])
        print('gate selection', [round(s, 3) for s, _ in scored], best, flush=True)
        methods = {'P5_16 tuned': est, 'Candidate': best, 'Candidate centred': best}
        args = [(n, sp, c, ti, i, controls) for n, sp in methods.items() for c, ti in S.CELLS for i in range(N_EVAL)]
        rows = list(ex.map(job, args, chunksize=1))
    L = lambda m: math.log1p(min(m['J'], 1000)); key = lambda r: (r['cond'], r['tissue'], r['i'])
    ref = {key(r): r['metrics'] for r in rows if r['name'] == 'P5_16 tuned'}
    groups = (('soft', lambda r: r['tissue'] == 'soft'), ('mid', lambda r: r['tissue'] == 'mid'), ('hard', lambda r: r['tissue'] == 'hard'),
              ('event', lambda r: r['cond'] in ('flinch', 'cough')), ('fast', lambda r: r['cond'] == 'shift_fast'), ('smooth', lambda r: r['cond'] in ('breath', 'shift_slow')),
              ('tremor', lambda r: r['cond'] == 'tremor'), ('rib+combo', lambda r: r['cond'] in ('rib', 'combo')))
    boot = lambda d: np.percentile(d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1), [2.5, 97.5])
    summary = {}
    for n in methods:
        rr = [r for r in rows if r['name'] == n]; s = stats(rr)
        d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(d)
        dm = np.array([r['metrics']['rmse']-ref[key(r)]['rmse'] for r in rr]); dm = dm[np.isfinite(dm)]; cm = boot(dm)
        calls = [r['metrics']['call_p50'] for r in rr]; calls99 = [r['metrics']['call_p99'] for r in rr]
        summary[n] = dict(s, vs_p5=[float(d.mean()), float(c[0]), float(c[1])], rmse_pair_median=float(np.median(dm)), win=float(np.mean(d < 0)), call_p50=float(np.nanmedian(calls)), call_p99=float(np.nanmedian(calls99)),
                          groups={g: stats([r for r in rr if f(r)]) for g, f in groups})
        print(f"{n:18s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.3f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs P5 {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] pair rmse {np.median(dm):+.4f} win {np.mean(d < 0):.2f} call {np.nanmedian(calls):.0f}/{np.nanmedian(calls99):.0f} us")
        print('   ', ' '.join(f"{g}:{summary[n]['groups'][g]['L']:.2f}/{summary[n]['groups'][g]['rmse_med']:.3f}/{summary[n]['groups'][g]['event_med'] or 0:.3f}/{summary[n]['groups'][g]['fail']}" for g, _ in groups))
    (HERE/'results'/'ridge_final.json').write_text(json.dumps(dict(best=best, gate=[s for s, _ in scored], summary=summary, rows=rows), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
