"""Signed error of the measurement gain k_good on hard tissue, band-pass against projection (same rows as main_rows.json)."""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
from kgood_intervention import BASE, HARD

HERE = Path(__file__).resolve().parent


def job(a):
    name, c, t, i, controls, est = a
    sc = S.cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    pol = S.policy(name, controls, est); wrap = S.Measured(pol); good = []

    class Trace:
        def reset(self, platform):
            wrap.reset(platform)

        def __call__(self, o):
            cmd = wrap(o); good.append(pol.k_good); return cmd

    log = S.PP.episode(sc, Trace(), seed)
    if log is None:
        return dict(name=name, cond=c, i=i, failed=True, err=[], write=[])
    g = np.array(good, float)[:len(log)]; kt = log[:len(g), 8]; t_ = log[:len(g), 0]
    # score only before the force first leaves the safe range, so a runaway's own stiffness rise is not counted
    over = np.flatnonzero(log[:len(g), 3] > 20.)
    stop = over[0] if len(over) else len(g)
    use = np.isfinite(g) & (kt > 0) & (t_ >= 1.) & (np.arange(len(g)) < stop)
    err = np.log(g[use]/kt[use])
    wr = np.flatnonzero(np.isfinite(g) & (kt > 0) & (np.r_[True, g[1:] != g[:-1]]) & (np.arange(len(g)) < stop))
    return dict(name=name, cond=c, i=i, failed=bool(log[:, 3].max() > 50), err=err.tolist(), write=np.log(g[wr]/kt[wr]).tolist())


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']
    args = [(n, c, ti, i, controls, est[n]['best']['params']) for n in ('Projection', 'BandPass') for c, ti in HARD for i in range(15)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, args, chunksize=1))
    out = {}
    for n in ('Projection', 'BandPass'):
        rr = [r for r in rows if r['name'] == n]
        e = np.concatenate([r['err'] for r in rr if r['err']]); w = np.concatenate([r['write'] for r in rr if r['write']])
        out[n] = dict(time_p10=float(np.percentile(e, 10)), time_p50=float(np.percentile(e, 50)), time_p90=float(np.percentile(e, 90)),
                      time_below_075=float(np.mean(e < np.log(.75))), time_above_12=float(np.mean(e > np.log(1.2))),
                      write_p10=float(np.percentile(w, 10)), write_p50=float(np.percentile(w, 50)), write_p90=float(np.percentile(w, 90)),
                      write_below_075=float(np.mean(w < np.log(.75))), writes=len(w))
        print(n, {k: round(v, 3) for k, v in out[n].items()})
    (HERE/'results'/'kgood_error.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
