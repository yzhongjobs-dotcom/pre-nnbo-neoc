"""Where is the remaining error? The D⁴ host with the feed-forward given the true tissue displacement (body position plus surface height under the probe),
extrapolated with its true velocity over the same horizon. Everything else unchanged. Diagnostic only."""
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


def _job(a):
    cfg, cond, tissue, dyn, phase, i, mode = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=phase))
    pol = D.make('D4', cfg, sc)
    grid = sc['field_x0']+np.arange(len(sc['h_s']))*sc['dx']
    if mode == 'true_d':
        original = pol.drift_pred; memory = {}
        def drift_pred(o, n, original=original):
            out = original(o, n)
            if out is None or not hasattr(o, 'truth_body'):
                return out
            s_now = o.truth_body+float(np.interp(o.log[n-1, 1], grid, sc['h_s']))
            prev = memory.get('s'); memory['s'] = s_now
            v = 0. if prev is None else (s_now-prev)/pol.Tc
            h = (pol.ndel+.5)*pol.Tc
            return (None, s_now+v*(h+pol.tau))
        pol.drift_pred = drift_pred
    log = P3.episode(sc, pol, 900+i, diagnostic=mode != 'own')
    m = M.metric(log, pol, sc); m.update(cond=cond, tissue=tissue, dynamics=dyn, phase=phase, i=i, mode=mode)
    return m


def main():
    chosen = pickle.loads(Path(sys.argv[1]).read_bytes())['chosen']['D4']
    jobs = [(chosen, c, t, dyn, ph, i, mode) for mode in ('own', 'true_d') for c, t in D.TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    base = {key(r): r for r in rows if r['mode'] == 'own'}
    for mode in ('own', 'true_d'):
        q = [r for r in rows if r['mode'] == mode]
        d = np.array([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in q])
        print(f"  {mode:6s} {D.line(q)}  vs own {d.mean():+.3f}")
        for c, t in D.TEST:
            qq = [r for r in q if r['cond'] == c]
            print(f"       {c:10s} rmse {np.median([r['rmse'] for r in qq]):.3f} event {np.nanmedian([r['event_rmse'] for r in qq]) if any(np.isfinite(r['event_rmse']) for r in qq) else float('nan'):.3f} peak max {max(r['peak'] for r in qq):.1f} L diff {np.mean([np.log1p(r['J'])-np.log1p(base[key(r)]['J']) for r in qq]):+.3f}")
    (HERE/'results'/'d_motion_ceiling.json').write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
