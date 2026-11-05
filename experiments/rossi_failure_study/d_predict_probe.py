"""Does the motion prediction help or hurt in each kind of event? Same D⁴ host; only the use of the tracked d in the feed-forward changes."""
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
CELLS = [('shift_fast', 'relaxed'), ('flinch', 'active'), ('cough', 'relaxed'), ('combo', 'mixed')]


def _job(a):
    cfg, cond, tissue, dyn, phase, i, mode = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=phase))
    pol = D.make('D4', cfg, sc)
    if mode == 'off':
        pol.params['feedforward'] = False
    elif mode in ('filtered', 'raw'):
        original = pol.drift_pred
        def drift_pred(o, n, original=original):
            out = original(o, n)
            if out is None:
                return out
            return (None, float(pol.tracker.x[0]) if mode == 'filtered' else float(pol.Z))
        pol.drift_pred = drift_pred
    log = P3.episode(sc, pol, 900+i)
    m = M.metric(log, pol, sc); m.update(cond=cond, tissue=tissue, dynamics=dyn, phase=phase, i=i, mode=mode)
    return m


def main():
    chosen = pickle.loads(Path(sys.argv[1]).read_bytes())['chosen']['D4']
    modes = ('predicted', 'filtered', 'raw', 'off')
    jobs = [(chosen, c, t, dyn, ph, i, mode) for mode in modes for c, t in CELLS for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    for c, _ in CELLS:
        print(c)
        for mode in modes:
            q = [r for r in rows if r['cond'] == c and r['mode'] == mode]
            print(f"   {mode:9s} rmse {np.median([r['rmse'] for r in q]):.3f}/{np.mean([r['rmse'] for r in q]):.3f} event {np.nanmedian([r['event_rmse'] for r in q]):.3f} peak max {max(r['peak'] for r in q):.1f} L {np.mean([np.log1p(r['J']) for r in q]):.3f}")
    (HERE/'results'/'d_predict_probe.json').write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
