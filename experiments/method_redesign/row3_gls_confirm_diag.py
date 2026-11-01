"""E4-A (GLS value in the measurement channel, projection acceptance and feedback unchanged) vs the final method:
diagnosis cohort 50 per family + push-ramp and ridge-ratio sweeps (20 per dose), paired bootstrap of mean RMSE.
Run: python -X utf8 row3_gls_confirm_diag.py
"""
import os, sys, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, HERE)
import common, cohort_eval as CE, sweep as SW
from row3_gls_meas_only import build

def _ep(a):
    cfg, (tag, fam, sid, p, T) = a
    try:
        log, _ = common.run(build(cfg), p, T, sid)
        m = common.task_metrics(log, 20.8)
        return dict(cfg=cfg, tag=tag, family=fam, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'], loss_max=m['loss_max'])
    except Exception as ex:
        return dict(cfg=cfg, tag=tag, family=fam, sid=sid, rmse=np.nan, peak=np.inf, ok=False, err=str(ex))

def boot(d, B=4000, seed=1):
    r = np.random.default_rng(seed); d = np.asarray(d)
    m = [r.choice(d, len(d)).mean() for _ in range(B)]
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

if __name__ == '__main__':
    sc = [('diag', f, sid, p, T) for (f, sid, p, T) in CE.cohort('diag', 50)]
    sc += [(kind, '%s:%g' % (kind, dose), sid, p, T) for (kind, dose, sid, p, T) in SW.scenes()]
    CFG = ('Proposed', 'E4-A')
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    json.dump(rows, open(os.path.join(HERE, 'results', 'row3_gls_confirm_diag.json'), 'w'), indent=1)
    groups = list(dict.fromkeys(r['family'] for r in rows))
    for g in groups:
        R = {c: {r['sid']: r for r in rows if r['cfg'] == c and r['family'] == g} for c in CFG}
        sids = sorted(set(R['Proposed']) & set(R['E4-A']))
        a = np.array([R['Proposed'][s]['rmse'] if R['Proposed'][s]['ok'] else 1e3 for s in sids]); b = np.array([R['E4-A'][s]['rmse'] if R['E4-A'][s]['ok'] else 1e3 for s in sids])
        bad = [sum((not R[c][s]['ok']) or R[c][s]['peak'] > 50 for s in sids) for c in CFG]
        lost = [sum(R[c][s].get('loss_max', 0) > .5 for s in sids) for c in CFG]
        d = boot(b-a)
        print('%-12s n=%2d  Proposed %.3f/%.3f  E4-A %.3f/%.3f  runaway %d/%d  longloss %d/%d  diff %+.3f [%+.3f,%+.3f]' % (
            g, len(sids), np.median(a), a.mean(), np.median(b), b.mean(), bad[0], bad[1], lost[0], lost[1], *d), flush=True)
