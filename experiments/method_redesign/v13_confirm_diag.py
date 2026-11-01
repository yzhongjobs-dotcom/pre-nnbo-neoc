"""V1+V3 candidates on the diagnosis cohort (50 per family) and the dose sweeps (20 per dose), paired with the final method.
  Proposed, V1 (servo-lag/delay identification), V1+E2-A, V1+E4-A. Run: python -X utf8 v13_confirm_diag.py
"""
import os, sys, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan')); sys.path.insert(0, HERE)
import common, variants, cohort_eval as CE, sweep as SW
from v1_servo_id import ServoID
from v3_temporal_estimators import MeasChannel

class SIDMeas(ServoID, MeasChannel):
    pass

def build(cfg):
    if cfg == 'Proposed': return variants.build('Proposed-IMM', dict(CE.FINAL))
    if cfg == 'V1': return ServoID(single_model=False, **CE.FINAL)
    return SIDMeas(cfg.split('+')[1], single_model=False, **CE.FINAL)

def _ep(a):
    cfg, (group, sid, p, T) = a
    try:
        log, _ = common.run(build(cfg), p, T, sid)
        m = common.task_metrics(log, 20.8)
        return dict(cfg=cfg, group=group, sid=sid, rmse=m['rmse'], peak=m['peak'], ok=not m['diverged'], loss_max=m['loss_max'])
    except Exception as ex:
        return dict(cfg=cfg, group=group, sid=sid, rmse=np.nan, peak=np.inf, ok=False, loss_max=0., err=repr(ex))

def boot(d, B=4000, seed=1):
    r = np.random.default_rng(seed); d = np.asarray(d)
    m = [r.choice(d, len(d)).mean() for _ in range(B)]
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

if __name__ == '__main__':
    sc = [(f, sid, p, T) for (f, sid, p, T) in CE.cohort('diag', 50)]
    sc += [('%s:%g' % (k, d), sid, p, T) for (k, d, sid, p, T) in SW.scenes()]
    CFG = ('Proposed', 'V1', 'V1+E2-A', 'V1+E4-A')
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    json.dump(rows, open(os.path.join(HERE, 'results', 'v13_confirm_diag.json'), 'w'), indent=1)
    for g in list(dict.fromkeys(r['group'] for r in rows)):
        R = {c: {r['sid']: r for r in rows if r['cfg'] == c and r['group'] == g} for c in CFG}
        sids = sorted(R['Proposed'])
        line = '%-12s' % g
        for c in CFG:
            v = np.array([R[c][s]['rmse'] if R[c][s]['ok'] else 1e3 for s in sids])
            bad = sum((not R[c][s]['ok']) or R[c][s]['peak'] > 50 for s in sids); ll = sum(R[c][s]['loss_max'] > .5 for s in sids)
            if c == 'Proposed': base = v; line += ' | %s %.3f/%.3f r%d L%d' % (c, np.median(v), v.mean(), bad, ll)
            else:
                d = boot(np.minimum(v, 1e3)-np.minimum(base, 1e3)); line += ' | %s %.3f/%.3f r%d L%d d%+.3f[%+.3f,%+.3f]' % (c, np.median(v), v.mean(), bad, ll, *d)
        print(line, flush=True)
