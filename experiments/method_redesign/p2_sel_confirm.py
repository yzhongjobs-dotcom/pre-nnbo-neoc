"""50-per-family confirmation (diagnosis cohort, all five families) of supervisory selection among servo-lag lead
filters {raw, 10, 6, 3 Hz}: monitoring signal = exponentially weighted (0.3 s) squared error between each candidate's
feed-forward passed through the identified servo model and the measured motion measurement; scale-independent
hysteresis h. Variants h = 1, 3, 7. Paired bootstrap (4000) of the per-scene RMSE difference to the final method.
Run: python -X utf8 p2_sel_confirm.py -> results/p2_sel_confirm.txt
"""
import os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p2_hf_attribution as H
from p2_hf_attribution import common, CE
CFG = {'Proposed': dict(), 'SelH1': dict(sched='select', servo_crit=True, hyst=1.),
       'SelH3': dict(sched='select', servo_crit=True, hyst=3.), 'SelH7': dict(sched='select', servo_crit=True, hyst=7.)}
if len(sys.argv) > 1 and sys.argv[1] == 'asym':
    CFG = {'Proposed': dict(), 'SelA1': dict(sched='select', servo_crit=True, hyst=1., asym=True),
           'SelA2': dict(sched='select', servo_crit=True, hyst=2., asym=True), 'SelA3': dict(sched='select', servo_crit=True, hyst=3., asym=True)}
def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(H.mk(**CFG[cfg]), p, T, sid); m = common.task_metrics(log, 20.8)
    return cfg, fam, sid, (m['rmse'] if not m['diverged'] else np.nan), m['peak'], m['loss_max']
if __name__ == '__main__':
    sc = CE.cohort('diag', 50)
    with Pool(16) as po: rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    rng = np.random.default_rng(0); out = []
    for fam in CE.FAMS:
        base = {r[2]: r for r in rows if r[0] == 'Proposed' and r[1] == fam}
        for c in CFG:
            R = [r for r in rows if r[0] == c and r[1] == fam]; v = np.array([r[3] for r in R])
            line = '%-9s %-8s %.3f/%.3f r%d L%d' % (fam, c, np.nanmedian(v), np.nanmean(v), sum((not np.isfinite(r[3])) or r[4] > 50 for r in R), sum(r[5] > .5 for r in R))
            if c != 'Proposed':
                d = np.array([r[3]-base[r[2]][3] for r in R]); d = d[np.isfinite(d)]
                bs = [np.mean(rng.choice(d, len(d))) for _ in range(4000)]
                line += ' | mean diff %+.3f [%+.3f, %+.3f]' % (d.mean(), *np.percentile(bs, [2.5, 97.5]))
            out.append(line); print(line, flush=True)
    open(os.path.join(HERE, 'results', 'p2_sel_confirm%s.txt' % ('_'+sys.argv[1] if len(sys.argv) > 1 else '')), 'w', encoding='utf-8').write('\n'.join(out)+'\n')
