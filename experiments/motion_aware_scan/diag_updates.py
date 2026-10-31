"""Which accepted windows give a wrong measurement gain? (diagnosis)

Records every update of k_good (precise and consistent window) with the true stiffness and window features, on a cohort,
and reports the error distribution ln(k_hat/k_true) conditioned on: relative precision, normalised residual, force level,
force rate over the window, time since the last push/retract ramp ended, and whether the window overlaps a ramp.
Run: COHORT=diag N=50 FAMS=T2m-stiff python -X utf8 diag_updates.py
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import cohort_eval as CE
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']


class Rec(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.upd = []; self._prev = np.nan

    def step(self, o):
        y = super().step(o)
        if self.k_good != self._prev and np.isfinite(self.k_good):
            n = o.n; f = self.fit or {}
            w = o.log[max(0, n-1-self.params['nw']):n, 5]
            self.upd.append((n, self.k_good, float(o.log[n-1, 8]), self.rel, f.get('chi2', np.nan),
                             float(o.log[n-1, 5]), float(np.max(w)-np.min(w)) if len(w) else np.nan))
            self._prev = self.k_good
        return y


def ramp_info(t, p):
    """(inside ramp, time since last ramp end)"""
    ends = []; inside = False
    for k in range(int(p[15])):
        t0, tr, hold = p[16+4*k], p[18+4*k], p[19+4*k]
        for a in (t0, t0+tr+hold):
            if a-0.13 <= t <= a+tr: inside = True
            if a+tr <= t: ends.append(t-(a+tr))
    return inside, (min(ends) if ends else np.inf)


def _ep(sc):
    fam, sid, p, T = sc
    pol = Rec(single_model=False, **BASE)
    log, pol = common.run(pol, p, T, sid)
    out = []
    for (n, kh, kt, rel, chi, F, dF) in pol.upd:
        t = log[n-1, 0]
        if t < 1 or kt <= 0: continue
        ins, since = ramp_info(t, p)
        out.append((float(math.log(kh/kt)), rel, chi, F, dF, ins, since))
    return out


if __name__ == '__main__':
    fams = os.environ.get('FAMS', 'T2m-stiff').split(',')
    sc = [s for s in CE.cohort(os.environ.get('COHORT', 'diag'), int(os.environ.get('N', '50'))) if s[0] in fams]
    with Pool(16) as po:
        rows = [r for res in po.map(_ep, sc) for r in res]
    A = np.array(rows, float)
    json.dump(dict(note='updates of k_good by the earlier version (no acceptance rule, model-reconstructed displacement) on the diag cohort; '
                   'err = ln(k_hat/k_true); too-low beyond band: err < -0.18 (= -ln 1.2)', columns=['ln_err', 'rel', 'chi2', 'F', 'window_dF', 'inside_ramp', 'since_ramp_end'],
                   rows=[[float(x) for x in r] for r in rows]), open(os.path.join(H, 'results_v2', 'diag_updates.json'), 'w'))
    e = A[:, 0]; big = e < -0.18                   # gain too low by more than 1/1.2 (dangerous side)
    print('updates %d; |err| q50 %.3f q95 %.3f; fraction too low beyond band %.3f, too high beyond band %.3f'
          % (len(e), np.median(abs(e)), np.percentile(abs(e), 95), big.mean(), (e > 0.29).mean()))
    def cond(name, m):
        if m.sum() == 0: return
        print('  %-34s n=%5d  too-low rate %.3f  median err %+.3f  q05 %+.3f' % (name, m.sum(), big[m].mean(), np.median(e[m]), np.percentile(e[m], 5)))
    rel, chi, F, dF, ins, since = A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5] > 0, A[:, 6]
    for lo, hi in ((0, .05), (.05, .1), (.1, .2), (.2, .3)): cond('rel in [%.2f,%.2f)' % (lo, hi), (rel >= lo) & (rel < hi))
    for lo, hi in ((0, 1), (1, 2), (2, 3)): cond('chi2 in [%d,%d)' % (lo, hi), (chi >= lo) & (chi < hi))
    for lo, hi in ((0, 3), (3, 5), (5, 8), (8, 1e9)): cond('force in [%g,%g)' % (lo, hi), (F >= lo) & (F < hi))
    for lo, hi in ((0, .3), (.3, .6), (.6, 1.2), (1.2, 1e9)): cond('window force range [%g,%g)' % (lo, hi), (dF >= lo) & (dF < hi))
    cond('inside ramp', ins); cond('outside ramp', ~ins)
    for lo, hi in ((0, .2), (.2, .5), (.5, 1.), (1., 1e9)): cond('since ramp end [%g,%g) s' % (lo, hi), (~ins) & (since >= lo) & (since < hi))
    cond('never a ramp before', (~ins) & ~np.isfinite(since))
