"""Open-loop causal replay: measurement gain of the final method versus the path-indexed modulus map (vp_map.py).

Prototype check, reconstructed 2026-11-02 (see vp_map.py). Closed-loop logs of the final method (diagnosis cohort, S2
and T3m, N per family; the original check used N=2). On each log, every UPD ticks the map estimator gets the same
16-tick window regression as the projection (common._build, encoder endpoints) and predicts the stiffness at the
current probe position and depth. Compared with the final method's measurement gain k_good (recorded per tick)
against the true tangent stiffness (log column 8): |ln(k_est/k_true)| in transition ticks (true |d ln k| over the
last 128 ms > 0.1) and elsewhere; share of ticks outside the tolerated band ln 1.2, and the estimator's own posterior
sd of ln k (calibration check). Run: python -X utf8 p3_replay.py [N]
"""
import math
import os
import sys
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
from vp_map import MapEstimator

UPD = 4


class Rec(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform); self.kg = []

    def step(self, o):
        out = super().step(o); self.kg.append((o.n, self.k_good)); return out


def replay(log, pol):
    n_all = log.shape[0]; t, x, F, ktrue = log[:, 0], log[:, 1], log[:, 5], np.maximum(log[:, 8], 1.)
    est = MapEstimator(x.min()-.05, x.max()+.05, upd=UPD)
    kg = dict(pol.kg); out = []
    for n in range(40, n_all, UPD):
        if t[n-1] < .5 or F[n-1] < 1.:
            continue
        j0 = max(pol.ndel+1, n-1-16); m = n-2-j0+1
        X, Y, Z, keep = common._build(log, j0, m, 3, pol.ns, pol.eb, pol.afc, pol.ndel, t[n-1], log[n-1, 2], 16*pol.ns*1e-3, True)
        rows = np.arange(j0, j0+m)[keep]
        if keep.sum() < 8:
            continue
        fit = common.window_fit(log, n, 16, 3, .04, pol.ns, pol.eb, pol.afc, pol.ndel, True)
        res = est.update(log, n, X[keep], Y[keep], Z[keep], rows, k_init=fit['k'] if fit['info'] > 0 else None)
        if res is None or t[n-1] < 1.:
            continue
        k_map = est.k_at(res['x'], res['d']); kt = ktrue[n-1]
        trans = abs(math.log(ktrue[n-1]/ktrue[max(n-17, 0)])) > .1
        sd = est.sd_lnk_at(res['x'], res['d']) if (n//UPD) % 10 == 0 else np.nan
        kgood = kg.get(n, np.nan)
        out.append((trans, math.log(k_map/kt), math.log(kgood/kt) if np.isfinite(kgood) else np.nan, sd))
    return np.array(out, float)


def _ep(a):
    fam, sid, p, T = a
    pol = variants.build('Proposed-IMM', dict(CE.FINAL)); pol.__class__ = Rec
    log, _ = common.run(pol, p, T, sid)
    return fam, replay(log, pol)


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    sc = [s for s in CE.cohort('diag', N) if s[0] in ('S2', 'T3m')]
    with Pool(16) as po:
        rows = po.map(_ep, sc, chunksize=1)
    for fam in ('S2', 'T3m'):
        A = np.concatenate([r[1] for r in rows if r[0] == fam and len(r[1])]); tr = A[:, 0] > .5
        for lab, m in (('transition', tr), ('elsewhere', ~tr)):
            em, eg = np.abs(A[m, 1]), np.abs(A[m, 2]); eg = eg[np.isfinite(eg)]
            print('%-4s %-10s n=%5d | map |ln err| median %.3f p90 %.3f out-of-band %.1f%% | final k_good median %.3f p90 %.3f out-of-band %.1f%%'
                  % (fam, lab, m.sum(), np.median(em), np.percentile(em, 90), 100*np.mean(em > math.log(1.2)),
                     np.median(eg), np.percentile(eg, 90), 100*np.mean(eg > math.log(1.2))), flush=True)
        sd = A[np.isfinite(A[:, 3]), 3]
        print('     map posterior sd(ln k) median %.3f; median |ln err| at those ticks %.3f'
              % (np.median(sd), np.median(np.abs(A[np.isfinite(A[:, 3]), 1]))), flush=True)
