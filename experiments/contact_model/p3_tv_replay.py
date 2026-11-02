"""Open-loop causal replay of measurement-gain filters (tv_projection.py) against the final method.

Closed-loop logs of the final method (diagnosis cohort, S2 / T3m / stiff push / fast push). At every tick: the final
method's measurement gain k_good (recorded), the constant-window estimate, and Kalman/IMM filters fed by every window,
each against the true tangent stiffness (log column 8). Reported per family: |ln(k/k_true)| median / p90 and the share
outside ln 1.2 (the band the feed-forward loop tolerates), split into transition ticks (true |d ln k| over the window
> 0.1) and the rest. Window-end variance of tv_fit and spline_fit: p3_evidence_checks.py.
Run: python -X utf8 p3_tv_replay.py [N]
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
from tv_projection import tv_fit, GainFilter, GainIMM

FAMS = ('S2', 'T3m', 'T2m-stiff', 'T2f-stiff')
CHI, REL = CE.FINAL['chi2_max'], CE.FINAL['rel_max']


class Rec(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform); self.kg = {}

    def step(self, o):
        out = super().step(o); self.kg[o.n] = self.k_good; return out


FILT = {'KF level': dict(q_level=1e-3), 'KF level+rate': dict(q_level=1e-4, q_rate=1.),
        'IMM .98': ('imm', dict()), 'IMM .9': ('imm', dict(p_stay=.9)), 'IMM q2r=10': ('imm', dict(q2r=10.))}


def replay(log, pol):
    t, F, kt = log[:, 0], log[:, 5], np.maximum(log[:, 8], 1.)
    out = []; filt = {k: (GainIMM(**v[1]) if isinstance(v, tuple) else GainFilter(**v)) for k, v in FILT.items()}
    for n in range(40, log.shape[0]):
        if t[n-1] < 1. or F[n-1] < 1.:
            continue
        j0 = max(pol.ndel+1, n-17); m = n-2-j0+1
        X, Y, Z, keep = common._build(log, j0, m, 3, pol.ns, pol.eb, pol.afc, pol.ndel, t[n-1], log[n-1, 2],
                                      16*pol.ns*1e-3, True)
        if keep.sum() < 10:
            continue
        f = tv_fit(X[keep], Y[keep], Z[keep], .04)
        if not np.isfinite(f['k_e']) or f['k_e'] <= 0 or f['k_const'] <= 0:
            continue
        ktrue = kt[n-1]; trans = abs(math.log(kt[n-1]/kt[max(n-17, 0)])) > .1
        kg = pol.kg.get(n, np.nan)
        acc_c = f['sd_const']/f['k_const'] <= REL and f['chi2_const'] <= CHI
        acc_t = f['sd_e']/f['k_e'] <= REL and f['chi2'] <= CHI
        lks = []
        for nm, fl in filt.items():
            v = fl.step(f['k_const'], f['sd_const']/f['k_const'], f['chi2_const'], F[n-1])
            lks.append(v-math.log(ktrue) if v is not None else np.nan)
        out.append((trans, math.log(kg/ktrue) if np.isfinite(kg) and kg > 0 else np.nan,
                    math.log(f['k_const']/ktrue), *lks))
    return np.array(out, float)


def _ep(a):
    fam, sid, p, T = a
    pol = variants.build('Proposed-IMM', dict(CE.FINAL)); pol.__class__ = Rec
    log, _ = common.run(pol, p, T, sid)
    return fam, replay(log, pol)


def stat(e):
    e = np.abs(e[np.isfinite(e)])
    return '%.3f/%.3f %4.1f%%' % (np.median(e), np.percentile(e, 90), 100*np.mean(e > math.log(1.2)))


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    sc = [s for s in CE.cohort('diag', N) if s[0] in FAMS]
    with Pool(16) as po:
        rows = po.map(_ep, sc, chunksize=1)
    print('|ln err| median/p90, share outside ln1.2')
    for fam in FAMS:
        A = np.concatenate([r[1] for r in rows if r[0] == fam and len(r[1])])
        for lab, msk in (('transition', A[:, 0] > .5), ('elsewhere', A[:, 0] < .5)):
            B = A[msk]
            print('%-9s %-10s n=%6d | k_good %s | window %s | ' % (fam, lab, len(B), stat(B[:, 1]), stat(B[:, 2]))
                  + ' | '.join('%s %s' % (nm, stat(B[:, 3+i])) for i, nm in enumerate(FILT)), flush=True)
