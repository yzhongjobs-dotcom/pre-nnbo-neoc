"""Saved evidence for statements in project_records/handoff/03 that had no raw output.

(a) Relative sd of the stiffness at the window end: constant window, first-order time-varying window (tv_fit) and the
    multi-block spline without arrival cost (spline_fit, prior=None), in closed-loop logs of the final method
    (diagnosis cohort, 10 per family).
(b) Mode probability of the "changing" model of GainIMM in closed loop (GainIMM-0.98+edge, a_f = 0.10; 5 per family).
Run: python -X utf8 p3_evidence_checks.py -> results/p3_evidence_checks.txt
"""
import math
import os
import numpy as np
from multiprocessing import Pool

import p3_gain_closed_loop as G
from p3_gain_closed_loop import CE, common, variants
from tv_projection import tv_fit, spline_fit

HERE = os.path.dirname(os.path.abspath(__file__))
FAMS = ('S2', 'T3m', 'T2m-stiff', 'T2f-stiff')


def _sd(a):
    fam, sid, p, T = a
    pol = variants.build('Proposed-IMM', dict(CE.FINAL))
    log, _ = common.run(pol, p, T, sid); t, F = log[:, 0], log[:, 5]; out = []
    for n in range(60, log.shape[0], 2):
        if t[n-1] < 1. or F[n-1] < 1.:
            continue
        j0 = max(pol.ndel+1, n-17); m = n-2-j0+1
        X, Y, Z, keep = common._build(log, j0, m, 3, pol.ns, pol.eb, pol.afc, pol.ndel, t[n-1], log[n-1, 2], 16*pol.ns*1e-3, True)
        if keep.sum() < 10:
            continue
        f = tv_fit(X[keep], Y[keep], Z[keep], .04); sp = spline_fit(log, n, pol.ns, pol.eb, pol.afc, pol.ndel, nb=3)
        if sp is None or not np.isfinite(f['k_e']) or min(f['k_e'], f['k_const'], sp['k_e']) <= 0:
            continue
        out.append((f['sd_const']/f['k_const'], f['sd_e']/f['k_e'], sp['sd_e']/sp['k_e']))
    return fam, np.array(out)


class RecMu(G.GainIMMMethod):
    def reset(self, pl):
        super().reset(pl); self.mus = []

    def step(self, o):
        r = super().step(o)
        if self.gimm.x is not None and o.log[o.n-1, 0] > 1.:
            self.mus.append(float(self.gimm.mu[1]))
        return r


def _mu(a):
    fam, sid, p, T = a
    pol = G.build('GainIMM-0.98+edge@a0.1'); pol.__class__ = RecMu
    common.run(pol, p, T, sid)
    return fam, np.array(pol.mus)


if __name__ == '__main__':
    with Pool(16) as po:
        sd = po.map(_sd, [s for s in CE.cohort('diag', 10) if s[0] in FAMS], chunksize=1)
        mu = po.map(_mu, CE.cohort('diag', 5), chunksize=1)
    out = ['median relative sd of the window-end stiffness: constant window | first-order TV window | 3-block spline without arrival cost']
    for fam in FAMS:
        A = np.concatenate([r[1] for r in sd if r[0] == fam and len(r[1])])
        out.append('%-9s %.3f | %.3f | %.3f (ticks %d)' % (fam, *np.median(A, 0), len(A)))
    out.append('GainIMM mode probability of the changing model: median, p5, p95')
    for fam in CE.FAMS:
        A = np.concatenate([r[1] for r in mu if r[0] == fam])
        out.append('%-9s %.3f %.3f %.3f' % (fam, np.median(A), np.percentile(A, 5), np.percentile(A, 95)))
    print('\n'.join(out))
    open(os.path.join(HERE, 'results', 'p3_evidence_checks.txt'), 'w', encoding='utf-8').write('\n'.join(out)+'\n')
