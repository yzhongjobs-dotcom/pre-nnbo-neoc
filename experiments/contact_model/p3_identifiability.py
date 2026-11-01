"""Identifiability of the local modulus field from window stiffness, normal force and moment (methods 21.5(a)).

Trajectories: the final method in closed loop on diagnosis scenes (S2 static span, T3m spine scan with breathing and
pushes; 10 each). Per non-overlapping 16-tick window j the estimator sees the projected stiffness k_j (relative sd
s_k), the mean normal force F_j (relative sd 15%: relaxation and damping are not modelled) and, optionally, the mean
centroid offset c_j = My/F (equal to m1/fe in the plant; sd from the moment-sensor grade). Unknowns: log modulus on
a 1 mm grid, the depth d_j per window, a dynamic factor log eta (instantaneous/static stiffness) and a moment bias.
Prior: Gaussian increments of log E between neighbouring cells (sd 0.5 per mm), weak absolute prior.
Fisher information at the truth, causal (windows up to the current one). Reported: posterior sd of ln k at the
current probe position and 5 / 10 mm ahead, in transition windows (true |d ln k| over the window > 0.1) and
elsewhere; and the within-window change of ln k, which the window-constant model ignores.
Caveat (recorded in methods 21.6): the kernels share the plant's elastic form and s_k = 0.05 is optimistic
(replayed windows have relative sd ~0.2); the result is qualitative only.
Run: python -X utf8 p3_identifiability.py
"""
import os
import sys
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
sys.path.insert(0, os.path.join(HERE, '..', 'cobot_scan_benchmark'))
import common, variants, cohort_eval as CE
import plant as PL
import footprint_kernels as FK

NW = 16
CASES = {'k only': dict(moment=None, sk=.05),
         'k + moment (Nano25 grade)': dict(moment=6e-5, sk=.05, bias=1e-3),
         'k + moment (RM65-6F grade)': dict(moment=3e-4, sk=.05, bias=5e-3),
         'k only, s_k 0.10': dict(moment=None, sk=.10)}


def window_rows(log, p):
    t, x, F, d, ktrue = log[:, 0], log[:, 1], log[:, 5], np.maximum(log[:, 6], 1e-6), log[:, 8]
    sel = np.where((t >= 1.) & (F > 1.))[0]; nwin = len(sel)//NW
    if nwin < 10:
        return None
    g = FK.Grid(x[sel].min()-.045, x[sel].max()+.045)
    E = np.array([PL.modulus(xx, p) for xx in g.x])
    rows = []
    h = 2e-6
    for j in range(nwin):
        idx = sel[j*NW:(j+1)*NW]; xc, dc = float(x[idx].mean()), float(d[idx].mean())
        Fe, M1, K = FK.kernels(g, xc, dc)
        Fp, Mp, Kp = FK.kernels(g, xc, dc+h); Fm, Mm, Km = FK.kernels(g, xc, max(dc-h, 1e-9))
        rows.append(dict(xc=xc, d=dc, Fe=Fe, M1=M1, K=K, fe=Fe@E, m1=M1@E, ks=K@E,
                         dfe=(Fp@E-Fm@E)/(2*h), dk=(Kp@E-Km@E)/(2*h),
                         dc=((Mp@E)/(Fp@E)-(Mm@E)/(Fm@E))/(2*h),
                         dlk=abs(np.log(max(ktrue[idx[-1]], 1.))-np.log(max(ktrue[idx[0]], 1.)))))
    return g, E, rows


def analyse(log, p):
    w = window_rows(log, p)
    if w is None:
        return None
    g, E, rows = w; n = g.n; nwin = len(rows); out = {}
    for case, cf in CASES.items():
        nu = n+nwin+2; Jr = []; wr = []
        for j, r in enumerate(rows):
            a = np.zeros(nu); a[:n] = r['K']*E/r['ks']; a[n+j] = r['dk']/r['ks']; a[n+nwin] = 1.
            Jr.append(a); wr.append(1/cf['sk']**2)
            a = np.zeros(nu); a[:n] = r['Fe']*E/r['fe']; a[n+j] = r['dfe']/r['fe']
            Jr.append(a); wr.append(1/.15**2)
            if cf['moment']:
                c = r['m1']/r['fe']; sc = cf['moment']/4.
                a = np.zeros(nu); a[:n] = (r['M1']*E-c*r['Fe']*E)/r['fe']; a[n+j] = r['dc']; a[n+nwin+1] = 1./4.
                Jr.append(a); wr.append(1/sc**2)
        Jr = np.array(Jr); wr = np.array(wr)
        D = np.zeros((n-1, nu)); D[np.arange(n-1), np.arange(n-1)] = -1; D[np.arange(n-1), np.arange(1, n)] = 1
        P0 = D.T@D/.5**2+np.diag(np.r_[np.full(n, 1/3.**2), np.full(nwin, 1/.005**2), 1/.5**2,
                                       1/(cf.get('bias', 1.))**2])
        nobs = 3 if cf['moment'] else 2; res = []
        for j in range(5, nwin):
            Jm = Jr[:nobs*(j+1)]; wm = wr[:nobs*(j+1)]
            C = np.linalg.inv(Jm.T@(Jm*wm[:, None])+P0)
            r = rows[j]; sds = []
            for ahead in (0., .005, .010):
                _, _, K = FK.kernels(g, r['xc']+ahead, r['d'])
                a = np.zeros(nu); a[:n] = K*E/(K@E); a[n+j] = r['dk']/r['ks']; a[n+nwin] = 1.
                sds.append(float(np.sqrt(a@C@a)))
            res.append((r['dlk'] > .1, *sds))
        out[case] = np.array(res, float)
    dl = np.array([r['dlk'] for r in rows[5:]])
    return out, dl


def _ep(a):
    fam, sid, p, T = a
    log, _ = common.run(variants.build('Proposed-IMM', dict(CE.FINAL)), p, T, sid)
    return fam, analyse(log, p)


if __name__ == '__main__':
    sc = [s for s in CE.cohort('diag', 10) if s[0] in ('S2', 'T3m')]
    with Pool(16) as po:
        rows = po.map(_ep, sc, chunksize=1)
    for fam in ('S2', 'T3m'):
        R = [r[1] for r in rows if r[0] == fam and r[1] is not None]
        dl = np.concatenate([r[1] for r in R]); tr = dl > .1
        print('%s: %d windows, %.0f%% in transition; within-window |d ln k| in transitions median %.2f, p90 %.2f'
              % (fam, len(dl), 100*tr.mean(), np.median(dl[tr]), np.percentile(dl[tr], 90)))
        for case in CASES:
            A = np.concatenate([r[0][case] for r in R]); T_ = A[:, 0] > .5
            print('  %-28s sd(ln k) transition: now %.3f +5mm %.3f +10mm %.3f | elsewhere: now %.3f +5mm %.3f +10mm %.3f'
                  % (case, *np.median(A[T_, 1:], 0), *np.median(A[~T_, 1:], 0)), flush=True)
