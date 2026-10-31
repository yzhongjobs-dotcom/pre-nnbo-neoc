"""Where do the residual peaks and long contact losses on stiff pushes come from? (diagnosis, paper draft §五)

Runs the final Proposed-IMM on the confirmation scenes of T2m-stiff (S3c) and T2f-stiff (S3d) with per-tick traces
of the feedback stiffness k_bar, the measurement stiffness k_good, the IMM maneuver probability and the supervisor
mode, and classifies each peak event (true force > 10 N) and each contact-loss episode by the state in the 0.2 s
before it: feedback ratio kappa = k/k_bar, measurement ratio kappa_g = k/k_good, push phase (ramp / hold / retract),
and whether it follows a re-contact.
Run: python -X utf8 diag_stiff.py -> results_v2/diag_stiff.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import scenes_v2 as SV
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
TC = 0.008


class Traced(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.tr = []

    def step(self, o):
        y = super().step(o)
        mu = float(self.tracker.mu[1]) if self.tracker.initialized else np.nan
        self.tr.append((o.n, self.kbar, self.k_good, mu, self.mode))
        return y


def phase(t, p):
    for k in range(int(p[15])):
        t0, tr, hold = p[16+4*k], p[18+4*k], p[19+4*k]
        if t0-0.05 <= t <= t0+tr: return 'push-ramp'
        if t0+tr < t <= t0+tr+hold: return 'hold'
        if t0+tr+hold < t <= t0+2*tr+hold+0.05: return 'retract-ramp'
    return 'quiet'


def _ep(sc):
    fam, sid, p, T = sc
    pol = Traced(single_model=False, **BASE)
    log, pol = common.run(pol, p, T, sid)
    n = len(log); kb = np.full(n, np.nan); kg = np.full(n, np.nan); mu = np.full(n, np.nan); md = np.zeros(n)
    for (i, a, b, c, d) in pol.tr:
        if i < n: kb[i], kg[i], mu[i], md[i] = a, b, c, d
    t, F, k = log[:, 0], log[:, 3], log[:, 8]
    ok = t >= 1.0
    events = []
    # peak events: local maxima above 10 N separated by 0.3 s
    idx = [i for i in range(1, n-1) if ok[i] and F[i] > 10 and F[i] >= F[i-1] and F[i] >= F[i+1]]
    last = -1e9
    for i in idx:
        if t[i]-last < 0.3: continue
        last = t[i]; w = slice(max(0, i-25), i+1)
        kap = k[w]/np.maximum(kb[w], 1); kapg = k[w]/np.maximum(kg[w], 1)
        lost_before = bool(np.any(F[max(0, i-63):i] < 0.3))
        events.append(dict(kind='peak', t=float(t[i]), F=float(F[i]), phase=phase(t[i], p), recontact=lost_before,
                           kappa_max=float(np.nanmax(kap)), kappag_max=float(np.nanmax(kapg)), kappag_min=float(np.nanmin(kapg)),
                           mu=float(np.nanmean(mu[w]))))
    # loss episodes
    loss = ok & (F < 0.5); i = 0
    while i < n:
        if loss[i]:
            j = i
            while j < n and loss[j]: j += 1
            w = slice(max(0, i-25), i+1)
            kap = k[w]/np.maximum(kb[w], 1); kapg = k[w]/np.maximum(kg[w], 1)
            events.append(dict(kind='loss', t=float(t[i]), dur=float((j-i)*TC), phase=phase(t[i], p),
                               kappa_max=float(np.nanmax(kap)), kappag_max=float(np.nanmax(kapg)), kappag_min=float(np.nanmin(kapg)),
                               mu=float(np.nanmean(mu[w])), supervisor_time=float(np.mean(md[i:j] == 1)) if j > i else 0.))
            i = j
        else:
            i += 1
    return dict(family=fam, sid=sid, events=events)


if __name__ == '__main__':
    sc = [s for s in SV.scenes('confirm', ('T2m-stiff', 'T2f-stiff'))]
    with Pool(16) as po:
        res = po.map(_ep, sc)
    json.dump(res, open(os.path.join(H, 'results_v2', 'diag_stiff.json'), 'w'), indent=1)
    for fam in ('T2m-stiff', 'T2f-stiff'):
        ev = [e for r in res if r['family'] == fam for e in r['events']]
        pk = [e for e in ev if e['kind'] == 'peak']; ls = [e for e in ev if e['kind'] == 'loss']
        print('==', fam, 'peaks', len(pk), 'losses', len(ls), 'loss>0.2s', sum(e['dur'] > 0.2 for e in ls))
        for e in sorted(pk, key=lambda e: -e['F'])[:10]:
            print(' peak %5.1f N %-12s recontact=%s kappa_max %.2f kappag [%.2f,%.2f] mu %.2f' % (
                e['F'], e['phase'], e['recontact'], e['kappa_max'], e['kappag_min'], e['kappag_max'], e['mu']))
        for e in sorted(ls, key=lambda e: -e['dur'])[:10]:
            print(' loss %.2f s %-12s kappa_max %.2f kappag [%.2f,%.2f] mu %.2f sup %.2f' % (
                e['dur'], e['phase'], e['kappa_max'], e['kappag_min'], e['kappag_max'], e['mu'], e['supervisor_time']))
        from collections import Counter
        print(' peak phases', Counter(e['phase'] for e in pk), ' loss phases', Counter(e['phase'] for e in ls))
