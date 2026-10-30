"""Does the force-derived equivalent motion need a switching tracker? (paper draft §二 difficulty 4, §五.6)

Equivalent motion seen through force is z = f/k - y; with f = k(t)(y - s), dz = -ds + (y - s) d ln k. It is quiet
during breathing and abrupt during pushes AND during ridge crossings (the d ln k term), so a single-bandwidth
tracker trades noise in quiet segments against lag at events. Two-ends test with every non-tracker parameter
fixed at the final values: single-model Kalman at the IMM's low (q_smooth) and high (q_maneuver) process noise,
the matched-tuned single Kalman (mechanism_checks.json), and the IMM. Force error is split into event segments
(push ramps +-0.1 s; ridge crossings where ln modulus changes by > ln 2 within +-0.1 s) and quiet segments; the
IMM's maneuver-mode probability is reported for both.
Run: python -X utf8 tracker_dilemma.py -> results_v2/tracker_dilemma.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import scenes_v2 as SV
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
KF_TUNED = json.load(open(os.path.join(H, 'results_v2', 'mechanism_checks.json')))['kf_matched']['best']
COMMON = {k: v for k, v in BASE.items() if k not in ('q_smooth', 'q_maneuver', 'p_stay', 'r_scale')}
CONF = {'KF-low': dict(q_kf=BASE['q_smooth'], r_scale=BASE['r_scale']),
        'KF-high': dict(q_kf=BASE['q_maneuver'], r_scale=BASE['r_scale']),
        'KF-tuned': KF_TUNED, 'IMM': None}
FD, WARM, TC = 4.0, 1.0, 0.008


class Logged(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.mu_log = {}

    def drift_pred(self, o, n):
        r = super().drift_pred(o, n)
        if not self.single_model and self.tracker.initialized:
            self.mu_log[n] = float(self.tracker.mu[1])
        return r


def event_mask(log, p):
    t = log[:, 0]; m = np.zeros(len(t), bool)
    for k in range(int(p[15])):
        t0, tr, hold = p[16+4*k], p[18+4*k], p[19+4*k]
        for a in (t0, t0+tr+hold):
            m |= (t >= a-0.1) & (t <= a+tr+0.1)
    lnE = np.log(np.maximum(log[:, 7], 1.))
    w = int(round(0.1/TC))
    for i in range(len(t)):
        seg = lnE[max(0, i-w):i+w+1]
        if seg.max()-seg.min() > math.log(2):
            m[i] = True
    return m


def _ep(a):
    name, (fam, sid, p, T) = a
    if name == 'IMM':
        pol = Logged(single_model=False, **BASE)
    else:
        pol = Logged(single_model=True, **dict(COMMON, **CONF[name]))
    log, pol = common.run(pol, p, T, sid)
    ev = event_mask(log, p); ok = log[:, 0] >= WARM
    e = log[:, 3]-FD
    rms = lambda mk: float(math.sqrt(np.mean(e[mk]**2))) if mk.any() else float('nan')
    mu = np.array([pol.mu_log.get(i, np.nan) for i in range(len(log))])
    mmu = lambda mk: float(np.nanmean(mu[mk])) if np.isfinite(mu[mk]).any() else float('nan')
    return dict(name=name, family=fam, sid=sid, ev=rms(ok & ev), quiet=rms(ok & ~ev), frac_ev=float(np.mean(ev[ok])),
                mu_ev=mmu(ok & ev), mu_quiet=mmu(ok & ~ev), peak=float(np.max(log[ok, 3])))


if __name__ == '__main__':
    sc = SV.scenes('confirm')
    with Pool(16) as po:
        rows = po.map(_ep, [(n, s) for n in CONF for s in sc], chunksize=1)
    out = {}
    for fam in SV.FAMILIES:
        for n in CONF:
            R = [r for r in rows if r['name'] == n and r['family'] == fam]
            g = lambda k: float(np.nanmean([r[k] for r in R]))
            out[fam+'|'+n] = {k: g(k) for k in ('ev', 'quiet', 'frac_ev', 'mu_ev', 'mu_quiet', 'peak')}
            v = out[fam+'|'+n]
            print('%-9s %-8s event %.3f quiet %.3f (event frac %.2f) mu_maneuver event %.2f quiet %.2f peak %.1f'
                  % (fam, n, v['ev'], v['quiet'], v['frac_ev'], v['mu_ev'], v['mu_quiet'], v['peak']), flush=True)
    json.dump(dict(conf={k: v for k, v in CONF.items()}, summary=out, rows=rows),
              open(os.path.join(H, 'results_v2', 'tracker_dilemma.json'), 'w'), indent=1)
