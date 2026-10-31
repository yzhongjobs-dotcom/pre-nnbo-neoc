"""Onset state of every runaway (true force first exceeding 15 N after t = 1 s) on a cohort (diagnosis).

For each scene and configuration: at the onset tick, the measurement ratio kappa_g = k/k_good, the feedback ratio
kappa = k/k_bar, the IMM maneuver probability, the force level, the age of k_good (time since it was last set), and the
precision/consistency of the window that set it, plus kappa_g just before onset (0.1 s earlier).
Run: CONFIGS=final,A2 COHORT=diag FAMS=T2m-stiff python -X utf8 diag_onset.py
"""
import json, os
import numpy as np
from multiprocessing import Pool
import common
import cohort_eval as CE
import variants

H = common.HERE
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']


class Traced(variants.MotionAwareIMM):
    def reset(self, platform):
        super().reset(platform)
        self.tr = {}; self.kg_set = (np.nan, np.nan, np.nan, -1)
        self._kg_last = np.nan

    def step(self, o):
        y = super().step(o)
        if self.k_good != self._kg_last:
            f = self.fit or {}
            self.kg_set = (self.k_good, self.rel, f.get('chi2', np.nan), o.n); self._kg_last = self.k_good
        mu = float(self.tracker.mu[1]) if self.tracker.initialized else np.nan
        self.tr[o.n] = (self.kbar, self.k_good, mu, self.kg_set)
        return y


def _ep(a):
    cfg, (fam, sid, p, T) = a
    extra = dict(CE.extra_configs()).get(cfg, {})
    pol = Traced(single_model=False, **dict(BASE, **extra))
    log, pol = common.run(pol, p, T, sid)
    t, F, k = log[:, 0], log[:, 3], log[:, 8]
    on = np.where((t > 1) & (F > 15))[0]
    if len(on) == 0:
        return None
    i = int(on[0]); j = max(0, i-12)
    kb, kg, mu, (kgv, rel, chi, n_set) = pol.tr.get(i, (np.nan,)*3+((np.nan,)*4,))
    _, kg0, mu0, _ = pol.tr.get(j, (np.nan,)*4)
    return dict(cfg=cfg, sid=sid, t=float(t[i]), F_prev=float(F[j]), kappa=float(k[i]/kb), kappag=float(k[i]/kg),
                kappag_before=float(k[j]/kg0), mu=mu, mu_before=mu0, age=float((i-n_set)*0.008), set_rel=rel, set_chi2=chi,
                k_ratio_since_set=float(k[i]/k[n_set]) if n_set >= 0 else np.nan, peak=float(F[t > 1].max()))


if __name__ == '__main__':
    cfgs = os.environ.get('CONFIGS', 'final').split(',')
    fams = os.environ.get('FAMS', 'T2m-stiff').split(',')
    sc = [s for s in CE.cohort(os.environ.get('COHORT', 'diag'), int(os.environ.get('N', '50'))) if s[0] in fams]
    with Pool(16) as po:
        rows = [r for r in po.map(_ep, [(c, s) for c in cfgs for s in sc]) if r]
    for r in sorted(rows, key=lambda r: (r['cfg'], -r['peak'])):
        print('%-6s %d peak %6.1f t %.2f | kappa_g %.2f (before %.2f) kappa %.2f mu %.2f (before %.2f) | k_good age %.2f s set rel %.2f chi2 %.1f, k changed x%.2f since'
              % (r['cfg'], r['sid'], r['peak'], r['t'], r['kappag'], r['kappag_before'], r['kappa'], r['mu'], r['mu_before'],
                 r['age'], r['set_rel'], r['set_chi2'], r['k_ratio_since_set']))
