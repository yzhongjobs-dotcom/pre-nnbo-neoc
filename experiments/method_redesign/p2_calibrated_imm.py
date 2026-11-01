"""Calibrated motion tracker: the IMM measurement variance is estimated online from the motion measurement itself
(causal second difference d2 = Z_n - 2 Z_{n-1} + Z_{n-2}, var = 6 R for white noise; short-time power, noise floor by
minimum statistics with bias factor B) instead of the tuned R = r (sigma_F/k)^2; the consistency-gate inflation is
kept. Process noise re-searched on a 3x3 grid on the tuning seeds (190000+, 10 per family), objective = mean
log(1+J) as in the tuning protocol. Run: python -X utf8 p2_calibrated_imm.py [grid|diag NAME]
"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE
CE.SEEDS['tune'] = 190000
F = dict(CE.FINAL)
class CalIMM(variants.MotionAwareIMM):
    cal_T, cal_win, cal_B = .1, 1.5, 2.
    def reset(self, platform):
        super().reset(platform)
        tr = self.tracker; up, rs = tr.update, tr.reset; st = self._cal = dict(z=[], P=None, N=np.inf, k=0)
        def update_(z):
            st['z'] = (st['z']+[z])[-3:]
            if len(st['z']) == 3:
                d2 = st['z'][2]-2*st['z'][1]+st['z'][0]; a = 1-np.exp(-self.Tc/self.cal_T)
                st['P'] = d2*d2/6 if st['P'] is None else st['P']+a*(d2*d2/6-st['P']); st['k'] += 1
                st['N'] = st['P'] if st['k'] < 12 else min(st['P'], st['N']*np.exp(self.Tc/self.cal_win*np.log(2.)))
            if st['k'] >= 12:
                gate = 1.
                if self.gate and self.fit is not None and np.isfinite(self.fit.get('chi2', np.nan)):
                    gate = max(1., self.fit['chi2']/self.params['chi2_max'])
                tr.R = max(self.cal_B*st['N']*gate, np.finfo(float).tiny)
            return up(z)
        def reset_(z0):
            st['z'] = []; st['P'] = None; st['N'] = np.inf; st['k'] = 0; return rs(z0)
        tr.update, tr.reset = update_, reset_
def mk(name):
    if name == 'Proposed': return variants.build('Proposed-IMM', F)
    a, b = name.split('_')[1:]
    pol = variants.build('Proposed-IMM', dict(F, q_smooth=F['q_smooth']*float(a), q_maneuver=F['q_maneuver']*float(b)))
    pol.__class__ = CalIMM; return pol
GRID = ['Proposed']+['Cal_%g_%g' % (a, b) for a in (1, 10, 100) for b in (1, 10, 100)]
def _ep(a):
    name, (fam, sid, p, T) = a
    log, _ = common.run(mk(name), p, T, sid); m = common.task_metrics(log, 20.8)
    return name, fam, m['rmse'], m['peak'], m['loss_max'], m['J'], not m['diverged']
def table(rows, names):
    for nm in names:
        line = '%-14s' % nm; J = []
        for fam in CE.FAMS:
            R = [r for r in rows if r[0] == nm and r[1] == fam]; ok = [r for r in R if r[6]]
            rm = [r[2] for r in ok]; J += [np.log1p(r[5]) if r[6] else np.log1p(1e3) for r in R]
            line += ' | %s %.3f/%.3f r%d L%d' % (fam, np.median(rm), np.mean(rm), len(R)-len(ok)+sum(r[3] > 50 for r in ok), sum(r[4] > .5 for r in ok))
        print('%s | mean log(1+J) %.3f' % (line, np.mean(J)), flush=True)
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'grid'
    if mode == 'grid':
        sc = CE.cohort('tune', 10); names = GRID
    else:
        sc = CE.cohort('diag', 20); names = ['Proposed']+sys.argv[2:]
    with Pool(16) as po: rows = po.map(_ep, [(n, s) for n in names for s in sc], chunksize=1)
    table(rows, names)
