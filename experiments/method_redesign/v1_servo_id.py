"""V1: online identification of the position-servo lag and command delay from encoder data.
Model (per control tick, command constant within a tick, acting ndel ticks after it is sent):
    y_{j+1} - c_j = rho (y_j - c_j),   c_j = command sent at tick j - ndel,   tau_hat = -Tc / ln(rho).
rho is a windowed least-squares estimate over the last W ticks; ndel in {1,2,3} is the one with the smallest residual.
Every place where the controller uses the servo model is updated: intra-period displacement interpolation (window
fit and motion measurement), the servo-lag velocity term of the feed-forward, and the delay-dependent horizons.
The force filter is taken from the sensor specification (assumption stated; its mismatch is still evaluated).
Run: python -X utf8 v1_servo_id.py
"""
import math, os, sys, json
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
import common, variants, cohort_eval as CE, robustness as RB
from common import P

class ServoID(variants.MotionAwareIMM):
    W = 150                       # 1.2 s of data
    def reset(self, platform):
        super().reset(platform); self._id_n = 0; self.tau_hat = self.tau; self.ndel_hat = self.ndel; self.id_hist = []
    def _apply(self, tau, ndel):
        self.tau = tau; self.eb = math.exp(-1e-3/tau) if tau > 0 else 0.; self.ndel = int(ndel)
        m = np.arange(self.ns)
        frac = -np.expm1(-m*.001/tau)/(-math.expm1(-self.Tc/tau)) if tau > 0 else (m > 0).astype(float)
        w = self.afc*(1-self.afc)**(self.ns-1-m); self._encoder_fraction = float(w@frac/w.sum())
    def _identify(self, log, n):
        best = None
        for nd in (1, 2, 3):
            j = np.arange(max(nd+1, n-1-self.W), n-1)
            if len(j) < 40: return
            y0, y1, c = log[j, 2], log[j+1, 2], log[j-nd, 11]
            u = y0-c; v = y1-c
            den = float(u@u)
            if den < 1e-12: return                                 # no information (commands equal to positions)
            rho = float(u@v)/den; res = float(((v-rho*u)**2).mean())
            if best is None or res < best[0]: best = (res, rho, nd)
        res, rho, nd = best
        if 0. < rho < 1.:
            tau = -self.Tc/math.log(rho)
            tau = min(max(tau, 1e-3), 0.1)
            self.tau_hat, self.ndel_hat = tau, nd
            self._apply(tau, nd)
    def step(self, o):
        n = o.n
        if n > 60 and n % 25 == 0: self._identify(o.log, n)
        if n % 125 == 0: self.id_hist.append((float(o.log[n-1, 0]), self.tau_hat, self.ndel_hat))
        return super().step(o)

COND = {'nominal': {}, 'lag x0.6': dict(tau_a=0.019*0.6), 'lag x1.5': dict(tau_a=0.019*1.5), 'delay +1': dict(ndel=2)}
def _ep(a):
    cond, cfg, (fam, sid, p, T) = a
    plat = dict(RB.NOMINAL); plat.update(COND[cond])
    if cfg == 'nominal-model': pol = RB.NominalModel(variants.build('Proposed-IMM', dict(CE.FINAL)), seed=sid)
    elif cfg == 'told': pol = variants.build('Proposed-IMM', dict(CE.FINAL))
    else: pol = RB.NominalModel(ServoID(single_model=False, **CE.FINAL), seed=sid)
    log = P.episode(p, pol, T=T, platform=plat, seed=sid)
    if log is None: return cond, cfg, fam, np.nan, np.inf, None
    m = common.task_metrics(log, 20.8)
    idh = getattr(getattr(pol, 'inner', pol), 'id_hist', None)
    return cond, cfg, fam, m['rmse'], m['peak'], (idh[-1] if idh else None)

if __name__ == '__main__':
    sc = CE.cohort('diag', 20); CFG = ('nominal-model', 'servo-ID', 'told')
    with Pool(16) as po: rows = po.map(_ep, [(c, g, s) for c in COND for g in CFG for s in sc], chunksize=1)
    json.dump([r[:5]+(r[5],) for r in rows], open(os.path.join(HERE, 'results', 'v1_servo_id.json'), 'w'), indent=1, default=float)
    for c in COND:
        ids = [r[5] for r in rows if r[0] == c and r[1] == 'servo-ID' and r[5]]
        print('%s: identified tau median %.4f s (true %.4f), ndel mode %d' % (c, np.median([i[1] for i in ids]), dict(RB.NOMINAL, **COND[c])['tau_a'],
              int(np.bincount([i[2] for i in ids]).argmax())))
        for fam in CE.FAMS:
            out = []
            for g in CFG:
                r = [x for x in rows if x[0] == c and x[1] == g and x[2] == fam]; v = [x[3] for x in r if np.isfinite(x[3])]
                out.append('%s %.3f/%.3f r%d' % (g, np.median(v), np.mean(v), sum((not np.isfinite(x[3])) or x[4] > 50 for x in r)))
            print('   %-10s %s' % (fam, ' | '.join(out)), flush=True)
