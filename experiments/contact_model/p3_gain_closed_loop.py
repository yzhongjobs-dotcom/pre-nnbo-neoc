"""Closed-loop check of the measurement-gain filter (methods 21.6): adopted method (final + V1 servo identification)
versus the same with the measurement gain from the two-model IMM gain filter (tv_projection.GainIMM) instead of the
hard accept-and-hold rule (b <= ln 1.2). Feedback stiffness, excitation, tracker and everything else unchanged.
Diagnosis cohort (110000+), 50 per family, nominal platform; paired bootstrap (4000) of per-scene RMSE.
Run: python -X utf8 p3_gain_closed_loop.py [N] -> results/p3_gain_closed_loop.txt
"""
import math
import os
import sys
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'motion_aware_scan'))
sys.path.insert(0, os.path.join(HERE, '..', 'method_redesign'))
import common, variants, cohort_eval as CE
from v1_servo_id import ServoID
from tv_projection import GainIMM


class GainIMMMethod(ServoID):
    imm_kw = {}; raise_rule = False; edge_gate = False; my_noise = 7.5e-4; my_extra = 0.

    def _edge(self, o, n):
        """pressure moment of the finite footprint: an edge under a symmetric footprint gives a first moment, normal
        motion on homogeneous tissue does not. Window-mean moment beyond 3 sd of the sensor's window noise (sensor
        spec: per-sample sd my_noise at 1 kHz, 30 Hz filter ~0.3, 16 ticks ~ 16 independent samples)."""
        My = o.log[max(n-16, 0):n, 13]
        if self.my_extra:
            My = My+self.my_extra*.3*self._rng.standard_normal(len(My))
        sd_w = math.hypot(self.my_noise, self.my_extra)*.3/4.
        return abs(float(np.mean(My))) > 3*sd_w

    def reset(self, platform):
        super().reset(platform); self.gimm = GainIMM(**self.imm_kw); self._rng = np.random.default_rng(1)

    def _stiffness_est(self, o, n):
        out = super()._stiffness_est(o, n)
        f = self.fit
        if f is not None and np.isfinite(f.get('k', np.nan)) and f['info'] > 0 and f['k'] > 0:
            rel = 1./math.sqrt(f['info'])/f['k']
            kh = min(max(f['k'], 200.), 3e5)
            v = self.gimm.step(kh, rel, f['chi2'])
            if self.raise_rule and v is not None and rel <= self.params['rel_max'] and math.log(kh)-v > 2*rel                     and (not self.edge_gate or self._edge(o, n)):
                # under-estimating the measurement gain is the dangerous direction (kappa_g > 1 tolerates only ln 1.2):
                # a precise window clearly above the fused level re-anchors it upward (same asymmetry as the feedback)
                for x in self.gimm.x: x[0] = math.log(kh)
                v = math.log(kh)
        else:
            v = self.gimm.step(np.nan, np.nan, np.nan)
        # the fused gain becomes the measurement gain once its own uncertainty is inside the band the feed-forward
        # loop tolerates (ln 1.2); before that k_good stays undefined as in the original method (startup protection)
        if v is not None and (np.isfinite(self.k_good) or self.gimm.sd() <= math.log(1.2)):
            self.k_good = math.exp(v); self.f_set = max(float(o.log[n-1, 5]), .3)
        return out


def build(cfg):
    base = dict(CE.FINAL)
    if '@a' in cfg:
        cfg, a = cfg.split('@a'); base['a_force'] = float(a)
    if cfg == 'Adopted':
        pol = variants.build('Proposed-IMM', base); pol.__class__ = ServoID
    else:
        pol = variants.build('Proposed-IMM', base); pol.__class__ = GainIMMMethod
        if cfg.endswith('+RM65'): cfg = cfg[:-5]; pol.my_extra = 4*7.5e-4
        if cfg.endswith('+edge'): cfg = cfg[:-5]; pol.raise_rule = True; pol.edge_gate = True
        if cfg.endswith('+raise'): cfg = cfg[:-6]; pol.raise_rule = True
        pol.imm_kw = dict(p_stay=float(cfg.split('-')[1])) if '-' in cfg else {}
    return pol


def _ep(a):
    cfg, (fam, sid, p, T) = a
    log, _ = common.run(build(cfg), p, T, sid); m = common.task_metrics(log, 20.8)
    return cfg, fam, sid, (m['rmse'] if not m['diverged'] else np.nan), m['peak'], m['loss_max'], m['J']


CE.SEEDS['tune'] = 190000

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'tune':
    CFG = [m+'@a%g' % a for m in (sys.argv[2:] or ['Adopted', 'GainIMM-0.98']) for a in (.10, .12, .14, .17)]
    sc = CE.cohort('tune', 10)
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    out = []
    for c in CFG:
        R = [r for r in rows if r[0] == c]; J = np.mean([np.log1p(r[6]) if np.isfinite(r[3]) else np.log1p(1e3) for r in R])
        fam_med = ' '.join('%s %.3f' % (f, np.nanmedian([r[3] for r in R if r[1] == f])) for f in CE.FAMS)
        out.append('%-18s mean log(1+J) %.4f | %s' % (c, J, fam_med)); print(out[-1], flush=True)
    open(os.path.join(HERE, 'results', 'p3_gain_tune_aforce.txt'), 'w', encoding='utf-8').write('\n'.join(out)+'\n')
    sys.exit()

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    CFG = tuple(sys.argv[2:]) if len(sys.argv) > 2 else ('Adopted', 'GainIMM-0.98')
    sc = CE.cohort('diag', N)
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in CFG for s in sc], chunksize=1)
    rng = np.random.default_rng(0); out = []
    for fam in CE.FAMS:
        base = {r[2]: r for r in rows if r[0] == CFG[0] and r[1] == fam}
        for c in CFG:
            R = [r for r in rows if r[0] == c and r[1] == fam]; v = np.array([r[3] for r in R])
            line = '%-9s %-13s %.3f/%.3f r%d L%d' % (fam, c, np.nanmedian(v), np.nanmean(v),
                                                    sum((not np.isfinite(r[3])) or r[4] > 50 for r in R), sum(r[5] > .5 for r in R))
            if c != CFG[0]:
                d = np.array([r[3]-base[r[2]][3] for r in R]); d = d[np.isfinite(d)]
                bs = [np.mean(rng.choice(d, len(d))) for _ in range(4000)]
                line += ' | mean diff %+.3f [%+.3f, %+.3f]' % (d.mean(), *np.percentile(bs, [2.5, 97.5]))
            out.append(line); print(line, flush=True)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    open(os.path.join(HERE, 'results', 'p3_gain_closed_loop.txt'), 'w', encoding='utf-8').write('\n'.join(out)+'\n')
