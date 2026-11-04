"""Which contact terms help when nothing is privileged: the deployed controller with its own k_good, population constants only.

  gamma     k_good follows the force, k = k_good (F/F_set)^gamma (static stiffening; repository option kg_gamma)
  memory    force-to-motion conversion through a weak power-law memory, beta = .1, three Prony branches that differ from the plant's
  both      the two together
Run on the nominal plant, then the same model on plants whose memory it does not match (robustness to a wrong contact model).
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import physics_plant as PP
import scenes as SC
import sweeps as S
from dynamic_gain_value import ConversionPolicy, stats
from failure_replay import Measured, cell_scene

HERE = Path(__file__).resolve().parent
BASE = 386000
N_SCENE = 5
from ridge_estimator import POP_ALPHA, POP_TAU
GAMMA = .25
VARIANTS = ('base', 'gamma', 'memory', 'both')
PLANTS = {'nominal': {}, 'low_memory': dict(alpha=[.1, .1, .1]), 'high_memory': dict(alpha=[.6, .6, .6]), 'one_branch': dict(alpha=[.5, 0., 0.], tau=[.04, .08, .8])}


class Deployed(ConversionPolicy):
    """Own k_good, population memory; no truth column is read."""
    def reset(self, platform):
        super().reset(platform)
        self.kr = 1e12; self.alpha = POP_ALPHA
        self.e = np.exp(-self.Tc/POP_TAU); self.g = np.exp(-self.Tc/(2*POP_TAU))
        z = np.exp(-1j*2*math.pi*float(self.params.get('f_exc', 20.8))*self.Tc)
        self.r_c = abs(1.+np.sum(POP_ALPHA*self.g*(1-z)/(1-self.e*z)))     # the response of the recursion below

    def _next_z(self):
        o, n = self._ctx; j = n-1
        decay = (1-self.afc)**self.ns; denom = 1-decay
        fin = (o.log[j, 5]-decay*o.log[j-1, 5])/denom
        y0, y1 = o.log[j-1, 2], o.log[j, 2]
        y = float(y0+self._encoder_fraction*(y1-y0))
        dF = fin-self._fin_prev; dy = y-self._y_prev
        if min(o.log[j, 5], o.log[j-1, 5]) <= .5:
            self.q = np.zeros(3)      # out of contact, by the estimators' test
        kc = self.k_good
        g = self.params.get('kg_gamma', 0.)
        if g and np.isfinite(kc) and np.isfinite(getattr(self, 'f_set', np.nan)):
            kc = kc*(max(o.log[j, 5], .3)/self.f_set)**g
        if not np.isfinite(kc) or kc <= 0:
            kc = self.kbar
        ke = kc/self.r_c
        dFe = (dF-float(np.sum((self.e-1)*self.q)))/(1+float(np.sum(self.alpha*self.g)))
        self.q = self.e*self.q+self.alpha*self.g*dFe
        self.Zd += dFe/ke-dy
        self.Z = self.Zd
        return self.Zd


def make(variant, controls, est):
    p = dict(S.FINAL_C); p.update(controls)
    if variant in ('gamma', 'both'):
        p['kg_gamma'] = GAMMA
    e = S.estimator('Projection', est)
    return Deployed(e, single_model=False, **p) if variant in ('memory', 'both') else S.SwapPolicy(e, single_model=False, **p)


def job(a):
    plant, variant, c, t, i, controls, est = a
    sc = cell_scene(c, t, i, BASE, PLANTS[plant]); seed = BASE+sum(ord(x) for x in c+t)*100+i
    pol = make(variant, controls, est); pol.sc = sc
    log = PP.episode(sc, Measured(pol), seed)
    return dict(plant=plant, variant=variant, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, controls.get('f_exc', 20.8)))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']['Projection']['best']['params']
    w = 2*math.pi*np.array([1., 20.8]); r = np.abs(1+np.sum(POP_ALPHA*1j*w[:, None]*POP_TAU/(1+1j*w[:, None]*POP_TAU), axis=1))
    print('population memory: gain 1 Hz / 20.8 Hz', round(float(r[0]/r[1]), 3))
    plants = sys.argv[1:] or list(PLANTS)
    args = [(pl, v, c, ti, i, controls, est) for pl in plants for v in VARIANTS for c, ti in S.CELLS for i in range(N_SCENE)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, args, chunksize=1))
    (HERE/'results'/'deployable_contact.json').write_text(json.dumps(rows, indent=0, allow_nan=True), encoding='utf-8')
    L = lambda m: math.log1p(min(m['J'], 1000)); key = lambda r: (r['cond'], r['tissue'], r['i'])
    for pl in plants:
        base = {key(r): r['metrics'] for r in rows if r['plant'] == pl and r['variant'] == 'base'}
        print(f'\n[{pl}]')
        for v in VARIANTS:
            rr = [r for r in rows if r['plant'] == pl and r['variant'] == v]; s = stats(rr)
            d = np.array([L(r['metrics'])-L(base[key(r)]) for r in rr]); m = d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1)
            print(f"{v:7s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.3f} event {s['event_med']:.3f} fail {s['fail']} long {s['long_loss']} | vs base {d.mean():+.3f} [{np.percentile(m,2.5):+.3f},{np.percentile(m,97.5):+.3f}]")
            print('        ', ' '.join(f"{g}:{stats([r for r in rr if r['tissue']==g])['L']:.2f}" for g in ('soft', 'mid', 'hard')), ' '.join(f"{n}:{stats([r for r in rr if r['cond'] in cs])['L']:.2f}" for n, cs in (('event', ('flinch', 'cough')), ('fast', ('shift_fast',)), ('smooth', ('breath', 'shift_slow')), ('rib+combo', ('rib', 'combo')))))


if __name__ == '__main__':
    main()
