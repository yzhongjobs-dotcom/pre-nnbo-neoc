"""Which contact model the estimator and the force-to-motion conversion should assume, and which parts of the candidate matter.

matrix    Local contact forms, each used in both channels, against plants whose contact law differs, so that no form wins just because
          the plant was written that way. These are local forms for this interface, not full implementations of the cited papers:
            elastic      F = k d, the static self-consistent normal law of the AuSoScan cylinder (d = d0, v = 0); its sliding term is
                         below .1% of F at scan speed. Ridge window of the candidate, not the existing P5 interface
            gamma 1/3    k rises with F as for a Hertz sphere (T-RO 2023 geometry), gain scheduling only
            hc form      Hunt-Crossley rate term with n = 1.26 from Moreira 2014, so gamma = .206; lam = .25 is an assumption, no value in the sources
            sls form     one relaxation branch, alpha = .6, tau = 50 ms, an assumption of the Kelvin-Boltzmann form
            weak power   three Prony branches approximating gain slope .1 (r = .125 in Moreira 2014), a recursive approximation of D^r
            power+gamma  weak power with gamma = .206
          Control parameters and gate scales are fixed for all, so this compares fixed configurations, not each form at its best.
ablation  The full candidate with one part removed at a time, on unused seeds, with every runaway classified.
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import physics_plant as PP
from dynamic_gain_value import Truth, stats
from estimator_swap import SwapPolicy
from failure_replay import cell_scene
from ridge_estimator import RidgeProjection, POP_ALPHA, POP_TAU
from root_cause import classify
from ridge_final import CENTRE

HERE = Path(__file__).resolve().parent
TUNING = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
CTL = TUNING['public']['C_neutral']['best']['params']
P5 = dict(nw=16, degree=5, sigma_scale=.796, chi_scale=.551)
GATE = dict(sigma_scale=.796, chi_scale=.551)
RIDGE = dict(nw=24, order=4, fc=14., **GATE)
MATRIX_BASE = 411000; N_MATRIX = 5
ABLATE_BASE = 413000; N_ABLATE = 10

MODELS = {
    'elastic': dict(),
    'gamma 1/3': dict(gamma=1/3),
    'hc form': dict(gamma=.206, lam=.25),
    'sls form': dict(alpha=np.array([.6]), tau=np.array([.05]), memory=True),
    'weak power': dict(alpha=POP_ALPHA, tau=POP_TAU, memory=True),
    'power+gamma': dict(alpha=POP_ALPHA, tau=POP_TAU, memory=True, gamma=.206),
    'P5 baseline': dict(baseline=True),
}
PLANTS = {
    'footprint': {},
    'power p=1': dict(law=1.),
    'power p=1.5': dict(law=1.5),
    'hc p=1.26 lam 1 no mem': dict(law=1.26, overrides=dict(alpha=[0., 0., 0.]), lam=1.),
    'high memory': dict(overrides=dict(alpha=[.6, .6, .6])),
    'one branch': dict(overrides=dict(alpha=[.5, 0., 0.], tau=[.04, .08, .8])),
}
# gamma in one channel only, on the plant where it is physically true and on the nominal plant
GAMMA_SPLIT = {'gamma regressor only': dict(reg=1/3, conv=0.), 'gamma conversion only': dict(reg=0., conv=1/3)}


class ModelPolicy(SwapPolicy):
    """Candidate controller: conversion through the given contact model; centring and first-gain acquisition can be removed."""
    centre = True; acquire = True; zalpha = None; ztau = None; gamma = 0.; lam = 0.

    def reset(self, platform):
        super().reset(platform)
        self.q = np.zeros(0 if self.zalpha is None else len(self.zalpha)); self.Zd = 0.; self._ctx = None; self.dc_prev = 0.
        if self.zalpha is not None:
            self.e = np.exp(-self.Tc/self.ztau); self.b = self.ztau/self.Tc*(1-self.e)    # exact for a ramp within the tick
            z = np.exp(-1j*2*math.pi*float(self.params.get('f_exc', 20.8))*self.Tc)
            self.r_c = abs(1.+np.sum(self.zalpha*self.b*(1-z)/(1-self.e*z)))
        else:
            self.r_c = 1.
        self._wrap()

    def _wrap(self):
        if getattr(self.tracker, '_conversion', None) is self:
            return
        original = self.tracker.update
        self.tracker.update = lambda _: original(self._next_z())
        self.tracker._conversion = self

    def _next_z(self):
        o, n = self._ctx; j = n-1
        decay = (1-self.afc)**self.ns; denom = 1-decay
        F = float(o.log[j, 5]); fin = (F-decay*o.log[j-1, 5])/denom
        y0, y1 = o.log[j-1, 2], o.log[j, 2]
        y = float(y0+self._encoder_fraction*(y1-y0))
        dF = fin-self._fin_prev; dy = y-self._y_prev
        if min(F, o.log[j-1, 5]) <= .5:
            self.q = np.zeros_like(self.q); self.dc_prev = 0.
        kc = self.k_good
        if self.gamma and np.isfinite(kc) and np.isfinite(getattr(self, 'f_set', np.nan)):
            kc = kc*(max(F, .3)/self.f_set)**self.gamma
        if not np.isfinite(kc) or kc <= 0:
            kc = self.kbar
        if self.zalpha is not None:
            dFe = (dF-float(np.sum((self.e-1)*self.q)))/(1+float(np.sum(self.zalpha*self.b)))
            self.q = self.e*self.q+self.zalpha*self.b*dFe
            dc = dFe*self.r_c/kc
        else:
            dc = dF/kc
        if self.lam:
            kt = kc
            for _ in range(2):      # kc is the magnitude at the probe frequency, the time constant needs the elastic tangent
                tau = self.lam*max(F, .3)/kt; kt = kc/math.sqrt(1+(2*math.pi*20.8*tau)**2)
            a = tau/self.Tc         # rate term: a first-order lag on the inverse
            dc = (dc+a*self.dc_prev)/(1+a); self.dc_prev = dc
        self.Zd += dc-dy
        self.Z = self.Zd
        return self.Zd

    def drift_pred(self, o, n):
        self._wrap()
        self._ctx = (o, n)
        out = super().drift_pred(o, n)
        if self.Z == 0. or self._fin_prev is None:
            self.Zd = 0.; self.q = np.zeros_like(self.q); self.dc_prev = 0.
        return out

    def _stiffness_est(self, o, n):
        p = self.params; before = self.k_good
        acquiring = self.acquire and not np.isfinite(self.k_good)
        saved = p['chi2_max'], p['meas_bound']
        if acquiring:
            p['chi2_max'] = np.inf; p['meas_bound'] = None
        out = super()._stiffness_est(o, n)
        p['chi2_max'], p['meas_bound'] = saved
        if acquiring and not np.isfinite(self.k_good) and np.isfinite(getattr(self, 'rel', np.nan)) and self.rel <= p['rel_max']:
            self.k_seek = self.k_seek/.85
        if self.centre and np.isfinite(self.k_good) and self.k_good != before:
            self.k_good *= CENTRE
        return out


def scene(c, t, i, base, plant):
    spec = PLANTS[plant]
    sc = cell_scene(c, t, i, base, spec.get('overrides'))
    if 'law' in spec:
        sc['law'] = spec['law']
    if 'lam' in spec:
        sc['lam'] = spec['lam']
    return sc


def build(model):
    p = dict(S.FINAL_C); p.update(CTL)
    if model == 'P5 baseline':
        pol = ModelPolicy(S.estimator('Projection', P5), single_model=False, **p)
        pol.centre = False; pol.acquire = False
        return pol
    if model in GAMMA_SPLIT:
        g = GAMMA_SPLIT[model]
        pol = ModelPolicy(RidgeProjection(gamma=g['reg'], **RIDGE), single_model=False, **p); pol.gamma = g['conv']
        return pol
    m = dict(MODELS[model])
    gamma = m.get('gamma', 0.); lam = m.get('lam', 0.)
    est = RidgeProjection(memory=m.get('memory', False), alpha=m.get('alpha'), tau=m.get('tau'), gamma=gamma, lam=lam, **RIDGE)
    pol = ModelPolicy(est, single_model=False, **p)
    pol.zalpha = m.get('alpha'); pol.ztau = m.get('tau'); pol.gamma = gamma; pol.lam = lam
    return pol


def matrix_job(a):
    plant, model, c, t, i = a
    sc = scene(c, t, i, MATRIX_BASE, plant); seed = MATRIX_BASE+sum(ord(x) for x in c+t)*100+i
    pol = build(model)
    log = PP.episode(sc, Truth(pol), seed)
    return dict(plant=plant, model=model, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, 20.8))


ABLATE = {
    'Full': dict(),
    'no relaxed prior (P5_16)': dict(estimator='p5'),
    'no y memory': dict(ymem=False),
    'no z memory': dict(zmem=False),
    'no centring': dict(centre=False),
    'no acquisition rule': dict(acquire=False),
    'D2 prior': dict(order=2),
    'relaxed prior only': dict(ymem=False, zmem=False, centre=False),
    'P5_16 only': dict(estimator='p5', zmem=False, centre=False, acquire=False),
}


def ablate_job(a):
    name, c, t, i = a
    f = ABLATE[name]
    sc = cell_scene(c, t, i, ABLATE_BASE); seed = ABLATE_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    if f.get('estimator') == 'p5':
        est = S.estimator('Projection', P5)
    else:
        est = RidgeProjection(memory=f.get('ymem', True), **dict(RIDGE, order=f.get('order', 4)))
    pol = ModelPolicy(est, single_model=False, **p)
    if f.get('zmem', True):
        pol.zalpha = POP_ALPHA; pol.ztau = POP_TAU
    pol.centre = f.get('centre', True); pol.acquire = f.get('acquire', True)
    w = Truth(pol); tr = []

    class Trace:
        def reset(self, platform):
            w.reset(platform)

        def __call__(self, o):
            u = w(o); tr.append((pol.kbar, pol.k_good)); return u
    log = PP.episode(sc, Trace(), seed)
    out = dict(name=name, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, 20.8), cause=None)
    out['cause'] = dict(kind='numeric') if log is None else classify(log, np.array(tr, float), sc)
    return out


L = lambda m: math.log1p(min(m['J'], 1000))
key = lambda r: (r['cond'], r['tissue'], r['i'])
boot = lambda d: np.percentile(d[np.random.default_rng(1).integers(0, len(d), (10000, len(d)))].mean(1), [2.5, 97.5])
GROUPS = (('soft', lambda r: r['tissue'] == 'soft'), ('mid', lambda r: r['tissue'] == 'mid'), ('hard', lambda r: r['tissue'] == 'hard'),
          ('event', lambda r: r['cond'] in ('flinch', 'cough')), ('fast', lambda r: r['cond'] == 'shift_fast'),
          ('smooth', lambda r: r['cond'] in ('breath', 'shift_slow')), ('rib+combo', lambda r: r['cond'] in ('rib', 'combo')))


def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else 'all'
    res = {}
    with ProcessPoolExecutor(max_workers=15) as ex:
        if phase in ('all', 'matrix'):
            extra = [(pl, m, c, t, i) for pl in ('footprint', 'power p=1.5') for m in GAMMA_SPLIT for c, t in S.CELLS for i in range(N_MATRIX)]
            rows = list(ex.map(matrix_job, [(pl, m, c, t, i) for pl in PLANTS for m in MODELS for c, t in S.CELLS for i in range(N_MATRIX)]+extra, chunksize=1))
            res['matrix'] = {}
            print('plant | model: L, rmse med, event, fail | vs elastic')
            for pl in PLANTS:
                ref = {key(r): r['metrics'] for r in rows if r['plant'] == pl and r['model'] == 'elastic'}
                res['matrix'][pl] = {}
                for m in MODELS:
                    rr = [r for r in rows if r['plant'] == pl and r['model'] == m]; s = stats(rr)
                    d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(d)
                    res['matrix'][pl][m] = dict(s, vs_elastic=[float(d.mean()), float(c[0]), float(c[1])])
                    print(f"{pl:20s} {m:11s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.2f} event {s['event_med']:.3f} fail {s['fail']:2d} long {s['long_loss']:2d} | {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]", flush=True)
            print('mean L over plants', {m: round(float(np.mean([res['matrix'][pl][m]['L'] for pl in PLANTS])), 3) for m in MODELS}, flush=True)
            for pl in ('footprint', 'power p=1.5'):
                ref = {key(r): r['metrics'] for r in rows if r['plant'] == pl and r['model'] == 'elastic'}
                for m in GAMMA_SPLIT:
                    rr = [r for r in rows if r['plant'] == pl and r['model'] == m]; s_ = stats(rr)
                    d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(d)
                    print(f"{pl:20s} {m:22s} L {s_['L']:.3f} event {s_['event_med']:.3f} fail {s_['fail']:2d} long {s_['long_loss']:2d} | {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}]", flush=True)
            (HERE/'results'/'contact_matrix.json').write_text(json.dumps(dict(res['matrix'], rows=rows), indent=0, allow_nan=True), encoding='utf-8')
        if phase in ('all', 'ablate'):
            rows = list(ex.map(ablate_job, [(n, c, t, i) for n in ABLATE for c, t in S.CELLS for i in range(N_ABLATE)], chunksize=1))
            ref = {key(r): r['metrics'] for r in rows if r['name'] == 'Full'}
            res['ablate'] = {}
            for n in ABLATE:
                rr = [r for r in rows if r['name'] == n]; s = stats(rr)
                d = np.array([L(r['metrics'])-L(ref[key(r)]) for r in rr]); c = boot(d)
                kinds = {}
                for r in rr:
                    if r['cause']:
                        kinds[r['cause']['kind']] = kinds.get(r['cause']['kind'], 0)+1
                g = {}
                for gn, f in GROUPS:
                    q = [r for r in rr if f(r)]; dg = np.array([L(r['metrics'])-L(ref[key(r)]) for r in q]); cg = boot(dg)
                    g[gn] = [float(dg.mean()), float(cg[0]), float(cg[1])]
                res['ablate'][n] = dict(s, vs_full=[float(d.mean()), float(c[0]), float(c[1])], causes=kinds, groups=g)
                print(f"{n:26s} L {s['L']:.3f} rmse {s['rmse_med']:.3f}/{s['rmse_p90']:.2f} event {s['event_med']:.3f} fail {s['fail']:2d} long {s['long_loss']:2d} | vs Full {d.mean():+.3f} [{c[0]:+.3f},{c[1]:+.3f}] | {kinds}")
                print('    '+' '.join(f"{gn}:{v[0]:+.2f}[{v[1]:+.2f},{v[2]:+.2f}]" for gn, v in g.items()), flush=True)
            (HERE/'results'/'ablation.json').write_text(json.dumps(dict(res['ablate'], rows=rows), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
