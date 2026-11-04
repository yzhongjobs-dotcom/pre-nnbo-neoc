"""Value ceiling of explaining probe motion correctly: same controller with true gains, only the force-to-motion conversion changes.

The tracker's measurement is Z += dF/k - dy, with k the gain at the probe frequency. The plant's response to motion at other
frequencies differs (gain at 1 Hz about .6-.7 of it), so part of the probe's own motion enters Z. Variants, all privileged:
  oracle    existing true-gain reference (rows in main_rows.json)
  dynamic   dF is inverted through the plant's true Prony memory and robot series stiffness, the exact elastic inverse
  low_gain  scalar conversion with the true gain at 1 Hz instead of the probe-frequency gain
If dynamic does not beat oracle, a better contact model cannot help the motion channel in closed loop.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import physics_plant as PP
import sweeps as S
from estimator_swap import SwapPolicy
from failure_replay import Measured

HERE = Path(__file__).resolve().parent
BASE = 386000
N_SCENE = 5


class Truth(Measured):
    """Exposes the true gains at the probe frequency and at 1 Hz to the privileged variants."""
    def __call__(self, o):
        a = np.full((o.n, 16), np.nan)
        a[:, [0, 1, 2, 5, 8, 9, 11, 12, 13]] = o.log[:o.n, [0, 1, 2, 5, 8, 9, 11, 12, 13]]
        return self.inner(PP.OLD.Obs(a, o.n))


class ConversionPolicy(SwapPolicy):
    conversion = 'dynamic'; sc = None

    def reset(self, platform):
        super().reset(platform)
        sc = self.sc; alpha, tau = sc['alpha'], sc['tau']
        self.kr = float(sc['kr']); self.alpha = alpha
        self.e = np.exp(-self.Tc/tau); self.g = np.exp(-self.Tc/(2*tau))
        w = 2*math.pi*20.8
        self.r_c = 1.+float(np.sum(alpha*(w*tau)**2/(1+(w*tau)**2)))
        self.q = np.zeros(3); self.Zd = 0.; self._ctx = None
        self._wrap()

    def _wrap(self):
        if getattr(self.tracker, '_conversion', None) is self:
            return
        original = self.tracker.update

        def update(_):
            return original(self._next_z())
        self.tracker.update = update
        self.tracker._conversion = self

    def _next_z(self):
        o, n = self._ctx; j = n-1
        decay = (1-self.afc)**self.ns; denom = 1-decay
        fin = (o.log[j, 5]-decay*o.log[j-1, 5])/denom
        y0, y1 = o.log[j-1, 2], o.log[j, 2]
        y = float(y0+self._encoder_fraction*(y1-y0))
        dF = fin-self._fin_prev; dy = y-self._y_prev
        if self.conversion == 'low_gain':
            dc = dF/max(float(o.log[j, 9]), 1.)
        else:
            k8 = max(float(o.log[j, 8]), 1.)
            kc = 1./max(1./k8-1./self.kr, 1e-9); ke = kc/self.r_c
            dFe = (dF-float(np.sum((self.e-1)*self.q)))/(1+float(np.sum(self.alpha*self.g)))
            self.q = self.e*self.q+self.alpha*self.g*dFe
            dc = dFe/ke+dF/self.kr
        self.Zd += dc-dy
        self.Z = self.Zd
        return self.Zd

    def drift_pred(self, o, n):
        self._wrap()
        self._ctx = (o, n)
        out = super().drift_pred(o, n)
        if self.Z == 0. or self._fin_prev is None:
            self.Zd = 0.; self.q = np.zeros(3)
        return out


def job(a):
    conversion, c, t, i, controls, est = a
    sc = S.cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(controls); p.update(oracle_meas=True, oracle_fb=True)
    pol = ConversionPolicy(S.estimator('Projection', est), single_model=False, **p)
    pol.conversion = conversion; pol.sc = sc
    log = PP.episode(sc, Truth(pol), seed)
    return dict(variant=conversion, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, controls.get('f_exc', 20.8)))


def stats(rows):
    L = [math.log1p(min(r['metrics']['J'], 1000)) for r in rows]
    rm = np.array([r['metrics']['rmse'] for r in rows], float); ev = np.array([r['metrics']['event_rmse'] for r in rows], float)
    return dict(n=len(rows), L=float(np.mean(L)), rmse_med=float(np.nanmedian(rm)), rmse_p90=float(np.nanpercentile(rm, 90)),
                event_med=float(np.nanmedian(ev)) if np.isfinite(ev).any() else None,
                fail=sum(bool(r['metrics']['failed']) or r['metrics']['peak'] > 50 for r in rows),
                long_loss=sum(r['metrics']['long_loss'] for r in rows))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']['Projection']['best']['params']
    args = [(v, c, ti, i, controls, est) for v in ('dynamic', 'low_gain') for c, ti in S.CELLS for i in range(N_SCENE)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, args, chunksize=1))
    main_rows = json.loads((HERE/'results'/'main_rows.json').read_text(encoding='utf-8'))
    base = {m: [r for r in main_rows if r['config'] == 'C_neutral' and r['method'] == m and r['i'] < N_SCENE] for m in ('Projection', 'Oracle')}
    table = {m: stats(base[m]) for m in base}
    for v in ('dynamic', 'low_gain'):
        table[v] = stats([r for r in rows if r['variant'] == v])
    groups = {'soft': lambda r: r['tissue'] == 'soft', 'mid': lambda r: r['tissue'] == 'mid', 'hard': lambda r: r['tissue'] == 'hard',
              'event': lambda r: r['cond'] in ('flinch', 'cough'), 'fast': lambda r: r['cond'] == 'shift_fast',
              'smooth': lambda r: r['cond'] in ('breath', 'shift_slow'), 'rib+combo': lambda r: r['cond'] in ('rib', 'combo')}
    by = {}
    for g, pick in groups.items():
        by[g] = {m: stats([r for r in base[m] if pick(r)]) for m in base}
        for v in ('dynamic', 'low_gain'):
            by[g][v] = stats([r for r in rows if r['variant'] == v and pick(r)])
    key = lambda r: (r['cond'], r['tissue'], r['i'])
    O = {key(r): r['metrics'] for r in base['Oracle']}
    D = {key(r): r['metrics'] for r in rows if r['variant'] == 'dynamic'}
    d = np.array([math.log1p(min(D[k]['J'], 1000))-math.log1p(min(O[k]['J'], 1000)) for k in O])
    rng = np.random.default_rng(1); m = d[rng.integers(0, len(d), (10000, len(d)))].mean(1)
    out = dict(table=table, by=by, dynamic_minus_oracle=dict(mean=float(d.mean()), ci=[float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]), rows=rows)
    (HERE/'results'/'dynamic_gain_value.json').write_text(json.dumps(out, indent=1, allow_nan=True), encoding='utf-8')
    fmt = lambda v: f"L {v['L']:.3f} rmse {v['rmse_med']:.3f}/{v['rmse_p90']:.3f} event {v['event_med'] if v['event_med'] is None else round(v['event_med'], 3)} fail {v['fail']} long {v['long_loss']}"
    for k, v in table.items():
        print(f'{k:11s} n{v["n"]}', fmt(v))
    print('dynamic-oracle', out['dynamic_minus_oracle'])
    for g, v in by.items():
        print(g, '|', ' | '.join(f"{m[:4]} {x['L']:.3f}/{x['rmse_med']:.3f}/{x['fail']}" for m, x in v.items()))


if __name__ == '__main__':
    main()
