"""Which part of the band-pass estimate causes its hard-tissue runaways: causal interventions on one signal at a time.

Same scenes, noise, control and estimator parameters as main_rows.json (base 386000, C_neutral, hard tissue).
Each variant replaces one signal with a privileged value and keeps everything else from the band-pass estimator:
  base        band-pass unchanged; must reproduce main_rows.json
  true_good   k_good = true gain each tick        -> is the measurement gain the cause?
  true_bar    feedback stiffness = true gain      -> is the feedback stiffness the cause?
  lag_good    k_good = true gain delayed .25 s    -> is the lag alone enough?
  stale_good  k_good = true gain, sampled only when the band-pass gate accepts -> is the update timing enough?
  open_gate   k_good = band-pass estimate, chi gate and bound removed -> does the gate cause it?
  memory_good k_good = true gain through the band-pass delay and EWRLS memory, written at its acceptances -> is the
              value error explained by its memory alone, or by bias from interference?
  margin_x    each band-pass k_good write multiplied by exp(x); x near the projection's +.12 median offset
              -> does the projection's safe-side offset explain the gap?
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import signal

import sweeps as S
from estimator_swap import SwapPolicy
from failure_replay import OracleMeasured

HERE = Path(__file__).resolve().parent
BASE = 386000
LAG_TICKS = int(round(.25/.008))
HARD = [(c, 'hard') for c in ('rest', 'breath', 'shift_slow', 'shift_fast', 'flinch', 'cough', 'tremor')]
VARIANTS = ('base', 'true_good', 'true_bar', 'lag_good', 'stale_good', 'open_gate', 'memory_good', 'margin_0.06', 'margin_0.12', 'margin_0.2')


class Intervened(SwapPolicy):
    """SwapPolicy with one signal replaced; the first-acceptance time stays the estimator's own."""
    intervene = 'base'

    def reset(self, platform):
        super().reset(platform)
        e = self.estimator; b, a = signal.sos2tf(e.sos)
        f0 = float(self.params.get('f_exc', 20.8))
        self.delay_ticks = int(round(float(signal.group_delay((b, a), w=[f0], fs=1./self.Tc)[1][0])))
        self.k_avg = np.nan

    def _stiffness_est(self, o, n):
        before = self.k_good
        out = super()._stiffness_est(o, n)
        if self.intervene in ('true_good', 'lag_good') and np.isfinite(self.k_good):
            j = n-1 if self.intervene == 'true_good' else max(n-1-LAG_TICKS, 0)
            self.k_good = max(float(o.log[j, 8]), 1.)
        if self.intervene.startswith('margin_') and np.isfinite(self.k_good) and self.k_good != before:
            self.k_good *= math.exp(float(self.intervene.split('_')[1]))
        if self.intervene == 'memory_good':
            kj = float(o.log[max(n-1-self.delay_ticks, 0), 8])
            if kj > 0:
                lam = self.estimator.lam
                self.k_avg = kj if not np.isfinite(self.k_avg) else lam*self.k_avg+(1-lam)*kj
            if np.isfinite(self.k_good) and self.k_good != before and np.isfinite(self.k_avg):
                self.k_good = max(self.k_avg, 1.)
            elif np.isfinite(before):
                self.k_good = before
        return out


def make(variant, controls, est):
    p = dict(S.FINAL_C); p.update(controls)
    if variant == 'true_bar':
        p['oracle_fb'] = True
    if variant == 'stale_good':
        p['oracle_meas'] = 'stale'
    if variant == 'open_gate':
        p['kgood_rule'] = 'precise'; p['meas_bound'] = None
    pol = Intervened(S.estimator('BandPass', est), single_model=False, **p)
    pol.intervene = variant
    return pol


def job(a):
    variant, c, t, i, controls, est = a
    sc = S.cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    pol = make(variant, controls, est)
    log = S.PP.episode(sc, OracleMeasured(pol), seed)
    return dict(variant=variant, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, controls.get('f_exc', 20.8)))


def summary(rows):
    fail = sum(bool(r['metrics']['failed']) or r['metrics']['peak'] > 50 for r in rows)
    rmse = np.array([r['metrics']['rmse'] for r in rows], float)
    level = np.mean([math.log1p(min(r['metrics']['J'], 1000)) for r in rows])
    return dict(n=len(rows), L=float(level), rmse_med=float(np.nanmedian(rmse)), rmse_p90=float(np.nanpercentile(rmse, 90)),
                fail=fail, long_loss=sum(r['metrics']['long_loss'] for r in rows))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']
    est = t['estimators']['C_neutral']['BandPass']['best']['params']
    args = [(v, c, ti, i, controls, est) for v in VARIANTS for c, ti in HARD for i in range(15)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, args, chunksize=1))
    main_rows = json.loads((HERE/'results'/'main_rows.json').read_text(encoding='utf-8'))
    ref = {m: [r for r in main_rows if r['config'] == 'C_neutral' and r['method'] == m and r['tissue'] == 'hard']
           for m in ('Projection', 'BandPass', 'Oracle')}
    table = {m: summary(ref[m]) for m in ref}
    table.update({v: summary([r for r in rows if r['variant'] == v]) for v in VARIANTS})
    by_key = {(r['cond'], r['i']): r['metrics']['J'] for r in ref['BandPass']}
    gap = max(abs(r['metrics']['J']-by_key[(r['cond'], r['i'])]) for r in rows if r['variant'] == 'base')
    out = dict(max_base_gap_J=float(gap), table=table, rows=rows)
    (HERE/'results'/'kgood_intervention.json').write_text(json.dumps(out, indent=1, allow_nan=True), encoding='utf-8')
    print('base reproduces main rows, max J gap', gap)
    for k, v in table.items():
        print(f"{k:11s} n{v['n']} L {v['L']:.3f} rmse {v['rmse_med']:.3f}/{v['rmse_p90']:.3f} fail {v['fail']} long {v['long_loss']}")


if __name__ == '__main__':
    sys.exit(main())
