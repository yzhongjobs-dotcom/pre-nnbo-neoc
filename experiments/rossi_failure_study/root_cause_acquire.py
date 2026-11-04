"""Deployable test of the acquisition root cause: stop growing the probe once a window is precise, and take the first gain on precision alone.

The chi test rejected first windows whose k was already within a few percent, while the probe kept growing. Also checks the
implementation at a large but well conditioned rho: the relaxed estimator must then match the cubic projection on the same arrays.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import sweeps as S
import common
import physics_plant as PP
from dynamic_gain_value import Truth
from failure_replay import cell_scene
from ridge_estimator import RidgeProjection
from ridge_final import BASE as FINAL_BASE
from root_cause import Full, Plain, estimator
import root_cause_confirm as C

HERE = Path(__file__).resolve().parent
CTL = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))['public']['C_neutral']['best']['params']


class Acquire:
    def _stiffness_est(self, o, n):
        p = self.params
        acquiring = not np.isfinite(self.k_good)
        saved = p['chi2_max'], p['meas_bound']
        if acquiring:
            p['chi2_max'] = np.inf; p['meas_bound'] = None
        out = super()._stiffness_est(o, n)
        p['chi2_max'], p['meas_bound'] = saved
        if acquiring and not np.isfinite(self.k_good) and np.isfinite(getattr(self, 'rel', np.nan)) and self.rel <= p['rel_max']:
            self.k_seek = self.k_seek/.85     # undoes this tick's growth, applied after the estimate
        return out


class PA(Acquire, Plain):
    pass


class FA(Acquire, Full):
    pass


def job(a):
    name, c, t, i = a
    sc = cell_scene(c, t, i, FINAL_BASE); seed = FINAL_BASE+sum(ord(x) for x in c+t)*100+i
    p = dict(S.FINAL_C); p.update(CTL)
    if name.startswith('P5_16'):
        pol = PA(estimator('P5_16'), single_model=False, **p)
    else:
        pol = FA(RidgeProjection(**C.FULL), single_model=False, **p); pol.sc = sc; pol.centre = True
    pol.intervene = ()
    log = PP.episode(sc, Truth(pol), seed)
    return dict(name=name, cond=c, tissue=t, i=i, metrics=S._metric(log, pol, sc, 20.8), acq=bool(np.isfinite(pol.k_good)))


def implementation():
    sc = cell_scene('breath', 'mid', 0, FINAL_BASE)
    p = dict(S.FINAL_C); p.update(CTL)
    pol = Plain(estimator('P3_16'), single_model=False, **p); log = PP.episode(sc, Truth(pol), FINAL_BASE)
    a = pol.params; D = np.diff(np.eye(16), 4, axis=0)
    out = {}
    for rho in (1e4, 1e6, 1e8):
        M = np.eye(16)-np.linalg.inv(np.eye(16)+rho*D.T@D)
        rk, rr = [], []
        for n in range(40, len(log), 3):
            X, Y, Z, keep = common._build(log, n-17, 16, 3, pol.ns, pol.eb, pol.afc, pol.ndel, log[n-1, 0], log[n-1, 2], 16*pol.ns*1e-3, a['enc_y'])
            if not keep.all():
                continue
            Q, _ = np.linalg.qr(X); P = np.eye(16)-Q@Q.T
            kp = (Y@P@Z)/(Y@P@Y); km = (Y@M@Z)/(Y@M@Y)
            r = Z-kp*Y
            rk.append(abs(km/kp-1)); rr.append((r@M@r)/(r@P@r))
        out[rho] = dict(k_rel_p50=float(np.median(rk)), k_rel_p99=float(np.percentile(rk, 99)), residual_ratio_p50=float(np.median(rr)), cond=float(np.linalg.cond(np.eye(16)+rho*D.T@D)))
        print('rho', rho, out[rho], flush=True)
    return out


def main():
    res = dict(implementation=implementation())
    with ProcessPoolExecutor(15) as ex:
        new = list(ex.map(job, [(n, c, ti, i) for n in ('P5_16 acquire', 'Full acquire') for c, ti in S.CELLS for i in range(10)], chunksize=1))
        ref = list(ex.map(C.full_cell, [('P5_16', c, ti, i, CTL, C.FULL, False) for c, ti in S.CELLS for i in range(10)], chunksize=1))
    res['summary'] = C.report(ref+new, ['P5_16 acquire', 'Full acquire'], 'P5_16')
    print('never acquired', {n: sum(not r['acq'] for r in new if r['name'] == n) for n in ('P5_16 acquire', 'Full acquire')})
    (HERE/'results'/f'root_cause_acquire_{FINAL_BASE}.json').write_text(json.dumps(dict(res, rows=new), indent=0, allow_nan=True), encoding='utf-8')


if __name__ == '__main__':
    main()
