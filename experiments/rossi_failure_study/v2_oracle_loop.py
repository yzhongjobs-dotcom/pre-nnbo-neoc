"""Privileged closed-loop ceiling: does V2's slope have control value if the run start and slope scale were known?

Diagnostic only, not a method. Arms share one implementation and controller (Kalman host, own mechanism, as the tick filter):
  C24      constant k, last 24 rows, D4 (W = I): the parity baseline for this harness
  ORACLE   k and slope over rows since the true start of the current monotone run (24 to 96), ridge sv2 = v_true^2;
           run start and v_true read from the recorded Huber-D4 closed loop of the same scene (the gain is set by scan position,
           so its timing barely changes between loops). Everything else causal.
Predictions fixed before running: A_large L visibly lower; A and A_slow_thin about unchanged; cough unchanged or slightly worse (no Huber weights).
Writes results/v2_oracle_loop.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import corrected_comparison as CC
import huber_prior_check as H
import id_baselines as I
import memory_segment_eval as SE
import overlay_compare as O
import v2_history_ceiling as C
from estimators import _fit_record

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent
TC, SIG = .008, H.SIGMA
CASES = (0, 1, 2, 6, 7, 8, 24, 25, 26, 36, 37, 38, 39, 40, 41)


class RowEstimator(I.EquationEstimator):
    def __init__(self, mode, src=None, sigma_scale=.796, chi_scale=.551):
        self.mode = mode; self.src = src; self.sigma_scale = sigma_scale; self.chi_scale = chi_scale; self.d_order = 2

    def _restart(self):
        super()._restart(); self.ys = []; self.zs = []; self.ts = []

    def estimate(self, o, n, policy):
        pair = self.pair(o, n, policy)
        if pair is None:
            return np.nan, np.inf, np.inf
        y, z = pair
        self.ys.append(y); self.zs.append(z); self.ts.append(float(o.log[n-2, 0]))
        if len(self.ys) > 120:
            self.ys.pop(0); self.zs.pop(0); self.ts.pop(0)
        L = 24
        if self.mode == 'oracle':
            tnow = float(o.log[n-1, 0])
            start, v = self.src(tnow)
            L = int(np.clip(round((tnow-start)/TC), 24, 96))
        if len(self.ys) < L:
            return np.nan, np.inf, np.inf
        yy = np.asarray(self.ys[-L:]); zz = np.asarray(self.zs[-L:]); tt = np.asarray(self.ts[-L:])
        h0 = yy-yy[-1]; Mx = C.mat(L); mh0 = Mx@h0
        a = h0@mh0
        if self.mode == 'c24':
            q = mh0/a
        else:
            h1 = (tt+.0035-float(o.log[n-1, 0]))*h0; mh1 = Mx@h1
            b, c1 = h0@mh1, h1@mh1
            cp = c1+SIG**2/max(v*v, 1e-9)
            q = (mh0-(b/cp)*mh1)/(a-b*b/cp)
        k = float(q@zz)
        sig = self.sigma_scale*SIG*math.sqrt(float(q@q))
        self.fit = _fit_record(k, sig, 1., degree=0, tw=policy.Tc)
        return k, sig, 1.


def source(index):
    """Run start and slope (N/m/s) at time t from the recorded Huber-D4 closed loop of this scene."""
    data = H.load(H.WIN/f'case_{index:02d}.npz'); log, _ = H.masks_for(data)
    t = log[:, 0]; k = log[:, 26]
    s = uniform_filter1d(np.log(np.maximum(k, 1.)), 25, mode='nearest'); r = np.gradient(s, t)
    sign = np.where(r > .2, 1, np.where(r < -.2, -1, 0))
    start = np.zeros(len(t))
    cur = t[0]
    for j in range(1, len(t)):
        if sign[j] != sign[j-1]:
            cur = t[j]
        start[j] = cur

    def at(tq):
        j = int(np.clip(np.searchsorted(t, tq), 0, len(t)-1))
        st = start[j]; m = (t >= max(st, tq-96*TC)) & (t <= tq) & (k > 200.)
        if m.sum() < 5:
            return st, 0.
        return st, float(np.polyfit(t[m], k[m], 1)[0])
    return at


def job(task):
    arm, idx = task
    O.setup()
    case = CC.cases()[idx]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Kalman', O.chosen()['Kalman'], sc)
    pol.estimator = I.ContactAdapter(RowEstimator(arm, source(idx) if arm == 'oracle' else None))
    CC.own_mechanism(pol)
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    return dict(arm=arm, index=idx, group=case['group'], L=float(np.log1p(m['J'])), peak=float(m['peak']), rmse=float(m['rmse']))


def main():
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, i) for a in ('c24', 'oracle') for i in CASES], chunksize=1))
    uni = {(r['method'], r['index']): r for r in json.loads((HERE/'results/huber_unified_compare.json').read_text(encoding='utf-8'))['closed']}
    get = {(r['arm'], r['index']): r for r in rows}
    print('case group | C24 (this harness) L peak | ORACLE slope L peak | Huber-D4 nogate 24 L | D4 nogate 24 L')
    for i in CASES:
        a, o = get[('c24', i)], get[('oracle', i)]
        print(f"{i:2d} {a['group']:12s} | {a['L']:.3f} {a['peak']:5.1f} | {o['L']:.3f} {o['peak']:5.1f} | {uni[('Huber-D4_nogate', i)]['L']:.3f} | {uni[('D4_nogate', i)]['L']:.3f}")
    for g in ('A', 'A_slow_thin', 'cough', 'A_large'):
        cs = [i for i in CASES if get[('c24', i)]['group'] == g]
        print(f"{g:12s} mean L: C24 {np.mean([get[('c24', i)]['L'] for i in cs]):.3f}  ORACLE {np.mean([get[('oracle', i)]['L'] for i in cs]):.3f}  Huber {np.mean([uni[('Huber-D4_nogate', i)]['L'] for i in cs]):.3f}")
    (HERE/'results/v2_oracle_loop.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
