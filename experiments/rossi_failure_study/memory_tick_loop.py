"""Closed loop for the tick-level Kalman: same controller and own-mechanism definition as the ungated arms; 42 paired cases; no retuning.

Predictions fixed before running: A and A+B near existing methods, since plateau force error is insensitive to k accuracy; A_large clearly better than window + IMM,
which adds lag; B still limited by calibration; combo not robust, so no conclusion from it.
The harness is checked first: with the unchanged JointKalman it must return the stored Kalman row exactly.
Writes results/memory_tick_loop.json.
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

import corrected_comparison as CC
import id_baselines as I
import memory_segment_eval as SE
import memory_tick_kalman as T
import overlay_compare as O
from estimators import _fit_record

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent
GROUPS = ('A', 'A_slow_thin', 'B', 'A+B', 'cough', 'combo', 'A_large')
ARMS = dict(T2=T.LADDER['T2_lo'], K6=T.LADDER['K6_rel_q'])


class TickEstimator(I.EquationEstimator):
    """Reads each causal interval once through the same pair() as JointKalman."""
    def __init__(self, sigma_scale=.796, chi_scale=.551, **cfg):
        self.cfg = cfg; self.sigma_scale = sigma_scale; self.chi_scale = chi_scale; self.d_order = 2

    def _restart(self):
        super()._restart()
        last = getattr(self, 'filter', None)
        k0 = last.x[0][0] if last is not None and last.started else 30.
        self.filter = T.TickKalman(**self.cfg); self.filter.reset(k0)
        self.y_prev = None

    def estimate(self, o, n, policy):
        pair = self.pair(o, n, policy)
        if pair is None:
            return np.nan, np.inf, np.inf
        y, z = pair
        dy = 0. if self.y_prev is None else y-self.y_prev
        self.y_prev = y
        k = self.filter.step(dy, z)
        f = self.filter
        var = sum(f.w[j]*(f.P[j][0, 0]+(f.x[j][0]-k/T.UNIT)**2) for j in range(f.modes))
        sig = self.sigma_scale*math.sqrt(max(var, 0.))*T.UNIT
        self.fit = _fit_record(k, sig, 1., degree=0, tw=policy.Tc)
        return k, sig, 1.


def job(task):
    arm, case = task
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    cfg = O.chosen()['Kalman']
    pol = B.build('Kalman', cfg, sc)
    if arm != 'Kalman':
        pol.estimator = I.ContactAdapter(TickEstimator(**ARMS[arm]))
    CC.own_mechanism(pol)
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    return dict(arm=arm, index=case['index'], group=case['group'], L=float(np.log1p(m['J'])), peak=float(m['peak']), rmse=float(m['rmse']), call_p50=float(m['call_p50']), call_p99=float(m['call_p99']))


def main():
    cases = CC.cases()
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, c) for a in ('Kalman', *ARMS) for c in cases], chunksize=1))
    worst = max(abs(r['L']-stored[('Kalman', r['index'])]['L']) for r in rows if r['arm'] == 'Kalman')
    print('harness vs stored Kalman, max |dL|', worst, flush=True)
    assert worst < 1e-9, worst
    get = lambda m: {i: stored[(m, i)] for i in range(42)}
    table = {'D4 window': get('D4_nogate'), 'D4 + IMM': get('D4_imm'), 'Huber-D4 gated': get('Huber-D4'), 'Oracle ideal cal': get('oracle_keff_idealcal')}
    for a in ARMS:
        table[a] = {r['index']: r for r in rows if r['arm'] == a}
    print('method | L mean median | A slow B A+B cough combo A_large | RMSE med | peak | >15 | >36.5')
    out = {}
    for name, v in table.items():
        L = np.array([v[i]['L'] for i in range(42)]); peak = np.array([v[i]['peak'] for i in range(42)])
        rm = np.array([v[i].get('rmse', v[i].get('force_rmse')) for i in range(42)])
        gr = {g: float(np.mean([v[c['index']]['L'] for c in cases if c['group'] == g])) for g in GROUPS}
        out[name] = dict(mean=float(L.mean()), median=float(np.median(L)), groups=gr, rmse=float(np.median(rm)), peak=float(peak.max()), over15=float(np.mean(peak > 15)), over36=float(np.mean(peak > 36.5)))
        print(f"{name:18s} {L.mean():.3f} {np.median(L):.3f} | "+' '.join(f'{gr[g]:.3f}' for g in GROUPS)+f" | {np.median(rm):.3f} | {peak.max():.1f} | {np.mean(peak > 15)*100:.0f}% | {np.mean(peak > 36.5)*100:.0f}%")
    print('call ms p50/p99 T2', np.median([r['call_p50'] for r in rows if r['arm'] == 'T2'])/1e3, np.median([r['call_p99'] for r in rows if r['arm'] == 'T2'])/1e3)
    (HERE/'results/memory_tick_loop.json').write_text(json.dumps(dict(summary=out, rows=rows, config={k: {a: b for a, b in v.items()} for k, v in ARMS.items()}), indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
