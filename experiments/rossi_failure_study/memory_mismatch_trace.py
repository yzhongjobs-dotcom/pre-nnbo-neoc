"""Root cause of the two mismatches of the tick-level Kalman (research/kalman_idea_analysis_20261109.md section 7).

Mismatch 1, cough: replay error improves, closed-loop L worsens. Closed loop, tick by tick, on cough cases 24, 26, 27:
the controller divides the force change by k_good (kappa_g = k_eff/k_good) and feedback by k_bar = k_hat + .807 sigma (kappa = k_eff/k_bar).
Stage 1 compares T2 with D4 + IMM (24 ticks) around the force peak. Stage 2 (argument 'swap') holds the estimate and replaces only sigma.
Writes results/memory_mismatch_trace.json.
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
import huber_unified_compare as U
import id_baselines as I
import memory_segment_eval as SE
import memory_tick_loop as TL
import overlay_compare as O
import rib_transition_eval as R

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent
CASES = (24, 26, 27)


def build(arm, sc, sigma_mode=None):
    if arm == 'T2':
        pol = B.build('Kalman', O.chosen()['Kalman'], sc)
        inner = TL.TickEstimator(**TL.ARMS['T2'])
        pol.estimator = I.ContactAdapter(inner)
        CC.own_mechanism(pol)
        if sigma_mode:
            # Replace only the sigma passed to the controller; the estimate is unchanged.
            est = pol.estimator; orig = est.estimate; scale = dict(x3=3., x10=10.)[sigma_mode]
            def estimate(o, n, policy):
                k, s, c = orig(o, n, policy)
                return k, s*scale if np.isfinite(s) else s, c
            est.estimate = estimate
        return pol
    U.configure(24)
    pol = B.build('D4', O.chosen()['D4'], sc)
    CC.own_mechanism(pol)
    R.make_filter(pol, 'imm')
    return pol


def trace(task):
    arm, idx, sigma_mode = task
    O.setup()
    case = CC.cases()[idx]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = build(arm, sc, sigma_mode)
    rec = []; step = pol.step

    def capture(o):
        u = step(o)
        rec.append((pol.k_good, pol.kbar, getattr(pol, 'k_hat', np.nan), getattr(pol, 'rel', np.nan)))
        return u
    pol.step = capture
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    rec = np.asarray(rec)[:len(log)]
    t = log[:, 0]; k = log[:, 26]; F = log[:, 18]
    ok = k > 200.
    kg = np.where(ok, k/rec[:, 0], np.nan); kb = np.where(ok, k/rec[:, 1], np.nan)
    rows = []
    for e in sc['event_info']:
        w = (t >= e['start']-.05) & (t <= e['end']+.3)
        j = np.flatnonzero(w)
        p = j[np.argmax(np.abs(F[j]-8.))]
        before = j[(t[j] <= t[p]) & (t[j] >= t[p]-.08)]
        rows.append(dict(start=e['start'], end=e['end'], peak_t=float(t[p]), peak_F=float(F[p]),
                         kappa_g_at_peak=float(np.nanmax(kg[before])), kappa_g_min=float(np.nanmin(kg[before])),
                         kappa_max=float(np.nanmax(kb[before])), rel_sigma=float(np.nanmedian(rec[before, 3])),
                         k_rise=float(np.log(np.nanmax(k[j])/k[j[0]])), force_rms=float(np.sqrt(np.mean((F[j]-8.)**2)))))
    return dict(arm=arm, sigma=sigma_mode, index=idx, L=float(np.log1p(m['J'])), peak=float(m['peak']), events=rows,
                series=dict(t=t[::1].tolist(), k=k.tolist(), F=F.tolist(), kg=rec[:, 0].tolist(), kb=rec[:, 1].tolist()) if idx == 26 else None)


def main():
    swap = len(sys.argv) > 1 and sys.argv[1] == 'swap'
    tasks = [(a, i, None) for a in ('T2', 'D4_imm') for i in CASES]
    if swap:
        tasks = [('T2', i, s) for i in CASES for s in ('x3', 'x10')]
    with ProcessPoolExecutor(6) as ex:
        res = list(ex.map(trace, tasks, chunksize=1))
    print('arm sigma case | L peak | per event: peak F, t; max k_eff/k_good and min in .08 s before; max k_eff/k_bar; median sigma/k; ln k rise in event')
    for r in res:
        print(f"{r['arm']:6s} {str(r['sigma']):4s} {r['index']} | {r['L']:.3f} {r['peak']:5.1f}")
        for e in r['events']:
            print(f"     F {e['peak_F']:5.1f} t {e['peak_t']:.3f} | kg {e['kappa_g_at_peak']:.2f} min {e['kappa_g_min']:.2f} | kb {e['kappa_max']:.2f} | rel {e['rel_sigma']:.3f} | rise {e['k_rise']:.2f} | Frms {e['force_rms']:.2f}")
    name = 'memory_mismatch_trace_swap.json' if swap else 'memory_mismatch_trace.json'
    (HERE/'results'/name).write_text(json.dumps(res), encoding='utf-8')


if __name__ == '__main__':
    main()
