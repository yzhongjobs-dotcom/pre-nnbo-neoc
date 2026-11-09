"""Tick by tick: where Huber-D4 and same-config D4 first separate in a closed loop, and what each estimate was.

Usage: python -X utf8 huber_case_trace.py   (prints the cases discussed in research/huber_d4_theory_vs_experiment_20261109.md)
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True

import numpy as np

import corrected_comparison as CC
import huber_same_config_loop as L
import memory_segment_eval as SE
import overlay_compare as O

P3 = SE.P3


def traced(method, case):
    original = O.chosen
    cfg = original(); cfg['D4'] = dict(L.SAME)
    O.chosen = lambda: cfg
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    host = method.replace('_nogate', '')
    pol = O.B.build(host, cfg[host], sc)
    if method.endswith('_nogate'):
        CC.own_mechanism(pol)
    est = pol.estimator; inner = est.estimate
    rec = {}

    def wrap(o, n, p):
        out = inner(o, n, p)
        rec[n-1] = (float(out[0]), float(out[1]), float(out[2]))
        return out
    est.estimate = wrap
    good = []
    step = pol.step

    def capture(o):
        u = step(o); good.append((pol.k_good, pol.kbar)); return u
    pol.step = capture
    log = P3.episode(sc, pol, case['seed'])
    O.chosen = original
    n = len(log)
    raw = np.array([rec.get(j, (np.nan,)*3) for j in range(n)])
    return dict(log=log, raw=raw, good=np.array(good[:n]), sc=sc)


def show(index, gated):
    case = CC.cases()[index]
    a = traced('Huber-D4' if gated else 'Huber-D4_nogate', case); b = traced('D4' if gated else 'D4_nogate', case)
    la, lb = a['log'], b['log']; n = min(len(la), len(lb))
    F = 5; U = 11
    diff = np.abs(la[:n, U]-lb[:n, U]); first = int(np.argmax(diff > 1e-9))
    big = int(np.argmax(np.abs(la[:n, 18]-lb[:n, 18]) > 2.))
    print(f"\n=== case {index} {case['group']} {'gated' if gated else 'nogate'} events {[(e['kind'], round(e['start'],2), round(e['end'],2)) for e in a['sc']['event_info']]}")
    print('peak Huber', round(float(np.max(np.abs(la[:, 18]))), 1), 'D4', round(float(np.max(np.abs(lb[:, 18]))), 1), 'first command diff t', round(float(la[first, 0]), 3), 'first force gap>2N t', round(float(la[big, 0]), 3) if np.any(np.abs(la[:n, 18]-lb[:n, 18]) > 2.) else None)
    idx = sorted(set(list(range(max(first-2, 0), min(first+6, n)))+list(range(max(big-14, 0), min(big+6, n)))))
    print(' t      keff   | raw H  raw D | good H good D | bar H  bar D | F H   F D')
    for j in idx:
        print(f" {la[j,0]:.3f} {la[j,26]:7.0f} | {a['raw'][j,0]:7.0f} {b['raw'][j,0]:7.0f} | {a['good'][j,0]:7.0f} {b['good'][j,0]:7.0f} | {a['good'][j,1]:7.0f} {b['good'][j,1]:7.0f} | {la[j,18]:6.2f} {lb[j,18]:6.2f}")
    # Where the force peaks in the worse run
    w = la if np.max(np.abs(la[:, 18])) > np.max(np.abs(lb[:, 18])) else lb
    p = int(np.argmax(np.abs(w[:, 18]))); print('worse run peak at t', round(float(w[p, 0]), 3))
    # Within the same-observation replay: activation history on the Huber run
    return a, b


if __name__ == '__main__':
    for idx, gated in ((26, False), (26, True), (36, True), (35, False)):
        show(idx, gated)
