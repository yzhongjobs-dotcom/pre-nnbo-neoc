"""Premise test for redesigning the active excitation: is stiffness information energy-limited, and what does the ripple cost?

Fixed excitation amplitude a_force x {0.5, 1, 2, 3}, ungated Huber-D4 (own mechanism), same controller, same cases. Nothing else changes.
Predictions fixed before running, from the window statistics (sd of ln k ~ sigma / (k A_F sqrt(info))):
 P1  plateau raw estimate spread scales as 1/A (x2 amplitude -> half).
 P2  transition error falls with A: noise and disturbance leakage both scale as 1/A; in-window k change does not.
 P3  the probe ripple in the force grows in proportion to A, and in plateau-dominated scenes (A, B) L rises with A.
 P4  where transitions or events dominate (A_large, cough) L falls with A, at least from 0.5 to 2.
If P1 fails the information is not energy-limited and the whole redesign is unfounded.
Writes results/excitation_energy_check.json.
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
import memory_segment_eval as SE
import overlay_compare as O

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent
SCALES = (.5, 1., 2., 3.)
GROUPS = ('A', 'B', 'A_large', 'cough')


def job(task):
    scale, case = task
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    pol.params['a_force'] *= scale
    raw = {}
    inner = pol.estimator.estimate

    def rec(o, n, p):
        out = inner(o, n, p); raw[n-1] = out[0]
        return out

    pol.estimator.estimate = rec
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    _, masks = SE.segment_masks(log, dict(event_info=sc['event_info']))
    k = np.array([raw.get(j, np.nan) for j in range(len(log))])
    e = np.log(np.where(k > 0, k, np.nan))-np.log(np.where(log[:, 26] > 0, log[:, 26], np.nan))
    seg = {s: e[masks[s]] for s in ('plateau', 'transition', 'event')}
    return dict(scale=scale, index=case['index'], group=case['group'], L=float(np.log1p(m['J'])), rmse=float(m['rmse']), ripple=float(m['ripple']), osc=float(m['osc']), peak=float(m['peak']),
                seg={s: v[np.isfinite(v)].tolist() for s, v in seg.items()})


def rsd(x):
    q = np.percentile(x, [25, 75]); return float((q[1]-q[0])/1.349)


def main():
    cases = [c for c in CC.cases() if c['group'] in GROUPS]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(s, c) for s in SCALES for c in cases], chunksize=1))
    out = {}
    print('group scale | L mean | force rmse med | ripple rms med | peak max | plateau: within-case sd med, |e| med | transition |e| med/P90 | event |e| med/P90')
    for g in GROUPS:
        for s in SCALES:
            rr = [r for r in rows if r['group'] == g and r['scale'] == s]
            pl = [np.array(r['seg']['plateau']) for r in rr]
            within = [rsd(x) for x in pl if len(x) > 20]
            cat = lambda key: np.concatenate([np.array(r['seg'][key]) for r in rr]) if rr else np.array([])
            tr, ev, pa = cat('transition'), cat('event'), cat('plateau')
            f = lambda x: (float(np.median(np.abs(x))), float(np.percentile(np.abs(x), 90))) if len(x) > 20 else (float('nan'), float('nan'))
            row = dict(L=float(np.mean([r['L'] for r in rr])), rmse=float(np.median([r['rmse'] for r in rr])), ripple=float(np.median([r['ripple'] for r in rr])), peak=float(max(r['peak'] for r in rr)),
                       plateau_within=float(np.median(within)) if within else float('nan'), plateau=f(pa), transition=f(tr), event=f(ev))
            out[f'{g}|{s}'] = row
            print(f"{g:8s} x{s:<3} | {row['L']:.3f} | {row['rmse']:.3f} | {row['ripple']:.3f} | {row['peak']:5.1f} | {row['plateau_within']:.3f} {row['plateau'][0]:.3f} | {row['transition'][0]:.3f}/{row['transition'][1]:.3f} | {row['event'][0]:.3f}/{row['event'][1]:.3f}")
    (HERE/'results/excitation_energy_check.json').write_text(json.dumps(dict(summary=out, rows=[{k: v for k, v in r.items() if k != 'seg'} for r in rows]), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
