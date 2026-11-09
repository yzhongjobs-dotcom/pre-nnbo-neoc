"""Does log[:,26] equal the actual probe-frequency gain along real trajectories?

Replays the Huber-D4 closed-loop command sequence noise-free twice, with and
without a small 20.8 Hz command perturbation. On the 0.5 ms log, the difference
isolates the incremental response, so gain = |dF(w0)|/|dy(w0)| per 0.5 s window,
with motion, scan, nonlinearity and contact state all present. Plant unchanged.

python -X utf8 experiments/rossi_failure_study/keff_verify.py
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import overlay_compare as O
import memory_segment_eval as SE

P3 = SE.P3
HERE = Path(__file__).resolve().parent
W0 = 2*np.pi*20.8


def phasor(t, v):
    s = t-t.mean()
    A = np.c_[np.cos(W0*t), np.sin(W0*t), np.ones(len(t)), s]
    b = np.linalg.lstsq(A, v, rcond=None)[0]
    return complex(b[0], -b[1])


def job(case):
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = O.B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    log = P3.episode(sc, pol, case['seed'])
    commands = log[:, 11].copy()
    tick = np.arange(len(commands))*.008
    keff = float(np.median(log[log[:, 0] >= 1., 26]))
    amplitude = min(2e-4, .12/max(keff, 1.))
    platform = dict(P3.PLATFORM, log_dt=.0005)
    if 'tau_a_true' in sc:
        platform['tau_a'] = sc['tau_a_true']
    a = P3.replay(sc, commands, sensor=False, platform=platform)
    b = P3.replay(sc, commands+amplitude*np.sin(W0*tick), sensor=False, platform=platform)
    if a is None or b is None:
        return dict(index=case['index'], group=case['group'], failed=True)
    t = a[:, 0]
    rows = []
    for start in np.arange(1., t[-1]-.5, .5):
        m = (t >= start) & (t < start+.5)
        contact = (a[m, 3] > 0.).all() and (b[m, 3] > 0.).all()
        if not contact:
            continue
        dy = phasor(t[m], b[m, 2]-a[m, 2])
        dF = phasor(t[m], b[m, 18]-a[m, 18])
        rows.append(dict(t=float(start), measured=abs(dF/dy), measured_phase=float(np.angle(dF/dy)),
                         col26=float(np.median(a[m, 26])), col27=float(np.median(a[m, 27])), col21=float(np.median(a[m, 21])),
                         force=float(np.median(a[m, 18])), speed=float(np.median(a[m, 16])), theta=float(np.median(a[m, 17]))))
    ratio = np.array([r['measured']/r['col26'] for r in rows])
    return dict(index=case['index'], group=case['group'], cell=case['cell'], amplitude=amplitude, windows=rows,
                measured_over_col26=dict(median=float(np.median(ratio)), p10=float(np.percentile(ratio, 10)), p90=float(np.percentile(ratio, 90))),
                keff_median=keff, mu=float(sc.get('mu', np.nan)), tool_mass=float(sc['tool_mass']))


if __name__ == '__main__':
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, CC.cases(), chunksize=1))
    (HERE/'results/keff_verify.json').write_text(json.dumps(dict(
        note='measured = |dF/dy| at 20.8 Hz on the 0.5 ms log, from a paired noise-free replay of Huber-D4 commands with a small command perturbation.', rows=rows), default=float), encoding='utf-8')
    for r in rows:
        if r.get('failed'):
            print(r['index'], r['group'], 'failed')
            continue
        q = r['measured_over_col26']
        print(r['index'], r['group'], r['cell'][2], r['cell'][4], 'keff', round(r['keff_median']), 'measured/col26', round(q['median'], 3), round(q['p10'], 3), round(q['p90'], 3), 'mu', round(r['mu'], 3), flush=True)
