"""Where is the 20.8 Hz probe amplitude lost on falls? Closed loop, Huber-D4 ungated, BASE and EXC_TRUE (as v2_excitation_control.py).
Per tick, lock-in amplitude at 20.8 Hz over the last 24 ticks (linear trend removed) of: the excitation part of the command (running sum of dx/k_used),
the total command, the encoder, and k_eff times the encoder; by segment, relative to the plateau median of the same case.
Writes results/v2_probe_chain.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O
import rib_transition_eval as R
import v2_verify as V

B = O.B
HERE = Path(__file__).resolve().parent
CASES = (36, 37, 38, 39, 40, 41)
TT = np.arange(24)*.008
CAR = np.exp(-2j*np.pi*20.8*TT)


def lockin(x):
    out = np.full(len(x), np.nan)
    for j in range(24, len(x)):
        w = x[j-24:j]
        if not np.isfinite(w).all():
            continue
        w = w-np.polyval(np.polyfit(TT, w, 1), TT)
        out[j] = 2*abs(np.mean(w*CAR))
    return out


def job(task):
    arm, idx = task
    O.setup()
    case = CC.cases()[idx]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    sensing = pol.step; inner = sensing.step
    rec = []; acc = [0.]

    def step(o):
        x_before = pol.x_prev
        u = inner(o)
        dx = pol.x_prev-x_before
        k_un = pol.k_good if np.isfinite(pol.k_good) else pol.k_seek
        kt = float(getattr(o, 'truth_keff', np.nan))
        k_used = k_un
        if arm == 'exc_true' and kt > 0. and np.isfinite(k_un) and k_un > 0.:
            u = u+dx*(1./kt-1./k_un); k_used = kt
        acc[0] += dx/k_used if np.isfinite(k_used) and k_used > 0 else 0.
        rec.append((acc[0], u))
        return u
    sensing.step = step
    log = R.oracle_run(sc, pol, case['seed'])
    rec = np.asarray(rec)[:len(log)]
    t = log[:, 0]; k = log[:, 26]
    seg = V.segments(t, k, sc['event_info'], k > 200.)
    chans = dict(exc_cmd=lockin(rec[:, 0]), command=lockin(rec[:, 1]), encoder=lockin(log[:, 2]))
    chans['probe_force'] = k*chans['encoder']
    out = {}
    for c, v in chans.items():
        base = np.nanmedian(v[seg['plateau']]) if seg['plateau'].sum() > 20 else np.nanmedian(v[k > 200.])
        out[c] = {s: float(np.nanmedian(v[seg[s]]/base)) if seg[s].any() else None for s in ('plateau', 'rising', 'falling', 'reversal')}
    return dict(arm=arm, index=idx, chans=out)


def main():
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, i) for a in ('base', 'exc_true') for i in CASES], chunksize=1))
    rep = {}
    for a in ('base', 'exc_true'):
        rep[a] = {}
        print('---', a, '(median over cases of segment median / plateau median)')
        for c in ('exc_cmd', 'command', 'encoder', 'probe_force'):
            rep[a][c] = {s: float(np.nanmedian([r['chans'][c][s] for r in rows if r['arm'] == a and r['chans'][c][s] is not None])) for s in ('rising', 'falling', 'reversal')}
            print(f"  {c:12s}", {s: round(v, 3) for s, v in rep[a][c].items()})
    (HERE/'results/v2_probe_chain.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
