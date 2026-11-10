"""Controlled-variable test: how much of the falling-segment information deficit and error comes from excitation scaled by the lagging k_good?

Closed loop, Huber-D4 ungated (CC.own_mechanism), unchanged except one factor. The command is sent + du/k_bar + (dx - dq)/k_un with dx the excitation
increment and dq the force-based compensation, both divided by k_un = k_good. Arm EXC_TRUE divides only dx by the true k_eff of the current tick
(privileged, diagnostic), so the probe force is held at its design value; feed-forward and feedback still use k_good and k_bar.
Per tick: the estimator's own sd sigma divided by the true k_eff (sigma ∝ 1/sqrt(window information), the information measure), its raw error ln(k_hat/k_eff),
and the segments of v2_verify.segments on each run's own k_eff.
Predictions fixed before running: on falls, sigma/k of EXC_TRUE is about .78 of BASE (the measured probe-force deficit) and on rises about 1.12;
falling |error| falls only by the noise share (bias from lag unchanged); closed-loop L changes little.
Writes results/v2_excitation_control.json.
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
import rib_transition_eval as R
import v2_verify as V

B, M = O.B, SE.M
HERE = Path(__file__).resolve().parent
CASES = (0, 2, 4, 36, 37, 38, 39, 40, 41)


def job(task):
    arm, idx = task
    O.setup()
    case = CC.cases()[idx]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    sensing = pol.step
    inner = sensing.step
    rec = []

    def step(o):
        x_before = pol.x_prev
        u = inner(o)
        dx = pol.x_prev-x_before
        k_un = pol.k_good if np.isfinite(pol.k_good) else pol.k_seek
        kt = float(getattr(o, 'truth_keff', np.nan))
        if arm == 'exc_true' and kt > 0. and np.isfinite(k_un) and k_un > 0.:
            u = u+dx*(1./kt-1./k_un)
        rec.append((float(o.log[o.n-1, 0]), getattr(pol, 'k_hat', np.nan), getattr(pol, 'rel', np.nan), kt, float(pol.exc_gain), float(k_un)))
        return u
    sensing.step = step
    log = R.oracle_run(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    rec = np.asarray(rec)[:len(log)]
    t = log[:, 0]; k = log[:, 26]
    seg = V.segments(t, k, sc['event_info'], k > 200.)
    err = np.log(np.where(rec[:, 1] > 0, rec[:, 1], np.nan)/np.where(k > 200., k, np.nan))
    # Probe force actually applied: k_eff times the 20.8 Hz encoder amplitude over the last 24 ticks (linear trend removed).
    y = log[:, 2]; pf = np.full(len(t), np.nan); tt = np.arange(24)*.008
    for j in range(24, len(t)):
        w = y[j-24:j]; w = w-np.polyval(np.polyfit(tt, w, 1), tt)
        pf[j] = k[j]*2*abs(np.mean(w*np.exp(-2j*np.pi*20.8*tt)))
    rel_true = rec[:, 2]*rec[:, 1]/np.where(k > 200., k, np.nan)      # sigma / k_eff
    out = dict(arm=arm, index=idx, group=case['group'], L=float(np.log1p(m['J'])), peak=float(m['peak']), rmse=float(m['rmse']), seg={})
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        mm = seg[s] & np.isfinite(err) & np.isfinite(rec[:, 2])
        out['seg'][s] = dict(rel=rel_true[mm].tolist(), err=err[mm].tolist(), pf=pf[mm].tolist(), gain=rec[mm, 4].tolist(), kun=(rec[mm, 5]/k[mm]).tolist())
    return out


def main():
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, i) for a in ('base', 'exc_true') for i in CASES], chunksize=1))
    reg = max(abs(r['L']-stored[('Huber-D4_nogate', r['index'])]['L']) for r in rows if r['arm'] == 'base')
    print('harness: base vs stored Huber-D4 nogate, max |dL|', reg)
    rep = dict(regression=reg, seg={}, cases=[{k: v for k, v in r.items() if k != 'seg'} for r in rows])
    print('segment | sigma/k median base, exc_true, ratio | |e| median, rms, bias base -> exc_true')
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        cell = {}
        for a in ('base', 'exc_true'):
            rel = np.concatenate([r['seg'][s]['rel'] for r in rows if r['arm'] == a]); e = np.concatenate([r['seg'][s]['err'] for r in rows if r['arm'] == a])
            cell[a] = dict(n=len(e), rel=float(np.median(rel)), med=float(np.median(np.abs(e))), rms=float(np.sqrt(np.mean(e*e))), bias=float(np.median(e)))
            cell[a]['pf'] = float(np.nanmedian(np.concatenate([r['seg'][s]['pf'] for r in rows if r['arm'] == a])))
            cell[a]['exc_gain'] = float(np.nanmedian(np.concatenate([r['seg'][s]['gain'] for r in rows if r['arm'] == a])))
            cell[a]['kun_over_k'] = float(np.nanmedian(np.concatenate([r['seg'][s]['kun'] for r in rows if r['arm'] == a])))
        cell['rel_ratio'] = cell['exc_true']['rel']/cell['base']['rel']
        rep['seg'][s] = cell
        b, x = cell['base'], cell['exc_true']
        print(f"  {s:9s} n {b['n']:5d}/{x['n']:5d} | exc_gain {b['exc_gain']:.2f} {x['exc_gain']:.2f} k_un/k {b['kun_over_k']:.2f} | probe force {b['pf']:.4f} {x['pf']:.4f} | {b['rel']:.3f} {x['rel']:.3f} ratio {cell['rel_ratio']:.2f} | {b['med']:.3f}/{b['rms']:.3f}/{b['bias']:+.3f} -> {x['med']:.3f}/{x['rms']:.3f}/{x['bias']:+.3f}")
    print('closed loop L (peak): base -> exc_true')
    for i in CASES:
        b = next(r for r in rows if r['arm'] == 'base' and r['index'] == i); x = next(r for r in rows if r['arm'] == 'exc_true' and r['index'] == i)
        print(f"  {i:2d} {b['group']:8s} {b['L']:.3f} ({b['peak']:.1f}) -> {x['L']:.3f} ({x['peak']:.1f})")
    (HERE/'results/v2_excitation_control.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
