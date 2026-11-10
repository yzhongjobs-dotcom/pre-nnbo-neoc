"""Where does the bias shared by short and long history windows come from? Component-removal replays of the same commands.

Base: Huber-D4 ungated closed loop (huber_error_mechanism.closed_loop). The command sequence is replayed with one more component removed at each step,
cumulatively: R1 sensor noise; R2 calibration (true mass, delay, servo lag in the sensing stage); R3 scene-event motion (push, pull, cough, flinch);
R4 muscle activation; R5 breathing and cardiac motion; then plant counterfactuals, diagnostic only and not a change of the settled setting:
R6 friction (mu_k = mu_s = k_d = 0); R7 tool inertia (tool mass and its estimate 0); R8 viscoelastic memory (alpha = 0).
For every replay, local-linear endpoint estimates over 24 and 64 ticks (D4, W = I) are scored against that replay's own k_eff.
The shared bias of a window is e24 when the long span is harmful noise-free (|e64| > |e24| + .05 in R1); a component causes it if removing it shrinks e24 there.
What remains after R8 is split into the curvature of the replay's own k inside the window (q24' cumsum((k - k_T - v tau) dY)) and the rest.
Writes results/v2_shared_bias_source.json.
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
import huber_error_mechanism as E
import v2_history_ceiling as C
import v2_verify as V

HERE = Path(__file__).resolve().parent
CASES = (0, 1, 2, 33, 34, 36, 37)
STEPS = ('noise', 'calibration', 'event_motion', 'activation', 'breath_cardiac', 'friction', 'tool_inertia', 'viscoelastic')
TC = .008


def local(Y, z, t, j, L):
    """Local-linear endpoint weights and estimate from a stored window (rows j-L..j-1), Y relative to its end."""
    h0 = Y[-L:]; zz = z[-L:]
    if not (np.isfinite(h0).all() and np.isfinite(zz).all()):
        return None, None
    r = np.arange(j-L, j)
    h1 = (t[r]+.0035-t[j-1])*h0
    M = C.mat(L); mh0, mh1 = M@h0, M@h1
    a, b, c1 = h0@mh0, h0@mh1, h1@mh1
    det = a*c1-b*b
    if det <= 0:
        return None, None
    q = (c1*mh0-b*mh1)/det
    return q, float(q@zz)


def job(index):
    case = CC.cases()[index]
    sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, raw = E.closed_loop(case)
    n = len(rows)
    W = E.W
    # Windows of 64 ticks need a longer buffer than the stored 32: rebuild with NW = 64 for every replay.
    W.NW = 64
    mods = []
    s = dict(sc, noise_scale=0.); mods.append(dict(s))
    CC.perfect_calibration(s); mods.append(dict(s))
    comp = np.sum(sc['motion_components'], axis=0) if len(sc['motion_components']) else np.zeros_like(sc['u'])
    s['u'] = sc['u']-comp; mods.append(dict(s))
    s['act'] = np.zeros_like(sc['act']); mods.append(dict(s))
    s['u'] = np.full_like(sc['u'], sc['u'][0]); mods.append(dict(s))
    s = dict(s, mu_k=0., mu_s=0., k_d=0.); mods.append(dict(s))
    s = dict(s, tool_mass=0., tool_mass_hat=0.); mods.append(dict(s))
    s = dict(s, alpha=np.zeros_like(np.asarray(sc['alpha'], float))); s.pop('memory_start', None); mods.append(dict(s))
    reps = []
    for m in mods:
        try:
            reps.append(E.replay(m, pol, sent, case['seed'], n, tau_h, mag_h, eb, policy))
        except Exception as ex:      # a counterfactual plant that cannot run is reported, not hidden
            reps.append(None)
    t = log[:n, 0]
    k1 = reps[0][2]
    seg = V.segments(t, k1, sc['event_info'], k1 > 200.)
    out = dict(index=index, group=case['group'], failed=[STEPS[i] for i, r in enumerate(reps) if r is None], rows=[])
    for j in range(100, n):
        sname = next((q for q in ('rising', 'falling', 'reversal') if seg[q][j]), None)
        if sname is None:
            continue
        vals = []
        for r in reps:
            if r is None:
                vals.append((np.nan, np.nan, np.nan)); continue
            Y, z, kk = r
            if kk[j] <= 200.:
                vals.append((np.nan, np.nan, np.nan)); continue
            q24, k24 = local(Y[j], z[j], t, j, 24)
            _, k64 = local(Y[j], z[j], t, j, 64)
            if q24 is None or k64 is None:
                vals.append((np.nan, np.nan, np.nan)); continue
            rr = np.arange(j-24, j)
            kr = .5*(kk[rr]+kk[rr+1]); tau = t[rr]+.0035-t[j-1]
            v = np.polyfit(tau, kr, 1)[0]
            dev = kr-kk[j]-v*tau
            curv = float(q24@np.r_[0., np.cumsum(dev[1:]*np.diff(Y[j][-24:]))])/kk[j]
            vals.append((k24/kk[j]-1., k64/kk[j]-1., curv))
        out['rows'].append((sname, vals))
    W.NW = 32
    return out


def main():
    with ProcessPoolExecutor(7) as ex:
        res = list(ex.map(job, CASES, chunksize=1))
    rep = dict(failed={r['index']: r['failed'] for r in res if r['failed']}, groups={})
    print('failed counterfactuals', rep['failed'])
    names = ('none',)+STEPS
    for g in ('A', 'combo', 'A_large', 'all'):
        rows = [row for r in res if g == 'all' or r['group'] == g for row in r['rows']]
        if not rows:
            continue
        V0 = np.array([[v for v in vals] for _, vals in rows])      # windows x steps x (e24, e64, curv)
        harm = (np.abs(V0[:, 0, 1])-np.abs(V0[:, 0, 0]) > .05) & np.isfinite(V0[:, 0, 0])
        cell = dict(windows=len(rows), harmful=int(harm.sum()), steps={})
        print(f'--- {g}: {len(rows)} change windows, {harm.sum()} with a harmful long span after noise removal.')
        print('    paired step i -> i+1 on harmful windows valid in both: n | median e24 before -> after | median |e24| before -> after | median e64 before -> after | median curvature after')
        for i in range(1, len(STEPS)):
            A_, B_ = V0[harm, i-1], V0[harm, i]
            ok = np.isfinite(A_[:, 0]) & np.isfinite(B_[:, 0]) & np.isfinite(A_[:, 1]) & np.isfinite(B_[:, 1])
            if ok.sum() < 5:
                continue
            c = dict(n=int(ok.sum()), e24=[float(np.median(A_[ok, 0])), float(np.median(B_[ok, 0]))], abs24=[float(np.median(np.abs(A_[ok, 0]))), float(np.median(np.abs(B_[ok, 0])))],
                     e64=[float(np.median(A_[ok, 1])), float(np.median(B_[ok, 1]))], curv=float(np.median(B_[ok, 2])), lost=int((np.isfinite(A_[:, 0]) & ~np.isfinite(B_[:, 0])).sum()))
            cell['steps'][STEPS[i]] = c
            print(f"    -{STEPS[i]:15s} n {c['n']:4d} (lost {c['lost']:3d}) | {c['e24'][0]:+.3f} -> {c['e24'][1]:+.3f} | {c['abs24'][0]:.3f} -> {c['abs24'][1]:.3f} | {c['e64'][0]:+.3f} -> {c['e64'][1]:+.3f} | {c['curv']:+.3f}")
        rep['groups'][g] = cell
    (HERE/'results/v2_shared_bias_source.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
