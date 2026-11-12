"""Replay of the compliance-domain two-channel filter (dual_channel_v3) on the 42 stored closed loops.

Differences from the synthetic form: the probe output is a 24-tick window average, so its regressor is the 24-tick mean of y and the
background is taken 11.5 ticks back; the output uses the current encoder sample (the probe ripple is below 1% of the indentation).
Breathing frequency: 'peak' = spectral peak of y over the whole episode (diagnostic, not causal); 'fixed' = 0.25 Hz.
Arms A0 (raw window) and A1 (probe-only IMM, parameters selected in dual_channel_replay.py) are the comparisons.
"""
import itertools, json, math, sys
from pathlib import Path

import numpy as np

import dual_channel_replay as D

TC = .008
HERE = Path(__file__).resolve().parent


def run(case, pars, f0):
    t, ys6, fs, m, R, valid, finite = case['data']
    y = case['y']; yw = case['yw']
    qs, qo, qbeta = pars
    dt = TC; tau = 11.5*TC; w = 2*math.pi*f0
    A = np.eye(5); A[0, 1] = dt
    cw, sw = math.cos(w*dt), math.sin(w*dt); A[2, 2], A[2, 3], A[3, 2], A[3, 3] = cw, sw, -sw, cw
    ct, st = math.cos(w*tau), math.sin(w*tau)
    Js = np.array([[dt**3/3, dt**2/2], [dt**2/2, dt]])
    n = len(t); out = np.full(n, np.nan); sd = np.full(n, np.nan); g = np.full(n, np.nan)
    active = False; last = -1000
    for j in range(n):
        if finite[j]:
            last = j
        if fs[j] < 1. or j-last > 25 or not np.isfinite(fs[j]+y[j]):
            active = False
            if j:
                out[j], sd[j], g[j] = out[j-1], sd[j-1], g[j-1]
            continue
        if not active:
            if not valid[j]:
                if j:
                    out[j], sd[j], g[j] = out[j-1], sd[j-1], g[j-1]
                continue
            C0 = fs[j]/math.exp(m[j]); y0 = y[j]; b0 = 1./3.
            x = np.array([C0, 0., 0., 0., b0]); P = np.diag([C0*C0*R[j]/24, (C0*.5)**2, (C0*.5)**2, (C0*.5)**2, (b0*1.1)**2]); active = True
        else:
            c2 = max(x[0]+x[2]+x[4]*(y[j]-y0), 1e-6)**2
            x = A@x; P = A@P@A.T
            P[:2, :2] += qs*c2*Js; P[2, 2] += qo*c2*dt; P[3, 3] += qo*c2*dt; P[4, 4] += qbeta*x[4]*x[4]
            if valid[j] and j >= 24 and fs[j-12] >= 1.:
                Cobs = fs[j-12]/math.exp(m[j])
                H = np.array([1., -tau, ct, -st, yw[j]-y0]); Rc = Cobs*Cobs*R[j]
                PH = P@H; S = float(H@PH)+Rc; K = PH/S
                x = x+K*(Cobs-H@x); I = np.eye(5)-np.outer(K, H); P = I@P@I.T+Rc*np.outer(K, K)
        H0 = np.array([1., 0., 1., 0., y[j]-y0]); C = max(float(H0@x), 1e-7)
        out[j] = math.log(fs[j])-math.log(C); sd[j] = math.sqrt(max(float(H0@P@H0), 0.))/C; g[j] = 1./x[4] if x[4] > 0 else np.nan
    return out, sd, g


def peak_freq(y):
    yy = y-np.polyval(np.polyfit(np.arange(len(y)), y, 1), np.arange(len(y)))
    f = np.arange(.12, .5, .005); tt = np.arange(len(y))*TC
    p = [abs(np.sum(yy*np.exp(-2j*np.pi*ff*tt))) for ff in f]
    return float(f[int(np.argmax(p))])


if __name__ == '__main__':
    root = HERE.parents[1]
    cases, _ = D.load_cases(root)
    for c in cases:
        c['y'] = np.array(c['log'][:, 2], copy=True); c['yw'] = D.trailing(c['y'], 24); c['f_peak'] = peak_freq(c['y'])
        c['pred']['A1'] = D.probe_filter(c['data'], (0., .01, .99))
    base = [c for c in cases if c['group'] in D.BASE_GROUPS]
    dev = [c for c in base if c['index'] % 2 == 0]; test = [c for c in base if c['index'] % 2 == 1]
    grid = list(itertools.product((1e-3, 1e-2, 1e-1, 1.), (1e-3, 1e-2, 1e-1, 1., 10.), (1e-7, 1e-5)))
    report = {}
    for mode in ('peak', 'fixed'):
        best = None
        for pars in grid:
            obj = float(np.mean([D.objective(run(c, pars, c['f_peak'] if mode == 'peak' else .25)[0], c) for c in dev]))
            if best is None or obj < best[0]:
                best = (obj, pars)
        print(mode, 'selected', best, 'A1 dev', float(np.mean([D.objective(c['pred']['A1'], c) for c in dev])), 'A0 dev', float(np.mean([D.objective(c['pred']['A0'], c) for c in dev])), flush=True)
        for c in cases:
            c['pred']['F_'+mode], c['sd_'+mode], c['g_'+mode] = run(c, best[1], c['f_peak'] if mode == 'peak' else .25)
        report[mode] = best
    arms = ('A0', 'A1', 'F_peak', 'F_fixed')
    print('\nTest cases (odd), rms of ln(k_hat/k_eff) by segment; cough and combo: all cases')
    for grp in ('A', 'A_slow_thin', 'A_large', 'B', 'A+B', 'cough', 'combo'):
        cs = [c for c in (test if grp in D.BASE_GROUPS else cases) if c['group'] == grp]
        print(grp, [c['index'] for c in cs])
        for seg in ('plateau', 'slow', 'rise', 'fall', 'reversal', 'event'):
            row = []
            for a in arms:
                e = np.concatenate([(c['pred'][a]-c['truth'])[c['masks'][seg]] for c in cs])
                row.append(f'{a} {D.rms(e):.3f}' if len(e) > 20 else f'{a}  -  ')
            n = sum(int(c['masks'][seg].sum()) for c in cs)
            if n > 20:
                print(f'   {seg:9s} n={n:5d}  '+' | '.join(row))
        print('   sustained under-estimate, median/worst: '+' | '.join(
            f"{a} {np.median([D.sustained(c['pred'][a]-c['truth'], c['masks']['outside_event']) or 0. for c in cs]):.2f}/{max(D.sustained(c['pred'][a]-c['truth'], c['masks']['outside_event']) or 0. for c in cs):.2f}" for a in arms))
    print('\nPer case, A, A_slow_thin, A_large: rms on dynamic ticks (slow+rise+fall+reversal) and plateau; learned g at end; f_peak; mean self-reported sd')
    for c in cases:
        if c['group'] in ('A', 'A_slow_thin', 'A_large'):
            dyn = c['masks']['dynamic']; pl = c['masks']['plateau']
            print(f"{c['index']:2d} {c['group']:11s} {'dev ' if c['index'] % 2 == 0 else 'test'} dyn n={int(dyn.sum()):4d} "+' '.join(f"{a} {D.rms((c['pred'][a]-c['truth'])[dyn]):.3f}" for a in arms)
                  + ' | plateau '+' '.join(f"{D.rms((c['pred'][a]-c['truth'])[pl]):.3f}" for a in arms)+f" | g {c['g_peak'][-1]:.1f} f {c['f_peak']:.2f} sd {np.nanmean(c['sd_peak'][c['masks']['all']]):.3f}")
    json.dump({k: dict(objective=v[0], pars=v[1]) for k, v in report.items()}, open(HERE/'results/dual_channel_v3_replay.json', 'w'))
