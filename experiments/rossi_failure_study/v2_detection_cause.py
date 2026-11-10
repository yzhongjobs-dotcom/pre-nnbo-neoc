"""Why does adaptive history miss the moment it becomes invalid? Recorded windows, noisy and noise-free (same commands).

A. Span compatibility (the test behind ICI and multi-scale history). For local-linear endpoint estimates over 24 and over 64 ticks:
   harm  = |e_64| - |e_24| on noise-free data (the long span is worse, so history should have been cut);
   D0    = k_64 - k_24 noise-free, the part of the long-short difference that more information could reveal;
   s_d   = noise sd of (k_64 - k_24) from the exact weights;
   common = e_64 - D0 = e_24: the bias both spans share and no compatibility test can see.
   Each harmful window (harm > .05) is classed: detectable (|D0| >= 2 s_d), noise-limited (|D0| < 2 s_d but |D0| >= .5 |e_64|),
   or common bias (|D0| < .5 |e_64|): even with no noise the test would not flag it, because the short span carries the same error.
B. V2 trend statistic of eq. (11): noise-free v0 against the true slope (bias) and the noise sd; missed falls classed the same way.
Writes results/v2_detection_cause.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import v2_history_ceiling as C
import v2_verify as V

HERE = Path(__file__).resolve().parent
TC, SIG = .008, H.SIGMA


def weights(dy, t, j, L):
    r = np.arange(j-L, j)
    if r[0] < 1 or not np.isfinite(dy[r]).all():
        return None
    y = np.cumsum(dy[r]); h0 = y-y[-1]; h1 = (t[r]+.0035-t[j-1])*h0
    M = C.mat(L); mh0, mh1 = M@h0, M@h1
    a, b, c1 = h0@mh0, h0@mh1, h1@mh1
    det = a*c1-b*b
    return None if det <= 0 else (c1*mh0-b*mh1)/det


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    dy, z = C.rows(data, False); dy0, z0 = C.rows(data, True)
    rows = []
    for j in range(100, n):
        if k[j] <= 200.:
            continue
        s = next((q for q in ('rising', 'falling', 'reversal') if seg[q][j]), None)
        if s is None:
            continue
        q24, q64 = weights(dy, t, j, 24), weights(dy, t, j, 64)
        q24n, q64n = weights(dy0, t, j, 24), weights(dy0, t, j, 64)
        if any(q is None for q in (q24, q64, q24n, q64n)):
            continue
        r24, r64 = np.arange(j-24, j), np.arange(j-64, j)
        if not (np.isfinite(z[r64]).all() and np.isfinite(z0[r64]).all()):
            continue
        kT = k[j]
        e24n = float(q24n@z0[r24])/kT-1; e64n = float(q64n@z0[r64])/kT-1
        # Noise sd of k64 - k24: both are linear in the same tick noise; pad q24 to the 64-tick support.
        qd = q64.copy(); qd[-24:] -= q24
        sd = SIG*math.sqrt(float(qd@qd))/kT
        rows.append((s, meta['group'], e24n, e64n, e64n-e24n, sd))
    # B. trend statistic, eq. (11), from the 24-tick window with W = I.
    trend = []
    for j in range(NW_START, n):
        if k[j] <= 200.:
            continue
        s = next((q for q in ('plateau', 'rising', 'falling', 'reversal') if seg[q][j]), None)
        if s is None:
            continue
        Y, Y1, zz = data['Y'][j, -24:], data['H1'][j, -24:], data['z'][j, -24:]
        Y0, Y10, zz0 = data['Y0'][j, -24:], data['H10'][j, -24:], data['z0'][j, -24:]
        if not all(np.isfinite(a).all() for a in (Y, Y1, zz, Y0, Y10, zz0)):
            continue
        M = C.mat(24)
        def stat(h0, h1, zv):
            a0, b0 = h0@M@h0, h0@M@h1; g = h1-(b0/a0)*h0; lv = M@g/(g@M@g)
            return float(lv@zv), SIG*math.sqrt(float(lv@lv))
        vt, sv = stat(Y, Y1, zz); v0, _ = stat(Y0, Y10, zz0)
        rr = np.arange(j-24, j); kr = .5*(k[rr]+k[rr+1]); vtrue = float(np.polyfit(t[rr]-t[j-1], kr, 1)[0])
        kT = k[j]
        trend.append((s, meta['group'], vt/kT, v0/kT, vtrue/kT, sv/kT))
    return dict(rows=rows, trend=trend)


NW_START = 100


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rows = [x for r in res for x in r['rows']]; trend = [x for r in res for x in r['trend']]
    rep = dict(span={}, trend={})
    print('A. 64- vs 24-tick local linear, noise-free. Harmful windows (|e64| > |e24| + .05) and their class')
    for g in ('all', 'A', 'A_large', 'A_slow_thin', 'combo'):
        for s in ('rising', 'falling', 'reversal'):
            a = [x for x in rows if x[0] == s and (g == 'all' or x[1] == g)]
            if len(a) < 30:
                continue
            e24 = np.array([x[2] for x in a]); e64 = np.array([x[3] for x in a]); D0 = np.array([x[4] for x in a]); sd = np.array([x[5] for x in a])
            harm = np.abs(e64)-np.abs(e24) > .05
            det = harm & (np.abs(D0) >= 2*sd)
            common = harm & (np.abs(D0) < .5*np.abs(e64))
            noise = harm & ~det & ~common
            cell = dict(n=len(a), harmful=float(harm.mean()), detectable=float(det.sum()/max(harm.sum(), 1)), noise_limited=float(noise.sum()/max(harm.sum(), 1)),
                        common_bias=float(common.sum()/max(harm.sum(), 1)), sd_diff=float(np.median(sd)), D0_harm=float(np.median(np.abs(D0[harm]))) if harm.any() else None,
                        e64_harm=float(np.median(np.abs(e64[harm]))) if harm.any() else None, e24_harm=float(np.median(np.abs(e24[harm]))) if harm.any() else None,
                        helpful=float((np.abs(e24)-np.abs(e64) > .05).mean()))
            rep['span'][f'{g}|{s}'] = cell
            if not harm.any():
                continue
            print(f"  {g:12s} {s:9s} n {len(a):5d} harmful {cell['harmful']:.2f} helpful {cell['helpful']:.2f} | of harmful: detectable {cell['detectable']:.2f} noise-limited {cell['noise_limited']:.2f} common {cell['common_bias']:.2f} | median |D0| {cell['D0_harm']:.3f} vs noise sd {cell['sd_diff']:.3f}, |e64| {cell['e64_harm']:.3f}, |e24| {cell['e24_harm']:.3f}")
    print('B. trend statistic of eq. (11), relative slope (1/s): noise-free bias rms, noise sd, true |v|; missed changes (|v| > 1/s, noisy |v~| < 2 s_v) classed')
    for g in ('all', 'A', 'A_large', 'A_slow_thin', 'cough', 'combo'):
        for s in ('plateau', 'rising', 'falling', 'reversal'):
            a = np.array([x[2:] for x in trend if x[0] == s and (g == 'all' or x[1] == g)])
            if len(a) < 30:
                continue
            vt, v0, vtrue, sv = a.T
            bias = v0-vtrue
            big = np.abs(vtrue) > 1.
            missed = big & (np.abs(vt) < 2*sv)
            cont = missed & (np.abs(v0) < 2*sv) & (np.abs(bias) > .5*np.abs(vtrue))
            noisel = missed & (np.abs(v0) >= 2*sv)
            cell = dict(n=len(a), bias_rms=float(np.sqrt(np.mean(bias**2))), noise_sd=float(np.median(sv)), v_abs=float(np.median(np.abs(vtrue))), corr0=float(np.corrcoef(v0, vtrue)[0, 1]),
                        big=int(big.sum()), missed=float(missed.sum()/max(big.sum(), 1)), missed_contaminated=float(cont.sum()/max(missed.sum(), 1)), missed_noise=float(noisel.sum()/max(missed.sum(), 1)))
            rep['trend'][f'{g}|{s}'] = cell
            print(f"  {g:12s} {s:9s} n {len(a):5d} | noise-free bias rms {cell['bias_rms']:.2f}, corr(v0, v) {cell['corr0']:+.2f} | noise sd {cell['noise_sd']:.2f} | |v| {cell['v_abs']:.2f} | |v|>1: {cell['big']:5d}, missed {cell['missed']:.2f} (contaminated {cell['missed_contaminated']:.2f}, noise would reveal {cell['missed_noise']:.2f})")
    (HERE/'results/v2_detection_cause.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
