"""Replay diagnosis of the two mismatches, on recorded observations (results/huber_windows), one factor at a time.

Cough: does the tick filter over-estimate at cough onset, with or without sensor noise, and does it need the fast k mode?
  T2 noisy / T2 noise-free / T1 (quiet k mode only) / window D4. Signed ln(k/k_eff) in each event, first .15 s and whole event.
Combo collapse: is it noise, the disturbance, or a mismatch between the plant's gain for non-probe motion and k_eff?
  T2 on z_syn = cumulative k_eff dY + recorded noise (recorded Y, no d): if no collapse, the cause is in z, not in Y.
  Per collapse window, noise-free: least-squares gain of Mz on the probe band (20.8 Hz +- 4) and on the rest of MY.
Predictions fixed before running: cough over-estimate is noise-free and vanishes in T1 (the fast k mode is what takes the sharp d);
combo collapse remains noise-free, vanishes on z_syn, and the non-probe gain in those windows is far below k_eff.
Writes results/memory_mismatch_replay.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import memory_tick_kalman as T

HERE = Path(__file__).resolve().parent
CFG = dict(form='incr', d='d4', sigma=T.SIG_TRUE)
TC, F0 = .008, 20.8


def run_filter(f, dY, z):
    n = len(z); est = np.full(n, np.nan)
    for j in range(n):
        if not np.isfinite(z[j]):
            if f.started:
                f.reset(f.x[0][0])
            continue
        est[j] = f.step(dY[j], z[j])
    return est


def series(data, which):
    Y = data['Y0' if which else 'Y']; Z = data['z0' if which else 'z']
    n = len(Y); dY = np.full(n, np.nan); z = np.full(n, np.nan)
    for j in range(n):
        if np.isfinite(Y[j, -1]):
            z[j] = Z[j, -1]
            if np.isfinite(Y[j, -2]):
                dY[j] = Y[j, -1]-Y[j, -2]
    return dY, z


def band_gain(Y, z):
    """Gain of Mz on the probe-band part and on the rest of MY (noise-free window)."""
    S = H.Solver(); a = S.fit(Y, z, iters=1)
    my = a['my']; mz = z-a['sz']
    m = len(my); t = np.arange(m)*TC
    basis = np.c_[np.sin(2*np.pi*F0*t), np.cos(2*np.pi*F0*t), np.sin(2*np.pi*(F0-4)*t), np.cos(2*np.pi*(F0-4)*t), np.sin(2*np.pi*(F0+4)*t), np.cos(2*np.pi*(F0+4)*t)]
    P = basis@np.linalg.pinv(basis)
    yp = P@my; yr = my-yp
    gp = float(yp@mz/(yp@yp)) if yp@yp > 0 else np.nan
    gr = float(yr@mz/(yr@yr)) if yr@yr > 0 else np.nan
    return gp, gr, float(np.sqrt(yr@yr/max(yp@yp, 1e-30)))


def case(path):
    data = H.load(path); meta = data['meta']; log, masks = H.masks_for(data)
    n = len(data['Y']); k = log[:n, 26]; t = log[:n, 0]
    dY, z = series(data, False); dY0, z0 = series(data, True)
    rng = np.random.default_rng(meta['index'])
    noise = z-z0
    # z_syn: what an exact model z = k_eff Y would give on the recorded Y, plus the recorded noise.
    zs = np.full(n, np.nan); acc = None
    for j in range(n):
        if not np.isfinite(z[j]):
            acc = None; continue
        if acc is None or not np.isfinite(dY[j]):
            acc = z[j]
        else:
            acc = acc+.5*(k[j]+k[j-1])*dY[j]
        zs[j] = acc+noise[j]
    est = dict(T2=run_filter(T.TickKalman(q_rel=(.002, .5), **CFG), dY, z),
               T2_noisefree=run_filter(T.TickKalman(q_rel=(.002, .5), **CFG), dY0, z0),
               T1_quiet=run_filter(T.TickKalman(q_rel=.002, **CFG), dY, z),
               T2_exact_model=run_filter(T.TickKalman(q_rel=(.002, .5), **CFG), dY, zs))
    S = H.Solver(); win = np.full(n, np.nan)
    for j in range(n):
        Yw, zw = H.tail(data['Y'][j], 24), H.tail(data['z'][j], 24)
        if len(Yw) >= 8:
            win[j] = H.clip(S.fit(Yw, zw, iters=1)['k'])
    est['window'] = win
    lk = np.log(np.where(k > 200., k, np.nan))
    err = {q: np.log(np.where(v > 0, v, np.nan))-lk for q, v in est.items()}
    ev = []
    for e in meta['event_info']:
        if e['kind'] != 'cough':
            continue
        onset = (t >= e['start']) & (t <= e['start']+.15); whole = (t >= e['start']) & (t <= e['end']+.2)
        ev.append(dict(group=meta['group'], index=meta['index'], start=e['start'], **{f'{q}_onset': float(np.nanmedian(err[q][onset])) for q in err},
                       **{f'{q}_onset_max': float(np.nanmax(err[q][onset])) for q in err}, **{f'{q}_whole': float(np.nanmedian(err[q][whole])) for q in err}))
    collapse = []
    ok = (t >= 1.) & np.isfinite(err['T2'])
    for j in np.flatnonzero(ok & (err['T2'] < math.log(.5))):
        Y0w, z0w = H.tail(data['Y0'][j], 24), H.tail(data['z0'][j], 24)
        if len(Y0w) != 24:
            continue
        gp, gr, ratio = band_gain(Y0w, z0w)
        collapse.append((meta['group'], float(k[j]), gp/k[j], gr/k[j], ratio))
    # Same quantities on a sample of normal windows, for contrast.
    normal = []
    for j in np.flatnonzero(ok & (np.abs(err['T2']) < .1))[::25]:
        Y0w, z0w = H.tail(data['Y0'][j], 24), H.tail(data['z0'][j], 24)
        if len(Y0w) == 24:
            gp, gr, ratio = band_gain(Y0w, z0w)
            normal.append((meta['group'], float(k[j]), gp/k[j], gr/k[j], ratio))
    severe = {q: int(np.sum(ok & (err[q] < math.log(.5)))) for q in err}
    return dict(index=meta['index'], group=meta['group'], events=ev, collapse=collapse, normal=normal, severe=severe)


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    names = ('window', 'T2', 'T2_noisefree', 'T1_quiet', 'T2_exact_model')
    print('--- cough events: median signed ln(k/k_eff) in first .15 s, max in first .15 s, median over event')
    ev = [e for r in res for e in r['events']]
    out = dict(cough={})
    for q in names:
        on = np.array([e[f'{q}_onset'] for e in ev]); mx = np.array([e[f'{q}_onset_max'] for e in ev]); wh = np.array([e[f'{q}_whole'] for e in ev])
        out['cough'][q] = dict(onset_median=float(np.nanmedian(on)), onset_max_median=float(np.nanmedian(mx)), whole_median=float(np.nanmedian(wh)), n=len(ev))
        print(f'  {q:15s} onset {np.nanmedian(on):+.3f}  onset max {np.nanmedian(mx):+.3f}  whole {np.nanmedian(wh):+.3f}  (n {len(ev)})')
    print('--- severe ticks (below half of truth) by group')
    out['severe'] = {}
    for g in H.GROUPS:
        out['severe'][g] = {q: sum(r['severe'][q] for r in res if r['group'] == g) for q in names}
        print(f'  {g:12s}', out['severe'][g])
    print('--- noise-free gain of Mz, relative to k_eff: probe band, rest of MY; rest/probe amplitude')
    out['gain'] = {}
    for label, key in (('collapse', 'collapse'), ('normal', 'normal')):
        for g in ('combo', 'A_large', 'all'):
            a = np.array([x[1:] for r in res for x in r[key] if g == 'all' or x[0] == g])
            if not len(a):
                continue
            out['gain'][f'{label}|{g}'] = dict(n=len(a), k=float(np.median(a[:, 0])), probe=float(np.median(a[:, 1])), rest=float(np.median(a[:, 2])), rest_p10=float(np.percentile(a[:, 2], 10)), amp=float(np.median(a[:, 3])))
            s = out['gain'][f'{label}|{g}']
            print(f"  {label:8s} {g:8s} n {s['n']:4d} k {s['k']:7.0f} probe {s['probe']:.2f} rest {s['rest']:+.2f} (P10 {s['rest_p10']:+.2f}) rest/probe amp {s['amp']:.1f}")
    (HERE/'results/memory_mismatch_replay.json').write_text(json.dumps(out, indent=1), encoding='utf-8')


if __name__ == '__main__' and len(sys.argv) == 1:
    main()


def follows_u(path):
    """In T2 collapse windows and in normal windows: does the non-probe part of MY follow the tissue motion u (log col 14)?"""
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(data['Y']); k = log[:n, 26]; t = log[:n, 0]; u = log[:n, 14]
    dY, z = series(data, False)
    est = run_filter(T.TickKalman(q_rel=(.002, .5), **CFG), dY, z)
    e = np.log(np.where(est > 0, est, np.nan))-np.log(np.where(k > 200., k, np.nan))
    S = H.Solver(); out = []
    tt = np.arange(24)*TC
    basis = np.c_[np.sin(2*np.pi*F0*tt), np.cos(2*np.pi*F0*tt), np.sin(2*np.pi*(F0-4)*tt), np.cos(2*np.pi*(F0-4)*tt), np.sin(2*np.pi*(F0+4)*tt), np.cos(2*np.pi*(F0+4)*tt)]
    P = basis@np.linalg.pinv(basis)
    for j in range(30, n):
        if t[j] < 1. or not np.isfinite(e[j]):
            continue
        kind = 'collapse' if e[j] < math.log(.5) else 'normal' if abs(e[j]) < .1 and j % 25 == 0 else None
        if kind is None:
            continue
        Y0 = H.tail(data['Y0'][j], 24)
        if len(Y0) != 24:
            continue
        a = S.fit(Y0, H.tail(data['z0'][j], 24), iters=1)
        yr = a['my']-P@a['my']
        uw = u[j-24:j]; mu = uw-S.fit(Y0, uw, iters=1)['sz']; ur = mu-P@mu
        if yr@yr <= 0 or ur@ur <= 0:
            continue
        out.append((meta['group'], kind, float(yr@ur/math.sqrt((yr@yr)*(ur@ur))), float(yr@ur/(ur@ur))))
    return out


def main_u():
    with ProcessPoolExecutor(8) as ex:
        rows = [x for r in ex.map(follows_u, sorted(H.WIN.glob('case_*.npz'))) for x in r]
    res = {}
    for kind in ('collapse', 'normal'):
        for g in ('combo', 'all'):
            a = np.array([x[2:] for x in rows if x[1] == kind and (g == 'all' or x[0] == g)])
            if len(a):
                res[f'{kind}|{g}'] = dict(n=len(a), corr=float(np.median(a[:, 0])), slope=float(np.median(a[:, 1])))
                print(kind, g, res[f'{kind}|{g}'])
    p = HERE/'results/memory_mismatch_replay.json'
    d = json.loads(p.read_text(encoding='utf-8')); d['nonprobe_vs_u'] = res; p.write_text(json.dumps(d, indent=1), encoding='utf-8')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'u':
    main_u()


def spectrum(path):
    """Per DFT bin of the 24-tick window (5.2 Hz apart): share of non-probe MY energy, and the plant gain Re(Mz/MY)/k_eff, noise-free."""
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(data['Y']); k = log[:n, 26]; t = log[:n, 0]
    dY, z = series(data, False)
    est = run_filter(T.TickKalman(q_rel=(.002, .5), **CFG), dY, z)
    e = np.log(np.where(est > 0, est, np.nan))-np.log(np.where(k > 200., k, np.nan))
    S = H.Solver(); out = []
    for j in range(30, n):
        if t[j] < 1. or not np.isfinite(e[j]):
            continue
        kind = 'collapse' if e[j] < math.log(.5) else 'normal' if abs(e[j]) < .1 and j % 25 == 0 else None
        if kind is None:
            continue
        Y0, z0 = H.tail(data['Y0'][j], 24), H.tail(data['z0'][j], 24)
        if len(Y0) != 24:
            continue
        a = S.fit(Y0, z0, iters=1)
        fy = np.fft.rfft(a['my']); fz = np.fft.rfft(z0-a['sz'])
        out.append((meta['group'], kind, float(k[j]), (np.abs(fy)**2).tolist(), (fz/np.where(np.abs(fy) > 0, fy, np.nan)).real.tolist()))
    return out


def main_spec():
    with ProcessPoolExecutor(8) as ex:
        rows = [x for r in ex.map(spectrum, sorted(H.WIN.glob('case_*.npz'))) for x in r]
    f = np.fft.rfftfreq(24, TC)
    res = {}
    for kind in ('collapse', 'normal'):
        q = [x for x in rows if x[1] == kind]
        P = np.array([x[3] for x in q]); P = P/P.sum(axis=1, keepdims=True)
        G = np.array([x[4] for x in q])/np.array([x[2] for x in q])[:, None]
        res[kind] = dict(n=len(q), freq=f.tolist(), energy=np.median(P, axis=0).tolist(), gain=np.nanmedian(G, axis=0).tolist())
        print(kind, len(q), 'median k', round(float(np.median([x[2] for x in q]))))
        print('  Hz     ', ' '.join(f'{v:6.1f}' for v in f))
        print('  energy ', ' '.join(f'{v:6.3f}' for v in res[kind]['energy']))
        print('  gain/k ', ' '.join(f'{v:6.2f}' for v in res[kind]['gain']))
    p = HERE/'results/memory_mismatch_replay.json'
    d = json.loads(p.read_text(encoding='utf-8')); d['spectrum'] = res; p.write_text(json.dumps(d, indent=1), encoding='utf-8')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'spec':
    main_spec()
