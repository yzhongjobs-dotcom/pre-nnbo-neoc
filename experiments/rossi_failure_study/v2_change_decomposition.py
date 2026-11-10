"""Where do V2's rising, falling and reversal errors come from, and what does the slope actually remove? Recorded windows only, no plant.

For a window with frozen Huber weights M_W and any linear estimator q with q'h0 = 1 (and q'P = 0): k_hat - k_T = q'(z - k_T h0). Exact split, all relative to k_T:
  noise     q'[(z - z0) - k_T (h0 - Y0)]                 sensor noise, including encoder noise in the regressor
  calib     q'[(z0 - z0i) - k_T (Y0 - Y0i)]              calibration error of the sensing stage
  k_lin     q'(v h1_0)                                    linear part of the within-window stiffness change, v the true slope over the window
  k_curv    q'c - k_lin, c = cumulative (k_i - k_T) dY0   curvature of the change (reversals)
  other     remainder: disturbance not absorbed by D4, plant dynamics, heterogeneity
Estimators, same M_W: const (Huber-D4, sv2 = 0), oracle (sv2 = v_true^2), free (sv2 = inf).
Also per window: slope SNR |v|/sigma_v from the free 2x2 information; V2 trend statistic of eq. (11); window + IMM and tick Kalman outputs
on the same observations for the memory routes.
Writes results/v2_change_decomposition.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import memory_tick_kalman as TK
import v2_verify as V

HERE = Path(__file__).resolve().parent
NW, TC, SIG = 24, .008, H.SIGMA
TICKS = Path(r'D:/pip_env/temp/claude/8fbcb2674095e6be/4c5dc125-7da2-41f3-aa59-d2d3f2f2c904/scratchpad/unified_ticks/24')
PARTS = ('noise', 'calib', 'k_lin', 'k_curv', 'other')
EST = ('const', 'oracle', 'free')


def q_of(M, h0, h1, lam):
    mh0, mh1 = M@h0, M@h1
    a, b, c1 = h0@mh0, h0@mh1, h1@mh1
    if lam == 0.:
        return mh0/a
    if not np.isfinite(lam):
        return np.zeros_like(h0)+mh0*0.+(mh0-(b/(c1+lam))*mh1)/(a-b*b/(c1+lam))
    cp = c1+lam
    return (mh0-(b/cp)*mh1)/(a-b*b/cp)


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    names = [f'{e}_{p}' for e in EST for p in PARTS]+[f'{e}_err' for e in EST]+['snr', 'vrel', 'vt_stat', 'kT']
    rec = {q: np.full(n, np.nan) for q in names}
    for j in range(NW+1, n):
        Y, Y1, z = data['Y'][j, -NW:], data['H1'][j, -NW:], data['z'][j, -NW:]
        Y0, Y10, z0 = data['Y0'][j, -NW:], data['H10'][j, -NW:], data['z0'][j, -NW:]
        Y0i, z0i = data['Y0i'][j, -NW:], data['z0i'][j, -NW:]
        if not all(np.isfinite(a).all() for a in (Y, Y1, z, Y0, Y10, z0, Y0i, z0i)) or k[j] <= 200.:
            continue
        r = np.arange(j-NW, j); kr = .5*(k[r]+k[r+1])
        if np.any(kr <= 200.):
            continue
        kT = k[j]
        tt = .5*(t[r]+t[r+1])-t[j-1]
        v = float(np.polyfit(tt, kr, 1)[0])
        _, _, w, _ = V.solve_v2_window(Y, Y1, z, rho=V.RHO, sd=V.SD, sigma_z=SIG, sv2=0.)
        M = V.m_of(w)
        c = np.r_[0., np.cumsum((kr[1:]-kT)*np.diff(Y0))]
        resid = {'noise': (z-z0)-kT*(Y-Y0), 'calib': (z0-z0i)-kT*(Y0-Y0i)}
        for e, lam in (('const', 0.), ('oracle', SIG**2/max(v*v, 1e-9)), ('free', 0.0 if False else np.inf)):
            if e == 'free':
                mh0, mh1 = M@Y, M@Y1
                a, b, c1 = Y@mh0, Y@mh1, Y1@mh1
                q = (c1*mh0-b*mh1)/(a*c1-b*b)
            else:
                q = q_of(M, Y, Y1, lam)
            rec[f'{e}_err'][j] = float(q@(z-kT*Y))/kT
            rec[f'{e}_noise'][j] = float(q@resid['noise'])/kT
            rec[f'{e}_calib'][j] = float(q@resid['calib'])/kT
            lin = float(q@(v*Y10))/kT
            rec[f'{e}_k_lin'][j] = lin
            rec[f'{e}_k_curv'][j] = float(q@c)/kT-lin
            rec[f'{e}_other'][j] = rec[f'{e}_err'][j]-sum(rec[f'{e}_{p}'][j] for p in ('noise', 'calib', 'k_lin', 'k_curv'))
        mh0, mh1 = M@Y, M@Y1
        a, b, c1 = Y@mh0, Y@mh1, Y1@mh1
        var_v = SIG**2*a/(a*c1-b*b)
        rec['snr'][j] = abs(v)/math.sqrt(var_v)
        rec['vrel'][j] = v/kT
        g = Y1-(b/a)*Y; Jv = g@M@g; lv = M@g/Jv
        rec['vt_stat'][j] = float(lv@z)**2/(SIG**2*float(lv@lv))      # (v~/s_v)^2 of eq. (11)-(12) with beta_v = 0
        rec['kT'][j] = kT
    # Memory routes on the same observations: window + IMM (replay ticks), tick Kalman (T2).
    with np.load(TICKS/f"Huber-D4_{meta['index']:02d}.npz") as f:
        imm = f['D4_imm__raw'].astype(float)[:n]
    dY, zz = TK.ticks(data)
    kf = TK.TickKalman(q_rel=(.002, .5), form='incr', d='d4', sigma=TK.SIG_TRUE); tk = np.full(n, np.nan)
    for j in range(n):
        if not np.isfinite(zz[j]):
            if kf.started:
                kf.reset(kf.x[0][0])
            continue
        tk[j] = kf.step(dY[j], zz[j])
    rec['imm_err'] = imm/np.where(k > 200., k, np.nan)-1.
    rec['tick_err'] = tk/np.where(k > 200., k, np.nan)-1.
    return dict(index=meta['index'], group=meta['group'], seg={s: m.tolist() for s, m in seg.items()}, rec={q: v.tolist() for q, v in rec.items()})


def lnerr(x):
    return np.log(np.maximum(1.+np.asarray(x), 1e-3))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    out = {}
    groups = ('all', 'A', 'A_slow_thin', 'A_large', 'combo', 'cough')
    for g in groups:
        for s in ('plateau', 'rising', 'falling', 'reversal', 'event'):
            pick = [r for r in res if g == 'all' or r['group'] == g]
            m = [np.asarray(r['seg'][s]) & np.isfinite(r['rec']['const_err']) for r in pick]
            if sum(x.sum() for x in m) < 20:
                continue
            cat = lambda key: np.concatenate([np.asarray(r['rec'][key])[mm] for r, mm in zip(pick, m)])
            row = dict(n=int(sum(x.sum() for x in m)))
            for e in EST+('imm', 'tick'):
                le = lnerr(cat(f'{e}_err')); le = le[np.isfinite(le)]
                row[e] = dict(med=float(np.median(np.abs(le))), rms=float(np.sqrt(np.mean(le**2))), bias=float(np.median(le)), p90=float(np.percentile(np.abs(le), 90)))
            for e in EST:
                row[e]['parts'] = {p: dict(mean=float(np.mean(cat(f'{e}_{p}'))), rms=float(np.sqrt(np.mean(cat(f'{e}_{p}')**2)))) for p in PARTS}
            snr = cat('snr'); vrel = cat('vrel')
            row['snr_med'] = float(np.median(snr)); row['vrel_med'] = float(np.median(vrel))
            out[f'{g}|{s}'] = row
    # Gain of the slope by slope SNR, all falling and rising windows outside events.
    bins = ((0, 1), (1, 2), (2, 4), (4, 1e9))
    out['by_snr'] = {}
    for s in ('rising', 'falling', 'reversal', 'event'):
        for lo, hi in bins:
            m = [np.asarray(r['seg'][s]) & np.isfinite(r['rec']['const_err']) & (np.asarray(r['rec']['snr']) >= lo) & (np.asarray(r['rec']['snr']) < hi) for r in res]
            if sum(x.sum() for x in m) < 10:
                continue
            cat = lambda key: np.concatenate([np.asarray(r['rec'][key])[mm] for r, mm in zip(res, m)])
            cell = dict(n=int(sum(x.sum() for x in m)), vrel=float(np.median(np.abs(cat('vrel')))))
            for e in EST+('imm', 'tick'):
                le = lnerr(cat(f'{e}_err')); le = le[np.isfinite(le)]
                cell[e] = dict(med=float(np.median(np.abs(le))), rms=float(np.sqrt(np.mean(le**2))), bias=float(np.median(le)))
            out['by_snr'][f'{s}|{lo}-{hi}'] = cell
    # Detection: how often the eq. (11) statistic exceeds kappa = 4 and 16 by segment and by SNR.
    out['detect'] = {}
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        m = [np.asarray(r['seg'][s]) & np.isfinite(np.asarray(r['rec']['vt_stat'])) for r in res]
        x = np.concatenate([np.asarray(r['rec']['vt_stat'])[mm] for r, mm in zip(res, m)])
        sn = np.concatenate([np.asarray(r['rec']['snr'])[mm] for r, mm in zip(res, m)])
        out['detect'][s] = dict(n=len(x), over4=float(np.mean(x > 4)), over16=float(np.mean(x > 16)), over4_snr_gt2=float(np.mean(x[sn > 2] > 4)) if np.any(sn > 2) else None, share_snr_gt2=float(np.mean(sn > 2)))
    (HERE/'results/v2_change_decomposition.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
    print('ln error |e| median / rms / bias: const (Huber-D4), oracle slope scale, free slope, window+IMM, tick Kalman; then const parts rms')
    for key, row in out.items():
        if '|' not in key or key.startswith('by_') or key == 'detect':
            continue
        f = lambda e: f"{row[e]['med']:.3f}/{row[e]['rms']:.3f}/{row[e]['bias']:+.3f}"
        print(f"{key:22s} n {row['n']:5d} snr {row['snr_med']:.2f} v {row['vrel_med']:+.2f} | c {f('const')} o {f('oracle')} f {f('free')} imm {f('imm')} tick {f('tick')}")
        print(' '*24+'const parts mean(rms): '+' '.join(f"{p} {row['const']['parts'][p]['mean']:+.3f}({row['const']['parts'][p]['rms']:.3f})" for p in PARTS))
        print(' '*24+'free  parts mean(rms): '+' '.join(f"{p} {row['free']['parts'][p]['mean']:+.3f}({row['free']['parts'][p]['rms']:.3f})" for p in PARTS))
    print('--- by slope SNR |v|/sigma_v: |e| median/rms for const, oracle, free, imm, tick')
    for key, c in out['by_snr'].items():
        print(f"{key:18s} n {c['n']:5d} |v/k| {c['vrel']:.2f} | "+' '.join(f"{e} {c[e]['med']:.3f}/{c[e]['rms']:.3f}/{c[e]['bias']:+.3f}" for e in EST+('imm', 'tick')))
    print('--- trend statistic (v~/s_v)^2 > kappa', json.dumps(out['detect'], indent=0))


if __name__ == '__main__':
    main()
