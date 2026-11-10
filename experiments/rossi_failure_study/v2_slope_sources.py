"""Two checks on where slope information could come from. Recorded observations only; truth is used for scoring.

1. Excitation coupling. The command is a_force sin / k_good, so the probe force is about a_force k_eff / k_good and the window information
   scales with (k_eff/k_good)^2. Measured: probe-force amplitude of each noise-free window (k_eff |M Y0| at 20.8 Hz) against k_eff / k_hat_prev,
   k_hat_prev the Huber-D4 estimate 12 ticks earlier (the controller holds a lagged estimate).
2. Slow force channel. Quasi-static contact F = f(delta, x): dF/dt = k_t d(delta)/dt + F d ln k_s/dt. With delta = y - u - s(x) - F/k_r,
   r = dF/dt - k_hat dy/dt = F v_s - k (du/dt + s' dx/dt) - ..., so psi = r / k_hat = (F/k) v_s - (slow motion terms).
   Breathing and surface terms are slow (period seconds); rib-driven stiffness changes are near 1-2 Hz. Band-pass psi to [f_lo, 5] Hz and
   read v_hat = psi_bp k_hat / F. Compare with the true slope of ln k_eff (0.2 s smoothed), band-passed the same way.
Predictions fixed before running: probe force falls with k_eff/k_hat_prev; v_hat correlates strongly with the truth outside events, with an
error sd far below the probe-band window's 2.5-4 /s; events contaminate it.
Writes results/v2_slope_sources.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, sosfiltfilt, sosfilt

import huber_prior_check as H
import v2_verify as V

HERE = Path(__file__).resolve().parent
TC, F0 = .008, 20.8
FS = 1/TC


def band(x, lo, hi, causal=False):
    sos = butter(2, [lo, hi], btype='band', fs=FS, output='sos')
    x = np.where(np.isfinite(x), x, 0.)
    return sosfilt(sos, x) if causal else sosfiltfilt(sos, x)


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]; F = log[:, 5]; y = log[:, 2]
    seg = V.segments(t, k, meta['event_info'], k > 200.)
    raw = np.asarray(data['raw'], float)
    kh = np.where(np.isfinite(raw) & (raw > 0), raw, np.nan)
    good = np.isfinite(kh)
    kh = np.interp(np.arange(n), np.flatnonzero(good), kh[good]) if good.any() else np.full(n, np.nan)
    kh_s = uniform_filter1d(kh, 12, mode='nearest')
    # 1. Probe force amplitude per noise-free window.
    S = H.Solver(); amp = np.full(n, np.nan); ratio = np.full(n, np.nan)
    tt = np.arange(24)*TC
    for j in range(40, n):
        Y0 = H.tail(data['Y0'][j], 24)
        if len(Y0) != 24 or k[j] <= 200.:
            continue
        my = S.fit(Y0, H.tail(data['z0'][j], 24), iters=1)['my']
        amp[j] = k[j]*2*abs(np.mean(my*np.exp(-2j*np.pi*F0*tt)))
        ratio[j] = k[j]/kh[j-12]
    # 2. Slow force channel. Derivatives over 5 Hz low-passed signals.
    contact = k > 200.
    Fl = band(F, .05, 5.); yl = band(y, .05, 5.)
    dF = np.gradient(Fl, t); dy = np.gradient(yl, t)
    Fm = uniform_filter1d(F, 25, mode='nearest')
    ell = np.log(np.maximum(k, 1.)); vtrue = np.gradient(uniform_filter1d(ell, 25, mode='nearest'), t)
    out = dict(index=meta['index'], group=meta['group'], seg={q: m.tolist() for q, m in seg.items()}, amp=amp.tolist(), ratio=ratio.tolist(), vtrue=vtrue.tolist())
    for lo in (.3, .5):
        psi = (dF-kh_s*dy)/kh_s
        vb = band(vtrue, lo, 5.)
        vhat = band(psi, lo, 5.)*kh_s/np.maximum(Fm, .5)
        vhat_c = band(psi, lo, 5., causal=True)*kh_s/np.maximum(Fm, .5)
        out[f'vhat_{lo}'] = np.where(contact, vhat, np.nan).tolist(); out[f'vband_{lo}'] = np.where(contact, vb, np.nan).tolist()
        out[f'vhatc_{lo}'] = np.where(contact, vhat_c, np.nan).tolist()
    return out


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rep = {}
    print('--- 1. probe force amplitude relative to its plateau median in the same case, by segment; and by k_eff/k_hat_prev bins (outside events)')
    rows = []
    for r in res:
        a = np.asarray(r['amp']); pl = np.asarray(r['seg']['plateau'])
        base = np.nanmedian(a[pl]) if np.isfinite(a[pl]).sum() > 20 else np.nanmedian(a)
        for s in ('plateau', 'rising', 'falling', 'reversal'):
            m = np.asarray(r['seg'][s]) & np.isfinite(a)
            rows += [(r['group'], s, x/base, q) for x, q in zip(a[m], np.asarray(r['ratio'])[m])]
    rep['probe_amp'] = {}
    for s in ('plateau', 'rising', 'falling', 'reversal'):
        x = np.array([q[2] for q in rows if q[1] == s]); rr = np.array([q[3] for q in rows if q[1] == s])
        rep['probe_amp'][s] = dict(median=float(np.median(x)), ratio_median=float(np.nanmedian(rr)))
        print(f"  {s:9s} amp/plateau {np.median(x):.2f}  k_eff/k_hat_prev {np.nanmedian(rr):.2f}")
    allx = np.array([q[2] for q in rows]); allr = np.array([q[3] for q in rows])
    for lo, hi in ((0, .7), (.7, .9), (.9, 1.1), (1.1, 1.4), (1.4, 9)):
        m = (allr >= lo) & (allr < hi)
        if m.sum() > 20:
            rep['probe_amp'][f'ratio {lo}-{hi}'] = float(np.median(allx[m]))
            print(f"  k_eff/k_hat_prev {lo}-{hi}: amp/plateau {np.median(allx[m]):.2f} (n {m.sum()}), predicted ~ ratio {np.median(allr[m]):.2f}")
    print('--- 2. slow-channel slope vs truth (both band-passed), outside events and inside: corr, regression gain, error sd (1/s), |v_true| sd; zero-phase and causal')
    rep['slow'] = {}
    for lo in (.3, .5):
        for g in ('A', 'A_slow_thin', 'A_large', 'A+B', 'B', 'cough', 'combo', 'all'):
            for scope in ('changes', 'events'):
                for kind in ('vhat', 'vhatc'):
                    xs, ys = [], []
                    for r in res:
                        if g != 'all' and r['group'] != g:
                            continue
                        if scope == 'changes':
                            m = np.zeros(len(r['vtrue']), bool)
                            for s in ('rising', 'falling', 'reversal', 'plateau'):
                                m |= np.asarray(r['seg'][s])
                        else:
                            m = np.asarray(r['seg']['event'])
                        xs.append(np.asarray(r[f'{kind}_{lo}'])[m]); ys.append(np.asarray(r[f'vband_{lo}'])[m])
                    x, yv = np.concatenate(xs), np.concatenate(ys)
                    ok = np.isfinite(x) & np.isfinite(yv)
                    if ok.sum() < 50:
                        continue
                    x, yv = x[ok], yv[ok]
                    c = float(np.corrcoef(x, yv)[0, 1]); gain = float(x@yv/(yv@yv))
                    err = x-yv; q = np.percentile(err, [25, 75])
                    rep['slow'][f'{lo}|{g}|{scope}|{kind}'] = dict(n=int(ok.sum()), corr=c, gain=gain, err_sd=float((q[1]-q[0])/1.349), err_rms=float(np.sqrt(np.mean(err**2))), v_sd=float(np.std(yv)))
                    if kind == 'vhat' or g in ('all', 'A_large'):
                        s = rep['slow'][f'{lo}|{g}|{scope}|{kind}']
                        print(f"  f_lo {lo} {g:12s} {scope:7s} {kind:5s} n {s['n']:5d} corr {s['corr']:+.2f} gain {s['gain']:.2f} err sd {s['err_sd']:.2f} rms {s['err_rms']:.2f} | v sd {s['v_sd']:.2f}")
    (HERE/'results/v2_slope_sources.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
