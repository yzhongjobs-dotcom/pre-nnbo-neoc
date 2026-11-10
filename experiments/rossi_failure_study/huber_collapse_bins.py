import os, sys
os.chdir(r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
sys.path.insert(0, os.getcwd()); sys.dont_write_bytecode = True
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_k] = '1'
import numpy as np
import huber_error_mechanism as E
import huber_prior_check as H
import corrected_comparison as CC

NW = 24
S = H.Solver()
for index, (t0, t1) in ((7, (6.1, 6.95)), (32, (5.9, 6.4))):
    case = CC.cases()[index]
    sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, raw = E.closed_loop(case)
    n = len(rows); k = log[:n, 26]; t = log[:n, 0]
    Y, z, _ = E.W.windows(E.W.hc_log(rows, tau_h, mag_h, policy['Tc']), policy, policy['enc'], eb)
    print(f'== case {index}: t | raw/k_eff | motion amplitude x k at 20.8 Hz and at 41.7-62.5 Hz (N) | per-bin share of the estimate s_b and window gain Re(G_b)/k for bins 4 (20.8), 8 (41.7), 10 (52.1), 12 (62.5)')
    for j in range(n):
        if not (t0 <= t[j] <= t1) or j % 3:
            continue
        y, zz = Y[j][-NW:], z[j][-NW:]
        if not np.isfinite(y).all():
            print(f'  {t[j]:.3f} no window'); continue
        f = S.fit(y, zz); q = f['my']/f['info']
        Q, Yf, Zf = np.fft.rfft(q), np.fft.rfft(y), np.fft.rfft(zz)
        share = np.real(np.conj(Q)*Yf)*2/NW; share[0] /= 2; share[-1] /= 2      # sum = q'y = 1
        G = np.real(Zf/np.where(np.abs(Yf) > 0, Yf, np.nan))/k[j]
        a4 = k[j]*2*abs(Yf[4])/NW; ah = k[j]*2*np.sqrt(np.sum(np.abs(Yf[8:13])**2))/NW
        print(f'  {t[j]:.3f} | {raw[j]/k[j]:6.3f} | {a4:.3f} {ah:.3f} | '+' '.join(f'b{b}: s {share[b]:+.2f} G {G[b]:+.2f}' for b in (4, 8, 10, 12))+f' | huber {np.mean(f["w"] < 1):.2f}')
