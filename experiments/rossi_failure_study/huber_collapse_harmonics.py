import os, sys
os.chdir(r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
sys.path.insert(0, os.getcwd()); sys.dont_write_bytecode = True
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_k] = '1'
import numpy as np
from scipy.ndimage import uniform_filter1d
import huber_error_mechanism as E
import huber_prior_check as H
import corrected_comparison as CC

NW = 24
S = H.Solver()
for index in (7, 32, 0):
    case = CC.cases()[index]
    sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, raw = E.closed_loop(case)
    n = len(rows); k = log[:n, 26]; t = log[:n, 0]; F = log[:n, 18]
    m0 = dict(sc, noise_scale=0.); CC.perfect_calibration(m0)
    Y0, z0, k0 = E.replay(m0, pol, sent, case['seed'], n, tau_h, mag_h, eb, policy)
    rows_ok = [j for j in range(NW+1, n) if np.isfinite(Y0[j][-NW:]).all() and k0[j] > 200.]
    stats = []
    for j in rows_ok:
        y0, zz0 = Y0[j][-NW:], z0[j][-NW:]
        f0 = S.fit(y0, zz0)
        Yf, Zf = np.fft.rfft(y0-np.polyval(np.polyfit(np.arange(NW), y0, 3), np.arange(NW))), np.fft.rfft(zz0-np.polyval(np.polyfit(np.arange(NW), zz0, 3), np.arange(NW)))
        # Bins of 5.21 Hz: 4 = 20.8 Hz, 8 = 41.7 Hz, 12 = 62.5 Hz.
        hz = abs(Zf[8])+abs(Zf[12]); hy = (abs(Yf[8])+abs(Yf[12]))*k0[j]
        stats.append((j, f0['k']/k0[j], float(F[max(j-NW, 0):j].min()), float(F[max(j-NW, 0):j].mean()), hz/abs(Zf[4]), hy/(abs(Yf[4])*k0[j]), k0[j]))
    a = np.array(stats)
    bad = a[:, 1] < .5; good = (a[:, 1] > .85) & (a[:, 1] < 1.2)
    print(f'== case {index}: noise-free windows {len(a)}, severe (<.5) {bad.sum()}')
    for name, m in (('severe', bad), ('normal', good)):
        if m.any():
            print(f'  {name:7s} min force in window med {np.median(a[m, 2]):.2f} N, mean force {np.median(a[m, 3]):.2f} | harmonic/fundamental in force {np.median(a[m, 4]):.3f}, in motion*k {np.median(a[m, 5]):.3f} | k {np.median(a[m, 6]):.0f}')
    for row in a[bad][:6]:
        print('   t', round(float(t[int(row[0])]), 3), 'k/k', round(row[1], 3), 'Fmin', round(row[2], 2), 'harm F', round(row[4], 3), 'harm Y', round(row[5], 3))
