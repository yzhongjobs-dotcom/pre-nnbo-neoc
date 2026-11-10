import os, sys
os.chdir(r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
sys.path.insert(0, os.getcwd()); sys.dont_write_bytecode = True
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_k] = '1'
import numpy as np
import corrected_comparison as CC
import memory_segment_eval as SE
import overlay_compare as O

P3, B = SE.P3, O.B
NW = 24
for index, (t0, t1) in ((7, (6.0, 6.9)), (0, (3.0, 3.5))):
    case = CC.cases()[index]
    O.setup(); sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc); CC.own_mechanism(pol)
    probe, kg, ff = [], [], []
    step = pol.step

    def cap(o):
        u = step(o)
        # The probe component the policy added this tick, the stiffness it used, and the motion feed-forward.
        probe.append(float(getattr(pol, 'u_exc', np.nan)) if hasattr(pol, 'u_exc') else np.nan)
        kg.append(pol.k_good)
        return u
    pol.step = cap
    log = P3.episode(sc, pol, case['seed'])
    t = log[:, 0]; cmd = log[:, 11]; y = log[:, 2]; F = log[:, 18]; k = log[:, 26]
    n = min(len(t), len(kg))
    print(f'== case {index}; probe attribute present: {np.isfinite(probe).any()}')
    print('  t | amplitude at 20.8 Hz over the last 24 ticks: command, encoder (mm) | encoder/command | k_good/k_eff | force amplitude at 20.8 Hz (N)')
    w = np.exp(-2j*np.pi*20.8*np.arange(NW)*.008)
    for j in range(NW, n, 6):
        if not (t0 <= t[j] <= t1):
            continue
        sl = slice(j-NW, j)
        dc = lambda x: x-np.polyval(np.polyfit(np.arange(NW), x, 3), np.arange(NW))
        ac, ay, af = (2*abs(np.mean(dc(v[sl])*w)) for v in (cmd, y, F))
        print(f'  {t[j]:.3f} | {ac*1e3:.4f} {ay*1e3:.4f} | {ay/max(ac, 1e-12):.2f} | {kg[j]/k[j]:.2f} | {af:.3f}')
