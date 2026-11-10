import os, sys
os.chdir(r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
sys.path.insert(0, os.getcwd()); sys.dont_write_bytecode = True
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_k] = '1'
import numpy as np
import probe_band_loop as L
import corrected_comparison as CC

for index in (29, 26):
    case = CC.cases()[index]
    traces = {}
    for arm in ('huber', 'band_huber'):
        L.O.setup()
        sc = L.SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
        pol = L.B.build('Huber-D4', L.O.chosen()['Huber-D4'], sc); CC.own_mechanism(pol)
        est = pol.estimator
        valid = []
        if arm != 'huber':
            f = L.BandEstimate(huber=True)
            def wrapped(o, n, p, f=f, est=est):
                out = f(est, o, n, p); valid.append(np.isfinite(out[0])); return out
            est.estimate = wrapped
        kg = []
        step = pol.step
        def cap(o, step=step, pol=pol):
            u = step(o); kg.append((pol.k_good, pol.kbar)); return u
        pol.step = cap
        log = L.P3.episode(sc, pol, case['seed'])
        traces[arm] = (log, np.array(kg), np.array(valid))
    log_h = traces['huber'][0]
    print(f"== case {index} events", [(e['kind'], round(e['start'], 2), round(e['end'], 2)) for e in sc['event_info']])
    for arm, (log, kg, valid) in traces.items():
        t = log[:, 0]; F = log[:, 18]; k = log[:, 26]; n = min(len(t), len(kg))
        pk = int(np.argmax(np.abs(F[:n]-sc['FD'])))
        print(f'  {arm:10s} peak {abs(F[pk]-sc["FD"]):.1f} N at t {t[pk]:.3f}; invalid share {1-valid.mean() if len(valid) else 0:.2f}')
        for j in range(pk-30, min(pk+8, n), 4):
            print(f'     t {t[j]:.3f} F {F[j]:6.2f} k_eff {k[j]:6.0f} k_good/k_eff {kg[j, 0]/k[j]:.2f} k_bar/k_eff {kg[j, 1]/k[j]:.2f}' + (f' valid {bool(valid[j])}' if len(valid) > j else ''))
