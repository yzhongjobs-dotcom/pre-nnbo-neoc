"""Static map of the plant's normal tangent stiffness k_n(x, F): what drives a stiffness change, position or load?

For representative scenes, along the scan path and at several forces around the reference: k_n, the load exponent
p = d ln k / d ln F at fixed position, and how far p is a function of the stiffness level alone.
Run: python -X utf8 stiffness_field_map.py
"""
import os, sys, math
for k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'): os.environ[k] = '1'
sys.dont_write_bytecode = True; sys.path.insert(0, os.getcwd())
import numpy as np
import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
P3 = SE.P3; O.setup()

FS = (.6, .8, 1., 1.25, 1.6)
if __name__ == '__main__':
    cases = CC.cases()
    print('groups:', {g: [i for i, c in enumerate(cases) if c['group'] == g] for g in dict.fromkeys(c['group'] for c in cases)})
    print('case group law | FD | ln k range along x at FD | load exponent p at FD: min median max | corr(p, ln k) | p from the level alone, 2(1-sqrt(kmin/k)): rms misfit | '
          'ln k change for F x0.8 / x1.25: median, max | share of positions with p > .3')
    for idx in (0, 3, 7, 12, 18, 24, 30, 36, 39):
        if idx >= len(cases):
            continue
        case = cases[idx]; sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant']); FD = sc['FD']
        xs = np.linspace(sc['xcmd'][0], sc['xcmd'][-1], 161)
        K = np.empty((len(FS), len(xs)))
        for a, fs in enumerate(FS):
            for b, x in enumerate(xs):
                d = P3.solve_depth(sc, fs*FD, x, 0.); K[a, b] = P3.contact(sc, d, x, 0.)[1]
        lk = np.log(K); p = (lk[3]-lk[1])/math.log(FS[3]/FS[1])
        pred = 2.*(1.-np.sqrt(K[2].min()/K[2]))
        print(f"{idx:2d} {case['group']:12s} {sc['contact_law']:10s} | {FD:4.1f} | {lk[2].max()-lk[2].min():.2f} | {p.min():+.2f} {np.median(p):+.2f} {p.max():+.2f} | "
              f"{np.corrcoef(p, lk[2])[0, 1] if lk[2].std() > 1e-6 else float('nan'):+.2f} | {np.sqrt(np.mean((p-pred)**2)):.2f} | "
              f"{np.median(lk[1]-lk[2]):+.2f}/{np.median(lk[3]-lk[2]):+.2f}, {np.min(lk[1]-lk[2]):+.2f}/{np.max(lk[3]-lk[2]):+.2f} | {np.mean(p > .3):.2f}")
        np.save(os.path.join('results', f'stiffness_field_map_{idx:02d}.npy'), np.vstack([xs, K]))
