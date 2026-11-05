import sys, pickle
pass
from concurrent.futures import ProcessPoolExecutor
import numpy as np, d_compare as D, method3 as M, physics_plant3 as P3, scenes3 as S3
cfg = pickle.loads(open(sys.argv[1],'rb').read())['chosen']['D4']
def job(a):
    dyn, ph, i, mode = a
    sc = S3.cell_scene3('cough', 'relaxed', 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=ph))
    if mode == 'no activation': sc['act'] = np.zeros_like(sc['act'])
    if mode == 'no motion': sc['u'] = np.zeros_like(sc['u'])
    pol = D.make('D4', cfg, sc); log = P3.episode(sc, pol, 900+i); m = M.metric(log, pol, sc)
    return mode, m['rmse'], m['event_rmse'], m['peak']
if __name__ == '__main__':
    jobs = [(dyn, ph, i, mode) for mode in ('as is', 'no activation', 'no motion') for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    with ProcessPoolExecutor(15) as ex: r = list(ex.map(job, jobs))
    for mode in ('as is', 'no activation', 'no motion'):
        q = [x for x in r if x[0] == mode]
        print(f"cough, {mode:14s} rmse median {np.median([x[1] for x in q]):.3f} event median {np.nanmedian([x[2] for x in q]):.3f} peak max {max(x[3] for x in q):.1f}")
