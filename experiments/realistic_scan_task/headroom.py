"""Step 0 headroom check: is there value in environment-dependent impedance?
(a) best single fixed gains on the cohort, (b) per-scenario best fixed (oracle),
with and without the explicit force-PI outer loop."""
import numpy as np, itertools, json, sys
from multiprocessing import Pool
import plant as P, scenarios as S
GRID=list(itertools.product([1000.,3000.],[10.,100.,400.],[60.,150.,300.],[5.,20.,60.,150.]))
def run(args):
    task,i,p,mode,ctrl=args
    out=[]
    for g in GRID:
        log,_=P.episode(p,lambda t,o,gg,g=g:np.array(g),mode=mode,ctrl=ctrl,seed=i)
        out.append(S.score(log))
    return task,i,mode,out
if __name__=='__main__':
    n=int(sys.argv[1]) if len(sys.argv)>1 else 6
    jobs=[]
    for task in ['T1','T2','T3']:
        for i,p in enumerate(S.cohort(task,n,2026)):
            jobs.append((task,i,p,0,(0.,0.))); jobs.append((task,i,p,1,(1.5,10.)))
    with Pool(4) as pool: res=pool.map(run,jobs)
    summary={}
    for task in ['T1','T2','T3']:
        for mode in [0,1]:
            rows=[r[3] for r in res if r[0]==task and r[2]==mode]
            sc=np.array([[x['score'] for x in row] for row in rows])   # scen x grid
            gbest=int(np.argmin(sc.mean(0))); single=sc[:,gbest]; oracle=sc.min(1)
            arg=[GRID[int(j)] for j in sc.argmin(1)]
            summary[f'{task}_mode{mode}']={'best_single_gains':GRID[gbest],'single_mean':float(single.mean()),'oracle_mean':float(oracle.mean()),
                'oracle_gain_percent':float(100*np.mean((single-oracle)/single)),'per_scene_best':arg,
                'single_components':{k:float(np.mean([row[gbest][k] for row in rows])) for k in rows[0][0]}}
            print(task,'PI' if mode else 'noPI',summary[f'{task}_mode{mode}']['best_single_gains'],'single %.3f oracle %.3f gain %.1f%%'%(single.mean(),oracle.mean(),summary[f'{task}_mode{mode}']['oracle_gain_percent']),arg)
    json.dump(summary,open('headroom_results.json','w'),indent=1)
