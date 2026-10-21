"""Shared evaluation and tuning utilities (same budget and score for every method)."""
import itertools, numpy as np
from multiprocessing import Pool
import plant as P, scenarios as S

GRID_M=[0.3,1.,3.]; GRID_D=list(np.round(np.geomspace(20,2500,12),1)); GRID_K=[0.,300.]
GRID_KI=[0.,3.,10.,30.]; GRID_KPF=[0.,2e-5,5e-5,1e-4,2e-4,4e-4]
FIXED_GRID=list(itertools.product(GRID_M,GRID_D,GRID_K,GRID_KI,GRID_KPF))

def _eval(args):
    p,T,ctrl,plat,seed,update=args
    log=P.episode(p,ctrl,T=T,platform=plat,seed=seed,update=update)
    return S.metrics(log,Tc=P.PLATFORMS[plat]['Tc'] if isinstance(plat,str) else plat['Tc'])

_pool=None
def pool():
    global _pool
    if _pool is None: _pool=Pool(4)
    return _pool

def evaluate(scenes,ctrls,plat='PC125',seed0=0,update=None):
    """metrics[scene][ctrl] for every scene x controller."""
    jobs=[(p,T,c,plat,seed0+i,update) for i,(p,T) in enumerate(scenes) for c in ctrls]
    res=pool().map(_eval,jobs,chunksize=max(1,len(jobs)//64))
    nc=len(ctrls); return [res[i*nc:(i+1)*nc] for i in range(len(scenes))]

def J(res): return np.array([[m['J'] for m in row] for row in res])

def feasible(res,peak=P.FD+3.0,frac=0.95):
    pk=np.array([[m['peak'] if np.isfinite(m['peak']) else 1e9 for m in row] for row in res])
    return (pk<=peak).mean(0)>=frac

def best_fixed(train,plat='PC125',seed0=0,refine=200):
    """B1: grid search then Nelder-Mead in log space from the 5 best grid points; selection by
    mean training J among settings with peak <= FD+3 N in >= 95 % of training scenes (if none is
    feasible, plain minimum J is used and reported)."""
    res=evaluate(train,FIXED_GRID,plat,seed0); Jm=J(res).mean(0); feas=feasible(res)
    order=np.argsort(np.where(feas,Jm,Jm+1e6))
    lo=np.log(P.GMIN[:5]+np.array([0,0,1,0.1,1e-6])); hi=np.log(P.GMAX[:5]+np.array([0,0,1,0.1,1e-6]))
    def obj(z):
        g=np.exp(np.clip(z,lo,hi))-np.array([0,0,1,0.1,1e-6])
        r=evaluate(train,[tuple(g)],plat,seed0); j=float(np.mean([row[0]['J'] for row in r]))
        f=feasible(r)[0]
        return j+(0 if f else 1e3), g
    best=None
    from scipy.optimize import minimize
    for j in order[:5]:
        g0=np.array(FIXED_GRID[j]); z0=np.log(g0+np.array([0,0,1,0.1,1e-6]))
        cache={}
        def f(z):
            v,g=obj(z); cache[tuple(np.round(z,6))]=(v,g); return v
        r=minimize(f,z0,method='Nelder-Mead',options=dict(maxfev=refine,xatol=1e-3,fatol=1e-4,initial_simplex=None))
        v,g=obj(r.x)
        if best is None or v<best[0]: best=(v,tuple(float(x) for x in g))
    gfeas=bool(feas.any())
    return dict(gains=best[1],train_J=best[0],grid_best=FIXED_GRID[int(order[0])],any_feasible_on_grid=gfeas,
                grid_res=res)

def summary(rows):
    """cohort means of physical metrics for a list of metric dicts."""
    keys=['J','rmse','p95','maxerr','inband','peak','t_over','loss','loss_max','osc','shear','diverged']
    out={k:float(np.nanmean([r[k] for r in rows])) for k in keys}
    pk=np.array([r['peak'] for r in rows]); out['frac_peak_over']=float(np.mean(~(pk<=P.FD+3.0)))
    return out

def paired_ci(a,b,n=10000,seed=0):
    """mean(a-b) with a bootstrap 95 % CI over scenes."""
    d=np.asarray(a)-np.asarray(b); rng=np.random.default_rng(seed)
    bs=d[rng.integers(0,len(d),(n,len(d)))].mean(1)
    return float(d.mean()),float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))
