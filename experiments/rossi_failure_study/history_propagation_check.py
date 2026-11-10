"""Noise-risk calibration of joint history/model compatibility. No truth in selection.

The Gaussian null simulations calibrate a threshold, not validate a tissue model.
The selected rule is then evaluated on the unchanged physical-plant recordings.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import json
from pathlib import Path
import numpy as np
import adaptive_history_v2 as A
import v2_verify as V

OUT=A.OUT;NS=np.array(list(A.SIZES)*2);DEG=np.r_[np.ones(5,int),np.zeros(5,int)]
GRID=np.arange(.5,4.01,.25)

def statistics(theta,Q,valid):
    variance=A.SIG**2*np.sum(Q[:,0]**2,axis=1)
    variance[~valid]=np.inf
    # theta: [sample, candidate, scaled parameter]
    scores=np.zeros((len(theta),10))
    for c in range(10):
        if not valid[c]:scores[:,c]=np.inf;continue
        for r in range(5):
            if not valid[r] or NS[r]>NS[c] or r==c:continue
            D=Q[c]-Q[r];C=A.SIG**2*(D@D.T)
            ev,U=np.linalg.eigh(C);use=ev>max(ev[-1]*1e-10,1e-20)
            delta=theta[:,c]-theta[:,r]
            stat=np.sum((delta@U[:,use])**2/ev[use],axis=1)
            if (~use).any():
                residual=np.linalg.norm(delta@U[:,~use],axis=1)
                stat[residual>1e-8*np.maximum(np.linalg.norm(delta,axis=1),1.)]=np.inf
            scores[:,c]=np.maximum(scores[:,c],stat)
    return scores,variance

def chosen(scores,var,c):
    order=np.argsort(var,kind='stable');ordered=np.isfinite(scores[:,order])&(scores[:,order]<=c*c)
    has=ordered.any(axis=1);idx=order[np.argmax(ordered,axis=1)]
    idx[~has]=-1
    return idx

def calibrate_packet(Q,valid,eps,budget=.05):
    """Smallest cutoff whose paired conditional null-risk UCB is within 5% of constant24."""
    null=np.einsum('bn,cpn->bcp',eps,np.nan_to_num(Q))
    scores,var=statistics(null,np.nan_to_num(Q),valid)
    if not valid[5] or var[5]<=0:return 4.,np.nan
    reference=null[:,5,0]**2
    last=np.inf
    for c in GRID:
        ids=chosen(scores,var,c)
        if (ids<0).any():continue
        excess=(null[np.arange(len(null)),ids,0]**2-reference)/var[5]
        upper=float(excess.mean()+3*excess.std(ddof=1)/np.sqrt(len(null)))
        last=upper
        if upper<=budget:return float(c),upper
    ids=chosen(scores,var,np.inf)
    excess=(null[np.arange(len(null)),ids,0]**2-reference)/var[5]
    upper=float(excess.mean()+3*excess.std(ddof=1)/np.sqrt(len(null)))
    return float('inf'),upper # report the actual Monte Carlo bound, including if it misses the budget


def run(index,eps):
    f=np.load(OUT/f'fits_{index:02d}.npz');theta=f['theta'];Q=f['Q'];valid=f['valid'];ids=f['ids']
    src=json.loads((OUT/f'case_{index:02d}.json').read_text())
    stream=np.load(OUT/f'stream_{index:02d}.npz');meta=json.loads(str(stream['meta']))
    truth=stream['truth'];k=truth[:,26];t=truth[:,0];seg=V.segments(t,k,meta['event_info'],k>200)
    estimates=[];nsel=[];psel=[];cutoffs=[];risk_upper=[]
    for j in range(len(theta)):
        c,ub=calibrate_packet(Q[j],valid[j],eps)
        cutoffs.append(c);risk_upper.append(ub)
        sc,var=statistics(np.nan_to_num(theta[j])[None],np.nan_to_num(Q[j]),valid[j])
        idx=chosen(sc,var,c)[0]
        estimates.append(theta[j,idx,0] if idx>=0 else np.nan)
        nsel.append(int(NS[idx]) if idx>=0 else 0);psel.append(int(DEG[idx]) if idx>=0 else -1)
    masks={name:m[ids] for name,m in seg.items()}
    metrics=A.score(estimates,k[ids],masks)
    common=np.isfinite(estimates)&np.isfinite(src['trace']['constant24'])
    paired={name:A.score(np.where(common,values,np.nan),k[ids],masks) for name,values in
            [('baseline',src['trace']['constant24']),('adaptive',estimates)]}
    report=dict(index=index,group=meta['group'],threshold_median=float(np.median(cutoffs)),metrics=metrics,paired=paired,
                history={name:dict(median=float(np.median(np.asarray(nsel)[m])),linear_fraction=float(np.mean(np.asarray(psel)[m]==1)),availability=float(np.mean(np.isfinite(estimates)[m]))) for name,m in masks.items() if m.sum()},
                trace=dict(k=estimates,n=nsel,degree=psel,t=t[ids].tolist(),cutoff=cutoffs,null_risk_upper=risk_upper))
    return report

def main():
    from concurrent.futures import ProcessPoolExecutor
    rng=np.random.default_rng(20261110)
    eps=rng.normal(size=(4096,max(A.SIZES)))*A.SIG
    cases=[0,12,24,36,1,13,25,37]
    with ProcessPoolExecutor(4) as pool:rows=list(pool.map(run,cases,[eps]*len(cases)))
    for r in rows:print(r['index'],'c',round(r['threshold_median'],2),{s:{m:round(v[s]['rmse'],3) for m,v in r['paired'].items()} for s in r['metrics']},flush=True)
    cal=dict(mode='per-design conditional null-risk calibration',samples=4096,seed=20261110,             budget_excess_risk=.05,grid=GRID.tolist(),baseline='same-observation constant24',             truth_used=False,shared_random_numbers=True)
    (OUT/'propagation.json').write_text(json.dumps(dict(calibration=cal,cases=rows),indent=1,allow_nan=True))
if __name__=='__main__':main()
