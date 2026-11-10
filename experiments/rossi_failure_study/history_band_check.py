"""Fixed physical probe-band sensitivity of the same adaptive-history rule.

24-tick projector is exactly probe_band_check.PB['band']; longer histories keep
15--30 Hz, not the original bin indices. Original plant, force and raw streams.
No model parameter, threshold or selected scale reads the evaluation truth.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import json,time
import numpy as np
import adaptive_history_v2 as A
import history_propagation_check as P


def projector(n):
    tt=np.arange(n);poly=np.vander((tt-tt.mean())/n,4)
    Qp=np.linalg.qr(poly)[0]
    bins=[b for b in range(1,n//2) if 15<=b/(n*A.TC)<=30]
    X=np.column_stack([f(2*np.pi*b*tt/n) for b in bins for f in (np.cos,np.sin)])
    X-=Qp@(Qp.T@X);U,s,_=np.linalg.svd(X,full_matrices=False)
    U=U[:,s>1e-10*s[0]]
    return U@U.T


def run(index,eps):
    path=A.OUT/f'bandfits_{index:02d}.npz'
    stream=np.load(A.OUT/f'stream_{index:02d}.npz')
    d={key:stream[key] for key in stream.files};d['meta']=json.loads(str(d['meta']))
    ids=np.load(A.OUT/f'fits_{index:02d}.npz')['ids'];nupdate=len(ids)
    if path.exists():
        f=np.load(path);theta,Q,valid=(f[key] for key in ('theta','Q','valid'))
    else:
        theta=np.full((nupdate,10,2),np.nan);Q=np.zeros((nupdate,10,2,96));valid=np.zeros((nupdate,10),bool)
        matrices={n:projector(n) for n in A.SIZES}
        for j,tick in enumerate(ids):
            for i,n in enumerate(A.SIZES):
                h0,h1,z=A.window(d,tick,n);X=np.column_stack([h0,h1/A.TREF]);M=matrices[n]
                for degree,slot in ((1,i),(0,i+5)):
                    xx=X if degree else X[:,:1];gram=xx.T@M@xx
                    if np.linalg.eigvalsh(gram)[0]<=1e-24:continue
                    qq=np.linalg.solve(gram,xx.T@M);th=qq@z
                    Q[j,slot,:degree+1,-n:]=qq;theta[j,slot]=[th[0],th[1] if degree else 0.]
                    valid[j,slot]=200<th[0]<3e5
        np.savez_compressed(path,theta=theta,Q=Q,valid=valid,ids=ids)
    fixed=[];calibrated=[];nfixed=[];ncal=[];degreecal=[];cutoffs=[]
    for j in range(nupdate):
        sc,var=P.statistics(np.nan_to_num(theta[j])[None],Q[j],valid[j])
        # Fixed pre-existing .05 familywise comparison is kept as a sensitivity,
        # not retuned against physical replay errors.
        from scipy.stats import chi2
        cfixed=np.sqrt(chi2.ppf(1-.05/40,2))
        idx=P.chosen(sc,var,cfixed)[0]
        fixed.append(theta[j,idx,0] if idx>=0 else np.nan);nfixed.append(int(P.NS[idx]) if idx>=0 else 0)
        c,ub=P.calibrate_packet(Q[j],valid[j],eps);idx=P.chosen(sc,var,c)[0]
        calibrated.append(theta[j,idx,0] if idx>=0 else np.nan);ncal.append(int(P.NS[idx]) if idx>=0 else 0)
        degreecal.append(int(P.DEG[idx]) if idx>=0 else -1);cutoffs.append(c)
    truth=d['truth'];k=truth[:,26];t=truth[:,0]
    seg=A.V.segments(t,k,d['meta']['event_info'],k>200);masks={name:mask[ids] for name,mask in seg.items()}
    estimates=dict(constant24=theta[:,5,0],constant96=theta[:,9,0],linear24=theta[:,0,0],linear96=theta[:,4,0],fixed=fixed,calibrated=calibrated)
    common=np.logical_and.reduce([np.isfinite(estimates[nm]) for nm in ('constant24','fixed','calibrated')])
    metrics={name:A.score(np.where(common,values,np.nan),k[ids],masks) for name,values in estimates.items()}
    result=dict(index=index,group=d['meta']['group'],metrics=metrics,threshold_median=float(np.median(cutoffs)),coverage=float(common.mean()),
          histories={name:dict(fixed=float(np.median(np.asarray(nfixed)[m])),calibrated=float(np.median(np.asarray(ncal)[m])),linear_fraction=float(np.mean(np.asarray(degreecal)[m]==1))) for name,m in masks.items() if m.sum()},
          trace=dict(t=t[ids].tolist(),truth=k[ids].tolist(),**{k:np.asarray(v).tolist() for k,v in estimates.items()},n=ncal,degree=degreecal,cutoff=cutoffs))
    print(index,{s:{m:round(metrics[m][s]['rmse'],3) for m in estimates if s in metrics[m]} for s in metrics['constant24']},flush=True)
    return result


def main():
    from concurrent.futures import ProcessPoolExecutor
    eps=np.random.default_rng(20261110).normal(size=(4096,96))*A.SIG
    with ProcessPoolExecutor(4) as pool:rows=list(pool.map(run,[0,12,24,36,1,13,25,37],[eps]*8))
    (A.OUT/'band.json').write_text(json.dumps(dict(band_hz=[15,30],cases=rows),allow_nan=True))
if __name__=='__main__':main()
