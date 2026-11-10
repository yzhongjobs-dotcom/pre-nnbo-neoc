"""Correct Huber differential covariance: Hessian curvature, not IRLS surrogate.

Finite differences verify the local derivative. This is a local approximation
for finite sensor noise and data-dependent closed-loop regressors, not exact
Gaussian coverage. Compares corrected local-covariance history selection using
unchanged original physical-plant records.
"""
import json
import numpy as np
import adaptive_history_v2 as A


def run(index):
    stream=np.load(A.OUT/f'stream_{index:02d}.npz');d={key:stream[key] for key in stream.files};d['meta']=json.loads(str(d['meta']))
    ids=np.load(A.OUT/f'fits_{index:02d}.npz')['ids'];fitter=A.Fitter()
    est=[];ns=[];degree=[];ratios=[];sd24=[];sd96=[];fd=[];gaps=[]
    base=json.loads((A.OUT/f'case_{index:02d}.json').read_text())
    for tick in ids:
        refs=[];cands=[];lookup={}
        for n in A.SIZES:
            inp=A.window(d,tick,n)
            for linear in (True,False):
                f=fitter.fit(*inp,linear=linear)
                if f is None:continue
                f['degree']=int(linear)
                # Point theta stays the original optimizer; only its local sensitivity changes.
                f['Qpoint']=f['Q'];f['Q']=f['J'];f['C']=f['Cjac']
                f['valid']=f['valid'] and f['jacvalid']
                if f['valid']:
                    oldsd=A.SIG*np.linalg.norm(f['Qpoint'][0]);ratios.append(np.sqrt(f['C'][0,0])/oldsd)
                    gaps.append(f['kink_margin'])
                    if linear and n==24:sd24.append(np.sqrt(f['C'][1,1])/A.TREF)
                    if linear and n==96:sd96.append(np.sqrt(f['C'][1,1])/A.TREF)
                if linear:refs.append(f)
                cands.append(f);lookup[linear,n]=f
        sel=A.select_models(refs,cands)
        est.append(sel['theta'][0] if sel else np.nan);ns.append(sel['n'] if sel else 0);degree.append(sel['degree'] if sel else -1)
        if len(fd)<12:
            f=lookup.get((True,48))
            if f and f['valid']:
                h0,h1,z=A.window(d,tick,48)
                rng=np.random.default_rng(100+tick);direction=rng.normal(size=len(z));direction/=np.linalg.norm(direction)
                step=1e-5
                plus=fitter.fit(h0,h1,z+step*direction);minus=fitter.fit(h0,h1,z-step*direction)
                observed=(plus['theta']-minus['theta'])/(2*step);pred=f['J']@direction
                fd.append(float(np.linalg.norm(observed-pred)/max(np.linalg.norm(pred),1e-9)))
    truth=d['truth'];k=truth[:,26];t=truth[:,0];seg=A.V.segments(t,k,d['meta']['event_info'],k>200);masks={key:m[ids] for key,m in seg.items()}
    common=np.isfinite(est)&np.isfinite(base['trace']['constant24'])
    paired={key:A.score(np.where(common,v,np.nan),k[ids],masks) for key,v in [('constant24',base['trace']['constant24']),('adaptive',est)]}
    result=dict(index=index,metrics=paired,coverage=float(common.mean()),
                jacobian_sd_over_IRLS_sd_median=float(np.median(ratios)),finite_difference_relative_errors=fd,
                median_slope_sd24=float(np.median(sd24)),median_slope_sd96=float(np.median(sd96)),
                trace=dict(t=t[ids].tolist(),k=est,n=ns,degree=degree))
    print(index,'fd',np.median(fd),'sdRatio',result['jacobian_sd_over_IRLS_sd_median'],
          {s:{m:round(v[s]['rmse'],3) for m,v in paired.items()} for s in paired['constant24']},flush=True)
    return result


def main():
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(4) as p:rows=list(p.map(run,[0,12,24,36,1,13,25,37]))
    (A.OUT/'jacobian.json').write_text(json.dumps(rows,indent=1,allow_nan=True))
if __name__=='__main__':main()
