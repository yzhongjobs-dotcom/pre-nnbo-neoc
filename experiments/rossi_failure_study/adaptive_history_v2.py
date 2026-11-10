"""Joint endpoint/slope compatibility on unique H-domain observations.

Reproduce: python experiments/rossi_failure_study/adaptive_history_v2.py
Original plant/controller unchanged. Captures the full causal H-processed stream
by instrumenting huber_windows.windows (observation only), then replays estimates.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
import json,math,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.linalg import solveh_banded
from scipy.stats import chi2
import huber_windows as W
import huber_prior_check as H
import corrected_comparison as CC
import v2_verify as V

HERE=Path(__file__).resolve().parent
OUT=HERE/'results/adaptive_history'
SIZES=(24,32,48,64,96)
TC=.008; TREF=24*TC/2; SIG=H.SIGMA


def capture(index):
    OUT.mkdir(exist_ok=True);p=OUT/f'stream_{index:02d}.npz'
    if p.exists():return
    recorded=[];original=W.windows
    def observe(log,policy,enc,eb):
        recorded.append((log.copy(),np.asarray(eb).copy(),dict(policy)))
        return original(log,policy,enc,eb)
    W.OUT.mkdir(exist_ok=True);W.windows=observe
    try:meta=W.job(CC.cases()[index])
    finally:W.windows=original
    assert len(recorded)==3
    log,eb,policy=recorded[0];log0=recorded[1][0][:len(log)]
    data=H.load(H.WIN/f'case_{index:02d}.npz');truth,_=H.masks_for(data)
    np.savez_compressed(p,log=log,log0=log0,eb=eb,truth=truth[:len(log)],meta=np.asarray(json.dumps(meta)))


class Fitter:
    def __init__(self):self.solver=H.Solver()
    def fit(self,h0,h1,z,linear=True):
        n=len(z);D,basis=self.solver.bands(n)
        X=np.column_stack([h0,h1/TREF]) if linear else h0[:,None]
        rhs=np.column_stack([X,z]);w=np.ones(n-4);converged=False
        for iteration in range(80):
            ab=H.RHO*(basis@w).reshape(5,n);ab[0]+=1
            smooth=solveh_banded(ab,rhs,lower=True,check_finite=False)
            prof=rhs-smooth
            A=X.T@prof[:,:X.shape[1]];A=(A+A.T)/2
            if np.linalg.eigvalsh(A)[0]<=1e-24:return None
            Q=np.linalg.solve(A,prof[:,:X.shape[1]].T)
            theta=Q@z
            d=smooth[:,-1]-smooth[:,:X.shape[1]]@theta
            u=D@d/(.04/math.sqrt(H.RHO))
            wn=np.minimum(1.,1.345/np.maximum(abs(u),1e-15))
            if max(abs(wn-w))<1e-6:converged=True;break
            w=wn
        # The IRLS point weight Q is NOT the derivative of a Huber solution.
        # The local Hessian uses psi'(u)=1 inside the quadratic region, 0 outside.
        active=(abs(u)<1.345).astype(float)
        abj=H.RHO*(basis@active).reshape(5,n);abj[0]+=1
        px=X-solveh_banded(abj,X,lower=True,check_finite=False)
        Aj=X.T@px;Aj=(Aj+Aj.T)/2
        ej=np.linalg.eigvalsh(Aj)
        jacvalid=ej[0]>max(ej[-1]*1e-10,1e-24)
        J=np.linalg.solve(Aj,px.T) if jacvalid else np.full_like(Q,np.nan)
        if not linear:
            theta=np.r_[theta[0],0.];Q=np.vstack([Q,np.zeros(n)]);J=np.vstack([J,np.zeros(n)])
        return dict(theta=theta,Q=Q,C=SIG**2*(Q@Q.T),J=J,Cjac=SIG**2*(J@J.T),jacvalid=jacvalid,
                    n=n,converged=converged,valid=converged and 200.<theta[0]<3e5,
                    active=float(np.mean(w<1.)),kink_margin=float(np.min(abs(abs(u)-1.345))))


def window(stream,j,n,zero=False):
    log=stream['log0'] if zero else stream['log'];eb=stream['eb'][j]
    meta=stream['meta'];p=meta['policy'];j0=j-n
    if j0<max(p['ndel']+1,0):return None
    _,h0,z,keep=W.common._build(log,j0,n,0,p['ns'],eb,p['afc'],p['ndel'],log[j,0],log[j,2],n*TC,True)
    if not keep.all():return None
    h1=W.time_weighted(log,j0,n,p['ns'],eb,p['afc'],log[j,0],log[j,2],True)
    return h0,h1,z


def pairwise(fits):
    out={};nmax=max(f['n'] for f in fits)
    for b in range(len(fits)):
        for a in range(b):
            fa,fb=fits[a],fits[b]
            Qa=np.pad(fa['Q'],((0,0),(nmax-fa['n'],0)))
            Qb=np.pad(fb['Q'],((0,0),(nmax-fb['n'],0)))
            dif=Qb-Qa;C=SIG**2*(dif@dif.T)
            d=fb['theta']-fa['theta'];e,U=np.linalg.eigh(C)
            use=e>max(e[-1]*1e-10,1e-20)
            stat=float(np.sum((U[:,use].T@d)**2/e[use]))
            # A deterministic difference in a zero-noise direction is incompatible.
            if np.linalg.norm(U[:,~use].T@d)>1e-8*max(np.linalg.norm(d),1.):stat=np.inf
            sk=float(d[0]**2/max(C[0,0],1e-20))
            out[a,b]=(stat,int(use.sum()),sk,C)
    return out


def select(fits,alpha=.05,joint=True):
    pairs=pairwise(fits);m=max(len(pairs),1);acceptable=[]
    for b,f in enumerate(fits):
        if not f['valid']:continue
        okay=True
        for a in range(b):
            if not fits[a]['valid']:continue
            stat,rank,sk,C=pairs[a,b]
            threshold=chi2.ppf(1-alpha/m,max(rank,1) if joint else 1)
            if (stat if joint else sk)>threshold:okay=False;break
        if okay:acceptable.append(b)
    if not acceptable:return None,pairs
    best=min(acceptable,key=lambda b:(fits[b]['C'][0,0],fits[b]['n']))
    return fits[best],pairs


def select_models(references,candidates,alpha=.05):
    """Compare each model against unshrunk local-linear references, including its own scale."""
    m=max(sum(r['n']<=c['n'] for c in candidates for r in references),1)
    eligible=[]
    for candidate in candidates:
        if not candidate['valid']:continue
        okay=True
        for reference in references:
            if not reference['valid'] or reference['n']>candidate['n']:continue
            nmax=max(candidate['n'],reference['n'])
            qc=np.pad(candidate['Q'],((0,0),(nmax-candidate['n'],0)))
            qr=np.pad(reference['Q'],((0,0),(nmax-reference['n'],0)))
            diff=qc-qr;C=SIG**2*(diff@diff.T)
            delta=candidate['theta']-reference['theta']
            ev,U=np.linalg.eigh(C);use=ev>max(ev[-1]*1e-10,1e-20)
            rank=int(use.sum())
            stat=float(np.sum((U[:,use].T@delta)**2/ev[use]))
            if np.linalg.norm(U[:,~use].T@delta)>1e-8*max(np.linalg.norm(delta),1.):stat=np.inf
            if stat>chi2.ppf(1-alpha/m,max(rank,1)):okay=False;break
        if okay:eligible.append(candidate)
    if not eligible:return None
    return min(eligible,key=lambda c:(c['C'][0,0],c['n']))


def score(vals,truth,seg):
    ans={}
    for key,mask in seg.items():
        vv=np.asarray(vals);use=mask&np.isfinite(vv)&(truth>200)
        if use.sum()<4:continue
        err=np.log(np.maximum(vv[use],200)/truth[use])
        ans[key]=dict(n=int(use.sum()),rmse=float(np.sqrt(np.mean(err**2))),
            median_abs=float(np.median(abs(err))),bias=float(np.median(err)))
    return ans


def process(index):
    capture(index)
    with np.load(OUT/f'stream_{index:02d}.npz') as ff:d={k:ff[k] for k in ff.files}
    d['meta']=json.loads(str(d['meta']));truth=d['truth'];t=truth[:,0];k=truth[:,26]
    seg=V.segments(t,k,d['meta']['event_info'],k>200)
    ids=np.arange(102,len(k),6);fitter=Fitter()
    trace={name:[] for name in ['model','chosen_model','degree_model','v_model','sdv_model','constant96','constant24','linear24','linear96','joint','scalar','chosen_joint','chosen_scalar','v_joint','sdv_joint','v24','sdv24','noise_free_linear24','noise_free_linear96','noise_free_adaptive']}
    events=[];valid_counts=[];packets=[];ctime=time.perf_counter()
    for j in ids:
        fits=[];constants=[]
        for n in SIZES:
            inp=window(d,j,n)
            f=fitter.fit(*inp) if inp is not None else None
            if f is not None:
                f['degree']=1;fits.append(f)
            constant=fitter.fit(*inp,linear=False) if inp is not None else None
            if constant is not None:
                constant['degree']=0;constants.append(constant)
        packet=[]
        for f in fits+constants:
            packet.append(dict(n=f['n'],degree=f['degree'],theta=f['theta'],Q=np.pad(f['Q'],((0,0),(max(SIZES)-f['n'],0))),valid=f['valid']))
        packets.append(packet)
        model=select_models(fits,fits+constants)
        selected,pairs=select(fits,joint=True) if fits else (None,{})
        scalar,_=select(fits,joint=False) if fits else (None,{})
        by_n={f['n']:f for f in fits}
        const=next((c for c in constants if c['n']==24),None)
        const96=next((c for c in constants if c['n']==96),None)
        chosen=selected['n'] if selected else 0
        q=by_n.get(24);long=by_n.get(96)
        values=dict(model=model['theta'][0] if model else np.nan,
                    chosen_model=model['n'] if model else 0,degree_model=model['degree'] if model else -1,
                    v_model=model['theta'][1]/TREF if model else np.nan,
                    sdv_model=np.sqrt(model['C'][1,1])/TREF if model else np.nan,
                    constant96=const96['theta'][0] if const96 else np.nan,constant24=const['theta'][0] if const else np.nan,
                    linear24=q['theta'][0] if q else np.nan,linear96=long['theta'][0] if long else np.nan,
                    joint=selected['theta'][0] if selected else np.nan,scalar=scalar['theta'][0] if scalar else np.nan,
                    chosen_joint=chosen,chosen_scalar=scalar['n'] if scalar else 0,
                    v_joint=selected['theta'][1]/TREF if selected else np.nan,
                    sdv_joint=np.sqrt(selected['C'][1,1])/TREF if selected else np.nan,
                    v24=q['theta'][1]/TREF if q else np.nan,sdv24=np.sqrt(q['C'][1,1])/TREF if q else np.nan)
        for n,name in [(24,'noise_free_linear24'),(96,'noise_free_linear96'),(chosen,'noise_free_adaptive')]:
            inp=window(d,j,n,zero=True) if n else None
            if inp is not None:
                ff=fitter.fit(*inp);values[name]=ff['theta'][0] if ff else np.nan
            else:values[name]=np.nan
        for name in trace:trace[name].append(values[name])
        valid_counts.append(sum(f['valid'] for f in fits))
    kk=k[ids];sg={name:m[ids] for name,m in seg.items()}
    metrics={name:score(trace[name],kk,sg) for name in ['constant24','constant96','linear24','linear96','joint','model','scalar','noise_free_linear24','noise_free_linear96','noise_free_adaptive']}
    info={}
    for name,m in sg.items():
        if m.sum()<4:continue
        sizes=np.asarray(trace['chosen_joint'])[m]
        ratio=np.asarray(trace['sdv_joint'])[m]/np.asarray(trace['sdv24'])[m]
        info[name]=dict(model_history_median=float(np.median(np.asarray(trace['chosen_model'])[m])),
                       model_linear_fraction=float(np.mean(np.asarray(trace['degree_model'])[m]==1)),n=int(m.sum()),chosen_median=float(np.median(sizes)),fraction96=float(np.mean(sizes==96)),
                       fraction24=float(np.mean(sizes==24)),slope_noise_ratio=float(np.nanmedian(ratio)))
    theta_cache=np.full((len(ids),10,2),np.nan);Q_cache=np.full((len(ids),10,2,max(SIZES)),np.nan);valid_cache=np.zeros((len(ids),10),bool)
    for ii,packet in enumerate(packets):
        for f in packet:
            slot=(0 if f['degree']==1 else 5)+SIZES.index(f['n'])
            theta_cache[ii,slot]=f['theta'];Q_cache[ii,slot]=f['Q'];valid_cache[ii,slot]=f['valid']
    np.savez_compressed(OUT/f'fits_{index:02d}.npz',theta=theta_cache,Q=Q_cache,valid=valid_cache,ids=ids)
    out=dict(index=index,group=d['meta']['group'],threshold_alpha=.05,metrics=metrics,information=info,availability={name:{segname:float(np.mean(np.isfinite(np.asarray(trace[name]))[mask])) for segname,mask in sg.items() if mask.sum()} for name in ['constant24','linear24','linear96','joint','model']},
             runtime_seconds=time.perf_counter()-ctime,valid_candidates_median=float(np.median(valid_counts)),
             trace=dict(t=t[ids].tolist(),k=kk.tolist(),**{q:np.asarray(v).tolist() for q,v in trace.items()}))
    (OUT/f'case_{index:02d}.json').write_text(json.dumps(out,allow_nan=True))
    print(index,d['meta']['group'],{s:{q:round(metrics[q][s]['rmse'],3) for q in ['constant24','constant96','linear24','linear96','joint','model'] if s in metrics[q]} for s in sg if s in metrics['joint']},flush=True)
    return out


def main():
    cases=[0,12,24,36,1,13,25,37]
    with ProcessPoolExecutor(4) as pool:rows=list(pool.map(process,cases))
    slim=[{k:v for k,v in r.items() if k!='trace'} for r in rows]
    (OUT/'summary.json').write_text(json.dumps(slim,indent=1,allow_nan=True))
if __name__=='__main__':main()
