"""Lateral dither phase-1; new files only, unchanged controller/estimator.
Run with python -X utf8. --analyze-only never runs an episode.
L1-L5 in PREDICTIONS are fixed before episodes and saved before execution.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys
sys.dont_write_bytecode=True
import argparse, hashlib, json, tempfile, time, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
PROJECT=Path(r'E:\github_project\scan-impedance-update')
ROOT=Path(__file__).resolve().parents[2]
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'numba_cache'))
os.environ['TRAJECTORY_LOOP_PROJECT']=str(PROJECT)
os.environ['TRAJECTORY_LOOP_CACHE']=os.environ['NUMBA_CACHE_DIR']
sys.path.insert(0,str(PROJECT/'experiments/rossi_failure_study'))
import trajectory_loop as T
import huber_windows as HW
from d_robust import common
import numpy as np
CASES=(36,39,38,0,12)
AMPS=(0.,.05,.1,.2)
TC,FX,FP,SIGMA=.008,15.625,125/6,.0245
PREDICTIONS={
'L1':'36,39 at 0.1 mm: best-lag correlation >=0.8; gain 0.3..1.5',
'L2':'36,39 at 0.1 mm: rate noise <=1/s and <=V2 noise/3',
'L3':'All at 0.1 mm: added changing/plateau ripple <=0.08/0.02 N rms; delta L<=0.03',
'L4':'All at 0.1 mm: median abs log stiffness error <=no-dither k_only error+0.01',
'L5':'12 at 0.1 mm: event mean differs from outside mean by <=2 plateau sd'}
SKILL=Path('C:/Users/user/.agents/skills/nature-figure')

def dump_new(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:
        json.dump(T.CC.clean(obj),f,indent=2,allow_nan=False)

def run(task):
    index,amp,folder=task
    started=time.perf_counter()
    case=T.CC.cases()[index]
    T.O.setup()
    sc=T.SE.MC.scene_for(tuple(case['cell']),T.CC.BASE,case['variant'])
    grid=np.arange(len(sc['xcmd']))*.001
    q=np.clip(grid/.3,0.,1.)
    ramp=q**3*(10.-15.*q+6.*q*q)
    sc['xcmd']=sc['xcmd']+amp*1e-3*np.sin(2*np.pi*FX*grid)*ramp
    pol=T.B.build('Huber-D4',T.O.chosen()['Huber-D4'],sc)
    T.CC.own_mechanism(pol)
    est=pol.estimator
    inner,inner_log=est.estimate,est._log
    seen,pre,ebs,raw,calls=[],[],{},{},[]
    policy={}
    def record_log(o,n):
        first=est.count
        incoming=[o.log[j].copy() for j in range(first,n)]
        result=inner_log(o,n)
        seen.extend(incoming)
        pre.extend([result[j].copy() for j in range(first,n)])
        return result
    def record(o,n,p):
        policy.update(ndel=p.ndel,ns=p.ns,afc=p.afc,enc=bool(p.params['enc_y']),
                      Tc=p.Tc,noise=float(p.params['noise']),nw=est.nw)
        ebs[n-1]=float(p.eb)
        out=inner(o,n,p)
        raw[n-1]=float(out[0]); calls.append(n-1)
        return out
    est._log,est.estimate=record_log,record
    log=T.P3.episode(sc,pol,case['seed'])
    metric=T.M.metric(log,pol,sc)
    rows,pre=np.asarray(seen),np.asarray(pre)
    assert len(rows)==est.count
    assert np.array_equal(pre[:,2],est.private[:len(pre),2])
    eb=np.array([ebs.get(j,np.nan) for j in range(len(rows))])
    meta=dict(index=index,amp_mm=amp,seed=case['seed'],group=case['group'],
              event_info=sc['event_info'],policy=policy,L=float(np.log1p(metric['J'])),
              metrics=metric,v_scan=float(sc['v']),tau_scan_calibrated=float(pol.tau),tau_scan_true=float(sc['tau_a_true']),
              K_t=sc['K_t'],m_t=sc['m_t'],zeta_t=sc['zeta_t'],
              command_velocity_ratio=2*np.pi*FX*amp*1e-3/sc['v'],
              command_min_velocity_after_1p5=float(np.min(np.gradient(sc['xcmd'],.001)[grid>1.5])),
              arm_zeta=sc.get('arm_zeta',-1.),runtime_s=time.perf_counter()-started)
    dest=Path(folder)/f'case_{index:02d}_a{amp:.2f}.npz'
    with dest.open('xb') as f:
        np.savez_compressed(f,log=log,seen=rows,preprocessed=pre,eb=eb,calls=calls,
            raw=np.array([raw.get(j,np.nan) for j in range(len(rows))]),
            meta=np.asarray(json.dumps(T.CC.clean(meta))))
    print(f'Completed case {index}, a={amp:.2f} mm, L={meta["L"]:.6f}',flush=True)
    return dest.name

def avg(x,n):
    # Local sums keep an initial missing truth value from poisoning the run.
    x=np.asarray(x)
    return np.array([np.mean(x[max(0,j-n+1):j+1]) for j in range(len(x))])
def med(x):
    x=np.asarray(x); x=x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else np.nan

def sd(x):
    x=np.asarray(x); x=x[np.isfinite(x)]
    return float(np.std(x,ddof=1)) if len(x)>1 else np.nan

def moments(x):
    x=np.asarray(x); x=x[np.isfinite(x)]
    return dict(n=len(x),mean=float(np.mean(x)) if len(x) else np.nan,sd=sd(x))

def regress(truth,signal,mask,shift=0):
    # Positive shift means signal follows truth by shift*Tc.
    ids=np.flatnonzero(mask); jj=ids+shift
    ok=(jj>=0)&(jj<len(signal)); ids,jj=ids[ok],jj[ok]
    ok=mask[jj]&np.isfinite(truth[ids])&np.isfinite(signal[jj])
    x,y=truth[ids[ok]],signal[jj[ok]]
    if len(x)<3 or np.std(x)<1e-10 or np.std(y)<1e-10:
        return dict(n=len(x),correlation=np.nan,gain=np.nan,intercept=np.nan,lag_s=shift*TC)
    xc,yc=x-x.mean(),y-y.mean()
    gain=float(xc@yc/(xc@xc))
    return dict(n=len(x),correlation=float(np.corrcoef(x,y)[0,1]),gain=gain,
                intercept=float(y.mean()-gain*x.mean()),lag_s=shift*TC)

def local_harmonics(values):
    # Centered 24-tick OLS, quadratic trend and both known frequency pairs.
    # RMS across local amplitudes avoids cancellation on opposite gradients.
    offsets=np.arange(-12,12)*TC; u=offsets/(24*TC)
    basis=np.column_stack([np.ones(24),u,u*u,np.sin(2*np.pi*FX*offsets),
        np.cos(2*np.pi*FX*offsets),np.sin(2*np.pi*FP*offsets),np.cos(2*np.pi*FP*offsets)])
    inv=np.linalg.pinv(basis)
    out=np.full((len(values),2),np.nan)
    for j in range(12,len(values)-11):
        co=inv@values[j-12:j+12]
        out[j]=np.hypot(co[3::2],co[4::2])
    return out

def analyze(path):
    with np.load(path,allow_pickle=False) as f:
        log,rows,pre,eb,calls=[f[k] for k in ('log','seen','preprocessed','eb','calls')]
        meta=json.loads(str(f['meta']))
    p=meta['policy']
    assert p['enc'] and p['Tc']==TC and meta['arm_zeta']<0
    N=len(log)
    Y,X,Z,H1=[np.full((N,24),np.nan) for _ in range(4)]
    lateral=pre.copy(); lateral[:,2]=rows[:,1]
    D4=np.diff(np.eye(24),n=4,axis=0)
    rho=(2*np.sin(np.pi*14*TC))**-8
    # Spectral square root of the specified M. This avoids subtraction error
    # in I-inverse(...) manufacturing information about a linear scan.
    _,sing,Vt=np.linalg.svd(D4,full_matrices=False)
    W=np.sqrt(rho*sing**2/(1.+rho*sing**2))[:,None]*Vt
    k,c,ko,cs,trend,gram_corr=[np.full(N,np.nan) for _ in range(6)]
    for j in calls:
        n=j+1; j0=max(p['ndel']+1,n-1-p['nw']); m=n-1-j0
        if m<24: continue
        args=(j0,m,0,p['ns'],eb[j],p['afc'],p['ndel'],pre[j,0])
        _,y,z,keep=common._build(pre,*args,pre[j,2],p['nw']*p['ns']*.001,True)
        _,x,zx,keepx=common._build(lateral,*args,lateral[j,2],p['nw']*p['ns']*.001,True)
        assert np.array_equal(z,zx) and np.array_equal(keep,keepx)
        if not np.all(keep[-24:]): continue
        h=HW.time_weighted(pre,j0,m,p['ns'],eb[j],p['afc'],pre[j,0],pre[j,2],True)[-24:]
        y,z,x=y[-24:],z[-24:],x[-24:]
        Y[j],X[j],Z[j],H1[j]=y,x,z,h
        wy,wx,wz,wh=W@y,W@x,W@z,W@h
        R=np.column_stack([wy,wx]); gram=R.T@R
        a0,b0=wy@wy,wy@wh
        g=wh-b0/a0*wy; trend[j]=g@g
        ko[j]=wy@wz/a0
        k[j]=ko[j]  # identifiable axial coefficient if lateral column is null
        gram_corr[j]=gram[0,1]/np.sqrt(gram[0,0]*gram[1,1]) if np.all(np.diag(gram)>0) else np.nan
        if np.linalg.eigvalsh(gram)[0] > 24*np.finfo(float).eps*np.linalg.eigvalsh(gram)[-1]:
            k[j],c[j]=np.linalg.solve(gram,R.T@wz)
            cs[j]=SIGMA*np.sqrt(np.linalg.inv(gram)[1,1])
    t=log[:,0]; lam=np.log(np.where(log[:,26]>0,log[:,26],np.nan)); lams=avg(lam,8)
    rate=np.gradient(lams,TC); velocity=np.gradient(log[:,1],TC)
    v24,F24=avg(velocity,24),avg(log[:,5],24)
    G=np.divide(rate,v24,out=np.full(N,np.nan),where=v24>.002)
    inside=np.zeros(N,dtype=bool)
    for ev in meta['event_info']: inside|=(t>=ev['start'])&(t<=ev['end']+.3)
    score=(t>1.5)&~inside
    valid=score&np.isfinite(G)&np.isfinite(c/F24)
    plateau,changing=score&(abs(rate)<.2),score&(abs(rate)>1.)
    corr=regress(G,c/F24,valid)
    lagged=[regress(G,c/F24,valid,s) for s in range(-25,26)]
    best=max(lagged,key=lambda q:q['correlation'] if np.isfinite(q['correlation']) else -np.inf)
    kap=best['gain']; c_sd=sd(c[plateau])
    implied=abs(v24)*c_sd/(abs(kap)*abs(F24))
    v2=np.divide(SIGMA,log[:,26]*np.sqrt(np.maximum(trend,0)),out=np.full(N,np.nan),where=(trend>0)&(log[:,26]>0))
    xc=np.r_[0.,np.cumsum((log[1:,16]+log[:-1,16])*.5*TC)]
    hx,hc,hf=[local_harmonics(a) for a in (log[:,1],xc,log[:,18])]
    def ripple(mask):
        return {str(f):float(np.sqrt(np.nanmean(hf[mask,q]**2)/2)) if np.any(mask) else np.nan
                for q,f in enumerate((FX,FP))}
    def err(estimate,mask):
        finite=mask&np.isfinite(estimate); good=finite&(estimate>0)
        return dict(median_abs_log=med(abs(np.log(estimate[good]/log[good,26]))),
            n_selected=int(mask.sum()),n_valid_positive=int(good.sum()),n_nonpositive=int((finite&~good).sum()))
    meta.update(file=path.name,gradient=dict(zero_lag=corr,best_lag=best),
        counts=dict(plateau=int(plateau.sum()),changing=int(changing.sum()),gradient=int(valid.sum()),lateral_rank_deficient=int(np.sum(score & ~np.isfinite(c)))),
        dither=dict(probe_mm=1e3*med(hx[score,0]),contact_integrated_mm=1e3*med(hc[score,0]),
            probe_velocity_ratio=2*np.pi*FX*med(hx[score,0])/meta['v_scan'],
            contact_velocity_ratio=2*np.pi*FX*med(hc[score,0])/meta['v_scan'],
            min_probe_velocity=float(np.min(velocity[score])),min_contact_velocity=float(np.min(log[score,16])),
            contact_reverse_ticks=int(np.sum(log[score,16]<0))),
        ripple=dict(plateau=ripple(plateau),changing=ripple(changing)),
        noise=dict(lateral_identifiable=bool(np.any(score & np.isfinite(c))),c_sd=c_sd,c_sd_predicted_median=med(cs[plateau]),rate_noise_median=med(implied[plateau]),
            v2_rate_noise_median=med(v2[plateau]),rate_noise_all_scored_median=med(implied[score]),
            v2_all_scored_median=med(v2[score]),F24_plateau_median=med(F24[plateau]),
            v24_plateau_median=med(v24[plateau]),normalized_plateau_sd=sd((c/F24)[plateau]),
            regressor_correlation_abs_median=med(abs(gram_corr[score]))),
        stiffness={s:dict(two_regressor=err(k,mask),k_only=err(ko,mask))
                   for s,mask in (('all',score),('plateau',plateau),('changing',changing))})
    if meta['index']==12:
        mi,mo=moments((c/F24)[(t>1.5)&inside]),moments((c/F24)[score])
        ps=sd((c/F24)[plateau])
        meta['body_motion']=dict(inside=mi,outside=mo,plateau_sd=ps,mean_difference=mi['mean']-mo['mean'],
                               mean_difference_in_plateau_sd=abs(mi['mean']-mo['mean'])/ps)
    arrays=dict(t=t,G_true=G,rate_true=rate,lambda_true=lam,lambda_smoothed=lams,k_hat=k,k_only=ko,
        c_hat=c,F24=F24,v24=v24,cs_predicted=cs,v2_rate_noise=v2,Y=Y,X=X,z=Z,H1=H1,score=score,
        plateau=plateau,changing=changing,force_harmonic_amplitudes=hf,probe_harmonic_amplitudes=hx,
        contact_harmonic_amplitudes=hc)
    return meta,arrays

def verdicts(rows):
    get={(r['index'],r['amp_mm']):r for r in rows}
    details={key:{} for key in PREDICTIONS}
    for i in CASES:
        r,b=get[i,.1],get[i,0.]; g,no=r['gradient']['best_lag'],r['noise']
        if i in (36,39):
            details['L1'][i]=dict(**g,met=bool(g['correlation']>=.8 and .3<=g['gain']<=1.5))
            n,v=no['rate_noise_median'],no['v2_rate_noise_median']
            details['L2'][i]=dict(rate_noise=n,v2_noise=v,reduction_factor=v/n,met=bool(n<=1 and n<=v/3))
        adds={s:r['ripple'][s][str(FX)]-b['ripple'][s][str(FX)] for s in ('plateau','changing')}
        dl=r['L']-b['L']
        details['L3'][i]=dict(added_ripple=adds,delta_L=dl,
            met=bool(adds['plateau']<=.02 and (not np.isfinite(adds['changing']) or adds['changing']<=.08) and dl<=.03),
            no_changing_ticks=not np.isfinite(adds['changing']))
        errors={s:r['stiffness'][s]['two_regressor']['median_abs_log']-b['stiffness'][s]['k_only']['median_abs_log']
                  for s in ('all','plateau','changing')}
        details['L4'][i]=dict(error_increase=errors,met=bool(errors['all']<=.01),
                            segmentwise_met=all(v<=.01 for v in errors.values() if np.isfinite(v)))
    bm=get[12,.1]['body_motion']
    details['L5'][12]=dict(**bm,met=bool(bm['mean_difference_in_plateau_sd']<=2))
    return {k:dict(prediction=PREDICTIONS[k],met=all(r['met'] for r in d.values()),cases=d) for k,d in details.items()}

def plot(row,a,folder):
    # Figure contract: quantitative grid; gradient tracking (top) and axial fit
    # disturbance (bottom). One simulated episode, full resolution, no CI.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0,str(SKILL/'scripts'))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'pdf.fonttype':42,'svg.fonttype':'none',
        'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
    fig,axes=plt.subplots(2,1,figsize=(7.2,4.5),sharex=True)
    fig.subplots_adjust(left=.12,right=.98,bottom=.11,top=.86,hspace=.32)
    kap=row['gradient']['best_lag']['gain']; vis=a['t']>=1.5; t=a['t'][vis]
    axes[0].plot(t,a['G_true'][vis],color='.2',lw=1,label='True gradient')
    axes[0].plot(t,(a['c_hat']/(kap*a['F24']))[vis],color='#26799b',lw=.7,alpha=.8,label='Lateral estimate / fitted gain')
    axes[0].set_ylabel('Gradient (1/m)')
    axes[0].legend(loc='lower left',bbox_to_anchor=(0,1.01),ncol=2,fontsize=7)
    axes[1].plot(t,a['lambda_true'][vis],color='.2',lw=1,label='True effective stiffness')
    axes[1].plot(t,np.log(np.where(a['k_hat'][vis]>0,a['k_hat'][vis],np.nan)),color='#26799b',lw=.7,label='Joint fit')
    axes[1].set_ylabel('ln stiffness (N/m)'); axes[1].set_xlabel('Time (s)')
    axes[1].legend(loc='lower left',bbox_to_anchor=(0,1.01),ncol=2,fontsize=7)
    axes[1].set_xlim(1.5,10.)
    for ax in axes:
        for ev in row['event_info']: ax.axvspan(ev['start'],ev['end']+.3,color='.8',alpha=.4,zorder=-1)
    fig.suptitle(f'Case {row["index"]}, commanded dither 0.1 mm; fitted gain {kap:.3g}',y=.98,fontsize=10)
    dest=folder/f'case_{row["index"]:02d}_a0.10'
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=str(dest)+'.alignment.json',strict=True)
    for ext in ('png','pdf','svg','tiff'):
        target=Path(str(dest)+'.'+ext)
        if target.exists(): raise FileExistsError(target)
    fig.savefig(str(dest)+'.png',dpi=300)
    fig.savefig(str(dest)+'.pdf',dpi=300)
    fig.savefig(str(dest)+'.svg',dpi=300)
    fig.savefig(str(dest)+'.tiff',dpi=300)
    plt.close(fig)

def subprocess_run(task):
    i,amp,folder=task
    stdout=Path(folder)/f'case_{i:02d}_a{amp:.2f}.stdout'
    with stdout.open('x',encoding='utf-8') as stream:
        result=subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).resolve()),
            '--worker-index',str(i),'--worker-amp',str(amp),'--worker-folder',folder],
            stdin=subprocess.DEVNULL,stdout=stream,stderr=stream,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    if result.returncode: raise RuntimeError(stdout.read_text(encoding='utf-8'))
    print(stdout.read_text(encoding='utf-8').strip(),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path,default=ROOT)
    parser.add_argument('--analyze-only',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--worker-index',type=int)
    parser.add_argument('--worker-amp',type=float)
    parser.add_argument('--worker-folder')
    args=parser.parse_args()
    if args.worker_index is not None:
        run((args.worker_index,args.worker_amp,args.worker_folder)); return
    root=args.output_root.resolve()
    out=root/'experiments/rossi_failure_study/results'; runs=out/'lateral_dither'
    figs=root/'research/figures/lateral_dither'
    runs.mkdir(parents=True,exist_ok=True); figs.mkdir(parents=True,exist_ok=True)
    if not args.analyze_only:
        dump_new(runs/'predictions_before_episodes.json',dict(predictions=PREDICTIONS,unix_time=time.time(),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
        tasks=[(i,a,str(runs)) for i in CASES for a in AMPS]
        with ThreadPoolExecutor(args.workers) as pool:
            futures=[pool.submit(subprocess_run,t) for t in tasks]
            for f in as_completed(futures): f.result()
    summaries,arrays=[],{}
    for i in CASES:
        for amp in AMPS:
            r,a=analyze(runs/f'case_{i:02d}_a{amp:.2f}.npz')
            summaries.append(r); arrays[i,amp]=a
            with (runs/f'case_{i:02d}_a{amp:.2f}_analysis.npz').open('xb') as f: np.savez_compressed(f,**a)
    for r in summaries:
        b=next(q for q in summaries if q['index']==r['index'] and q['amp_mm']==0.)
        r['delta_L']=r['L']-b['L']
        r['added_ripple_15p625']={s:r['ripple'][s][str(FX)]-b['ripple'][s][str(FX)] for s in ('plateau','changing')}
    for i in (36,39,12):
        plot(next(r for r in summaries if r['index']==i and r['amp_mm']==.1),arrays[i,.1],figs)
    output=dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),predictions=verdicts(summaries),runs=summaries,conventions=dict(
        truth='ln(log[:,26]), trailing 8-tick mean, centered finite-difference derivative',
        masks='t>1.5 s; outside start..end+0.3 s; plateau abs(rate)<0.2; changing abs(rate)>1; gradient v24>0.002 m/s',
        windows='Exact 24-interval estimator _build; actual estimator preprocessed rows and eb; X uses raw measured x without H',
        covariance='sigma^2 inverse(R.T M R), surrogate; M evaluated through spectral square root from D4 SVD to preserve its nullspace; rank requires min eigenvalue>24*eps*max eigenvalue; rank-deficient lateral fit c/covariance undefined, axial k_hat=k_only',
        lag='Positive lag means signal follows truth; maximize positive Pearson correlation over -25..25 ticks; both endpoints valid',
        gain='OLS with intercept; best-lag gain used for implied rate noise and figure normalization',
        rate_noise='median abs(v24)*plateau sd(c)/(abs(best-lag gain)*abs(F24)); plateau median primary; all scored median also retained',
        harmonic='24-tick centered OLS: quadratic trend plus joint sin/cos at 15.625 and 20.833333 Hz; median displacement amplitude; sqrt(mean(local force amplitude squared)/2) for force',
        contact='Trapezoidal integration of logged contact velocity, without attenuation correction',
        added_ripple='Signed rms difference from no-dither case',
        L4='All scored ticks primary; plateau/changing differences also retained; nonpositive fits excluded from log and counted',
        limitations='One episode per cell, common seed within case, overlapping windows; truth k_eff includes depth/body/velocity dependence; gain/lag are offline'),
        figures=[str((figs/f'case_{i:02d}_a0.10.png').relative_to(root)) for i in (36,39,12)])
    dump_new(out/'lateral_dither_check.json',output)
    print(json.dumps(T.CC.clean(output['predictions']),indent=2),flush=True)
if __name__=='__main__': main()







