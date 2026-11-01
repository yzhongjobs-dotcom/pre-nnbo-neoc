"""Closed-loop evaluation with frozen open-loop selection, paired scenario bootstrap.
Run after step3b_predictors.py. Temporary transport/checkpoints stay in --cache.
"""
import argparse,json,time
from copy import deepcopy
from step3b_core import *
import sweep
_MODELS=None

def episode(job):
    global _MODELS
    name,gain,kind,dose,scene=job
    if _MODELS is None: _,_MODELS=load_models()
    fam,sid,p,T,code=scene;use_plant(code)
    t0=time.perf_counter()
    c=build(name,p,deepcopy(_MODELS),gain)
    log,c=common.run(c,p,T,sid)
    mt=common.task_metrics(log,20.8)
    calls=np.array(c.t_call[100:])*1e3
    return dict(model=name,gain=gain,kind=kind,dose=dose,family=fam,sid=sid,code=code,
                ok=not bool(mt['diverged']),rmse=float(mt['rmse']) if np.isfinite(mt['rmse']) else None,
                peak=float(mt['peak']) if np.isfinite(mt['peak']) else None,
                loss_max=float(mt['loss_max']) if np.isfinite(mt['loss_max']) else None,
                loss=float(mt['loss']),p50_ms=float(np.percentile(calls,50)),p99_ms=float(np.percentile(calls,99)),elapsed_s=time.perf_counter()-t0)

def aggregate(rows):
    good=[r for r in rows if r['ok']]
    rm=[r['rmse'] for r in good]
    return dict(n=len(rows),diverged=len(rows)-len(good),runaway=sum((not r['ok']) or r['peak']>50 for r in rows),
                long_loss_05=sum((not r['ok']) or r['loss_max']>.5 for r in rows),
                long_loss_1=sum((not r['ok']) or r['loss_max']>1 for r in rows),
                median=float(np.median(rm)) if rm else None,mean=float(np.mean(rm)) if rm else None,
                timing_p50_ms=float(np.median([r['p50_ms'] for r in rows])),timing_p99_ms=float(np.percentile([r['p99_ms'] for r in rows],99)))

def bootstrap(rows,name,seed=219901):
    a={r['sid']:r for r in rows if r['model']==name};b={r['sid']:r for r in rows if r['model']=='Proposed'}
    ids=sorted(set(a)&set(b));valid=[i for i in ids if a[i]['ok'] and b[i]['ok']]
    delta=np.array([a[i]['rmse']-b[i]['rmse'] for i in valid])
    if not len(delta): return dict(n=0)
    rng=np.random.default_rng(seed)
    boots=delta[rng.integers(0,len(delta),(10000,len(delta)))].mean(axis=1)
    # Missing/diverged are explicitly counted; no arbitrary RMSE imputation.
    return dict(n=len(delta),excluded_divergent_pairs=len(ids)-len(valid),mean_delta=float(delta.mean()),ci95=np.percentile(boots,[2.5,97.5]).tolist(),draws=10000)

def summary(rows,names):
    out={}
    for fam in FAMS:
        rr=[r for r in rows if r['family']==fam]
        out[fam]={n:aggregate([r for r in rr if r['model']==n]) for n in names}
        for n in names:
            if n!='Proposed': out[fam][n]['paired']=bootstrap(rr,n)
    return out

def ss_radius(kappa,kg,As,Bs,Cs):
    import ff_loop as ff
    # Same physical force-loop and filter inversion as ff_loop.tick, FINAL gains.
    ff.M1,ff.D1,ff.KPF1=(CE.FINAL[k] for k in ('M1','D1','KPF1'))
    d=len(Bs);nx=9+d+1
    def tick(x):
        y,cp,Ff,I,wd,cmd,mp=x[:7];Fp,yp=x[7:9];xh=x[9:9+d];spp=x[-1]
        e,ep=-Ff,-mp
        ed=math.exp(-ff.D1*DT/ff.M1);wd2=wd*ed+(1-ed)*e/ff.D1
        fin=(Ff-ff.DEC*Fp)/(1-ff.DEC);ybar=yp+ff.FRAC*(y-yp)
        z=kg/kappa*fin-ybar
        xh2=As@xh+Bs*z;sp=Cs@xh2
        cmd2=cmd+.5*DT*(wd+wd2)+ff.KPF1*(e-ep)-(sp-spp)
        cp2=cmd
        for _ in range(ff.NS):
            Ff+=ff.ALPHA*(kappa*y-Ff);y+=ff.BETA*(cp2-y)
        return np.r_[[y,cp2,Ff,I,wd2,cmd2,-e],[x[2],x[0]],xh2,sp]
    eye=np.eye(nx);A=np.column_stack([tick(eye[:,i]) for i in range(nx)])
    keep=[i for i in range(nx) if i!=3]
    ev=np.linalg.eigvals(A[np.ix_(keep,keep)]);ev=ev[np.abs(ev-1)>1e-7]
    return float(np.max(abs(ev)))

def filter_ss(A,Q,R,H):
    from scipy.linalg import solve_discrete_are
    d=len(H)
    # DARE is the steady predicted covariance; avoids slow scalar random-walk
    # convergence when the stationary mode has tiny process variance.
    h=np.zeros((1,d));h[0,0]=1.
    pp=solve_discrete_are(A.T,h.T,Q,np.array([[R]]))
    K=pp[:,0]/(pp[0,0]+R)
    return (np.eye(d)-np.outer(K,np.eye(d)[0]))@A,K,H

def linear_check(names,models):
    kggrid=np.unique(np.r_[np.arange(.5,1.501,.025),1/1.2,1.2,.75,.82,1.18,1.25]);out={}
    for name in names:
        model=models.get(name)
        if name=='P3': continue
        records=[]
        for k in (1e3,3e3,1e4,3e4,1e5):
            for kap in (.5,1.,2.):
                if name=='P2': modes=[('FIR',None)]
                else:
                    if name=='Proposed':
                        mm=MultiIMM(CE.FINAL['q_smooth'],CE.FINAL['q_maneuver'],stay=CE.FINAL['p_stay'],third=False)
                    else: mm=model
                    modes=[(str(i),i) for i in range(mm.m)]
                for mode,index in modes:
                    radii=[]
                    for kg in kggrid:
                        if name=='P2':
                            A=np.zeros((W,W));A[:-1,1:]=np.eye(W-1);B=np.zeros(W);B[-1]=1.;C=model.weights
                        else:
                            d=1 if mm.m==3 and index==0 else (4 if mm.m==3 and index==2 else 3)
                            A,B,C=filter_ss(mm.A[index,:d,:d],mm.Q[index,:d,:d],(CE.FINAL.get('noise',.04)/(k/kg))**2*CE.FINAL['r_scale'],mm.H[index,:d])
                        radii.append(ss_radius(kap,float(kg),A,B,C))
                    i=int(np.argmin(abs(kggrid-1)));lo=hi=i
                    if radii[i]<1-1e-9:
                        while lo>0 and radii[lo-1]<1-1e-9: lo-=1
                        while hi<len(kggrid)-1 and radii[hi+1]<1-1e-9: hi+=1
                        interval=[float(kggrid[lo]),float(kggrid[hi])]
                    else: interval=None
                    within=(kggrid>=1/1.2-1e-9)&(kggrid<=1.2+1e-9)
                    records.append(dict(k=k,kappa=kap,mode=mode,interval=interval,max_radius_accepted=float(np.max(np.array(radii)[within])),grid=kggrid.tolist(),radius=radii))
        out[name]=dict(records=records,all_accepted_stable=all(r['max_radius_accepted']<1-1e-9 for r in records),
                       common_interval=[max(r['interval'][0] for r in records),min(r['interval'][1] for r in records)] if all(r['interval'] for r in records) else None)
    return out

def write_md(out):
    names=out['models'];lines=['# Step 3b: closed-loop replacement','',
        'Candidates were frozen using only validation prediction MSE before ID/OOD or closed-loop evaluation. All other controller parameters remain FINAL. No confidence-dependent feed-forward attenuation was introduced. All errors are true force RMSE in N.', '',
        '| Family | Model | Median / mean | Peak>50 N (incl divergence) | Long loss >0.5s | Paired mean Δ [95% CI] |', '|---|---|---:|---:|---:|---|']
    for fam,g in out['diagnostic'].items():
        for n,r in g.items():
            p=r.get('paired');ci='' if not p else f"{p['mean_delta']:+.3f} [{p['ci95'][0]:+.3f}, {p['ci95'][1]:+.3f}]"
            lines.append(f"| {fam} | {n} | {r['median']:.3f}/{r['mean']:.3f} | {r['runaway']}/{r['n']} | {r['long_loss_05']}/{r['n']} | {ci} |")
    lines+=['','## Dose response','', '| Kind / dose | '+' | '.join(names)+' |', '|---|'+'---|'*len(names)]
    for key,g in out['sweep'].items():
        lines.append('| '+key.replace('|',' / ')+' | '+' | '.join(f"{g[n]['median']:.3f}/{g[n]['mean']:.3f} ({g[n]['runaway']}/{g[n]['n']})" for n in names)+' |')
    lines+=['','## Measurement-gain tolerance','',
        'Linear analysis reproduces ff_loop.py substep dynamics, inverse force filter, encoder interpolation and one-command delay, with FINAL force gains. P1 freezes each IMM mode and recomputes its steady Kalman gain for each κ_g; arbitrary switching is outside this analysis. P2 uses its complete 63-tap FIR realization. The accepted log-gain band b≤ln1.2 corresponds to κ_g∈[1/1.2,1.2]. Intervals are grid observations, not analytic stability guarantees.', '',
        '| Predictor | Common grid interval | Entire accepted band stable? |','|---|---|---|']
    for n,r in out['linear_tolerance'].items(): lines.append(f"| {n} | {r['common_interval']} | {r['all_accepted_stable']} |")
    lines+=['','Gain stress scales inverse measurement stiffness ONLY inside drift prediction; feedback stiffness, excitation normalization and gain acceptance stay FINAL. It multiplies the existing gain error rather than imposing known true stiffness. P3 is assessed empirically at 0.75/0.82/1.18/1.25; this cannot certify a continuum or convert an estimator bias bound into a stability proof.', '',
        '| Gain factor | Model | Runaway | Long loss>0.5s | Median / mean |','|---|---|---:|---:|---:|']
    for gain,g in out['gain_stress'].items():
        for n,r in g.items(): lines.append(f"| {gain} | {n} | {r['runaway']}/{r['n']} | {r['long_loss_05']}/{r['n']} | {r['median']:.3f}/{r['mean']:.3f} |")
    lines+=['','## Computation','', 'Whole-controller hot serial timing (ms), one diagnostic scene/family, includes stiffness estimation, predictor and control, excludes plant and cold compilation:', '', '| Model | Median P50 | P99 of scene P99 |','|---|---:|---:|']
    for n,r in out['serial_timing'].items(): lines.append(f"| {n} | {r['timing_p50_ms']:.3f} | {r['timing_p99_ms']:.3f} |")
    lines+=['','## Statistical boundaries','',
        'Diagnostic seeds110000+,50/family; sweeps exactly reuse sweep.py ramp doses0.6/0.45/0.3/0.2/0.15/0.1s and ridge ratios3/10/20/30/45/60,20/dose, seeds160000+/170000+. No final cohort accessed. Paired bootstrap resamples scenes within family,10000 replicates, seed219901; CI is percentile95%, exploratory and without multiplicity correction. Means/medians include finite runaway trajectories. Divergent pairs would be excluded only from RMSE/CI and explicitly counted, never from safety counts. TrueK/Oracle+exc are privileged references, not feasible competitors. Long loss means maximum contiguous true force<0.5N duration>0.5s. Runtime in parallel workers is retained for auditing but serial hot timing is used for interpretation. No evaluated outcome feeds selection or retuning.']
    (RES/'step3b_closed_loop.md').write_text('\n'.join(lines)+'\n',encoding='utf8')

def cached_jobs(jobs,path,workers):
    import pickle,hashlib
    model_hash=hashlib.sha256(json.dumps(json.loads((RES/'step3b_predictors.json').read_text(encoding='utf8'))['selected'],sort_keys=True).encode()).hexdigest()
    if path.exists():
        saved=pickle.loads(path.read_bytes())
        assert saved.get('model_hash')==model_hash, 'stale model cache: use a fresh --cache directory'
        assert saved['keys']==[(j[0],j[1],j[2],j[3],j[4][1]) for j in jobs]
        return saved['rows']
    rows=file_pool('step3b_closed_loop','episode',jobs,path.parent/'closed_workers',workers)
    path.write_bytes(pickle.dumps(dict(model_hash=model_hash,keys=[(j[0],j[1],j[2],j[3],j[4][1]) for j in jobs],rows=rows)))
    return rows

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8);parser.add_argument('--cache',type=Path,required=True);args=parser.parse_args();args.cache.mkdir(exist_ok=True,parents=True)
    cfg,models=load_models();chosen=cfg['closed_loop_candidates'];names=['Proposed',*chosen,'TrueK','Oracle+exc'];started=time.time()
    scenes=make_scenes(110000,50)
    jobs=[(n,1.,'diagnostic',None,s) for n in names for s in scenes]
    print('diagnostic jobs',len(jobs),names,flush=True)
    rows=cached_jobs(jobs,args.cache/'closed_diag.pkl',args.workers)
    print('diagnostic done',round(time.time()-started),flush=True)
    sc=[]
    for kind,dose,sid,p,T in sweep.scenes():
        sc.append((kind,dose,('T2f-stiff' if kind=='ramp' else 'S2',sid,p,T,0)))
    jobs=[(n,1.,kind,dose,s) for n in names for kind,dose,s in sc]
    sw=cached_jobs(jobs,args.cache/'closed_sweep.pkl',args.workers);print('sweep done',round(time.time()-started),flush=True)
    stressnames=['Proposed',*chosen]
    stress=[]
    for gain in (.75,.82,1.18,1.25):
        jobs=[(n,gain,'gain_stress',None,s) for n in stressnames for s in scenes]
        stress+=cached_jobs(jobs,args.cache/f'closed_gain{gain}.pkl',args.workers)
        print('gain done',gain,round(time.time()-started),flush=True)
    print('linear tolerance',flush=True);linear=linear_check(['Proposed','P1','P2'],models)
    # Timing serial in main process: independent of worker congestion.
    serial=[episode((n,1.,'timing',None,s)) for n in names for s in make_scenes(110000,1)]
    out=dict(models=names,diagnostic=summary(rows,names),
             sweep={f'{kind}|{dose}':{n:dict(aggregate([r for r in sw if r['kind']==kind and r['dose']==dose and r['model']==n]),paired=bootstrap([r for r in sw if r['kind']==kind and r['dose']==dose],n)) for n in names} for kind,doses in [('ramp',sweep.RAMPS),('ratio',sweep.RATIOS)] for dose in doses},
             gain_stress={str(g):{n:aggregate([r for r in stress if r['gain']==g and r['model']==n]) for n in stressnames} for g in (.75,.82,1.18,1.25)},
             linear_tolerance=linear,serial_timing={n:aggregate([r for r in serial if r['model']==n]) for n in names},
             scene_summaries=rows+sw+stress,
             elapsed_s=time.time()-started)
    (RES/'step3b_closed_loop.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf8');write_md(out)
    print('completed closed-loop',round(out['elapsed_s']),flush=True)

if __name__=='__main__': main()
