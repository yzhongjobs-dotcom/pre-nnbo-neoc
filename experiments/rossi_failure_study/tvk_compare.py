"""TVK validation. Reuses baseline host, causal replay, scoring and bootstrap.

python -X utf8 experiments/rossi_failure_study/tvk_compare.py --stage pilot
Stages: unit, pilot, dev, test, replay, timing, report (restartable JSON).
All development choices are frozen before looking at diagnostic test results.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import math
import time
from pathlib import Path
from types import SimpleNamespace
from itertools import product
import numpy as np
from scipy.linalg import solveh_banded
import baseline_compare as B
import baseline_report as R
import tvk_mhe as T

REPORT=Path('D:/pip_env/temp/claude/8fbcb2674095e6be/841f52cd-14a8-4db4-b56f-36be6ffcc337/scratchpad/tvk_report.md')
ORIGINAL_BUILD=B.build


def build(method,cfg,sc):
    return T.build(sc,cfg) if method.startswith('TVK') else ORIGINAL_BUILD(method,cfg,sc)


def configs():
    return [dict(q=q,order=p,lam_k=lam,delta_k=.002 if q==2 else .02,gate=g)
            for q,p,lam,g in product((2,1),(4,3),(.01,.1,1.),(.551,1.))]


def execute(job):
    B.build=build
    return B.execute(job)


def solver(q=2,p=4,lam=.1,delta=None):
    e=T.TVKRidge(q=q,order=p,nw=24,lam_k=lam,delta_k=delta or (.002 if q==2 else .02),sigma_scale=1.,chi_scale=1.)
    e.reset(SimpleNamespace(Tc=.008,params={}))
    return e


def constant(e,Y,z,noise):
    m=len(Y);Y=Y-Y.mean();D,b=e._difference_bands(m);w=np.ones(m-e.order)
    for it in range(e.iters):
        ab=(b@w).reshape(e.order+1,m);ab[0]+=1.
        sy,sz=solveh_banded(ab,np.column_stack((Y,z)),lower=True,check_finite=False).T
        my=Y-sy;info=Y@my;k=my@z/info
        if it<e.iters-1:w=np.minimum(1.,e.c/np.maximum(abs(D@(sz-k*sy))/(noise/math.sqrt(e.rho)),1e-12))
    res=z-k*Y;num=res@((z-sz)-k*my)
    return dict(k=float(k),sigma=float(noise/math.sqrt(info)),chi=float(num/noise**2/(m-e.order-1)),d=sz-k*sy)


def unit():
    rng=np.random.default_rng(7863);dt=.008;m=24;t=np.arange(m)*dt;noise=.04
    Y=.00012*np.sin(2*np.pi*20.8*t)+.0002*t
    z=3000*np.exp(2*t)*Y+.5*np.sin(2*np.pi*.3*t)+rng.normal(0,noise,m)
    audits=[]
    for q,p in product((2,1),(4,3)):
        e=solver(q,p);a=e.solve_window(Y,z,noise,True);b=e.solve_window(Y+3.7,z,noise)
        inv=float(np.max(abs(a['k']-b['k']))/np.max(abs(a['k'])))
        assert inv<1e-7,inv
        assert a['dense_var_relative_error']<1e-6,a
        assert a['normal_residual']<1e-8,a
        audits.append(dict(q=q,p=p,reference_relative_error=inv,**{k:a[k] for k in ('dense_theta_error','dense_var_relative_error','dense_trace_error','normal_residual')}))
    limits=[];ref=constant(solver(),Y,z,noise)
    policy=SimpleNamespace(Tc=dt,params=dict(noise=noise,enc_y=0.),ndel=0,ns=8,eb=.5,afc=.2)
    original=T.common._build
    T.common._build=lambda *args:(None,Y.copy(),z.copy(),np.ones(m,bool))
    try:
        from d_robust import RobustRidge
        baseline=RobustRidge(nw=m,order=4,fc=14.);baseline.reset(policy)
        log=np.zeros((m+10,26));log[:,0]=np.arange(len(log))*dt
        real_baseline=baseline.estimate(SimpleNamespace(log=log),len(log),policy)
    finally:T.common._build=original
    const_audit={key:float(abs(ref[key]-value)/max(abs(value),1e-15)) for key,value in zip(('k','sigma','chi'),real_baseline)}
    assert max(const_audit.values())<1e-6,const_audit
    for q,lam in product((1,2),(1e2,1e4,1e6,1e7)):
        e=solver(q,lam=lam);a=e.solve_window(Y,z,noise)
        limits.append(dict(q=q,lam_k=lam,k=float(a['k'][-1]),sigma=a['sigma'],chi=a['chi'],legacy_chi=a['legacy_chi'],
                           relative_to_constant={k:float(abs((a['k'][-1] if k=='k' else a[k])-ref[k])/abs(ref[k])) for k in ('k','sigma','chi')}))
    variants={'q2':solver(2),'q1':solver(1),'constant':solver(1)}
    ramp={}
    for label,e in variants.items():
        rng=np.random.default_rng(16305)  # common phases AND noise for all variants
        hats=[];sigs=[]
        for j in range(160):
            phase=rng.uniform(0,2*np.pi);yy=.00012*np.sin(2*np.pi*20.8*t+phase)
            kk=3000*np.exp(2*(t-t[-1]));dd=.8*np.sin(2*np.pi*.3*t+.4)
            # Generate with the specified reference-invariant observation model.
            D,b=e._difference_bands(m);ab=(b@np.ones(m-e.order)).reshape(e.order+1,m);ab[0]+=1
            sy=solveh_banded(ab,yy,lower=True);yt=yy-sy
            zz=kk.mean()*yy+(kk-kk.mean())*yt+dd+rng.normal(0,noise,m)
            a=constant(e,yy,zz,noise) if label=='constant' else e.solve_window(yy,zz,noise)
            hats.append(a['k'] if label=='constant' else a['k'][-1]);sigs.append(a['sigma'])
        ramp[label]=dict(n=160,endpoint_bias=float(np.mean(hats)-3000),empirical_sd=float(np.std(hats,ddof=1)),mean_sigma=float(np.mean(sigs)))
    steps=[];tt=np.arange(250)*dt;yy=.00012*np.sin(2*np.pi*20.8*tt)
    # Raw k*Y is physical here: report mismatch to the proposed model honestly.
    eps=rng.normal(0,noise,len(tt));pulse=1.2*np.exp(-((tt-.85)/.035)**2)
    for ratio,label in product((2,4),variants):
        e=variants[label];kk=np.where(tt<.8,3000.,3000.*ratio);zz=kk*yy+.8*np.sin(2*np.pi*.3*tt)+pulse+eps
        hats=[];sigs=[];chi=[];ends=[]
        for j in range(m,len(tt)):
            a=constant(e,yy[j-m:j],zz[j-m:j],noise) if label=='constant' else e.solve_window(yy[j-m:j],zz[j-m:j],noise)
            hats.append(a['k'] if label=='constant' else a['k'][-1]);sigs.append(a['sigma']);chi.append(a['chi']);ends.append(tt[j-1])
        hats=np.array(hats);sigs=np.array(sigs);ends=np.array(ends);truth=np.where(ends<.8,3000.,3000.*ratio)
        under=(hats<.9*truth)&(ends>=.8);post=(ends>=.8)&(ends<=1.3)
        recover=None
        for j in np.flatnonzero(ends>=.8):
            if j+5<=len(hats) and np.all(abs(hats[j:j+5]/truth[j:j+5]-1)<.2):recover=float(ends[j]-.8);break
        worst=np.flatnonzero(post)[int(np.argmin(hats[post]/truth[post]))]
        steps.append(dict(ratio=ratio,method=label,recovery_s=recover,max_underestimate=float(max(0.,1-np.min(hats[post]/truth[post]))),
                          sigma_pre=float(np.mean(sigs[ends<.75])),sigma_under=float(np.mean(sigs[under])) if under.any() else None,
                          sigma_at_max_under=float(sigs[worst]),max_under_time_s=float(ends[worst]-.8),
                          sigma_first_post=float(sigs[np.flatnonzero(ends>=.8)[0]]),
                          under_2sigma_coverage=float(np.mean(abs(hats[under]-truth[under])<=2*sigs[under])) if under.any() else None))
    modulation=[]
    # Leakage is the L2 projection of d_hat-d onto the true modulated response;
    # signed fractions and ratios >1 are retained, never clipped.
    for f,label in product((3.,6.,9.),variants):
        e=variants[label];leaks=[];errs=[]
        for phase in np.linspace(0,2*np.pi,32,endpoint=False):
            kk=3000*(1+.3*np.sin(2*np.pi*f*t+phase));yy0=.00012*np.sin(2*np.pi*20.8*t+.3)
            dd=.8*np.sin(2*np.pi*.3*t+.4)
            z0=kk*yy0+dd
            a=constant(e,yy0,z0,noise) if label=='constant' else e.solve_window(yy0,z0,noise)
            # d has an arbitrary constant reference: remove DC from both
            # projection vectors, so the leakage metric is reference invariant.
            target=(kk-kk.mean())*yy0
            target=target-target.mean()
            dh=a['d']-dd;dh=dh-dh.mean()
            leaks.append(float(target@dh/(target@target)));errs.append(float((a['k'] if label=='constant' else a['k'][-1])/kk[-1]-1))
        modulation.append(dict(f_Hz=f,method=label,leakage_mean=float(np.mean(leaks)),leakage_p10_p90=np.percentile(leaks,[10,90]).tolist(),endpoint_relative_rmse=float(np.sqrt(np.mean(np.array(errs)**2)))))
    ramp_shape=[]
    kk=3000*np.exp(2*(t-t[-1]));yy=.00012*np.sin(2*np.pi*20.8*t+.3)
    e=solver();_,bb=e._difference_bands(m);ab=(bb@np.ones(m-e.order)).reshape(e.order+1,m);ab[0]+=1
    sy=solveh_banded(ab,yy,lower=True);zz=kk.mean()*yy+(kk-kk.mean())*(yy-sy)+.8*np.sin(2*np.pi*.3*t+.4)
    for label,e in variants.items():
        a=constant(e,yy,zz,noise) if label=='constant' else e.solve_window(yy,zz,noise)
        traj=np.full(m,a['k']) if label=='constant' else a['k']
        ramp_shape.append(dict(method=label,trajectory=traj.tolist(),truth=kk.tolist(),endpoint_relative_bias=float(traj[-1]/kk[-1]-1),
                               near_flat_fraction=float(np.mean(abs(np.diff(traj))<.1*np.diff(kk))),
                               relative_slope_median=float(np.median(np.diff(traj)/np.diff(kk)))))
    B.dump('tvk_unit.json',dict(audits=audits,constant_actual_RobustRidge_relative_errors=const_audit,constant_reference={k:ref[k] for k in ('k','sigma','chi')},limits=limits,ramp=ramp,ramp_shape=ramp_shape,steps=steps,modulation=modulation,
                              caveat='q2 infinite penalty retains slope; effective-df chi differs from baseline. Synthetic ramp is model matched; steps/modulation are physical raw kY.'))
    print('UNIT',audits,'RAMP',ramp,flush=True)


def pilot():
    jobs=[('TVK',configs()[j],('breath','relaxed','neural_white',21.,0),B.D.BASE,0.,False,True) for j in (0,12)]
    jobs.append(('Huber-current',{},('breath','relaxed','neural_white',21.,0),B.D.BASE,0.,False,True))
    tic=time.perf_counter();rows=B.run(jobs,worker=execute)
    B.dump('tvk_pilot.json',dict(rows=rows,seconds=time.perf_counter()-tic,workers=15,blas_threads=1))
    print([(r['method'],r['seconds'],r['J']) for r in rows],flush=True)


def dev():
    cfgs=configs();B.dump('tvk_plan.json',dict(configs=cfgs,criterion='mean log1p(J), original 8 DEV scenes only',test_base=B.TEST_BASE,workers=15,blas_threads=1))
    cells=[(c,t,'neural_white',21.,i) for c,t in B.D.DEV for i in range(2)]
    tic=time.perf_counter();rows=B.run([('TVK',cfg,cell,B.D.BASE,0.,False,True) for cfg in cfgs for cell in cells],worker=execute)
    scores=[B.L(rows[j*8:(j+1)*8]) for j in range(len(cfgs))];j=int(np.argmin(scores))
    # q1 gets its own best development setting; no diagnostic test selection.
    j1=min((j for j,c in enumerate(cfgs) if c['q']==1),key=lambda j:scores[j])
    j2=min((j for j,c in enumerate(cfgs) if c['q']==2),key=lambda j:scores[j])
    B.dump('tvk_dev.json',dict(configs=cfgs,scores=scores,index=j,q1_index=j1,q2_index=j2,chosen={'TVK-selected':cfgs[j],'TVK-q1':cfgs[j1],'TVK-q2':cfgs[j2]},rows=rows,seconds=time.perf_counter()-tic))
    print('SELECTED',j,scores[j],cfgs[j],'; Q1',j1,scores[j1],cfgs[j1],flush=True)


def test():
    dev=B.load('tvk_dev.json');chosen=dev['chosen'];saved=dict(chosen=chosen,rows=[],seconds=0.,base=B.TEST_BASE,note='Previously seen grid; diagnosis only, no retuning')
    if (B.OUT/'tvk_test.json').exists():saved=B.load('tvk_test.json')
    for label,cfg in chosen.items():
        if any(r['method']==label for r in saved['rows']):continue
        same=next((m for m,c in chosen.items() if m!=label and c==cfg and any(r['method']==m for r in saved['rows'])),None)
        if same:
            saved['rows'] += [dict(r,method=label) for r in saved['rows'] if r['method']==same]
            B.dump('tvk_test.json',saved);continue
        tic=time.perf_counter();rr=B.run([(label,cfg,c,B.TEST_BASE,0.,False,False) for c in B.cells()],worker=execute)
        saved['rows']+=rr;saved['seconds']+=time.perf_counter()-tic;B.dump('tvk_test.json',saved)
        print(label,'TEST',B.L(rr),flush=True)


def refresh_p3():
    """Recompute affected development rows after enforcing the specified S=D4.

    No new configurations and no test-based selection. A normal fresh `dev`
    already uses the correct solver; this stage only avoids rerunning p=4.
    """
    d=B.load('tvk_dev.json');cfgs=d['configs'];indices=[j for j,c in enumerate(cfgs) if c['order']==3]
    cells=[(c,t,'neural_white',21.,i) for c,t in B.D.DEV for i in range(2)]
    tic=time.perf_counter();rows=B.run([('TVK',cfgs[j],cell,B.D.BASE,0.,False,True) for j in indices for cell in cells],worker=execute)
    for offset,j in enumerate(indices):
        d['rows'][j*8:(j+1)*8]=rows[offset*8:(offset+1)*8]
        d['scores'][j]=B.L(rows[offset*8:(offset+1)*8])
    d['seconds']+=time.perf_counter()-tic
    j=int(np.argmin(d['scores']));j1=min((j for j,c in enumerate(cfgs) if c['q']==1),key=lambda j:d['scores'][j]);j2=min((j for j,c in enumerate(cfgs) if c['q']==2),key=lambda j:d['scores'][j])
    old=d['chosen'];d.update(index=j,q1_index=j1,q2_index=j2,chosen={'TVK-selected':cfgs[j],'TVK-q1':cfgs[j1],'TVK-q2':cfgs[j2]})
    B.dump('tvk_dev.json',d)
    assert old==d['chosen'],('development choice changed: diagnostic test/replay must be invalidated',old,d['chosen'])
    m=B.load('tvk_module.json');m['rows']=[r for r in m['rows'] if not (r['method'].startswith('TVK') and m['chosen'][r['method']]['order']==3)]
    m['completed']=[label for label in m['completed'] if not (label.startswith('TVK') and m['chosen'][label]['order']==3)]
    B.dump('tvk_module.json',m)
    print('P3 refreshed; choices unchanged',d['chosen'],flush=True)


def replay_job(job, log_folder='tvk_logs'):
    source,label,cfg=job;B.build=build
    raw=[];holder={}
    def instrument(method,cfg,sc):
        pol=build(method,cfg,sc);holder['pol']=pol;orig=pol.estimator.estimate
        def estimate(o,n,p):
            k,s,c=orig(o,n,p);raw.append((n-1,k,s,c));return k,s,c
        pol.estimator.estimate=estimate;return pol
    B.build=instrument
    r=B.replay_job((source,label,cfg,0.))
    with np.load(B.OUT/source['log_file']) as zz:
        log=np.full((len(zz['log']),26),np.nan);log[:,zz['cols']]=zz['log']
    sc=B.scene(B.key(source),B.TEST_BASE,0.);t=log[:,0];truth=log[:,21]
    ev=np.zeros(len(t),bool)
    for e in sc['event_info']:ev|=(t>=e['start'])&(t<=e['end']+.3)
    grad=np.gradient(np.log(np.maximum(truth,1.)),t)
    change=np.convolve(abs(grad)>2.,np.ones(30),'same')>0
    n=np.array([a[0] for a in raw]);a=np.array([a[1:] for a in raw]);valid=np.isfinite(a[:,0])&(a[:,0]>0)&np.isfinite(a[:,1])
    err=np.full(len(n),np.nan);err[valid]=np.log(a[valid,0]/np.maximum(truth[n[valid]],1.))
    # B's accepted trace is already preserved in the per-group accepted errors.
    # Recover per-window host decisions by capture of the trace after replay.
    # instrument stores the policy so the exact same gate is consumed once.
    accept=holder['pol'].accept_trace
    accepted={round(x[0],7):bool(x[4]) for x in accept}
    ok=np.array([accepted.get(round(float(t[j]),7),False) for j in n])&valid
    eligible=(t[n]>=1.)&(truth[n]>200.)
    group=np.where(ev[n],'event',np.where(change[n],'changing','quiet'))
    extra={}
    for kind in ('quiet','event','changing','all'):
        mask=eligible&(group==kind if kind!='all' else True)
        vv=mask&valid;under=vv&(err<math.log(.9));rising=vv&(grad[n]>2.)
        extra[kind]=dict(n=int(mask.sum()),valid=int(vv.sum()),signed_sum=float(np.nansum(err[vv])),
                         under10=int(under.sum()),under_pass=int((under&ok).sum()),rising_valid=int(rising.sum()),
                         rising_under10=int((rising&(err<math.log(.9))).sum()),
                         covered_all=int(np.sum(abs(err[vv])<=2*a[vv,1]/a[vv,0])),
                         lag_samples_s=(-err[rising]/grad[n[rising]]).tolist())
    r['safety']=extra
    folder=B.OUT/log_folder;folder.mkdir(exist_ok=True)
    name=f"{label}_{r['source']}_{r['cond']}_{r['dynamics']}_{int(r['phase'])}_{r['i']}.npz"
    data=dict(index=n,estimate=a,truth=truth[n],accepted=ok,group=group)
    for kind,g in r['groups'].items():
        for field in ('error_samples','accepted_errors','times'):data[kind+'_'+field]=np.array(g.pop(field))
    for kind,g in extra.items():data[kind+'_lag_samples_s']=np.array(g.pop('lag_samples_s'))
    np.savez_compressed(folder/name,**data);r['samples_file']=log_folder+'/'+name
    return r


def replay():
    dev=B.load('tvk_dev.json');selected=dev['chosen']['TVK-selected']
    chosen=dict(dev['chosen'])
    # Matched p/q variants hold lambda/delta/gate of the corresponding selected q.
    for q,p in product((2,1),(4,3)):
        root=next((c for c in dev['chosen'].values() if c['q']==q),dict(selected,q=q,delta_k=.002 if q==2 else .02))
        chosen[f'TVK-q{q}-p{p}']=dict(root,q=q,order=p)
    base=B.load('baseline_dev_r2.json')['chosen']
    for m in ('Huber-D4','D4','Kalman','MHE'):chosen[m]=base[m]
    chosen['Huber-current']={}
    sources=[r for r in B.load('baseline_test.json')['rows'] if r['method'] in ('Huber-D4','D4') and 'log_file' in r]
    saved=B.load('tvk_module.json') if (B.OUT/'tvk_module.json').exists() else dict(chosen=chosen,rows=[],completed=[],seconds=0.)
    # Deduplicate identical TVK configurations while retaining labels in final results.
    import json
    cache={}
    for label,cfg in chosen.items():
        if label in saved['completed']:cache[json.dumps(cfg,sort_keys=True)]=label;continue
        signature=json.dumps(cfg,sort_keys=True)
        if label.startswith('TVK') and signature in cache:
            rr=[dict(r,method=label) for r in saved['rows'] if r['method']==cache[signature]]
        else:
            tic=time.perf_counter();rr=B.run([(s,label,cfg) for s in sources],worker=replay_job)
            saved['seconds']+=time.perf_counter()-tic
        saved['rows']+=rr;saved['completed'].append(label)
        if label.startswith('TVK'):cache[signature]=label
        B.dump('tvk_module.json',saved);print('REPLAY',label,len(rr),flush=True)


def timing():
    source=next(r for r in B.load('baseline_test.json')['rows'] if r['method']=='Huber-D4' and r['cond']=='breath')
    with np.load(B.OUT/source['log_file']) as zz:
        log=np.full((len(zz['log']),26),np.nan);log[:,zz['cols']]=zz['log'];observed=zz['observed'].copy()
    sc=B.scene(B.key(source),B.TEST_BASE,0.);rows=[]
    for method,cfg in dict(B.load('tvk_dev.json')['chosen'],**{'Huber-current':{}}).items():
        pol=build(method,cfg,sc);pol.reset(dict(B.P3.PLATFORM,tau_a=sc['tau_a_hat']));est=[];orig=pol.estimator.estimate
        def timed(o,n,p):
            tic=time.perf_counter();a=orig(o,n,p)
            if log[n-1,0]>=1.:est.append(time.perf_counter()-tic)
            return a
        pol.estimator.estimate=timed
        for n in range(1,len(log)+1):pol(B.source_observation(log,observed,n))
        rows.append(dict(method=method,estimator_p50_p99_us=(np.percentile(est,[50,99])*1e6).tolist(),
                         policy_p50_p99_us=(np.percentile(np.array(pol.t_call)[125:],[50,99])*1e6).tolist(),estimator_calls=len(est)))
    try:
        from threadpoolctl import threadpool_info
        pools=threadpool_info()
        assert all(x['num_threads']==1 for x in pools if x['user_api']=='blas'),pools
    except ImportError:pools=[]
    B.dump('tvk_timing.json',dict(rows=rows,processes=1,blas_threads=1,threadpools=pools,note='local prototype, not a real-time certificate'))
    print(rows,flush=True)


def aggregate_module(rows):
    """One pass through primary replay arrays; preserve all unfavorable tails."""
    pools={}
    for r in rows:
        with np.load(B.OUT/r['samples_file']) as data:
            for kind,g in r['groups'].items():
                for source,cond in ((r['source'],'all'),('all',r['cond']),('all','all'),(r['source'],r['cond'])):
                    key=(r['method'],source,cond,kind)
                    out=pools.setdefault(key,dict(n=0,valid=0,accepted=0,covered=0,under10=0,under_pass=0,rising_valid=0,rising_under10=0,signed_sum=0.,covered_all=0,errors=[],accepted_errors=[],lag=[]))
                    for f in ('n','valid','accepted','covered'):out[f]+=g[f]
                    for f in ('under10','under_pass','rising_valid','rising_under10','signed_sum','covered_all'):out[f]+=r['safety'][kind][f]
                    out['errors'].append(data[kind+'_error_samples'])
                    out['accepted_errors'].append(data[kind+'_accepted_errors'])
                    out['lag'].append(data[kind+'_lag_samples_s'])
    summary=[]
    for (method,source,cond,kind),a in pools.items():
        errors=np.concatenate(a.pop('errors'));accepted_errors=np.concatenate(a.pop('accepted_errors'));lag=np.concatenate(a.pop('lag'))
        ratio=lambda x,y:float(x/y) if y else None
        summary.append(dict(method=method,source=source,cond=cond,kind=kind,**a,
                            available=ratio(a['accepted'],a['n']),abslog_med=R.median(errors),abslog_p90=float(np.percentile(errors,90)) if len(errors) else None,
                            signed_log_mean=ratio(a['signed_sum'],a['valid']),accepted_abslog_med=R.median(accepted_errors),
                            coverage_accepted=ratio(a['covered'],a['accepted']),coverage_all=ratio(a['covered_all'],a['valid']),
                            under10_rate=ratio(a['under10'],a['valid']),rising_under10_rate=ratio(a['rising_under10'],a['rising_valid']),
                            under_pass_rate=ratio(a['under_pass'],a['under10']),lag_p50_p90_s=np.percentile(lag,[50,90]).tolist() if len(lag) else None))
    return summary


def report():
    u=B.load('tvk_unit.json');dev=B.load('tvk_dev.json');test=B.load('tvk_test.json');module=B.load('tvk_module.json');timing=B.load('tvk_timing.json')
    baseline_dev=B.load('baseline_dev_r2.json')
    primary=B.load('baseline_test_r2.json');base={B.key(r):r for r in primary['rows'] if r['method']=='Huber-current'}
    methods=list(dev['chosen'])+['Huber-current','Huber-D4','D4','Kalman','MHE']
    allrows=test['rows']+[r for r in primary['rows'] if r['method'] in methods]
    stats={m:R.stat([r for r in allrows if r['method']==m]) for m in methods}
    differences={};condition_stats={}
    conds=list(dict.fromkeys(c[0] for c in B.cells()))
    for m in methods:
        rr=[r for r in allrows if r['method']==m]
        differences[m]={'all':R.paired(rr,base)}
        condition_stats[m]={}
        for c in conds:
            sub=[r for r in rr if r['cond']==c]
            differences[m][c]=R.paired(sub,base);condition_stats[m][c]=R.stat(sub)
    ms=aggregate_module(module['rows'])
    pairwise={}
    for m in dev['chosen']:
        rr=[r for r in allrows if r['method']==m]
        pairwise[m]={other:R.paired(rr,{B.key(r):r for r in allrows if r['method']==other}) for other in ('Huber-current','Huber-D4','D4','Kalman','MHE')}
    import hashlib
    filenames=['tvk_mhe.py','tvk_compare.py','baseline_compare.py','baseline_round2.py','d_robust.py']
    provenance={name:hashlib.sha256((B.HERE/name).read_bytes()).hexdigest() for name in filenames}
    runtime=dict(pilot=B.load('tvk_pilot.json')['seconds'],development=dev['seconds'],
                 diagnostic_test=test['seconds'],module_replay=module['seconds'])
    runtime['accounted_total']=sum(runtime.values())
    identities=[r['source_command_max_error'] for r in module['rows'] if r['source_command_max_error'] is not None]
    summary=dict(stats=stats,differences=differences,pairwise=pairwise,condition_stats=condition_stats,module=ms,code_sha256=provenance,
                 runtime_seconds=runtime,source_reproduction=dict(n=len(identities),max_command_error=max(identities)),chosen=dev['chosen'],diagnostic_only=True,workers=15,blas_threads=1)
    B.dump('tvk_summary.json',summary)
    f=R.fmt;tab=R.table
    cn={'breath':'呼吸','shift_fast':'快推拉','flinch':'惊跳','rib':'跨肋骨','cough':'咳嗽','combo':'综合'}
    # Scene identifiers come from the baseline, never assumed in aggregation.
    lines=['# TVK：时变刚度的鲁棒窗口 MHE 验证报告',
           '本次只新增 `tvk_mhe.py`、`tvk_compare.py` 和 `results/tvk_*.json`；逐窗口数组在已忽略的 `results/tvk_logs/`。现有算法、宿主、基线结果未修改。`.gitignore` 的 TVK 忽略项和本次 prompt 记录在开始时已存在，予以保留。未提交、未推送、未修改 HANDOFF 或 decision_log。',
           '## A. 实现与单元检查',
           '**先说明验收条件中的数学冲突：单元检查 1 按原表述不成立，不能宣称全部验收通过。** D²k=0 的零空间是 a+bi；λₖ→∞ 只排除曲率，不能排除斜率。若要常数极限，必须另加趋强的斜率先验，改变本次指定目标函数。q=1 才具有常数极限。另外当前 Huber-D⁴ 用 N−p−1；要求 TVK 用 N−tr(H) 时，χ 不能在常数极限与它完全相同。实现保持指定目标，不增加隐含斜率惩罚，也不对 infinity 特判转调旧估计器。',
           '变量 x=k/K，K 为每窗冷启动二次常数估计的截断尺度（200–100000 N/m），只用观测、全程冻结。先中心化 Y；常数参考由 d 的未惩罚常数分量吸收。令 v=K(Y−SY)，u=KSY，则 A=diag(v)+u1ᵀ/N。目标为 ‖z−Ax−d‖²+2ρs²Σh₁.₃₄₅(Dᵖd/s)+2λₖnoise²Σh₁(Dᑫx/δₖ)，s=noise/√ρ。δₖ 是相对于 K 的每拍增量（q=1）或每拍斜率增量（q=2），分别固定 0.02、0.002；Tc=8 ms。Huber 标准定义在小量处为 u²/2，因此曲率权重为 λₖnoise²/δₖ²。p=3、4 均按 14 Hz 计算各自 ρ。',
           '联合变量按 (xᵢ,dᵢ) 交错，带宽 2max(p,q)。k̄ 使观测矩阵成为带状项加秩一项，正规矩阵是带状项加秩二项；Woodbury 小系统求解。六次冷启动 IRLS，与 RobustRidge 同步使用第六次解的冻结权重，不热启动。σ 来自完整冻结 Hessian 逆的末端 k 对角元，保留所有斜率和 d 的不确定度。χ 分子延续原实现的剖面二次目标（测量残差加冻结先验项），分母用 noise²[N−tr(H)]；另暴露 legacy_chi 便于核对。它们是冻结权重代理统计，不是严格的鲁棒频率学覆盖保证。',
           '窗口、断联截取、H 预处理、相位估计矩阵、噪声尺度和 estimate(o,n,policy) 接口沿用原代码。输出末端 k；负值或越界值沿用原 _clip_k 与宿主拒绝规则；控制器门控逻辑、保守量、前馈和跟踪器均未改。σ 的 .796 和 χ 的 gate 是与 baseline 一样的宿主尺度。',
           '### 1–2. 极限、带状核验与参考不变性',
           tab(['q','p','参考不变相对差','稠密方差相对差','正规方程相对残差'],[[a['q'],a['p'],f(a['reference_relative_error']),f(a['dense_var_relative_error']),f(a['normal_residual'])] for a in u['audits']]),
           '上表参考偏移为 Y+3.7（原探测幅值约 0.00012），逐点 k 轨迹比较；稠密核验还检查帽子矩阵迹。精确数值保存在 tvk_unit.json。',
           '合成常数基准同时通过注入相同回归数组，直接调用实际 RobustRidge.estimate 核验；k/σ/χ 的相对差为 '+str(u['constant_actual_RobustRidge_relative_errors'])+'，没有以复制出的错误基准检验极限。',
           tab(['q','λₖ','末端 k','σ','χ（有效 df）','χ（旧 df）','k 相对常数差','σ 相对差','χ 相对差'],[[a['q'],a['lam_k'],f(a['k']),f(a['sigma']),f(a['chi']),f(a['legacy_chi']),f(a['relative_to_constant']['k']),f(a['relative_to_constant']['sigma']),f(a['relative_to_constant']['chi'])] for a in u['limits']]),
           '常数基准为同一窗口、同一六次 Huber-D⁴ IRLS；基准 k/σ/χ 分别为 '+str(u['constant_reference'])+'。q=1 的 k、σ 收敛，χ 因自由度约定不收敛到旧口径；q=2 的差异是模型零空间造成，非数值误差。',
           '### 3. 对数斜率 2/s，呼吸和噪声',
           '24 拍、噪声 0.04 N、真实末端 3000 N/m；160 组共同探测相位和共同噪声，q1/q2 使用 λₖ=.1，避免以开发结果回调单元检查。生成式使用指定的参考不变模型，属于模型匹配试验。',
           tab(['方法','端点偏差 N/m','经验标准差 N/m','平均 σ N/m'],[[m,f(a['endpoint_bias']),f(a['empirical_sd']),f(a['mean_sigma'])] for m,a in u['ramp'].items()]),
           tab(['无噪声斜坡','末端相对偏差','近平台比例','估计/真实中位斜率'],[[a['method'],f(a['endpoint_relative_bias']),f(a['near_flat_fraction']),f(a['relative_slope_median'])] for a in u['ramp_shape']]),
           '近平台定义为 |Δk̂|<0.1Δk真值。此处 q=1 是一阶差分 Huber，有限 δ 下小增量仍是二次惩罚，不能等同于纯 L1-TV；若没有形成平台，不能强行把端点偏差称为“台阶化”。完整无噪声轨迹保留在 tvk_unit.json。',
           'q=2 的斜坡偏差更小，同时方差显著增加；约四倍端点方差只是特定探测与窗口下的近似，当前有限窗口不能把它当常数定理。',
           '### 4. ×2/×4 台阶与咳嗽式 d 脉冲',
           '在 0.8 s 令 k 跳变，叠加中心 0.85 s、宽 35 ms、1.2 N 高斯脉冲和呼吸、共同噪声。恢复定义为连续五拍相对误差 <20%；最大低估在跳变后 0.5 s 内统计。此项是单个共同噪声实例，不能证明尾部安全。',
           tab(['倍数','方法','恢复 s','最大低估比例','前段 σ','低估时平均 σ','最大低估当拍 σ','跳变首拍 σ','低估时物理 2σ 覆盖'],[[a['ratio'],a['method'],f(a['recovery_s']),f(a['max_underestimate']),f(a['sigma_pre']),f(a['sigma_under']),f(a['sigma_at_max_under']),f(a['sigma_first_post']),f(a['under_2sigma_coverage'])] for a in u['steps']]),
           'σ 增大不等于安全：见 ×4 的低估覆盖；冻结 Hessian 无法完整反映突变模型偏差。',
           '### 5. 3/6/9 Hz 调制与 d 泄漏',
           '24 拍、30% 正弦调制，32 个共同调制相位，无噪声隔离结构泄漏。泄漏定义为 (d̂−d) 在真实调制响应 (k−mean k)Y 上的有符号最小二乘投影系数，两向量均去均值使其不依赖 d 的任意常数参考；不截断到 [0,1]。生成式为物理 kY，因而也检验指定模型的近似误差。',
           tab(['调制 Hz','方法','平均泄漏比例','泄漏 p10/p90','末端相对 RMSE'],[[a['f_Hz'],a['method'],f(a['leakage_mean']),str([round(x,3) for x in a['leakage_p10_p90']]),f(a['endpoint_relative_rmse'])] for a in u['modulation']]),
           '**未验证出硬性的约 7 Hz 可辨识界。** 14 Hz 是软截止，不是带限；0.192 s 窗口有明显端点效应，Huber 权重还依赖结构。20.8−fₖ<14 只描述下边带进入 d 的弱惩罚区，上边带仍在。即使 fₖ<7，也不能从截止差直接推出无泄漏或唯一可辨识。表中 3 Hz 已有泄漏，6→9 Hz 只是部分方法泄漏上升。稀疏 d 与周期幅值变化可被先验偏好区分，不是普适唯一性保证。',
           '### 6. 单进程计时',
           tab(['方法','估计器 p50/p99 μs','整个策略 p50/p99 μs','估计器调用数'],[[a['method'],str([round(x,1) for x in a['estimator_p50_p99_us']]),str([round(x,1) for x in a['policy_p50_p99_us']]),a['estimator_calls']] for a in timing['rows']]),
           '同一呼吸源日志，先前已完成导入/编译，去除首秒；估计器计时含 H 预处理和窗口构造，策略为完整调用。此为本机单进程原型计时，不是实时认证。',
           '## B. 模块级结果',
           '复用 baseline_compare.replay_job/source_observation：原 Huber-D4、D4 两份源日志各 168 场，共 336 场/配置。先使用已保存的真实 callback 观测恢复因果输入，再运行同一宿主的原门控；统计包含所有尝试和无效窗口。t≥1 s、真值>200；event 优先于 changing，changing 为 |dlogk/dt|>2/s 的 30 拍扩张，与 baseline 完全相同。误差用原真值列 21，不把它当控制输入。',
           '同源同方法的 '+str(len(identities))+' 个回放均核验命令复现，最大差为 '+str(max(identities))+' m。',
           'Huber-current 是当前默认；Huber-D4、D4、Kalman、MHE 使用第二轮选定配置。q/p 对照固定对应 q 的开发最佳 λ/δ/gate；改变 p 后未重新择优。这些是诊断变体，不参与闭环重新选参。相同配置的标签复用同一批回放，不伪装独立试验。',
           '实际回放配置（未列出的 TVK 通用参数采用上节默认）：\n\n```json\n'+__import__('json').dumps(module['chosen'],ensure_ascii=False,indent=1)+'\n```',
           '可用率=宿主接受数/全部符合分组的时间拍；绝对对数误差中位/p90 使用所有有限估计；2σ 覆盖沿用 baseline 的 |log(k̂/k)|≤2σ/k̂，分别给出接受窗口和全部有限窗口。带符号误差为 mean log(k̂/k)。低估>10% 用 k̂<0.9k；低估通过率=这些低估窗口中通过原门控的比例。上升低估率只在未扩张的实际 dlogk/dt>2/s 统计。',
           '滞后指标用上升段等效局部滞后 −log(k̂/k)/(dlogk/dt)，保留负值（超前/高估）；不把它当真实纯时间延迟。完整样本保留在 npz。自然突变簇恢复时间与删失沿用 baseline recovery，另表给出。',
           '### 双源分别统计：quiet / event / changing']
    heads=['源','方法','组','可用率','绝对log误差中位/p90','接受/全部 2σ覆盖','带符号误差','低估率','上升低估率','低估通过率','滞后p50/p90 s']
    def mrow(a):
        return [a['source'],a['method'],a['kind'],f(a['available']),f(a['abslog_med'])+'/'+f(a['abslog_p90']),f(a['coverage_accepted'])+'/'+f(a['coverage_all']),f(a['signed_log_mean']),f(a['under10_rate']),f(a['rising_under10_rate']),f(a['under_pass_rate']),str([round(x,3) for x in a['lag_p50_p90_s']]) if a['lag_p50_p90_s'] else '—']
    lines[-1:-1]=['### 合并双源的主要对比',tab(heads,[mrow(a) for a in ms if a['cond']=='all' and a['source']=='all' and a['method'] in ('Huber-current','TVK-selected','TVK-q2','D4','Kalman','MHE')])]
    lines.append(tab(heads,[mrow(a) for a in ms if a['cond']=='all' and a['source']!='all']))
    lines+=['### 合并双源、逐工况与分组',tab(['工况']+heads[1:],[[cn.get(a['cond'],a['cond'])]+mrow(a)[1:] for a in ms if a['source']=='all' and a['cond']!='all'])]
    recover=[]
    for m in module['chosen']:
        rr=[r for r in module['rows'] if r['method']==m];rec=[x for r in rr for x in r['recovery']];delay=[x['delay'] for x in rec if x['delay'] is not None]
        recover.append([m,len(rec),len(delay),len(rec)-len(delay),f(R.median(delay))])
    lines += ['### 自然突变后的接受估计恢复（连续五拍在20%以内）',tab(['方法','突变簇数','恢复数','右删失数','已恢复延迟中位 s'],recover),
              '完整源×工况×分组的计数、比例以及不接受窗口的误差都在 tvk_summary.json；未以海量相关窗口做独立样本显著性检验。',
              '## C. 闭环开发选参',
              '严格使用 d_compare.DEV 的四种(工况,组织)×两个 i=0/1，共八场，原开发场景/传感器种子与原评分 L=mean(log1p(J))。24 组在 tvk_plan.json 中于运行前固定：q∈{2,1}、p∈{4,3}、λₖ∈{.01,.1,1}、gate∈{.551,1}，δ 随 q 固定；σ_scale=.796，nw=24、fc=14 Hz，六次 IRLS。仅开发集按 L 最小选择。',
              tab(['编号','q','p','λₖ','δₖ','gate','平均 L'],[[j,c['q'],c['order'],c['lam_k'],c['delta_k'],c['gate'],f(dev['scores'][j])] for j,c in enumerate(dev['configs'])]),
              '相同八场开发集的既有基准：Huber-current L='+f(baseline_dev['scores']['Huber-D4']['default'])+'；重调Huber-D4='+f(baseline_dev['scores']['Huber-D4']['best'])+'；Kalman='+f(baseline_dev['scores']['Kalman']['best'])+'；MHE='+f(baseline_dev['scores']['MHE']['best'])+'。这些直接复用基线开发记录，不另增基线运行。',
              '运行预算：小样 '+f(runtime['pilot'])+' s，开发 '+f(runtime['development'])+' s，诊断闭环 '+f(runtime['diagnostic_test'])+' s，模块回放 '+f(runtime['module_replay'])+' s；四项累计 '+f(runtime['accounted_total']/60)+' min（不含编码、单元检查和单进程计时）。正式批次使用15个worker，BLAS单线程；单进程计时在批次结束后运行。',
              '选中配置：`'+str(dev['chosen'])+'`。总体最优与 q=1 最优可能为同一配置，只计算一次实际测试再复用标签；额外报告 q=2 的开发最佳，使默认物理先验也接受闭环检验。λₖ 最优在上边界，未扩展预算；不能声称已经找到 TVK 全局最优。',
              '## D. 诊断测试与配对差',
              '使用已看过的 base=574500、原 SENSOR_BASE=857000、原 168 场格子。**只作诊断，不作论文数字；未据此改参数。** Huber-current 与所有外部/内部基线逐行复用 baseline_test_r2.json；TVK 使用相同场景和传感器种子。簇 bootstrap 原封不动调用 baseline_report.paired：5000 次，按(工况,组织,i)聚类，整体42簇，每簇四个动态/相位版本；分工况3或9簇。',
              tab(['方法','L','ΔL对当前 [95%簇CI]','RMSE 中位/均值','事件 RMSE 均值','峰值>15N','最大力 N','长失联','failed'],[[m,f(stats[m]['L']),R.difftext(differences[m]['all']),f(stats[m]['rmse_med'])+'/'+f(stats[m]['rmse_mean']),f(stats[m]['event_mean']),stats[m]['over15'],f(stats[m]['peak']),stats[m]['long'],stats[m]['fail']] for m in methods]),
              tab(['TVK配置','对照','ΔL [95%簇CI]'],[[m,other,R.difftext(d)] for m,pp in pairwise.items() for other,d in pp.items()]),
              tab(['工况','方法','L','ΔL对当前 [95%簇CI]','峰值>15N','最大力 N','长失联'],[[cn.get(c,c),m,f(condition_stats[m][c]['L']),R.difftext(differences[m][c]),condition_stats[m][c]['over15'],f(condition_stats[m][c]['peak']),condition_stats[m][c]['long']] for c in conds for m in methods]),
              '饱和命令、非有限命令、>36.5 N、全部失败计数和非聚类/分层区间保存在 tvk_summary.json，未省去不利回合。',
              '## E. 结论',
              '当前诊断中，TVK-selected 相对 Huber-current 的 ΔL 为 '+R.difftext(differences['TVK-selected']['all'])+'，q=2 为 '+R.difftext(differences['TVK-q2']['all'])+'。与 Kalman/MHE 的整体 L 对照分别为 '+f(stats['TVK-selected']['L'])+' / '+f(stats['Kalman']['L'])+' / '+f(stats['MHE']['L'])+'；不能仅据均值宣布普适优越，尤其需一起看事件覆盖、上升低估通过率和尾部力。',
              'q=2 在模型匹配的平滑斜坡上确实减小端点偏差，符合分段线性先验的预期；但它无惩罚的斜率也放大静态窗口方差和噪声斜率，有限窗口中并不保证闭环更合适。q=1 引入斜率收缩，可能付出斜坡偏差，开发集却选中了它，说明当前控制使用方式与噪声环境需要偏差—方差平衡。不能把跨肋骨的物理斜坡直接推成 q=2 综合最优。',
              '在此已见诊断格，q=1 相对 Kalman 的配对改善区间不跨零；相对 MHE 的簇区间跨零，相对当前 Huber-D⁴ 也未显示改善。因此不能得出 TVK 更适合替代 Huber-D⁴；也不能把对 MHE 的非聚类区间替代预定的簇区间来宣布显著优势。q=2 在完整宿主上明显不适合当前设置。',
              'TVK 允许 k 在窗内变化，可降低常数模型的平均化偏差；与此同时，χ 可把模型失配解释掉，σ 的局部曲率也不覆盖所有跳变偏差。×4 台阶已直接展示风险，闭环仍受前馈对刚度低估的敏感性限制。本次保留为候选，不据模块精度或单一工况收益自动替换当前默认。',
              '## F. 不确定或可能不公平之处',
              '1. 验收要求存在两处不可同时满足的数学条件（q=2常数极限、χ自由度）；本次原样报告失败，不改变目标掩盖。若 Claude 要求修订模型，应明确新增斜率惩罚或统一χ口径，再开展新的开发选择，不能回用此诊断格调参。',
              '2. 六次 IRLS 是固定计算预算的近似解，稠密核验验证的是冻结权重线性系统，而非所有窗口已达到鲁棒目标最小值。Huber 对噪声和突变的区分是先验偏好；模型失配、非正估计和局部曲率校准问题仍在。',
              '3. TVK 的相对尺度 K 来自当前窗二次估计，固定后目标凸；其数值及 δ 影响先验强度。q1/q2 的 δ 不同以对应增量/曲率不同物理单位，不能解释成相同数学惩罚强度下的纯阶数消融。p3变体只在回放中固定其他设置，S严格固定为二次D4、14Hz平滑器。最终复核纠正了p3平滑器误随d阶数改变的实现，受影响的12组开发配置与两组双源回放均已重跑，旧p3数值不再作为证据。',
              '4. 模型匹配斜坡、原始 kY 台阶和调制回答不同问题。仅有160组斜坡噪声、一个台阶噪声实例、32个调制相位；约7Hz只是频谱分离的启发，不是本次通过的可辨识定理。',
              '5. baseline 的对手此前有38–72组开发预算，本次TVK24组；当前默认固定，另列旧24组重调Huber作参照。TVK的λ最优触边，未扩大搜索。双源为闭环产生的同分布日志，不能等价于任意控制轨迹或真实人体。',
              '6. 原真值列21未完全计入编码器—接触增益中的探头惯性与伺服近似偏差，因此共享约正偏差仍影响log误差/覆盖。未修正真值以偏袒TVK。门控后的coverage有选择效应，因此同时给出全部有限窗口coverage。沿用原 _clip_k：负或越界的有限端点会截断至既有上下界，σ仍是未约束Hessian的曲率，非约束后严格分布。',
              '7. 低估率与滞后只在有限估计上计算；可用率分母包括原回放全部符合条件的时间拍。局部等效滞后可为负，不代表实际预测超前。自然恢复有右删失，表中中位只对已恢复簇，不隐藏删失数。',
              '8. common._build 的最后可用回归行与当前 callback 相差一拍（8 ms），所有方法沿用同一口径；TVK输出该窗口末端，不额外外推到callback时刻。168格子此前已见且各方法复用相同物理场景，置信区间仅反映此仿真分布。报告不作论文推广结论，也不把局部Hessian σ当安全证书。',
              '复现命令：依次运行 `python -X utf8 experiments/rossi_failure_study/tvk_compare.py --stage unit|pilot|dev|test|replay|timing|report`（每次指定一个stage）。test/replay可续跑；dev显式重跑会覆盖其同名结果。所有主要配置、统计和耗时均在 results/tvk_*.json。']
    lookup={(a['method'],a['kind']):a for a in ms if a['source']=='all' and a['cond']=='all'}
    a,b=lookup['Huber-current','event'],lookup['TVK-selected','event']
    c,d=lookup['Huber-current','changing'],lookup['TVK-selected','changing']
    at=lines.index('## E. 结论')+1
    lines.insert(at,'模块主结果：q=1 事件段绝对log误差中位 '+f(a['abslog_med'])+'→'+f(b['abslog_med'])+'，接受窗口2σ覆盖 '+f(a['coverage_accepted'])+'→'+f(b['coverage_accepted'])+'，但可用率 '+f(a['available'])+'→'+f(b['available'])+'。刚度变化段误差 '+f(c['abslog_med'])+'→'+f(d['abslog_med'])+'，上升段低估>10%比例 '+f(c['rising_under10_rate'])+'→'+f(d['rising_under10_rate'])+'，等效滞后p90 '+f(c['lag_p50_p90_s'][1])+'→'+f(d['lag_p50_p90_s'][1])+' s。改进与退化同时存在。')
    REPORT.parent.mkdir(parents=True,exist_ok=True);REPORT.write_text('\n\n'.join(lines)+'\n',encoding='utf-8')
    print('REPORT',REPORT,'SUMMARY',[(m,stats[m]['L'],differences[m]['all']) for m in dev['chosen']],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True,choices=['unit','pilot','dev','test','replay','timing','report','refresh_p3']);args=ap.parse_args()
    globals()[args.stage]()
