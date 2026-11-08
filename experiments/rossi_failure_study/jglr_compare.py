"""JGLR checks; reuse baseline_compare and tvk_compare, no controller changes.

Run --stage unit|pilot|dev|test|replay|timing|report. Workers=15, BLAS=1.
Synthetic seeds/design and 12 development configurations are fixed in code.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_k]='1'
import argparse
import time
import math
import json
from pathlib import Path
from itertools import product
from types import SimpleNamespace
import numpy as np
from scipy.stats import chi2, poisson
import tvk_compare as V
import baseline_compare as B
import baseline_report as R
import jglr as J

REPORT=Path('D:/pip_env/temp/claude/8fbcb2674095e6be/841f52cd-14a8-4db4-b56f-36be6ffcc337/scratchpad/jglr_report.md')
BASE_BUILD=V.ORIGINAL_BUILD
CONDITIONS=('rib','cough','combo','breath','flinch')


def build(method,cfg,sc):
    if method.startswith('JGLR'):return J.build(sc,cfg)
    if method.startswith('TVK'):return V.T.build(sc,cfg)
    return BASE_BUILD(method,cfg,sc)


def configs():
    a=[dict(q_k=q,L=l,alpha=alpha,impulse=True)
       for q,l,alpha in product((0.,.002),(24,96),(1e-4,1e-3))]
    a += [dict(q_k=q,L=48,alpha=1e-4,impulse=imp)
          for q,imp in product((.0005,.002),(False,True))]
    return a


def execute(job):
    B.build=build
    return B.execute(job)


def solver(**kw):
    e=J.JGLRRidge(nw=24,fc=14.,**kw)
    e.reset(SimpleNamespace(Tc=.008,params={}))
    return e


def simulate_case(task):
    ratio,seed=task;rng=np.random.default_rng(seed);n=384;theta=128
    noise=.04;A=.167;k0=3000.;dt=.008;t=np.arange(n)*dt
    phase=rng.uniform(0,2*np.pi);y=A/k0*np.sin(2*np.pi*20.8*t+phase)
    k=np.where(np.arange(n)<theta,k0,k0*ratio)
    rho=solver().rho;f,g=J.dynamics();dstate=np.zeros(5);dd=[]
    for _ in range(n):
        dstate=f@dstate+g*rng.normal(0,noise/math.sqrt(rho));dd.append(dstate[1])
    # Exact D4 process plus a slow breathing mean; a k jump accompanies d step.
    z=k*y+np.array(dd)+.4*np.sin(2*np.pi*.3*t)+.25*(np.arange(n)>=theta)+rng.normal(0,noise,n)
    blind=J.JumpFilter(q_k=0.,L=96,alpha=1e-4)
    oracle=J.JumpFilter(q_k=0.,detect=False)
    ep=np.zeros((5,2));ep[0,0]=1.;ep[1,1]=1.
    info=np.zeros((2,2));b=np.zeros(2);snap=[];blind_snap=[];cov_min=0.
    for i in range(n):
        out=blind.step(y[i],z[i]);oracle.step(y[i],z[i])
        if i<theta:continue
        h=np.array([3000*y[i],1.,0.,0.,0.])
        # Reconstruct nominal gain at this sample from P after update.
        gain=oracle.P@h/noise**2
        if i>theta:ep=f@ep
        sig=h@ep;ep=ep-np.outer(gain,sig)
        info+=np.outer(sig,sig)/oracle.last['variance']
        b+=sig*oracle.last['innovation']/oracle.last['variance']
        m=i-theta+1
        if m in (6,12,24,48,96):
            inv=np.linalg.inv(info);xc=oracle.x+ep@inv@b;pc=oracle.P+ep@inv@ep.T
            snap.append(dict(m=m,error=float(3000*xc[0]-k[i]),sigma=float(3000*math.sqrt(pc[0,0]))))
            blind_snap.append(dict(m=m,error=float(out[0]-k[i]),sigma=out[1],
                                   covered=bool(abs(out[0]-k[i])<=2*out[1])))
            cov_min=min(cov_min,float(np.linalg.eigvalsh(blind.P).min()))
    detections=[e for e in blind.events if e['index']>=theta and e['kind']=='k_d_step']
    false=[e for e in blind.events if e['index']<theta]
    return dict(ratio=ratio,seed=seed,delay=None if not detections else detections[0]['index']-theta+1,
                first_detection=None if not detections else detections[0],false_pre=len(false),
                oracle=snap,blind=blind_snap,events=blind.events,cov_min=cov_min)


def null_case(seed):
    rng=np.random.default_rng(seed);f=J.JumpFilter(q_k=0.,L=96,alpha=1e-4)
    F,g=J.dynamics();x=np.array([1.,0.,0.,0.,0.]);n=5000
    for i in range(n):
        x=F@x+g*rng.normal(0,.04/math.sqrt(f.rho))
        y=.167/3000*math.sin(2*math.pi*20.8*.008*i+.31)
        z=3000*x[0]*y+x[1]+rng.normal(0,.04)
        f.step(y,z)
    return dict(seed=seed,n=n-8,alarms=len(f.events),events=f.events)


def signature_check():
    rows=[];rng=np.random.default_rng(941)
    for phase in np.linspace(0,2*np.pi,12,endpoint=False):
        f=J.JumpFilter(q_k=0.,detect=False);F,_=J.dynamics()
        ep=np.zeros((5,3));ep[0,0]=1.;ep[1,1]=1.;info=np.zeros((3,3));raw=[]
        for i in range(224):
            yy=.167/3000*np.sin(2*np.pi*20.8*.008*i+phase)
            f.step(yy,3000*yy)
            if i<128:continue
            h=np.array([3000*yy,1.,0.,0.,0.]);gain=f.P@h/.04**2
            if i>128:ep=F@ep
            sig=h@ep
            if i==128:sig[2]+=1.
            ep-=np.outer(gain,sig);info+=np.outer(sig,sig)/f.last['variance']
            raw.append([3000*yy,1.,float(i==128)])
            m=i-128+1
            if m not in (1,2,3,4,6,12,24,48,96):continue
            rawgram=np.array(raw).T@np.array(raw)
            angle=lambda gram,a,b:float(np.degrees(np.arccos(np.clip(abs(gram[a,b])/math.sqrt(max(gram[a,a]*gram[b,b],1e-300)),0,1))))
            out=dict(m=m,phase=float(phase),angle_raw_step=angle(rawgram,0,1),
                     angle_white_step=angle(info,0,1),angle_white_impulse=angle(info,0,2))
            if m>=2:
                vals,vec=np.linalg.eigh(info);root=vec@np.diag(np.sqrt(np.maximum(vals,0)))
                rand=rng.normal(size=(5000,3))@root.T
                schur=info[0,0]-info[0,1]**2/info[1,1]
                for label,amp in [('k30_d',np.array([.3,.25,0.])),('k2_d',np.array([1.,.25,0.])),('d',np.array([0.,.25,0.])),('impulse',np.array([0.,0.,.25]))]:
                    b=rand+info@amp
                    sd=b[:,1]**2/info[1,1]
                    sk=(b[:,0]-info[0,1]*b[:,1]/info[1,1])**2/max(schur,1e-300)
                    si=b[:,2]**2/max(info[2,2],1e-300)
                    scores=np.column_stack((np.zeros(5000),sd-f.h1,si-f.h1,
                                            np.where(sk>f.h1,sd+sk-f.h2,-np.inf)))
                    classification=np.argmax(scores,axis=1)
                    out[label]={key:float(np.mean(classification==j)) for j,key in enumerate(('none','d_step','d_impulse','k_d_step'))}
            rows.append(out)
    return rows


def ramp_case(task):
    rate,seed=task;rng=np.random.default_rng(seed);n=256;t=np.arange(n)*.008
    theta=128;k=np.where(np.arange(n)<theta,3000.,3000*np.exp(rate*(t-t[theta])))
    y=.167/3000*np.sin(2*np.pi*20.8*t+rng.uniform(0,2*np.pi))
    z=k*y+.4*np.sin(2*np.pi*.3*t+.4)+rng.normal(0,.04,n)
    e=J.JumpFilter(q_k=0.,L=96,alpha=1e-4)
    for yy,zz in zip(y,z):a=e.step(yy,zz)
    q1=V.solver(q=1,lam=1.,delta=.02);base=V.constant(q1,y[-24:],z[-24:],.04)
    tvk=q1.solve_window(y[-24:],z[-24:],.04)
    return dict(rate=rate,seed=seed,truth=float(k[-1]),
                estimates={'JGLR':dict(k=a[0],sigma=a[1]),
                           'Huber-current':dict(k=base['k'],sigma=.796*base['sigma']),
                           'TVK-q1':dict(k=float(tvk['k'][-1]),sigma=.796*tvk['sigma'])})


def unit():
    tic=time.perf_counter();rng=np.random.default_rng(1651);m=24;dt=.008;t=np.arange(m)*dt
    y=.167/3000*np.sin(2*np.pi*20.8*t+.3);z=3000*y+.4*np.sin(2*np.pi*.3*t)+rng.normal(0,.04,m)
    e=solver();ref=J.batch_arrival(y,z,.04,e.rho)
    f=J.JumpFilter(q_k=0.,detect=False)
    for yy,zz in zip(y,z):a=f.step(yy,zz)
    from ridge_estimator import RidgeProjection
    rp=RidgeProjection(nw=24);rp.reset(SimpleNamespace(Tc=.008,params={}))
    metric,_=rp._matrix(m);kk=float(y@metric@z/(y@metric@y));ss=.04/math.sqrt(y@metric@y)
    equivalence=dict(k_relative=abs(a[0]-kk)/abs(kk),sigma_relative=abs(a[1]-ss)/ss,
                     state_max=float(np.max(abs(f.x-ref[0]))),cov_max=float(np.max(abs(f.P-ref[1]))))
    assert max(equivalence['k_relative'],equivalence['sigma_relative'])<1e-7,equivalence
    weights=e._probe_weights(24);offset_error=abs(weights@(y+3.7-y[-1]-3.7)-weights@(y-y[-1]))
    carrier=weights@np.exp(2j*np.pi*20.8*dt*(np.arange(m)-m+1))
    assert offset_error<1e-14 and abs(carrier-1)<1e-12
    # Host reference invariance: invoke the actual estimator with injected rows.
    policy=SimpleNamespace(Tc=.008,params=dict(noise=.04,enc_y=0.),ndel=0,ns=8,eb=.5,afc=.2)
    host=[];orig=J.common._build
    try:
        for offset in (0.,3.7):
            est=solver(q_k=0.);log=np.zeros((160,26));log[:,0]=np.arange(160)*dt
            def injected(log,j0,mm,*args):
                tt=np.arange(j0,j0+mm)*dt;yy=.167/3000*np.sin(2*np.pi*20.8*tt)
                return None,yy+offset,3000*yy+.4*np.sin(2*np.pi*.3*tt),np.ones(mm,bool)
            J.common._build=injected
            outputs=[est.estimate(SimpleNamespace(log=log),n,policy) for n in range(10,160)]
            host.append(np.array(outputs));
    finally:J.common._build=orig
    mask=np.isfinite(host[0])&np.isfinite(host[1]);host_error=float(np.max(abs(host[0][mask]-host[1][mask])))
    assert host_error<1e-4,host_error
    B.dump('jglr_unit_basic.json',dict(equivalence=equivalence,reference_regressor_abs=offset_error,
           host_reference_max_abs=host_error,carrier_real=float(carrier.real),carrier_imag=float(carrier.imag),
           thresholds=dict(h1=f.h1,h2=f.h2),seconds=time.perf_counter()-tic))
    steps=B.run([(ratio,17710+s) for ratio in (1.3,2.,4.) for s in range(120)],worker=simulate_case)
    B.dump('jglr_synthetic_steps.json',dict(rows=steps,n_per_ratio=120,memory=96,alpha=1e-4,q_k=0.,seconds=time.perf_counter()-tic))
    null=B.run(list(range(19610,19630)),worker=null_case)
    B.dump('jglr_synthetic_null.json',dict(rows=null,total_samples=sum(r['n'] for r in null),alarms=sum(r['alarms'] for r in null),alpha=1e-4))
    angles=signature_check();B.dump('jglr_signatures.json',dict(rows=angles,replicates_per_phase=5000,n_phases=12))
    ramps=B.run([(rate,21710+s) for rate in (.5,1.,2.) for s in range(80)],worker=ramp_case)
    B.dump('jglr_synthetic_ramps.json',dict(rows=ramps,n_per_rate=80))
    B.dump('jglr_unit_runtime.json',dict(seconds=time.perf_counter()-tic))
    print('UNIT',equivalence,'host_reference',host_error,'null',sum(r['alarms'] for r in null),flush=True)


def pilot():
    cell=('breath','relaxed','neural_white',21.,0)
    tic=time.perf_counter()
    jobs=[('JGLR',configs()[j],cell,B.D.BASE,0.,False,True) for j in (0,7)]
    jobs.append(('Huber-current',{},cell,B.D.BASE,0.,False,True))
    rows=B.run(jobs,worker=execute)
    B.dump('jglr_pilot.json',dict(rows=rows,seconds=time.perf_counter()-tic,workers=15,blas_threads=1))
    print('PILOT',[(r['method'],r['seconds'],r['J']) for r in rows],flush=True)


def dev():
    cfgs=configs();rowspec=[(c,t,'neural_white',21.,i) for c,t in B.D.DEV for i in range(2)]
    thresholds=[dict(c,h1=J.JumpFilter(**c).h1,h2=J.JumpFilter(**c).h2) for c in cfgs]
    B.dump('jglr_plan.json',dict(configs=thresholds,criterion='original DEV8 mean log1p(J)',workers=15,blas_threads=1))
    tic=time.perf_counter();rows=B.run([('JGLR',c,cell,B.D.BASE,0.,False,True) for c in cfgs for cell in rowspec],worker=execute)
    scores=[B.L(rows[j*8:j*8+8]) for j in range(len(cfgs))];best=int(np.argmin(scores))
    B.dump('jglr_dev.json',dict(configs=cfgs,scores=scores,index=best,chosen={'JGLR':cfgs[best]},rows=rows,seconds=time.perf_counter()-tic))
    print('DEV',best,scores[best],cfgs[best],flush=True)


def test():
    cfg=B.load('jglr_dev.json')['chosen']['JGLR'];tic=time.perf_counter()
    rows=B.run([('JGLR',cfg,cell,B.TEST_BASE,0.,False,False) for cell in B.cells()],worker=execute)
    B.dump('jglr_test.json',dict(chosen=cfg,rows=rows,seconds=time.perf_counter()-tic,diagnostic_only=True))
    print('TEST',B.L(rows),flush=True)


def replay_job(job):
    holder={}
    def capture(method,cfg,sc):
        pol=build(method,cfg,sc);holder['pol']=pol;return pol
    V.build=capture
    result=V.replay_job(job,log_folder='jglr_logs')
    est=holder['pol'].estimator;events=est.events
    result['detections']=events
    # Match to the SAME natural transition clusters used for recovery.
    changes=result['recovery'];matches=[]
    for event in changes:
        valid=[d for d in events if d['kind']=='k_d_step' and event['start']<=d['callback_time']<=event['end']+1.5]
        matches.append(dict(**event,detection_delay_from_start=None if not valid else valid[0]['callback_time']-event['start'],
                            detection_delay_from_end=None if not valid else valid[0]['callback_time']-event['end']))
    result['detection_matches']=matches
    result['detection_counts']={kind:sum(e['kind']==kind for e in events) for kind in ('k_d_step','d_step','d_impulse')}
    # Preserve actual causal regressors for information-limit log checks.
    with np.load(B.OUT/result['samples_file']) as z:data={k:z[k] for k in z.files}
    data['regressor_trace']=np.array(est.regressor_trace);np.savez_compressed(B.OUT/result['samples_file'],**data)
    return result


def replay():
    cfg=B.load('jglr_dev.json')['chosen']['JGLR']
    sources=[r for r in B.load('baseline_test.json')['rows'] if r['method'] in ('Huber-D4','D4') and r['cond'] in CONDITIONS and 'log_file' in r]
    # These exact source files and baseline configurations have already been replayed.
    old=B.load('tvk_module.json');reuse=[r for r in old['rows'] if r['method'] in ('Huber-current','TVK-q1') and r['cond'] in CONDITIONS]
    tic=time.perf_counter();rows=B.run([(s,'JGLR',cfg) for s in sources],worker=replay_job)
    B.dump('jglr_module.json',dict(chosen={'JGLR':cfg,'Huber-current':{},'TVK-q1':old['chosen']['TVK-q1']},
                                 rows=rows+reuse,seconds=time.perf_counter()-tic,reused_rows=len(reuse),new_rows=len(rows)))
    print('REPLAY',len(rows),'new',len(reuse),'reused',flush=True)


def timing():
    old_load=V.B.load;old_build=V.build;old_dump=V.B.dump
    chosen={'JGLR':B.load('jglr_dev.json')['chosen']['JGLR']}
    V.build=build
    V.B.load=lambda name:dict(chosen=chosen) if name=='tvk_dev.json' else old_load(name)
    V.B.dump=lambda name,obj:old_dump('jglr_timing.json' if name=='tvk_timing.json' else name,obj)
    try:V.timing()
    finally:V.B.load=old_load;V.B.dump=old_dump;V.build=old_build


def report():
    import jglr_report
    jglr_report.main()


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True,choices=['unit','pilot','dev','test','replay','timing','report'])
    globals()[ap.parse_args().stage]()
