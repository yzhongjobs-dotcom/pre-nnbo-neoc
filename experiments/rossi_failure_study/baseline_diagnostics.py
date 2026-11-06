"""Post-selection diagnostics: synthetic stiffness steps and warm serial timing.

These results never select parameters. Natural scan transitions are reported
separately. The synthetic equation is generated independently at 1 kHz.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
import math
import time
from types import SimpleNamespace

import numpy as np

import baseline_compare as B
import id_baselines as I
import physics_plant3 as P3


def synthetic(seed=18,moving=False):
    dt=.008;t=np.arange(751)*dt;step=3.
    k=np.where(t<step,2000.,8000.)
    y=.004+.0001*np.sin(2*np.pi*20.8*t)
    surface=.001*np.sin(2*np.pi*.3*t) if moving else np.zeros(len(t))
    a=1-math.exp(-2*np.pi*30*.001);eb=math.exp(-.001/.019);f=8.;rng=np.random.default_rng(seed)
    log=np.full((len(t),26),np.nan);log[:,0]=t;log[:,1]=0.;log[:,2]=y;log[:,11]=y
    log[:,12:14]=0.;log[0,5]=f
    for j in range(len(t)-1):
        for mm in range(1,9):
            fr=(1-eb**(mm-1))/(1-eb**8)
            yy=y[j]+fr*(y[j+1]-y[j]);tt=t[j]+(mm-1)*.001
            kk=2000. if tt<step else 8000.
            # Inverse-filter equation noise is small but nonzero; it is not
            # drawn from the same prior used by any candidate estimator.
            ss=surface[j]+(mm-1)/8*(surface[j+1]-surface[j])
            f+=a*(kk*(yy-ss)+rng.normal(0,.06)-f)
        log[j+1,5]=f
    return log,k,step


def diagnostic(method,cfg,moving):
    log,truth,step=synthetic(moving=moving)
    # Only known settings below are used; no plant state, geometry or truth.
    sc=dict(FD=8.,geometry={'radius':.03},tool_mass=0.,sensor_delay_ms=0.,initial_contact_observed=0.)
    pol=B.build(method,cfg,sc);pol.reset(dict(P3.PLATFORM))
    raw=[];original=pol.estimator.estimate
    def record(o,n,p):
        tic=time.perf_counter();kh,sig,chi=original(o,n,p);took=time.perf_counter()-tic
        raw.append((n-1,kh,sig,chi,took));return kh,sig,chi
    pol.estimator.estimate=record
    for n in range(1,len(log)+1):
        # Estimator+host gate/phase only: fictitious fixture commands should
        # not update the robot-servo identifier in this ideal equation check.
        pol._stiffness_est(P3.Obs(log,n),n)
    accepted={round(v[0],7):v[-1] for v in pol.accept_trace}
    kh=np.full(len(log),np.nan);ok=np.zeros(len(log),bool);times=[]
    for j,k,sig,chi,took in raw:
        kh[j]=k;ok[j]=accepted.get(round(float(log[j,0]),7),False)
        if log[j,0]>=1.:times.append(took)
    match=np.isfinite(kh)&(np.abs(np.log(np.maximum(kh,1.)/truth))<=math.log(1.2))
    delays={}
    for label,valid in [('finite',np.isfinite(kh)),('accepted',ok)]:
        delay=None
        for j in np.flatnonzero(log[:,0]>=step):
            if j+5<=len(log) and np.all(match[j:j+5]&valid[j:j+5]):
                delay=float(log[j,0]-step);break
        delays[label]=delay
    warm=log[:,0]>=1.;pre=warm&(log[:,0]<step)
    before_step=np.flatnonzero(log[:,0]<step)[-5:]
    return dict(method=method,moving=moving,step_time=step,step_ratio=4.,delay=delays,censored=delays['accepted'] is None,
                pre_step_accurate=bool(np.all(match[before_step])),pre_step_accepted=bool(np.all(ok[before_step])),
                before_error=float(np.nanmedian(np.abs(np.log(kh[pre]/truth[pre])))),
                after_error=float(np.nanmedian(np.abs(np.log(kh[log[:,0]>5.]/truth[log[:,0]>5.])))),
                available_before=float(np.mean(ok[pre])),available_after=float(np.mean(ok[log[:,0]>5.])),
                serial_estimator_p50_us=float(np.median(times)*1e6),serial_estimator_p99_us=float(np.percentile(times,99)*1e6))


def selected_cases(saved):
    """Frozen-config illustrative failures, never a new selection criterion."""
    current={B.key(r):r for r in saved['rows'] if r['method']=='Huber-current'}
    rossi=[r for r in saved['rows'] if r['method']=='Rossi-control' and r['cond']=='rib']
    best=min(rossi,key=lambda r:math.log1p(r['J'])-math.log1p(current[B.key(r)]['J']))
    worst=max(current.values(),key=lambda r:r['J'])
    requests=[best,current[B.key(best)],worst,
              next(r for r in saved['rows'] if r['method']=='Rossi-control' and B.key(r)==B.key(worst))]
    for method in ('Huber-D4','DO-MPAC','Xue-const'):
        requests.append(max((r for r in saved['rows'] if r['method']==method),key=lambda r:r['peak'] or 0.))
    folder=B.OUT/'baseline_cases';folder.mkdir(exist_ok=True);rows=[];seen=set()
    columns=['time','compensated_force','requested_command','k_good','k_est','sigma','chi','accepted','phase',
             'eso_e','eso_edot','eso_disturbance','adaptive_compensation','command_velocity','rossi_khat']
    for source in requests:
        method=source['method'];cell=B.key(source)
        if (method,cell) in seen:continue
        seen.add((method,cell));sc=B.scene(cell,B.TEST_BASE,0.)
        if method.startswith('Xue'):
            sc['initial_contact_observed']=B.load('baseline_calibration.json')['positions'][B.calkey(cell,False)]
        pol=B.build(method,saved['chosen'][method],sc);states=[];original=pol.step
        def recorded(o):
            u=original(o);t=float(o.log[o.n-1,0]);trace=getattr(pol,'accept_trace',[])
            last=trace[-1][1:] if trace and abs(trace[-1][0]-t)<1e-7 else [np.nan]*4
            z=np.asarray(getattr(pol,'z',None));z=z if z.shape==(3,) else np.full(3,np.nan)
            states.append([t,float(o.log[o.n-1,5]),u,getattr(pol,'k_good',np.nan),*last,
                           getattr(pol,'phi',np.nan),*z,getattr(pol,'comp',np.nan),
                           getattr(pol,'velocity',np.nan),getattr(pol,'k_hat',np.nan)])
            return u
        pol.step=recorded
        log=P3.episode(sc,pol,B.sensor_seed(cell));metric=B.M.metric(log,pol,sc)
        assert abs(metric['J']-source['J'])<1e-9,(method,cell,metric['J'],source['J'])
        state=np.asarray(states,float);warm=log[:,0]>=1.;j=int(np.argmax(np.where(warm,log[:,18],-np.inf)))
        near=(log[:,0]>=log[j,0]-.25)&(log[:,0]<=log[j,0]);ratios=state[near,3]/np.maximum(log[near,21],1.)
        ratio=float(np.nanmedian(ratios)) if np.isfinite(ratios).any() else None
        path=folder/f"{method}_{cell[0]}_{cell[2]}_{int(cell[3])}_{cell[4]}.npz"
        np.savez_compressed(path,log=log,state=state,state_columns=columns)
        inner=getattr(getattr(pol,'estimator',None),'inner',None);extra={}
        if isinstance(inner,I.XueHC):extra=dict(invalid_depth_or_contact=inner.invalid_depth,updates=inner.updates,final_theta=inner.theta.tolist())
        rows.append(dict(method=method,cell=cell,J=metric['J'],peak=metric['peak'],peak_time=float(log[j,0]),
                         peak_true_gain=float(log[j,21]),preceding_k_good_to_true_gain=ratio,
                         J_reproduction_error=float(metric['J']-source['J']),diagnostic=extra,
                         file=str(path.relative_to(B.OUT)).replace('\\','/')))
    return rows


if __name__=='__main__':
    started=time.perf_counter()
    chosen=B.load('baseline_dev.json')['chosen'];rows=[]
    # One untimed warm-up removes compilation; then deterministic diagnostic
    # repeats are for timing only, no parameter search or statistical claims.
    diagnostic('D4',chosen['D4'],False)
    for method in B.ID_METHODS:
        for moving in (False,True):rows.append(diagnostic(method,chosen[method],moving))
    serial=[]
    if (B.OUT/'baseline_test.json').exists():
        saved=B.load('baseline_test.json');source=next(r for r in saved['rows'] if r['method']=='Huber-D4' and r['cond']=='breath')
        z=np.load(B.OUT/source['log_file']);log=np.full((len(z['log']),26),np.nan);log[:,z['cols']]=z['log'];observed=z['observed'];z.close()
        sc=B.scene(B.key(source),B.TEST_BASE,0.)
        sc['initial_contact_observed']=B.load('baseline_calibration.json')['positions'][B.calkey(B.key(source),False)]
        for method,cfg in saved['chosen'].items():
            pol=B.build(method,cfg,sc);pol.reset(dict(P3.PLATFORM,tau_a=sc['tau_a_hat']))
            for n in range(1,len(log)+1):pol(B.source_observation(log,observed,n))
            calls=np.array(pol.t_call)[125:]*1e6
            serial.append(dict(method=method,p50_us=float(np.median(calls)),p99_us=float(np.percentile(calls,99))))
    cases=selected_cases(saved) if (B.OUT/'baseline_test.json').exists() else []
    B.dump('baseline_diagnostics.json',dict(rows=rows,serial_full_calls=serial,cases=cases,seconds=time.perf_counter()-started,
                                         note='post-selection independent equation, warm serial fixed-source replay and illustrative frozen-config failures; not tuning or real-time certification'))
    for r in rows:print(r)
