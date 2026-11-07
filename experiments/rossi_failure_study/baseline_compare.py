"""Equal-budget development, fresh paired test grid and two-source replay.

Run --pilot first, then --stage dev/test/replay/report. Workers=15, BLAS=1.
Every stage is restartable from its own JSON. No git operations.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_key]='1'

import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import id_baselines as I
import ctrl_baselines as K
import method3 as M
import d_compare as D
import physics_plant3 as P3
import scenes3 as S3

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
TEST_BASE=574500  # fresh scene AND sensor seeds; old checks used 520500/900+i
SENSOR_BASE=857000
ID_METHODS=('Huber-D4','D4','P5','Rossi-ID','Kalman','MHE','Xue-const','Xue-chi')
CTRL_METHODS=('Rossi-control','DO-MPAC','Xing-DAHIC')
METHODS=ID_METHODS+CTRL_METHODS
OBS_COLS=[0,1,2,5,11,12,13]
SAVE_COLS=OBS_COLS+[3,6,18,20,21,22]


def configurations(method):
    if method in ('Huber-D4','D4','P5'):
        structs=([dict(nw=24,fc=14.),dict(nw=32,fc=14.),dict(nw=24,fc=10.)]
                 if method!='P5' else [dict(nw=24,degree=5),dict(nw=16,degree=5),dict(nw=16,degree=3)])
        return [dict(s,gate=g,sigma_scale=a) for s,g,a in product(structs,(.551,.2,1.,2.),(.796,1.5))]
    if method=='Rossi-ID':
        ss=[dict(order=o,half=h,memory=m) for o,h,m in product((1,2),(4.,8.),(.1,.4))]
        ss += [dict(filter='r2',Q=q,memory=m) for q,m in product((2.,4.),(.1,.4))]
        # Put the documented existing default first, without adding a 25th trial.
        ss.sort(key=lambda s:s!=dict(order=2,half=4.,memory=.4))
        return [dict(s,gate=g) for s,g in product(ss,(.551,1.))]
    if method=='Kalman':
        return [dict(q_k=k,q_d=d,d_order=o,gate=g) for k,d,o,g in product((300.,3000.),(.03,.3,3.),(2,1),(.551,1.))]
    if method=='MHE':
        ss=[(3000.,.3),(0.,0.),(300.,.03),(300.,.3),(3000.,.03),(3000.,3.)]
        return [dict(N=n,q_k=k,q_d=d,gate=g) for n,(k,d),g in product((24,16),ss,(.551,1.))]
    if method.startswith('Xue'):
        return [dict(delta=d,gain=g,initial_k=k,initial_n=n) for d,g,k,n in product((.05,.005,.5),(1.,.3),(3000.,30000.),(1.26,1.))]
    if method=='Rossi-control':
        return [dict(K_I=g,A_E=a,lam=math.exp(-.008/m)) for g,a,m in product((5.,15.,30.,60.),(2e-5,5e-5,1e-4),(.1,.4))]
    if method=='DO-MPAC':
        return [dict(eso_bandwidth=w,q_over_r=q,k_h0=k,Np=n) for w,q,k,n in product((8.,20.,50.),(1e4,1e6),(1500.,15000.),(15,40))]
    if method=='Xing-DAHIC':
        return [dict(alpha=a,m=m,b=b,k1=k1,k2=k2) for a,(m,b),(k1,k2) in product((1.,.99,.9),((1.,100.),(3.5,1000.),(1.,300.),(3.5,3000.)),((1.,1.),(.1,.01)))]
    raise ValueError(method)


def build(method,cfg,sc):
    if method=='Huber-current':method='Huber-D4'
    if method.startswith('Xue'):
        cfg=dict(cfg,contact_position=sc.get('initial_contact_observed'))
    return I.build(method,cfg,sc) if method in ID_METHODS else K.build(method,cfg,sc)


def cells():
    # Exactly the union of d_prior_check standard + "more": slow cells 3
    # seeds, fast-change cells 9 seeds -> 24+144=168 per method.
    return [(c,t,dy,ph,i) for c,t in D.TEST for dy in ('neural_white','viscous')
            for ph in (14.,28.) for i in range(3 if c in ('breath','shift_fast') else 9)]


def scene(cell,base,start_gap):
    c,t,dy,ph,i=cell
    sc=S3.cell_scene3(c,t,'nominal',i,base=base,overrides=dict(slide='long',force=8.,dynamics=dy,phase=ph))
    if start_gap:sc['start_gap']=start_gap
    return sc


def sensor_seed(cell,dev=False):
    c,t,dy,ph,i=cell
    if dev:return 900+i  # d_compare development sensor seeds, same scene seeds
    return SENSOR_BASE+1000*next(j for j,ct in enumerate(D.TEST) if ct==(c,t))+100*(dy=='viscous')+10*(ph==28.)+i


def key(r):return (r['cond'],r['tissue'],r['dynamics'],r['phase'],r['i'])


class CommandAudit:
    def __init__(self,pol):self.pol=pol;self.bad=0;self.clipped=0
    def reset(self,platform):self.pol.reset(platform);self.limit=platform['VMAX']*platform['Tc']
    def __call__(self,o):
        u=self.pol(o)
        self.bad+=int(not np.isfinite(u))
        self.clipped+=int(np.isfinite(u) and abs(u-o.log[o.n-1,11])>self.limit+1e-12)
        return u


def execute(a):
    method,cfg,cell,base,gap,keep,dev=a
    sc=scene(cell,base,gap)
    if method.startswith('Xue') and gap==0.:
        calibration=load('baseline_calibration.json')
        sc['initial_contact_observed']=calibration['positions'][calkey(cell,dev)]
    pol=build(method,cfg,sc)
    audit=CommandAudit(pol)
    start=time.perf_counter();log=P3.episode(sc,audit,sensor_seed(cell,dev));elapsed=time.perf_counter()-start
    r=M.metric(log,pol,sc)
    c,t,dy,ph,i=cell
    r.update(method=method,cfg=cfg,cond=c,tissue=t,dynamics=dy,phase=ph,i=i,seconds=elapsed,
             nonfinite_commands=audit.bad,saturated_commands=audit.clipped)
    if audit.bad:r['failed']=True;r['J']=1000.
    if keep and log is not None:
        folder=OUT/'baseline_logs';folder.mkdir(exist_ok=True)
        path=folder/f'{method}_{c}_{t}_{dy}_{int(ph)}_{i}.npz'
        np.savez_compressed(path,log=log[:,SAVE_COLS],cols=SAVE_COLS)
        r['log_file']=str(path.relative_to(OUT)).replace('\\','/')
    return r


class CaptureCommands:
    """Echo saved commands and capture the pre-tick observation, not truth.

    P3 writes the NEW command and post-substep x/Ft/m into the current log
    row after the callback. Historical rows and the newest row therefore
    need different treatment during exact causal fixed-source replay.
    """
    def __init__(self,commands):self.commands=commands;self.observed=[]
    def reset(self,platform):self.observed=[]
    def __call__(self,o):
        self.observed.append(o.log[o.n-1,OBS_COLS].copy())
        return self.commands[o.n-1]


def capture_job(source):
    path=OUT/source['log_file']
    with np.load(path) as z:original=z['log'].copy();cols=z['cols'].copy()
    sc=scene(key(source),TEST_BASE,0.)
    policy=CaptureCommands(original[:,list(cols).index(11)])
    log=P3.episode(sc,policy,sensor_seed(key(source)))
    if log is None or not np.array_equal(log[:,cols],original):
        raise AssertionError(('scripted physical replay differs',key(source),source['method']))
    np.savez_compressed(path,log=original,cols=cols,observed=np.asarray(policy.observed),observed_cols=OBS_COLS)
    return dict(source=source['method'],cell=key(source),physical_log_bitwise_equal=True,callbacks=len(policy.observed))


def capture(gap):
    assert gap==0.
    sources=[r for r in load('baseline_test.json')['rows'] if r['method'] in ('Huber-D4','D4') and 'log_file' in r]
    start=time.perf_counter();rows=run(sources,worker=capture_job)
    dump('baseline_capture.json',dict(rows=rows,seconds=time.perf_counter()-start))


def source_observation(log,observed,n):
    """Historical finalized log plus exactly the measured newest callback row."""
    obs=P3.Obs(log,n);obs.log[-1,OBS_COLS]=observed[n-1]
    return obs


def run(jobs,worker=execute,workers=15):
    start=time.perf_counter();rows=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for j,r in enumerate(ex.map(worker,jobs,chunksize=1),1):
            rows.append(r)
            if j%48==0 or j==len(jobs):
                print(f'{j}/{len(jobs)}, elapsed {time.perf_counter()-start:.1f}s',flush=True)
    return rows


def dump(name,obj):
    # RFC8259 JSON: non-finite measurements are explicit nulls.
    def clean(v):
        if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
        if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
        if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
        if isinstance(v,np.integer):return int(v)
        if isinstance(v,np.bool_):return bool(v)
        return v
    (OUT/name).write_text(json.dumps(clean(obj),ensure_ascii=False,indent=1,allow_nan=False),encoding='utf-8')


def load(name):return json.loads((OUT/name).read_text(encoding='utf-8'))
def L(rows):return float(np.mean([math.log1p(r['J']) for r in rows]))


def calkey(cell,dev):return ('dev:' if dev else 'test:')+':'.join(map(str,cell))


class ContactCalibration(K.common.DirectPolicy):
    """Known first-contact position required by Xue Assumption 1, sensor only."""
    def __init__(self):self.params={}
    def reset(self,platform):super().reset(platform);self.position=None;self.detected_at=None
    def step(self,o):
        row=o.log[o.n-1]
        if self.position is None and row[5]>.5:
            self.position=float(row[2]);self.detected_at=float(row[0])
        return float(row[11]+(.005*self.Tc if self.position is None else 0.))


def calibration_job(a):
    cell,dev=a;sc=scene(cell,D.BASE if dev else TEST_BASE,0.)
    # A copy of the same object, stationary at the initial scan position.
    # No force/depth/geometry truth is passed to the acquisition controller.
    sc['T']=2.;sc['start_gap']=.001
    for k in ('u','act','fext','rip'):sc[k]=np.full_like(sc[k],sc[k][0])
    sc['xcmd']=np.full_like(sc['xcmd'],sc['xcmd'][0])
    pol=M.attach(ContactCalibration(),sc);tic=time.perf_counter()
    log=P3.episode(sc,pol,sensor_seed(cell,dev)+50000)
    if log is None or pol.position is None:raise RuntimeError(f'first contact not acquired: {cell}')
    return dict(key=calkey(cell,dev),position=pol.position,detected_at=pol.detected_at,seconds=time.perf_counter()-tic)


def calibrate():
    dev=[(c,t,'neural_white',21.,i) for c,t in D.DEV for i in range(2)]
    rows=run([(c,True) for c in dev]+[(c,False) for c in cells()],worker=calibration_job)
    dump('baseline_calibration.json',dict(positions={r['key']:r['position'] for r in rows},rows=rows,
                                        threshold_N=.5,approach_speed=.005,extra_simulated_seconds=2.))


def develop(gap,methods=METHODS):
    rows=[];chosen={};scores={};start=time.perf_counter()
    for method in methods:
        cfgs=configurations(method);assert len(cfgs)==24
        dev=[(c,t,'neural_white',21.,i) for c,t in D.DEV for i in range(2)]
        rr=run([(method,cfg,cell,D.BASE,gap,False,True) for cfg in cfgs for cell in dev])
        val=[L(rr[j*8:(j+1)*8]) for j in range(24)];idx=int(np.argmin(val))
        chosen[method]=cfgs[idx];scores[method]=dict(default=val[0],best=val[idx],index=idx,all=val,configs=cfgs)
        rows+=rr
        print(method,scores[method]['default'],scores[method]['best'],chosen[method],flush=True)
        dump('baseline_dev.json',dict(start_gap=gap,base=D.BASE,chosen=chosen,scores=scores,rows=rows,seconds=time.perf_counter()-start))


def test(gap):
    dev=load('baseline_dev.json');assert gap==dev['start_gap']
    rows=[];start=time.perf_counter()
    configs=dict(dev['chosen'])
    if configs['Huber-D4']!=configurations('Huber-D4')[0]:
        configs['Huber-current']=configurations('Huber-D4')[0]
    for method,cfg in configs.items():
        rr=run([(method,cfg,cell,TEST_BASE,gap,method in ('Huber-D4','D4'),False) for cell in cells()])
        rows+=rr;print(method,L(rr),flush=True)
        dump('baseline_test.json',dict(start_gap=gap,base=TEST_BASE,sensor_base=SENSOR_BASE,chosen=configs,rows=rows,seconds=time.perf_counter()-start))


def recovery(t,k,kh,valid):
    """Natural sharp-transition clusters; recovery within 20% for 5 ticks.

    No ideal imposed steps exist in this scan grid. Return censored records,
    not a fabricated stiffness-step response. Also report synthetic steps in QA.
    """
    changing=(np.abs(np.gradient(np.log(np.maximum(k,1.)),t))>2.)&(t>1.)&(k>200)
    inds=np.flatnonzero(changing);clusters=[]
    if not len(inds):return []
    cuts=np.r_[0,np.flatnonzero(np.diff(inds)>25)+1,len(inds)]
    for a,b in zip(cuts[:-1],cuts[1:]):
        start,end=inds[a],inds[b-1]
        if abs(math.log(k[end]/k[start]))<math.log(1.2):continue
        limit=min(len(t)-5,end+int(1.5/.008));delay=None
        for j in range(end,limit):
            if np.all(valid[j:j+5] & (np.abs(np.log(np.maximum(kh[j:j+5],1.)/np.maximum(k[j:j+5],1.)))<=math.log(1.2))):
                delay=float(t[j]-t[end]);break
        clusters.append(dict(start=float(t[start]),end=float(t[end]),ratio=float(k[end]/k[start]),delay=delay,censor_seconds=float(t[limit]-t[end])))
    return clusters


def replay_job(a):
    source,method,cfg,gap=a;cell=key(source);sc=scene(cell,TEST_BASE,gap)
    if method.startswith('Xue') and gap==0.:
        sc['initial_contact_observed']=load('baseline_calibration.json')['positions'][calkey(cell,False)]
    saved=np.load(OUT/source['log_file']);log=np.full((len(saved['log']),26),np.nan);log[:,saved['cols']]=saved['log']
    observed=saved['observed']
    pol=build(method,cfg,sc);platform=dict(P3.PLATFORM,tau_a=sc['tau_a_hat']);pol.reset(platform)
    # Instrument every attempted estimator call, including invalid ones. Host
    # accept_trace omits invalid windows; using it as denominator inflates uptime.
    raw=[];orig=pol.estimator.estimate
    def record(o,n,p):
        a=time.perf_counter();k,sig,chi=orig(o,n,p);elapsed=time.perf_counter()-a
        raw.append((n-1,k,sig,chi,elapsed));return k,sig,chi
    pol.estimator.estimate=record
    command_error=0.
    for n in range(1,len(log)+1):
        requested=pol(source_observation(log,observed,n))
        if method==source['method']:
            sent=observed[n-1,OBS_COLS.index(11)];limit=P3.PLATFORM['VMAX']*P3.PLATFORM['Tc']
            command_error=max(command_error,abs(float(np.clip(requested,sent-limit,sent+limit)-log[n-1,11])))
    if method==source['method'] and command_error>=1e-10:
        raise AssertionError(('causal source command replay differs',method,cell,command_error))
    accepted={round(x[0],7):bool(x[4]) for x in pol.accept_trace}
    t=log[:,0];k=log[:,21];ev=np.zeros(len(t),bool)
    for e in sc['event_info']:ev|=(t>=e['start'])&(t<=e['end']+.3)
    changing=np.convolve(np.abs(np.gradient(np.log(np.maximum(k,1.)),t))>2.,np.ones(30),'same')>0
    # quiet means no detected event/change, not a physically motionless surface.
    groups={'quiet':[],'event':[],'changing':[]};kh=np.full(len(log),np.nan);okseq=np.zeros(len(log),bool)
    indexed={r[0]:r[1:] for r in raw}
    for j in range(len(log)):
        hat,sig,chi,elapsed=indexed.get(j,(np.nan,np.inf,np.inf,np.nan))
        ok=accepted.get(round(float(t[j]),7),False)
        kh[j]=hat;okseq[j]=ok
        if t[j]<1. or k[j]<=200.:continue
        kind='event' if ev[j] else ('changing' if changing[j] else 'quiet')
        finite=np.isfinite(hat) and hat>0 and np.isfinite(sig)
        err=math.log(hat/k[j]) if finite else np.nan
        groups[kind].append((err,ok,2*sig/hat if finite else np.nan,elapsed))
    result={}
    for kind,w in groups.items():
        arr=np.array(w,float).reshape(-1,4);finite=np.isfinite(arr[:,0]);ok=arr[:,1].astype(bool)&finite
        errs=np.abs(arr[finite,0]);ea=np.abs(arr[ok,0])
        result[kind]=dict(n=len(w),valid=int(finite.sum()),accepted=int(ok.sum()),
                          error_sum=float(errs.sum()),error_samples=errs.tolist(),accepted_errors=ea.tolist(),
                          covered=int(np.sum(ea<=arr[ok,2])),times=arr[np.isfinite(arr[:,3]),3].tolist())
    extra={}
    inner=getattr(pol.estimator,'inner',None)
    if isinstance(inner,I.XueHC):
        extra=dict(invalid_depth_or_contact=inner.invalid_depth,updates=inner.updates,
                   final_theta=inner.theta.tolist(),initial_contact_position=inner.contact)
    record=dict(source=source['method'],method=method,cond=cell[0],tissue=cell[1],dynamics=cell[2],phase=cell[3],i=cell[4],
                groups=result,recovery=recovery(t,k,kh,okseq),diagnostic=extra,
                source_command_max_error=command_error if method==source['method'] else None)
    if source.get('retain_series', False):
        record['series']=dict(t=t,truth=k,estimate=kh,accepted=okseq,
                              sigma=np.array([indexed.get(j,(np.nan,np.inf,np.inf,np.nan))[1] for j in range(len(log))]))
    return record


def replay(gap):
    test=load('baseline_test.json');sources=[r for r in test['rows'] if r['method'] in ('Huber-D4','D4') and 'log_file' in r]
    start=time.perf_counter();rows=[]
    for method in ID_METHODS:
        rr=run([(s,method,test['chosen'][method],gap) for s in sources],worker=replay_job)
        rows+=rr
        # Compact per-scenario summaries; pooled samples remain in a compressed
        # NPZ so report aggregation doesn't require huge JSON arrays.
        folder=OUT/'baseline_replay';folder.mkdir(exist_ok=True)
        for r in rr:
            data={}
            for kind,g in r['groups'].items():
                for field in ('error_samples','accepted_errors','times'):
                    data[kind+'_'+field]=np.array(g.pop(field))
            path=folder/f"{method}_{r['source']}_{r['cond']}_{r['dynamics']}_{int(r['phase'])}_{r['i']}.npz"
            np.savez_compressed(path,**data);r['samples_file']=str(path.relative_to(OUT)).replace('\\','/')
        dump('baseline_module.json',dict(start_gap=gap,rows=rows,seconds=time.perf_counter()-start))


def pilot(gap):
    start=time.perf_counter();jobs=[]
    for method in METHODS:
        jobs.append((method,configurations(method)[0],('breath','relaxed','neural_white',21.,0),D.BASE,gap,False,True))
    rows=run(jobs,workers=15);dump('baseline_pilot.json',dict(rows=rows,seconds=time.perf_counter()-start))
    for r in rows:print(r['method'],round(r['seconds'],2),'s',round(math.log1p(r['J']),3),r['failed'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['calibrate','pilot','dev','test','capture','replay','report'],default='pilot');parser.add_argument('--start-gap',type=float,default=0.)
    args=parser.parse_args();OUT.mkdir(exist_ok=True)
    if args.stage=='report':
        from baseline_report import generate
        generate()
    elif args.stage=='calibrate':calibrate()
    else:globals()[args.stage if args.stage!='dev' else 'develop'](args.start_gap)
