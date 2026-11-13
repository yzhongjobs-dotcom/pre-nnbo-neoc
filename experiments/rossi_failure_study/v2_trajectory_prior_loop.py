"""Prespecified V2 trajectory prior: 66 closed loops, 8 workers, read-only references."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.environ['TRAJECTORY_LOOP_PROJECT']=r'E:\github_project\scan-impedance-update'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
PROJECT=Path(os.environ['TRAJECTORY_LOOP_PROJECT'])
sys.path.insert(0,str(PROJECT/'experiments/rossi_failure_study'))
import trajectory_loop as TL
import trajectory_robust as TR
import huber_windows as HW
from d_robust import common
from scipy.linalg import solveh_banded
import numpy as np
import argparse, hashlib, json, math, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor,as_completed
BASE=PROJECT/'experiments/rossi_failure_study/results'
OUT=Path(__file__).resolve().parent/'results/v2_trajectory_prior'
ARMS=('huber','delag2','v2t')

class PriorEstimate:
    def __init__(self,est,zero=False):
        self.est,self.original,self.zero=est,est.estimate,zero
        self.filter=TL.DelagCorrection(force_ratio_factor=False)
        self.arrays={k:[] for k in ('t_est','k_w','k_ret','g','vt','vp','sv','tau_w','beta')}
        self.checks=dict(finite_calls=0,m0_checks=0,m0_violations=0,m0_max_relative=0.,clipped_calls=0,fallback=0,beta_zero_ticks=0,beta_zero_mismatches=0)
    def statistics(self,o,n,p,kw):
        est=self.est; log=est._log(o,n)
        j0=max(p.ndel+1,n-1-est.nw); m=n-2-j0+1
        te,ye=log[n-1,0],log[n-1,2]
        _,Y,z,keep=common._build(log,j0,m,0,p.ns,p.eb,p.afc,p.ndel,te,ye,est.nw*p.ns*1e-3,p.params['enc_y'])
        h1=HW.time_weighted(log,j0,m,p.ns,p.eb,p.afc,te,ye,p.params['enc_y'])
        broken=np.flatnonzero(~keep)
        if len(broken):
            cut=broken[-1]+1; Y,z,h1=Y[cut:],z[cut:],h1[cut:]
        m=len(Y)
        assert est.order==4 and len(est.weights)==m-4
        # Exact same band representation of I+rho D' diag(w) D.
        _,basis=est._difference_bands(m)
        ab=(basis@est.weights).reshape((5,m)); ab[0]+=1.
        rhs=np.column_stack((Y,z,h1))
        my,mz,mh=(rhs-solveh_banded(ab,rhs,lower=True,check_finite=False)).T
        a,b,c1=float(Y@my),float(Y@mh),float(h1@mh)
        p0,p1=float(my@z),float(h1@mz)
        m0=p0/a
        if 200.<kw<3e5:
            rel=abs(m0-kw)/abs(kw)
            self.checks['m0_checks']+=1
            self.checks['m0_max_relative']=max(self.checks['m0_max_relative'],rel)
            self.checks['m0_violations']+=int(not np.isfinite(rel) or rel>1e-8)
        else: self.checks['clipped_calls']+=1
        Jv,tau=c1-b*b/a,-b/a
        if Jv<=0 or tau<=0 or not np.all(np.isfinite((a,b,c1,p0,p1,m0,Jv,tau))): return None
        vt=(p1-b*m0)/Jv; sv=est.sigma_scale*p.params['noise']/math.sqrt(Jv)
        return (vt,sv,tau) if np.all(np.isfinite((vt,sv,tau))) else None
    def __call__(self,o,n,p):
        kw,sigma,chi=self.original(o,n,p); stats=None
        if np.isfinite(kw):
            self.checks['finite_calls']+=1; stats=self.statistics(o,n,p,kw)
        if self.zero: self.filter.beta=0.
        kd,beta=self.filter.step(n,float(o.log[n-1,2]),float(o.log[n-1,5]),kw,sigma)
        if self.zero: self.filter.beta=0.
        k,sig=kw,sigma; g,vt,vp,sv,tau=(np.nan,)*5
        if np.isfinite(kw):
            if stats is None:
                self.checks['fallback']+=1; k=kd; sig=sigma if kd==kw else sigma*kd/kw
            else:
                vt,sv,tau=stats; f=self.filter; j=len(f.forces)-1
                F24=(f.f_sums[j+1]-f.f_sums[max(0,j-23)])/(j+1-max(0,j-23))
                if j>=23 and F24>=1. and np.isfinite(sigma):
                    y6=(f.y_sums[j+1]-f.y_sums[max(0,j-5)])/(j+1-max(0,j-5))
                    r=float(np.clip(1.+beta*(y6-f.y24[j])/f.cw[j],.5,2.))
                    vp=kw*(1./r-1.)/tau; Sp=.5*abs(vp)
                    g=Sp**2/(Sp**2+sv**2) if Sp else 0.
                    vh=vp+g*(vt-vp)
                    if np.all(np.isfinite((vp,Sp,g,vh))):
                        if beta!=0.:
                            k=float(np.clip(np.clip(kw+tau*vh,.5*kw,2.*kw),200.,3e5))
                            sig=math.sqrt(sigma**2+tau**2*(1.-g)*Sp**2)
                    else:
                        self.checks['fallback']+=1; k=kd; sig=sigma if kd==kw else sigma*kd/kw
        if beta==0.:
            self.checks['beta_zero_ticks']+=1
            self.checks['beta_zero_mismatches']+=int(np.float64(k).view(np.uint64)!=np.float64(kw).view(np.uint64))
        for key,val in zip(self.arrays,(float(o.log[n-1,0]),kw,k,g,vt,vp,sv,tau,beta)): self.arrays[key].append(val)
        return k,sig,chi

def event_stats(r,z):
    te,t,kr=z['t_est'],z['t'],z['k_ret']; ev=r['event_info'][0]
    truth=np.interp(te,t,z['k_eff']); fax=np.interp(te,t,z['fax'])
    mask=(te>=ev['start'])&(te<=ev['end']+.6)
    good=mask&(te>1.5)&(truth>0)&(fax>.5)&np.isfinite(kr)&(kr>0)
    e=np.log(kr[good]/truth[good]); fm=(t>=ev['start'])&(t<=ev['end']+.6)
    return dict(force_rms=TL.D.rms(z['fax'][fm]-r['FD']),mean_log_error=float(e.mean()),min_log_error=float(e.min()),beta_before=float(z['beta'][te<ev['start']][-1]),mean_g=float(np.nanmean(z['g'][mask])) if 'g' in z else None)

def job(task):
    TL.O.setup(); sc=TL.SE.MC.scene_for(tuple(task['cell']),TL.CC.BASE,task['variant'])
    if task['block']==5: TR.shift_event(sc,task['sign'])
    pol=TL.B.build('Huber-D4',TL.O.chosen()['Huber-D4'],sc); TL.CC.own_mechanism(pol)
    rec=PriorEstimate(pol.estimator,task.get('fixed_zero',False)); pol.estimator.estimate=rec
    log=TL.P3.episode(sc,pol,task['seed']); metric=TL.M.metric(log,pol,sc)
    row=dict(task,FD=float(sc['FD']),speed_mm_s=float(sc['v']*1000),event_info=sc['event_info'],L=float(np.log1p(metric['J'])),**{k:float(metric[k]) for k in ('rmse','peak','over_15','loss_time')},checks=rec.checks)
    arrays={k:np.asarray(v) for k,v in rec.arrays.items()}; arrays.update(t=log[:,0],fax=log[:,18],k_eff=log[:,26])
    filename=f"run_{task['run_id']:03d}_v2t.npz"
    with (OUT/filename).open('xb') as stream: np.savez_compressed(stream,**arrays)
    row['file']=filename
    if task['block']==5: row['event']=event_stats(row,arrays)
    return row

def worker(task,spool):
    stem=str(task['run_id']); spec,result,stdout=(spool/(stem+s) for s in ('.task.json','.json','.stdout'))
    spec.write_text(json.dumps(task),encoding='utf-8')
    cmd=[sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--worker-task',str(spec),'--worker-result',str(result)]
    with stdout.open('x',encoding='utf-8') as stream:
        done=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=stream,stderr=stream,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
    if done.returncode: raise RuntimeError(stdout.read_text(encoding='utf-8'))
    return json.loads(result.read_text(encoding='utf-8'))

REFERENCE_MIRROR=Path(r'D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/results')
def folder(r):
    if r['arm']=='v2t': return OUT
    name='trajectory_loop' if r['block']==0 else 'trajectory_robust'
    if (BASE/name/r['file']).exists(): return BASE/name
    for source in ('trajectory_loop.json','trajectory_robust.json'):
        assert (BASE/source).read_bytes()==(REFERENCE_MIRROR/source).read_bytes(), 'Reference mirror differs from repository'
    return REFERENCE_MIRROR/name
def load(r):
    with np.load(folder(r)/r['file'],allow_pickle=False) as z: return {k:z[k] for k in z.files}

def summarize(rows):
    groups={g:{a:dict(mean=float(np.mean([r['L'] for r in rows if r['group']==g and r['arm']==a])),median=float(np.median([r['L'] for r in rows if r['group']==g and r['arm']==a]))) for a in ARMS} for g in ('A','A_slow_thin','B','A+B','A_large original','A_large fresh','fresh A')}
    ident={}
    for label,rr in (('original_fast',[r for r in rows if r['block']==0 and r['index'] in TL.FAST]),('fresh_large_all',[r for r in rows if r['block']==1]),('fresh_large_fast',[r for r in rows if r['block']==1 and r['speed_mm_s']>15.])):
        ident[label]={}
        for arm in ('delag2','v2t'):
            scored=[TL.scored_arrays(r,folder(r)) for r in rr if r['arm']==arm]
            e=np.concatenate([d[0] for d in scored]); vals=[d[2] for d in scored if d[2] is not None]
            ident[label][arm]=dict(segments={s:TL.error_stats(e,np.concatenate([d[1][s] for d in scored])) for s in ('plateau','rise','fall','reversal')},sustained=dict(median=float(np.median(vals)),worst=float(max(vals))))
    events=[]
    for v in [r for r in rows if r['block']==5 and r['arm']=='v2t']:
        pair={a:next(r for r in rows if r['block']==5 and r['scene']==v['scene'] and r['sign']==v['sign'] and r['arm']==a) for a in ARMS}
        entry=dict(scene=v['scene'],sign=v['sign'],direction=pair['huber']['direction'],arms={})
        for a,r in pair.items(): entry['arms'][a]=dict(L=r['L'],peak=r['peak'],**(r['event'] if a=='huber' else event_stats(r,load(r))))
        events.append(entry)
    change=[r for r in rows if r['arm']=='v2t' and ((r['block']==0 and r['index'] in (36,37)) or (r['block']==1 and r['scene'] in ('neural_white_28_i5','viscous_14_i5')))]
    event=[r for r in rows if r['arm']=='v2t' and r['block']==5 and r['event']['beta_before']>0.]
    mechanism={}
    for label,rr in (('traced_large_change',change),('learned_event',event)):
        gg,tt,ss=[],[],[]
        for r in rr:
            z=load(r); te=z['t_est']
            if label=='traced_large_change':
                _,m,_,_=TL.scored_arrays(r,folder(r)); mask=m['rise']|m['fall']|m['reversal']
            else:
                ev=r['event_info'][0]; mask=(te>=ev['start'])&(te<=ev['end']+.6)
            gg.extend(z['g'][mask&np.isfinite(z['g'])]); tt.extend(z['tau_w'][mask&np.isfinite(z['tau_w'])])
            ratio=z['sv']/z['k_w']; ss.extend(ratio[mask&np.isfinite(ratio)])
        mechanism[label]=dict(runs=len(rr),g_ticks=len(gg),mean_g=float(np.mean(gg)),median_tau_w=float(np.median(tt)),median_sv_over_kw=float(np.median(ss)))
    cases=[dict(block=r['block'],scene=r['scene'],index=r['index'],L={a:next(q['L'] for q in rows if q['block']==r['block'] and q['scene']==r['scene'] and q['arm']==a) for a in ARMS}) for r in rows if r['arm']=='v2t' and r['group'] in ('A_large original','A_large fresh')]
    gd={g:groups[g]['v2t']['mean']-groups[g]['delag2']['mean'] for g in ('A_large original','A_large fresh')}
    di={label:{s:dict(mean_delta=ident[label]['v2t']['segments'][s]['mean']-ident[label]['delag2']['segments'][s]['mean'],median_abs_delta=ident[label]['v2t']['segments'][s]['median_abs']-ident[label]['delag2']['segments'][s]['median_abs']) for s in ('rise','fall')} for label in ('original_fast','fresh_large_fast')}
    learned=[e for e in events if e['arms']['delag2']['beta_before']>0.]
    ad={e['scene']:e['arms']['v2t']['L']-e['arms']['delag2']['L'] for e in learned if e['direction']=='away from probe'}
    n6=[e['arms']['v2t']['L'] for e in learned if e['scene']=='neural_white_28_i6']
    rex={e['scene']+':'+str(e['sign']):e['arms']['v2t']['force_rms']-(e['arms']['huber']['force_rms']+e['arms']['delag2']['force_rms'])/2 for e in learned}
    zd=[r['L']-next(q['L'] for q in rows if q['block']==r['block'] and q['scene']==r['scene'] and q['arm']=='huber') for r in rows if r['arm']=='v2t' and ((r['block']==0 and r['index']<24) or r['block']==2)]
    predictions=dict(V1=dict(met=bool(abs(gd['A_large original'])<=.02 and abs(gd['A_large fresh'])<=.03 and all(abs(d['mean_delta'])<=.03 and d['median_abs_delta']<=.02 for x in di.values() for d in x.values())),group_deltas=gd,identification_deltas=di),V2=dict(met=bool(max(ad.values())<=0 and max(n6)<=.40 and max(rex.values())<=0),away_L_deltas=ad,N6_max_L=max(n6),event_rms_excess=rex),V3=dict(met=bool(max(map(abs,zd))<=.003),max_abs_L_delta=max(map(abs,zd))),V4=dict(met=bool(mechanism['traced_large_change']['mean_g']<=.3 and mechanism['learned_event']['mean_g']>=.5)))
    per_case={}
    for r in rows:
        if r['arm'] in ('delag2','v2t') and ((r['block']==0 and r['index'] in TL.FAST) or r['block']==1):
            e,mm,su,_=TL.scored_arrays(r,folder(r))
            per_case.setdefault(r['scene'],{})[r['arm']]=dict(segments={s:TL.error_stats(e,mm[s]) for s in ('plateau','rise','fall','reversal')},sustained=su)
    return dict(groups=groups,large_cases=cases,identification=ident,per_case_identification=per_case,late_events=events,mechanism=mechanism,predictions=predictions)

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--worker-task',type=Path); parser.add_argument('--worker-result',type=Path); parser.add_argument('--summarize',action='store_true'); args=parser.parse_args()
    if args.worker_task:
        row=job(json.loads(args.worker_task.read_text(encoding='utf-8')))
        with args.worker_result.open('x',encoding='utf-8') as stream: json.dump(row,stream,allow_nan=False)
        return
    dest=OUT.parent/'v2_trajectory_prior_loop.json'
    if args.summarize:
        obj=json.loads(dest.read_text(encoding='utf-8'))
        obj['summary']=summarize(obj['rows'])
        obj['conventions']['fresh_labels']='N/V labels use actual scene i: fresh i1,i3,i4,i5,i7,i8; late i0,i2,i6'
        obj['reference_trace_mirror']=str(REFERENCE_MIRROR)
        obj['reference_mirror_json_byte_identity']=all((BASE/n).read_bytes()==(REFERENCE_MIRROR/n).read_bytes() for n in ('trajectory_loop.json','trajectory_robust.json'))
        candidate=[r for r in obj['rows'] if r['arm']=='v2t']
        obj['checks']={k:sum(r['checks'][k] for r in candidate) for k in candidate[0]['checks'] if k!='m0_max_relative'}
        obj['checks']['m0_max_relative']=max(r['checks']['m0_max_relative'] for r in candidate)
        required=('t','fax','k_eff','t_est','k_w','k_ret','g','vt','vp','sv','tau_w','beta')
        for r in candidate:
            z=load(r)
            assert all(k in z for k in required)
            assert all(len(z[k])==len(z['t_est']) for k in required[3:])
        assert len(obj['rows'])==198 and len(candidate)==66
        assert len({(r['block'],r['scene'],r['sign'],r['arm']) for r in obj['rows']})==198
        fixed=np.load(OUT/'run_066_v2t.npz')
        ref=np.load(BASE/'trajectory_loop/case_36_huber.npz')
        obj['beta_zero_validation']['force_trace_bitwise_equal']=bool(np.array_equal(fixed['fax'].view(np.uint64),ref['fax'].view(np.uint64)))
        obj['beta_zero_validation']['truth_trace_bitwise_equal']=bool(np.array_equal(fixed['k_eff'].view(np.uint64),ref['k_eff'].view(np.uint64)))
        fixed.close(); ref.close()
        for name,digest in obj['sources'].items(): assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
        obj['validation']='198 unique rows; all66 candidate trace schemas/lengths and source hashes verified; fixed-zero force/truth bitwise checked'
        dest.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(dict(checks=obj['checks'],beta_zero=obj['beta_zero_validation'],validation=obj['validation'])))
        return
    if dest.exists() or OUT.exists(): raise FileExistsError('Refusing to overwrite outputs')
    OUT.mkdir(parents=True)
    loop=json.loads((BASE/'trajectory_loop.json').read_text(encoding='utf-8')); robust=json.loads((BASE/'trajectory_robust.json').read_text(encoding='utf-8')); rows,tasks=[],[]
    for r in loop['rows']:
        if r['arm'] in ('huber','delag2'):
            case=TL.CC.cases()[r['index']]; group='A_large original' if r['group']=='A_large' else r['group']
            rows.append(dict(r,block=0,scene=f"case_{r['index']}",sign=1,group=group))
            if r['arm']=='huber': tasks.append(dict(case,index=r['index'],block=0,scene=f"case_{r['index']}",sign=1,arm='v2t',group=group))
    for r in robust['rows']:
        if r['arm'] in ('huber','delag2') and r['block'] in (1,2,5):
            group={1:'A_large fresh',2:'fresh A',5:'late push/pull'}[r['block']]; rows.append(dict(r,group=group))
            if r['arm']=='huber': tasks.append(dict({k:r[k] for k in ('block','scene','cell','variant','seed','index','sign')},arm='v2t',group=group))
    assert len(tasks)==66 and len(rows)==132
    for i,task in enumerate(tasks): task['run_id']=i
    check=dict(tasks[24],run_id=66,fixed_zero=True); assert check['index']==36
    spool=Path(tempfile.mkdtemp(prefix='v2_prior_workers_')); print(f'66 runs + zero check, 8 workers; IPC {spool}',flush=True); start=time.perf_counter(); zero=None
    with ThreadPoolExecutor(8) as pool:
        pending=[pool.submit(worker,t,spool) for t in tasks+[check]]
        for number,future in enumerate(as_completed(pending),1):
            row=future.result()
            if row.get('fixed_zero'): zero=row
            else: rows.append(row)
            print(f"{number}/67 {row['scene']} sign={row['sign']} L={row['L']:.9f} checks={row['checks']}",flush=True)
    ref=next(r for r in rows if r['block']==0 and r['index']==36 and r['arm']=='huber')
    zero_check=dict(checks=zero['checks'],beta_all_zero=bool(np.all(load(zero)['beta']==0.)),metric_deltas={k:zero[k]-ref[k] for k in ('L','rmse','peak','over_15','loss_time')})
    rows.sort(key=lambda r:(r['block'],r['scene'],r['sign'],r['arm']))
    obj=dict(status='complete',workers=8,new_episodes=66,validation_episodes=1,rows=rows,beta_zero_validation=zero_check,runtime_s=time.perf_counter()-start,sources={n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in ('trajectory_loop.json','trajectory_robust.json')},conventions=dict(error='ln(k_ret/k_eff), original masks/validity and sustained 13 contiguous ticks',beta='Imported delag2 history and gain learner; pre-update beta',mechanism='Pooled finite ticks, cases36/37+fresh N5/V5 change masks; learned events with pre-event beta>0',fresh_labels='N/V labels use actual scene i: fresh i1,i3,i4,i5,i7,i8; late i0,i2,i6',matrix='Exact banded M action, final frozen Huber weights',traces='All66 + fixed-zero validation saved; requested subset included'))
    # Preserve primary completed rows before read-only analysis.
    with dest.open('x',encoding='utf-8') as stream: json.dump(obj,stream,indent=2,allow_nan=False)
    obj['summary']=summarize(rows); dest.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('SUMMARY '+json.dumps(obj['summary']),flush=True); print('ZERO '+json.dumps(zero_check),flush=True); print(f'Written {dest}',flush=True)
if __name__=='__main__': main()




