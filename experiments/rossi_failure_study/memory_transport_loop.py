"""Fixed log-memory closed loop. Eight subprocess workers; source repository read only."""
import os, sys
sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
from pathlib import Path
PROJECT = Path(sys.argv[sys.argv.index('--project-root')+1])
os.environ['TRAJECTORY_LOOP_PROJECT'] = str(PROJECT)
sys.path.insert(0,str(PROJECT/'experiments/rossi_failure_study'))
import trajectory_loop as TL
import probe_band_loop as PB
import argparse, hashlib, json, math, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
ARMS = ('huber','delag2','kf','tmem24','tmem48','mem24')
SEG = ('plateau','rise','fall','reversal')
PREDICTIONS = dict(
 K1='kf original A_large mean L<=.195; fresh A_large mean L<=.265.',
 K2='kf case7 L<.5.',
 K3='kf-huber L within +/-.02 case by case, originals 0..23 except7.',
 K4='Combined original/fresh A_large mean L(kf_noT)-L(kf)>=.03.',
 K5='Every run of every arm: total time k_ret/k_eff>3 after .5s <=.5s.',
 K6='band_kf case7 L<=kf; elsewhere group means within +/-.02 of kf.',
 K7='Pooled 30-original plateau median |e| reduction >=15% kf, >=35% kf_long; kf sustained-underestimate median below huber in every group.')

class LogMemory:
 def __init__(self,arm):
  self.q = 1e-4 if arm=='kf_long' else 5e-4
  self.noT = arm=='kf_noT'
  self.resets=0
  self.reset()
 def reset(self):
  self.gain=TL.DelagCorrection(force_ratio_factor=False)
  self.n=None
  self.lam=self.P=np.nan
  self.active=False
  self.missing=0
  self.y6_prev=None
 def step(self,n,y,force,kw,sigma,chi):
  if self.n is not None and n!=self.n+1:
   self.resets+=1
   self.reset()
  self.n=n
  g=self.gain
  j=len(g.forces)
  ys,fs=g.y_sums[-1]+y,g.f_sums[-1]+force
  lo6,lo24=max(0,j-5),max(0,j-23)
  y6=(ys-g.y_sums[lo6])/(j+1-lo6)
  y24=(ys-g.y_sums[lo24])/(j+1-lo24)
  F24=(fs-g.f_sums[lo24])/(j+1-lo24)
  beta=0. if self.noT else g.beta
  finite=np.isfinite(kw)
  rel=sigma/kw if finite else np.inf
  self.missing=0 if finite else self.missing+1
  if self.missing>25: self.active=False
  was_active=self.active
  if was_active:
   u=0.
   if self.y6_prev is not None and F24>=1.:
    u=-math.log(float(np.clip(1.+beta*(y6-self.y6_prev)/(F24/math.exp(self.lam)),.8,1.25)))
    self.lam+=u
    self.P+=self.q+(.5*u)**2
   if finite and np.isfinite(sigma) and rel<=3.:
    Cw=F24/kw
    d=math.log(float(np.clip(1.+beta*(y6-y24)/Cw,.5,2.))) if Cw!=0. else 0.
    z=math.log(kw)-d
    R=24.*rel**2+(.5*d)**2
    K=self.P/(self.P+R)
    self.lam+=K*(z-self.lam)
    self.P*=1.-K
  if not was_active and finite and np.isfinite(sigma) and rel<=1.:
   self.lam,self.P=math.log(kw),(2.*rel)**2
   self.active=True
  # Imported, unchanged delag2 beta update, after using the old beta twice.
  g.step(n,y,force,kw,sigma)
  self.y6_prev=y6
  if not self.active: return (kw,sigma,chi),beta
  k=float(np.clip(math.exp(self.lam),200.,3e5))
  return (k,k*math.sqrt(self.P),chi),beta

class FiniteMemory:
 def __init__(self,arm):
  self.H={'tmem24':24,'tmem48':48,'mem24':24,'h0':0}[arm]
  self.noT=arm=='mem24'; self.resets=0; self.reset()
 def reset(self):
  self.gain=TL.DelagCorrection(force_ratio_factor=False); self.n=None
  self.windows=[]; self.P=np.nan; self.active=False
 def step(self,n,y,force,kw,sigma,chi):
  if self.n is not None and n!=self.n+1:
   self.resets+=1; self.reset()
  self.n=n; g=self.gain; j=len(g.forces)
  ys,fs=g.y_sums[-1]+y,g.f_sums[-1]+force
  lo6,lo24=max(0,j-5),max(0,j-23)
  y6=(ys-g.y_sums[lo6])/(j+1-lo6); y24=(ys-g.y_sums[lo24])/(j+1-lo24)
  F24=(fs-g.f_sums[lo24])/(j+1-lo24); beta=0. if self.noT else g.beta
  finite=np.isfinite(kw) and kw>0
  self.windows.append((math.log(kw) if finite else np.nan,F24/kw if finite else np.nan,
                       y24,sigma/kw if finite else np.nan))
  terms=[]
  for i in range(max(0,j-self.H),j+1):
   lk,cw,yi,rel=self.windows[i]
   if not (np.isfinite(lk) and np.isfinite(rel) and rel<=3. and np.isfinite(cw) and cw!=0.): continue
   d=math.log(float(np.clip(1.+beta*(y6-yi)/cw,.5,2.)))
   v=rel**2+(.5*d)**2
   if v<=0.: continue
   terms.append((i,1./v,lk-d))
  self.active=bool(terms); self.P=np.nan; triple=(kw,sigma,chi)
  if terms:
   sw=sum(x[1] for x in terms); lam=sum(x[1]*x[2] for x in terms)/sw
   self.P=(1./sw)*len(terms)/((j-terms[0][0]+24)/24)
   k=float(np.clip(math.exp(lam),200.,3e5)); triple=(k,k*math.sqrt(self.P),chi)
  g.step(n,y,force,kw,sigma)
  return triple,beta

class Recorder:
 def __init__(self,original,arm):
  self.original,self.memory=original,FiniteMemory(arm)
  self.t,self.kw,self.kr,self.p,self.beta=[],[],[],[],[]
 def __call__(self,o,n,policy):
  kw,sigma,chi=self.original(o,n,policy)
  triple,beta=self.memory.step(n,float(o.log[n-1,2]),float(o.log[n-1,5]),kw,sigma,chi)
  self.t.append(float(o.log[n-1,0])); self.kw.append(kw); self.kr.append(triple[0])
  self.p.append(self.memory.P if self.memory.active else np.nan); self.beta.append(beta)
  return triple

def force_metrics(t,fax):
 warm=t>=1.
 f,tt=fax[warm],t[warm]
 w=12
 hp=(f-np.convolve(f,np.ones(w)/w,mode='same'))[w:-w]
 tw=tt[w:-w]
 A=np.c_[np.cos(2*np.pi*20.8*tw),np.sin(2*np.pi*20.8*tw)]
 coef=np.linalg.lstsq(A,hp,rcond=None)[0]
 return dict(osc=TL.D.rms(hp-A@coef),ripple=float(np.hypot(*coef)/math.sqrt(2.)))

def score(row,arrays):
 t,fax,kt,te,kr=[arrays[s] for s in ('t','fax','k_eff','t_est','k_ret')]
 ids=np.rint(te/TL.V3.TC).astype(int)
 assert np.max(np.abs(t[ids]-te),initial=0.)<1e-9
 proxy=np.zeros((len(t),10)); proxy[:,0],proxy[:,7],proxy[:,9]=t,fax,kt
 _,masks=TL.D.masks_for(proxy,row)
 e=np.full(len(te),np.nan)
 valid=np.isfinite(kr)&(kr>0)&(kt[ids]>0)
 e[valid]=np.log(kr[valid]/kt[ids][valid])
 full_e=np.full(len(t),np.nan); full_e[ids]=e
 row['sustained_underestimate']=TL.D.sustained(full_e,masks['outside_event'])
 row['identification']={s:TL.error_stats(e,masks[s][ids]) for s in SEG}
 # Hold last finite return on control timeline for accumulated K5 duration.
 held=np.full(len(t),np.nan); held[ids]=kr
 ix=np.maximum.accumulate(np.where(np.isfinite(held),np.arange(len(t)),-1))
 ok=(ix>=0)&(t>.5)&(kt>0)
 ratio=held[ix[ok]]/kt[ok]
 row['over3_seconds']=float(np.count_nonzero(ratio>3.)*TL.V3.TC)
 row['max_ratio']=float(np.max(ratio)) if len(ratio) else None
 for k,v in force_metrics(t,fax).items():
  if k in row: assert abs(row[k]-v)<1e-12
  row[k]=v
 crossing=fax[masks['fast']]-row['FD']
 row['crossing_force_rms']=TL.D.rms(crossing) if len(crossing) else None
 row['_score']=dict(errors={s:e[masks[s][ids]&np.isfinite(e)].tolist() for s in SEG},crossing=crossing.tolist())

def job(task,output):
 started=time.perf_counter()
 TL.O.setup()
 sc=TL.SE.MC.scene_for(tuple(task['cell']),TL.CC.BASE,task['variant'])
 pol=TL.B.build('Huber-D4',TL.O.chosen()['Huber-D4'],sc)
 TL.CC.own_mechanism(pol)
 original=pol.estimator.estimate
 if task['arm']=='band_kf':
  # Disable only REL_MAX rejection within this isolated process.
  PB.REL_MAX=math.inf
  band=PB.BandEstimate(huber=False); est=pol.estimator
  original=lambda o,n,p: band(est,o,n,p)
 rec=Recorder(original,task['arm']); pol.estimator.estimate=rec
 log=TL.P3.episode(sc,pol,task['seed'])
 if log is None: raise RuntimeError(str(task))
 metric=TL.M.metric(log,pol,sc)
 row=dict(task,FD=float(sc['FD']),speed_mm_s=float(sc['v']*1000),event_info=sc['event_info'],L=float(np.log1p(metric['J'])),
          **{k:float(metric[k]) for k in TL.METRICS},filter_resets=rec.memory.resets)
 arrays=dict(t=log[:,0],fax=log[:,18],k_eff=log[:,26],t_est=np.asarray(rec.t),k_w=np.asarray(rec.kw),
             k_ret=np.asarray(rec.kr),P=np.asarray(rec.p),beta=np.asarray(rec.beta))
 score(row,arrays)
 if task['id'] in ('7','36','37','N4','N5','V4','V5'):
  row['file']=f"{task['id']}_{task['arm']}.npz"
  with (Path(output)/row['file']).open('xb') as stream: np.savez_compressed(stream,**arrays)
 row['runtime_s']=time.perf_counter()-started
 return row

def subprocess_job(task,output,spool):
 stem=f"{task['id']}_{task['arm']}"
 spec,result,stdout=[spool/(stem+s) for s in ('.task.json','.json','.stdout')]
 spec.write_text(json.dumps(task),encoding='utf-8')
 cmd=[sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--project-root',str(PROJECT),
      '--worker-task',str(spec),'--worker-output',str(output),'--worker-result',str(result)]
 with stdout.open('x',encoding='utf-8') as stream:
  done=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=stream,stderr=stream,
                      creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
 if done.returncode: raise RuntimeError(stdout.read_text(encoding='utf-8'))
 row=json.loads(result.read_text(encoding='utf-8')); result.unlink()
 return row

def summarize(rows):
 groups=('A','A_slow_thin','B','A+B','A_large_original','A_large_fresh')
 def agg(rr):
  out={}
  for arm in ARMS:
   aa=[r for r in rr if r['arm']==arm]
   ss=[r['sustained_underestimate'] for r in aa if r['sustained_underestimate'] is not None]
   cf=np.concatenate([r['_score']['crossing'] for r in aa])
   out[arm]=dict(mean_L=float(np.mean([r['L'] for r in aa])),median_L=float(np.median([r['L'] for r in aa])),
    identification={s:TL.D.stats(np.concatenate([r['_score']['errors'][s] for r in aa])) for s in SEG},
    sustained=dict(median=float(np.median(ss)) if ss else None,worst=max(ss) if ss else None),
    mean_ripple=float(np.mean([r['ripple'] for r in aa])),crossing_force_rms=TL.D.rms(cf) if len(cf) else None)
  return out
 summary={g:agg([r for r in rows if r['group']==g]) for g in groups}
 for g,indices in (('A_large_fast',TL.FAST),('A_large_slow',(38,41))):
  summary[g]=agg([r for r in rows if r['block']==0 and r['index'] in indices])
 summary['original_pooled']=agg([r for r in rows if r['block']==0])
 by={(r['id'],r['arm']):r for r in rows}
 k1={g:summary[g]['kf']['mean_L'] for g in groups[-2:]}
 d3={str(i):by[str(i),'kf']['L']-by[str(i),'huber']['L'] for i in range(24) if i!=7}
 large=[r for r in rows if r['group'] in groups[-2:]]
 d4=float(np.mean([r['L'] for r in large if r['arm']=='kf_noT'])-np.mean([r['L'] for r in large if r['arm']=='kf']))
 worst={a:max((r for r in rows if r['arm']==a),key=lambda r:(r['over3_seconds'],r['max_ratio'])) for a in ARMS}
 d6={g:float(np.mean([r['L'] for r in rows if r['group']==g and r['id']!='7' and r['arm']=='band_kf'])-
             np.mean([r['L'] for r in rows if r['group']==g and r['id']!='7' and r['arm']=='kf'])) for g in groups}
 pool=summary['original_pooled']
 reductions={a:1.-pool[a]['identification']['plateau']['median_abs']/pool['huber']['identification']['plateau']['median_abs'] for a in ('kf','kf_long')}
 under={g:summary[g]['kf']['sustained']['median']-summary[g]['huber']['sustained']['median'] for g in (*groups,'A_large_fast','A_large_slow')}
 verdicts=dict(
 K1=dict(met=k1[groups[-2]]<=.195 and k1[groups[-1]]<=.265,mean_L=k1),
 K2=dict(met=by['7','kf']['L']<.5,L=by['7','kf']['L']),
 K3=dict(met=all(abs(v)<=.02 for v in d3.values()),per_case_delta=d3),
 K4=dict(met=d4>=.03,mean_L_delta=d4),
 K5=dict(met=all(r['over3_seconds']<=.5 for r in rows),worst={a:{k:r[k] for k in ('id','over3_seconds','max_ratio')} for a,r in worst.items()}),
 K6=dict(met=by['7','band_kf']['L']<=by['7','kf']['L'] and all(abs(v)<=.02 for v in d6.values()),
         case7={a:by['7',a]['L'] for a in ('kf','band_kf')},elsewhere_group_delta=d6),
 K7=dict(met=reductions['kf']>=.15 and reductions['kf_long']>=.35 and all(v<0 for v in under.values()),
         plateau_reduction=reductions,sustained_median_delta=under))
 return summary,verdicts

def main():
 parser=argparse.ArgumentParser(); parser.add_argument('--project-root',type=Path,required=True)
 parser.add_argument('--output-root',type=Path); parser.add_argument('--worker-task',type=Path)
 parser.add_argument('--worker-output',type=Path); parser.add_argument('--worker-result',type=Path)
 args=parser.parse_args()
 if args.worker_task:
  row=job(json.loads(args.worker_task.read_text(encoding='utf-8')),args.worker_output)
  with args.worker_result.open('x',encoding='utf-8') as f: json.dump(row,f,allow_nan=False)
  return
 outbase=args.output_root/'experiments/rossi_failure_study/results'
 dest=outbase/'memory_transport_loop.json'; output=outbase/'memory_transport'
 data=json.loads(dest.read_text(encoding='utf-8')); old=data['rows']
 assert data['status']=='complete' and not any(r['arm']=='tmem24' for r in old)
 predictions={
 'M1':'tmem24 original A_large mean L <= 0.195, fresh <= 0.265, each A_large L - delag2 <= 0.05',
 'M2':'tmem24 case 7 and fresh V4 L < 0.5',
 'M3':'tmem24 within +/-0.02 of huber for each original 0-23 except 7',
 'M4':'combined original and fresh A_large mean L(mem24)-L(tmem24) >= 0.03',
 'M5':'each new arm and run: total time k_ret/k_eff > 3 after 0.5 s <= 0.5 s',
 'M6':'pooled original plateau median |e| reduction >=15% tmem24 and >=25% tmem48; original fast A_large fall mean e <=+0.15 tmem24'}
 data['round2']={'status':'running','prespecified_predictions':predictions,'workers':8}
 dest.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
 print('Predictions BEFORE RUN '+json.dumps(predictions),flush=True)
 tasks=[]
 for r in old:
  if r['arm']!='huber': continue
  task={k:r[k] for k in ('block','id','index','group','seed','cell','variant')}
  if 'scene' in r: task['scene']=r['scene']
  tasks.extend(dict(task,arm=a) for a in ARMS[3:])
 assert len(tasks)==126
 h0check=compare_h0(dict(tasks[0],arm='h0'))
 print('H=0 check '+json.dumps(h0check),flush=True)
 spool=Path(tempfile.mkdtemp(prefix='finite_memory_workers_')); print('Worker IPC: '+str(spool),flush=True)
 rows=[]; started=time.perf_counter()
 with ThreadPoolExecutor(8) as pool:
  futures=[pool.submit(subprocess_job,t,output,spool) for t in tasks]
  for i,future in enumerate(as_completed(futures),1):
   row=future.result(); rows.append(row); print(f"{i}/126 {row['id']} {row['arm']} L={row['L']:.9f}",flush=True)
 groups=('A','A_slow_thin','B','A+B','A_large_original','A_large_fresh','A_large_fast','A_large_slow','original_pooled')
 summary={g:{a:data['summary'][g][a] for a in ARMS[:3]} for g in groups}
 for g in groups:
  rr=[r for r in rows if (r['group']==g if g in groups[:6] else r['block']==0 and (g=='original_pooled' or r['index'] in (TL.FAST if g=='A_large_fast' else (38,41))))]
  for arm in ARMS[3:]:
   aa=[r for r in rr if r['arm']==arm]; ss=[r['sustained_underestimate'] for r in aa if r['sustained_underestimate'] is not None]
   cf=np.concatenate([r['_score']['crossing'] for r in aa])
   summary[g][arm]=dict(mean_L=float(np.mean([r['L'] for r in aa])),median_L=float(np.median([r['L'] for r in aa])),
    identification={s:TL.D.stats(np.concatenate([r['_score']['errors'][s] for r in aa])) for s in SEG},
    sustained=dict(median=float(np.median(ss)) if ss else None,worst=max(ss) if ss else None),
    mean_ripple=float(np.mean([r['ripple'] for r in aa])),crossing_force_rms=TL.D.rms(cf) if len(cf) else None)
 by={(r['id'],r['arm']):r for r in old+rows}
 large=[r for r in rows if r['group'] in ('A_large_original','A_large_fresh')]
 deltas={r['id']:r['L']-by[r['id'],'delag2']['L'] for r in large if r['arm']=='tmem24'}
 d3={str(i):by[str(i),'tmem24']['L']-by[str(i),'huber']['L'] for i in range(24) if i!=7}
 d4=float(np.mean([r['L'] for r in large if r['arm']=='mem24'])-np.mean([r['L'] for r in large if r['arm']=='tmem24']))
 worst={a:max((r for r in rows if r['arm']==a),key=lambda r:(r['over3_seconds'],r['max_ratio'])) for a in ARMS[3:]}
 red={a:1-summary['original_pooled'][a]['identification']['plateau']['median_abs']/summary['original_pooled']['huber']['identification']['plateau']['median_abs'] for a in ('tmem24','tmem48')}
 fall=summary['A_large_fast']['tmem24']['identification']['fall']['mean']
 means={g:summary[g]['tmem24']['mean_L'] for g in ('A_large_original','A_large_fresh')}
 verdicts={
 'M1':dict(met=means['A_large_original']<=.195 and means['A_large_fresh']<=.265 and max(deltas.values())<=.05,mean_L=means,max_case_delta=max(deltas.items(),key=lambda x:x[1]),per_case_delta=deltas),
 'M2':dict(met=by['7','tmem24']['L']<.5 and by['V4','tmem24']['L']<.5,L={i:by[i,'tmem24']['L'] for i in ('7','V4')}),
 'M3':dict(met=all(abs(v)<=.02 for v in d3.values()),per_case_delta=d3),
 'M4':dict(met=d4>=.03,mean_L_delta=d4),
 'M5':dict(met=all(r['over3_seconds']<=.5 for r in rows),worst={a:{k:r[k] for k in ('id','over3_seconds','max_ratio')} for a,r in worst.items()}),
 'M6':dict(met=red['tmem24']>=.15 and red['tmem48']>=.25 and fall<=.15,plateau_reduction=red,fast_fall_mean=fall)}
 for row in rows: del row['_score']
 data['rows'].extend(rows)
 data['round2'].update(status='complete',new_episodes=126,reused_rows=126,summary=summary,predictions=verdicts,h0_check=h0check,runtime_s=time.perf_counter()-started)
 dest.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
 print('Verdicts '+json.dumps(verdicts),flush=True)

def compare_h0(task):
 TL.O.setup(); sc=TL.SE.MC.scene_for(tuple(task['cell']),TL.CC.BASE,task['variant'])
 pol=TL.B.build('Huber-D4',TL.O.chosen()['Huber-D4'],sc); TL.CC.own_mechanism(pol)
 rec=Recorder(pol.estimator.estimate,'h0'); pol.estimator.estimate=rec
 TL.P3.episode(sc,pol,task['seed'])
 path=PROJECT/'experiments/rossi_failure_study/results/trajectory_loop/case_00_delag2.npz'
 if not path.exists(): path=Path('D:/pip_env/temp/trajectory_robust_round4/experiments/rossi_failure_study/results/trajectory_loop/case_00_delag2.npz')
 with np.load(path,allow_pickle=False) as f: expected=f['k_ret']; times=f['t_est']
 assert np.array_equal(np.asarray(rec.t),times)
 actual=np.asarray(rec.kr); same_nan=np.array_equal(np.isnan(actual),np.isnan(expected))
 delta=np.abs(actual-expected); finite=np.isfinite(delta); max_abs=float(np.max(delta[finite]))
 return dict(case=task['id'],n=len(actual),same_nan=same_nan,max_abs=max_abs,met=bool(same_nan and max_abs<=1e-12),trace=str(path))
if __name__=='__main__': main()
