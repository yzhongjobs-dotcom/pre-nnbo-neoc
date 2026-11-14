import os,sys,json,math,tempfile,subprocess,argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
for x in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[x]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'; sys.dont_write_bytecode=True
sys.path.insert(0,r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
import v2_trajectory_prior_loop as V
import numpy as np
TL=V.TL; BASE=V.BASE; DEST=Path(__file__).parent/'results/quantile_estimand_loop.json'
MET=('L','rmse','peak','over_15','loss_time')
def score(row,log,te,kr):
 t,fax,kt=log[:,0],log[:,18],log[:,26]; proxy=np.zeros((len(t),10)); proxy[:,0],proxy[:,7],proxy[:,9]=t,fax,kt
 _,m=TL.D.masks_for(proxy,row); ids=np.rint(te/TL.V3.TC).astype(int); assert np.max(np.abs(t[ids]-te),initial=0)<1e-9
 truth=np.interp(te,t,kt); force=np.interp(te,t,fax); valid=(te>1.5)&(truth>0)&(force>.5)
 e=np.full(len(te),np.nan); ok=np.isfinite(kr)&(kr>0)&(truth>0); e[ok]=np.log(kr[ok]/truth[ok]); seg={}
 for s in ('plateau','rise','fall'):
  v=e[m[s][ids]&valid&np.isfinite(e)]; seg[s]=dict(sum=float(v.sum()),n=len(v))
 full=np.full(len(t),np.nan); full[ids]=e
 return dict(segments=seg,sustained=TL.D.sustained(full,m['outside_event']))
def job(task):
 TL.O.setup(); sc=TL.SE.MC.scene_for(tuple(task['cell']),TL.CC.BASE,task['variant']); pol=TL.B.build('Huber-D4',TL.O.chosen()['Huber-D4'],sc); TL.CC.own_mechanism(pol)
 
 if task['arm'] in ('huber','delag2'):
  recorder=TL.ReturnedEstimate(pol.estimator.estimate,task['arm'],task['index']); pol.estimator.estimate=recorder
  log=TL.P3.episode(sc,pol,task['seed']); m=TL.M.metric(log,pol,sc)
  row={k:task[k] for k in ('block','scene','index','group','arm')}; row.update({k:float(m[k]) for k in ('rmse','peak','over_15','loss_time','ripple')}); row['L']=float(np.log1p(m['J'])); return row
 prior=V.PriorEstimate(pol.estimator); ratios=[]; returned=[]
 def wrap(o,n,p):
  k,sig,chi=prior(o,n,p)
  if np.isfinite(k) and k>0 and np.isfinite(sig):
   arm=task['arm']; x=sig/k;
   if task.get('collect'): ratios.append(math.exp(x))
   if arm=='v2q05': k=3e5 if .5*x>=math.log(3e5/k) else float(np.clip(k*math.exp(.5*x),200.,3e5))
   elif arm=='v2q10': k=3e5 if x>=math.log(3e5/k) else float(np.clip(k*math.exp(x),200.,3e5))
   elif arm=='v2c': k=float(np.clip(k*task['C'],200.,3e5))
   elif arm=='zero': k=float(np.clip(k*math.exp(0.*sig/k),200.,3e5))
  returned.append(k); return k,sig,chi
 pol.estimator.estimate=wrap; log=TL.P3.episode(sc,pol,task['seed']); m=TL.M.metric(log,pol,sc)
 row={k:task[k] for k in ('block','scene','index','group','arm')}; row.update({k:float(m[k]) for k in ('rmse','peak','over_15','loss_time','ripple')}); row['L']=float(np.log1p(m['J']))
 row['identification']=score(dict(row,event_info=sc['event_info'],FD=float(sc['FD'])),log,np.asarray(prior.arrays['t_est']),np.asarray(returned))
 if task.get('collect'): row['ratios']=ratios
 return row
def worker(t,spool):
 name=f"{t['block']}_{t['scene']}_{t['arm']}"; spec=spool/(name+'.in'); out=spool/(name+'.out'); err=spool/(name+'.err'); 
 if out.exists(): return json.loads(out.read_text(encoding='utf8'))
 spec.write_text(json.dumps(t),encoding='utf8')
 with err.open('w',encoding='utf8') as f: done=subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--worker',str(spec),'--out',str(out)],stdin=subprocess.DEVNULL,stdout=f,stderr=f,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
 if done.returncode: raise RuntimeError(name+' '+err.read_text(encoding='utf8'))
 return json.loads(out.read_text(encoding='utf8'))
def run(tasks,spool):
 rows=[]
 with ThreadPoolExecutor(8) as pool:
  for f in as_completed([pool.submit(worker,t,spool) for t in tasks]):
   r=f.result(); rows.append(r); print(f"{len(rows)}/{len(tasks)} {r['scene']} {r['arm']} {r['L']:.5f}",flush=True)
 return rows
def main():
 p=argparse.ArgumentParser(); p.add_argument('--fill-ripple',action='store_true'); p.add_argument('--resume-spool',type=Path); p.add_argument('--worker',type=Path); p.add_argument('--out',type=Path); a=p.parse_args()
 if a.worker: a.out.write_text(json.dumps(job(json.loads(a.worker.read_text(encoding='utf8'))),allow_nan=False),encoding='utf8'); return
 
 if a.fill_ripple:
  obj=json.loads(DEST.read_text(encoding='utf8')); stored=json.loads((BASE/'v2_trajectory_prior_loop.json').read_text(encoding='utf8'))['rows']; missing=[r for r in obj['rows'] if r['ripple'] is None]
  tasks=[next(x for x in stored if x['block']==r['block'] and x['scene']==r['scene'] and x['arm']==r['arm']) for r in missing]; spool=Path(tempfile.mkdtemp(prefix='quantile_ripple_')); print('ripple IPC',spool,flush=True)
  for rr in run(tasks,spool):
   r=next(x for x in missing if x['block']==rr['block'] and x['scene']==rr['scene'] and x['arm']==rr['arm']); assert all(r[k]==rr[k] for k in MET); r['ripple']=rr['ripple']
  obj['conventions']+=' Missing reference ripple recovered by exact-metric matched reruns.'; DEST.write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf8'); print('Filled',len(missing),flush=True); return
 if DEST.exists(): raise FileExistsError(DEST)
 stored=json.loads((BASE/'v2_trajectory_prior_loop.json').read_text(encoding='utf8'))['rows']; refs=[r for r in stored if r['block'] in (0,1,2) and r['arm'] in ('huber','delag2','v2t')]; assert len(refs)==162
 templates=[r for r in refs if r['arm']=='v2t']; assert len(templates)==54
 spool=a.resume_spool or Path(tempfile.mkdtemp(prefix='quantile_estimand_')); print('IPC',spool,flush=True)
 rerun=run([dict(r,collect=True) for r in templates if r['block']==0],spool)
 def verify(rows):
  for r in rows:
   ref=next(x for x in templates if x['block']==r['block'] and x['scene']==r['scene']); assert all(r[k]==ref[k] for k in MET),(r['scene'],r['arm'])
 verify(rerun); rerun.sort(key=lambda r:r['index']); ratios=np.concatenate([r.pop('ratios') for r in rerun]); C=float(ratios.mean()); print('C',C,'median',np.median(ratios),'p90',np.percentile(ratios,90),'ticks',len(ratios),flush=True)
 tasks=[dict(r,arm=arm,C=C) for r in templates for arm in ('v2q05','v2q10','v2c')]+[dict(r,arm='v2t') for r in templates if r['block'] in (1,2)]+[dict(r,arm='zero') for r in templates if r['block']==0 and r['index'] in (0,7,36)]
 new=run(tasks,spool); zeros=[r for r in new if r['arm']=='zero']; vr=rerun+[r for r in new if r['arm']=='v2t']; verify(vr)
 output=[]
 for r in refs:
  rr={k:r[k] for k in ('block','scene','index','group','arm','L','rmse','peak','over_15','loss_time')}
  if r['arm']=='v2t':
   extra=next(x for x in vr if x['block']==r['block'] and x['scene']==r['scene']); rr.update(ripple=extra['ripple'],identification=extra['identification'])
  elif 'ripple' in r: rr['ripple']=r['ripple']
  else:
   try:
    z=V.load(r); t,f=z['t'],z['fax']; warm=t>=1.; tt=t[warm]; f=f[warm]; w=12; hp=(f-np.convolve(f,np.ones(w)/w,mode='same'))[w:-w]; tt=tt[w:-w]; A=np.c_[np.cos(2*np.pi*20.8*tt),np.sin(2*np.pi*20.8*tt)]; rr['ripple']=float(np.hypot(*np.linalg.lstsq(A,hp,rcond=None)[0])/math.sqrt(2))
   except (FileNotFoundError,KeyError): rr['ripple']=None
  output.append(rr)
 output += [r for r in new if r['arm'] in ('v2q05','v2q10','v2c')]; assert len(output)==324
 check={str(r['index']):r['L']-next(x['L'] for x in templates if x['block']==0 and x['index']==r['index']) for r in zeros}
 obj=dict(status='complete',workers=8,C=C,margin_median=float(np.median(ratios)),margin_p90=float(np.percentile(ratios,90)),margin_ticks=len(ratios),zero_dL=check,rows=output,conventions='30 original v2t reruns for C, 24 fresh v2t reruns for identification/ripple, checked exactly against stored metrics; e=ln(k_ret/k_eff); original masks, pooled segment ticks and 13-tick sustained windows; no traces saved.')
 with DEST.open('x',encoding='utf8') as f: json.dump(obj,f,indent=2,allow_nan=False)
 print('Written',DEST,'zero',check,flush=True)
if __name__=='__main__': main()




