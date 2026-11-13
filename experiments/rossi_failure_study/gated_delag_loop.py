import os,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'; sys.dont_write_bytecode=True
from pathlib import Path
import json,tempfile,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
ROOT=Path(r'E:\github_project\scan-impedance-update')
BASE=ROOT/'experiments/rossi_failure_study/results'
DEST=Path(__file__).parent/'results/gated_delag_loop.json'
os.environ['TRAJECTORY_LOOP_PROJECT']=str(ROOT)
sys.path.insert(0,str(ROOT/'experiments/rossi_failure_study'))
import trajectory_loop as T
P={'G1':'gated L < 0.5 on case 7 and fresh V4','G2':'gated_delag2 improves mean L >= 0.05 in each A_large group, and no A_large scene worsens > 0.05','G3':'each case 0-23 gated_delag2 within +/-0.01 of gated','G4':'gated_delag2 lowest mean L over all 42 scenes'}
SEG=('plateau','rise','fall','reversal')
def job(task):
 T.O.setup(); sc=T.SE.MC.scene_for(tuple(task['cell']),T.CC.BASE,task['variant']); pol=T.B.build('Huber-D4',T.O.chosen()['Huber-D4'],sc)
 if task['arm']=='gated_delag2': pol.estimator.estimate=T.ReturnedEstimate(pol.estimator.estimate,'delag2',task['index'])
 step=pol.step; rec=[]
 def capture(o):
  u=step(o); rec.append((float(o.log[o.n-1,0]),pol.k_good)); return u
 pol.step=capture; log=T.P3.episode(sc,pol,task['seed'])
 if log is None: raise RuntimeError(str(task))
 m=T.M.metric(log,pol,sc); row=dict(task,L=float(np.log1p(m['J'])),**{k:float(m[k]) for k in ('rmse','peak','over_15','loss_time','accept')})
 rec=np.asarray(rec); assert len(rec)==len(log) and np.allclose(rec[:,0],log[:,0],rtol=0,atol=1e-9)
 proxy=np.zeros((len(log),10)); proxy[:,0],proxy[:,7],proxy[:,9]=log[:,0],log[:,18],log[:,26]
 _,masks=T.D.masks_for(proxy,{'event_info':sc['event_info']}); valid=np.isfinite(rec[:,1])&(rec[:,1]>0)&(log[:,26]>0)
 e=np.full(len(log),np.nan); e[valid]=np.log(rec[valid,1]/log[valid,26])
 row['identification']={s:T.error_stats(e,masks[s]) for s in SEG}; row['sustained_underestimate']=T.D.sustained(e,masks['outside_event'])
 row['_errors']={s:e[masks[s]&np.isfinite(e)].tolist() for s in SEG}; return row
def subjob(t,spool):
 name=t['id']+'_'+t['arm']; a,b,c=[spool/(name+x) for x in ('.task','.result','.stdout')]; a.write_text(json.dumps(t),encoding='utf8')
 with c.open('w',encoding='utf8') as f: p=subprocess.run([sys.executable,'-X','utf8',__file__,'--worker',str(a),str(b)],stdout=f,stderr=f,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
 if p.returncode: raise RuntimeError(c.read_text(encoding='utf8'))
 return json.loads(b.read_text(encoding='utf8'))
def run(tasks,spool):
 rows=[]
 with ThreadPoolExecutor(8) as ex:
  fs=[ex.submit(subjob,t,spool) for t in tasks]
  for f in as_completed(fs):
   r=f.result(); rows.append(r); print(f"{len(rows)}/{len(tasks)} {r['id']} {r['arm']} L={r['L']:.12f}",flush=True)
 return rows
def main():
 if '--worker' in sys.argv:
  a,b=map(Path,sys.argv[2:4]); b.write_text(json.dumps(job(json.loads(a.read_text(encoding='utf8'))),allow_nan=False),encoding='utf8'); return
 if DEST.exists(): raise FileExistsError(DEST)
 source=json.loads((BASE/'memory_transport_loop.json').read_text(encoding='utf8'))['rows']; orig=[r for r in source if r['arm']=='huber']; assert len(orig)==42
 tasks=[{**{k:r[k] for k in ('block','id','index','group','seed','cell','variant')},**({'scene':r['scene']} if 'scene' in r else {}),'arm':'gated'} for r in orig]
 out={'status':'running','prespecified_predictions':P,'workers':8,'rows':[]}; DEST.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8'); print('PREDICTIONS BEFORE RUN '+json.dumps(P),flush=True)
 spool=Path(tempfile.mkdtemp(prefix='gated_delag_')); first=run([t for t in tasks if t['block']==0],spool)
 ref={r['index']:r for r in json.loads((BASE/'completion_closed_loop.json').read_text(encoding='utf8'))['rows'] if r['method']=='Huber-D4'}
 diff={r['id']:r['L']-ref[r['index']]['L'] for r in first}; check={'max_abs_dL':max(abs(v) for v in diff.values()),'per_case_dL':diff}; print('HARNESS '+json.dumps(check),flush=True)
 if check['max_abs_dL']!=0:
  out.update(status='harness_failed',harness=check,rows=[{k:v for k,v in r.items() if k!='_errors'} for r in first]); DEST.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf8'); return
 new=first+run([t for t in tasks if t['block']!=0]+[dict(t,arm='gated_delag2') for t in tasks],spool); ident={}
 for group,selected in [('A_large_fast',lambda r:r['block']==0 and r['index'] in T.FAST),('A_large_fresh',lambda r:r['group']=='A_large_fresh')]:
  ident[group]={}
  for arm in ('gated','gated_delag2'):
   rr=[r for r in new if r['arm']==arm and selected(r)]; sust=[r['sustained_underestimate'] for r in rr if r['sustained_underestimate'] is not None]
   ident[group][arm]={'segments':{s:T.D.stats(np.concatenate([r['_errors'][s] for r in rr])) for s in SEG},'sustained':{'median':float(np.median(sust)) if sust else None,'worst':max(sust) if sust else None}}
 reused=[{k:r[k] for k in ('arm','id','index','group','block','seed','cell','variant','L','rmse','peak','over_15','loss_time')}|({'scene':r['scene']} if 'scene' in r else {})|{'accept':None} for r in source if r['arm'] in ('huber','delag2')]; assert len(reused)==84
 for r in new: del r['_errors']
 out.update(status='complete',harness=check,identification=ident,rows=reused+new); DEST.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf8'); print('COMPLETE',flush=True)
if __name__=='__main__': main()
