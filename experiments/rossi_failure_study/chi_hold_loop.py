import os,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import json,tempfile,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
SOURCE=Path(r'E:\github_project\scan-impedance-update'); HERE=Path(__file__).parent; BASE=SOURCE/'experiments/rossi_failure_study/results'; os.environ['TRAJECTORY_LOOP_PROJECT']=str(SOURCE); sys.path.insert(0,str(BASE.parent))
import trajectory_loop as T
P={'C1':'case 7 and V4 L < 0.5','C2':'original large mean <= 0.20; fresh excluding V4 <= 0.27; each large delta <= 0.05','C3':'each original 0-23 except 7 and fresh A within +/-0.02 huber','C4':'42-scene mean below 0.220 and 0.228','C5':'each original 0-23 except 7 held fraction <5% and longest <0.5 s'}
class Est:
 def __init__(self,original,arm,limit): self.original,self.arm,self.limit=original,arm,limit; self.g=T.DelagCorrection(False); self.last=None; self.prev=None; self.data=[]
 def __call__(self,o,n,p):
  if self.prev is not None and n!=self.prev+1: self.last=None
  self.prev=n; kw,sig,chi=self.original(o,n,p); ok=bool(np.isfinite(kw) and np.isfinite(sig) and np.isfinite(chi) and kw>0 and chi<=self.limit and sig/kw<=1); beta=0.
  if self.arm=='delag2_chi':
   kc,beta=self.g.step(n,float(o.log[n-1,2]),float(o.log[n-1,5]),kw if ok else np.nan,sig if ok else np.nan)
   triple=(kc,sig*kc/kw if ok and kc!=kw and np.isfinite(sig) else sig,chi) if ok else (kw,sig,chi)
  else: triple=(kw,sig,chi)
  held=self.arm!='huber' and not ok and self.last is not None
  if ok: self.last=triple
  elif held: triple=self.last
  self.data.append((float(o.log[n-1,0]),kw,sig,chi,triple[0],held,beta)); return triple

def job(task):
 T.O.setup(); sc=T.SE.MC.scene_for(tuple(task['cell']),T.CC.BASE,task['variant']); pol=T.B.build('Huber-D4',T.O.chosen()['Huber-D4'],sc); limit=float(pol.params['chi2_max']); T.CC.own_mechanism(pol); rec=Est(pol.estimator.estimate,task['arm'],limit); pol.estimator.estimate=rec
 log=T.P3.episode(sc,pol,task['seed']); m=T.M.metric(log,pol,sc); a=np.array(rec.data); h=a[:,5].astype(bool); count=longest=0
 for v in h:
  count=count+1 if v else 0; longest=max(longest,count)
 row=dict(task,L=float(np.log1p(m['J'])),**{k:float(m[k]) for k in ('rmse','peak','over_15','loss_time')},held_fraction=float(h.mean()),longest_hold_s=float(longest*T.V3.TC),FD=float(sc['FD']),event_info=sc['event_info'],chi_max=limit)
 arrays=dict(t=log[:,0],fax=log[:,18],k_eff=log[:,26],t_est=a[:,0],k_w=a[:,1],sigma=a[:,2],chi=a[:,3],k_ret=a[:,4],held=h,beta=a[:,6]); ids=np.rint(a[:,0]/T.V3.TC).astype(int); proxy=np.zeros((len(log),10)); proxy[:,0],proxy[:,7],proxy[:,9]=log[:,0],log[:,18],log[:,26]; _,masks=T.D.masks_for(proxy,row); e=np.full(len(a),np.nan); valid=np.isfinite(a[:,4])&(a[:,4]>0)&(log[ids,26]>0); e[valid]=np.log(a[valid,4]/log[ids,26][valid]); score=(a[:,0]>1.5)&(log[ids,26]>0)&(log[ids,18]>.5)
 row['_errors']={s:e[masks[s][ids]&score&np.isfinite(e)].tolist() for s in ('plateau','rise','fall','reversal')}; full=np.full(len(log),np.nan); full[ids]=e; row['sustained_underestimate']=T.D.sustained(full,masks['outside_event'])
 if task['arm']=='delag2_chi' and task['id'] in ('7','36','37','V4','N4','N5'):
  path=HERE/'results/chi_hold'/f"{task['id']}_delag2_chi.npz"; np.savez_compressed(path,**arrays); row['file']=path.name
 return row

def sub(t,spool):
 stem=t['id']+'_'+t['arm']; spec=spool/(stem+'.task'); result=spool/(stem+'.json'); out=spool/(stem+'.log'); spec.write_text(json.dumps(t),encoding='utf8')
 with out.open('w',encoding='utf8') as f: p=subprocess.run([sys.executable,'-X','utf8',__file__,'--worker',str(spec),str(result)],stdout=f,stderr=f,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
 if p.returncode: raise RuntimeError(out.read_text(encoding='utf8'))
 return json.loads(result.read_text(encoding='utf8'))

def main():
 if '--worker' in sys.argv:
  a,b=map(Path,sys.argv[2:]); b.write_text(json.dumps(job(json.loads(a.read_text(encoding='utf8'))),allow_nan=False),encoding='utf8'); return
 base=HERE/'results'; dest=base/'chi_hold_loop.json'; out=base/'chi_hold'; assert not dest.exists() and not out.exists(); out.mkdir(parents=True)
 dest.write_text(json.dumps(dict(status='running',prespecified_predictions=P,workers=8),indent=2)+'\n',encoding='utf8'); print('PREDICTIONS BEFORE RUN '+json.dumps(P),flush=True)
 old=json.loads((BASE/'gated_delag_loop.json').read_text(encoding='utf8'))['rows']; ref=[r for r in old if r['arm']=='huber']; assert len(ref)==42
 tasks=[{k:r[k] for k in ('block','id','index','group','seed','cell','variant','scene') if k in r} for r in ref]; robust=json.loads((BASE/'trajectory_robust.json').read_text(encoding='utf8'))['rows']; fresh=[r for r in robust if r['block']==2 and r['arm']=='huber']; assert len(fresh)==12
 for i,r in enumerate(fresh): tasks.append(dict(block=2,id='A'+str(i),index=None,group='A_fresh',seed=r['seed'],cell=r['cell'],variant=r['variant'],scene=r['scene']))
 runs=[dict(t,arm=a) for t in tasks for a in ('huber_chi','delag2_chi')]+[dict(t,arm='huber') for t in tasks if t['block']==2]; spool=Path(tempfile.mkdtemp(prefix='chi_hold_workers_')); rows=[]
 with ThreadPoolExecutor(8) as pool:
  for f in as_completed([pool.submit(sub,t,spool) for t in runs]):
   row=f.result(); rows.append(row); print(f"{len(rows)}/{len(runs)} {row['id']} {row['arm']} L={row['L']:.6f}",flush=True)
 for i,h in enumerate(fresh):
  ident='A'+str(i); rr=next(r for r in rows if r['id']==ident and r['arm']=='huber'); assert rr['L']==h['L']; rr['source']='stored trajectory_robust; rerun matched'; rr['held_fraction']=rr['longest_hold_s']=None
  x=next(r for r in robust if r['block']==2 and r['arm']=='delag2' and r['scene']==h['scene'])
  rows.append(dict(arm='delag2',id=ident,block=2,group='A_fresh',index=None,scene=x['scene'],seed=x['seed'],cell=x['cell'],variant=x['variant'],**{k:x[k] for k in ('L','rmse','peak','over_15','loss_time')},held_fraction=None,longest_hold_s=None,source='stored trajectory_robust'))
 for r in old: r.setdefault('held_fraction',None); r.setdefault('longest_hold_s',None)
 dest.write_text(json.dumps(dict(status='complete',prespecified_predictions=P,workers=8,rows=old+rows),indent=2,allow_nan=False)+'\n',encoding='utf8')
if __name__=='__main__': main()
