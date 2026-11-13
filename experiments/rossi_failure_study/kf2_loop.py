import os,sys,json,math,tempfile,subprocess,time
from pathlib import Path
for v in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[v]='1'
sys.dont_write_bytecode=True
ROOT=Path(r'E:\github_project\scan-impedance-update'); HERE=Path(__file__).parent; OUT=HERE/'results'; BASE=ROOT/'experiments/rossi_failure_study/results'; sys.path.insert(0,str(BASE.parent)); sys.argv.extend(['--project-root',str(ROOT)])
import memory_transport_loop as M
import numpy as np
from concurrent.futures import ThreadPoolExecutor,as_completed
T=M.TL
PRED={'D1':'case 7 and V4 L < 0.5','D2':'original and fresh A_large means <= 0.20 and 0.25; each delta vs delag2 <= 0.05','D3':'each original 0-23 except 7 and fresh A within +/-0.02 huber','D4':'42-scene mean L < 0.18','D5':'every run ratio > 3 after 0.5 s totals <= 0.5 s','D6':'each of seven scenes at +/-0.005 within 0.05 L of unperturbed'}
class Filter:
 def __init__(self): self.resets=0; self.reset()
 def reset(self):
  self.g=T.DelagCorrection(force_ratio_factor=False); self.n=None; self.active=False; self.missing=0; self.prev=None; self.x=np.zeros(2); self.P=np.full((2,2),np.nan)
 def step(self,n,y,f,kw,sig,chi):
  if self.n is not None and n!=self.n+1: self.resets+=1; self.reset()
  self.n=n; g=self.g; j=len(g.forces); ys=g.y_sums[-1]+y; fs=g.f_sums[-1]+f; a=max(0,j-5); b=max(0,j-23)
  y6=(ys-g.y_sums[a])/(j+1-a); y24=(ys-g.y_sums[b])/(j+1-b); F=(fs-g.f_sums[b])/(j+1-b); beta=g.beta
  finite=np.isfinite(kw); rel=sig/kw if finite and kw>0 else np.inf; self.missing=0 if finite else self.missing+1
  if self.missing>25: self.active=False
  if self.active:
   if self.prev is not None and F>=1:
    u=-math.log(float(np.clip(1+beta*(y6-self.prev)/(F/math.exp(self.x[0])),.8,1.25))); s=float(np.clip(1+self.x[1],0,2)); self.x[0]+=s*u
    fm=np.array([[1.,u],[0.,1.]]); self.P=fm@self.P@fm.T+np.diag([5e-4,1e-5])
   if finite and kw>0 and np.isfinite(sig) and rel<=3:
    cw=F/kw; d=math.log(float(np.clip(1+beta*(y6-y24)/cw,.5,2))) if cw!=0 else 0.; s=float(np.clip(1+self.x[1],0,2)); nu=math.log(kw)-self.x[0]-s*d; h=np.array([1.,d]); R=24*rel**2; S=h@self.P@h+R
    if S>0:
     K=self.P@h/S; self.x+=K*nu; a=np.eye(2)-np.outer(K,h); self.P=a@self.P@a.T+R*np.outer(K,K); self.x[1]=np.clip(self.x[1],-1,1)
  elif finite and kw>0 and np.isfinite(sig) and rel<=1:
   self.x=np.array([math.log(kw),0.]); self.P=np.diag([(2*rel)**2,.5**2]); self.active=True
  g.step(n,y,f,kw,sig); self.prev=y6
  if not self.active: return (kw,sig,chi),beta
  k=float(np.clip(math.exp(self.x[0]),200,3e5)); return (k,k*math.sqrt(max(0,self.P[0,0])),chi),beta
class Recorder:
 def __init__(self,orig,e): self.orig=orig; self.e=e; self.memory=Filter(); self.data=[]
 def __call__(self,o,n,p):
  kw,sig,chi=self.orig(o,n,p); triple,beta=self.memory.step(n,float(o.log[n-1,2]),float(o.log[n-1,5]),kw,sig,chi); triple=(triple[0]*(1+self.e),triple[1],triple[2]); q=self.memory
  self.data.append((float(o.log[n-1,0]),kw,triple[0],q.P[0,0] if q.active else np.nan,q.x[1] if q.active else np.nan,beta)); return triple
def job(task):
 T.O.setup(); sc=T.SE.MC.scene_for(tuple(task['cell']),T.CC.BASE,task['variant']); pol=T.B.build('Huber-D4',T.O.chosen()['Huber-D4'],sc); T.CC.own_mechanism(pol); rec=Recorder(pol.estimator.estimate,task['perturbation']); pol.estimator.estimate=rec
 log=T.P3.episode(sc,pol,task['seed']); m=T.M.metric(log,pol,sc); a=np.asarray(rec.data); row=dict(task,FD=float(sc['FD']),event_info=sc['event_info'],L=float(np.log1p(m['J'])),**{k:float(m[k]) for k in ('rmse','peak','over_15','loss_time')},filter_resets=rec.memory.resets)
 arrays=dict(t=log[:,0],fax=log[:,18],k_eff=log[:,26],t_est=a[:,0],k_w=a[:,1],k_ret=a[:,2],P00=a[:,3],eps=a[:,4],beta=a[:,5]); M.score(row,arrays); vals=a[(a[:,0]>3)&np.isfinite(a[:,4]),4]; row['eps_median_after_3s']=float(np.median(vals)) if len(vals) else None
 if task['id'] in ('7','36','37','N4','N5','V4','V5') and task['perturbation']==0:
  path=OUT/'kf2'/f"{task['id']}_kf2.npz"; np.savez_compressed(path,**arrays); row['file']=str(path.relative_to(OUT))
 return row
def sub(task,spool):
 stem=task['id']+'_'+str(task['perturbation']); spec=spool/(stem+'.task'); result=spool/(stem+'.json'); output=spool/(stem+'.log'); spec.write_text(json.dumps(task),encoding='utf8')
 with output.open('w',encoding='utf8') as f: p=subprocess.run([sys.executable,'-X','utf8',__file__,'--worker',str(spec),str(result)],stdout=f,stderr=f,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
 if p.returncode: raise RuntimeError(output.read_text(encoding='utf8'))
 return json.loads(result.read_text(encoding='utf8'))
def main():
 if '--worker' in sys.argv:
  i=sys.argv.index('--worker'); task=json.loads(Path(sys.argv[i+1]).read_text(encoding='utf8')); Path(sys.argv[i+2]).write_text(json.dumps(job(task),allow_nan=False),encoding='utf8'); return
 OUT.mkdir(exist_ok=True); dest=OUT/'kf2_loop.json'; traces=OUT/'kf2'; assert not dest.exists() and not traces.exists(); traces.mkdir(); dest.write_text(json.dumps({'status':'running','prespecified_predictions':PRED,'rows':[]},indent=2),encoding='utf8'); print('PRE-RUN PREDICTIONS',json.dumps(PRED),flush=True)
 old=json.loads((BASE/'memory_transport_loop.json').read_text(encoding='utf8'))['rows']; fresh=json.loads((BASE/'chi_hold_loop.json').read_text(encoding='utf8'))['rows']; ref=[r for r in old if r['arm']=='huber']+[r for r in fresh if r['arm']=='huber' and r['group']=='A_fresh']; assert len(ref)==54
 tasks=[dict({k:r[k] for k in ('block','id','index','group','seed','cell','variant','scene') if k in r},arm='kf2',perturbation=0.) for r in ref]; tasks += [dict(t,perturbation=e) for t in tasks if t['id'] in ('7','36','37','N4','N5','V4','V5') for e in (.005,-.005)]; assert len(tasks)==68
 spool=Path(tempfile.mkdtemp(prefix='kf2_workers_')); rows=[]; start=time.time()
 with ThreadPoolExecutor(8) as pool:
  for f in as_completed([pool.submit(sub,t,spool) for t in tasks]):
   row=f.result(); rows.append(row); print(len(rows),'/68',row['id'],row['perturbation'],round(row['L'],6),flush=True)
 dest.write_text(json.dumps({'status':'complete','prespecified_predictions':PRED,'runtime_s':time.time()-start,'rows':rows},indent=2,allow_nan=False),encoding='utf8')
if __name__=='__main__': main()
