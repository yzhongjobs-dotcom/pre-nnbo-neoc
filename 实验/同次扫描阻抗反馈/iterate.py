"""Second design round; joint force/pose is original-task cost, not a test retune.

Keeps force-only and force+normal-pose tasks separately. All model/weight choices
are exploratory on train/validation; independent tests only after protocol freeze.
"""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import model as m
from run import noises,batch,NAMES
HERE=Path(__file__).resolve().parent
ALL=NAMES+['transport_fixed','task_matched','task_matched_transport','temporal_fixed_ff']


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--task',choices=['force','joint'],required=True);args=ap.parse_args()
    q=float(args.task=='joint')
    ps=m.sample(18,8101);ns=noises(ps,8102)
    configs=[]
    for kt in [300.,1000.]:
      for kn in [0.,20.,80.,240.]:
       for dn in [2.,6.,18.,54.]:
        for kp,ki in [(0.,0.),(.5,0.),(1.5,0.),(.5,8.)]:
         configs.append(np.array([kt,kn,2*.8*np.sqrt(m.MASS[0]*kt),dn,kp,ki,.6,.008,q,1.]))
    m.episode(ps[0],ns[0],0,configs[0]);start=time.perf_counter();rows=[]
    for i,cfg in enumerate(configs):
        values,_,_,_=batch(ps,ns,0,cfg);rows.append((float(values[:,0].mean()),cfg))
        if i%48==0:print(args.task,'fixed tuning',i,min(v[0] for v in rows),flush=True)
    cfg=min(rows,key=lambda v:v[0])[1];cfgs=[]
    for method in range(len(ALL)):
        cc=cfg.copy()
        if method in (8,10,11):
            tries=[]
            for scale in [0.,.25,.5,.75,1.,1.25]:
                cf=cfg.copy();cf[9]=scale
                values,_,_,_=batch(ps,ns,method,cf);tries.append((float(values[:,0].mean()),scale))
            cc[9]=min(tries)[1]
        cfgs.append(cc)
    vp=m.sample(18,8201);vn=noises(vp,8202)
    summary={};arrays={'train_params':ps,'validation_params':vp,'cfgs':np.array(cfgs)}
    for method,name in enumerate(ALL):
        values,traces,x,cost=batch(vp,vn,method,cfgs[method],collect=method==6)
        arrays[name]=values
        pose=np.array([np.sqrt(np.mean((tr[:,2]-tr[0,2])**2)) for tr in traces])
        summary[name]={'mean_J':float(values[:,0].mean()),'force_rms_N':float(values[:,1].mean()),
          'normal_pose_rms_m':float(pose.mean()),'tangent_rms_m':float(values[:,2].mean()),
          'worst_peak_force_N':float(values[:,3].max()),'energy_active_fraction':float(values[:,5].mean())}
        if method==6:
            arrays['validation_features']=x;arrays['validation_costs']=cost
            for i in range(3):arrays[f'trace_{i}']=traces[i]
        print(args.task,name,summary[name],flush=True)
    payload={'stage':'Second exploratory design iteration, no independent test',
      'task':args.task,'normal_pose_weight':q,'normal_pose_scale_m':.0015,'fixed_configs':len(configs),
      'cfgs':{name:cfgs[i].tolist() for i,name in enumerate(ALL)},'train_seed':8101,'validation_seed':8201,
      'summary':summary,'seconds':time.perf_counter()-start,
      'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','iterate.py','run.py']}}
    (HERE/f'{args.task}_validation.json').write_text(json.dumps(payload,indent=2)+'\n')
    np.savez_compressed(HERE/f'{args.task}_validation.npz',**arrays)
if __name__=='__main__':main()
