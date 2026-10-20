"""Tune on training scans, inspect validation, then freeze an independent test."""
from pathlib import Path
import argparse,hashlib,json,time,platform
import numpy as np
import model as m
HERE=Path(__file__).resolve().parent
NAMES=['fixed_PI','observer_fixed','critical','transport_ff','temporal_joint',
       'spatial_gain','spatial_joint','spatial_ff_only']


def noises(ps,seed):
    rng=np.random.default_rng(seed)
    return [rng.normal(size=(int(round(.064/p[9]/m.CONTROL_DT))+1,3))*[.03,.03,.0001] for p in ps]


def batch(ps,ns,method,cfg,collect=False):
    values=[];traces=[];xs=[];qs=[]
    for p,n in zip(ps,ns):
        val,tr,x,q=m.episode(p,n,method,cfg,collect=collect)
        values.append(val);traces.append(tr)
        if collect:xs.append(x);qs.append(q)
    return np.array(values),traces,np.concatenate(xs) if xs else np.empty((0,16)),np.concatenate(qs) if qs else np.empty((0,25))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--test',action='store_true');args=ap.parse_args()
    ps=m.sample(18,8101);ns=noises(ps,8102)
    if args.test:
        config=json.loads((HERE/'validation.json').read_text());cfg=np.array(config['selected_cfg'])
        raise NotImplementedError('Freeze final protocol after the pilot; no test run yet')
    configs=[]
    # Common PD tangent is deliberately strong; normal direct force PI may use Kn=0.
    for kt in [300.,1000.]:
      for kn in [0.,150.]:
       for dn in [12.,35.,80.]:
        for kp,ki in [(0.,0.),(.5,8.),(1.5,15.),(3.,30.)]:
         dt=2*.8*np.sqrt(m.MASS[0]*kt)
         configs.append(np.array([kt,kn,dt,dn,kp,ki,.6,.008]))
    # Compile before timing or comparing configurations.
    m.episode(ps[0],ns[0],0,configs[0])
    start=time.perf_counter();rows=[]
    for j,cfg in enumerate(configs):
        values,_,_,_=batch(ps,ns,0,cfg)
        score=float(values[:,0].mean());rows.append((score,cfg))
        if j%24==0:print('tuning',j,'best',min(x[0] for x in rows),flush=True)
    _,cfg=min(rows,key=lambda r:r[0]);tune_seconds=time.perf_counter()-start
    valps=m.sample(18,8201);valns=noises(valps,8202)
    selected={};summary={};arrays={'train_params':ps,'validation_params':valps,'selected_cfg':cfg}
    for method,name in enumerate(NAMES):
        values,traces,x,q=batch(valps,valns,method,cfg,collect=method==6)
        arrays[name]=values
        summary[name]={'mean_J':float(values[:,0].mean()),'mean_force_rms_N':float(values[:,1].mean()),
            'mean_tangent_rms_m':float(values[:,2].mean()),'max_peak_force_N':float(values[:,3].max()),
            'mean_contact_loss_fraction':float(values[:,4].mean()),'mean_energy_active_fraction':float(values[:,5].mean()),
            'min_energy_margin':float(values[:,8].min()),'max_gain_violation':float(values[:,9].max())}
        print(name,summary[name],flush=True)
        if method==6:
            arrays['validation_features']=x;arrays['validation_model_costs']=q
            arrays['representative_trace']=traces[0]
    out={'scope':'Exploratory train/validation, not independent test or full B0/prior replication',
      'selected_cfg':cfg.tolist(),'configuration_order':['Kt','Kn','Dt','Dn','force_kp','force_ki','measurement_alpha','speed_overwritten_per_scan'],
      'train_seed':8101,'train_noise_seed':8102,'validation_seed':8201,'validation_noise_seed':8202,
      'train_scans':len(ps),'validation_scans':len(valps),'fixed_configs':len(configs),'tune_seconds':tune_seconds,
      'summary':summary,'model_sha256':hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
      'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime':platform.python_version()}
    (HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    np.savez_compressed(HERE/'validation_arrays.npz',**arrays)
    print('saved exploratory validation',flush=True)
if __name__=='__main__':main()
