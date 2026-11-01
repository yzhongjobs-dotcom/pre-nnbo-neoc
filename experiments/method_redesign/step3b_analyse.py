"""Post-selection diagnostics; never changes fitted weights or candidate ranking."""
import argparse,json,pickle
from step3b_core import *
from step3b_predictors import eval_pred,write_md

def analyse(data,models):
    grouped={n:{} for n in models};truths=[]
    for e in data:
        assert e['p'][12]==0, 'benchmark has zero ridge bump; x does not affect surface'
        plant=use_plant(e['code']);p=e['p'];ds=[]
        for t in e['t']:
            h=plant.surface(BASE_PLANT.X0,t,p)[0]
            hf,ht,_=plant.surface(BASE_PLANT.X0,t+L,p)
            ds.append(hf+TAU*ht-h)
        truths.append(np.array(ds))
    for n,model in models.items():
        groups={}
        for e,dtruth in zip(data,truths):
            y,v=eval_pred(n,model,e);mask=e['t']>=1
            horizon=(y-e['z'])-dtruth
            inc=np.r_[np.nan,np.diff(y)-np.diff(e['target'])]
            imask=mask&~e['reset'];imask[0]=False
            err=y-e['target'];centered=err.copy()
            for s in np.unique(e['segment']):
                m=(e['segment']==s)&mask
                if np.any(m):centered[m]-=err[m].mean()
            for fam,phase in [('all','all'),(e['family'],'all')]+[('all',p) for p in ('breathing','pre50ms','ramp','hold','retract')]+[(e['family'],p) for p in ('breathing','pre50ms','ramp','hold','retract')]:
                m=mask if phase=='all' else mask&(e['phase']==phase);mi=m&imask
                if np.any(m):
                    g=groups.setdefault(fam+'|'+phase,[[],[],[]]);g[0].append(horizon[m]);g[1].append(inc[mi]);g[2].append(centered[m])
        for key,arrays in groups.items():
            result={}
            for metric,a in zip(('local_horizon','oracle_increment','segment_centered'),arrays):
                values=np.concatenate(a)
                result[metric]=dict(n_ticks=len(values),rmse_um=float(np.sqrt(np.mean(values**2))*1e6),p95_um=float(np.percentile(abs(values),95)*1e6))
            grouped[n][key]=result
    return grouped

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);args=ap.parse_args()
    cfg,models=load_models();models={'P0':None,**models,'P0r':MultiIMM(**cfg['selected']['P0r'],third=False)}
    diag=pickle.loads((args.cache/'diag.pkl').read_bytes());ood=[]
    for c in (5,6,7):ood+=pickle.loads((args.cache/f'ood{c}.pkl').read_bytes())
    cfg['offset_diagnostics']={'id':analyse(diag,models),'ood':analyse(ood,models)}
    cfg['offset_diagnostics_definition']='Local horizon error = (prediction-current measured eta) - (surface(t+12ms)+tau*surface_dot(t+12ms)-surface(t)). Oracle increment error = diff(prediction)-diff(physical target), excluding tracker resets. Segment-centered error removes one best-fit constant per reset segment for diagnosis only (uses future labels, never control or selection). None changes the original absolute target or validation criterion.'
    (RES/'step3b_predictors.json').write_text(json.dumps(cfg,indent=2,allow_nan=False),encoding='utf8');write_md(cfg)
    path=RES/'step3b_predictors.md';s=path.read_text(encoding='utf8')
    s+='\n## Offset-free diagnostics (post-selection only)\n\n'+cfg['offset_diagnostics_definition']+'\n\n| Split / family | Model | Local horizon RMSE μm | Oracle increment RMSE μm | Centered RMSE μm |\n|---|---|---:|---:|---:|\n'
    for split in ('id','ood'):
        for fam in ('all',*FAMS):
            for n in models:
                g=cfg['offset_diagnostics'][split][n][fam+'|all']
                s+=f"| {split} / {fam} | {n} | {g['local_horizon']['rmse_um']:.2f} | {g['oracle_increment']['rmse_um']:.2f} | {g['segment_centered']['rmse_um']:.2f} |\n"
    s+='\nThis separates cumulative measurement/contact bias from short-horizon extrapolation, but does not provide privileged bias correction to any predictor. A lower pooled absolute error need not imply better motion compensation. GP fitting reached its61-evaluation cap without convergence; its result is a bounded prototype, not a GP performance upper bound.\n'
    path.write_text(s,encoding='utf8');print('post-selection diagnostics completed')

if __name__=='__main__':main()
