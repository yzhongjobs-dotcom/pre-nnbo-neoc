import json
from pathlib import Path
import numpy as np

SOURCE = Path(r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study\results\huber_windows')
OUT = Path(__file__).parent/'results'/'force_term_replay.json'
G = {'A_large fast':[36,37,39,40], 'A_large slow':[38,41], 'A':list(range(6)), 'A_slow_thin':list(range(6,12))}
E = {'B':list(range(12,18)), 'A+B':list(range(18,24))}
def trail(a,n):
    s=np.r_[0.,np.cumsum(a)]; j=np.arange(len(a)); lo=np.maximum(0,j-n+1)
    return (s[j+1]-s[lo])/(j+1-lo)
def load(i):
    with np.load(SOURCE/f'case_{i:02d}.npz',allow_pickle=False) as f:
        log,raw,meta=f['log'],f['raw'],json.loads(str(f['meta']))
    t,k,fax=log[:,0],log[:,9],log[:,7]; lam=np.log(np.maximum(k,np.finfo(float).tiny)); l=trail(lam,6)
    dl=trail(np.gradient(l,.008),5); event=np.zeros(len(t),bool)
    for ev in meta['event_info']: event|=(t>=ev['start']-.05)&(t<=ev['end']+.6)
    score=(t>1.5)&(k>0)&(fax>.5); contact=(k>0)&(fax>.5)
    runs=np.cumsum(contact & ~np.r_[False,contact[:-1]])
    ks=np.full(len(t),np.nan); state=np.nan; alpha=1-np.exp(-.008/.3)
    for j,v in enumerate(raw):
        if np.isfinite(v) and v>0: state=v if not np.isfinite(state) else state+alpha*(v-state)
        ks[j]=state
    return dict(y=trail(log[:,2],6),f=trail(log[:,4],6),l=l,ks=ks,contact=contact,runs=runs,score=score,event=event,rise=(dl>1)&score&~event,fall=(dl< -1)&score&~event)
def variables(c,h):
    n=len(c['y']); j=np.arange(h,n); p=j-h; f=c['f']; l=c['l']; y=c['y']
    ok=c['contact'][j]&c['contact'][p]&(c['runs'][j]==c['runs'][p])&(f[j]>0)&(f[p]>0)
    T=np.full(n,np.nan); x=T.copy(); xm=T.copy(); z=T.copy(); valid=np.zeros(n,bool)
    T[j]=l[j]-l[p]; x[j]=-(y[j]-y[p])*np.exp(l[p])/f[p]; xm[j]=-(y[j]-y[p])*c['ks'][p]/f[p]; z[j]=np.log(f[j]/f[p]); valid[j]=ok
    return T,x,xm,z,valid
def corr(a,b,m):
    m=m&np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 and np.std(a[m])>0 and np.std(b[m])>0 else None
def fit(T,X,m):
    X=np.column_stack(X); m=m&np.isfinite(T)&np.all(np.isfinite(X),axis=1)
    a=np.linalg.lstsq(X[m],T[m],rcond=None)[0]; d=T[m]-X[m]@a; s=np.sum((T[m]-np.mean(T[m]))**2)
    return dict(coef=a.tolist(),r2=float(1-np.sum(d*d)/s),n=int(m.sum()))
def med(v):
    v=[x for x in v if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else None
def lag(x,y,m):
    dx=np.r_[np.full(5,np.nan),np.diff(x,prepend=np.nan)[5:]] if False else np.r_[np.full(5,np.nan),x[5:]-x[:-5]]
    dy=np.r_[np.full(5,np.nan),y[5:]-y[:-5]]; best=None
    for k in range(-25,26):
        a=np.arange(max(0,-k),min(len(x),len(x)-k)); b=a+k
        r=corr(dy[a],dx[b],m[a]&m[b])
        if r is not None and (best is None or r>best[1]): best=(k*8,r)
    return best

def main():
    cases={i:load(i) for i in range(42)}
    out={'definitions':{'fit':'No-intercept per-case OLS, outside events, centered-target R²; group medians.', 'contact':'Contiguous runs with positive true stiffness and fax_true > 0.5; both endpoints in same run.', 'smooth':'Positive finite raw, exponential linear-stiffness smooth, tau 0.3 s, first valid initialization, held through gaps.', 'timing':'5-tick differences; maximum Pearson correlation within ±200 ms, positive means predictor lags truth; outside-event valid ticks, per-case fits.', 'aggregation':'Group medians of per-case values; residuals are per-case segment means.'},'groups':{},'events':{}}
    for name,ids in G.items():
        out['groups'][name]={}
        for h in (9,24,48):
            rows=[]
            for i in ids:
                c=cases[i]; T,x,xm,z,v=variables(c,h); m=v&c['score']&~c['event']
                fs={k:fit(T,X,m) for k,X in {'x1':[x],'x2':[z],'x1_x2':[x,z],'x1m':[xm],'x1m_x2':[xm,z]}.items()}
                residual={}
                for key,X in [('x1',[x]),('x1_x2',[x,z])]:
                    pred=np.column_stack(X)@fs[key]['coef']; residual[key]={s:float(np.mean((T-pred)[m&c[s]])) if np.any(m&c[s]) else None for s in ('rise','fall')}
                a=fs['x1_x2']['coef']; combo=a[0]*x+a[1]*z
                timing={key:lag(xx,c['l'],m) for key,xx in [('x1',x),('x1_x2',combo),('x2',z)]}
                rows.append(dict(case=i,fits=fs,correlation=corr(x,z,m),correlation_measured=corr(xm,z,m),residual=residual,timing=timing))
            summary={'fits':{},'correlation':med([r['correlation'] for r in rows]),'correlation_measured':med([r['correlation_measured'] for r in rows]),'residual':{},'timing_ms':{},'cases':rows}
            for key in rows[0]['fits']:
                fs=[r['fits'][key] for r in rows]; summary['fits'][key]={'r2':med([f['r2'] for f in fs]),'coef':[med([f['coef'][q] for f in fs]) for q in range(len(fs[0]['coef']))]}
            for key in ('x1','x1_x2'): summary['residual'][key]={s:med([r['residual'][key][s] for r in rows]) for s in ('rise','fall')}
            for key in ('x1','x1_x2','x2'): summary['timing_ms'][key]=med([r['timing'][key][0] if r['timing'][key] else None for r in rows])
            out['groups'][name][str(h)]=summary
    coef=out['groups']['A_large fast']['9']['fits']['x1_x2']['coef']
    for name,ids in E.items():
        rows=[]
        for i in ids:
            c=cases[i]; T,x,xm,z,v=variables(c,9); m=v&c['score']&c['event']; combo=coef[0]*x+coef[1]*z
            rms=lambda a:float(np.sqrt(np.mean(a[m]**2))) if m.any() else None
            rows.append(dict(case=i,n=int(m.sum()),x1_rms=rms(x),x2_rms=rms(z),combo_rms=rms(combo),correlation=corr(x,z,m)))
        out['events'][name]={k:med([r[k] for r in rows]) for k in ('x1_rms','x2_rms','combo_rms','correlation')}; out['events'][name]['cases']=rows
    a=out['groups']['A_large fast']['9']; fs=a['fits']; r=a['residual']; t=a['timing_ms']
    out['predictions']={'F1':bool(fs['x1_x2']['r2']-fs['x1']['r2']>=.15 and .3<=fs['x1_x2']['coef'][1]<=2),'F2':bool((r['x1']['rise']<-.03 or r['x1']['fall']>.03) and abs(r['x1_x2']['rise'])<=.03 and abs(r['x1_x2']['fall'])<=.03),'F3':bool(t['x1']>=24 and abs(t['x1_x2'])<=8),'F4':all(out['events'][g]['combo_rms']<=out['events'][g]['x1_rms'] for g in E)}
    OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf8'); print(OUT); print(out['predictions'])
if __name__=='__main__': main()
