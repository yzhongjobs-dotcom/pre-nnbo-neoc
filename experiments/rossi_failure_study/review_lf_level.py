"""Verification/feasibility on the 42 stored closed loops (no new simulation).
1) physics: is ln k_eff a function of the indentation change at constant force?  psi = -int (k/F)(dy + s du)
2) replay: probe window output m = ln raw; covariate estimate = EMA[m - b*psibar] + b*psi, b by causal recursive LS.
   psi from truth (y and u: mechanism ceiling) or from the encoder alone, high-passed (naive deployable form)."""
import numpy as np, json, math, sys
TC=.008
def ma(v,w=6):
    out=np.convolve(v,np.ones(w)/w)[:len(v)]; out[:w-1]=np.nan; return out
def ema(v,tau):
    a=1-math.exp(-TC/tau); o=np.full(len(v),np.nan); s=None
    for i,x in enumerate(v):
        if np.isfinite(x): s=x if s is None else s+a*(x-s)
        if s is not None: o[i]=s
    return o
def hp(v,tau): return v-ema(v,tau)
def cov_est(m,psi,tau_a=1.5,tau_h=2.,lam=.9995):
    """a(t)+b(t)*psi(t); b from recursive LS of high-passed m on high-passed window-averaged psi (forgetting lam per tick)."""
    pb=ma(psi,24); hm=hp(m,tau_h); hpb=hp(pb,tau_h)
    sxx=1e-6; sxy=0.; out=np.full(len(m),np.nan); bs=np.full(len(m),np.nan); a=None; al=1-math.exp(-TC/tau_a)
    for i in range(len(m)):
        if not (np.isfinite(m[i]) and np.isfinite(pb[i]) and np.isfinite(hm[i]) and np.isfinite(hpb[i])): continue
        b=min(max(sxy/sxx,0.),3.) if sxx>1e-3 else 0.   # coefficient from data strictly before tick i
        r=m[i]-b*pb[i]; a=r if a is None else a+al*(r-a)
        out[i]=a+b*psi[i]; bs[i]=b
        sxx=lam*sxx+hpb[i]**2; sxy=lam*sxy+hpb[i]*hm[i]
    return out,bs
groups={}
for i in range(42):
    d=np.load(f'results/huber_windows/case_{i:02d}.npz'); meta=json.loads(str(d['meta'])); log=d['log']
    t=log[:,0]; y=ma(log[:,2]); u=ma(log[:,6]); F=np.maximum(ma(log[:,7]),1.); k=np.maximum(log[:,9],1.); raw=d['raw']
    lnk=np.log(k); ok=(t>1.5)&(log[:,9]>0)&(log[:,7]>.5)&np.isfinite(raw)&(raw>0)
    m=np.where(np.isfinite(raw)&(raw>0),np.log(np.clip(raw,200,3e5)),np.nan)
    kc=np.exp(ema(m,.5))                       # causal slow scale for the conversion (no truth)
    dy=np.r_[0,np.diff(np.nan_to_num(y))]; du=np.r_[0,np.diff(np.nan_to_num(u))]
    # sign of u: tracking means dy ~ -s du
    g=np.isfinite(y)&np.isfinite(u)&(t>1.5); s=-np.sign(np.sum(dy[g]*du[g])) or 1.
    psi_true=-np.cumsum(np.nan_to_num(kc/8.)*(dy+s*du))      # indentation change in units of F/k (F set to 8 N nominal)
    psi_enc=hp(-np.cumsum(np.nan_to_num(kc/8.)*dy),1.0)     # encoder only, breathing-band removed causally
    dl=ma(np.gradient(ma(lnk),TC),5)
    ev=np.zeros(len(t),bool)
    for e in meta['event_info']: ev|=(t>=e['start']-.05)&(t<=e['end']+.6)
    # 1) physics fit on 2 s high-passed quantities (non-causal diagnostic)
    gg=ok&np.isfinite(psi_true)&~ev
    def fit(psi):
        a=hp(lnk,2.)[gg]; b=hp(psi,2.)[gg]
        if len(a)<50 or b.std()<1e-9: return np.nan,np.nan
        sl=np.polyfit(b,a,1)[0]; return sl,np.corrcoef(a,b)[0,1]**2
    sT,rT=fit(psi_true); sE,rE=fit(psi_enc)
    # 2) replay
    eT,bT=cov_est(m,psi_true); eE,bE=cov_est(m,psi_enc); eM=ema(m,1.5)
    err=lambda est: est-lnk
    seg={'plateau':ok&~ev&(np.abs(dl)<.2),'rise':ok&~ev&(dl>1.),'fall':ok&~ev&(dl<-1.),'slow':ok&~ev&(np.abs(dl)>.4)&(np.abs(dl)<=1.),'event':ok&ev}
    G=groups.setdefault(meta['group'],dict(sT=[],rT=[],sE=[],rE=[],seg={sname:{n:[] for n in ('raw','mem','covT','covE')} for sname in seg}))
    G['sT'].append(sT);G['rT'].append(rT);G['sE'].append(sE);G['rE'].append(rE)
    for sname,mask in seg.items():
        for n,est in (('raw',m),('mem',eM),('covT',eT),('covE',eE)):
            e=err(est)[mask&np.isfinite(est)]; G['seg'][sname][n].append(e)
print('group | physics: slope/R2 with true motion | encoder only || segment: n, then for raw / memory 1.5 s / covariate(true u) / covariate(encoder): rms error and P05 (most negative 5%)')
for gname,G in groups.items():
    print(f"{gname:12s} | {np.nanmedian(G['sT']):.2f}/{np.nanmedian(G['rT']):.2f} | {np.nanmedian(G['sE']):.2f}/{np.nanmedian(G['rE']):.2f}")
    for sname,dd in G['seg'].items():
        n=sum(len(x) for x in dd['raw'])
        if n<40: continue
        line=f"     {sname:8s} n={n:5d} "
        for nme in ('raw','mem','covT','covE'):
            e=np.concatenate(dd[nme]) if len(dd[nme]) else np.array([np.nan])
            line+=f"| {nme} {np.sqrt(np.mean(e**2)):.3f} {np.percentile(e,5):+.2f} "
        print(line)
