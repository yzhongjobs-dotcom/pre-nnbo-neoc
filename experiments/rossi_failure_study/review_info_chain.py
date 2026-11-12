"""Verification on the stored closed loops (results/huber_windows): where the probe information goes between the force ripple and the window estimate. Run: python -X utf8 review_info_chain.py 0 3 6 12 15 18 24 30 36 38"""
import numpy as np, json, math, sys
TC=.008; F0=20.8; N=24
def Mmat(m, fc, order=4):
    rho=(2*math.sin(math.pi*fc*TC))**(-2*order); D=np.diff(np.eye(m),order,axis=0)
    return np.eye(m)-np.linalg.inv(np.eye(m)+rho*D.T@D)
def proj_perp(B):
    Q,_=np.linalg.qr(B); return np.eye(len(B))-Q@Q.T
M14=Mmat(N,14.); M35=Mmat(N,3.5); M7=Mmat(N,7.)
tt=np.arange(N)*TC
P3=proj_perp(np.vander(np.linspace(-1,1,N),4))
P0=proj_perp(np.ones((N,1)))
SC=np.c_[np.sin(2*np.pi*F0*tt),np.cos(2*np.pi*F0*tt)]
def band(v):  # energy of v (after cubic removal) inside span{sin,cos} after cubic removal
    B=P3@SC; Q,_=np.linalg.qr(B); c=Q.T@(P3@v); return float(c@c)
rows=[]
for i in [int(a) for a in sys.argv[1:]]:
    d=np.load(f'results/huber_windows/case_{i:02d}.npz'); meta=json.loads(str(d['meta'])); log=d['log']
    Y=d['Y'][:,-N:]; z=d['z'][:,-N:]; Y0=d['Y0'][:,-N:]; z0=d['z0'][:,-N:]
    ok=np.isfinite(Y).all(1)&np.isfinite(Y0).all(1)&(log[:,0]>1.)&(log[:,9]>0)
    for _,a0,a1 in []: pass
    k=log[:,9]
    eps=(z-z0)[ok]; sig=eps[:,-1].std()   # newest element: one fresh sample per window
    lag1=np.corrcoef(eps[:,-1],eps[:,-2])[0,1]
    out=dict(case=i,group=meta['group'],sig=sig,lag1=lag1,k_med=np.median(k[ok]))
    idx=np.flatnonzero(ok)[::6]
    acc={n:[] for n in ('none','mean','poly3','D4_3.5','D4_7','D4_14','band','probeF_rms','allF_hp_rms')}
    kn=[];k0=[]
    for n in idx:
        y=Y0[n]; kk=k[n]
        acc['none'].append(kk*kk*(y@y)/sig**2/N); acc['mean'].append(kk*kk*(y@P0@y)/sig**2/N)
        acc['poly3'].append(kk*kk*(y@P3@y)/sig**2/N)
        acc['D4_3.5'].append(kk*kk*(y@M35@y)/sig**2/N); acc['D4_7'].append(kk*kk*(y@M7@y)/sig**2/N); acc['D4_14'].append(kk*kk*(y@M14@y)/sig**2/N)
        acc['band'].append(kk*kk*band(y)/sig**2/N)
        acc['probeF_rms'].append(kk*math.sqrt(band(y)/N)); acc['allF_hp_rms'].append(kk*math.sqrt((y@P3@y)/N))
        my=M14@Y[n]; kn.append((my@z[n])/(Y[n]@my)); my0=M14@Y0[n]; k0.append((my0@z0[n])/(Y0[n]@my0))
    for a,v in acc.items(): out[a]=float(np.median(v))
    kn=np.array(kn);k0=np.array(k0); e=(kn-k0)/k[idx]
    out['noise_sd_rel_actual']=float(1.4826*np.median(np.abs(e-np.median(e))))
    out['noise_sd_rel_pred']=float(np.median(1/np.sqrt(N*np.array(acc['D4_14']))))
    # whole-episode 20.8 Hz ripple in true axial force and in measured force
    t=log[ok,0]; 
    for name,c in (('fax',7),('Fmeas',4)):
        f=log[ok,c]; w=12; hp=f-np.convolve(f,np.ones(w)/w,mode='same'); hp=hp[w:-w]; tw=t[w:-w]
        A=np.c_[np.cos(2*np.pi*F0*tw),np.sin(2*np.pi*F0*tw)]; coef=np.linalg.lstsq(A,hp,rcond=None)[0]
        out['ripple_'+name]=float(np.hypot(*coef)/math.sqrt(2)); out['hp_rms_'+name]=float(hp.std())
    rows.append(out)
keys=[k for k in rows[0] if k not in('case','group')]
print('case group '+' '.join(keys))
for r in rows: print(r['case'],r['group'],' '.join(f"{r[k]:.4g}" for k in keys))
