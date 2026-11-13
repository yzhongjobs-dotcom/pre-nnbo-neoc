import os, sys, math
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[k]='1'
sys.dont_write_bytecode=True; sys.path.insert(0, os.getcwd())
import numpy as np
import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
P3,B=SE.P3,O.B; O.setup(); TC=.008
def ma(v,w=6):
    o=np.convolve(v,np.ones(w)/w)[:len(v)]; o[:w-1]=np.nan; return o
def hp(v,T):
    w=int(T/TC)|1; p=w//2; vv=np.pad(v,p,mode='edge'); return v-np.convolve(vv,np.ones(w)/w,mode='valid')
def r2(a,b):
    m=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[m],b[m])[0,1]**2) if m.sum()>50 and np.std(b[m])>0 else float('nan')
print('case group | theta range(deg) | R2 of ln k_eff on ln|M_true| whole, 2s-detrended | same with measured M | with -M_true/tan(theta) (truth pose) whole | measured moment noise: hf rms (mN m) vs |M| median | events: rms change of ln|M_meas| inside vs outside')
for idx in (36,39,38,0,3,12,24):
    case=CC.cases()[idx]; sc=SE.MC.scene_for(tuple(case['cell']),CC.BASE,case['variant'])
    pol=B.build('Huber-D4',O.chosen()['Huber-D4'],sc); CC.own_mechanism(pol)
    log=P3.episode(sc,pol,case['seed'])
    t=log[:,0]; ok=(t>1.5)&(log[:,26]>0)&(log[:,18]>.5)
    lk=np.log(np.maximum(ma(log[:,26]),1.)); Mt=ma(log[:,3]*log[:,24]); Mm=ma(log[:,13]); th=log[:,17]
    lMt=np.log(np.maximum(np.abs(Mt),1e-6)); lMm=np.log(np.maximum(np.abs(Mm),1e-6))
    prox=np.log(np.maximum(np.abs(Mt/np.tan(np.where(np.abs(th)>1e-4,th,np.nan))),1e-9))
    ev=np.zeros(len(t),bool)
    for e in sc.get('event_info',[]): ev|=(t>=e['start']-.05)&(t<=e['end']+.6)
    g=ok&~ev&np.isfinite(lk)&np.isfinite(lMt)
    hf=log[:,13]-ma(log[:,13],6)
    dM=np.abs(np.gradient(lMm,TC))
    print(f"{idx:2d} {case['group']:8s} | {math.degrees(np.nanmin(th[ok])):+5.1f}..{math.degrees(np.nanmax(th[ok])):+5.1f} | {r2(lk[g],lMt[g]):.2f} {r2(hp(np.nan_to_num(lk),2.)[g],hp(np.nan_to_num(lMt),2.)[g]):.2f} | {r2(lk[g],lMm[g]):.2f} {r2(hp(np.nan_to_num(lk),2.)[g],hp(np.nan_to_num(lMm),2.)[g]):.2f} | {r2(lk[g],prox[g]):.2f} | {1e3*np.nanstd(hf[ok]):.2f} vs {1e3*np.nanmedian(np.abs(Mm[ok])):.1f} | {np.sqrt(np.nanmean(dM[ok&ev]**2)) if (ok&ev).sum()>5 else float('nan'):.2f} / {np.sqrt(np.nanmean(dM[g]**2)):.2f}")
