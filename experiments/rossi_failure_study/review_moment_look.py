"""Feasibility look only: does the moment channel the plant already gives the Method mark where stiffness changes along the scan?
Prediction written before running: for a smooth gradient CoP ~ (L^2/12) dln(kappa)/dx, so dln(k_eff)/dt ~ 12 v CoP / L^2 on rib crossings;
body-motion events (push/pull, cough) should move the CoP far less than rib crossings do."""
import os, sys, json, math
for _k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[_k]='1'
sys.dont_write_bytecode=True; sys.path.insert(0, os.getcwd())
import numpy as np
import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
P3,B=SE.P3,O.B
TC=.008
def ma(v,w=6):
    out=np.convolve(v,np.ones(w)/w)[:len(v)]; out[:w-1]=np.nan; return out
cases={c['index']:c for c in CC.cases()}
O.setup()
print('case group | L(m) v(m/s) | corr(true CoP, dlnk/dt) lead-lag best(s) | corr(measured M/F, dlnk/dt) | slope check 12 v CoP/L^2 vs dlnk (fit gain) | CoP rms mm: plateau / change / event | measured M/F rms mm same')
for i in [int(a) for a in sys.argv[1:]]:
    case=cases[i]; sc=SE.MC.scene_for(tuple(case['cell']),CC.BASE,case['variant'])
    pol=B.build('Huber-D4',O.chosen()['Huber-D4'],sc)
    log=P3.episode(sc,pol,case['seed'],diagnostic=True) if 'diagnostic' in P3.episode.__code__.co_varnames else P3.episode(sc,pol,case['seed'])
    t=log[:,0]; ok=(t>1.5)&(log[:,26]>0)&(log[:,18]>.5)
    k=np.maximum(log[:,26],1.); dl=ma(np.gradient(ma(np.log(k)),TC),5)
    cop=ma(log[:,24]); F=np.maximum(ma(log[:,5]),.5); M=ma(log[:,13]); ft=ma(log[:,12]); L=np.nanmedian(log[ok,23]); v=np.gradient(ma(log[:,1]),TC)
    mcop=M/F
    ev=np.zeros(len(t),bool)
    for e in sc.get('event_info',[]): ev|=(t>=e['start']-.05)&(t<=e['end']+.5)
    g=ok&~ev&np.isfinite(dl)&np.isfinite(cop)&np.isfinite(mcop)
    def lagcorr(a,b):
        best=(0,0.)
        for s in range(-40,41,2):
            aa=np.roll(a,s); m=g.copy(); m[:45]=False; m[-45:]=False
            c=np.corrcoef(aa[m],b[m])[0,1]
            if abs(c)>abs(best[1]): best=(s*TC,c)
        return best
    c0=np.corrcoef(cop[g],dl[g])[0,1]; lg=lagcorr(cop,dl); c1=np.corrcoef(mcop[g],dl[g])[0,1]
    pred=12*np.nanmedian(np.abs(v[g]))*cop/L**2; gain=np.polyfit(dl[g],pred[g],1)[0]
    pl=g&(np.abs(dl)<.2); ch=g&(np.abs(dl)>1.); e=ok&ev&np.isfinite(cop)
    r=lambda x,m: 1e3*float(np.sqrt(np.mean((x[m]-np.nanmedian(x[g]))**2))) if m.sum()>3 else float('nan')
    print(f"{i:2d} {case['group']:10s} | {L:.4f} {np.nanmedian(np.abs(v[g])):.4f} | {c0:+.2f}  best {lg[1]:+.2f} at {lg[0]:+.2f}s | {c1:+.2f} | {gain:+.2f} | {r(cop,pl):.2f} / {r(cop,ch):.2f} / {r(cop,e):.2f} | {r(mcop,pl):.2f} / {r(mcop,ch):.2f} / {r(mcop,e):.2f}",flush=True)
