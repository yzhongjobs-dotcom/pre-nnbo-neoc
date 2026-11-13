"""Round 2. Four unchanged round-1 loops; logging-only instrumentation.
All outputs new, exclusively created. python -X utf8 lateral_dither_truth.py
"""
import os, sys, tempfile
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.environ['NUMBA_CACHE_DIR']=tempfile.mkdtemp(prefix='lateral_truth_jit_')
sys.path.insert(0,r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
import lateral_dither_check as R
import inspect, json, time, hashlib
import numpy as np
P=R.T.P3
PREDICTIONS=dict(T1='Best lag corr>=0.85 in 36,38,0',T2='Relief variance / total variance>=0.30 in 36,0',T3='Case39 rms(c_static*a_contact)<0.02 N',T4='All corr(dlnk_x,f_x/F)>=0.9 and gain 0.4..1.1')

def instrument():
    source=inspect.getsource(P._tick.py_func).replace('@njit(cache=True)','@njit(cache=False)')
    needle='log[nlog, 26] = abs(keff); log[nlog, 27] = math.atan2(keff.imag, keff.real)'
    assert source.count(needle)==1
    source=source.replace(needle,needle+'\n            log[nlog, 28] = xc; log[nlog, 29] = activation\n            log[nlog, 30] = fe; log[nlog, 31] = hc_factor\n            log[nlog, 32] = hereditary_force')
    exec(compile(source,'<logging_only_tick>','exec'),P.__dict__)
    P.LOGC=33

def rms(x):
    return float(np.sqrt(np.nanmean(np.asarray(x)**2)))

def analyze(path,out):
    meta,a=R.analyze(path)
    with np.load(path) as z: log=z['log']
    case=R.T.CC.cases()[meta['index']]
    sc=R.T.SE.MC.scene_for(tuple(case['cell']),R.T.CC.BASE,case['variant'])
    geom=P._geometry(sc)
    def force(x,d,th,vc,act):
        return P._contact(x,d,th,vc,geom,sc['q_field'],sc['E_field'],sc['ribs'],act)[0]
    ex,ed=.0002,1e-7
    n=len(log)
    fs,fx,ks,hx,dl,direct,kr=[np.full(n,np.nan) for _ in range(7)]
    for j,row in enumerate(log):
        x,d,th,vc,act=row[[28,6,17,16,29]]
        def f(xx,dd): return force(xx,dd,th,vc,act)
        def k(xx): return (f(xx,d+ed)-f(xx,d-ed))/(2*ed)
        def h(xx): return P._lookup(sc['h_s'],sc['dx'],xx-sc['field_x0'])
        fs[j]=f(x,d)
        fx[j]=(f(x+ex,d)-f(x-ex,d))/(2*ex)
        ks[j]=k(x)
        hx[j]=(h(x+ex)-h(x-ex))/(2*ex)
        kp,km=k(x+ex),k(x-ex)
        if kp>0 and km>0: dl[j]=(np.log(kp)-np.log(km))/(2*ex)
        direct[j]=(f(x+ex,d+h(x+ex)-h(x))-f(x-ex,d+h(x-ex)-h(x)))/(2*ex)
        kr[j]=P._contact(x,d,th,vc,geom,sc['q_field'],sc['E_field'],sc['ribs'],act)[1]
    relief=ks*hx; cs=fx+relief
    mask=a['score']&(log[:,3]>0)&(log[:,23]>0)
    def reg(x,y,m=mask,s=0): return R.regress(x,y,m,s)
    zero=reg(cs,a['c_hat'])
    best=max([reg(cs,a['c_hat'],mask,s) for s in range(-25,26)],key=lambda r:r['correlation'] if np.isfinite(r['correlation']) else -np.inf)
    def residual(r):
        ids=np.flatnonzero(mask); jj=ids+round(r['lag_s']/R.TC)
        ok=(jj>=0)&(jj<n); ids,jj=ids[ok],jj[ok]
        ok=mask[jj]&np.isfinite(cs[ids])&np.isfinite(a['c_hat'][jj]); ids,jj=ids[ok],jj[ok]
        return dict(rms=rms(a['c_hat'][jj]-r['gain']*cs[ids]),predicted_sd_median=R.med(a['cs_predicted'][jj]),predicted_sd_rms=rms(a['cs_predicted'][jj]))
    def shares(x,y,m=mask):
        ok=m&np.isfinite(x)&np.isfinite(y); x,y=x[ok],y[ok]; total=np.var(x+y)
        return dict(first=float(np.var(x)/total),second=float(np.var(y)/total),cross=float(2*np.mean((x-x.mean())*(y-y.mean()))/total),correlation=float(np.corrcoef(x,y)[0,1]))
    spatial=dl*a['v24']; other=a['rate_true']-spatial
    amp=a['contact_harmonic_amplitudes'][:,0]
    # Column21 is evaluated at 20.8 Hz. Reconstruct its 15.625Hz analogue.
    const=P._constants(sc,P.PLATFORM); alpha,tau,_=P.memory_arrays(sc)
    w=2*np.pi*R.FX; k15=np.full(n,np.nan); k20=np.full(n,np.nan)
    for j,row in enumerate(log):
        th,vc=row[[17,16]]; ct,st=np.cos(th),np.sin(th)
        tangent=kr[j]*row[31]
        if sc['contact_law']=='halfspace' and row[3]>0:
            tangent=(force(row[28],row[6]+1e-6,th,vc,row[29])-force(row[28],row[6]-1e-6,th,vc,row[29]))/2e-6*row[31]
        friction=P._friction(1.,vc,row[23],const,geom)-const[17]*vc
        axial=ct+friction*st; lateral=friction*ct-st
        series=ct*(axial/sc['kr']+const[26]*lateral)
        for freq,dest in ((R.FX,k15),(20.8,k20)):
            w=2*np.pi*freq
            kc=tangent*(1+np.sum(alpha*(1j*w*tau)/(1+1j*w*tau)))+1j*w*const[4]*max(row[32],0)
            dest[j]=abs(kc/(1+kc*series)*axial*ct)
    pick=int(np.flatnonzero(mask)[len(np.flatnonzero(mask))//2])
    row=dict(index=meta['index'],n_scored=int(mask.sum()),measurement=dict(zero=zero,best=best,residual_zero=residual(zero),residual_best=residual(best),K20p8_over_ks=R.med((log[:,21]/ks)[mask]),K15p625_over_ks_reconstructed=R.med((k15/ks)[mask])),physics=dict(variance=shares(fx,relief),fx_over_F=reg(dl,fx/a['F24']),fx_over_static_F=reg(dl,fx/fs),c_over_F=reg(a['G_true'],cs/a['F24']),spatial_rate=reg(a['rate_true'],spatial),rate_variance=shares(spatial,other),spatial_rms=rms(spatial[mask]),other_rms=rms(other[mask]),total_rate_rms=rms(a['rate_true'][mask])),signal=dict(amplitude_mm_median=1e3*R.med(amp[mask]),rms_N=rms((cs*amp)[mask]),noise_N=R.SIGMA),checks=dict(force_max_error=float(np.max(abs(fs-log[:,30]))),theta_max=float(np.max(abs(log[:,17]))),direct_tick=pick,direct_t=float(log[pick,0]),c_static=float(cs[pick]),direct=float(direct[pick]),direct_rms_error=rms((direct-cs)[mask]),K20_reconstruction_max_relative_error=float(np.nanmax(abs(k20/log[:,21]-1))),k_derivative_vs_returned_relative_median=R.med((ks/kr-1)[mask]),contact_law=sc['contact_law']))
    a.update(log=log,f_static=fs,f_x=fx,k_s=ks,h_x=hx,c_static=cs,dlnk_x=dl,c_direct=direct,spatial_rate=spatial,other_rate=other,score_contact=mask,K15p625=k15)
    with (out/f'case_{meta["index"]:02d}_truth.npz').open('xb') as f: np.savez_compressed(f,**a)
    print(json.dumps(R.T.CC.clean(row)),flush=True)
    return row

def main():
    out=Path(__file__).resolve().parent/'results'/time.strftime('lateral_truth_%Y%m%d_%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    R.dump_new(out/'predictions_before_episodes.json',dict(predictions=PREDICTIONS,unix_time=time.time(),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    instrument(); rows=[]
    for i in (36,39,38,0):
        name=R.run((i,.1,str(out)))
        rows.append(analyze(out/name,out))
    get={r['index']:r for r in rows}
    verdicts=dict(T1=all(get[i]['measurement']['best']['correlation']>=.85 for i in (36,38,0)),T2=all(get[i]['physics']['variance']['second']>=.3 for i in (36,0)),T3=get[39]['signal']['rms_N']<.02,T4=all(r['physics']['fx_over_F']['correlation']>=.9 and .4<=r['physics']['fx_over_F']['gain']<=1.1 for r in rows))
    R.dump_new(out/'summary.json',dict(predictions=PREDICTIONS,verdicts=verdicts,runs=rows,conventions=dict(derivatives='dx=0.2mm; dd=0.1um; fixed logged theta,vc,activation',columns='0:28 unchanged;28 xc;29 activation;30 elastic force;31 HC factor;32 hereditary force',height='positive y and h increase depth; requested c_static=fx+ks*hx; theta checked',F='round1 F24; static F alternative also saved',variance='Var(first)/Var(sum),Var(second)/Var(sum),2Cov/Var(sum); sum=1',regression='OLS with intercept; x=truth,y=signal; residual requested omits intercept; positive lag means c_hat delayed',K='col21=20.8Hz, not15.625; 15.625 reconstructed with plant formula')))
    print('OUTPUT',out,flush=True)
if __name__=='__main__': main()