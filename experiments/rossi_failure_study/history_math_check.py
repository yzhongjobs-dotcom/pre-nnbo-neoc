"""Algebra and bias accounting on frozen original physical-plant observations.

The logged k_eff is used ONLY for offline attribution. A k_eff path does not
identify true body motion; residual below combines d, dynamic mismatch, and
errors of this reference contact representation. No replacement plant.
"""
import json
import numpy as np
import adaptive_history_v2 as A


def quantities(d,j,n,T=None,yT=None,kpath=None):
    log=d['log'];T=log[j,0] if T is None else T;yT=log[j,2] if yT is None else yT
    p=d['meta']['policy'];ns=p['ns'];eb=d['eb'][j];afc=p['afc']
    a=1-afc;ebn=eb**ns
    shape=(n,ns);tt=np.zeros(shape);pp=np.zeros(shape);ww=np.zeros(shape)
    for row,tick in enumerate(range(j-n,j)):
        for m in range(ns):
            w=afc*a**(ns-1-m);fr=(1-eb**m)/(1-ebn) if ebn<1 else m/ns
            tt[row,m]=log[tick,0]+m*.001
            pp[row,m]=log[tick,2]+fr*(log[tick+1,2]-log[tick,2])-yT;ww[row,m]=w
    h0=(ww*pp).sum(axis=1);h1=(ww*pp*(tt-T)).sum(axis=1)
    c0=ww.sum(axis=1);c1=(ww*tt).sum(axis=1)
    out=dict(h0=h0,h1=h1,c0=c0,c1=c1)
    if kpath is not None:
        kval=np.interp(tt.ravel(),d['truth'][:,0],kpath).reshape(shape)
        out['Ak']=(ww*pp*kval).sum(axis=1)
    return out


def main():
    rows=[];anchor=0.;moment=0.;decomp=0.;mc=0.;audit_done=False;fit_transport=0.
    for index in (0,24,36,37):
        f=np.load(A.OUT/f'fits_{index:02d}.npz');stream=np.load(A.OUT/f'stream_{index:02d}.npz')
        d={k:stream[k] for k in stream.files};d['meta']=json.loads(str(d['meta']))
        true=d['truth'];k=true[:,26];t=true[:,0];v=np.gradient(k,t)
        seg=A.V.segments(t,k,d['meta']['event_info'],k>200);ids=f['ids'];theta,Q,valid=(f[key] for key in ('theta','Q','valid'))
        components={name:np.full((len(ids),4),np.nan) for name in ('linear24','linear96','constant24','constant96')}
        for u,j in enumerate(ids):
            for n,slot,name in ((24,0,'linear24'),(96,4,'linear96'),(24,5,'constant24'),(96,9,'constant96')):
                if not valid[u,slot]:continue
                g=quantities(d,j,n,kpath=k);h0,h1,z=A.window(d,j,n);z0=A.window(d,j,n,zero=True)[2]
                assert np.allclose(h0,g['h0'],rtol=1e-10,atol=1e-15)
                assert np.allclose(h1,g['h1'],rtol=1e-10,atol=1e-15)
                qq=Q[u,slot,0,-n:]
                model=k[j]*h0+v[j]*h1
                parts=np.array([qq@(z-z0),qq@(g['Ak']-model),qq@(z0-g['Ak']),v[j]*(qq@h1)])
                if slot<5:
                    moment=max(moment,abs(qq@h0-1),abs(qq@h1))
                else:moment=max(moment,abs(qq@h0-1))
                error=theta[u,slot,0]-k[j]
                decomp=max(decomp,abs(parts.sum()-error)/max(k[j],1))
                components[name][u]=parts/max(k[j],1)
                if index==36 and not audit_done and n==96 and slot==4:
                    audit_done=True
                    T=t[j];delta=.017;dy=.00023
                    shifted=quantities(d,j,n,T+delta,d['log'][j,2]+dy)
                    r0=shifted['h0']-(g['h0']-dy*g['c0'])
                    r1=shifted['h1']-(g['h1']-delta*g['h0']-dy*(g['c1']-(T+delta)*g['c0']))
                    anchor=max(abs(r0).max(),abs(r1).max())
                    eps=np.random.default_rng(7).normal(size=(30000,96))*A.SIG
                    D=Q[u,4]-Q[u,0];pred=A.SIG**2*(D@D.T)
                    empirical=np.cov((eps@D.T).T)
                    mc=np.max(abs(empirical-pred)/np.sqrt(np.outer(np.diag(pred),np.diag(pred))))
                    fitter=A.Fitter();old=fitter.fit(h0,h1,z);new=fitter.fit(shifted['h0'],shifted['h1'],z)
                    target=np.array([old['theta'][0]+delta*old['theta'][1]/A.TREF,old['theta'][1]])
                    fit_transport=float(np.max(abs(new['theta']-target))/max(abs(target[0]),1))
        result=dict(index=index,segments={})
        for name,mask in seg.items():
            use=mask[ids]
            if use.sum()<4:continue
            result['segments'][name]={nm:dict(force_noise_rms=float(np.nanmean(c[use,0]**2)**.5),
                curvature_rms=float(np.nanmean(c[use,1]**2)**.5),
                remaining_dynamics_rms=float(np.nanmean(c[use,2]**2)**.5),
                linear_trend_omission_rms=float(np.nanmean(c[use,3]**2)**.5)) for nm,c in components.items()}
        rows.append(result)
    assert audit_done
    # Independent polynomial benchmark for the noise/curvature calculation.
    tt=np.linspace(-.4,0.,1001);X=np.column_stack([np.ones(len(tt)),tt])
    a=3.;curve=.5*a*tt**2;coef=np.linalg.lstsq(X,curve,rcond=None)[0]
    expected=np.array([-a*.4**2/12,-a*.4/2])
    quadratic_error=float(np.max(abs((coef-expected)/expected)))
    result=dict(checks=dict(uniform_quadratic_relative_discretization_error=quadratic_error,audit_executed=audit_done,fit_transport_relative=fit_transport,anchor_absolute=anchor,moment_max=moment,relative_decomposition_max=decomp,covariance_mc_normalized_error=mc),
                attribution='relative to logged k_eff path; remaining dynamics is NOT a uniquely identified body-motion component',cases=rows)
    (A.OUT/'math_check.json').write_text(json.dumps(result,indent=1))
    print(json.dumps(result,indent=1))
if __name__=='__main__':main()
