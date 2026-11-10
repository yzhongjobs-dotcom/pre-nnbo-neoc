"""Conditional slope information and equal-physical-force design; no controller changes.

Reproduce: python experiments/rossi_failure_study/v2_information_design.py
Regenerates only four existing fixtures when absent. Truth is used only for
labels/offline normalization; no proposed online rule consumes it. The design
calculation uses the frozen H operator and the exact original observation
sampling; its encoder trajectories are NOT executable robot commands.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
import json,math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
import huber_prior_check as H
import v2_verify as V

HERE=Path(__file__).resolve().parent
OUT=HERE/'results/v2_information_design.json'
N=24; TC=.008; NS=8; SIG=H.SIGMA
D=np.diff(np.eye(N),4,axis=0)
M=np.eye(N)-np.linalg.solve(np.eye(N)+H.RHO*D.T@D,np.eye(N))
T=np.arange(N)*TC
P=np.vander(np.linspace(-1,1,N),4)
Qp=np.linalg.qr(P)[0]
R=np.eye(N)-Qp@Qp.T  # optimistic: all non-polynomial disturbances known/absent

def weights(h0,h1,G):
    a=float(h0@G@h0);b=float(h0@G@h1)
    g=h1-b/a*h0;J=float(g@G@g)
    lv=G@g/J;q0=G@h0/a;qf=q0-(b/a)*lv
    return dict(a=a,b=b,J=J,lv=lv,q0=q0,qf=qf)

def quantile(x):
    x=np.asarray(x);x=x[np.isfinite(x)]
    return dict(median=float(np.median(x)),p10=float(np.percentile(x,10)),p90=float(np.percentile(x,90)),n=len(x)) if len(x) else None

def recorded_case(idx):
    p=H.WIN/f'case_{idx:02d}.npz'
    if not p.exists():
        import huber_windows as W,corrected_comparison as C
        H.WIN.mkdir(exist_ok=True);W.job(C.cases()[idx])
    d=H.load(p);log,_=H.masks_for(d);k=log[:,26];t=log[:,0]
    seg=V.segments(t,k,d['meta']['event_info'],k>200)
    # Covariance from the timeline's new observation, not from repeated overlapping rows.
    nz=d['z'][:,-1]-d['z0'][:,-1];ok=np.isfinite(nz)
    zz=nz[ok];zz=zz-np.mean(zz)
    ac=np.array([np.mean(zz[:len(zz)-lag]*zz[lag:]) for lag in range(N)])
    cov=ac[np.abs(np.arange(N)[:,None]-np.arange(N)[None,:])]
    ev,evec=np.linalg.eigh(cov);floor=max(ac[0]*1e-5,1e-14)
    cov=(evec*np.maximum(ev,floor))@evec.T
    inv=np.linalg.solve(cov,np.eye(N))
    Rg=inv-inv@P@np.linalg.solve(P.T@inv@P,P.T@inv)
    ew,uw=np.linalg.eigh(cov);white=(uw/np.sqrt(ew))@uw.T
    qp=np.linalg.qr(white@P)[0];proj=np.eye(N)-qp@qp.T
    tlocal=(np.arange(N)-(N-1))*TC
    carrier=np.c_[np.sin(2*np.pi*20.8*tlocal),np.cos(2*np.pi*20.8*tlocal)]
    def match_matrix(B):
        U,sv,_=np.linalg.svd(proj@white@B,full_matrices=False)
        U=U[:,sv>sv[0]*1e-10]
        return white.T@U@U.T@white
    R2=match_matrix(carrier)
    R4=match_matrix(np.column_stack([carrier,tlocal[:,None]*carrier]))
    # Existing probe-band implementation: bins 3--5 after the D4 nuisance space.
    band=np.column_stack([f(2*np.pi*b*np.arange(N)/N) for b in [3,4,5] for f in [np.cos,np.sin]])
    band=band-Qp@(Qp.T@band);qb=np.linalg.qr(band)[0];PB=qb@qb.T
    result={s:[] for s in ['all','plateau','rising','falling','reversal','event']}
    for j in range(N+1,len(k),6):
        h0=d['Y'][j,-N:];h1=d['H1'][j,-N:]
        if not np.isfinite(np.r_[h0,h1]).all() or k[j]<=200:continue
        _,_,w0,_=V.solve_v2_window(h0,h1,d['z'][j,-N:],rho=H.RHO,sd=.04/math.sqrt(H.RHO),sigma_z=SIG,sv2=0.)
        Mw=V.m_of(w0)
        v=weights(h0,h1,Mw);g=weights(h0,h1,Rg)
        b2=weights(h0,h1,R2);b4=weights(h0,h1,R4)
        Wiv=np.vstack([PB@M@h0,PB@M@h1]);Aiv=Wiv@np.column_stack([h0,h1]);qiv=np.linalg.solve(Aiv,Wiv)[1]
        I00=g['a'];I01=g['b'];I11=float(h1@Rg@h1)
        rho=I01/math.sqrt(I00*I11)
        sv=math.sqrt(float(v['lv']@cov@v['lv']))/k[j]
        opt=math.sqrt(float(g['lv']@cov@g['lv']))/k[j]
        s0=math.sqrt(float(v['q0']@cov@v['q0']))/k[j]
        sf=math.sqrt(float(v['qf']@cov@v['qf']))/k[j]
        noise_y=d['Y'][j,-N:]-d['Y0'][j,-N:]
        vtrue=float(np.polyfit(t[j-N:j]-t[j-1],.5*(k[j-N:j]+k[j-N+1:j+1]),1)[0])
        zz0=d['z0'][j,-N:]
        row=dict(existing_d4_probe_band_slope_sd=math.sqrt(float(qiv@cov@qiv))/k[j],huber_active_fraction=float(np.mean(w0<1.)),reference_rate=vtrue/k[j],
                 slope_snr=abs(vtrue/k[j])/max(sv,1e-15),
                 matched_slope_snr=abs(vtrue/k[j])/(math.sqrt(float(b4['lv']@cov@b4['lv']))/k[j]),
                 information_needed_snr1=(sv/max(abs(vtrue/k[j]),1e-12))**2,
                 template2_slope_sd=math.sqrt(float(b2['lv']@cov@b2['lv']))/k[j],
                 template4_slope_sd=math.sqrt(float(b4['lv']@cov@b4['lv']))/k[j],
                 template4_noisefree_slope_error=float((b4['lv']@zz0-vtrue)/k[j]),
                 effective_time_lever=-v['b']/v['a'],v2_slope_sd=sv,efficient_slope_sd=opt,efficiency_information_gain=(sv/opt)**2,
                 v2_noise_free_free_k_error=float(v['qf']@zz0/k[j]-1),
                 gls_noise_free_free_k_error=float(g['qf']@zz0/k[j]-1),
                 v2_noise_free_slope_error=float((v['lv']@zz0-vtrue)/k[j]),
                 gls_noise_free_slope_error=float((g['lv']@zz0-vtrue)/k[j]),
                 free_k_noise_inflation=sf/s0,gls_correlation=rho,gls_free_k_inflation=1/math.sqrt(max(1-rho*rho,1e-15)),
                 sigma_z_proxy_ratio=sv/(SIG*np.linalg.norm(v['lv'])/k[j]),
                 noisy_regressor_to_force_noise=float(np.std(k[j]*noise_y))/max(np.std(zz),1e-15))
        result['all'].append(row)
        for s in seg:
            if seg[s][j]:result[s].append(row)
    return dict(index=idx,group=d['meta']['group'],noise_sd=float(np.std(zz)),noise_ac1=float(ac[1]/ac[0]),
                covariance_note='Toeplitz from unique new observations; eigenvalue floor 1e-5 of variance; descriptive single fixture',
                segments={s:{a:quantile([r[a] for r in rr]) for a in rr[0]} for s,rr in result.items() if rr})

def basis(phase,phi=.4,tau_a=.02):
    w0=2*np.pi*20.8;dd=(1-np.exp(-1j*w0*TC))/TC
    th=np.tan(phi)/(dd.imag-np.tan(phi)*dd.real);mh=abs(1+th*dd);a=th/TC
    tt=np.arange(-NS,N*NS+1)*.001;xi=tt/(N*TC)
    win=np.where((xi>=0)&(xi<=1),64*xi**3*(1-xi)**3,0.)
    B=np.column_stack([np.ones(len(tt))]+[win*np.polynomial.legendre.legval(2*xi-1,[0]*i+[1]) for i in range(4)])
    Y=np.sin(w0*tt+phase)[:,None]*B
    # Same H lag, on the fine grid; only current/future force enters cost.
    F=((1+a)*Y[NS:]-a*Y[:-NS])/mh
    coarse=Y[NS::NS];previous=Y[0:1]
    yh=((1+a)*coarse-a*np.vstack([previous,coarse[:-1]]))/mh
    frac=(1-np.exp(-np.arange(NS)*.001/tau_a))/(1-np.exp(-TC/tau_a))
    af=1-np.exp(-2*np.pi*30*.001);ww=af*(1-af)**np.arange(NS-1,-1,-1)
    pp=yh[:-1,None,:]+frac[None,:,None]*(yh[1:]-yh[:-1])[:,None,:]-yh[-1]
    h0=np.sum(ww[None,:,None]*pp,axis=1)
    st=(np.arange(N)[:,None]*TC+np.arange(NS)[None,:]*.001-N*TC)
    h1=np.sum(ww[None,:,None]*st[:,:,None]*pp,axis=1)
    env=B[NS:]
    return h0,h1,F[1:],env[1:]

def design(phase,peak_ratio,seed):
    B0,B1,F,env=basis(phase);base=np.r_[1.,np.zeros(4)]
    forceQ=F.T@F/len(F);C0=float(base@forceQ@base);peak0=np.max(abs(F@base))
    def calc(c,G=R):
        z0=B0@c;z1=B1@c;return weights(z0,z1,G)
    wbase=calc(base);J0=wbase['J']
    def objective(x):return -calc(np.r_[1.,x])['J']/J0
    freq=np.fft.rfftfreq(len(F),.001)
    Ffft=np.fft.rfft(F,axis=0);outband=(freq<15)|(freq>30)
    spectral_weight=np.full(len(freq),2.);spectral_weight[0]=1.
    if len(F)%2==0:spectral_weight[-1]=1.
    def constraints(x):
        c=np.r_[1.,x];power=spectral_weight*abs(Ffft@c)**2
        leak=np.sum(power[outband])/max(np.sum(power),1e-15)
        return np.r_[1-c@forceQ@c/C0,env@c,peak_ratio-abs(F@c)/peak0,.01-leak]
    rng=np.random.default_rng(seed);best=None
    for start in [np.zeros(4)]+[rng.normal(0,.15,4) for _ in range(5)]:
        res=minimize(objective,start,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                     options=dict(maxiter=300,ftol=1e-10))
        if min(constraints(res.x))>=-1e-6 and (best is None or res.fun<best.fun):best=res
    c=np.r_[1.,best.x];wp=calc(c);before=weights(B0@base,B1@base,M);after=weights(B0@c,B1@c,M)
    BB=B1+(N*TC/2)*B0
    # Global upper bound in the enlarged five-dimensional linear span; no peak/endpoint restrictions.
    eig=eigh(BB.T@R@BB,forceQ,eigvals_only=True)
    bound=float(eig[-1]*C0/J0)
    vals=dict(phase=phase,peak_ratio=peak_ratio,force_outband_limit=.01,coefficients=best.x.tolist(),
              fisher_information_gain=float(wp['J']/J0),force_rms_ratio=float(np.sqrt(c@forceQ@c/C0)),
              peak_force_ratio=float(np.max(abs(F@c))/peak0),
              v2_actual_noise_information_gain=float((before['lv']@before['lv'])/(after['lv']@after['lv'])),
              endpoint_k_variance_ratio=float((wp['qf']@wp['qf'])/(wbase['qf']@wbase['qf'])),
              optimistic_span_information_bound=bound,converged=bool(best.success))
    freq=np.fft.rfftfreq(len(F),.001)
    def outband(c):
        ff=spectral_weight*abs(np.fft.rfft(F@c))**2
        return float(np.sum(ff[(freq<15)|(freq>30)])/max(np.sum(ff),1e-15))
    vals.update(outside_15_30_force_power=outband(c),baseline_outside_15_30_force_power=outband(base))
    return vals

def main():
    records=[recorded_case(i) for i in [0,12,24,36]]
    for r in records:print(r['group'],r['segments']['all'],flush=True)
    designs=[design(phase,cap,17) for phase in [0.,np.pi/3] for cap in [1.,1.25]]
    for r in designs:print('design',r,flush=True)
    old=json.loads((HERE/'results/v2_verify.json').read_text())['trend']
    requirements={g:dict(sd=x['noise_sd'],reference_rate=x['true_p90'],
                    information_factor_snr1=(x['noise_sd']/x['true_p90'])**2,
                    information_factor_mean_1p96sd=(1.96*x['noise_sd']/x['true_p90'])**2)
                  for g,x in old.items()}
    out=dict(source_commit='a2d1b977f92451f1b96ad588ccffb758591fa79d',records=records,designs=designs,
             requirements=requirements,conditions='Gaussian fixed-regressor conditional noise; design uses exact sampled H and polynomial nuisance. No new closed loop or control changes.')
    OUT.write_text(json.dumps(out,indent=1))
if __name__=='__main__':main()
