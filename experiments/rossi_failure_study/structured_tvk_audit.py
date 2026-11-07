"""Causal, same-trajectory audit of time-varying stiffness versus motion.

This is an identification study, not a new controller. The unchanged Huber-D4
policy generates every physical trajectory. Window solvers receive only the
existing preprocessed Y,z and a slope column constructed through the same
within-tick measurement operator. Truth is used only for scoring/attribution.
Run with --collect first (Numba compilation may dominate the first episode),
then --analyse. Results and paired physical logs are reproducible from seeds.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from numba import njit
from scipy.linalg import solveh_banded

import sweeps  # installs the repository's shared import paths
import common
import method3 as M
import scenes3 as S3
import physics_plant3 as P3

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'/'structured_tvk_audit'
CONDITIONS = ('rest', 'breath', 'shift_fast', 'rib', 'cough', 'combo')
BASE = 611700
REFERENCE_COMMIT = '858f3c128589e26c1682b5f55874df12491c535c'


@njit(cache=True)
def slope_column(log, j0, m, ns, eb, afc, te, ye, tw):
    """H[(t-T)(y-yref)/tw], using _build's exact encoder interpolation.

    This makes the local multiplicative-gain model measurement-consistent.
    It does NOT make a memoryless gain an exact model of viscoelastic tissue.
    """
    out = np.zeros(m)
    aa = 1.-afc
    ebn = eb**ns
    for jj in range(m):
        j = j0+jj
        for mm in range(1, ns+1):
            w = afc*aa**(ns-mm)
            fr = (1.-eb**(mm-1))/(1.-ebn) if ebn < 1. else (mm-1)/ns
            yy = log[j,2]+fr*(log[j+1,2]-log[j,2])-ye
            tau = (log[j,0]+(mm-1)*.001-te)/tw
            out[jj] += w*tau*yy
    return out


@njit(cache=True)
def two_pole_moments(log,j0,m,ns,afc,ndel,t1,t2,te,ye,tw):
    """Endpoint-conditioned two-pole interpolation, using calibration only.

    For a held command c: y(s)=c+(y0-c)h0(s)+v0*h1(s).
    y1 determines v0. Both measured endpoints are retained, so a small
    calibrated pole error does not accumulate over the full episode.
    """
    Y=np.zeros(m);G=np.zeros(m);Tc=ns*.001;aa=1.-afc
    h0T=(t1*math.exp(-Tc/t1)-t2*math.exp(-Tc/t2))/(t1-t2)
    h1T=t1*t2*(math.exp(-Tc/t1)-math.exp(-Tc/t2))/(t1-t2)
    for jj in range(m):
        j=j0+jj;y0=log[j,2];y1=log[j+1,2];c=log[j-ndel,11]
        v0=(y1-c-(y0-c)*h0T)/h1T
        for mm in range(ns):
            s=mm*.001;h0=(t1*math.exp(-s/t1)-t2*math.exp(-s/t2))/(t1-t2)
            h1=t1*t2*(math.exp(-s/t1)-math.exp(-s/t2))/(t1-t2)
            yy=c+(y0-c)*h0+v0*h1-ye;w=afc*aa**(ns-1-mm)
            Y[jj]+=w*yy;G[jj]+=w*(log[j,0]+s-te)/tw*yy
    return Y,G


def collect(cond, replica, duration=6., stride=4):
    tissue = 'mixed' if cond in ('rib', 'combo') else 'relaxed'
    sc = S3.cell_scene3(cond, tissue, 'nominal', replica, base=BASE, T=duration,
                       overrides=dict(force=4., slide='long', dynamics='neural_white', phase=14.))
    pol = M.build('H', sc)
    est = pol.estimator
    original = est.estimate
    windows = []

    def record(o, n, policy):
        answer = original(o, n, policy)
        j0 = max(policy.ndel+1, n-1-est.nw)
        count = n-1-j0
        if count != est.nw or o.log[n-1,0] < 1. or not np.isfinite(answer[0]):
            return answer
        reg = est._log(o,n)
        te, ye = reg[n-1,0], reg[n-1,2]
        _, Y, z, keep = common._build(reg,j0,count,0,policy.ns,policy.eb,policy.afc,
                                    policy.ndel,te,ye,est.nw*policy.Tc,policy.params['enc_y'])
        if not keep.all():
            return answer
        G = slope_column(reg,j0,count,policy.ns,policy.eb,policy.afc,te,ye,est.nw*policy.Tc)
        # Sensing calibration is already available to the original Method;
        # no tissue, body, true stiffness, or future observations enter here.
        sensing=policy.step
        t1=getattr(sensing,'t1',policy.tau-P3.ACTUATOR_T2)
        t2=getattr(sensing,'t2',P3.ACTUATOR_T2)
        rawY,rawG=two_pole_moments(o.log,j0-1,count+1,policy.ns,policy.afc,
                                  policy.ndel,t1,t2,te,o.log[n-1,2],est.nw*policy.Tc)
        metric,_=est._matrix(count)
        carrier=np.exp(-2j*np.pi*20.8*o.log[j0:j0+count,0])
        _,old_rawY,_,_=common._build(o.log,j0,count,0,policy.ns,policy.eb,policy.afc,
                                    policy.ndel,te,o.log[n-1,2],est.nw*policy.Tc,policy.params['enc_y'])
        yc=carrier@(metric@rawY[1:]);yo=carrier@(metric@old_rawY)
        # Preserve the baseline's causal phase tracker; transform its reference
        # when the displacement operator changes. Do not relearn phase from
        # force residuals during coughs. This is a measured-displacement-only
        # reference correction, with no access to contact phase truth.
        phase_correction=float(np.angle(yo/yc)) if abs(yc)>1e-15 else 0.
        phi=float(np.clip(policy.phi+phase_correction,0.,.5))
        diff=(1.-np.exp(-2j*np.pi*20.8*policy.Tc))/policy.Tc
        tan=math.tan(phi);tau_h=tan/max(diff.imag-tan*diff.real,1.)
        mag=abs(1.+tau_h*diff)
        correctedY=(rawY[1:]+tau_h/policy.Tc*np.diff(rawY))/mag
        correctedG=(rawG[1:]+tau_h/policy.Tc*(np.diff(rawG)-policy.Tc/(est.nw*policy.Tc)*rawY[:-1]))/mag
        if n%stride:return answer
        windows.append((n-1,j0,Y.copy(),z.copy(),G,float(policy.params['noise']),
                        float(est.rho),np.array(answer),correctedY,correctedG,phi))
        return answer

    est.estimate = record
    start = time.perf_counter()
    log = P3.episode(sc,pol,9100+replica)
    if log is None or not windows:
        raise RuntimeError(f'No usable trajectory: {cond} {replica}')
    OUT.mkdir(parents=True,exist_ok=True)
    path = OUT/f'{cond}_{replica}.npz'
    np.savez_compressed(path,log=log,idx=np.array([w[0] for w in windows]),
                        j0=np.array([w[1] for w in windows]),Y=np.array([w[2] for w in windows]),
                        z=np.array([w[3] for w in windows]),G=np.array([w[4] for w in windows]),
                        noise=np.array([w[5] for w in windows]),rho=np.array([w[6] for w in windows]),
                        baseline=np.array([w[7] for w in windows]),
                        correctedY=np.array([w[8] for w in windows]),correctedG=np.array([w[9] for w in windows]),
                        correctedPhi=np.array([w[10] for w in windows]),
                        events=np.array([[e[1],e[2]] for e in sc['events']]).reshape(-1,2),
                        sc_act=sc['act'],sc_u=sc['u'])
    result=dict(cond=cond,replica=replica,seconds=time.perf_counter()-start,
                windows=len(windows),file=path.name,source_commit=REFERENCE_COMMIT,
                controller=M.metric(log,pol,sc),
                scene={k:sc[k] for k in ('region','subject','probe','contact_law','m','k','c','v','cover','kr')})
    (OUT/f'{cond}_{replica}.json').write_text(json.dumps(result,indent=2,allow_nan=True),encoding='utf-8')
    print(json.dumps(result),flush=True)


def huber_fit(Y,z,G,noise,rho,threshold=np.inf,iters=6):
    """Same IRLS as baseline; optional affine gain with L1 slope penalty.

    threshold is a fixed multiple of noise*sqrt(J0), calculated once using the
    quadratic D4 metric; it is an optimization penalty, not a valid calibrated
    hypothesis test. G uses normalized time, so its coefficient is a gain change.
    """
    m=len(Y);D=np.diff(np.eye(m),4,axis=0)
    w=np.ones(m-4); scale=noise/math.sqrt(rho)
    basis=np.zeros((5,m,len(D)))
    for q in range(5):basis[q,:m-q]=rho*(D[:,q:]*D[:,:m-q]).T
    basis=basis.reshape(-1,len(D))
    rhs=np.c_[Y,G,z]
    eta=None
    for it in range(iters):
        ab=(basis@w).reshape(5,m);ab[0]+=1.
        smooth=solveh_banded(ab,rhs,lower=True,check_finite=False)
        filtered=rhs-smooth
        info=float(Y@filtered[:,0]); cross=float(Y@filtered[:,1])
        J=max(float(G@filtered[:,1])-cross*cross/info,1e-30)
        score=float(G@filtered[:,2])-cross*float(Y@filtered[:,2])/info
        if eta is None:eta=threshold*noise*math.sqrt(J)
        change=0. if math.isinf(threshold) else math.copysign(max(abs(score)-eta,0.),score)/J
        k=(float(Y@filtered[:,2])-cross*change)/info
        d=smooth[:,2]-k*smooth[:,0]-change*smooth[:,1]
        if it != iters-1:
            w=np.minimum(1.,1.345/np.maximum(np.abs(D@d)/scale,1e-12))
    return k,change,score/(noise*math.sqrt(J)),d


def summary(error,mask):
    e=error[mask]
    if not len(e):return None
    return dict(n=len(e),rrmse=float(np.sqrt(np.mean(e*e))),
                bias=float(np.mean(e)),p90_abs=float(np.percentile(np.abs(e),90)))


def analyse():
    rows=[];max_equiv=0.
    thresholds=(0.,1.,2.,3.,5.,8.,np.inf)
    for cond in CONDITIONS:
        for replica in range(3):
            f=OUT/f'{cond}_{replica}.npz'
            if not f.exists():continue
            a=np.load(f);log=a['log'];idx=a['idx'];truth=log[idx,21]
            event=np.zeros(len(idx),bool)
            for lo,hi in a['events']:
                event|=(log[idx,0]>=lo)&(log[idx,0]<=hi+.2)
            kval=np.array([log[j:j+24,21] for j in a['j0']])
            relspan=np.ptp(kval,axis=1)/truth
            rr=np.linspace(-1.,0.,24)
            design=np.c_[np.ones(24),rr]
            affine=design@np.linalg.lstsq(design,kval.T,rcond=None)[0]
            curve=np.max(np.abs(kval-affine.T),axis=1)/truth
            base=a['baseline'][:,0]
            rows.append(dict(cond=cond,replica=replica,method='Huber-current',
                             all=summary((base-truth)/truth,np.ones(len(truth),bool)),
                             event=summary((base-truth)/truth,event),
                             non_event=summary((base-truth)/truth,~event),
                             change_p50_p90=np.percentile(relspan,[50,90]).tolist(),
                             affine_remainder_p50_p90=np.percentile(curve,[50,90]).tolist()))
            for label,Yset,Gset in [('affine',a['Y'],a['G']),('matched',a['correctedY'],a['correctedG'])]:
              for tau in thresholds:
                est=[];change=[];scores=[]
                for Y,z,G,sigma,rho in zip(Yset,a['z'],Gset,a['noise'],a['rho']):
                    kh,delta,score,d=huber_fit(Y,z,G,sigma,rho,tau)
                    est.append(np.clip(kh,200.,3e5));change.append(delta);scores.append(score)
                est=np.array(est)
                if math.isinf(tau) and label=='affine':max_equiv=max(max_equiv,float(np.max(np.abs(est-base)/np.maximum(base,1.))))
                rows.append(dict(cond=cond,replica=replica,method=label+'-'+str(tau),
                                 all=summary((est-truth)/truth,np.ones(len(truth),bool)),
                                 event=summary((est-truth)/truth,event),
                                 non_event=summary((est-truth)/truth,~event),
                                 fraction_nonzero=float(np.mean(np.abs(change)>1e-8))))
            print(f'analysed {cond} {replica}',flush=True)
    result=dict(source_commit=REFERENCE_COMMIT,baseline_reconstruction_max_relative=max_equiv,
                protocol='Exploratory same-physical-trajectory replay; replica 0 is used for parameter selection. Replicas 1/2 were also inspected during measurement-model revision and are not final independent validation. No controller replacement.',rows=rows)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(baseline_reconstruction_max_relative=max_equiv,rows=len(rows))),flush=True)


def mechanisms():
    """Truth is used here only to characterize the physical test distribution."""
    rows=[];all_gains=[]
    for cond in CONDITIONS:
        spans=[];curves=[];event_curves=[];cal=[]
        for replica in range(3):
            a=np.load(OUT/f'{cond}_{replica}.npz');log=a['log'];idx=a['idx'];truth=log[idx,21]
            good=(log[:,0]>=1.)&(log[:,18]>.5)&(log[:,21]>200.)
            all_gains.extend(log[good,21])
            kval=np.array([log[j:j+24,21] for j in a['j0']])
            X=np.c_[np.ones(24),np.linspace(-1.,0.,24)]
            fit=X@np.linalg.lstsq(X,kval.T,rcond=None)[0]
            curvature=np.max(np.abs(kval-fit.T),axis=1)/truth
            spans.extend(np.ptp(kval,axis=1)/truth);curves.extend(curvature)
            ev=np.zeros(len(idx),bool)
            for lo,hi in a['events']:ev|=(log[idx,0]>=lo)&(log[idx,0]<=hi+.2)
            event_curves.extend(curvature[ev])
            tissue='mixed' if cond in ('rib','combo') else 'relaxed'
            sc=S3.cell_scene3(cond,tissue,'nominal',replica,base=BASE,T=6.,
                             overrides=dict(force=4.,slide='long',dynamics='neural_white',phase=14.))
            omega=2*np.pi*20.8;kt=float(np.median(log[:,21]))
            cal.append(dict(replica=replica,mass_true=sc['tool_mass'],mass_calibrated=sc['tool_mass_hat'],
                            approximate_real_gain_bias=-(sc['tool_mass']-sc['tool_mass_hat'])*omega**2,
                            approximate_relative_calibration_halfwidth=.05*sc['tool_mass_hat']*omega**2/kt))
        rows.append(dict(cond=cond,window_relative_span_p50_p90_p99=np.percentile(spans,[50,90,99]).tolist(),
                         best_affine_max_remainder_p50_p90_p99=np.percentile(curves,[50,90,99]).tolist(),
                         event_affine_remainder_p50_p90=np.percentile(event_curves,[50,90]).tolist() if len(event_curves) else None,
                         mass_diagnostic=cal))
    omega=2*np.pi*20.8;wb=2*np.pi*10.
    minimum_impedance=15.*np.sqrt((omega**2-wb**2)**2+(.5*wb*omega)**2)
    gains=np.asarray(all_gains);q=gains/minimum_impedance
    body=dict(assumptions='Frozen local LTI contact and passive body; m=15..40 kg, natural frequency 7..10 Hz, damping ratio .25...5 from the simulation, not measured human bounds.',
              minimum_impedance_N_per_m=float(minimum_impedance),
              displacement_per_force_um_per_N=float(1e6/minimum_impedance),
              observed_gain_p50_p90_p99=np.percentile(gains,[50,90,99]).tolist(),
              relative_passive_gain_error_bound_p50_p90_p99=np.percentile(q/(1.-q),[50,90,99]).tolist())
    result=dict(caution='Simulated dynamics and activation trajectories; not dynamic human ground-truth validation. Mass formulas are single-line approximations, not full-plant bias attribution.',passive_body_line_bound=body,rows=rows)
    (OUT/'mechanisms.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result),flush=True)


def verify_recording():
    rows=[]
    for cond in ('rest','rib','cough'):
        tissue='mixed' if cond=='rib' else 'relaxed'
        sc=S3.cell_scene3(cond,tissue,'nominal',0,base=BASE,T=6.,
                         overrides=dict(force=4.,slide='long',dynamics='neural_white',phase=14.))
        pol=M.build('H',sc);log=P3.episode(sc,pol,9100)
        stored=np.load(OUT/f'{cond}_0.npz')['log']
        assert np.array_equal(log,stored),cond
        finite=np.isfinite(log)&np.isfinite(stored)
        err=float(np.max(np.abs(log[finite]-stored[finite])))
        rows.append(dict(cond=cond,bitwise_identical=True,max_absolute_difference=err,
                         comparison='np.array_equal, including equal infinity patterns; finite values identical'))
    (OUT/'recording_equivalence.json').write_text(json.dumps(rows,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(rows),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--collect',action='store_true');p.add_argument('--analyse',action='store_true')
    p.add_argument('--mechanisms',action='store_true')
    p.add_argument('--verify-recording',action='store_true')
    p.add_argument('--cond',choices=CONDITIONS);p.add_argument('--replica',type=int);p.add_argument('--duration',type=float,default=6.)
    a=p.parse_args()
    if a.collect:
        for cond in (CONDITIONS if a.cond is None else (a.cond,)):
            for replica in (range(3) if a.replica is None else (a.replica,)):
                collect(cond,replica,a.duration)
    if a.analyse:analyse()
    if a.mechanisms:mechanisms()
    if a.verify_recording:verify_recording()


if __name__=='__main__':main()
