"""Identical-observation, simulation-only estimator forensics. No source edits.

python -X utf8 experiments/rossi_failure_study/id_layer_root_cause.py --workers 6
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import ast
import copy
import hashlib
import json
import math
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from types import MethodType
import numpy as np
from numba import njit
from scipy.signal import hilbert
from scipy.stats import spearmanr
import low_estimate_forensics as L
import corrected_comparison as C
import rib_transition_eval as R
import memory_segment_eval as SE
import memory_filters2 as F

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE/'results/id_layer_root_cause.json'
FIG = ROOT/'research/figures/id_layer'
MEM = C.MEMORY
METHODS = ('raw', 'gated')+MEM
COLORS = dict(zip(METHODS, ('#777777','#000000','#0072B2','#56B4E9','#009E73','#E69F00','#CC79A7','#D55E00','#882255')))


def stats(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x): return dict(n=0, bias=None, sd=None, rmse=None)
    m = np.median(x)
    return dict(n=len(x), bias=float(m), sd=float(1.4826*np.median(abs(x-m))),
                rmse=float(np.sqrt(np.mean(x*x))), mean=float(np.mean(x)),
                min=float(x.min()), max=float(x.max()), below_1=float(np.mean(x<0)),
                below_085=float(np.mean(x<np.log(.85))))


def factory(kind):
    f = F.FACTORIES[kind]()
    if kind in ('imm','ewa'):
        f.update = SE.bind_source(f, F.IMM.update if kind=='imm' else F.EWA.update, [
            ('innovation_var = prior+R','innovation_var = prior+R/4.'),
            ('-.5*(', '-.5*.25*('), ('K = prior/innovation_var','K = prior/(prior+R)')])
    return f


def filter_replay(kind, windows, count, white=None):
    f = factory(kind); out = np.full(count,np.nan); pred=out.copy()
    last=None; m=None
    for w in windows:
        n=w['n']; k=w['k_hat']; sig=w['sigma']
        valid=np.isfinite(k) and k>0 and np.isfinite(sig) and sig>0 and w['YMY']>0
        if valid and (last is None or n-last>=6):
            y=math.log(k) if white is None else float(white[n-1])
            m,p=f.update(y,4*(sig/k)**2,0.,.048 if last is None else (n-last)*.008)
            last=n
        if m is not None:
            out[n-1]=np.clip(m,math.log(200),math.log(300000)) if white is None else m
            pred[n-1]=p
    # Causal hold also at ticks when the host did not request a window.
    for j in range(1,count):
        if not np.isfinite(out[j]): out[j]=out[j-1]; pred[j]=pred[j-1]
    return out,pred


def allan(loghat, mask):
    runs=R.spans(mask); longest=max((b-a for a,b in runs),default=0)
    scales=np.unique(np.rint(np.geomspace(25,max(25,longest),16)).astype(int))
    rows=[]
    for m in scales:
        differences=[]
        for a,b in runs:
            x=loghat[a:b]
            if len(x)<2*m or not np.isfinite(x).all(): continue
            means=(np.r_[0.,np.cumsum(x)][m:]-np.r_[0.,np.cumsum(x)][:-m])/m
            d=means[m:]-means[:-m]
            if m==25:d=d[::m]  # shortest scale: non-overlapping adjacent blocks
            differences.extend(d.tolist())
        rows.append(dict(T=m*.008,pairs=len(differences),adev=float(np.sqrt(np.mean(np.square(differences))/2)) if differences else None))
    return rows


def event_scores(log, hats, events):
    t=log[:,0]; truth=np.log(np.maximum(log[:,21],1e-30)); rows=[]
    for ev in events:
        a=ev['start']; b=np.searchsorted(t,ev['t1']+.5+1e-9)
        ix=np.arange(a,b); pre=(t>=ev['t0']-.5-1e-9)&(t<ev['t0']-1e-9)
        taus=np.r_[np.arange(0,1.5,.008),1.5]
        shifted=np.array([np.interp(t[ix]-tau,t,truth,left=np.nan,right=np.nan) for tau in taus])
        row=dict(ev,channels={})
        for kind,h in hats.items():
            e=h-truth; offset=float(np.nanmedian(e[pre])); local=e-offset
            def lag(off):
                loss=np.nanmean(abs(h[ix][None,:]-shifted-off),axis=1)
                return float(taus[np.nanargmin(loss)])
            row['channels'][kind]=dict(offset=offset,lag=lag(offset),lag_without_offset_removal=lag(0),
                settling=R.hold_time(local,t,ev['t1'],.1),settling_unadjusted=R.hold_time(e,t,ev['t1'],.1),
                adjusted=stats(local[ix]),unadjusted=stats(e[ix]),
                transition_only=stats(e[a:ev['end']+1]),interval_end=float(t[b-1]))
        rows.append(row)
    return rows


def fit(t,x):
    s=t-t.mean(); A=np.c_[np.cos(2*np.pi*20.8*t),np.sin(2*np.pi*20.8*t),np.ones(len(t)),s,s*s,s**3]
    v=np.linalg.lstsq(A,x,rcond=None)[0]
    return complex(v[0],-v[1])


def sampled_gain(case,commands,log):
    """Central-difference noise-free sampled y->true axial F at 20.8 Hz.

    Diagnostic only: a +/-1 um held-command sinusoid around the recorded
    commands. This does not feed any replayed estimator or change its data.
    Local fits need not have a stationary plateau; their truth span is retained.
    """
    L.O.setup();sc=L.MC.scene_for(tuple(case['cell']),C.BASE,case['variant'])
    log=np.asarray(log,float);command=np.asarray(commands);tt=log[:,0]
    delta=1e-6*np.sin(2*np.pi*20.8*tt)
    plus=L.P3.replay(sc,command+delta,case['seed'],sensor=False)
    minus=L.P3.replay(sc,command-delta,case['seed'],sensor=False)
    records=[]
    for start in range(125,len(log)-62,31):
        ix=np.arange(start,start+62)
        if min(plus[ix,21].min(),minus[ix,21].min())<=0:continue
        dy=plus[ix,2]-minus[ix,2];df=plus[ix,18]-minus[ix,18]
        Y=fit(tt[ix],dy);Z=fit(tt[ix],df)
        if abs(Y)<1e-10:continue
        k=float(np.median(log[ix,21]));g=abs(Z/Y)
        records.append(dict(t=float(tt[start]),sampled_gain=float(g),logged_k=k,
            log_gap=float(np.log(g/k)),truth_log_span=float(np.log(log[ix,21].max()/log[ix,21].min())),
            encoder_harmonic_amplitude_m=float(abs(Y)/2),gain_phase_deg=float(np.angle(Z/Y)*180/np.pi),
            plateau_only=False))
    return records


# A frozen local affine constitutive law, not just a frozen spatial lookup.
# c[36:] = enabled, depth reference, elastic intercept, tangent, viscous
# coefficient, theta, friction velocity, cop, contact length, cover.
@njit(cache=False)
def frozen_contact(x,d,theta,vc,geom,qfield,efield,ribs,activation,c):
    if c[36]: return c[38]+c[39]*(d-c[37]),c[39],c[43],c[44],c[45]
    return L.P3._contact(x,d,theta,vc,geom,qfield,efield,ribs,activation)


def freeze_tick():
    def transform(fn):
        class Edit(ast.NodeTransformer):
            def visit_Call(self,node):
                if isinstance(node.func,ast.Name) and node.func.id=='_contact':
                    node.func=ast.Name('frozen_contact',ast.Load());node.args.append(ast.Name('c',ast.Load()))
                return self.generic_visit(node)
        Edit().visit(fn)
        class Assign(ast.NodeTransformer):
            def visit_Assign(self,node):
                self.generic_visit(node)
                target=ast.unparse(node.targets[0]); value=ast.unparse(node.value)
                replacement=None
                if target=='theta':replacement='c[41] if c[36] else '+value
                elif target=='rate':replacement='1. if c[36] else '+value
                elif target=='fn' and 'fe + qsum' in value:
                    replacement='max(fe+qsum+c[40]*(3*d-4*s[13]+s[22])/(2*dt),0.) if c[36] else '+value
                elif target=='fn' and 'hereditary_force' in value:
                    replacement='max(hereditary_force+c[40]*normal_velocity,0.) if c[36] else '+value
                elif target=='hc_factor':replacement='1. if c[36] else '+value
                elif target=='kc':replacement='complex(tangent,w*c[40]) if c[36] else '+value
                elif target=='derivative':replacement='ke*(1.+np.sum(weights))+1.5*c[40]/dt if c[36] else '+value
                if replacement:node.value=ast.parse(replacement,mode='eval').body
                return node
        Assign().visit(fn)
        # Constant friction slope as well as angle is needed for constant axial |K|.
        class Friction(ast.NodeTransformer):
            def visit_BinOp(self,node):
                if ast.unparse(node)=='c[17] * vc':
                    return ast.parse('c[17]*(c[42] if c[36] else vc)',mode='eval').body
                return self.generic_visit(node)
            def visit_Call(self,node):
                if isinstance(node.func,ast.Name) and node.func.id=='_friction':
                    node.args[1]=ast.parse('c[42] if c[36] else vc',mode='eval').body
                return self.generic_visit(node)
        Friction().visit(fn)
    # private_function retains the original globals; add dependency in a copied namespace.
    fn=L.private_function(L.P3._tick.py_func,transform)
    import types
    fn=types.FunctionType(fn.__code__,dict(fn.__globals__,frozen_contact=frozen_contact),fn.__name__)
    return njit(cache=False)(fn)


def plant_states(sc,commands,seed):
    P=L.P3;pl=dict(P.PLATFORM,tau_a=sc.get('tau_a_true',P.PLATFORM['tau_a']))
    s,buf=P._init(sc,pl);c=P._constants(sc,pl);alpha,tau,_=P.memory_arrays(sc)
    noise=P._noise(np.random.default_rng(seed),int(sc['T']*1000)+10,False)
    geom=P._geometry(sc); states=[];buffers=[];ks=[]
    log=np.zeros((len(commands)+2,26)); n=kk=0
    for command in commands:
        states.append(s.copy());buffers.append(buf.copy());ks.append(kk)
        seen=s[5]; n0=n
        n,kk=P._tick(s,buf,command,c,geom,sc['q_field'],sc['E_field'],sc['ribs'],sc['h_s'],sc['xcmd'],sc['u'],sc['fext'],sc['act'],sc['rip'],alpha,tau,noise,kk,log,n)
        assert n>=0
        log[n0,5]=seen
    return log[:n],(states,buffers,ks,c,geom,alpha,tau,noise)


def frozen_window(sc,commands,log,hi,states,w,tick):
    P=L.P3;n=w['n'];a=w['j0'];end=n-1
    ss,bb,kk,c0,geom,alpha,tau,noise=states
    s=ss[a].copy();buf=bb[a].copy();c=np.r_[c0,np.zeros(10)]
    ref=hi[end*16];x=ref[26];d=ref[6];theta=ref[17];vc=ref[16]
    fe,ke,cop,length,cover=P._contact(x,d,theta,vc,geom,sc['q_field'],sc['E_field'],sc['ribs'],0.)
    ct=np.cos(theta);st=np.sin(theta)
    friction=P._friction(1.,vc,length,c0,geom)-c0[17]*vc
    axial=ct+friction*st;series=ct*(axial/c0[3]+c0[26]*(friction*ct-st))
    kt=ref[21]*np.exp(1j*ref[22]);kc=kt/(axial*ct-kt*series)
    ww=2*np.pi*20.8;A=1+np.sum(alpha*(1j*ww*tau)/(1+1j*ww*tau))
    tangent=kc.real/A.real;eta=(kc.imag-tangent*A.imag)/ww
    # Preserve the initial force-memory state; change the slope at the start of the window.
    intercept=fe*tangent/max(ke,1e-30)
    c[36:]=[1,d,intercept,tangent,max(eta,0),theta,vc,cop,length,cover]
    s[9]=intercept+tangent*(s[13]-d)
    part=np.zeros((n-a+2,26));p=0;k=kk[a]
    for j in range(a,n):
        seen=s[5];p0=p
        p,k=tick(s,buf,commands[j],c,geom,sc['q_field'],sc['E_field'],sc['ribs'],sc['h_s'],sc['xcmd'],sc['u'],sc['fext'],sc['act'],sc['rip'],alpha,tau,noise,k,part,p)
        assert p>=0
        part[p0,5]=seen
    frozen=log.copy();frozen[a:n]=part[:n-a]
    return frozen,dict(target_k=float(ref[21]),realized_k=stats(np.log(np.maximum(part[:n-a,21],1e-30)/ref[21])),
                       linearized_constitutive=True)


def causes(values):
    if not values:return dict(n=0)
    a=np.asarray(values,float); valid=np.isfinite(a).all(axis=1);a=a[valid]
    if not len(a):return dict(n=0)
    # Signed telescope and absolute-size fractions; these are not independent variance shares.
    size=np.mean(abs(a),axis=0);den=size.sum()
    return dict(n=len(a),names=['noise','disturbance','within_window_k_change','systematic_remainder'],
                mean=np.mean(a,axis=0),median=np.median(a,axis=0),rms=np.sqrt(np.mean(a*a,axis=0)),
                absolute_share=size/den,telescoping_error=float(np.max(abs(a.sum(axis=1)-a.sum(axis=1)))))


def plot_episode(label,log,hats,events,breath,probe,limits=None,suffix=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=log[:,0];truth=np.log(np.maximum(log[:,21],1e-30))
    fig,ax=plt.subplots(3,1,figsize=(15,10),sharex=True,gridspec_kw={'height_ratios':[1.3,1,1]})
    ax[0].plot(t,log[:,21]/1000,color='#AA0000',lw=2,label='True |K(20.8)|')
    for name,h in hats.items():
        ax[0].plot(t,np.exp(h)/1000,color=COLORS[name],lw=.9,label=name,alpha=.85)
        ax[1].plot(t,h-truth,color=COLORS[name],lw=.9)
    for ev in events:
        for aa in ax:aa.axvspan(ev['t0'],ev['t1'],color='purple',alpha=.12)
    ax[1].axhline(0,color='k',lw=.6);ax[1].axhline(np.log(.85),color='red',ls=':',lw=.8)
    ax[0].set_ylabel('Stiffness (kN/m)');ax[1].set_ylabel('ln(estimate / true K)')
    ax[2].plot(t,1000*log[:,14],label='Exogenous u (mm)',color='#0072B2')
    ax[2].plot(t,1000*breath,label='Breathing (mm)',color='#009E73',ls='--')
    aa=ax[2].twinx();aa.plot(t,log[:,18],label='True force (N)',color='#D55E00')
    aa.set_ylabel('True force (N)')
    pa=ax[2].twinx();pa.spines['right'].set_position(('outward',55))
    pa.plot(t,probe,label='Probe force amplitude (N)',color='#882255',lw=1.1)
    pa.set_ylabel('Probe amplitude (N)',color='#882255');pa.tick_params(axis='y',colors='#882255')
    ax[2].set_ylabel('Motion (mm)');ax[2].set_xlabel('Time (s)')
    ax[0].legend(ncol=5,fontsize=9);ax[2].legend(loc='upper left',fontsize=9)
    handles1,labels1=aa.get_legend_handles_labels();handles2,labels2=pa.get_legend_handles_labels()
    aa.legend(handles1+handles2,labels1+labels2,loc='upper right',fontsize=9)
    for a in ax:a.grid(alpha=.2)
    if limits:ax[0].set_xlim(*limits)
    else:ax[0].set_xlim(t[0],t[-1])
    ax[0].set_title(label+suffix+' — identical recorded observations; shaded transitions')
    fig.tight_layout();path=FIG/(label+suffix+'.png');fig.savefig(path,dpi=150);plt.close(fig)
    return str(path.relative_to(ROOT))


def job(case):
    started=time.perf_counter();label=case['group']+'_'+case['cell'][2]+'_i'+str(case['cell'][-1])
    cfg=dict(label=label,cell=case['cell'],base=C.BASE,variant=case['variant'],seed=case['seed'])
    _,sc,log,windows,commands,steps,contexts=L.episode(cfg,full=True)
    print(label,'recorded',len(log),flush=True)
    count=len(log);t=log[:,0];truth=np.log(np.maximum(log[:,21],1e-30))
    targets=[w['n'] for w in windows];wm={w['n']:w for w in windows}
    plain=L.P3.replay(sc,np.array(commands),case['seed'],sensor=True)
    assert L.max_difference(plain,log)<1e-10,'STOP: command replay differs'
    repro,_,_=L.recorded_estimates(sc,plain,contexts,targets)
    delta=max(abs(repro[n]['k_hat']-wm[n]['k_hat']) for n in targets if np.isfinite(wm[n]['k_hat']))
    assert delta<1e-7,('STOP: raw replay differs',delta)
    # Rerun the actual gates and phase state causally on identical observations.
    p=L.B.build('Huber-D4',L.O.chosen()['Huber-D4'],sc)
    pl=dict(L.P3.PLATFORM,tau_a=sc.get('tau_a_hat',L.P3.PLATFORM['tau_a']));p.reset(pl)
    gated=[]
    for n in range(1,count+1):p.step(L.P3.Obs(plain,n));gated.append(p.k_good)
    g=np.array([s['kgood'] for s in steps]);gg=np.array(gated)
    assert np.array_equal(g,gg,equal_nan=True),'STOP: gated replay differs'
    print(label,'replay checks passed',flush=True)
    raw=np.full(count,np.nan);sig=raw.copy();chi=raw.copy()
    for w in windows:raw[w['n']-1]=w['k_hat'];sig[w['n']-1]=w['sigma'];chi[w['n']-1]=w['chi']
    hats=dict(raw=np.log(raw),gated=np.log(g));pred={}
    for kind in MEM:hats[kind],pred[kind]=filter_replay(kind,windows,count)
    _,masks=SE.segment_masks(log,sc);plateau=masks['plateau'];trans=masks['transition']
    _,_,events,excluded=R.transitions(log)
    breath=np.interp(t,np.arange(len(sc['motion_breath']))*.001,sc['motion_breath'])
    phase=np.angle(hilbert(sc['motion_breath']-np.mean(sc['motion_breath'])))
    phase=np.interp(t,np.arange(len(phase))*.001,np.unwrap(phase))%(2*np.pi)
    probe=np.array([s['probe_amplitude'] for s in steps])
    plateau_stats={kind:stats((h-truth)[plateau]) for kind,h in hats.items()}
    # Conditional white-noise null: fixed observed sigma schedule, iid window errors;
    # filter parameters/likelihoods unchanged. Pool eight deterministic draws.
    rng=np.random.default_rng(case['seed']+770000);white={k:[] for k in METHODS};overlap={k:[] for k in METHODS}
    for rep in range(8):
        z=rng.normal(size=count)*sig/raw
        white['raw'].extend(z[plateau & np.isfinite(z)])
        for kind in MEM:
            zz,_=filter_replay(kind,windows,count,white=z)
            white[kind].extend(zz[plateau & np.isfinite(zz)])
        correlated=np.convolve(rng.normal(size=count+23),np.ones(24)/np.sqrt(24),mode='valid')*sig/raw
        overlap['raw'].extend(correlated[plateau & np.isfinite(correlated)])
        for kind in MEM:
            zz,_=filter_replay(kind,windows,count,white=correlated)
            overlap[kind].extend(zz[plateau & np.isfinite(zz)])
    for kind in METHODS:
        plateau_stats[kind]['white_prediction_sd']=stats(white[kind])['sd'] if kind!='gated' else None
        plateau_stats[kind]['white_overlap_prediction_sd']=stats(overlap[kind])['sd'] if kind!='gated' else None
    phase_stats={k:[dict(phase_bin=j,**stats((h-truth)[plateau&(phase>=j*np.pi/4)&(phase<(j+1)*np.pi/4)])) for j in range(8)] for k,h in hats.items()}
    no,states=plant_states(sc,commands,case['seed'])
    tick=L.diagnostic_tick();hi=L.local_replay(sc,commands,case['seed'],tick,highrate=True,sensor=False)
    assert L.max_difference(hi[::16,:26],no)<1e-10
    selected=[n for n in targets if trans[n-1] or (plateau[n-1] and (n-1)%25==0)]
    noest,nocal,nopol=L.recorded_estimates(sc,no,contexts,selected)
    gauge,gcheck=L.motion_gauge(sc,hi,nocal,nopol,'u')
    motionest,_,_=L.recorded_estimates(sc,gauge,contexts,selected,corrected=True)
    ft=freeze_tick();cf=[];gaps=[]
    print(label,'counterfactual windows',len(selected),flush=True)
    for n in selected:
        w=wm[n]; frozen,fcheck=frozen_window(sc,commands,no,hi,states,w,ft)
        # Same held commands imply identical actuator/encoder trajectory. Thus
        # sensing compensation, H history, and the depth-preserving gauge are
        # identical; only Z changes. Reusing Y is exact, not a fixed-weight fit:
        # array_estimate recomputes all six IRLS iterations on the changed Z.
        a0=w['j0'];ctx=contexts[n]
        assert np.max(abs(frozen[a0:n,2]-no[a0:n,2]))<1e-12
        dz=frozen[:,5]-no[:,5];decay=(1-ctx['afc'])**ctx['ns']
        delta_z=dz[a0+1:n]-decay*dz[a0:n-1]
        Y=np.asarray(motionest[n]['Y']);z=np.asarray(motionest[n]['z'])
        assert len(Y)==len(delta_z)
        fk=L.array_estimate(Y,z+delta_z,nopol)['k_hat']
        a=math.log(w['k_hat']);b=math.log(noest[n]['k_hat']);c=math.log(motionest[n]['k_hat']);d=math.log(fk)
        cf.append(dict(n=n,t=t[n-1],segment='transition' if trans[n-1] else 'plateau',
                       values=[a,b,c,d],parts=[a-b,b-c,c-d,d-truth[n-1]],freeze=fcheck))
    # Alternative order on a deterministic 1/10 subsample: motion before noise.
    fullhi=L.local_replay(sc,commands,case['seed'],tick,highrate=True,sensor=True)
    _,fullcal,fullpol=L.recorded_estimates(sc,plain,contexts,selected)
    fullg,_=L.motion_gauge(sc,fullhi,fullcal,fullpol,'u')
    sub=selected[::10];alternative,_,_=L.recorded_estimates(sc,fullg,contexts,sub,corrected=True)
    for r in cf:
        if r['n'] in alternative:
            v=math.log(alternative[r['n']]['k_hat']);a,b,c,d=r['values']
            r['alternative_motion_first']=[v-c,a-v,c-d,d-truth[r['n']-1]]
    gain=sampled_gain(case,commands,log)
    shares={s:causes([r['parts'] for r in cf if r['segment']==s]) for s in ('plateau','transition')}
    method_shares={}
    for kind,h in hats.items():
        method_shares[kind]={}
        for seg in ('plateau','transition'):
            selected_rows=[r for r in cf if r['segment']==seg]
            arrays=[]
            for r in selected_rows:
                j=r['n']-1;center=math.log(SE.CC.CENTRE) if kind=='gated' else 0.
                arrays.append(r['parts']+[float(h[j]-hats['raw'][j]-center),center])
            aa=np.asarray(arrays)
            if aa.size:
                aa=aa[np.isfinite(aa).all(axis=1)];sz=np.mean(abs(aa),axis=0)
                method_shares[kind][seg]=dict(names=shares[seg]['names']+['memory_or_gate_selection','centering'],
                    mean=np.mean(aa,axis=0),absolute_share=sz/sz.sum(),n=len(aa))
    memory_delta={kind:{s:stats((h-hats['raw'])[masks[s]]) for s in ('plateau','transition')} for kind,h in hats.items() if kind!='raw'}
    excitation={s:dict(n=int((masks[s]&np.isfinite(raw)).sum()),
        error_probe_spearman=float(spearmanr(abs((hats['raw']-truth)[masks[s]]),probe[masks[s]],nan_policy='omit').statistic),
        error_relative_sigma_spearman=float(spearmanr(abs((hats['raw']-truth)[masks[s]]),(sig/raw)[masks[s]],nan_policy='omit').statistic)) for s in ('plateau','transition')}
    figures=[plot_episode(label,log,hats,events,breath,probe)]
    # Two largest transitions and two largest method/tick underestimates per episode.
    for j,ev in enumerate(sorted(events,key=lambda e:e['size'],reverse=True)[:2]):
        figures.append(plot_episode(label,log,hats,events,breath,probe,(max(0,ev['t0']-.6),min(t[-1],ev['t1']+1.5)),f'_transition{j+1}'))
    candidates=[]
    eligible=(t>=1)&(log[:,21]>0)
    for kind,h in hats.items():
        e=np.where(eligible,h-truth,np.nan)
        candidates.extend((e[idx],idx,kind) for idx in np.flatnonzero(np.isfinite(e)))
    centers=[]
    for _,ix,kind in sorted(candidates):
        if all(abs(t[ix]-tt)>.5 for tt in centers):
            centers.append(t[ix]);figures.append(plot_episode(label,log,hats,events,breath,probe,(max(0,t[ix]-.8),min(t[-1],t[ix]+.8)),f'_under{len(centers)}_{kind}'))
        if len(centers)==2:break
    result=dict(case=case,label=label,replay=dict(raw_max_abs_difference=delta,gated_exact=True,commands_max_abs_difference=L.max_difference(plain,log)),
        plateau=plateau_stats,allan=allan(hats['raw'],plateau),phase=phase_stats,
        transitions=event_scores(log,hats,events),excluded_transitions=excluded,causes=shares,method_cause_shares=method_shares,counterfactuals=cf,
        memory_delta=memory_delta,excitation=excitation,sampled_gain=gain,sampled_gap=stats([v['log_gap'] for v in gain]),
        motion_depth_check=gcheck,figures=figures,seconds=time.perf_counter()-started,
        plateau_runs=[dict(t0=float(t[a]),duration=float((b-a)*.008)) for a,b in R.spans(plateau)],
        records=dict(log=log,commands=commands,windows=windows,steps=steps,contexts=contexts,
                     breathing=breath,phase=phase,probe=probe,estimates={k:np.exp(v) for k,v in hats.items()},
                     masks=masks,filter_covariance=pred))
    return L.clean(result)


def summary_figures(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(15,5.5))
    offsets=dict(raw=(-30,-5),gated=(5,5),kal_long=(5,5),kal_short=(5,-15),
                 imm=(-8,22),lepski=(6,-28),ewa=(8,7),forget_time=(5,7),ard_fixed=(-6,-36))
    for kind in METHODS:
        sd=[r['plateau'][kind]['sd'] for r in rows if r['plateau'][kind]['sd'] is not None]
        ev=[e['channels'][kind] for r in rows for e in r['transitions']]
        if sd and ev:
            xx=np.mean(sd);lag=np.median([e['lag'] for e in ev]);under=-min(e['unadjusted']['min'] for e in ev)
            for j,(a,yy) in enumerate(zip(ax,[lag,under])):
                a.scatter(xx,yy,color=COLORS[kind],s=55)
                offset=offsets[kind] if j==1 else ((-30,4) if kind=='raw' else (4,4))
                a.annotate(kind,(xx,yy),xytext=offset,textcoords='offset points',fontsize=9,
                           arrowprops=dict(arrowstyle='-',color=COLORS[kind],lw=.6) if j==1 else None)
    for a in ax:a.set_xlabel('Plateau robust SD (log units; episode mean)');a.grid(alpha=.2)
    ax[0].set_ylabel('Transition lag (s; event median, offset removed)');ax[1].set_ylabel('Largest underestimate (-min log error)')
    ax[1].set_yscale('log');ax[1].set_ylim(.18,5.5)
    for a in ax:a.set_xlim(.02,.15)
    fig.tight_layout();fig.savefig(FIG/'noise_tracking_tradeoff.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(2,3,figsize=(15,8));fig2,bx=plt.subplots(2,3,figsize=(15,8))
    for a,b,r in zip(ax.flat,bx.flat,rows):
        rr=[v for v in r['allan'] if v['adev'] is not None]
        if rr:
            T=np.array([v['T'] for v in rr]);v=np.array([v['adev'] for v in rr]);a.loglog(T,v,'o-',label='Raw Allan deviation');a.loglog(T,v[0]*np.sqrt(T[0]/T),'--',label='1/sqrt(T)')
        else:a.text(.1,.5,'No plateau long enough for two 0.2 s blocks',transform=a.transAxes,wrap=True)
        a.set_title(r['label'],fontsize=9);a.set_xlabel('Averaging time (s)');a.set_ylabel('Allan deviation (log k)');a.grid(alpha=.2)
        for kind in METHODS:
            pp=r['phase'][kind];b.plot(np.arange(8)*45+22.5,[v['bias'] if v['bias'] is not None else np.nan for v in pp],color=COLORS[kind],label=kind)
        b.set_title(r['label'],fontsize=9);b.set_xlabel('Breathing phase (deg)');b.set_ylabel('Median log error');b.grid(alpha=.2)
    for a in ax.flat:
        if a.get_legend_handles_labels()[0]:a.legend(fontsize=8);break
    bx.flat[0].legend(fontsize=7,ncol=3)
    fig.tight_layout();fig.savefig(FIG/'plateau_allan.png',dpi=150);plt.close(fig)
    fig2.tight_layout();fig2.savefig(FIG/'breathing_phase_bias.png',dpi=150);plt.close(fig2)


def finalize(result):
    rows=result['episodes'];summary={}
    for row in rows:
        sample_bias={kind:[] for kind in METHODS}
        log=np.asarray(row['records']['log'],float);t=log[:,0]
        for gain in row['sampled_gain']:
            ix=(t>=gain['t']-1e-9)&(t<gain['t']+62*.008-1e-9)&np.asarray(row['records']['masks']['plateau'])
            for kind in METHODS:
                h=np.asarray(row['records']['estimates'][kind],float)
                sample_bias[kind].extend((np.log(h[ix]/log[ix,21])-gain['log_gap']).tolist())
        row['plateau_relative_sampled_gain']={k:stats(v) for k,v in sample_bias.items()}
        row['plateau_noise_off']=stats([v['values'][1]-math.log(log[v['n']-1,21]) for v in row['counterfactuals'] if v['segment']=='plateau'])
        row['freeze_max_abs_log_gap']=max((max(abs(v['freeze']['realized_k']['min']),abs(v['freeze']['realized_k']['max'])) for v in row['counterfactuals']),default=None)
        # Numerical identity of reused-regressor replay is checked against the
        # original window before applying the constitutive force change.
        row['alternative_order']={s:causes([v['alternative_motion_first'] for v in row['counterfactuals'] if v['segment']==s and 'alternative_motion_first' in v]) for s in ('plateau','transition')}
        row['original_order_same_subsample']={s:causes([v['parts'] for v in row['counterfactuals'] if v['segment']==s and 'alternative_motion_first' in v]) for s in ('plateau','transition')}
    for kind in METHODS:
        pp=[r['plateau'][kind] for r in rows if r['plateau'][kind]['n']]
        ev=[e['channels'][kind] for r in rows for e in r['transitions']]
        def avg(key):
            vv=[p.get(key) for p in pp if p.get(key) is not None]
            return float(np.mean(vv)) if vv else None
        summary[kind]=dict(plateau_episodes=len(pp),plateau={k:avg(k) for k in ('bias','sd','rmse','white_prediction_sd','white_overlap_prediction_sd')},
            transition=dict(n=len(ev),lag=R.dist([v['lag'] for v in ev]),
                settling=R.dist([v['settling']['time'] for v in ev]),settling_censored=sum(v['settling']['censored'] for v in ev),
                mean_error=R.dist([v['unadjusted']['mean'] for v in ev]),
                mean_error_adjusted=R.dist([v['adjusted']['mean'] for v in ev]),
                worst_min=min((v['unadjusted']['min'] for v in ev),default=None),worst_min_adjusted=min((v['adjusted']['min'] for v in ev),default=None),
                below_1=R.dist([v['unadjusted']['below_1'] for v in ev]),below_085=R.dist([v['unadjusted']['below_085'] for v in ev])))
    result['summary']=summary
    result['transition_counts']={size+'_'+direction:sum(e['size_class']==size and e['direction']==direction for r in rows for e in r['transitions']) for size in ('large','small') for direction in ('rising','falling')}
    result['metadata']['white_null_overlap']='Also 24-tick moving-average white noise with the same marginal sigma/k_hat, preserving nominal overlap correlation. Eight fixed draws, finite-sample reference rather than exact Huber noise model.'
    result['metadata']['systematic_gap']='All requested e metrics use logged K. Separate per-episode measured sampled F/y gap and plateau errors relative to that measured gain are supplied; the gap is target-definition/systematic, not estimator error.'
    result['metadata']['plateau_absence']='Empty plateaus, phase bins and Allan scales remain missing; no substituted segments.'
    result['metadata']['sampled_gain']='Central difference of noise-free replays with +/-1um held-command20.8Hz sinusoid around original commands; 0.496s local harmonic+cubic fits, step0.248s. True axial F/encoder y. Does not change estimator observations. Nonstationary truth spans retained; not broadband control gain.'
    result['metadata']['code_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(L.__file__),Path(R.__file__),Path(SE.__file__),Path(F.__file__),Path(C.__file__)]}
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=6,choices=range(1,9));args=parser.parse_args()
    assert not OUT.exists(),'Refuse to overwrite existing results'
    FIG.mkdir(parents=True,exist_ok=True)
    cases=[c for c in C.cases() if (c['group']=='A' and ((c['cell'][2]=='neural_white' and c['cell'][4] in (2,6)) or (c['cell'][2]=='viscous' and c['cell'][4]==6))) or
           (c['group']=='A_slow_thin' and ((c['cell'][2]=='neural_white' and c['cell'][4] in (1,5)) or (c['cell'][2]=='viscous' and c['cell'][4]==5)))]
    assert len(cases)==6
    start=time.perf_counter();rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(job,cases):
            rows.append(row);print('completed',row['label'],round(row['seconds'],1),flush=True)
    summary_figures(rows)
    result=dict(metadata=dict(base=C.BASE,methods=METHODS,centering=SE.CC.CENTRE,
        centering_log=math.log(SE.CC.CENTRE),workers=args.workers,
        replay='All methods share recorded measurements. Estimator-dependent excitation and loop effects are omitted.',
        memory='All valid raw windows, six-tick updates; R=4*(sigma/k_hat)^2. R.make_filter constants and tempered IMM/EWA likelihoods; no tuning, no host centering or gates.',
        white_null='8 deterministic independent Gaussian window-noise draws, conditional observed sigma/k_hat schedule, same filters; ignores overlap. Not posterior covariance.',
        counterfactual_order='noise off, depth-preserving u gauge, frozen endpoint local affine constitutive law, systematic remainder',
        freeze='Local affine tangent and Prony/viscous coefficients plus geometry/friction frozen to endpoint; original commands, initial state/history retained. A linearized intervention, not spatial freeze. Realized log-K range checked.',
        shares='mean absolute log contribution divided by sum; signed mean and RMS also stored. Interactions/order-dependent, not variance explained.',
        scoring='Exact col21 truth. SE masks; R transitions including small events. Event interval t0:t1+0.5, offsets previous0.5s, lag0:1.5s, settling tolerance0.1 for0.2s cap2s.',
        sampled_gain='Noise-off original commands, contiguous plateau 0.496s harmonic+cubic fits at20.8Hz, true force col18 / encoder col2. Conditional spectral diagnostic; not broadband loop gain.'),episodes=rows,seconds=time.perf_counter()-start)
    finalize(result)
    OUT.write_text(json.dumps(L.clean(result),ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
    print('saved',OUT,'seconds',result['seconds'],flush=True)


if __name__=='__main__':main()
