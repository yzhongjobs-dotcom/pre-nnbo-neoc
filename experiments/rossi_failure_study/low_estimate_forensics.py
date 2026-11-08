"""Fixed-case causal forensics. Only this file and low_estimate_forensics.json.

python -X utf8 experiments/rossi_failure_study/low_estimate_forensics.py --workers 6
Instance-bound source instrumentation retains the adopted execution path.
No class/module monkey patches, policy attachments, new scene variants or tuning.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import ast
import copy
import inspect
import json
import math
import textwrap
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import MethodType, SimpleNamespace

import overlay_compare as O
import numpy as np
from scipy.special import ndtr
from scipy.signal import periodogram, find_peaks
from scipy.linalg import solve_banded, solveh_banded
from numba import njit
import tail_channel_diagnosis as T
import memory_closed_loop as MC
import method3 as M
import physics_plant3 as P3
import proposed
import common
from estimator_swap import SwapPolicy
from d_robust import RobustRidge
from contact_compare import ModelPolicy
from ridge_final import CENTRE

HERE = Path(__file__).resolve().parent
OUT = HERE/'results/low_estimate_forensics.json'
B = O.B


def selected_cases():
    cases = list(T.cases())
    rows = json.loads((HERE/'results/memory_closed_loop_true_base.json').read_text())
    grid = [c for c in O.cells() if c[0] == 'breath_rib']
    for i in (1, 4):
        cell = ('breath_rib', 'rib', 'neural_white', 28., i)
        r = next(r for r in rows if r['name'] == 'Huber-true dev' and
                 r['variant'] == 'slow_thin' and r['dynamics'] == cell[2] and r['phase'] == 28. and r['i'] == i)
        cases.append(dict(label=f'dev{i}', cell=cell, base=MC.DEV_BASE,
                          variant='slow_thin', seed=MC.SEED+grid.index(cell),
                          stored_peak=r['peak'], rounded=False))
    cell = ('breath_rib', 'rib', 'neural_white', 28., 2)
    rows = json.loads((HERE/'results/overlay_test.json').read_text())['rows']
    r = next(r for r in rows if r['method'] == 'Huber-D4' and B.key(r) == cell)
    cases.append(dict(label='control', cell=cell, base=674500, variant='nominal',
                      seed=O.sensor_seed(cell), stored_peak=r['peak'], rounded=False))
    return cases


def private_function(fn, transform):
    tree = ast.parse(textwrap.dedent(inspect.getsource(fn)))
    tree.body[0].decorator_list=[]
    transform(tree.body[0]); ast.fix_missing_locations(tree)
    ns = {}
    exec(compile(tree, '<instance low estimate forensics>', 'exec'), fn.__globals__, ns)
    return ns[tree.body[0].name]


class CaptureReturns(ast.NodeTransformer):
    def visit_Return(self, node):
        return ast.copy_location(ast.Return(ast.Call(
            ast.Attribute(ast.Name('self', ast.Load()), '_gate_capture', ast.Load()),
            [ast.Call(ast.Name('locals', ast.Load()), [], []), node.value], [])), node)


def attach_instrumentation(pol, run, windows, commands, steps, contexts):
    def gate_capture(v, result):
        if 'k_hat' not in v:
            n=v['n'];j0=max(pol.ndel+1,n-1-pol.estimator.nw)
            windows.append(dict(n=n,t=float(v['o'].log[n-1,0]),k_hat=np.nan,sigma=np.inf,chi=np.inf,
                                prior_kgood=contexts[n]['kgood'],prior_kbar=contexts[n]['kbar'],
                                precise=None,consistent=None,bound=None,within=None,better=None,ok_meas=None,
                                accepted=False,path='no_valid_estimate',k_ref=None,precision_k_ref=None,z=np.nan,
                                kgood_after=float(pol.k_good),kbar_return=float(result[0]),raw_k=np.nan,YMY=0.,
                                j0=j0,m=max(n-1-j0,0),preproc_tau_h=contexts[n]['tau_h'],preproc_mag_h=contexts[n]['mag_h']))
            return result
        n = v['n']; prior = contexts[n]['kgood']; kh = float(v['k_hat']); sig = float(v['sig'])
        accepted = ('ok_meas' in v and bool(v['ok_meas']))
        r = dict(n=n, t=float(v['o'].log[n-1, 0]), k_hat=kh, sigma=sig,
                 chi=float(v['chi_est']), prior_kgood=prior, prior_kbar=contexts[n]['kbar'],
                 precise=bool(v['precise']), consistent=bool(v['consistent']),
                 bound=v.get('bound'), within=v.get('within'), better=v.get('better'),
                 ok_meas=v.get('ok_meas'), accepted=accepted,
                 path=('measurement:within' if v.get('within') else 'measurement:better' if v.get('better') else 'measurement:no_bound_acquisition') if accepted else
                      ('reject:precision' if not v['precise'] else 'reject:consistency' if not v['consistent'] else 'reject:measurement'),
                 k_ref=float(v['k_ref']), precision_k_ref=kh if run=='gate_prior_free' else float(v['k_ref']), z=(kh-prior)/sig,
                 kgood_after=float(pol.k_good), kbar_return=float(result[0]),
                 raw_k=pol._raw_fit['raw_k'], YMY=pol._raw_fit['YMY'],
                 j0=pol._raw_fit['j0'], m=pol._raw_fit['m'],
                 preproc_tau_h=contexts[n]['tau_h'], preproc_mag_h=contexts[n]['mag_h'])
        windows.append(r)
        return result
    pol._gate_capture = gate_capture
    def swap_transform(fn):
        if run=='gate_prior_free':
            for i,node in enumerate(fn.body):
                if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='precise' for t in node.targets):
                    fn.body[i]=ast.parse("precise = sig <= p['rel_max']*k_hat and (not self.gate_on_prior or k_est > 200.)").body[0]
        CaptureReturns().visit(fn)
    swap = private_function(SwapPolicy._stiffness_est, swap_transform)
    swap = MethodType(swap, pol)
    original_est = pol._stiffness_est
    pol._private_swap = swap
    def model_transform(fn):
        class ReplaceSuper(ast.NodeTransformer):
            def visit_Call(self,node):
                if ast.unparse(node.func)=='super()._stiffness_est':
                    node.func=ast.Attribute(ast.Name('self',ast.Load()),'_private_swap',ast.Load())
                return self.generic_visit(node)
        ReplaceSuper().visit(fn)
    model = MethodType(private_function(ModelPolicy._stiffness_est,model_transform),pol)
    def estimate(o, n):
        contexts[n] = dict(kgood=float(pol.k_good), kbar=float(pol.kbar),
                           tau_h=float(pol.tau_h), mag_h=float(pol.mag_h),
                           ndel=int(pol.ndel), eb=float(pol.eb), ns=int(pol.ns), afc=float(pol.afc))
        # CandidatePolicy's only extra path around SwapPolicy is phase updating.
        result = model(o, n)
        if run in ('kgood_limit','both') and np.isfinite(contexts[n]['kgood']):
            pol.k_good=max(pol.k_good,(1-pol.params['rho_d'])*contexts[n]['kgood'])
        if windows and windows[-1]['n']==n:
            windows[-1]['kgood_after']=float(pol.k_good)
        if np.isfinite(getattr(pol, 'rel', np.nan)) and pol.rel <= pol.params['rel_max']:
            pol.update_phase(o, n)
        return result
    assert 'out = super()._stiffness_est(o, n)' in inspect.getsource(original_est)
    pol._stiffness_est = estimate
    raw_est = RobustRidge.estimate
    def raw_transform(fn):
        # Save raw, unclipped k and the exact final-weight information.
        pos = next(i for i, x in enumerate(fn.body) if isinstance(x, ast.Assign) and
                   ast.unparse(x) == 'k = _clip_k(k)')
        fn.body.insert(pos, ast.parse("policy._raw_fit = dict(raw_k=float(k), YMY=float(info), j0=int(j0)+int(broken[-1]+1 if len(broken) else 0), m=int(m))").body[0])
    pol.estimator.estimate = MethodType(private_function(raw_est, raw_transform), pol.estimator)

    def transform_step(fn):
        final = fn.body[-1]
        assert ast.unparse(final.value) == 'sent + du / self.kbar + (dx - dq) / k_un'
        val = final.value
        if run in ('exc_from_kbar', 'both'):
            val = ast.parse('sent + du/self.kbar + dx/(self.kbar if np.isfinite(self.k_good) else k_un) - dq/k_un', mode='eval').body
        fn.body[-1] = ast.Return(ast.Call(ast.Attribute(ast.Name('self', ast.Load()), '_command_capture', ast.Load()),
                                        [ast.Call(ast.Name('locals', ast.Load()), [], []), val], []))
    step = MethodType(private_function(proposed.MotionAware.step, transform_step), pol)
    def command_capture(v, command):
        k = float(getattr(v['o'], 'truth_k', np.nan)); n = v['n']
        excden = pol.kbar if run in ('exc_from_kbar', 'both') and np.isfinite(pol.k_good) else v['k_un']
        steps.append(dict(n=n, t=float(v['t']), k=k, kbar=float(pol.kbar), kgood=float(pol.k_good),
                          fb=float(v['du']/pol.kbar), ff=float(-v['dq']/v['k_un']),
                          exc=float(v['dx']/excden), k_un=float(v['k_un']),
                          probe_amplitude=float(pol.params['a_force']*pol.exc_gain*k/excden),
                          command=float(command), mode=int(pol.mode)))
        return command
    pol._command_capture = command_capture
    def inherited(o):
        n = o.n
        if n > 60 and n % 25 == 0:
            pol._identify(o.log, n)
        if n % 125 == 0:
            pol.id_hist.append((float(o.log[n-1, 0]), pol.tau_hat, pol.ndel_hat))
        if pol.mode == 1 or o.log[n-1, 5] < .3:
            pol._clear_motion()
        return step(o)
    if isinstance(pol.step, M.ForceSensing):
        pol.step.step = inherited
    else:
        pol.step = inherited
    saved = pol.step
    def outer(o):
        c = saved(o); commands.append(float(c))
        if not steps or steps[-1]['n']!=o.n:
            steps.append(dict(n=o.n,t=float(o.log[o.n-1,0]),k=float(getattr(o,'truth_k',np.nan)),
                              kbar=float(pol.kbar),kgood=float(pol.k_good),fb=0.,ff=0.,exc=0.,
                              k_un=float(pol.k_good),probe_amplitude=0.,command=float(c),mode=int(pol.mode),
                              approach_increment=float(c-o.log[o.n-1,11])))
        return c
    pol.step = outer


def finish_windows(windows, log, steps, nw):
    sm = {r['n']: r for r in steps}
    for r in windows:
        n = r['n']; j = n-1
        # Diagnostic truth seen by the command is the preceding plant sample.
        r['k'] = float(log[max(0,j-1),21]); r['force'] = float(log[j,18])
        truth = log[r['j0']:j,21]
        r['true_k_min'] = float(truth.min()) if len(truth) else np.nan
        r['true_k_max'] = float(truth.max()) if len(truth) else np.nan
        r['ln_k_span'] = float(np.log(truth.max()/truth.min())) if len(truth) and truth.min()>0 else math.inf
        r['rel_sigma'] = r['sigma']/r['k'] if r['k']>0 else math.inf
        if n in sm:
            s = sm[n]; r['probe_force_amplitude'] = s['probe_amplitude']
            r['k_over_kgood'] = r['k']/s['kgood']; r['k_over_kbar'] = r['k']/s['kbar']
            r['k_over_prior_kgood'] = r['k']/r['prior_kgood']
            r['probe_force_amplitude_prior']=r['probe_force_amplitude']*s['kgood']/r['prior_kgood']
            amplitudes=[sm[v]['probe_amplitude'] for v in range(r['j0']+1,n) if v in sm]
            r['probe_force_amplitude_in_window']=distribution(amplitudes)


def event_summary(windows, steps, log, peak_t):
    # All exits are kept, including overestimation (low-side) exits. Initial
    # acquisition is excluded until each channel has entered the stated band.
    exits = []; prev = {}; entered = set()
    for s in steps:
        if s['k'] <= 0 or s['t'] > peak_t: continue
        for channel in ('kgood','kbar'):
            if not np.isfinite(s[channel]): continue
            r = s['k']/s[channel]; inside = .75 <= r <= 1.3
            if inside: entered.add(channel)
            if channel in entered and prev.get(channel,True) and not inside:
                exits.append(dict(n=s['n'], t=s['t'], channel=channel, ratio=r,
                                  direction='underestimate' if r>1.3 else 'overestimate'))
            prev[channel] = inside
    # First upper-band exit precedes the peak, with 10 full predecessor windows.
    eligible = [e for e in exits if e['direction']=='underestimate']
    first = eligible[0] if eligible else None
    if first:
        wi = max((i for i,r in enumerate(windows) if r['n']<=first['n']), default=0)
        history = windows[max(0,wi-10):wi+1]
        first = dict(first, window_n=windows[wi]['n'], window_t=windows[wi]['t'],
                     force_at_event=float(log[first['n']-1,18]))
    else: history=[]
    abnormal=np.flatnonzero((log[:,18]>11.) | (log[:,18]<.5))
    warm_abnormal=np.flatnonzero((log[:,0]>=1.)&((log[:,18]>11.) | (log[:,18]<.5)))
    return dict(first_upper_exit=first, first_either_exit=exits[0] if exits else None,
                first_abnormal_t=float(log[abnormal[0],0]) if len(abnormal) else None,
                first_abnormal_after_1s_t=float(log[warm_abnormal[0],0]) if len(warm_abnormal) else None,
                abnormal_definition='true axial force >11 N or <0.5 N, including before 1 s',
                exits=exits, history=history)


def response(steps, log, event, peak_t):
    if not event: return None
    t0=event['t']; mask=(log[:,0]>=t0)&(log[:,0]<=peak_t)
    # Spectrum of the final second leading to the peak; frequency resolution is retained.
    fs=125.; w=(log[:,0]>=max(t0,peak_t-1.))&(log[:,0]<=peak_t)
    freq,power=periodogram(log[w,18],fs=fs,detrend='linear')
    q=(freq>=5)&(freq<=55)
    hz=float(freq[q][np.argmax(power[q])]) if q.any() else None
    times=log[w,0];forces=log[w,18]
    peaks,_=find_peaks(forces,distance=2,prominence=max(.3,.1*float(np.ptp(forces))))
    intervals=np.diff(times[peaks]);cv=float(np.std(intervals)/np.mean(intervals)) if len(intervals)>1 else None
    oscillation_hz=float(1/np.median(intervals)) if len(intervals)>=2 and cv<.35 else None
    selected=[s for s in steps if t0<=s['t']<=peak_t]
    peak_step=min(steps,key=lambda s:abs(s['t']-peak_t))
    around=[]
    for s in steps:
        if max(t0,peak_t-.3)<=s['t']<=peak_t:
            j=s['n']-1
            around.append(dict(s,force=float(log[j,18]),depth=float(log[j,6]),encoder=float(log[j,2]),
                               applied_command=float(log[j,11]),command_clipped=abs(s['command']-log[j,11])>1e-10))
    return dict(event_to_peak_s=peak_t-t0, frequency_hz=hz,oscillation_frequency_hz=oscillation_hz,
                force_peak_times=times[peaks].tolist(),inter_peak_interval_cv=cv,
                spectral_window_s=float(w.sum()/fs),
                max_increment_mm={c:1000*max(abs(s[c]) for s in selected) for c in ('fb','ff','exc')},
                rms_increment_mm={c:1000*float(np.sqrt(np.mean([s[c]**2 for s in selected]))) for c in ('fb','ff','exc')},
                peak_step=peak_step, force_min=float(log[mask,18].min()),
                approach_ticks=sum(s['mode']==1 for s in selected),
                max_approach_increment_mm=1000*max((abs(s.get('approach_increment',0.)) for s in selected),default=0.),
                max_k_over_kgood=max(s['k']/s['kgood'] for s in selected if s['kgood']>0),
                ff_upper_band=[1.23,1.29],around_peak=around)


def episode(case, run='own', full=False):
    O.setup(); sc=MC.scene_for(tuple(case['cell']),case['base'],case['variant'])
    pol=B.build('Huber-D4',O.chosen()['Huber-D4'],sc)
    windows=[]; commands=[]; steps=[]; contexts={}
    attach_instrumentation(pol,run,windows,commands,steps,contexts)
    log=P3.episode(sc,pol,case['seed'],diagnostic=True)
    assert log is not None
    metric=M.metric(log,pol,sc); peak=float(metric['peak'])
    if run=='own' and 'stored_peak' in case:
        assert abs(peak-case['stored_peak'])<(.05 if case['rounded'] else 1e-7),(case,peak)
    j=int(np.argmax(np.where(log[:,0]>=1.,log[:,18],-np.inf))); tp=float(log[j,0])
    finish_windows(windows,log,steps,pol.estimator.nw)
    ev=event_summary(windows,steps,log,tp)
    result=dict(case=case['label'],run=run,peak=peak,L=float(np.log1p(metric['J'])),
                peak_t=tp,event=ev,response=response(steps,log,ev['first_upper_exit'],tp),
                downward_jumps=[dict(t=r['t'],prior=r['prior_kgood'],after=r['kgood_after'],
                                     sigma_units=(r['prior_kgood']-r['kgood_after'])/r['sigma'])
                                for r in windows if r['accepted'] and r['kgood_after']<r['prior_kgood']])
    if full: return result,sc,log,windows,commands,steps,contexts
    return result


def clean(x):
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,(np.bool_,)): return bool(x)
    if isinstance(x,(np.integer,)): return int(x)
    if isinstance(x,(float,np.floating)): return float(x) if np.isfinite(x) else None
    return x


STAT_FIELDS=('k','k_hat','raw_k','sigma','accepted','prior_kgood','kgood_after',
             'probe_force_amplitude','YMY','t','ln_k_span')


def stats_job(case):
    result,sc,log,windows,commands,steps,contexts=episode(case,full=True)
    rows=[[r.get(k,np.nan) for k in STAT_FIELDS] for r in windows if r['k']>0 and np.isfinite(r['k_hat']) and np.isfinite(r['sigma'])]
    return dict(case=case,peak=result['peak'],L=result['L'],rows=rows)


def distribution(a):
    a=np.asarray(a,float);a=a[np.isfinite(a)]
    return dict(n=len(a),min=float(a.min()),max=float(a.max()),quantiles=dict(zip(('p01','p05','p25','p50','p75','p95','p99'),
                np.quantile(a,[.01,.05,.25,.5,.75,.95,.99]).tolist()))) if len(a) else dict(n=0)


def population_statistics(scenes):
    a=np.array([r for s in scenes for r in s['rows']],float)
    col={k:a[:,i] for i,k in enumerate(STAT_FIELDS)}
    k=col['k'];sig=col['sigma'];kh=col['k_hat'];raw=col['raw_k'];rel=sig/k
    accepted=col['accepted'].astype(bool);prior=col['prior_kgood'];after=col['kgood_after']
    bins=[0,.05,.1,.15,.2,.3,.5,1.,math.inf];rows=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        q=(rel>=lo)&(rel<hi);n=int(q.sum())
        if not n:continue
        c=float(np.median(kh[q]/k[q]));low=kh<.5*c*k;ql=q&low
        jump=q&accepted&(after<prior)&np.isfinite(prior)
        rows.append(dict(rel_sigma_bin=[lo,hi],n=n,c=c,
                         centered=distribution((kh[q]-c*k[q])/sig[q]),
                         raw_centered=distribution((raw[q]-c*k[q])/sig[q]),
                         low_count=int(ql.sum()),low_fraction=float(low[q].mean()),
                         gaussian_prediction=float(ndtr(-.5*c*k[q]/sig[q]).mean()),
                         low_accepted_count=int((ql&accepted).sum()),
                         low_accept_fraction=float(accepted[ql].mean()) if ql.any() else None,
                         accepted_downward_jumps=distribution((prior[jump]-after[jump])/sig[jump]),
                         accepted_low_downward_jumps=distribution((prior[ql&jump]-after[ql&jump])/sig[ql&jump])))
    def conditional(x,edges):
        rr=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            q=(x>=lo)&(x<hi)
            if q.any():rr.append(dict(bin=[lo,hi],relative_sigma=distribution(rel[q]),
                                       probe_amplitude=distribution(col['probe_force_amplitude'][q]),
                                       k=distribution(k[q])))
        return rr
    ratio=k/prior
    by_ratio=conditional(ratio,[0,.25,.5,.75,1.,1.3,2.,math.inf])
    by_after_ratio=conditional(k/after,[0,.25,.5,.75,1.,1.3,2.,math.inf])
    by_k=conditional(k,[0,2000,5000,10000,15000,25000,50000,math.inf])
    joint=[]
    for lo,hi in zip([0,5000,10000,15000,25000],[5000,10000,15000,25000,math.inf]):
        for rlo,rhi in zip([0,.5,.75,1.,1.3],[.5,.75,1.,1.3,math.inf]):
            q=(k>=lo)&(k<hi)&(ratio>=rlo)&(ratio<rhi)
            if q.sum():joint.append(dict(k_bin=[lo,hi],ratio_bin=[rlo,rhi],n=int(q.sum()),
                                        relative_sigma=distribution(rel[q])))
    # Within-scene, within-true-k-bin association controls the principal hardness confound.
    out=dict(windows=len(a),relative_sigma_bins=rows,by_k_over_prior_kgood=by_ratio,by_k_over_updated_kgood=by_after_ratio,
             by_k=by_k,joint_k_and_ratio=joint,
             conventions=dict(all_windows=True,includes_warmup=True,zero_truth_excluded=True,
                              overlapping_windows_not_independent=True,
                              invalid_windows_excluded=True,
                              center='median k_hat/k separately in each relative sigma bin',
                              gaussian='mean Phi(-0.5*c*k/sigma) over the SAME windows',
                              accepted='actual measurement assignment, including ModelPolicy centring',
                              sigma='existing calibrated frozen-IRLS surrogate, not a fitted standard deviation'),
             fields=STAT_FIELDS,scene_samples=scenes)
    return out


def pool_map(fn,jobs,workers):
    with ProcessPoolExecutor(max_workers=workers) as pool:
        out=[]
        for r in pool.map(fn,jobs):
            out.append(r)
            print('completed',fn.__name__,len(out),'/',len(jobs),flush=True)
    return out


def intervention_job(job):
    case,run,reference=job
    result,sc,log,windows,commands,steps,contexts=episode(case,run,full=True)
    e=reference['event']['first_upper_exit']
    if e:
        n=e['n'];s=next((s for s in steps if s['n']==n),None)
        crossing=event_summary(windows,steps,log,float(log[-1,0]))
        abnormal=crossing['first_abnormal_t']
        exits=crossing['exits']
        recurrence=dict(reference_t=e['t'],reference_n=n,
                        before_any_abnormal=any(x['direction']=='underestimate' and (abnormal is None or x['t']<abnormal) for x in exits),
                        first_upper_anywhere=crossing['first_upper_exit'],
                        at_reference=None if s is None else dict(t=s['t'],mode=s['mode'],k_over_kgood=s['k']/s['kgood'],k_over_kbar=s['k']/s['kbar'],
                                                                upper_outside=s['k']/s['kgood']>1.3 or s['k']/s['kbar']>1.3),
                        upper_exit_near_reference=any(x['direction']=='underestimate' and abs(x['t']-e['t'])<=.1 for x in exits),
                        definition='Same-clock upper-band status and any upper departure within +/-0.1 s; not a claim of identical endogenous trajectories.')
        result['event_recurrence']=recurrence
    return result


def save(data):
    OUT.write_text(json.dumps(clean(data),ensure_ascii=False,indent=2),encoding='utf-8')


def recorded_estimates(sc,log,contexts,targets,corrected=False):
    """Estimator replay with the original causal calibration/phase schedule.

    No gates/control feedback are rerun: only the factor's measurement path is
    changed. The H preprocessing state is rebuilt chronologically, NOT reset
    per window. For the plain replay every window must reproduce.
    """
    pol=B.build('Huber-D4',O.chosen()['Huber-D4'],sc)
    pl=dict(P3.PLATFORM)
    if 'tau_a_hat' in sc:pl['tau_a']=sc['tau_a_hat']
    pol.reset(pl)
    sensing=pol.step if isinstance(pol.step,M.ForceSensing) else None
    estimates={}; corrected_log=log.copy()
    def read(o):
        n=o.n
        # The live code calls estimate only at the recorded contact ticks.
        if n in contexts:
            ctx=contexts[n]
            for k in ('tau_h','mag_h','ndel','eb','ns','afc'):
                setattr(pol,k,ctx[k])
            kh,sigma,chi=pol.estimator.estimate(o,n,pol)
            estimates[n]=dict(k_hat=float(kh),sigma=float(sigma),chi=float(chi))
            if n in targets and np.isfinite(kh):
                est=pol.estimator;ll=est.private
                j0=max(pol.ndel+1,n-1-est.nw);m=n-1-j0
                _,Y,z,keep=common._build(ll,j0,m,0,pol.ns,pol.eb,pol.afc,pol.ndel,ll[n-1,0],ll[n-1,2],est.nw*pol.ns*1e-3,pol.params['enc_y'])
                broken=np.flatnonzero(~keep)
                if len(broken):Y,z=Y[broken[-1]+1:],z[broken[-1]+1:]
                estimates[n]['Y']=Y.tolist();estimates[n]['z']=z.tolist()
        corrected_log[:n]=o.log[:n]
        return 0.
    if sensing:sensing.step=read
    for j in range(max(targets)):
        n=j+1
        # Set before sensing as well: its command-to-acceleration model uses ndel.
        if n in contexts:
            for k in ('tau_h','mag_h','ndel','eb','ns','afc'):
                setattr(pol,k,contexts[n][k])
        o=P3.Obs(log,n)
        if sensing and not corrected:sensing(o)
        else:read(o)
    return {n:estimates[n] for n in targets},corrected_log,pol


_DIAG_TICK=None


def diagnostic_tick():
    global _DIAG_TICK
    if _DIAG_TICK is not None:return _DIAG_TICK
    fn=P3._tick.py_func
    def transform(f):
        class Freeze(ast.NodeTransformer):
            def visit_Call(self,node):
                if isinstance(node.func,ast.Name) and node.func.id=='_contact':
                    node.args[0]=ast.parse('c[9] if c[10] else xc',mode='eval').body
                return self.generic_visit(node)
        Freeze().visit(f)
        for node in ast.walk(f):
            if isinstance(node,ast.If) and ast.unparse(node.test)=='j % int(round(c[32] / dt)) == 0':
                idx=next(i for i,x in enumerate(node.body) if isinstance(x,ast.AugAssign) and ast.unparse(x.target)=='nlog')
                for off,expr in reversed(list(enumerate(('xc','s[24]','s[2]','s[3]'),26))):
                    node.body.insert(idx,ast.parse(f'log[nlog, {off}] = {expr}').body[0])
    _DIAG_TICK=njit(cache=False)(private_function(fn,transform))
    return _DIAG_TICK


def local_replay(sc,commands,seed,tick,highrate=False,sensor=True,frozen_x=None,freeze_from_tick=0):
    """P3.replay's numerical path, with private extra diagnostic columns."""
    pl=dict(P3.PLATFORM)
    if 'tau_a_true' in sc:pl['tau_a']=sc['tau_a_true']
    if highrate:pl['log_dt']=.0005
    s,buf=P3._init(sc,pl);c=P3._constants(sc,pl)
    if frozen_x is not None:c[9]=float(frozen_x)
    a,tau,_=P3.memory_arrays(sc);geom=P3._geometry(sc)
    noise=P3._noise(np.random.default_rng(seed),int(sc['T']*1000)+10,sensor)
    stride=int(round(pl['Tc']/pl.get('log_dt',pl['Tc'])))
    log=np.zeros((len(commands)*stride+2,30));n=0;kk=0
    for index,command in enumerate(commands):
        if frozen_x is not None:c[10]=float(index>=freeze_from_tick)
        seen=s[5];n0=n
        if not np.isfinite(command):command=buf[-1]
        n,kk=tick(s,buf,float(command),c,geom,sc['q_field'],sc['E_field'],sc['ribs'],sc['h_s'],
                  sc['xcmd'],sc['u'],sc['fext'],sc['act'],sc['rip'],a,tau,noise,kk,log,n)
        assert n>=0
        log[n0,5]=seen
    return log[:n]


def motion_gauge(sc,hi,raw,pol,source='u'):
    """Depth-preserving offline coordinate replay, including probe inertia.

    w is the exact linear support response to the removed u, in the same
    semi-implicit 0.5 ms scheme as _tick. Preserve d by solving
          delta_y - (tool_mass/kr) delta_a = w.
    Contact/lateral forces and depth are then identical. This is a diagnostic
    coordinate replay, not a claim that its inverse-actuator input is held
    constant over each original 8 ms command interval.
    """
    dt=.0005;n=len(hi);w=np.zeros(n);pos=float(sc['u'][0]);vel=0.
    for j,t in enumerate(hi[:,0]):
        it=min(int(round(t/dt))//2,len(sc['u'])-2);fr=(t-it*.001)/.001
        uu=sc['u'][it]*(1-fr)+sc['u'][it+1]*fr
        ud=(sc['u'][it+1]-sc['u'][it])*1000
        vel+=dt*(sc['k']*(uu-pos)+sc['c']*(ud-vel))/sc['m'];pos+=dt*vel
        w[j]=pos
    w-=float(sc['u'][0])
    body_initial=float(P3._init(sc,dict(P3.PLATFORM))[0][3])
    surface=np.interp(hi[:,26],sc['field_x0']+np.arange(len(sc['h_s']))*sc['dx'],sc['h_s'])
    surface-=float(np.interp(P3._init(sc,dict(P3.PLATFORM))[0][7],
                           sc['field_x0']+np.arange(len(sc['h_s']))*sc['dx'],sc['h_s']))
    if source=='surface':w=surface
    elif source=='body_load':w=hi[:,29]-body_initial-w
    elif source=='body_all':w=hi[:,29]-body_initial
    elif source=='all':w=hi[:,29]-body_initial+surface
    else:assert source=='u'
    mass=sc.get('tool_mass',0.);lam=mass/sc['kr']/dt**2
    ab=np.zeros((3,n));ab[0,1:]=-lam;ab[1]=1+2*lam;ab[2,:-1]=-lam
    ab[1,0]=ab[1,-1]=1.;ab[0,1]=0.;ab[2,-2]=0.
    dy=solve_banded((1,1),ab,w)
    da=np.zeros(n);da[1:-1]=(dy[2:]-2*dy[1:-1]+dy[:-2])/dt**2
    depth_error=(dy-w-mass/sc['kr']*da)*np.cos(hi[:,17])
    # Remove u but transform the axis, body and raw sensor together.
    # The true and calibrated mass/delay streams remain distinct.
    dv=np.r_[0.,np.diff(dy)/dt];tau=sc.get('tau_a_true',P3.PLATFORM['tau_a'])
    inverse=dy+tau*dv+(tau-P3.ACTUATOR_T2)*P3.ACTUATOR_T2*da
    hat_tau=sc.get('tau_a_hat',P3.PLATFORM['tau_a']);hat_t2=P3.ACTUATOR_T2
    hat_y=float(dy[0]);hat_v=float(dv[0]);hat_a=np.zeros(n)
    for j in range(n):
        acc=(inverse[j]-hat_y-hat_tau*hat_v)/((hat_tau-hat_t2)*hat_t2)
        hat_a[j]=acc;hat_v+=dt*acc;hat_y+=dt*hat_v
    def filtered(mass,delay,afc,acceleration):
        state=0.;hist=np.zeros(4);out=np.zeros(n);lag=min(int(delay),2);frac=delay-lag
        for j in range(n):
            before=state
            if j%2==0:
                hist[1:]=hist[:-1];hist[0]=mass*acceleration[j]
                state+=afc*((1-frac)*hist[lag]+frac*hist[lag+1]-state)
            out[j]=before if j%16==0 else state
        return out
    actual_delay=sc.get('sensor_delay_ms',0.);hat_delay=sc.get('sensor_delay_hat',0.)
    afc=pol.afc
    d_raw=filtered(mass,actual_delay,afc,da)
    d_cal=filtered(sc.get('tool_mass_hat',0.),hat_delay,afc,hat_a)
    out=raw.copy();ix=np.arange(len(raw))*16
    # raw here is already the adopted, calibrated sensing output.
    shift=hat_delay*.001/.008
    yshift=dy[ix].copy();yshift[1:]-=shift*(dy[ix[1:]]-dy[ix[:-1]])
    out[:,2]+=yshift;out[:,5]+=(d_raw-d_cal)[ix]
    out[:,4]-=w[ix]
    # A command change is irrelevant to _build's selected enc_y=True branch;
    # retain its inverse-actuator diagnostic value, without using it as a
    # deployable held-command claim.
    out[:,11]+=inverse[ix]
    return out,dict(source=source,max_abs_depth_error_m=float(np.max(np.abs(depth_error))),
                    max_axis_shift_m=float(np.max(np.abs(dy))),
                    max_axis_shift_minus_support_m=float(np.max(np.abs(dy-w))),
                    max_raw_sensor_shift_N=float(np.max(np.abs(d_raw))),
                    max_calibration_residual_shift_N=float(np.max(np.abs(d_raw-d_cal))),
                    equivalent='offline depth-preserving coordinate replay; full force and geometry unchanged; original calibration/phase schedule',
                    boundary='zero extra acceleration at both ends; inverse actuator input need not be held at 8 ms')


def array_estimate(Y,z,pol):
    Y=np.asarray(Y);z=np.asarray(z);est=pol.estimator;m=len(Y)
    D,basis=est._difference_bands(m);rhs=np.column_stack((Y,z));w=np.ones(m-est.order)
    scale=pol.params['noise']/math.sqrt(est.rho)
    for it in range(est.iters):
        ab=(basis@w).reshape(est.order+1,m);ab[0]+=1
        smooth=solveh_banded(ab,rhs,lower=True,overwrite_ab=True,check_finite=False)
        sy,sz=smooth[:,0],smooth[:,1];my=Y-sy;info=float(Y@my);k=float(my@z)/info
        if it==est.iters-1:break
        u=np.abs(D@(sz-k*sy))/scale;w=np.minimum(1.,est.c/np.maximum(u,1e-12))
    res=z-k*Y
    return dict(k_hat=float(np.clip(k,200,3e5)),raw_k=k,
                sigma=est.sigma_scale*pol.params['noise']/math.sqrt(info),
                chi=est.chi_scale*float(res@((z-sz)-k*my))/pol.params['noise']**2/max(m-est.order-1,1))


def target_windows(result,windows,log):
    targets={}
    ev=result['event'];first=ev['first_upper_exit']
    if first:targets[first['window_n']]=['first_upper_exit']
    if ev['first_either_exit']:
        n=ev['first_either_exit']['n'];r=next((r for r in windows if r['n']==n),None)
        if r:targets.setdefault(n,[]).append('first_either_exit')
    # Peak precursor: the largest accepted downward change in the preceding
    # second, including windows BEFORE the crossing. This selection is stated,
    # not silently labelled the first episode-wide crossing.
    tp=result['peak_t'];near=[r for r in windows if tp-1.<=r['t']<=tp and r['accepted'] and
                            np.isfinite(r['prior_kgood']) and r['kgood_after']<r['prior_kgood']]
    if near:
        r=min(near,key=lambda r:r['kgood_after']/r['prior_kgood'])
        targets.setdefault(r['n'],[]).append('largest_fractional_drop_in_prepeak_second')
    # Keep the last pre-peak upper exit's first member within that final second.
    ee=[e for e in ev['exits'] if e['direction']=='underestimate' and tp-1.<=e['t']<=tp]
    if ee:targets.setdefault(ee[0]['n'],[]).append('first_upper_exit_in_prepeak_second')
    return targets


def forensics_job(case):
    result,sc,log,windows,commands,steps,contexts=episode(case,full=True)
    target_map=target_windows(result,windows,log);targets=sorted(target_map)
    nt=max(targets);cmd=np.array(commands[:nt])
    plain=P3.replay(sc,cmd,case['seed'],sensor=True)
    assert plain is not None
    plain_error=max_difference(plain,log[:nt])
    assert plain_error<1e-10,plain_error
    repro,corrected,pol=recorded_estimates(sc,plain,contexts,targets)
    wm={r['n']:r for r in windows}
    for n in targets:
        assert abs(repro[n]['k_hat']-wm[n]['k_hat'])<1e-7,(case['label'],n,repro[n],wm[n])
    tick=diagnostic_tick();hi=local_replay(sc,cmd,case['seed'],tick,highrate=True)
    assert max_difference(hi[::16,:26],plain)<1e-10
    gauges={}
    for source in ('u','surface','body_load','body_all','all'):
        gauge,gauge_check=motion_gauge(sc,hi,corrected,pol,source)
        gauge_est,_,_=recorded_estimates(sc,gauge,contexts,targets,corrected=True)
        gauges[source]=(gauge_est,gauge_check)
    variants={}
    for removal in ('noise_off','mass_off','sampling_offset_off','event_off'):
        ss=copy.deepcopy(sc)
        if removal=='mass_off':ss['tool_mass']=ss['tool_mass_hat']=0.
        if removal=='sampling_offset_off':ss['sensor_delay_ms']=ss['sensor_delay_hat']=0.
        if removal=='event_off':
            if not ss.get('event_info'):continue
            ss['u']=ss['u']-sum(ss['motion_components'],np.zeros_like(ss['u']))
            ss['act']=np.zeros_like(ss['act']);ss['fext']=np.zeros_like(ss['fext'])
        ll=P3.replay(ss,cmd,case['seed'],sensor=removal!='noise_off')
        assert ll is not None
        est,_,_=recorded_estimates(ss,ll,contexts,targets)
        variants[removal]=(est,ll)
    # A literal u-added command replay is reported as a check, alongside the
    # depth-exact equivalent; it is NOT called depth-matched when it isn't.
    ss=copy.deepcopy(sc);ss['u']=np.zeros_like(ss['u'])
    delta=np.interp(np.arange(nt)*.008,np.arange(len(sc['u']))*.001,sc['u'])
    ll=P3.replay(ss,cmd+delta,case['seed'],sensor=True)
    literal_est,_,_=recorded_estimates(ss,ll,contexts,targets)
    result['reproduction']=dict(max_abs_log_difference=plain_error,
                               max_abs_khat_difference=max(abs(repro[n]['k_hat']-wm[n]['k_hat']) for n in targets))
    result['forensic_windows']=[]
    for n in targets:
        r=wm[n];index=windows.index(r);xf=float(hi[(n-1)*16-1,26])
        frozen=local_replay(sc,cmd[:n],case['seed'],tick,frozen_x=xf,freeze_from_tick=r['j0'])
        fest,_,_=recorded_estimates(sc,frozen[:,:26],contexts,[n])
        values=dict(plain=repro[n],noise_off=variants['noise_off'][0][n],
                    motion_off_depth_matched=gauges['u'][0][n],
                    spatial_frozen=fest[n],mass_off=variants['mass_off'][0][n],
                    sampling_offset_off=variants['sampling_offset_off'][0][n],
                    literal_u_added=literal_est[n])
        for source in ('surface','body_load','body_all','all'):
            values[source+'_off_depth_matched']=gauges[source][0][n]
        # Does the realised sensor-noise residual alone suffice at the ORIGINAL
        # regressor? This replaces all noiseless mismatch by k_end*Y; it is a
        # conditional regression diagnostic, not a new closed-loop scene.
        Y=np.asarray(repro[n]['Y']);z=np.asarray(repro[n]['z'])
        zn=np.asarray(variants['noise_off'][0][n]['z'])
        check=array_estimate(Y,z,pol)
        assert abs(check['k_hat']-r['k_hat'])<1e-7
        values['noise_only_fixed_regressor']=array_estimate(Y,r['k']*Y+(z-zn),pol)
        if 'event_off' in variants:values['event_off']=variants['event_off'][0][n]
        for key,val in values.items():
            val['ratio_to_original_truth']=val['k_hat']/r['k']
            val['restored_original_band']=.75<=r['k']/(CENTRE*val['k_hat'])<=1.3 if val['k_hat']>0 else False
        jj=slice(r['j0'],n-1)
        checks={key:dict(max_abs_depth_difference_m=float(np.max(abs(l[jj,6]-plain[jj,6]))),
                         truth_k_end=float(l[n-2,21])) for key,(_,l) in variants.items()}
        checks['spatial_frozen']=dict(frozen_contact_x=xf,truth_k_end=float(frozen[n-2,21]),
                                      true_ln_k_span=float(np.log(frozen[jj,21].max()/frozen[jj,21].min())) if frozen[jj,21].min()>0 else None,
                                      max_abs_depth_difference_m=float(np.max(abs(frozen[jj,6]-plain[jj,6]))))
        checks['literal_u_added']=dict(max_abs_depth_difference_m=float(np.max(abs(ll[jj,6]-plain[jj,6]))))
        row=dict(reasons=target_map[n],n=n,window=r,history=windows[max(0,index-10):index+1],
                 counterfactuals=values,checks=checks,motion_depth_check={k:v[1] for k,v in gauges.items()},
                 response=response(steps,log,dict(t=r['t']),result['peak_t']))
        result['forensic_windows'].append(row)
    result['all_windows']=windows
    result['limitations']=[
        'Single realised noise sequence and original endogenous commands: noise-off establishes realised contribution, not probability or noise alone sufficiency.',
        'Phase and calibration schedule frozen to the own run; the estimator H preprocessing and IRLS weights are recomputed.',
        'Spatial freeze holds the material lookup position at window end, not the load-dependent true tangent; changed force also changes depth and history.',
        'Motion depth match is an offline coordinate replay with exact depth equation; its inverse actuator input is not constrained to 8 ms holds. Literal command+u depth error is separately measured.',
        'Mass removal changes both sensor inertia and series-compliance depth; sampling-offset removal does not remove encoder intra-tick interpolation error.',
        'One-factor removals cannot uniquely allocate interactions or distinguish a probe-band residual aligned with Y from a gain change.',
        'A band exit alone does not establish a destabilizing cause; both early harmless exits and later peak precursors are retained.'
    ]
    return result


def max_difference(a,b):
    finite=np.isfinite(a)&np.isfinite(b)
    assert np.array_equal(np.isfinite(a),np.isfinite(b))
    assert np.array_equal(a[~finite],b[~finite],equal_nan=True)
    return float(np.max(abs(a[finite]-b[finite])))


def preprocessing_tests(n,corrected,contexts,pol,original):
    """Remove individual parts of the sampled H coordinate transform only.

    The adopted transform applies changing mag_h to (y-y_at_episode_start),
    caching each past row permanently. These are conditional estimator tests,
    not policy interventions or proposed changes to the Method.
    """
    ctx=contexts[n];keys=sorted(contexts);y=corrected[:n,2].copy();dy=np.r_[0.,np.diff(y)]
    row_ctx=[contexts[keys[min(int(np.searchsorted(keys,j+1)),len(keys)-1)]] for j in range(n)]
    tau=np.array([c['tau_h'] for c in row_ctx]);mag=np.array([c['mag_h'] for c in row_ctx])
    variants=dict(original_rebuilt=y[0]+((y-y[0])+tau*dy/.008)/mag,
                  H_off=y,
                  H_mag_fixed=y[0]+((y-y[0])+tau*dy/.008)/ctx['mag_h'],
                  H_tau_fixed=y[0]+((y-y[0])+ctx['tau_h']*dy/.008)/mag,
                  H_both_fixed=y[0]+((y-y[0])+ctx['tau_h']*dy/.008)/ctx['mag_h'],
                  H_reference_at_window_end=y[-1]+((y-y[-1])+tau*dy/.008)/mag)
    out={};j0=max(ctx['ndel']+1,n-1-pol.estimator.nw);m=n-1-j0
    for name,yp in variants.items():
        ll=corrected[:n].copy();ll[:,2]=yp
        _,Y,z,keep=common._build(ll,j0,m,0,ctx['ns'],ctx['eb'],ctx['afc'],ctx['ndel'],ll[-1,0],ll[-1,2],pol.estimator.nw*ctx['ns']*.001,pol.params['enc_y'])
        broken=np.flatnonzero(~keep)
        if len(broken):Y,z=Y[broken[-1]+1:],z[broken[-1]+1:]
        out[name]=array_estimate(Y,z,pol)
    assert abs(out['original_rebuilt']['k_hat']-original['k_hat'])<1e-6,(n,out['original_rebuilt'],original)
    return out,dict(y_minus_episode_origin_m=float(y[-1]-y[0]),
                    tau_h_range=[float(tau[j0:].min()),float(tau[j0:].max())],
                    mag_h_range=[float(mag[j0:].min()),float(mag[j0:].max())])


def audit_job(job):
    case,prior=job
    current,sc,log,windows,commands,steps,contexts=episode(case,full=True)
    targets=[r['n'] for r in prior['forensic_windows']];nt=max(targets)
    plain=P3.replay(sc,commands[:nt],case['seed'],sensor=True)
    estimates,corrected,pol=recorded_estimates(sc,plain,contexts,targets)
    wm={r['n']:r for r in windows}
    prior['event']=current['event'];prior['response']=current['response'];prior['all_windows']=windows
    for row in prior['forensic_windows']:
        n=row['n'];r=wm[n];index=windows.index(r)
        extra,diagnostics=preprocessing_tests(n,corrected,contexts,pol,estimates[n])
        row['counterfactuals'].update(extra);row['preprocessing_diagnostics']=diagnostics
        row['window']=r;row['history']=windows[max(0,index-10):index+1]
        row['response']=response(steps,log,dict(t=r['t']),current['peak_t'])
    prior['response_from_first_abnormal']=response(steps,log,dict(t=current['event']['first_abnormal_t']),current['peak_t']) if current['event']['first_abnormal_t'] is not None else None
    prior['limitations'].append('H preprocessing tests change only its cached encoder-coordinate transform, conditional on original commands, phase and calibration; H_off also removes contact-model phase compensation.')
    return prior


def validate_results(data):
    """Cross-check all delivered numbers, without rerunning completed episodes."""
    own={r['case']:r for r in data['own']}
    assert len(own)==7 and len(data['interventions'])==28
    for case in data['cases']:
        r=own[case['label']]
        assert abs(r['peak']-case['stored_peak'])<(.05 if case['rounded'] else 1e-7)
        assert r['reproduction']['max_abs_log_difference']==0.
        assert r['reproduction']['max_abs_khat_difference']==0.
        for w in r['forensic_windows']:
            assert len(w['history'])==11,(r['case'],w['n'],len(w['history']))
            assert w['counterfactuals']['plain']['k_hat']==w['window']['k_hat']
            assert abs(w['counterfactuals']['original_rebuilt']['k_hat']-w['window']['k_hat'])<1e-6
            assert w['checks']['noise_off']['max_abs_depth_difference_m']==0.
            for check in w['motion_depth_check'].values():assert check['max_abs_depth_error_m']<1e-12
            k=w['window']['k']
            for name,cf in w['counterfactuals'].items():
                truth=w['checks'].get(name,{}).get('truth_k_end',k)
                if truth is not None and truth>0 and cf.get('k_hat') is not None:
                    cf['counterfactual_truth_k']=truth
                    cf['ratio_to_counterfactual_truth']=cf['k_hat']/truth
    assert len(data['statistics']['scene_samples'])==72
    assert {v:sum(s['case']['variant']==v for s in data['statistics']['scene_samples']) for v in ('nominal','slow_thin')}==dict(nominal=36,slow_thin=36)
    scenes=data['statistics']['scene_samples']
    for scene in scenes:
        c=scene['case'];assert c['base']==520000
        grid=[cell for cell in O.cells() if cell[0]=='breath_rib']
        assert c['seed']==MC.SEED+grid.index(tuple(c['cell']))
    # Recompute delivered aggregates from saved primary window samples.
    data['statistics']=population_statistics(scenes)
    by_variant={}
    for variant in ('nominal','slow_thin'):
        s=population_statistics([r for r in scenes if r['case']['variant']==variant]);s.pop('scene_samples')
        by_variant[variant]=s
    data['statistics']['by_variant']=by_variant
    for r in data['interventions']:
        assert 'event_recurrence' in r
        assert r['event_recurrence']['reference_n']==own[r['case']]['event']['first_upper_exit']['n']
    data['definitions']=dict(
        peak_and_L='method3.metric: true axial force peak after t>=1 s, L=log1p(J); startup is included in event search and window statistics',
        events='First departure after a channel has first entered [0.75,1.3]. All departures, both directions, and the first upper departure are retained; no crossing is automatically called destabilizing.',
        low_window_targets='First either-direction exit, first upper exit, first upper exit in the prepeak second, and largest fractional accepted downward update in that second; no force-peak causality inferred solely from order.',
        spatial_freeze='Material lookup location fixed at the window-end contact location, only within the target regression window, with original prefix state. Load/activation dependence of true tangent is retained.',
        gate_prior_free='Only the precision comparison uses sigma<=rel_max*k_hat. The original floor guard and prior-based measurement bound are retained.',
        kgood_limit='Limit the final, centred k_good decrease to (1-rho_d) times prior k_good per estimator call; rising is unchanged.',
        exc_from_kbar='Only the excitation denominator changes after k_good is finite. Initial k_seek acquisition, feed-forward, and feedback remain the adopted paths.',
        both='kgood_limit plus exc_from_kbar, and no other policy changes',
        probe_amplitude='a_force*exc_gain*k/k_un; current, prior, and in-window amplitudes are recorded separately for each forensic window.',
        H_tests='Conditional regression re-estimation at original commands, noise, calibration and phase schedule. H_mag_fixed changes only mag_h to its window-end value; original_rebuilt verifies the sampled transform exactly.',
        motion_check='Depth equation is checked on every 0.5 ms substep, with mass and the calibrated sensing path transformed. Literal command+u replay depth errors are separately retained; the depth-exact inverse command need not respect 8 ms holds.',
        null='Unavailable/unevaluated/invalid scalar, or unbounded histogram endpoint; gate path distinguishes no valid estimate from rejection.',
        resources='At most 6 simulation worker processes; BLAS/Numba thread limits 1; no external services or hardware; no git mutations')
    data['validation']=dict(passed=True,own_peaks=7,replay_logs_exact=7,window_estimates_exact=7,
                            original_H_transform_checks=sum(len(r['forensic_windows']) for r in own.values()),
                            histories_have_10_predecessor_attempts=True,scene_count=72,intervention_count=28)
    print(json.dumps(data['validation']),flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6)
    ap.add_argument('--stage',choices=('cases','stats','interventions','forensics','audit','validate','all'),default='all');args=ap.parse_args()
    assert 1<=args.workers<=6
    O.setup(); cases=selected_cases()
    data=json.loads(OUT.read_text(encoding='utf-8')) if OUT.exists() else dict(cases=cases)
    if args.stage in ('cases','all'):
        own=pool_map(episode,cases,args.workers);data['own']=own;save(data)
        for r in own:
            print(json.dumps(clean({k:r[k] for k in ('case','run','peak','L','peak_t')})),flush=True)
    if args.stage in ('stats','all'):
        grid=[c for c in O.cells() if c[0]=='breath_rib']
        jobs=[dict(label=f'{variant}_{q}',cell=c,base=MC.DEV_BASE,variant=variant,seed=MC.SEED+q)
              for variant in ('nominal','slow_thin') for q,c in enumerate(grid)]
        assert len(jobs)==72
        scenes=pool_map(stats_job,jobs,args.workers)
        data['statistics']=population_statistics(scenes);save(data)
    if args.stage in ('interventions','all'):
        refs={r['case']:r for r in data['own']}
        jobs=[(case,run,refs[case['label']]) for run in ('kgood_limit','exc_from_kbar','gate_prior_free','both') for case in cases]
        data['interventions']=pool_map(intervention_job,jobs,args.workers);save(data)
    if args.stage in ('forensics','all'):
        data['own']=pool_map(forensics_job,cases,args.workers);save(data)
    if args.stage in ('audit','all'):
        refs={r['case']:r for r in data['own']}
        data['own']=pool_map(audit_job,[(c,refs[c['label']]) for c in cases],args.workers);save(data)
    if args.stage in ('validate','all'):
        validate_results(data);save(data)


if __name__=='__main__': main()
