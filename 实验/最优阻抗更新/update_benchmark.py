"""Same target and stopping residual: warm Newton vs chord/sensitivity.

Parameters change between solves; this is not a time-integrated robot trial.
Sensitivity construction, Hessian refresh and fallback are all charged.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import model as m

HERE=Path(__file__).resolve().parent
METHODS=['cold_Newton_analytic_start','warm_Newton','warm_chord','sensitivity_chord']
TOL=1e-6
ANCHOR_EVERY=8


def prepare(x,g,colored,sensed,sensitivity):
    p=m.Problem(x,colored,sensed);_,gr,h=p.evaluate(g,True)
    derivative=np.zeros((4,len(x)))
    if sensitivity:
        epsilon=1e-4
        for i in range(len(x)):
            shift=np.eye(len(x))[i]*epsilon
            plus=m.Problem(x+shift,colored,sensed);minus=m.Problem(x-shift,colored,sensed)
            derivative[:,i]=(plus.evaluate(g)[1]-minus.evaluate(g)[1])/(2*epsilon)
            p.calls+=plus.calls+minus.calls
            p.lyapunov_solves+=plus.lyapunov_solves+minus.lyapunov_solves
    active=((g<=m.LOW+1e-8)&(gr>0))|((g>=m.HIGH-1e-8)&(gr<0))
    free=~active;gain_derivative=np.zeros_like(derivative)
    if sensitivity and free.any():
        gain_derivative[free]=-np.linalg.solve(h[np.ix_(free,free)],derivative[free])
    return h,gain_derivative,p.calls,p.lyapunov_solves


def path(rng,dim,kind):
    centers=rng.uniform(-.5,.5,dim);amplitudes=rng.uniform(.05,.25,dim)
    phases=rng.uniform(-np.pi,np.pi,dim)
    phase=np.linspace(0,2*np.pi,33)
    xx=centers+amplitudes*np.sin(phase[:,None]+phases)
    if kind=='with_model_revision':
        xx[12:]+=rng.uniform(-.35,.35,dim)
        xx[24:]+=rng.uniform(-.35,.35,dim)
    return np.clip(xx,-1,1)


def main():
    protocol=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
        seed=3901,paths_per_family=12,updates_per_path=32,anchor_every=ANCHOR_EVERY,
        stopping_residual=TOL,methods=METHODS,
        sequence_kinds=['smooth','with_model_revision'],
        timing='Random method order per path; preparation+monitoring+line search+fallback charged; construction of common problem data excluded',
        scope='Parameterized stationary four-gain problem, NOT physical gain integration or B0 finite-root certification')
    (HERE/'update_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    records=[]
    for family_index,(label,colored,sensed) in enumerate([('white',False,False),('colored',True,False),('colored_observation_noise',True,True)]):
        rng=np.random.default_rng(3901+family_index);dim=4 if colored else 3
        for case in range(12):
            kind='smooth' if case<6 else 'with_model_revision'
            xx=path(rng,dim,kind)
            initial=m.Problem(xx[0],colored,sensed).solve().x
            refs=[]
            for x in xx[1:]:
                p=m.Problem(x,colored,sensed);refs.append(p.solve([p.analytic(),initial]).fun)
            for method in rng.permutation(len(METHODS)):
                previous=initial.copy();anchorx=xx[0].copy();anchorg=initial.copy()
                for i,x in enumerate(xx[1:]):
                    p=m.Problem(x,colored,sensed);start=time.perf_counter();prep_calls=prep_lyap=0
                    if method>=2 and i%ANCHOR_EVERY==0:
                        anchorx=xx[i].copy();anchorg=previous.copy()
                        h,derivative,prep_calls,prep_lyap=prepare(anchorx,anchorg,colored,sensed,method==3)
                    guess=p.analytic() if method==0 else previous.copy()
                    if method==3:guess=np.clip(anchorg+derivative@(x-anchorx),m.LOW,m.HIGH)
                    g,info=m.newton(p,guess,frozen_hessian=h if method>=2 else None,maxiter=15,tol=TOL)
                    fallback=False
                    if not info['converged']:
                        fallback=True
                        g,info=m.newton(p,g,maxiter=30,tol=TOL)
                    if not info['converged']:
                        g=p.solve([g,p.analytic()]).x
                        _,gr=p.evaluate(g);info['projected_gradient']=float(np.max(np.abs(m.projected_gradient(g,gr))))
                        info['converged']=info['projected_gradient']<=TOL
                    elapsed=time.perf_counter()-start
                    value=p.evaluate(g)[0] # evaluation for research is outside reported online time/counters
                    records.append(dict(family=label,case=case,kind=kind,update=i,method=METHODS[method],
                        time_s=elapsed,online_calls=p.calls-1+prep_calls,
                        lyapunov_solves=p.lyapunov_solves-2+prep_lyap,fallback=fallback,
                        relative_gap=float(value/refs[i]-1),J=value,reference_J=refs[i],**info))
                    previous=g
            print(label,case,flush=True)
    summary={}
    for family in ['white','colored','colored_observation_noise']:
        summary[family]={}
        for kind in ['smooth','with_model_revision']:
            summary[family][kind]={}
            for method in METHODS:
                rr=[r for r in records if r['family']==family and r['kind']==kind and r['method']==method]
                times=np.array([r['time_s'] for r in rr])
                summary[family][kind][method]=dict(mean_ms=float(times.mean()*1000),
                    p95_ms=float(np.quantile(times,.95)*1000),
                    mean_evaluations=float(np.mean([r['online_calls'] for r in rr])),
                    mean_lyapunov_solves=float(np.mean([r['lyapunov_solves'] for r in rr])),
                    max_projected_gradient=max(r['projected_gradient'] for r in rr),
                    max_relative_gap=max(r['relative_gap'] for r in rr),
                    fallback_count=sum(r['fallback'] for r in rr),
                    converged_count=sum(r['converged'] for r in rr),cases=len(rr))
    (HERE/'update_results.json').write_text(json.dumps(dict(protocol=protocol,summary=summary,updates=records),indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
