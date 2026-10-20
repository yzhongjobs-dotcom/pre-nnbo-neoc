"""Strong-rival audit; test reduced-critic NEOC against structured chord etc.

Same minimized root and information. Closed-form solution is included since the
scalar problem permits it. Timings are diagnostic only, not deployment evidence.
"""
from pathlib import Path
import time,json
import numpy as np
import multinode as v

P=Path(__file__).resolve().parent
XS=v.XS;UB=v.UB
M=np.maximum(2*XS/(v.LO*v.LO+1)**1.5,2*XS**2*(XS**2+UB**2)/(v.LO*XS+UB)**3)
RADIUS=float(np.sqrt(2*.002/max(M)))

def F(w,th):
    u=np.clip(-w/2,-UB,UB)
    return XS**2-th*XS*w+u*u+w*u

def J(w,th):return -th*XS+np.clip(-w/2,-UB,UB)

def run(path,mode):
    anchor=path[0];wa,ca=v.exact_derivatives(anchor);w=wa.copy()
    chord_coefficient=XS*wa/J(wa,anchor)
    roots=1;root_steps=0;f_calls=0;event_updates=0;rows=[];references=[]
    switches=XS/(2*UB)-UB/(2*XS)
    for th in path:
        if mode=='closed_form':
            w,_=v.exact_derivatives(th);roots+=1
        elif mode=='newton_full':
            roots+=1
            for _ in range(12):
                residual=F(w,th);f_calls+=1
                if np.max(abs(residual))<1e-12:break
                w-=residual/J(w,th);root_steps+=1
            else:raise AssertionError('root did not converge')
        elif mode=='one_newton':
            w-=F(w,th)/J(w,th);f_calls+=1;root_steps+=1
        else:
            if mode=='segmented_value':
                crossed=switches[(switches>min(anchor,th)+1e-12)&(switches<max(anchor,th)-1e-12)]
                if len(crossed):
                    # Reinitialize at the last event passed along the direction.
                    anchor=float(max(crossed) if th>anchor else min(crossed))
                    wa,ca=v.exact_derivatives(anchor);event_updates+=len(crossed);roots+=1
            if abs(th-anchor)>RADIUS:
                anchor=th;wa,ca=v.exact_derivatives(anchor);roots+=1
                chord_coefficient=XS*wa/J(wa,anchor)
            if mode in ['neoc','segmented_value']:w=wa+ca*(th-anchor)
            elif mode=='structured_chord':
                # F(wa,th)-F(wa,anchor)=-XS*wa*(th-anchor), exactly.
                # Reuse the same derivative; no redundant Hamiltonian evaluations.
                w=wa+chord_coefficient*(th-anchor)
            else:raise ValueError(mode)
        rows.append(w.copy())
    return np.array(rows),dict(anchor_or_exact_roots=roots,Newton_steps=root_steps,
                              residual_evaluations=f_calls,event_updates=event_updates)

def main():
    modes=['closed_form','newton_full','one_newton','structured_chord','segmented_value','neoc']
    paths=dict(release=np.linspace(.9,1.2,121),activate=np.linspace(1.2,.9,121),
               shifted_release=np.linspace(.905,1.2,113))
    output={}
    for name,path in paths.items():
        exact=np.array([v.exact_derivatives(t)[0] for t in path]);out={};pred={}
        for mode in modes:
            values,counts=run(path,mode);pred[mode]=values
            timings=[]
            for _ in range(11):
                tic=time.perf_counter();run(path,mode);timings.append(time.perf_counter()-tic)
            error=float(np.max(abs(values-exact)))
            out[mode]=dict(maximum_gradient_error=error,within_common_budget=error<=.002,
                median_family_seconds=float(np.median(timings)),min_family_seconds=min(timings),
                max_family_seconds=max(timings),**counts)
        out['equivalence']=dict(neoc_structured_chord_max_abs=float(np.max(abs(pred['neoc']-pred['structured_chord']))),
            interpretation='For this parameter-affine reduced root, NEOC is algebraically the same structured cached-J chord predictor.')
        assert out['equivalence']['neoc_structured_chord_max_abs']<1e-12
        assert all(out[m]['within_common_budget'] for m in modes)
        output[name]=out
    output['decision']='Keep elimination as a regularity/conditioning repair; this test establishes neither a distinct updater nor NEOC necessity against strong competitors.'
    output['scope']='scalar constrained force-action family, not a reproduction of published ENE nor R15 gain-rate HJB; segmented_value is an event-aware diagnostic'
    (P/'strong_comparator_results.json').write_text(json.dumps(output,indent=2),encoding='utf8')
    print(json.dumps(output,indent=2),flush=True)

if __name__=='__main__':main()
