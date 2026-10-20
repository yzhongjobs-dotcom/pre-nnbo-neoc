"""Few structured first-order planning steps on unchanged transient instances.

Sequential linear programming is a mature alternative; this tests its actual
cost/quality, not a claim of a new optimizer or certified real-time deployment.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.optimize import linprog
import model as m
import transient as t
import adjoint as ad

HERE=Path(__file__).resolve().parent


def main():
    source=json.loads((HERE/'transient_results.json').read_text())
    raw=np.load(HERE/'transient_paired.npz');rows=[]
    for row in source['cases']:
        label,case=row['family'],row['case'];colored=label!='white';dim=4 if colored else 3;n=10 if colored else 4
        arr=raw[label+f'_{case}_initial'];g0=arr[2*dim:2*dim+4]
        mu=arr[2*dim+4:2*dim+4+n];p0=arr[-n*n:].reshape(n,n)
        problem=m.Problem(row['new_features'],colored,colored)
        p=ad.AdjointTransient(problem,g0,mu,p0,blocks=2)
        p.evaluate(np.zeros(8));p.cached=None # exclude JIT compilation, not optimization setup
        start=time.perf_counter()
        zero=np.zeros(8)
        analytic=t.ramp(g0,problem.analytic()).reshape(2,6,4).mean(axis=1).ravel()
        candidates=[zero,analytic]
        values=[float(p.run(z)[0][0]) for z in candidates]
        selected=int(np.argmin(values));z=candidates[selected].copy();value=values[selected]
        aub=np.vstack((p.gmap,-p.gmap))
        bub=np.r_[np.tile(m.HIGH-g0,t.SEGMENTS),np.tile(g0-m.LOW,t.SEGMENTS)]
        history=[]
        for iteration in range(5):
            history.append(dict(iteration=iteration,J=value,
                gap_to_12_block_reference=value/row['methods'][-1]['J']-1,
                elapsed_ms=(time.perf_counter()-start)*1000))
            if iteration==4:break
            _,grad=p.evaluate(z)
            bounds=list(zip(np.maximum(-1.,z-1.),np.minimum(1.,z+1.)))
            linear=linprog(grad,A_ub=aub,b_ub=bub,bounds=bounds,method='highs')
            if not linear.success:break
            direction=linear.x-z
            if grad@direction>=-1e-12:break
            accepted=False
            for fraction in [1.,.5,.25,.125,.0625]:
                trial=z+fraction*direction
                jt=float(p.run(trial)[0][0])
                if jt<=value+1e-4*fraction*(grad@direction):
                    z,value=trial,jt;accepted=True;break
            if not accepted:break
        rows.append(dict(family=label,case=case,initialization='zero' if selected==0 else 'analytic-rate block average',
            history=history,rate_blocks=z.tolist(),
            minimum_gain_margin=float(min((p.gmap@z-np.tile(m.LOW-g0,t.SEGMENTS)).min(),
                                         (np.tile(m.HIGH-g0,t.SEGMENTS)-p.gmap@z).min()))))
    out=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        adjoint_sha256=hashlib.sha256((HERE/'adjoint.py').read_bytes()).hexdigest(),
        protocol=dict(task='Same 8 transient cases; no cost, model, bounds or samples changed',
            blocks=2,maximum_iterations=4,trust_radius_in_normalized_rate=1.,
            initialization='Choose zero vs averaged analytic rate by their same-horizon cost; no stationary optimizer required',
            line_search=[1.,.5,.25,.125,.0625],
            timing='Includes initialization cost evaluations, LP/line search and gradients; excludes JIT and common mechanical data construction',
            status='Exploratory algorithm iteration after inspecting full solves, not a preregistered new test set'),cases=rows)
    (HERE/'bounded_results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
