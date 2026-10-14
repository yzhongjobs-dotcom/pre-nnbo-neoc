"""Numerical implementation checks of the same fixed common root.

No new contact scenario, no NEOC execution or candidate correction.
Run after run_experiment.py: python -X utf8 diagnose_anchor.py
"""
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import run_experiment as r

out=Path(__file__).resolve().parent
p=r.RootProblem(); eta=np.array([1.,0.])
w=np.load(out/'anchor.npz')['w']
res,L,u,lam,active=p.evaluate(w,eta,True)
h=1e-6
fd=np.column_stack([(p.evaluate(w+np.eye(28)[j]*h,eta)[0]-p.evaluate(w-np.eye(28)[j]*h,eta)[0])/(2*h) for j in range(28)])
f,G,gphi,ell,A,b=p.terms(eta)
grad=gphi+np.einsum('nwj,w->nj',p.px,w)
c=np.einsum('nji,nj->ni',G,grad)
station=2*r.CFG['action_weight']*u+c+np.einsum('nki,nk->ni',A,lam)
slack=b-np.einsum('nki,ni->nk',A,u)
delta=r.DREF*(1+p.y[:,0]); a=r.halfwidth(delta)
result=dict(scope='same fixed 28-weight/336-variable anchor, no extra simulation',
            jacobian_relative_error=float(np.linalg.norm(fd-L)/np.linalg.norm(L)),
            jacobian_condition=float(np.linalg.cond(L)),
            basis_rank=int(np.linalg.matrix_rank(p.ps)),
            node_stationarity_residual=float(np.max(abs(station))),
            complementarity_residual=float(np.max(abs(lam*slack))),
            minimum_slack=float(slack.min()),minimum_multiplier=float(lam.min()),
            geometry_inversion_residual=float(np.max(abs(r.psi(a)-delta))),
            anchor_bellman_residual=float(np.max(abs(res))))
assert result['jacobian_relative_error']<1e-6
assert result['node_stationarity_residual']<1e-8
assert result['complementarity_residual']<1e-8
assert result['minimum_slack']>=-1e-8 and result['minimum_multiplier']>=-1e-8
assert result['geometry_inversion_residual']<1e-12
def fun(z):return p.evaluate(z*p.weight_scale,eta)[0]*p.row_scale
def jac(z):return p.evaluate(z*p.weight_scale,eta,True)[1]*p.row_scale[:,None]*p.weight_scale
fits={}
for name,initial,limit in [('from_newton_iterate',w/p.weight_scale,300),('extended_same_iterate',None,3000),('from_zero',np.zeros(28),1000)]:
    if initial is None:initial=last.x
    st=time.perf_counter()
    fit=least_squares(fun,initial,jac=jac,xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=limit)
    ww=fit.x*p.weight_scale; residual=p.evaluate(ww,eta)[0]
    fits[name]=dict(root_tolerance_passed=bool(np.max(abs(residual))<r.CFG['newton_tolerance']),
                   max_residual=float(np.max(abs(residual))),nfev=int(fit.nfev),
                   message=fit.message,seconds=time.perf_counter()-st,
                   jacobian_condition=float(np.linalg.cond(p.evaluate(ww,eta,True)[1])))
    np.savez(out/('diagnostic_'+name+'.npz'),w=ww,eta=eta)
    last=fit
result['same_equation_solver_diagnostics']=fits
(out/'implementation_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
