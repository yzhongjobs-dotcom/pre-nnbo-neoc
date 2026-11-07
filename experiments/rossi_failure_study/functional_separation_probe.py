"""Construct observationally identical endpoint-gain alternatives.

An adversarial diagnostic, NOT an estimator or a physiological calibration.
The candidate gain correction is monotone, starts at zero, ends at +/-20%,
and has bounded first/second differences. The force-fit residual is held
EXACTLY unchanged; only the original disturbance-prior cost is compared.
"""
import json
import numpy as np
import cvxpy as cp
import structured_tvk_audit as A


def huber_cost(v,c=1.345):
    av=np.abs(v)
    return float(np.sum(np.where(av<=c,.5*v*v,c*av-.5*c*c)))


def main():
    rows=[]
    for cond in ('rest','rib','cough','combo'):
      for replica in (0,1):
        a=np.load(A.OUT/f'{cond}_{replica}.npz');log=a['log'];idx=a['idx']
        if cond=='rest':j=len(idx)//2
        else:
            spans=np.array([np.ptp(log[jj:jj+24,21])/log[ii,21] for ii,jj in zip(idx,a['j0'])])
            j=int(np.argmax(spans))  # truth chooses a diagnostic window only
        Y,z,G=a['Y'][j],a['z'][j],a['G'][j]
        sigma,rho=float(a['noise'][j]),float(a['rho'][j])
        K,_,_,d0=A.huber_fit(Y,z,G,sigma,rho)
        D=np.diff(np.eye(24),4,axis=0);D1=np.diff(np.eye(24),axis=0);D2=np.diff(np.eye(24),2,axis=0)
        B0=huber_cost(np.sqrt(rho)*D@d0/sigma)
        for sign in (-1.,1.):
            h=cp.Variable(24)
            dalternative=d0/sigma-cp.multiply(K*Y/sigma,h)
            objective=.5*cp.sum(cp.huber(np.sqrt(rho)*D@dalternative,1.345))
            constraints=[h[0]==0.,h[-1]==sign*.2,sign*(D1@h)>=0.,cp.abs(D1@h)<=.05,cp.abs(D2@h)<=.02]
            p=cp.Problem(cp.Minimize(objective),constraints)
            p.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
            if p.status!='optimal':raise RuntimeError(p.status)
            hh=h.value;dd=d0-K*Y*hh
            identity=float(np.max(np.abs(K*Y+d0-(K*(1.+hh)*Y+dd))))
            rows.append(dict(cond=cond,replica=replica,time_s=float(log[idx[j],0]),
                             sign=sign,baseline_huber_cost=B0,alternative_huber_cost=float(p.value),
                             cost_difference=float(p.value-B0),same_prediction_max_error_N=identity,
                             correction_fraction=hh.tolist(),Y=Y.tolist(),z=z.tolist(),
                             original_d=d0.tolist(),alternative_d=dd.tolist(),K=K,
                             caution='The rate/curvature limits define an illustrative class, not a validated human class.'))
    result=dict(protocol='Identical fit residual, no noise manipulation; monotone smooth endpoint alternatives on selected physical windows.',rows=rows)
    (A.OUT/'separation_witnesses.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps([{k:r[k] for k in ('cond','replica','sign','baseline_huber_cost','alternative_huber_cost','same_prediction_max_error_N')} for r in rows],indent=2))


if __name__=='__main__':main()
