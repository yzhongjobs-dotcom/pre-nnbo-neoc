"""A direct, single-convex-objective TV-gain prototype on physical replay.

min_{x,d} .5||z/noise-diag(K*Y/noise)x-d||^2
          + sum huber_1.345(sqrt(rho) D4 d) + lam ||D1 x||_1.

K is a measured scale from the original estimate, frozen within each solve.
Unlike tvk_mhe.py, the forward gain column is the direct multiplicative model,
not the high-pass/mean operator. No scene labels or truth enter the solver.
This is a replay prototype; post-hoc output clipping matches existing practice.
"""
import json
import time
import argparse
import numpy as np
from numba import njit
import structured_tvk_audit as A


@njit(cache=True)
def admm(invH,rhs,Bd,Bk,x,d,lam,maxiter=3000,tol=1e-6):
    n=len(x);p=Bd@d;q=Bk@x;up=np.zeros(len(p));uq=np.zeros(len(q))
    primal=1e9;dual=1e9
    for it in range(maxiter):
        rr=rhs.copy();rr[:n]+=Bk.T@(q-uq);rr[n:]+=Bd.T@(p-up)
        theta=invH@rr;x=theta[:n];d=theta[n:]
        vd=Bd@d;vk=Bk@x
        p0=p.copy();q0=q.copy()
        for j in range(len(p)):
            v=vd[j]+up[j]
            p[j]=.5*v if abs(v)<=2.*1.345 else v-1.345*np.sign(v)
        for j in range(len(q)):
            v=vk[j]+uq[j];q[j]=np.sign(v)*max(abs(v)-lam,0.)
        up+=vd-p;uq+=vk-q
        primal=np.sqrt(np.sum((vd-p)**2)+np.sum((vk-q)**2))
        dual=np.sqrt(np.sum((Bd.T@(p-p0))**2)+np.sum((Bk.T@(q-q0))**2))
        dual_scale=np.sqrt(np.sum((Bd.T@up)**2)+np.sum((Bk.T@uq)**2))
        if primal<tol*(1.+np.sqrt(np.sum(vd*vd)+np.sum(vk*vk))) and dual<tol*(1.+dual_scale):
            break
    return x,d,it+1,primal,dual


def fit(Y,z,noise,rho,K,lam,maxiter=3000,tol=1e-6):
    n=len(Y);v=K*Y/noise;zz=z/noise
    Bk=np.diff(np.eye(n),axis=0);Bd=np.sqrt(rho)*np.diff(np.eye(n),4,axis=0)
    V=np.diag(v)
    H=np.block([[np.diag(v*v)+Bk.T@Bk,V],[V,np.eye(n)+Bd.T@Bd]])
    rhs=np.r_[v*zz,zz]
    # A linear solve is used once; ADMM then reuses this inverse. An optimized
    # deployment would use the banded factorization, not certify this timing.
    invH=np.linalg.solve(H,np.eye(2*n))
    x,d,it,pr,du=admm(invH,rhs,Bd,Bk,np.ones(n),zz-v,lam,maxiter,tol)
    return K*x[-1],dict(iterations=it,primal=pr,dual=du,trajectory=K*x)


def main():
    rows=[];lam_grid=(1.,3.,10.,30.,100.)
    timings=[];maxiter_count=0
    for replica in range(3):
      for cond in A.CONDITIONS:
        a=np.load(A.OUT/f'{cond}_{replica}.npz');truth=a['log'][a['idx'],21]
        take=np.arange(0,len(truth),3)  # targeted replay, not every overlapping window
        for lam in lam_grid:
            vals=[];iters=[]
            for j in take:
                t0=time.perf_counter()
                kh,diag=fit(a['Y'][j],a['z'][j],a['noise'][j],a['rho'][j],a['baseline'][j,0],lam)
                timings.append(time.perf_counter()-t0);iters.append(diag['iterations'])
                maxiter_count+=diag['iterations']==3000
                vals.append(np.clip(kh,200.,3e5))
            err=(np.array(vals)-truth[take])/truth[take]
            base=(a['baseline'][take,0]-truth[take])/truth[take]
            rows.append(dict(cond=cond,replica=replica,lam=lam,score=A.summary(err,np.ones(len(err),bool)),
                             baseline=A.summary(base,np.ones(len(base),bool)),
                             max_iterations=max(iters),mean_iterations=float(np.mean(iters))))
        print(f'TV replay {cond} {replica}',flush=True)
    dev=[r for r in rows if r['replica']==0]
    feasible=[]
    for lam in lam_grid:
        rr=[r for r in dev if r['lam']==lam]
        safe=all(r['score']['rrmse']<=1.05*r['baseline']['rrmse'] for r in rr if r['cond'] in ('rest','breath','shift_fast'))
        bad=[r['score']['rrmse']/r['baseline']['rrmse'] for r in rr if r['cond'] in ('rib','cough','combo')]
        if safe:feasible.append((float(np.mean(bad)),lam))
    chosen=min(feasible)[1] if feasible else None
    result=dict(protocol='Dev replica 0; replicas 1/2 are held-out for lambda selection. Comparisons use identical physical data.',
                selected_lambda=chosen,selection_feasible=feasible,rows=rows,
                timings_ms=(1000*np.percentile(timings,[50,99])).tolist(),
                maxiter_count=maxiter_count,n_solves=len(timings))
    result['caution']='Early ADMM diagnostic with a finite iteration budget; use direct_tv.json from the independent conic solve for scientific conclusions.'
    (A.OUT/'direct_tv_admm_diagnostic.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


def conic_main():
    """Independent converged solver, including an exact constant-k ablation.

    All trajectory comparisons are exploratory: some were inspected during
    measurement-model revision. No independent final-validation claim is made.
    """
    import cvxpy as cp
    n=24;Dk=np.diff(np.eye(n),axis=0)
    first=np.load(A.OUT/'rest_0.npz');rho=float(first['rho'][0])
    Bd=np.sqrt(rho)*np.diff(np.eye(n),4,axis=0)
    v=cp.Parameter(n);zz=cp.Parameter(n);lam=cp.Parameter(nonneg=True)
    x=cp.Variable(n);d=cp.Variable(n)
    prob=cp.Problem(cp.Minimize(.5*cp.sum_squares(zz-cp.multiply(v,x)-d)
                                +.5*cp.sum(cp.huber(Bd@d,1.345))+lam*cp.norm1(Dk@x)))
    k0=cp.Variable();d0=cp.Variable(n)
    constant=cp.Problem(cp.Minimize(.5*cp.sum_squares(zz-v*k0-d0)+.5*cp.sum(cp.huber(Bd@d0,1.345))))
    opts=dict(solver='CLARABEL',tol_gap_abs=1e-8,tol_gap_rel=1e-8,tol_feas=1e-8,max_iter=300)
    rows=[];times=[];grid=(1.,3.,10.,30.,100.);audits=[]
    for replica in range(3):
      for cond in A.CONDITIONS:
        a=np.load(A.OUT/f'{cond}_{replica}.npz');truth=a['log'][a['idx'],21]
        take=np.arange(0,len(truth),3)
        estimates={str(l):[] for l in grid};estimates['constant']=[]
        for j in take:
            K=a['baseline'][j,0];sigma=a['noise'][j]
            v.value=K*a['Y'][j]/sigma;zz.value=a['z'][j]/sigma
            constant.solve(**opts)
            if constant.status!='optimal':raise RuntimeError(constant.status)
            estimates['constant'].append(np.clip(K*k0.value,200.,3e5))
            for l in grid:
                lam.value=l;t0=time.perf_counter();prob.solve(**opts);times.append(time.perf_counter()-t0)
                if prob.status!='optimal':raise RuntimeError(prob.status)
                estimates[str(l)].append(np.clip(K*x.value[-1],200.,3e5))
                if replica==0 and j==take[len(take)//2] and l==30.:
                    ka,diag=fit(a['Y'][j],a['z'][j],sigma,rho,K,l,maxiter=50000,tol=1e-8)
                    audits.append(dict(cond=cond,endpoint_relative_difference=float(abs(ka-K*x.value[-1])/K),
                                       iterations=diag['iterations'],primal=diag['primal'],dual=diag['dual']))
        base=(a['baseline'][take,0]-truth[take])/truth[take]
        for l,vals in estimates.items():
            err=(np.array(vals)-truth[take])/truth[take]
            rows.append(dict(cond=cond,replica=replica,lam=l,score=A.summary(err,np.ones(len(err),bool)),
                             baseline=A.summary(base,np.ones(len(base),bool))))
        print(f'conic replay {cond} {replica}',flush=True)
    feasible=[]
    for l in grid:
        rr=[r for r in rows if r['replica']==0 and r['lam']==str(l)]
        good=all(r['score']['rrmse']<=1.05*r['baseline']['rrmse'] for r in rr if r['cond'] in ('rest','breath','shift_fast'))
        bad=np.mean([r['score']['rrmse']/r['baseline']['rrmse'] for r in rr if r['cond'] in ('rib','cough','combo')])
        if good:feasible.append((float(bad),l))
    result=dict(protocol='Exploratory same-physical-trajectory replay; replica 0 selects lambda. Replicas 1/2 are not final independent validation after model exploration.',
                selected_lambda=min(feasible)[1] if feasible else None,selection_feasible=feasible,rows=rows,
                independent_ADMM_audit=audits,conic_solver='CLARABEL',cvxpy_version=cp.__version__,
                timings_ms=(1000*np.percentile(times,[50,99])).tolist(),n_solves=len(times))
    (A.OUT/'direct_tv.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--conic',action='store_true');p.add_argument('--admm',action='store_true');a=p.parse_args()
    main() if a.admm else conic_main()
