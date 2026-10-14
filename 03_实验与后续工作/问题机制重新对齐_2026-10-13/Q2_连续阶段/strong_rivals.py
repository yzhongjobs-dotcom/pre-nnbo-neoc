"""Tune and strengthen SAME-information cheap predictors before NEOC testing."""
import argparse,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize,minimize_scalar
from contact_stages import Problem,ROOT,dump,rollout,step,margins,safe_rate


@njit(cache=True)
def extension(first,y0,index,targets,h,sub,dt,p,pulses,lo,hi):
    y=y0.copy();cost=0.;minimum=1e99;n=len(targets)
    for j in range(index,n):
        z=first.copy() if j==index else safe_rate(y,targets[j],j*h,h,sub,dt,p,pulses,lo,hi)
        for s in range(sub):
            t=j*h+s*dt;yn,ell,ya,yb,yc=step(y,z,t,dt,p,pulses)
            for state,now in ((y,t),(ya,t+dt/2),(yb,t+dt/2),(yc,t+dt)):
                for v in margins(state,z,now,p,pulses):minimum=min(minimum,v)
            y=yn;cost+=ell
        for v in margins(y,z,(j+1)*h,p,pulses):minimum=min(minimum,v)
    # Surface is zero at the planned end; use the same terminal comparison.
    from_theta=p[4]+p[5]
    cost+=p[29]*((y[0]/p[14])**2+p[20]*(from_theta*y[0]/p[16])**2)
    return cost,minimum


def improve(p,targets):
    y=p.y0.copy();zs=[];timings=[];details=[]
    extension(np.zeros(2),y,0,targets,p.h,p.sub,p.dt,p.p,p.pulses,p.lo,p.hi)
    for j in range(p.n):
        tic=time.perf_counter()
        incumbent=safe_rate(y,targets[j],j*p.h,p.h,p.sub,p.dt,p.p,p.pulses,p.lo,p.hi)
        bounds=list(zip(np.maximum(-1,(p.lo-y[2:4])/(p.h*p.rate)),np.minimum(1,(p.hi-y[2:4])/(p.h*p.rate))))
        cache={}
        def evaluate(x):
            key=tuple(x)
            if key not in cache:cache[key]=extension(x,y,j,targets,p.h,p.sub,p.dt,p.p,p.pulses,p.lo,p.hi)
            return cache[key]
        value,margin=evaluate(incumbent)
        # Grid seeds and local current-action optimization receive the same
        # future models as the full solve; only the tail policy is cheaper.
        starts=[incumbent]
        grid=np.linspace(-1,1,5)
        for a in grid:
            for b in grid:
                q=np.clip([a,b],np.array(bounds)[:,0],np.array(bounds)[:,1])
                v,m=evaluate(q)
                if m>=-1e-8 and v<value:incumbent=q;value=v
        starts.append(incumbent)
        successes=[]
        for start in starts:
            res=minimize(lambda x:evaluate(x)[0],start,method='SLSQP',bounds=bounds,
               constraints=[dict(type='ineq',fun=lambda x:evaluate(x)[1])],
               options=dict(maxiter=45,ftol=1e-8,eps=3e-5))
            v,m=evaluate(res.x)
            if res.success and m>=-1e-8 and v<value:incumbent=res.x;value=v
            successes.append(dict(success=bool(res.success),iterations=int(res.nit),message=str(res.message)))
        if evaluate(incumbent)[1]<-1e-8:raise RuntimeError('No feasible predictive incumbent')
        zs.append(incumbent);timings.append(time.perf_counter()-tic);details.append(successes)
        for s in range(p.sub):y=step(y,incumbent,j*p.h+s*p.dt,p.dt,p.p,p.pulses)[0]
    return np.array(zs).ravel(),dict(total_seconds=sum(timings),worst_seconds=max(timings),
             decision_seconds=timings,local_solve_diagnostics=details)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reserve',type=float,default=.00005)
    parser.add_argument('--name',default='strong_rivals');args=parser.parse_args()
    out=ROOT/args.name;out.mkdir(exist_ok=True);p=Problem(reserve=args.reserve)
    p.evaluate(np.zeros(2*p.n));p.policy(p.targets('fixed',[50,3]))
    rows=[]
    # Fixed oracle and independent phase-local solutions are tuned, not just
    # three arbitrarily chosen gains or damping settings.
    def fixed_cost(g):
        z,_=p.policy(p.targets('fixed',g));J,c=p.run(z)
        return float(J[0]+1e5*np.sum(np.minimum(c,0)**2))
    fixed=[]
    for seed in [[50,3],[200,8],[500,12]]:
        res=minimize(fixed_cost,seed,method='Powell',bounds=list(zip(p.lo,p.hi)),
                     options=dict(maxiter=50,xtol=1e-5,ftol=1e-9))
        z,secs=p.policy(p.targets('fixed',res.x));fixed.append(dict(gains=res.x.tolist(),success=bool(res.success),metrics=p.metrics(z),rates=z.tolist(),seconds=secs))
    best_fixed=min(fixed,key=lambda r:r['metrics']['cost'])
    for lead in [0.,.125,.25,.375,.5,.75,1.,1.5,4.]:
        kind='quasi' if lead==0 else 'preview'
        def objective(zeta):
            target=p.targets(kind,[zeta,lead]);z,_=p.policy(target);J,c=p.run(z)
            return float(J[0]+1e5*np.sum(np.minimum(c,0)**2))
        fits=[]
        for lower,upper in [(.03,.2),(.2,.7)]:
            res=minimize_scalar(objective,bounds=(lower,upper),method='bounded',options=dict(xatol=1e-5))
            fits.append(res)
        best=min(fits,key=lambda r:r.fun);target=p.targets(kind,[best.x,lead]);z,seconds=p.policy(target)
        row=dict(kind=kind,zeta=float(best.x),lead=lead,metrics=p.metrics(z),rates=z.tolist(),seconds=seconds)
        rows.append(row)
        print(kind,lead,'cost',row['metrics']['cost'],'margin',row['metrics']['minimum_margin'],flush=True)
    feasible=[r for r in rows if r['metrics']['minimum_margin']>=-1e-8]
    best=min(feasible,key=lambda r:r['metrics']['cost'])
    result=dict(config=p.cfg,fixed_fits=fixed,best_fixed=best_fixed,preview_fits=rows,best_preview=best)
    dump(out/'results.json',result)
    target=p.targets(best['kind'],[best['zeta'],best['lead']])
    z,timing=improve(p,target);m=p.metrics(z)
    result['rollout']=dict(metrics=m,timing=timing,rates=z.tolist());dump(out/'results.json',result)
    np.save(out/'rollout_rates.npy',z);np.savetxt(out/'rollout_trace.csv',p.run(z,True)[2],delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
    print('rollout',m,'worst seconds',timing['worst_seconds'],flush=True)
    res,seconds=p.solve(z)
    result['independent_shooting_from_rollout']=dict(success=bool(res.success),message=str(res.message),
         seconds=seconds,iterations=int(res.nit),metrics=p.metrics(res.x),rates=res.x.tolist())
    dump(out/'results.json',result);np.save(out/'independent_opt_rates.npy',res.x)
    np.savetxt(out/'independent_opt_trace.csv',p.run(res.x,True)[2],delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
    print('independent optimization',result['independent_shooting_from_rollout']['metrics'],flush=True)


if __name__=='__main__':main()
