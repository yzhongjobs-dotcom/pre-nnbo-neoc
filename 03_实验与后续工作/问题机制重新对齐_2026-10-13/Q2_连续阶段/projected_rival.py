"""Strong continuous low-level feasible-action rival, same physical task.

Analytic gain/rate/tank projection at each RK stage replaces expensive phase
acceptance. The reference constrained held-rate shooting is a subset of these
controls; a negative excess reveals a reference-parameterization limitation.
"""
import argparse,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize
from contact_stages import Problem,ROOT,dump,flow,margins,safe_rate,environment


@njit(cache=True)
def action(y,desired,p,lo,hi):
    e,v,k,d,E=y;u=np.empty(2)
    for j in range(2):
        lower=max(-p[6+j],-20*(y[2+j]-lo[j]))
        upper=min(p[6+j],20*(hi[j]-y[2+j]))
        if j==0 and e*e>1e-25:
            upper=min(upper,(p[10]*(p[1]+d)*v*v+p[11]*(E-p[8]))/(.5*e*e))
        u[j]=max(lower,min(desired[j]*p[6+j],upper))/p[6+j]
    return u


@njit(cache=True)
def step_projected(y,desired,t,h,p,pulses,lo,hi):
    u=action(y,desired,p,lo,hi);a,la=flow(y,u,t,p,pulses);ya=y+h*a/2
    ua=action(ya,desired,p,lo,hi);b,lb=flow(ya,ua,t+h/2,p,pulses);yb=y+h*b/2
    ub=action(yb,desired,p,lo,hi);c,lc=flow(yb,ub,t+h/2,p,pulses);yc=y+h*c
    uc=action(yc,desired,p,lo,hi);d,ld=flow(yc,uc,t+h,p,pulses)
    return y+h*(a+2*b+2*c+d)/6,h*(la+2*lb+2*lc+ld)/6,ya,yb,yc,u,ua,ub,uc


@njit(cache=True)
def target_rate(y,target,h,p):
    return np.array([max(-1.,min((target[0]-y[2])/(h*p[6]),1.)),
                     max(-1.,min((target[1]-y[3])/(h*p[7]),1.))])


@njit(cache=True)
def future(first,y0,index,targets,h,sub,dt,p,pulses,lo,hi):
    y=y0.copy();J=0.;minimum=1e99;power_min=1e99
    for j in range(index,len(targets)):
        desired=first if j==index else target_rate(y,targets[j],h,p)
        for s in range(sub):
            t=j*h+s*dt;yn,ell,ya,yb,yc,u,ua,ub,uc=step_projected(y,desired,t,dt,p,pulses,lo,hi)
            for state,rate,phase in ((y,u,t),(ya,ua,t+dt/2),(yb,ub,t+dt/2),(yc,uc,t+dt)):
                mar=margins(state,rate,phase,p,pulses)
                for k in range(4):minimum=min(minimum,mar[k])
                power_min=min(power_min,mar[4])
            y=yn;J+=ell
    th,gam=environment(len(targets)*h,p,pulses)
    J+=p[29]*((y[0]/p[14])**2+p[20]*(th*(gam-y[0])/p[16])**2)
    return J,minimum,power_min


def policy(p,targets,improve=True,causal=False,tail_parameters=None,frozen=False,
           maxiter=30,tolerance=1e-8,use_grid=True,warm_start=False,trend=False):
    y=p.y0.copy();rows=[[0.,*y,0.,0.,*environment(0.,p.p,p.pulses)]]
    total=0.;timings=[];desireds=[];diagnostics=[]
    # Exclude JIT only from service time, keep its measured preparation field.
    warm_model=p.p
    if causal or frozen:
        forecast=Problem(reserve=p.cfg['reserve'],segments=p.n,theta_span=0.)
        forecast.p[4]=p.p[4];forecast.p[5]=0.
        warm_model=forecast.p;targets=forecast.targets('preview',tail_parameters)
    tic=time.perf_counter();future(np.zeros(2),y,0,targets,p.h,p.sub,p.dt,warm_model,p.pulses,p.lo,p.hi)
    setup=time.perf_counter()-tic
    observed=[]
    for j in range(p.n):
        tic=time.perf_counter();incumbent=target_rate(y,targets[j],p.h,p.p)
        planning=p.p
        if causal or frozen:
            # Sensor information is only theta at the current physical time.
            # The future material path is absent from this forecast copy.
            theta_now=p.p[4] if frozen else environment(j*p.h,p.p,p.pulses)[0]
            forecast=Problem(reserve=p.cfg['reserve'],segments=p.n,theta_span=0.)
            forecast.p[4]=theta_now;forecast.p[5]=0.
            slope=0.
            if trend:
                if observed:slope=(theta_now-observed[-1]['observed_theta'])/p.h
                forecast.p[31]=1.;forecast.p[32]=j*p.h;forecast.p[33]=slope
            planning=forecast.p
            targets=forecast.targets('preview',tail_parameters)
            incumbent=target_rate(y,targets[j],p.h,planning)
            observed.append(dict(time=j*p.h,observed_theta=float(theta_now),
                                 forecast_theta0=float(planning[4]),forecast_theta_span=float(planning[5]),
                                 estimated_slope=float(slope)))
        cache={}
        def evaluate(z):
            key=tuple(z)
            if key not in cache:cache[key]=future(z,y,j,targets,p.h,p.sub,p.dt,planning,p.pulses,p.lo,p.hi)
            return cache[key]
        value,margin,_=evaluate(incumbent);successes=[]
        if improve:
            starts=[incumbent]
            if warm_start and desireds:
                old=desireds[-1];v,m,_=evaluate(old)
                if m>=-1e-8 and v<value:incumbent=old;value=v
            if use_grid:
                for a in [-1.,0.,1.]:
                    for b in [-1.,0.,1.]:
                        z=np.array([a,b]);v,m,_=evaluate(z)
                        if m>=-1e-8 and v<value:incumbent=z;value=v
            if np.linalg.norm(incumbent-starts[0])>1e-10:starts.append(incumbent)
            for start in starts:
                res=minimize(lambda z:evaluate(z)[0],start,method='SLSQP',bounds=[(-1,1)]*2,
                   constraints=[dict(type='ineq',fun=lambda z:evaluate(z)[1])],
                   options=dict(maxiter=maxiter,ftol=tolerance,eps=1e-5))
                v,m,_=evaluate(res.x)
                if res.success and m>=-1e-8 and v<value:incumbent=res.x;value=v
                successes.append(dict(success=bool(res.success),iterations=int(res.nit),evaluations=int(res.nfev)))
        if evaluate(incumbent)[1]<-1e-8:raise RuntimeError('Physical state constraints not met by incumbent')
        timings.append(time.perf_counter()-tic);diagnostics.append(successes);desireds.append(incumbent)
        for s in range(p.sub):
            t=j*p.h+s*p.dt;yn,ell,*_=step_projected(y,incumbent,t,p.dt,p.p,p.pulses,p.lo,p.hi)
            y=yn;total+=ell
            u=action(y,incumbent,p.p,p.lo,p.hi)*p.rate
            rows.append([t+p.dt,*y,*u,*environment(t+p.dt,p.p,p.pulses)])
    th,gam=environment(p.cfg['horizon'],p.p,p.pulses)
    total+=p.p[29]*((y[0]/p.p[14])**2+p.p[20]*(th*(gam-y[0])/p.p[16])**2)
    tr=np.array(rows);fn=p.cfg['preload']+tr[:,8]*(tr[:,9]-tr[:,1]);fr=-tr[:,3]*tr[:,1]-tr[:,4]*tr[:,2]
    margin=[];power=[]
    for row in tr:
        u=row[6:8]/p.rate;mar=margins(row[1:6],u,row[0],p.p,p.pulses);margin.append(mar[:4]);power.append(mar[4])
    metrics=dict(cost=total,position_rms_mm=float(np.sqrt(np.trapz(tr[:,1]**2,tr[:,0])/p.cfg['horizon'])*1e3),
         force_rms_N=float(np.sqrt(np.trapz((fn-p.cfg['preload'])**2,tr[:,0])/p.cfg['horizon'])),
         minimum_margin=float(np.min(margin)),min_power_margin=float(np.min(power)),min_energy=float(tr[:,5].min()),
         max_increment_force=float(abs(fr).max()),min_contact=float(fn.min()),max_contact=float(fn.max()),
         max_K=float(tr[:,3].max()),max_D=float(tr[:,4].max()))
    timing=dict(setup_seconds=setup,total_seconds=sum(timings),worst_seconds=max(timings),
                 decision_seconds=timings,local_diagnostics=diagnostics,material_observations=observed)
    return metrics,timing,tr,np.array(desireds)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reserve',type=float,default=.00005)
    parser.add_argument('--name',default='projected_rival');parser.add_argument('--causal',action='store_true');args=parser.parse_args()
    out=ROOT/args.name;out.mkdir(exist_ok=True);p=Problem(reserve=args.reserve)
    calibration=json.loads((ROOT/'strong_rivals'/'results.json').read_text())
    best=calibration['best_preview'];targets=p.targets(best['kind'],[best['zeta'],best['lead']])
    results=dict(config=p.cfg,causal_current_material_only=args.causal,tail_parameters=dict(zeta=best['zeta'],lead=best['lead']),
            scope='Low-level analytic viable-action projection; held-rate optimum is a restricted reference, not global optimum.')
    for label,improvement in [('preview',False),('rollout',True)]:
        m,t,tr,z=policy(p,targets,improvement,args.causal,[best['zeta'],best['lead']]);results[label]=dict(metrics=m,timing=t)
        dump(out/'results.json',results);np.save(out/(label+'_desired_rates.npy'),z)
        np.savetxt(out/(label+'_trace.csv'),tr,delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
        print(label,m,'worst seconds',t['worst_seconds'],flush=True)


if __name__=='__main__':main()
