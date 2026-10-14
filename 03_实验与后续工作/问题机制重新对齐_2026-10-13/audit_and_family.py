"""Independent ODE/value/gradient checks and material/state hold-out cases."""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm,solve_continuous_lyapunov
from scipy.optimize import minimize,differential_evolution
from dynamic_tradeoff import Problem,tune,dump,ROOT


def independent_integral(p,z):
    y=p.y0.copy();cost=0.;tr=[y.copy()]
    for u in np.asarray(z).reshape(p.n,2):
        # Independently implement the physical equations, not call RK4 rhs.
        def rhs(t,x):
            e,v,k,d,E=x[:5];uk,ud=u*p.rate;c=p.cfg;w=c['weights']
            fr=-k*e-d*v
            ell=w['position']*(e/c['error_scale'])**2+w['force']*(p.theta*e/c['force_scale'])**2
            ell+=w['velocity']*(v/c['velocity_scale'])**2+w['command']*(fr/c['command_scale'])**2
            ell+=w['compliance']*(((k-p.pref[0])/p.gscale[0])**2+((d-p.pref[1])/p.gscale[1])**2)
            ell+=w['rate']*np.sum(u*u)
            return [v,(-(p.theta+k)*e-(c['robot_damping']+d)*v)/c['mass'],
                    uk,ud,c['recharge']*(c['robot_damping']+d)*v*v-.5*e*e*uk,ell]
        sol=solve_ivp(rhs,[0,p.h],np.r_[y,0.],method='DOP853',rtol=2e-11,atol=1e-13)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:5,-1];cost+=sol.y[5,-1];tr.append(y.copy())
    cost+=p.cfg['weights']['terminal']*((y[0]/p.cfg['error_scale'])**2+
                  p.cfg['weights']['force']*(p.theta*y[0]/p.cfg['force_scale'])**2)
    return cost,np.array(tr)


def constant_oracle(p):
    x=p.y0[:2];c=p.cfg;w=c['weights']
    def cost(g):
        k,d=g;A=np.array([[0.,1.],[-(p.theta+k)/c['mass'],-(c['robot_damping']+d)/c['mass']]])
        Q=np.diag([w['position']/c['error_scale']**2+w['force']*(p.theta/c['force_scale'])**2,
                    w['velocity']/c['velocity_scale']**2])
        Q+=w['command']/c['command_scale']**2*np.outer(g,g)
        P=solve_continuous_lyapunov(A.T,-Q);final=expm(A*c['horizon'])@x
        return float(x@P@x-final@P@final+c['horizon']*w['compliance']*np.sum(((g-p.pref)/p.gscale)**2)+
                     w['terminal']*(1/c['error_scale']**2+w['force']*(p.theta/c['force_scale'])**2)*final[0]**2)
    runs=[minimize(cost,start,method='L-BFGS-B',bounds=list(zip(p.lo,p.hi)),
              options=dict(maxiter=200,ftol=1e-14,gtol=1e-9)) for start in [p.pref,[200,8],[600,15]]]
    best=min(runs,key=lambda r:r.fun)
    saved=p.y0.copy();p.y0[2:4]=best.x
    metrics=p.metrics(np.zeros(2*p.n));p.y0=saved
    return dict(gains=best.x.tolist(),cost=float(best.fun),optimizer_success=bool(best.success),
         integral_relative_error=abs(metrics['cost']-best.fun)/best.fun,
         metrics=metrics,scope='Generous fixed-gain oracle may preselect actual initial gains; no preparation charge.')


def repair_state(p,initial,global_search=False):
    low=np.array([0.,0.,-4.,-3.,-4.,.025]);high=np.array([1180.,27.,-2.,-1.,4.,.25])
    def phys(q):
        a=low+(high-low)*q;a[2:4]=10**a[2:4];return a
    def norm(a):
        a=np.array(a);a[2:4]=np.log10(a[2:4]);return (a-low)/(high-low)
    def objective(q):
        z=p.policy(phys(q),'state');J,c=p.rollout(z)
        return float(J[0]+1e4*np.sum(np.minimum(c,0)**2))
    if global_search:
        res=differential_evolution(objective,[(0,1)]*6,seed=15,popsize=8,maxiter=60,tol=1e-7,polish=False)
        start=res.x
    else:start=norm(initial)
    res=minimize(objective,start,method='Powell',bounds=[(0,1)]*6,
                options=dict(maxiter=100,xtol=1e-8,ftol=1e-10))
    z=p.policy(phys(res.x),'state')
    return dict(parameters=phys(res.x).tolist(),optimizer_success=bool(res.success),
                message=str(res.message),metrics=p.metrics(z)),z


def audit():
    folder=ROOT/'pilot_compiled';r=json.loads((folder/'results.json').read_text())
    p=Problem();idx=r['best_optimization']['index'];z=np.load(folder/f'opt_{idx}_rates.npy')
    J,c,tr=p.rollout(z,True);ind,itr=independent_integral(p,z)
    fine=Problem(dt=p.dt/2);Jfine,cfine,tfine=fine.rollout(z,True)
    rng=np.random.default_rng(15);gradient=[]
    for i in range(5):
        point=z*.8+rng.normal(0,.002,len(z));d=rng.normal(size=len(z));d/=np.linalg.norm(d)
        value,g,margin,A=p.evaluate(point);eps=1e-6
        jp,cp=p.rollout(point+eps*d);jm,cm=p.rollout(point-eps*d)
        gradient.append(dict(objective_abs=float(abs(g@d-(jp[0]-jm[0])/2/eps)),
                      constraint_max_abs=float(np.max(abs(A@d-(cp-cm)[0]/2/eps)))))
    K=p.y0[2]+p.gmap@z
    # Structural storage identity, independently derived from robot force.
    defects=[]
    for _ in range(50):
        e=rng.uniform(-.002,.002);v=rng.uniform(-.02,.02);k=rng.uniform(20,1200);d=rng.uniform(1,30)
        E=.002;uk=rng.uniform(-3000,3000);ud=rng.uniform(-80,80);fr=-k*e-d*v
        acc=(-(p.theta+k)*e-(p.cfg['robot_damping']+d)*v)/p.cfg['mass']
        edot=p.cfg['recharge']*(p.cfg['robot_damping']+d)*v*v-.5*e*e*uk
        storage=p.cfg['mass']*v*acc+.5*uk*e*e+k*e*v+edot
        supply=v*(-p.theta*e)-(1-p.cfg['recharge'])*(p.cfg['robot_damping']+d)*v*v
        defects.append(abs(storage-supply))
    const=constant_oracle(p)
    repaired,zstate=repair_state(p,r['state']['parameters'],True)
    np.save(folder/'state_repaired_rates.npy',zstate)
    report=dict(independent_DOP853_cost=ind,relative_cost_difference=abs(ind-J[0])/J[0],
      max_segment_state_difference=np.max(abs(itr-tr[::p.sub,1:6]),axis=0).tolist(),
      half_step_relative_cost_difference=float(abs(Jfine[0]-J[0])/J[0]),
      half_step_minimum_margin=float(cfine.min()),gradient_checks=gradient,
      max_energy_identity_defect=float(max(defects)),constant_initial_gain_oracle=const,
      state_global_and_Powell=repaired,
      repaired_state_excess=repaired['metrics']['cost']/r['best_optimization']['metrics']['cost']-1,
      status='diagnostic validation, no global nonlinear optimality or real tissue claim')
    dump(folder/'audit.json',report);print(json.dumps(report,indent=2),flush=True)


def family():
    folder=ROOT/'pilot_compiled';pilot=json.loads((folder/'results.json').read_text())
    audit=json.loads((folder/'audit.json').read_text());pars=audit['state_global_and_Powell']['parameters']
    cases=[(100.,.0015,.008),(600.,.0015,.008),(900.,.0015,.008),
           (300.,.0008,.004),(300.,.0025,-.008),(600.,.002,-.012)]
    report={}
    out=ROOT/'material_and_state_holdout';out.mkdir(exist_ok=True)
    for i,(theta,e0,v0) in enumerate(cases):
        p=Problem(theta,initial_error=e0,initial_velocity=v0);p.evaluate(np.zeros(2*p.n))
        z=p.policy(pars,'state');state=dict(parameters=pars,metrics=p.metrics(z))
        fixed=p.policy(pilot['fixed']['parameters'],'fixed')
        local,zlocal=tune(p,'state',maxiter=100)
        local2,zlocal2=repair_state(p,local['parameters'])
        if local2['metrics']['cost']<local['metrics']['cost']:local,zlocal=local2,zlocal2
        sols=[]
        for start in [z,zlocal,np.zeros(2*p.n)]:
            res,seconds=p.solve(start,maxiter=180);m=p.metrics(res.x)
            sols.append(dict(success=bool(res.success),iterations=int(res.nit),seconds=seconds,
                              metrics=m,rates=res.x.tolist()))
        okay=[r for r in sols if r['success'] and r['metrics']['minimum_margin']>=-1e-6]
        if not okay:raise RuntimeError('No accepted optimum '+str(i))
        opt=min(okay,key=lambda r:r['metrics']['cost']);best=opt['metrics']['cost']
        report[str(i)]=dict(theta=theta,e0=e0,v0=v0,frozen_state=state,
            frozen_fixed=p.metrics(fixed),retuned_state=local,constant_oracle=constant_oracle(p),
            optimization=opt,optimization_runs=sols,
            frozen_state_excess=state['metrics']['cost']/best-1,
            retuned_state_excess=local['metrics']['cost']/best-1)
        dump(out/'results.json',report)
        tr=p.rollout(np.array(opt['rates']),True)[2]
        np.savetxt(out/f'case{i}_opt_trace.csv',tr,delimiter=',',header='t,e,v,K,D,E,uK,uD',comments='')
        print(i,theta,'frozen',report[str(i)]['frozen_state_excess'],
              'retuned',report[str(i)]['retuned_state_excess'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['audit','family']);args=parser.parse_args()
    audit() if args.mode=='audit' else family()
