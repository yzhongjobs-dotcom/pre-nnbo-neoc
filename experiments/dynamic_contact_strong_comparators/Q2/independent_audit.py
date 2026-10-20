"""Independent equations, adaptive integration, storage and static-target check."""
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from contact_stages import Problem,ROOT,dump,environment


def env(t,c):
    a=np.clip(t/c['horizon'],0,1);theta=c['theta0']+c['theta_span']*(1-np.cos(np.pi*a))/2
    gamma=sum(A*np.sin(np.pi*(t-t0)/duration)**2 for t0,duration,A in c['pulses'] if t0<t<t0+duration)
    return theta,gamma


def independent(p,desired,projected=False):
    c=p.cfg;y=p.y0.copy();cost=0.;minimum_energy=y[4];power_min=1e99;segment_states=[y.copy()]
    def actual(x,z):
        u=z*p.rate
        if projected:
            lower=np.maximum(-p.rate,-20*(x[2:4]-p.lo))
            upper=np.minimum(p.rate,20*(p.hi-x[2:4]))
            if x[0]**2>1e-25:
                upper[0]=min(upper[0],(c['recharge']*(c['robot_damping']+x[3])*x[1]**2+
                               c['energy_alpha']*(x[4]-c['energy_min']))/(.5*x[0]**2))
            u=np.clip(u,lower,upper)
        return u
    for j,z in enumerate(desired.reshape(p.n,2)):
        t0=j*p.h;t1=(j+1)*p.h
        def rhs(t,x):
            e,v,k,d,E=x[:5];th,gam=env(t,c);uk,ud=actual(x[:5],z)
            fr=-k*e-d*v;ferr=th*(gam-e)
            ell=c['w_position']*(e/c['error_scale'])**2+c['w_velocity']*(v/c['velocity_scale'])**2
            ell+=c['w_force']*(ferr/c['force_scale'])**2+c['w_command']*(fr/c['command_scale'])**2
            ell+=c['w_compliance']*(((k-p.pref[0])/c['gain_scale'][0])**2+((d-p.pref[1])/c['gain_scale'][1])**2)
            ell+=c['w_rate']*((uk/p.rate[0])**2+(ud/p.rate[1])**2)
            return [v,(fr+ferr-c['robot_damping']*v)/c['mass'],uk,ud,
                    c['recharge']*(c['robot_damping']+d)*v*v-.5*e*e*uk,ell]
        sol=solve_ivp(rhs,[t0,t1],np.r_[y,0.],method='DOP853',rtol=2e-10,atol=1e-13,
                      max_step=p.dt,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        xs=sol.sol(np.linspace(t0,t1,101));minimum_energy=min(minimum_energy,float(xs[4].min()))
        for x in xs[:5].T:
            uk,ud=actual(x,z);e,v,k,d,E=x
            power_min=min(power_min,c['recharge']*(c['robot_damping']+d)*v*v+
                c['energy_alpha']*(E-c['energy_min'])-.5*e*e*uk)
        y=sol.y[:5,-1];cost+=sol.y[5,-1];segment_states.append(y.copy())
    th,gam=env(c['horizon'],c)
    cost+=c['terminal']*((y[0]/c['error_scale'])**2+c['w_force']*(th*(gam-y[0])/c['force_scale'])**2)
    return dict(cost=cost,min_energy=minimum_energy,min_power_W=power_min,segment_states=np.array(segment_states))


def main():
    p=Problem();r=json.loads((ROOT/'strong_rivals'/'results.json').read_text());rng=np.random.default_rng(29)
    z=np.array(r['independent_shooting_from_rollout']['rates']);J,c,tr=p.run(z,True)
    ind=independent(p,z);fine=Problem(dt=p.dt/4);mf=fine.metrics(z)
    gradients=[]
    for _ in range(3):
        point=z*.8;direction=rng.normal(size=len(z));direction/=np.linalg.norm(direction)
        v,g,c,A=p.evaluate(point);eps=1e-6
        ja,ca=p.run(point+eps*direction);jb,cb=p.run(point-eps*direction)
        gradients.append(dict(cost_error=float(abs(g@direction-(ja[0]-jb[0])/2/eps)),
                      constraint_error=float(np.max(abs(A@direction-(ca-cb)[0]/2/eps)))))
    storage=[]
    for _ in range(100):
        e=rng.uniform(-.002,.002);v=rng.uniform(-.03,.03);k=rng.uniform(20,1200);d=rng.uniform(1,30)
        t=rng.uniform(0,4);E=.00105;uk=rng.uniform(-400,400);th,gam=env(t,p.cfg)
        acc=(-(th+k)*e-(p.cfg['robot_damping']+d)*v+th*gam)/p.cfg['mass']
        Edot=p.cfg['recharge']*(p.cfg['robot_damping']+d)*v*v-.5*e*e*uk
        deriv=p.cfg['mass']*v*acc+.5*uk*e*e+k*e*v+Edot
        supply=th*(gam-e)*v-(1-p.cfg['recharge'])*(p.cfg['robot_damping']+d)*v*v
        storage.append(abs(deriv-supply))
    target_errors=[];targets=p.targets('quasi',[.15,0])
    for j in [0,6,7,8,22,23,24]:
        th,gam=env(j*p.h,p.cfg)
        def cost(k):
            e=th*gam/(th+k)
            return p.p[18]*(e/p.p[14])**2+(p.p[20]/p.p[16]**2+p.p[21]/p.p[17]**2)*(k*e)**2+p.p[22]*((k-p.pref[0])/p.p[26])**2
        fit=minimize_scalar(cost,bounds=(p.lo[0],p.hi[0]),method='bounded',options=dict(xatol=1e-9))
        target_errors.append(abs(fit.x-targets[j,0]))
    proj=json.loads((ROOT/'projected_rival'/'results.json').read_text())
    zp=np.load(ROOT/'projected_rival'/'rollout_desired_rates.npy');pi=independent(p,zp,True)
    pi_states=pi.pop('segment_states');ind_states=ind.pop('segment_states')
    proj_tr=np.loadtxt(ROOT/'projected_rival'/'rollout_trace.csv',delimiter=',',skiprows=1)
    out=dict(direct_adaptive_integral=ind,direct_cost_relative_difference=abs(ind['cost']-J[0])/J[0],
        direct_step_divide4=mf,direct_segment_state_error=np.max(abs(ind_states-tr[::p.sub,1:6]),axis=0).tolist(),
        derivative_checks=gradients,storage_identity_error=max(storage),
        static_target_K_max_error=max(target_errors),projected_adaptive_integral=pi,
        projected_cost_relative_difference=abs(pi['cost']-proj['rollout']['metrics']['cost'])/pi['cost'],
        projected_segment_state_error=np.max(abs(pi_states-proj_tr[::p.sub,1:6]),axis=0).tolist(),
        scope='Independent numerical and sampled checks; not a global optimum, deployment or interval certificate.')
    dump(ROOT/'audit.json',out);print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':main()
