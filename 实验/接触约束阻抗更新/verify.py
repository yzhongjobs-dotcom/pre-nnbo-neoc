"""Independent projection, energy identity and adaptive-ODE checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize
import model as m

HERE=Path(__file__).resolve().parent


def rhs(y,target,p):
    # Rebuild from physical scalars, without calling model.flow/action/contact.
    physical=y[4:8]*np.array([500.,500.,40.,40.])
    lo=np.maximum(-np.array([1.6,1.6,2.,2.]),-40*(y[4:8]-np.array([.01,.01,.05,.05])))
    hi=np.minimum(np.array([1.6,1.6,2.,2.]),40*(np.array([3.2,3.2,3.,3.])-y[4:8]))
    raw=(target-y[4:8])/.025
    a=.5*np.r_[y[:2]**2*np.array([500.,500.]),0.,0.]
    dissipated=physical[2]*y[2]**2+physical[3]*y[3]**2+.5*((y[2]+.008)**2+y[3]**2)
    budget=.5*dissipated+40*(y[8]-1e-5)
    u=np.clip(raw,lo,hi)
    if a@u > budget:
        def residual(lam):return a@np.clip(raw-lam*a,lo,hi)-budget
        upper=1.
        while residual(upper)>0:upper*=10.
        u=np.clip(raw-brentq(residual,0,upper,xtol=1e-10)*a,lo,hi)
    depth=max(y[1]+4.5/p[0],0.)
    force=3*(depth/(4.5/p[0]))**1.5
    friction=p[1]*force*np.tanh((y[2]+.008)/.002)
    command=np.array([3*p[1]-physical[0]*y[0]-physical[2]*y[2],
                      3-physical[1]*y[1]-physical[3]*y[3]])
    executed=np.clip(command,-12.,12.)
    acceleration=(executed-np.array([friction,force])-.5*np.array([y[2]+.008,y[3]]))/np.array([.8,.5])
    f=np.r_[y[2:4],acceleration,u,.5*dissipated-a@u]
    cost=(y[0]/.002)**2+(y[1]/.006)**2+np.sum((y[2:4]/.03)**2)+((force-3)/.7)**2
    cost+=np.sum(((executed-np.array([3*p[1],3.]))/4.)**2)
    return f,cost,force,u,a,budget,lo,hi,raw,executed,np.array([friction,force])


def independent(targets,state,param,delay=0.):
    value=np.r_[state,0.]
    switches=sorted(set([0.,min(delay,.3),min(delay+.15,.3),.3]))
    min_energy=state[8]-1e-5
    for left,right in zip(switches[:-1],switches[1:]):
        if right<=left:continue
        target=state[4:8] if left<delay-1e-12 else targets[min(int(max(left-delay,0.)/.15+1e-8),1)]
        def f(t,y):
            dy,cost,*_=rhs(y[:9],target,param)
            return np.r_[dy,cost]
        sol=solve_ivp(f,(left,right),value,method='DOP853',rtol=2e-10,atol=2e-12,max_step=.001)
        assert sol.success
        min_energy=min(min_energy,float(sol.y[8].min()-1e-5))
        value=sol.y[:,-1]
    return float(value[-1]),value[:9],min_energy


def main():
    result=json.loads((HERE/'results.json').read_text());raw=np.load(HERE/'paired.npz')
    params,states=m.sample(64,7401);rng=np.random.default_rng(7402)
    # Explicitly exercise active halfspaces; these are algebra checks, not task trials.
    extra=states[:16].copy();extra[:,:2]=rng.uniform(.003,.01,(16,2))
    extra[:,2:4]=0.;extra[:,8]=m.E_MIN+1e-8
    params=np.vstack((params,params[:16]));states=np.vstack((states,extra))
    projection_errors=[];energy_errors=[];vector_errors=[];active_count=0
    for i,(y,p) in enumerate(zip(states,params)):
        target=m.HIGH.copy() if i>=64 else m.normalized_targets(rng.uniform(0,1,4))[0]
        f,c,force,u,a,budget,lo,hi,rawu,executed,fc=rhs(y,target,p)
        actual=m.flow(y,target,p)
        vector_errors.append(float(max(np.max(abs(f-actual[0])),abs(c-actual[1]),abs(force-actual[2]))))
        norm=max(float(np.linalg.norm(a)),1e-10)
        active_count+=int(a@np.clip(rawu,lo,hi)>budget)
        constraints=[dict(type='ineq',fun=lambda v:(budget-a@v)/norm,
                          jac=lambda v:-a/norm)]
        objective_scale=max(float(np.linalg.norm(rawu)),1.)
        sol=minimize(lambda v:.5*np.sum((v-rawu)**2)/objective_scale,np.zeros(4),
                     jac=lambda v:(v-rawu)/objective_scale,
                     bounds=list(zip(lo,hi)),constraints=constraints,method='SLSQP',
                     options=dict(ftol=1e-11,maxiter=100))
        assert sol.success,(i,sol.message)
        projection_errors.append(float(np.max(abs(sol.x-u))))
        # H = kinetic error + virtual spring + tank; reference/feedforward ports retained.
        g=y[4:8]*np.array([500.,500.,40.,40.]);v=np.array([y[2]+.008,y[3]])
        hdot=np.sum(np.array([.8,.5])*y[2:4]*f[2:4])+np.sum(g[:2]*y[:2]*f[:2])
        hdot+=.5*np.sum(y[:2]**2*f[4:6]*500.)+f[8]
        port=y[2:4]@(executed+g[:2]*y[:2]+g[2:]*y[2:4]-fc)
        expected=port-np.sum(g[2:]*y[2:4]**2)-.5*y[2:4]@v
        expected+=.5*(np.sum(g[2:]*y[2:4]**2)+.5*np.sum(v*v))
        energy_errors.append(float(abs(hdot-expected)))
    checks=[]
    # Preselected first two test contexts; no choosing by the observed ranking.
    for case in [0,1]:
        s=raw['test_states'][case];p=raw['test_params'][case]
        for name in ['hold','fixed','critical','library','teacher','actor_11','actor_library']:
            target=raw['test_'+name+'_targets'][case]
            jr=float(raw['test_'+name+'_J'][case]);j,final,emin=independent(target,s,p)
            fine=m.rollout(target,s,p,dt=.0005)
            checks.append(dict(case=case,method=name,adaptive_J=j,
                rk4_relative_error=abs(jr/j-1),twofold_refinement_relative_change=abs(float(fine[0])/jr-1),
                independent_minimum_energy_margin_J=emin))
        for name in ['library','teacher']:
            delay=min(float(raw['test_'+name+'_ms'][case])/1000,.3)
            j,_,_=independent(raw['test_'+name+'_targets'][case],s,p,delay=delay)
            stored=float(raw['test_'+name+'_delayed_J'][case])
            checks.append(dict(case=case,method=name+'_measured_release',delay_s=delay,
                               rk4_relative_error=abs(stored/j-1)))
    import torch
    reload_errors=[]
    for seed in [11,22,33]:
        data=torch.load(HERE/f'actor_{seed}.pt',weights_only=False)
        net=torch.nn.Sequential(torch.nn.Linear(11,64),torch.nn.Tanh(),torch.nn.Linear(64,64),
                                torch.nn.Tanh(),torch.nn.Linear(64,8),torch.nn.Sigmoid()).double()
        net.load_state_dict(data['weights']);net.eval()
        x=m.features(raw['test_params'],raw['test_states'])
        with torch.no_grad():out=net(torch.as_tensor((x-data['mean'])/data['std'],dtype=torch.float64)).numpy()
        target=m.normalized_targets(out.reshape(-1,4)).reshape(-1,2,4)
        reload_errors.append(float(np.max(abs(target-raw[f'test_actor_{seed}_targets']))))
    out=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
             independent_states=80,constructed_active_halfspace_checks=16,active_halfspace_cases=active_count,
             maximum_vector_or_cost_difference=max(vector_errors),
             maximum_projection_difference=max(projection_errors),maximum_energy_identity_error=max(energy_errors),
             max_actor_reload_target_error=max(reload_errors),adaptive_ODE=checks,
             max_adaptive_relative_cost_error=max(c['rk4_relative_error'] for c in checks),
             scope='Local model, numerical projection and energy accounting only; not robot validation or force/passivity certification')
    assert active_count>=16
    assert max(vector_errors)<1e-7 and max(projection_errors)<1e-6 and max(energy_errors)<1e-10
    assert max(reload_errors)<1e-12
    assert out['max_adaptive_relative_cost_error']<.005
    (HERE/'independent_checks.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
