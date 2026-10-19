"""Independent contact integral, execution, task gradient and value checks."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
import torch
import physics as p
from local_hjb import LocalHJB,CENTRE
import learning

def main():
    rng=np.random.default_rng(16701);results={};errs=[];qerrs=[]
    eta=np.array([1e5,1.5e5,.15,.22,2.,3.,0.,0.,0.,-.001,.001,.01]);x,eta=p.initial(eta)
    for _ in range(12):
        xx=x.copy();xx[:4]=rng.uniform(-1,1,4)*np.array([.001,.00015,.002,.002]);a=float(p.halfwidth(eta[10]-xx[1]));s=eta[9]+xx[0];v=eta[11]+xx[2]
        def pressure(u):return p.L/(2*p.R)*(eta[0] if s+u<eta[6] else eta[1])*np.sqrt(max(0,a*a-u*u))
        def traction(u):
            j=0 if s+u<eta[6] else 1
            return -np.sign(v)*eta[2+j]*pressure(u)-v*eta[4+j]*np.sqrt(max(0,a*a-u*u))/(np.pi*a*a/2)
        points=[r for r in [eta[6]-s] if -a<r<a]
        ff=np.array([quad(pressure,-a,a,points=points,epsabs=1e-10)[0],quad(traction,-a,a,points=points,epsabs=1e-10)[0],
                     quad(lambda u:u*u/(2*p.R)*traction(u)-u*pressure(u),-a,a,points=points,epsabs=1e-10)[0]])
        errs.append(float(np.max(abs(ff-p.wrench(0,xx,eta)))))
    for _ in range(50):
        xx=x.copy();xx[:2]=rng.uniform(-.003,.003,2);xx[8]=p.EMIN+rng.uniform(0,.00001)
        _,_,lo,hi,power,budget,_=p.terms(0,xx,eta);linear=rng.normal(0,.1,4);aa=p.qp(linear,lo,hi,power,budget)
        opt=minimize(lambda z:np.sum(p.RA*z*z)+linear@z,np.clip(np.zeros(4),lo,hi),jac=lambda z:2*p.RA*z+linear,
                     bounds=list(zip(lo,hi)),constraints=[{'type':'ineq','fun':lambda z:budget-power@z[:2],
                     'jac':lambda z:-np.r_[power,0.,0.]}],method='SLSQP',options={'ftol':1e-12,'maxiter':100})
        qerrs.append(float(abs(np.sum(p.RA*aa*aa)+linear@aa-opt.fun)))
        assert np.all(aa>=lo-1e-8) and np.all(aa<=hi+1e-8) and power@aa[:2]<=budget+1e-8
    # Independent torch task objective versus numpy midpoint plus projection.
    action=np.array([.05,-.07,.02,-.01]);xx=x.copy();val=0.;dt=.01
    for i in range(round(p.HORIZON/dt)):
        xx,c,_=p.step(i*dt,xx,eta,action,dt);val+=dt*c
    val+=p.terminal(xx,eta)
    tx=torch.tensor(x[None],dtype=torch.float64);te=torch.tensor(eta[None],dtype=torch.float64);ta=torch.tensor(action[None],dtype=torch.float64,requires_grad=True)
    tv=learning.task_loss(tx,te,ta);tv.sum().backward();direction=rng.normal(size=4);h=1e-5
    vv=[]
    for sign in (-1,1):
        v=learning.task_loss(tx,te,torch.tensor((action+sign*h*direction)[None],dtype=torch.float64));vv.append(float(v.detach()))
    fd=(vv[1]-vv[0])/(2*h);ad=float(ta.grad.numpy()[0]@direction)
    solver=LocalHJB();policy,info=solver.solve(x,eta)
    data=solver.fun(0.,np.zeros(9),np.zeros(4),eta)
    point,_,_=p.step(0.,CENTRE.copy(),eta,np.zeros(4),.01)
    results={'contact_max_absolute_error':max(errs),'QP_max_cost_difference':max(qerrs),
             'torch_numpy_cost_difference':float(tv.detach())-float(val),
             'task_directional_gradient_relative_error':abs(ad-fd)/max(abs(fd),1e-8),
             'symbolic_numpy_midpoint_state_difference':float(np.max(abs(np.asarray(data[0]).ravel()-(point-CENTRE)/p.XSCALE))),
             'value_solver':info}
    active_x=x.copy();active_x[:2]=[.001,.00015];active_x[8]=p.EMIN+1e-6
    request=np.array([.6,.5,.1,-.1])
    nd,nc,na=p.flow(0.,active_x,eta,request)
    taction=torch.tensor(request[None],dtype=torch.float64,requires_grad=True)
    td,tc=learning.flow(torch.zeros(1,dtype=torch.float64),torch.tensor(active_x[None],dtype=torch.float64),te,taction)
    tc.sum().backward()
    results['active_energy_torch_numpy_flow_difference']=float(np.max(abs(td.detach().numpy()[0]-nd)))
    results['active_energy_gradient_finite']=bool(torch.isfinite(taction.grad).all())
    equivalent_errors=[]
    for _ in range(30):
        ex=eta.copy();ex[6]=rng.uniform(-.004,.004);xx=x.copy();xx[0]=rng.uniform(-.001,.001)
        equivalent_errors.append(float(np.max(abs(p.wrench(0,xx,ex)[:2]-p.wrench(0,xx,p.current_equivalent(xx,ex))[:2]))))
    results['current_equivalent_FnFt_max_error']=max(equivalent_errors)
    assert max(errs)<1e-8 and max(qerrs)<1e-7
    assert abs(float(tv.detach())-val)<1e-7 and results['task_directional_gradient_relative_error']<1e-4
    assert info['feasible'] and info['cost']<info['history'][0]
    assert results['active_energy_torch_numpy_flow_difference']<1e-7 and results['active_energy_gradient_finite']
    assert results['current_equivalent_FnFt_max_error']<1e-10
    print(json.dumps(results,indent=2));Path(__file__).with_name('checks.json').write_text(json.dumps(results,indent=2))

if __name__=='__main__':main()
