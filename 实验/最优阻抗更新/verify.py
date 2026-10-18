"""Independent covariance/frequency, derivative, actor and transient checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad,solve_ivp
import model as m

HERE=Path(__file__).resolve().parent


def independent(x,g,colored,sensed,with_covariance=True):
    # Rebuild from physical matrices, not Problem.matrices/evaluate.
    kn=300*3**((x[0]+1)/2);ratio=.20+.15*(x[1]+1);cross=.60*x[2]
    ke=np.array([[kn,cross*kn*np.sqrt(ratio)],[cross*kn*np.sqrt(ratio),ratio*kn]])
    mass=np.diag([1.,.25]);be=np.diag([2.,1.]);base=6 if colored else 4;n=base+4*sensed
    k,d=np.diag(g[:2]*[500,500]),np.diag(g[2:]*[40,40])
    a=np.zeros((n,n));a[:2,2:4]=np.eye(2)
    a[2:4,:2]=np.linalg.solve(mass,-ke-k);a[2:4,2:4]=np.linalg.solve(mass,-be-d)
    l=np.zeros((2,n));l[:,:4]=np.column_stack((k,d));diff=np.zeros((n,n))
    if colored:
        rate=4*10**((x[3]+1)/2)
        a[2:4,4:6]=np.linalg.inv(mass);a[4:6,4:6]=-rate*np.eye(2)
        diff[4:6,4:6]=2*rate*np.array([[.04,.012],[.012,.04]])
    else:
        diff[2:4,2:4]=np.linalg.solve(mass,np.array([[.0009,.00027],[.00027,.0009]]))@np.linalg.inv(mass)
    if sensed:
        rate=60*np.pi;a[base:,base:]=-rate*np.eye(4)
        a[2:4,base:]=-np.linalg.solve(mass,np.column_stack((k,d)))
        l[:,base:]=np.column_stack((k,d))
        diff[base:,base:]=2*rate*np.diag(np.array([50e-6,50e-6,.002,.002])**2)
    cf=np.column_stack((ke,be));c=l.T@l
    c[:4,:4]+=cf.T@cf/.7**2
    c[:2,:2]+=np.eye(2)/.003**2;c[2:4,2:4]+=np.eye(2)/.03**2
    if with_covariance:
        operator=np.kron(np.eye(n),a)+np.kron(a,np.eye(n))
        p=np.linalg.solve(operator,-diff.ravel(order='F')).reshape(n,n,order='F')
        return float(np.trace(c@p)),a,c,diff,p
    return None,a,c,diff,None


def main():
    rng=np.random.default_rng(5901);derivatives=[];covariance=[];frequency=[]
    for colored,sensed in [(False,False),(True,False),(True,True)]:
        for i in range(6):
            x=rng.uniform(-.8,.8,4 if colored else 3);p=m.Problem(x,colored,sensed)
            g=p.analytic()*.9+m.LOW*.1
            j,gr,h=p.evaluate(g,True);eps=1e-5;eye=np.eye(4)
            fd=np.array([(p.evaluate(g+eps*v)[0]-p.evaluate(g-eps*v)[0])/(2*eps) for v in eye])
            hf=np.column_stack([(p.evaluate(g+eps*v)[1]-p.evaluate(g-eps*v)[1])/(2*eps) for v in eye])
            derivatives.append(dict(colored=colored,sensed=sensed,gradient_max_abs=float(np.max(abs(gr-fd))),hessian_max_abs=float(np.max(abs(h-hf)))))
            jj,a,c,q,cov=independent(x,g,colored,sensed)
            covariance.append(dict(colored=colored,sensed=sensed,relative_J_difference=abs(jj/j-1),minimum_covariance_eigenvalue=float(np.linalg.eigvalsh(cov).min())))
            if i==0:
                def spectrum(w):
                    transfer=np.linalg.inv(1j*w*np.eye(len(a))-a)
                    return float(np.trace(c@transfer@q@transfer.conj().T).real/np.pi)
                integral,err=quad(spectrum,0,np.inf,epsabs=1e-10,epsrel=1e-9,limit=200)
                frequency.append(dict(colored=colored,sensed=sensed,relative_J_difference=abs(integral/j-1),quadrature_error=err))
    exact=[]
    for knfeature in np.linspace(-1,1,9):
        p=m.Problem([knfeature,.1,0.]);g=p.analytic();j,gr=p.evaluate(g)
        exact.append(dict(feature=float(knfeature),scalar_formula_vs_numerical_J_difference=abs(j-p.solve().fun),
                          scalar_formula_vs_CARE_J_difference=abs(j-p.care_reference()[0]),
                          projected_gradient=float(np.max(abs(m.projected_gradient(g,gr))))))
    import torch
    from torch import nn
    data=np.load(HERE/'paired.npz');actors=[]
    for label in ['white','colored','colored_observation_noise']:
        for seed,method in [(11,3),(22,4),(33,5)]:
            saved=torch.load(HERE/f'actor_{label}_{seed}.pt',weights_only=True)
            net=nn.Sequential(nn.Linear(saved['input_dim'],64),nn.Tanh(),nn.Linear(64,64),nn.Tanh(),nn.Linear(64,4)).double()
            net.load_state_dict(saved['state_dict']);xx=data[label+'_features'][-96:]
            with torch.no_grad():pred=net(torch.tensor(xx[:3])).numpy()
            gg=np.clip(np.exp(pred*np.array(saved['sd'])+np.array(saved['mean'])),m.LOW,m.HIGH)
            actors.append(dict(family=label,seed=seed,max_gain_difference=float(np.max(abs(gg-data[label+'_test_gains'][:3,method])))))
    results=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        covariance=covariance,frequency=frequency,derivatives=derivatives,scalar_exact_cases=exact,actor_reload=actors)
    if (HERE/'transient_results.json').exists():
        import transient as t
        import adjoint as ad
        raw=json.loads((HERE/'transient_results.json').read_text());arrs=np.load(HERE/'transient_paired.npz')
        checks=[];gradient_checks=[]
        for row in [r for r in raw['cases'] if r['case']==0]:
            label=row['family'];colored=label!='white';dim=4 if colored else 3;n=10 if colored else 4
            arr=arrs[label+'_0_initial'];g0=arr[2*dim:2*dim+4];mu0=arr[2*dim+4:2*dim+4+n];p0=arr[-n*n:].reshape(n,n)
            x=row['new_features'];pr=m.Problem(x,colored,colored)
            tp=t.Transient(pr,g0,mu0,p0);ap=ad.AdjointTransient(pr,g0,mu0,p0)
            z=arrs[label+'_0_rates'][-1];value,grad=tp.evaluate(z);adj,ag=ap.evaluate(z)
            gradient_checks.append(dict(family=label,relative_value_difference=abs(adj/value-1),max_abs_gradient_difference=float(np.max(abs(grad-ag)))))
            pairs=list(zip(row['methods'],arrs[label+'_0_rates']))
            if (HERE/'bounded_results.json').exists():
                bounded=json.loads((HERE/'bounded_results.json').read_text())
                selected=next(r for r in bounded['cases'] if r['family']==label and r['case']==0)
                zz=np.repeat(np.array(selected['rate_blocks']).reshape(2,4),6,axis=0).ravel()
                pairs.append((dict(method='bounded_structured_steps',J=selected['history'][-1]['J']),zz))
            for method,z in pairs:
                state=np.r_[mu0,p0.ravel(),0.];horizon=t.HORIZON/t.SEGMENTS;g=g0.copy()
                for phase in range(t.SEGMENTS):
                    rate=z[phase*4:phase*4+4]*t.RATE
                    def rhs(tt,y):
                        _,a,c,q,_=independent(x,g+rate*tt,colored,colored,False)
                        mu=y[:n];cov=y[n:n+n*n].reshape(n,n)
                        return np.r_[a@mu,(a@cov+cov@a.T+q).ravel(),mu@c@mu+np.trace(c@cov)]
                    sol=solve_ivp(rhs,(0,horizon),state,method='DOP853',rtol=1e-10,atol=1e-13)
                    state=sol.y[:,-1];g+=horizon*rate
                checks.append(dict(family=label,method=method['method'],DOP853_J=float(state[-1]),
                    relative_cost_difference=abs(state[-1]/method['J']-1)))
        results['transient_DOP853']=checks;results['adjoint_vs_complex_step']=gradient_checks
        if (HERE/'bounded_results.json').exists():
            bounds=[]
            for row in json.loads((HERE/'bounded_results.json').read_text())['cases']:
                label=row['family'];dim=3 if label=='white' else 4
                arr=arrs[label+f'_{row["case"]}_initial'];g0=arr[2*dim:2*dim+4]
                zz=np.array(row['rate_blocks']).reshape(2,4)
                endpoints=g0+np.cumsum(zz*t.RATE*(t.HORIZON/2),axis=0)
                bounds.append(dict(family=label,case=row['case'],max_abs_rate=float(abs(zz).max()),
                    minimum_gain_margin=float(min((endpoints-m.LOW).min(),(m.HIGH-endpoints).min()))))
            results['bounded_gain_and_rate_checks']=bounds
    (HERE/'independent_checks.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(max_gradient_error=max(r['gradient_max_abs'] for r in derivatives),
        max_hessian_error=max(r['hessian_max_abs'] for r in derivatives),
        max_independent_covariance_relative_error=max(r['relative_J_difference'] for r in covariance),frequency=frequency,
        transient=results.get('transient_DOP853'),adjoint=results.get('adjoint_vs_complex_step')),indent=2))


if __name__=='__main__':main()
