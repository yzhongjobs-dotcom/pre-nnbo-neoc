"""Structured continuous-adjoint gradient for the finite-horizon screen.

Integrates second moments directly. Forward RK4 + backward RK4/Hermite moment
interpolation gives an O(dt^4) gradient approximation; verify against independent
complex-step and DOP853 rather than calling it an exact discrete derivative.
"""
import numpy as np
from numba import njit
from scipy.optimize import minimize
import transient as t


@njit(cache=True)
def matrices(g,a0,ais,lis,c0):
    a=a0+g[0]*ais[0]+g[1]*ais[1]+g[2]*ais[2]+g[3]*ais[3]
    l=g[0]*lis[0]+g[1]*lis[1]+g[2]*lis[2]+g[3]*lis[3]
    return a,c0+l.T@l,l


@njit(cache=True)
def forward(s,g,a0,ais,lis,c0,diff):
    a,c,_=matrices(g,a0,ais,lis,c0)
    return a@s+s@a.T+diff,np.trace(c@s)


@njit(cache=True)
def backward(v,s,g,a0,ais,lis,c0):
    a,c,l=matrices(g,a0,ais,lis,c0)
    gg=np.empty(4)
    for i in range(4):
        ci=lis[i].T@l+l.T@lis[i]
        gg[i]=np.trace(ci@s)+2*np.trace(v@ais[i]@s)
    return a.T@v+v@a+c,gg


@njit(cache=True)
def value_gradient(z,g0,s0,a0,ais,lis,c0,diff,sub):
    n=len(s0);steps=t.SEGMENTS*sub;dt=t.HORIZON/steps
    moments=np.empty((steps+1,n,n));gains=np.empty((steps+1,4))
    moments[0]=s0;gains[0]=g0;value=0.
    for step in range(steps):
        phase=step//sub;rate=z[4*phase:4*phase+4]*t.RATE
        s=moments[step];g=gains[step]
        a,ja=forward(s,g,a0,ais,lis,c0,diff)
        b,jb=forward(s+dt*a/2,g+dt*rate/2,a0,ais,lis,c0,diff)
        c,jc=forward(s+dt*b/2,g+dt*rate/2,a0,ais,lis,c0,diff)
        d,jd=forward(s+dt*c,g+dt*rate,a0,ais,lis,c0,diff)
        moments[step+1]=s+dt*(a+2*b+2*c+d)/6
        gains[step+1]=g+dt*rate
        value+=dt*(ja+2*jb+2*jc+jd)/6
    v=np.zeros((n,n));lam=np.zeros(4);gradient=np.zeros(len(z))
    for rev in range(steps):
        step=steps-rev-1;phase=step//sub
        s0i=moments[step];s1=moments[step+1];g0i=gains[step];g1=gains[step+1]
        f0,_=forward(s0i,g0i,a0,ais,lis,c0,diff)
        f1,_=forward(s1,g1,a0,ais,lis,c0,diff)
        smid=.5*(s0i+s1)+dt*(f0-f1)/8
        gmid=.5*(g0i+g1)
        av,ag=backward(v,s1,g1,a0,ais,lis,c0)
        bv,bg=backward(v+dt*av/2,smid,gmid,a0,ais,lis,c0)
        cv,cg=backward(v+dt*bv/2,smid,gmid,a0,ais,lis,c0)
        dv,dg=backward(v+dt*cv,s0i,g0i,a0,ais,lis,c0)
        # Integral of lambda_g over this gain-rate interval.
        gradient[4*phase:4*phase+4]+=dt*t.RATE*(lam+2*(lam+dt*ag/2)+2*(lam+dt*bg/2)+(lam+dt*cg))/6
        lam+=dt*(ag+2*bg+2*cg+dg)/6
        v+=dt*(av+2*bv+2*cv+dv)/6
    return value,gradient


class AdjointTransient(t.Transient):
    def evaluate(self,z):
        if self.cached is not None and np.array_equal(self.cached[0],z):return self.cached[1]
        value,gradient=value_gradient(self.expand@z,self.g0,self.p0+np.outer(self.mu0,self.mu0),
            self.a0,self.ais,self.lis,self.c0,self.diff,self.sub)
        out=value,self.expand.T@gradient
        self.cached=(z.copy(),out)
        return out


def main():
    from pathlib import Path
    import hashlib,json,time
    import model as m
    here=Path(__file__).resolve().parent
    source=json.loads((here/'transient_results.json').read_text())
    data=np.load(here/'transient_paired.npz');records=[]
    for row in source['cases']:
        label,case=row['family'],row['case'];colored=label!='white';dim=4 if colored else 3
        arr=data[label+f'_{case}_initial'];n=10 if colored else 4
        g0=arr[2*dim:2*dim+4];mu=arr[2*dim+4:2*dim+4+n];p0=arr[-n*n:].reshape(n,n)
        problem=m.Problem(row['new_features'],colored,colored)
        for blocks in [1,2,12]:
            p=AdjointTransient(problem,g0,mu,p0,blocks=blocks)
            # Identical starts to complex-step implementation, not the already
            # converged full solution from the earlier benchmark.
            stationary=data[label+f'_{case}_rates'][1]
            if blocks<12:guess=stationary.reshape(blocks,t.SEGMENTS//blocks,4).mean(axis=1).ravel()
            else:
                # Same best baseline selected in the original full solve.
                baseline=int(np.argmin([r['J'] for r in row['methods'][:-1]]))
                guess=data[label+f'_{case}_rates'][baseline]
            p.evaluate(guess);p.cached=None
            sol,elapsed=p.solve(guess,maxiter=150 if blocks==12 else 100)
            z=p.expand@sol.x
            physical=float(t.Transient(problem,g0,mu,p0).run(z)[0][0])
            reference=row['methods'][{1:3,2:4,12:5}[blocks]]
            records.append(dict(family=label,case=case,blocks=blocks,
                time_ms=elapsed*1000,complex_step_time_ms=reference['optimization_ms'],
                J_in_original_integrator=physical,original_reference_J=reference['J'],
                relative_cost_difference=physical/reference['J']-1,
                success=bool(sol.success),iterations=sol.nit))
        print(label,case,flush=True)
    (here/'adjoint_results.json').write_text(json.dumps(dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        transient_sha256=hashlib.sha256((here/'transient.py').read_bytes()).hexdigest(),
        scope='Same starts, bounds and objective; structured continuous-adjoint derivative, not B0 NEOC',
        cases=records),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
