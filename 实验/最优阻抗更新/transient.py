"""Finite-horizon impedance-rate screen after a model change.

Expected physical cost of mean/covariance dynamics, actual four gain states
and gain-rate limits. No preferred gain trajectory, realized future noise or
ultrasound proxy. Energy tank and hard-force constraints are NOT implemented.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from numba import njit
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import minimize,LinearConstraint
import model as m

HERE=Path(__file__).resolve().parent
SEGMENTS,HORIZON,SUB=12,.30,25
RATE=np.array([800.,800.,80.,80.])/m.SCALE


@njit(cache=True)
def flow(mu,p,g,a0,ais,lis,c0,diff):
    a=a0+g[0]*ais[0]+g[1]*ais[1]+g[2]*ais[2]+g[3]*ais[3]
    l=g[0]*lis[0]+g[1]*lis[1]+g[2]*lis[2]+g[3]*lis[3]
    c=c0+l.T@l
    return a@mu,a@p+p@a.T+diff,mu@c@mu+np.trace(c@p)


@njit(cache=True)
def rollout(batch,g0,mu0,p0,a0,ais,lis,c0,diff,sub):
    h=HORIZON/SEGMENTS;dt=h/sub
    values=np.empty(len(batch),dtype=batch.dtype)
    trace=np.empty((SEGMENTS+1,5),dtype=batch.dtype)
    for b in range(len(batch)):
        g=g0+batch[b,0]*0;mu=mu0+batch[b,0]*0;p=p0+batch[b,0]*0;value=batch[b,0]*0
        if b==0:trace[0,0]=0.;trace[0,1:]=g
        for phase in range(SEGMENTS):
            rate=batch[b,4*phase:4*phase+4]*RATE
            for step in range(sub):
                am,ap,aj=flow(mu,p,g,a0,ais,lis,c0,diff)
                bm,bp,bj=flow(mu+dt*am/2,p+dt*ap/2,g+dt*rate/2,a0,ais,lis,c0,diff)
                cm,cp,cj=flow(mu+dt*bm/2,p+dt*bp/2,g+dt*rate/2,a0,ais,lis,c0,diff)
                dm,dp,dj=flow(mu+dt*cm,p+dt*cp,g+dt*rate,a0,ais,lis,c0,diff)
                mu=mu+dt*(am+2*bm+2*cm+dm)/6
                p=p+dt*(ap+2*bp+2*cp+dp)/6
                value+=dt*(aj+2*bj+2*cj+dj)/6
                g=g+dt*rate
            if b==0:trace[phase+1,0]=(phase+1)*h;trace[phase+1,1:]=g
        values[b]=value
    return values,trace


class Transient:
    def __init__(self,problem,g0,mu0,p0,blocks=12,sub=SUB):
        self.problem=problem;self.g0=g0;self.mu0=mu0;self.p0=p0
        self.blocks=blocks;self.sub=sub;self.cached=None
        self.a0=problem.a0;self.ais=np.array(problem.ais);self.lis=np.array(problem.lis)
        self.c0=problem.basecost;self.diff=problem.diffusion
        self.expand=np.kron(np.repeat(np.eye(blocks),SEGMENTS//blocks,axis=0),np.eye(4))
        self.gmap=np.kron(np.tril(np.ones((SEGMENTS,SEGMENTS))),np.diag(RATE*HORIZON/SEGMENTS))@self.expand
        self.constraint=LinearConstraint(self.gmap,np.tile(m.LOW-g0,SEGMENTS),np.tile(m.HIGH-g0,SEGMENTS))

    def run(self,z,sub=None):
        zz=np.atleast_2d(z)@self.expand.T
        return rollout(np.ascontiguousarray(zz),self.g0,self.mu0,self.p0,self.a0,self.ais,self.lis,self.c0,self.diff,self.sub if sub is None else sub)

    def evaluate(self,z):
        if self.cached is not None and np.array_equal(self.cached[0],z):return self.cached[1]
        batch=np.tile(z.astype(complex),(len(z)+1,1));batch[1:]+=1j*1e-23*np.eye(len(z))
        values,_=self.run(batch)
        out=float(values[0].real),values[1:].imag/1e-23
        self.cached=(z.copy(),out)
        return out

    def solve(self,start=None,maxiter=100):
        if start is None:start=np.zeros(4*self.blocks)
        t=time.perf_counter()
        sol=minimize(self.evaluate,start,jac=True,method='SLSQP',bounds=[(-1,1)]*len(start),
            constraints=[self.constraint],options=dict(maxiter=maxiter,ftol=1e-10))
        return sol,time.perf_counter()-t


def ramp(g0,target):
    g=g0.copy();out=[];h=HORIZON/SEGMENTS
    for _ in range(SEGMENTS):
        rate=np.clip((target-g)/(h*RATE),-1,1)
        out.extend(rate);g+=h*RATE*rate
    return np.array(out)


def main():
    protocol=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256((HERE/'model.py').read_bytes()).hexdigest(),
        seed=4901,cases_per_family=4,families=['white','colored_observation_noise'],
        horizon_s=HORIZON,segments=SEGMENTS,integration_substeps=SUB,
        actual_gain_rate=(RATE*m.SCALE).tolist(),
        initial_gains='old-model stationary optimal gains, shared',
        initial_covariance='old-model stationary distribution, shared',
        current_mean_offset_m=.002,current_mean_velocity_offset_m_s=.01,
        information='Same known updated model, belief, noise statistics; future noise realizations unknown to all',
        goal='Expected integral of physical displacement/velocity/contact-force/actuator costs, no gain preference or terminal gain target',
        scope='8 exploratory mechanism cases, not an independent population benchmark or a deployed controller',
        exclusions=['energy state and hard-force constraints','MHE','full robot/image validation'])
    (HERE/'transient_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    records=[];arrays={}
    for family_index,(label,colored,sensed) in enumerate([('white',False,False),('colored_observation_noise',True,True)]):
        rng=np.random.default_rng(4901+family_index);dim=4 if colored else 3
        for case in range(4):
            oldx,newx=rng.uniform(-1,1,(2,dim))
            old=m.Problem(oldx,colored,sensed);problem=m.Problem(newx,colored,sensed)
            g0=old.solve().x;target=problem.solve().x
            a,_,_=old.matrices(g0);p0=solve_continuous_lyapunov(a,-old.diffusion)
            mu0=np.zeros(problem.n);mu0[:2]=rng.uniform(-.002,.002,2);mu0[2:4]=rng.uniform(-.01,.01,2)
            full=Transient(problem,g0,mu0,p0)
            # JIT compilation before timings. Both dtypes share the same code.
            full.run(np.zeros(48));full.evaluate(np.zeros(48))
            zfreeze=np.zeros(48);zstationary=ramp(g0,target);zanalytic=ramp(g0,problem.analytic())
            costs=[float(full.run(z)[0][0]) for z in [zfreeze,zstationary,zanalytic]]
            zs=[zfreeze,zstationary,zanalytic];times=[0.,0.,0.];solmeta=[]
            for blocks in [1,2]:
                low=Transient(problem,g0,mu0,p0,blocks=blocks)
                guess=zstationary.reshape(blocks,SEGMENTS//blocks,4).mean(axis=1).ravel()
                low.evaluate(guess) # exclude compile/cache initialization, not numerical solves
                low.cached=None
                sol,elapsed=low.solve(guess)
                z=low.expand@sol.x;zs.append(z);costs.append(float(full.run(z)[0][0]));times.append(elapsed)
                solmeta.append(dict(blocks=blocks,success=bool(sol.success),iterations=sol.nit,message=str(sol.message)))
            best=int(np.argmin(costs));full.cached=None
            sol,elapsed=full.solve(zs[best],maxiter=150)
            zref=sol.x if sol.fun<=costs[best] else zs[best]
            zs.append(zref);costs.append(float(full.run(zref)[0][0]));times.append(elapsed)
            names=['freeze_old_gains','slew_to_new_stationary_optimum','slew_to_analytic','one_rate_block','two_rate_blocks','12_rate_blocks_reference']
            detailed=[]
            for name,z,j,elapsed in zip(names,zs,costs,times):
                gains=g0+(full.gmap@z).reshape(SEGMENTS,4)
                detailed.append(dict(method=name,J=j,gap_to_reference=j/costs[-1]-1,
                    optimization_ms=elapsed*1000,
                    minimum_gain_margin=float(min((gains-m.LOW).min(),(m.HIGH-gains).min())),
                    maximum_abs_normalized_rate=float(np.abs(z).max())))
            arrays[label+f'_{case}_rates']=np.array(zs)
            arrays[label+f'_{case}_initial']=np.r_[oldx,newx,g0,mu0,p0.ravel()]
            # Two preselected case-0 paths get a 2x integration refinement.
            refine=None
            if case==0:
                fine=[float(full.run(z,sub=2*SUB)[0][0]) for z in zs]
                refine=dict(costs=fine,max_relative_change=max(abs(np.array(fine)/np.array(costs)-1)))
            records.append(dict(family=label,case=case,old_features=oldx.tolist(),new_features=newx.tolist(),
                initial_mean=mu0.tolist(),initial_gains=(g0*m.SCALE).tolist(),stationary_target=(target*m.SCALE).tolist(),
                methods=detailed,low_dimension_solvers=solmeta,
                reference_solver=dict(success=bool(sol.success),iterations=sol.nit,message=str(sol.message)),refinement=refine))
            print(label,case,json.dumps(detailed),flush=True)
    np.savez_compressed(HERE/'transient_paired.npz',**arrays)
    (HERE/'transient_results.json').write_text(json.dumps(dict(protocol=protocol,cases=records),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
