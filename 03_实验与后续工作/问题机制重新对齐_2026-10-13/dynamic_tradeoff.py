"""Physical transient contact diagnostic; no imposed optimal gain trajectory.

Direct gain-rate shooting is a demand diagnostic, not the R15 HJB algorithm.
Complex-step batch derivatives differentiate the same RK4 and stage quadrature.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.optimize import minimize, LinearConstraint
import compiled_flow

ROOT = Path(__file__).resolve().parent
CFG = dict(mass=.2, robot_damping=1., force_preload=3., horizon=.8,
           segments=16, integration_step=.005, initial_error=.0015,
           initial_velocity=.008, initial_gains=[50.,3.],
           gain_low=[20.,1.], gain_high=[1200.,30.],
           gain_rate=[3000.,80.], energy_min=.001, energy_reserve=.0015,
           recharge=.3, actuator_max=2., force_min=.2, force_max=5.,
           error_scale=.001, velocity_scale=.025, force_scale=.5,
           command_scale=1., gain_scale=[400.,12.],
           weights=dict(position=1., velocity=.05, force=.3, command=.05,
                        compliance=.3, rate=.002, terminal=.05))


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding='utf8')


class Problem:
    def __init__(self, theta=300., horizon=None, softness=None, dt=None,
                 initial_error=None, initial_velocity=None, overrides=None):
        self.cfg = dict(CFG)
        self.cfg['weights'] = dict(CFG['weights'])
        if overrides is not None: self.cfg.update(overrides)
        if horizon is not None: self.cfg['horizon'] = horizon
        if softness is not None: self.cfg['weights']['compliance'] = softness
        if initial_error is not None: self.cfg['initial_error'] = initial_error
        if initial_velocity is not None: self.cfg['initial_velocity'] = initial_velocity
        self.theta = float(theta)
        self.n = self.cfg['segments']; self.h = self.cfg['horizon']/self.n
        self.sub = round(self.h/(dt or self.cfg['integration_step']))
        self.dt = self.h/self.sub
        self.rate = np.array(self.cfg['gain_rate'])
        self.lo = np.array(self.cfg['gain_low']); self.hi = np.array(self.cfg['gain_high'])
        self.pref = np.array(self.cfg['initial_gains'])
        self.gscale = np.array(self.cfg['gain_scale'])
        self.y0 = np.r_[self.cfg['initial_error'],self.cfg['initial_velocity'],
                         self.pref,self.cfg['energy_min']+self.cfg['energy_reserve']]
        c=self.cfg;w=c['weights']
        self.params=np.array([c['mass'],c['robot_damping'],c['force_preload'],c['recharge'],
          c['energy_min'],c['energy_reserve'],c['actuator_max'],c['force_min'],c['force_max'],
          c['error_scale'],c['velocity_scale'],c['force_scale'],c['command_scale'],
          *self.pref,*self.gscale,w['position'],w['velocity'],w['force'],w['command'],
          w['compliance'],w['rate'],w['terminal'],*self.rate])
        mat = np.kron(np.tril(np.ones((self.n,self.n))),np.eye(2))
        self.gmap = mat*np.tile(self.h*self.rate,self.n)[None,:]
        self.linear = LinearConstraint(self.gmap,
            np.tile(self.lo-self.pref,self.n),np.tile(self.hi-self.pref,self.n))
        self.cached = None

    def rhs(self, y, z):
        c = self.cfg; e,v,k,d,energy = y.T; u = z*self.rate
        f = -k*e-d*v
        dy = np.column_stack((v,(f-self.theta*e-c['robot_damping']*v)/c['mass'],
                              u[:,0],u[:,1],
                              c['recharge']*(c['robot_damping']+d)*v*v-.5*e*e*u[:,0]))
        w = c['weights']
        ell = w['position']*(e/c['error_scale'])**2
        ell += w['force']*(self.theta*e/c['force_scale'])**2
        ell += w['velocity']*(v/c['velocity_scale'])**2
        ell += w['command']*(f/c['command_scale'])**2
        ell += w['compliance']*np.sum(((y[:,2:4]-self.pref)/self.gscale)**2,axis=1)
        ell += w['rate']*np.sum(z*z,axis=1)
        return dy, ell

    def constraints(self,y):
        c=self.cfg; e,v,k,d,E=y.T; fr=-k*e-d*v
        fn=c['force_preload']-self.theta*e
        return np.column_stack(((E-c['energy_min'])/c['energy_reserve'],
             (c['actuator_max']**2-fr*fr)/c['actuator_max']**2,
             (fn-c['force_min'])/c['force_preload'],
             (c['force_max']-fn)/c['force_preload']))

    def rollout(self,z,record=False):
        out=compiled_flow.rollout(np.atleast_2d(z),self.theta,self.params,self.y0,
                                  self.n,self.sub,self.dt)
        return out if record else out[:2]

    def rollout_reference(self,z,record=False):
        z=np.atleast_2d(z); batch=len(z); z=z.reshape(batch,self.n,2)
        y=np.tile(self.y0,(batch,1)).astype(z.dtype)
        J=np.zeros(batch,dtype=z.dtype); margins=[]; trace=[]
        if record: trace.append(np.r_[0.,y[0],0.,0.])
        for j in range(self.n):
            for s in range(self.sub):
                a,la=self.rhs(y,z[:,j]); ya=y+self.dt*a/2
                b,lb=self.rhs(ya,z[:,j]); yb=y+self.dt*b/2
                cc,lc=self.rhs(yb,z[:,j]); yc=y+self.dt*cc
                d,ld=self.rhs(yc,z[:,j])
                # Audit endpoints and RK stages; trajectory certificate not claimed.
                for stage in [y,ya,yb,yc]: margins.append(self.constraints(stage))
                y=y+self.dt*(a+2*b+2*cc+d)/6
                J+=self.dt*(la+2*lb+2*lc+ld)/6
                if record: trace.append(np.r_[(j*self.sub+s+1)*self.dt,y[0],z[0,j]*self.rate])
        margins.append(self.constraints(y))
        c=self.cfg; J+=c['weights']['terminal']*((y[:,0]/c['error_scale'])**2+
                         c['weights']['force']*(self.theta*y[:,0]/c['force_scale'])**2)
        out=(J,np.concatenate(margins,axis=1))
        return (*out,np.array(trace)) if record else out

    def evaluate(self,z):
        if self.cached is not None and np.array_equal(z,self.cached[0]): return self.cached[1]
        n=len(z); batch=np.tile(np.asarray(z,dtype=complex),(n+1,1))
        batch[1:]+=1j*1e-25*np.eye(n)
        J,c=self.rollout(batch)
        result=(float(J[0].real),J[1:].imag/1e-25,c[0].real,c[1:].imag.T/1e-25)
        self.cached=(z.copy(),result)
        return result

    def solve(self,start=None,maxiter=160):
        z=np.zeros(2*self.n) if start is None else np.asarray(start)
        t=time.perf_counter()
        res=minimize(lambda x:self.evaluate(x)[:2],z,jac=True,method='SLSQP',
          bounds=[(-1,1)]*len(z),constraints=[self.linear,dict(type='ineq',
          fun=lambda x:self.evaluate(x)[2],jac=lambda x:self.evaluate(x)[3])],
          options=dict(maxiter=maxiter,ftol=1e-10))
        elapsed=time.perf_counter()-t
        return res,elapsed

    def metrics(self,z):
        J,c,tr=self.rollout(z,True); e,v,k,d,E=tr[:,1:6].T
        fr=-k*e-d*v; fn=self.cfg['force_preload']-self.theta*e
        return dict(cost=float(J[0]),position_rmse_mm=float(np.sqrt(np.mean(e**2))*1e3),
             force_rmse_N=float(np.sqrt(np.mean((fn-self.cfg['force_preload'])**2))),
             peak_actuator_N=float(np.max(abs(fr))),min_contact_N=float(np.min(fn)),
             min_energy_J=float(np.min(E)),max_gain_K=float(k.max()),max_gain_D=float(d.max()),
             final_gains=tr[-1,3:5].tolist(),minimum_margin=float(c.min()),
             max_gain_violation=float(max(np.max(self.lo-tr[:,3:5]),np.max(tr[:,3:5]-self.hi),0.)))

    def policy(self,pars,mode):
        """Common physical interface; state feedback applies at each segment."""
        y=self.y0.copy(); zz=[]
        for j in range(self.n):
            if mode=='fixed': target=np.asarray(pars)
            elif mode=='lqr':
                q=1/self.cfg['error_scale']**2+self.cfg['weights']['force']*(self.theta/self.cfg['force_scale'])**2
                r=self.cfg['weights']['command']/self.cfg['command_scale']**2
                kt=np.sqrt(self.theta**2+q/r)-self.theta
                dd=np.sqrt(self.cfg['robot_damping']**2+2*self.cfg['mass']*kt)-self.cfg['robot_damping']
                target=np.array([kt,dd])*np.array(pars)
            elif mode=='state':
                # Six freely tuned parameters: max K/D deviations, e/v gates,
                # approach/recession term and gain target response time.
                ak,ad,se,sv,bias,tau=pars
                gate=(y[0]/se)**2+(y[1]/sv)**2
                gate=gate/(1+gate)
                sign=np.tanh(bias*y[0]*y[1]/(se*sv))
                target=self.pref+np.array([ak,ad])*np.clip(gate*(1+sign)/2,0,1)
            target=np.clip(target,self.lo,self.hi)
            tau=pars[-1] if mode=='state' else .04
            z=np.clip((target-y[2:4])/np.maximum(tau,self.h)/self.rate,-1,1)
            z=np.clip(z,(self.lo-y[2:4])/(self.h*self.rate),(self.hi-y[2:4])/(self.h*self.rate))
            # A nonempty rate interval always contains zero. Worst-case current
            # energy reserve is allocated to one segment; common viable cap.
            reserve=max(y[4]-self.cfg['energy_min'],0.)
            needed=.5*y[0]**2*max(z[0]*self.rate[0],0.)*self.h
            if needed>reserve and needed>0: z[0]*=reserve/needed
            zz.append(z)
            for s in range(self.sub):
                y=compiled_flow.step(y,z,self.theta,self.params,self.dt)[0]
        return np.array(zz).ravel()


def tune(p,mode,maxiter=75):
    if mode=='fixed': bounds=[(20.,1200.),(1.,30.)];start=[200.,8.]
    elif mode=='lqr': bounds=[(.01,1.5),(.01,1.5)];start=[.2,.3]
    else:
        # Transform to normalized coordinates; avoid false failure from scale.
        bounds=[(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)];start=[.5,.5,.5,.5,.5,.5]
    def physical(x):
        if mode!='state':return np.asarray(x)
        return np.array([1180*x[0],27*x[1],10**(-4+2*x[2]),
                         10**(-3+2*x[3]),-4+8*x[4],.025+.225*x[5]])
    def objective(x):
        z=p.policy(physical(x),mode);J,c=p.rollout(z)
        violation=np.minimum(c,0.)
        return float(J[0]+1e4*np.sum(violation**2))
    runs=[]
    starts=[start]
    if mode=='state':starts += [[.2,.8,.2,.4,.5,.1],[.8,.2,.4,.2,.5,.1]]
    if mode=='fixed':starts += [[50.,3.],[600.,15.]]
    if mode=='lqr': starts += [[.02,.12],[.1,.2]]
    for x in starts:
        res=minimize(objective,x,method='L-BFGS-B',bounds=bounds,
                     options=dict(maxiter=maxiter,ftol=1e-11))
        runs.append(res)
    best=min(runs,key=lambda r:r.fun);pars=physical(best.x);z=p.policy(pars,mode)
    return dict(parameters=pars.tolist(),optimizer_success=bool(best.success),
           message=str(best.message),metrics=p.metrics(z)),z


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--theta',type=float,default=300.)
    parser.add_argument('--softness',type=float,default=.3)
    parser.add_argument('--name',default='pilot');args=parser.parse_args()
    out=ROOT/args.name;out.mkdir(exist_ok=True)
    p=Problem(args.theta,softness=args.softness);dump(out/'config.json',dict(p.cfg,theta=p.theta))
    p.evaluate(np.zeros(2*p.n))  # compilation is setup, separate from solve timings
    results={};starts=[np.zeros(2*p.n)]
    for mode in ['fixed','lqr','state']:
        r,z=tune(p,mode);results[mode]=r;starts.append(z)
        np.save(out/(mode+'_rates.npy'),z)
        tr=p.rollout(z,True)[2];np.savetxt(out/(mode+'_trace.csv'),tr,delimiter=',',
               header='t,e,v,K,D,E,uK,uD',comments='')
        print(mode,json.dumps(r),flush=True);dump(out/'results.json',results)
    sols=[]
    for j,start in enumerate(starts):
        res,seconds=p.solve(start);r=dict(success=bool(res.success),message=str(res.message),
           iterations=int(res.nit),seconds=seconds,metrics=p.metrics(res.x))
        sols.append(r);np.save(out/f'opt_{j}_rates.npy',res.x)
        print('opt',j,json.dumps(r),flush=True)
        results['optimization_runs']=sols;dump(out/'results.json',results)
    feasible=[(i,r) for i,r in enumerate(sols) if r['success'] and r['metrics']['minimum_margin']>=-1e-6]
    if not feasible:raise RuntimeError('No accepted common optimizer result')
    i,r=min(feasible,key=lambda a:a[1]['metrics']['cost']);z=np.load(out/f'opt_{i}_rates.npy')
    tr=p.rollout(z,True)[2];np.savetxt(out/'optimal_trace.csv',tr,delimiter=',',
               header='t,e,v,K,D,E,uK,uD',comments='')
    results['best_optimization']=dict(index=i,**r)
    results['relative_simple_excess']={k:results[k]['metrics']['cost']/r['metrics']['cost']-1 for k in ['fixed','lqr','state']}
    dump(out/'results.json',results);print('relative_simple_excess',results['relative_simple_excess'],flush=True)


if __name__=='__main__': main()
