"""Slow material and repeated surface motion, actual gain-rate control.

Shooting reference and rival controllers use identical material/surface forecasts.
No imposed gain target: simple targets come from the physical equilibrium cost.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from numba import njit
from scipy.optimize import minimize,LinearConstraint

ROOT=Path(__file__).resolve().parent
CFG=dict(mass=.2,robot_damping=1.,preload=3.,horizon=4.,segments=32,dt=.005,
   theta0=300.,theta_span=250.,initial_gains=[50.,3.],gain_low=[20.,1.],gain_high=[1200.,30.],
   gain_rate=[400.,30.],energy_min=.001,reserve=.00005,recharge=.3,energy_alpha=10.,
   actuator_increment_max=2.,force_min=.2,force_max=5.,error_scale=.001,
   velocity_scale=.025,force_scale=.5,command_scale=1.,gain_scale=[400.,12.],
   w_position=1.,w_velocity=.05,w_force=.3,w_command=.05,w_compliance=.3,w_rate=.002,
   terminal=.05,planning_guard=3e-5,pulses=[[.7,.6,.0015],[2.7,.6,.003]])


def dump(path,data):path.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf8')


@njit(cache=True)
def environment(t,p,pulses):
    phase=max(0.,min(t/p[3],1.))
    th=p[4]+p[5]*(1-np.cos(np.pi*phase))/2
    if len(p)>31 and p[31]>0:
        th=max(50.,p[4]+p[33]*(t-p[32]))
    gam=0.
    for row in pulses:
        if row[0]<t<row[0]+row[1]:gam+=row[2]*np.sin(np.pi*(t-row[0])/row[1])**2
    return th,gam


@njit(cache=True)
def flow(y,z,t,p,pulses):
    e,v,k,d,E=y;uk=z[0]*p[6];ud=z[1]*p[7];th,gam=environment(t,p,pulses)
    force_error=th*(gam-e);fr=-k*e-d*v
    dy=np.array([v,(fr+force_error-p[1]*v)/p[0],uk,ud,
                  p[10]*(p[1]+d)*v*v-.5*e*e*uk])
    ell=p[18]*(e/p[14])**2+p[19]*(v/p[15])**2+p[20]*(force_error/p[16])**2
    ell+=p[21]*(fr/p[17])**2+p[22]*(((k-p[24])/p[26])**2+((d-p[25])/p[27])**2)
    ell+=p[23]*(z[0]*z[0]+z[1]*z[1])
    return dy,ell


@njit(cache=True)
def step(y,z,t,h,p,pulses):
    a,la=flow(y,z,t,p,pulses);ya=y+h*a/2
    b,lb=flow(ya,z,t+h/2,p,pulses);yb=y+h*b/2
    c,lc=flow(yb,z,t+h/2,p,pulses);yc=y+h*c
    d,ld=flow(yc,z,t+h,p,pulses)
    return y+h*(a+2*b+2*c+d)/6,h*(la+2*lb+2*lc+ld)/6,ya,yb,yc


@njit(cache=True)
def margins(y,z,t,p,pulses):
    e,v,k,d,E=y;th,gam=environment(t,p,pulses);fr=-k*e-d*v
    fn=p[2]+th*(gam-e)
    power=p[10]*(p[1]+d)*v*v+p[11]*(E-p[8])-.5*e*e*z[0]*p[6]
    scale=p[11]*p[9]+.5*p[14]**2*p[6]
    return np.array([(E-p[8])/p[9],(p[12]**2-fr*fr)/p[12]**2,
                 (fn-p[13])/p[2],(p[28]-fn)/p[2],power/scale])


@njit(cache=True)
def rollout(z,y0,t0,n,sub,dt,p,pulses):
    batch=z.shape[0];nsteps=n*sub
    J=np.zeros(batch,dtype=z.dtype);c=np.empty((batch,5*(4*nsteps+1)),dtype=z.dtype)
    tr=np.empty((nsteps+1,10),dtype=z.dtype)
    for b in range(batch):
        y=y0.astype(z.dtype);idx=0
        if b==0:
            th,gam=environment(t0,p,pulses)
            tr[0,0]=t0;tr[0,1:6]=y;tr[0,6:8]=0.;tr[0,8]=th;tr[0,9]=gam
        for j in range(n):
            u=z[b,2*j:2*j+2]
            for s in range(sub):
                t=t0+(j*sub+s)*dt
                yn,ell,ya,yb,yc=step(y,u,t,dt,p,pulses)
                c[b,idx:idx+5]=margins(y,u,t,p,pulses);idx+=5
                c[b,idx:idx+5]=margins(ya,u,t+dt/2,p,pulses);idx+=5
                c[b,idx:idx+5]=margins(yb,u,t+dt/2,p,pulses);idx+=5
                c[b,idx:idx+5]=margins(yc,u,t+dt,p,pulses);idx+=5
                y=yn;J[b]+=ell
                if b==0:
                    row=j*sub+s+1;th,gam=environment(t+dt,p,pulses)
                    tr[row,0]=t+dt;tr[row,1:6]=y;tr[row,6]=u[0]*p[6];tr[row,7]=u[1]*p[7]
                    tr[row,8]=th;tr[row,9]=gam
        c[b,idx:idx+5]=margins(y,z[b,-2:],t0+nsteps*dt,p,pulses)
        th,gam=environment(t0+nsteps*dt,p,pulses)
        J[b]+=p[29]*((y[0]/p[14])**2+p[20]*(th*(gam-y[0])/p[16])**2)
    return J,c,tr


@njit(cache=True)
def phase_minimum(y,z,t,sub,dt,p,pulses):
    m=1e99
    for s in range(sub):
        now=t+s*dt;yn,ell,ya,yb,yc=step(y,z,now,dt,p,pulses)
        for state,phase in ((y,now),(ya,now+dt/2),(yb,now+dt/2),(yc,now+dt)):
            for a in margins(state,z,phase,p,pulses):m=min(m,a)
        y=yn
    for a in margins(y,z,t+sub*dt,p,pulses):m=min(m,a)
    return m


@njit(cache=True)
def safe_rate(y,target,t,h,sub,dt,p,pulses,lo,hi):
    z=np.empty(2)
    for i in range(2):
        z[i]=max(-1.,min((target[i]-y[2+i])/(h*p[6+i]),1.))
        z[i]=max((lo[i]-y[2+i])/(h*p[6+i]),min(z[i],(hi[i]-y[2+i])/(h*p[6+i])))
    if phase_minimum(y,z,t,sub,dt,p,pulses)<p[30]:
        lower=z.copy();lower[0]=min(z[0],0.)
        if phase_minimum(y,lower,t,sub,dt,p,pulses)>=p[30]:
            high=z[0];low=lower[0]
            for _ in range(25):
                trial=z.copy();trial[0]=(high+low)/2
                if phase_minimum(y,trial,t,sub,dt,p,pulses)>=p[30]:low=trial[0]
                else:high=trial[0]
            z[0]=low
    return z


class Problem:
    def __init__(self,reserve=None,segments=None,theta_span=None,dt=None):
        self.cfg=dict(CFG)
        if reserve is not None:self.cfg['reserve']=reserve
        if segments is not None:self.cfg['segments']=segments
        if theta_span is not None:self.cfg['theta_span']=theta_span
        if dt is not None:self.cfg['dt']=dt
        c=self.cfg;self.n=c['segments'];self.h=c['horizon']/self.n
        self.sub=round(self.h/c['dt']);self.dt=self.h/self.sub
        self.rate=np.array(c['gain_rate']);self.lo=np.array(c['gain_low']);self.hi=np.array(c['gain_high'])
        self.pref=np.array(c['initial_gains']);self.pulses=np.array(c['pulses'])
        self.y0=np.r_[0.,0.,self.pref,c['energy_min']+c['reserve']]
        self.p=np.array([c['mass'],c['robot_damping'],c['preload'],c['horizon'],c['theta0'],c['theta_span'],
          *self.rate,c['energy_min'],c['reserve'],c['recharge'],c['energy_alpha'],
          c['actuator_increment_max'],c['force_min'],c['error_scale'],c['velocity_scale'],c['force_scale'],c['command_scale'],
          c['w_position'],c['w_velocity'],c['w_force'],c['w_command'],c['w_compliance'],c['w_rate'],
          *self.pref,*c['gain_scale'],c['force_max'],c['terminal'],c['planning_guard'],0.,0.,0.])
        self.gmap=np.kron(np.tril(np.ones((self.n,self.n))),np.eye(2))*np.tile(self.h*self.rate,self.n)[None,:]
        self.linear=LinearConstraint(self.gmap,np.tile(self.lo-self.pref,self.n),np.tile(self.hi-self.pref,self.n))
        self.cache=None

    def run(self,z,record=False):
        out=rollout(np.atleast_2d(z),self.y0,0.,self.n,self.sub,self.dt,self.p,self.pulses)
        return out if record else out[:2]

    def evaluate(self,z):
        if self.cache is not None and np.array_equal(z,self.cache[0]):return self.cache[1]
        n=len(z);b=np.tile(z.astype(complex),(n+1,1));b[1:]+=1j*1e-25*np.eye(n)
        J,c=self.run(b)
        r=(float(J[0].real),J[1:].imag/1e-25,c[0].real,c[1:].imag.T/1e-25)
        self.cache=(z.copy(),r);return r

    def solve(self,start=None,maxiter=240):
        start=np.zeros(2*self.n) if start is None else start
        tic=time.perf_counter()
        res=minimize(lambda z:self.evaluate(z)[:2],start,jac=True,method='SLSQP',bounds=[(-1,1)]*len(start),
           constraints=[self.linear,dict(type='ineq',fun=lambda z:self.evaluate(z)[2]-self.cfg['planning_guard'],jac=lambda z:self.evaluate(z)[3])],
           options=dict(maxiter=maxiter,ftol=2e-9))
        return res,time.perf_counter()-tic

    def metrics(self,z):
        J,c,tr=self.run(z,True);e,v,k,d,E=tr[:,1:6].T
        fn=self.cfg['preload']+tr[:,8]*(tr[:,9]-e);fr=-k*e-d*v
        return dict(cost=float(J[0]),position_rms_mm=float(np.sqrt(np.trapz(e*e,tr[:,0])/self.cfg['horizon'])*1e3),
                    force_rms_N=float(np.sqrt(np.trapz((fn-self.cfg['preload'])**2,tr[:,0])/self.cfg['horizon'])),
                    minimum_margin=float(c.min()),min_energy=float(E.min()),max_increment_force=float(abs(fr).max()),
                    min_contact=float(fn.min()),max_contact=float(fn.max()),max_K=float(k.max()),max_D=float(d.max()))

    def targets(self,kind,params):
        # Static targets minimize the same force/position/compliance cost.
        times=np.arange(self.n)*self.h
        instantaneous=[]
        for t in times:
            th,gam=environment(t,self.p,self.pulses)
            # Unique stationary derivative in this box, independently checked
            # against scalar optimization in the validation script.
            qp=self.p[18]/self.p[14]**2
            qf=self.p[20]/self.p[16]**2+self.p[21]/self.p[17]**2
            def gradient(k):
                return (th*gam)**2*(-qp+qf*th*k)/(th+k)**3+self.p[22]*(k-self.pref[0])/self.p[26]**2
            low=self.lo[0];high=self.hi[0]
            if gradient(low)>=0:k=low
            elif gradient(high)<=0:k=high
            else:
                for _ in range(45):
                    middle=(low+high)/2
                    if gradient(middle)>0:high=middle
                    else:low=middle
                k=(low+high)/2
            instantaneous.append(k)
        instantaneous=np.array(instantaneous)
        rows=[]
        for j,t in enumerate(times):
            if kind=='fixed':k,d=params[:2]
            else:
                k=instantaneous[j]
                if kind=='preview':
                    lead=params[1]
                    later=(times>=t)&(times<=t+lead)
                    k=max(k,float(np.max(instantaneous[later]-self.rate[0]*(times[later]-t))))
                th,gam=environment(t,self.p,self.pulses)
                d=2*params[0]*np.sqrt(self.cfg['mass']*(th+k))-self.cfg['robot_damping']
            rows.append(np.clip([k,d],self.lo,self.hi))
        return np.array(rows)

    def policy(self,targets):
        y=self.y0.copy();zs=[];control_seconds=0.
        for j in range(self.n):
            tic=time.perf_counter()
            z=safe_rate(y,targets[j],j*self.h,self.h,self.sub,self.dt,self.p,self.pulses,self.lo,self.hi)
            control_seconds+=time.perf_counter()-tic;zs.append(z)
            for s in range(self.sub):y=step(y,z,j*self.h+s*self.dt,self.dt,self.p,self.pulses)[0]
        return np.array(zs).ravel(),control_seconds


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reserve',type=float,default=.00005)
    parser.add_argument('--name',default='pilot');args=parser.parse_args()
    out=ROOT/args.name;out.mkdir(exist_ok=True);p=Problem(reserve=args.reserve)
    dump(out/'config.json',p.cfg);p.evaluate(np.zeros(2*p.n))
    records={};seeds=[np.zeros(2*p.n)]
    choices=[('fixed',[50.,3.]),('fixed',[200.,8.]),('fixed',[500.,12.])]
    choices+=[(k,[zeta,lead]) for k in ['quasi','preview'] for zeta in [.15,.35,.6] for lead in ([0.] if k=='quasi' else [.5,1.5,4.])]
    for idx,(kind,pars) in enumerate(choices):
        targets=p.targets(kind,pars);z,seconds=p.policy(targets);m=p.metrics(z)
        records[str(idx)]=dict(kind=kind,parameters=pars,metrics=m,decision_seconds=seconds,rates=z.tolist())
        print(kind,pars,m['cost'],m['minimum_margin'],flush=True)
    feasible=[r for r in records.values() if r['metrics']['minimum_margin']>=-1e-8]
    if not feasible:raise RuntimeError('No feasible rival incumbent')
    best=min(feasible,key=lambda r:r['metrics']['cost']);seeds.append(np.array(best['rates']))
    results=dict(rivals=records,best_rival=best,optimization_runs=[]);dump(out/'results.json',results)
    for i,start in enumerate(seeds):
        res,seconds=p.solve(start);r=dict(success=bool(res.success),message=str(res.message),iterations=int(res.nit),
            seconds=seconds,metrics=p.metrics(res.x),rates=res.x.tolist())
        results['optimization_runs'].append(r);dump(out/'results.json',results)
        np.save(out/f'opt_{i}_rates.npy',res.x)
        np.savetxt(out/f'opt_{i}_trace.csv',p.run(res.x,True)[2],delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
        print('optimization',i,json.dumps({k:v for k,v in r.items() if k!='rates'}),flush=True)


if __name__=='__main__':main()
