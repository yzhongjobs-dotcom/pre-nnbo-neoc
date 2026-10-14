"""Finite-horizon physical contact carrier with a nodal single value critic.

Numerical repair, not an exact reproduction of the manuscript collocation root.
State = indentation, normal velocity, actual normal stiffness. Damping fixed.
Energy is tracked and its inactive-envelope bound checked in the experiment.
"""
from pathlib import Path
import importlib.util
import numpy as np

OLD=Path(__file__).resolve().parent.parent/'最小闭环验证_2026-10-12/run_experiment.py'
spec=importlib.util.spec_from_file_location('r15_contact_source',OLD)
src=importlib.util.module_from_spec(spec);spec.loader.exec_module(src)
DREF=src.DREF;FD=src.FD
T=1.5; MASS=.2; DAMP=8.;BR=.2;RATE=500.;ALPHA=10.
KMIN=30.;KMAX=300.;KREF=150.;RU=.01/RATE**2

class Carrier:
    def __init__(self,nd=17,nv=17,steps=100,focused=False):
        self.d=np.unique(np.r_[np.linspace(.05*DREF,3*DREF,nd),DREF])
        self.v=np.linspace(-.008,.008,nv)
        self.k=np.linspace(KMIN,KMAX,10)
        if focused:
            # Resolve actual slow contact motion instead of wasting resolution
            # at unattained large velocities; retain the same outer envelope.
            self.d=np.unique(np.r_[np.linspace(.05*DREF,3*DREF,11),
                                  np.linspace(.45*DREF,1.3*DREF,nd),DREF])
            self.v=.008*np.sinh(np.linspace(-5,5,nv))/np.sinh(5)
            self.k=np.array([30.,60.,90.,120.,135.,150.,165.,180.,195.,210.,240.,270.,300.])
        self.shape=(len(self.d),len(self.v),len(self.k))
        self.N=steps;self.h=T/steps;self.ts=np.linspace(0,T,steps+1)
        dd,vv,kk=np.meshgrid(self.d,self.v,self.k,indexing='ij')
        self.x=np.c_[dd.ravel(),vv.ravel(),kk.ravel()]
        self.S=len(self.x);self.a=src.halfwidth(self.x[:,0])
        self.norm=(self.x[:,2]-KREF)/KREF
        self.model_cache={}

    def fn_nodes(self,t,eta):
        # eta=(Ef/E0,b/a0). Pressure geometry follows the original cylinder.
        r=eta[1]*.003+.006-.008*t
        z=np.clip(r/self.a,-1,1)
        hh=np.sqrt(np.maximum(self.a*self.a-r*r,0))
        area=.5*(self.a*self.a*np.arccos(z)-np.clip(r,-self.a,self.a)*hh)
        return .02/(2*.03)*30000*(np.pi*self.a*self.a/2+(eta[0]-1)*area)

    def model(self,eta):
        key=tuple(eta)
        if key in self.model_cache:return self.model_cache[key]
        x=self.x;d,v,k=x.T;rows=[]
        for t in self.ts[:-1]:
            f=self.fn_nodes(t,eta)
            acc=(FD-f-k*(d-DREF)-(DAMP+BR)*v)/MASS
            dn=d+self.h*v;vn=v+self.h*acc
            geometry=self.geometry(dn,vn)
            stage=((d-DREF)/DREF)**2+(v/.008)**2+(f/FD-1)**2+.1*self.norm**2
            rows.append((geometry,stage))
        ft=self.fn_nodes(T,eta)
        phi=((d-DREF)/DREF)**2+(ft/FD-1)**2
        result=(rows,phi);self.model_cache[key]=result
        return result

    def geometry(self,d,v):
        # Clamps the VALUE LOOKUP outside the grid, not the physical dynamics.
        # Grid expansion/refinement and actual rollout checks audit this truncation.
        id=np.clip(np.searchsorted(self.d,d)-1,0,len(self.d)-2)
        iv=np.clip(np.searchsorted(self.v,v)-1,0,len(self.v)-2)
        ad=np.clip((d-self.d[id])/(self.d[id+1]-self.d[id]),0,1)
        av=np.clip((v-self.v[iv])/(self.v[iv+1]-self.v[iv]),0,1)
        rows=np.array([id*len(self.v)+iv,id*len(self.v)+iv+1,
                       (id+1)*len(self.v)+iv,(id+1)*len(self.v)+iv+1],dtype=np.int32)
        weights=np.array([(1-ad)*(1-av),(1-ad)*av,ad*(1-av),ad*av])
        return rows,weights

    def lines(self,next_value,geom):
        rows,weights=geom;data=next_value.reshape(-1,len(self.k))
        return sum(data[rows[j]]*weights[j,:,None] for j in range(4))

    def minimize(self,lines,k):
        # Exact minimization on every piecewise-linear next-value gain interval.
        # Action is continuous gain RATE, not a hand-written gain schedule.
        k=np.atleast_1d(k)
        lo=np.maximum(-RATE,-ALPHA*(k-KMIN))
        hi=np.minimum(RATE,ALPHA*(KMAX-k))
        kl=self.k[:-1];kr=self.k[1:];slope=np.diff(lines,axis=1)/(kr-kl)
        ulo=np.maximum(lo[:,None],(kl-k[:,None])/self.h)
        uhi=np.minimum(hi[:,None],(kr-k[:,None])/self.h)
        u=np.minimum(np.maximum(-slope/(2*RU),ulo),uhi)
        q=self.h*RU*u*u+lines[:,:-1]+slope*(k[:,None]+self.h*u-kl)
        q=np.where(ulo<=uhi+1e-12,q,np.inf)
        best=np.argmin(q,axis=1);ids=np.arange(len(k));action=u[ids,best]
        beta=(k+self.h*action-self.k[best])/(self.k[best+1]-self.k[best])
        return q[ids,best],action,best.astype(np.int16),np.clip(beta,0,1)

    def bellman(self,next_value,index,eta,cache=False):
        (geom,stage)=self.model(eta)[0][index]
        lines=self.lines(next_value,geom)
        val,u,j,beta=self.minimize(lines,self.x[:,2])
        value=self.h*stage+val
        policy=(geom,j,beta) if cache else None
        return value,policy

    def solve(self,eta,cache=False):
        model=self.model(eta)
        V=np.empty((self.N+1,self.S));V[-1]=model[1]
        P=[None]*self.N
        for i in range(self.N-1,-1,-1):V[i],P[i]=self.bellman(V[i+1],i,eta,cache)
        return V,P

    def apply_policy(self,value,P):
        (rows,weights),j,beta=P;data=value.reshape(-1,len(self.k))
        return sum(weights[c]*((1-beta)*data[rows[c],j]+beta*data[rows[c],j+1])
                   for c in range(4))

    def residual(self,V,eta):
        F=np.empty_like(V);F[-1]=V[-1]-self.model(eta)[1]
        for i in range(self.N):F[i]=V[i]-self.bellman(V[i+1],i,eta)[0]
        return F

    def linear_solve(self,F,P):
        # J d = -F for J=(I, -P), terminal derivative I.
        d=np.empty_like(F);d[-1]=-F[-1]
        for i in range(self.N-1,-1,-1):d[i]=-F[i]+self.apply_policy(d[i+1],P[i])
        return d

    def sensitivity(self,V,P,eta,eps=1e-4):
        # Directional derivative on the soft->hard cone. Piecewise interpolation
        # is not claimed C2 on an open domain; convergence is checked separately.
        new=np.array(eta,dtype=float);new[0]+=eps
        F=(self.residual(V,new)-self.residual(V,eta))/eps
        Cv=self.linear_solve(F,P)
        phi0=self.model(eta)[1];dphi=(self.model(new)[1]-phi0)/eps
        return Cv-dphi[None,:],dphi

    def candidate(self,method,eta,V,P,Cw,dphi):
        phi0=self.model((1.,0.))[1];phi=self.model(eta)[1]
        # V=Phi(eta)+W, with terminal W=0. All methods update the same W root.
        base=V+(phi-phi0)[None,:]
        if method=='NEOC':
            result=base+Cw*(eta[0]-1.)
            result[-1]=phi
            return result
        if method=='cached_chord':
            return base+self.linear_solve(self.residual(base,eta),P)
        if method=='one_Newton':
            F=np.empty_like(base);newP=[None]*self.N
            F[-1]=base[-1]-phi
            for i in range(self.N):
                b,newP[i]=self.bellman(base[i+1],i,eta,True);F[i]=base[i]-b
            return base+self.linear_solve(F,newP)
        if method=='full_backward_solve':return self.solve(eta)[0]
        raise ValueError(method)

    def policy(self,V,t,state,eta):
        d,v,k=state[:3];i=min(int(t/self.h),self.N-1)
        fn=src.fast_force(t,d,eta)
        acc=(FD-fn-k*(d-DREF)-(DAMP+BR)*v)/MASS
        geom=self.geometry(np.array([d+self.h*v]),np.array([v+self.h*acc]))
        # Preserve the actual phase rather than holding a whole value time slice.
        fraction=min(max(t/self.h-i,0.),1.)
        line=(1-fraction)*self.lines(V[i+1],geom)+fraction*self.lines(V[min(i+2,self.N)],geom)
        return float(self.minimize(line,np.array([k]))[1][0])

def rollout(carrier,V,model_eta,true_eta=(1.5,0.),event_time=0.,oldV=None,
            delayed_publish=0.,fixed_gain=False,end=T,initial=None):
    h=.001;steps=round(end/h);x=np.array(initial or [DREF,0.,150.,.02],dtype=float)
    data=[];cost=0.
    for j in range(steps+1):
        t=j*h;f=src.fast_force(t,x[0],true_eta)
        old=t<event_time+delayed_publish and oldV is not None
        u=0. if fixed_gain else carrier.policy(oldV if old else V,t,x,(1.,0.) if old else model_eta)
        stage=((x[0]-DREF)/DREF)**2+(x[1]/.008)**2+(f/FD-1)**2+.1*((x[2]-150)/150)**2+RU*u*u
        data.append([t,*x,f,u,stage,cost])
        if j==steps:break
        def rhs(time,y):
            ff=src.fast_force(time,y[0],true_eta)
            vel=y[1];err=y[0]-DREF
            acc=(FD-ff-y[2]*err-(DAMP+BR)*vel)/MASS
            ed=.3*(DAMP+BR)*vel*vel-.5*err*err*u
            return np.array([vel,acc,u,ed])
        a=rhs(t,x);b=rhs(t+h/2,x+h*a/2);c=rhs(t+h/2,x+h*b/2);d=rhs(t+h,x+h*c)
        x=x+h*(a+2*b+2*c+d)/6;cost+=stage*h
        if x[0]<=0:raise FloatingPointError('Physical rollout lost contact')
    f=src.fast_force(end,x[0],true_eta)
    terminal=((x[0]-DREF)/DREF)**2+(f/FD-1)**2
    data=np.array(data);force=data[:,5]
    valid=bool(np.all((data[:,1]>=carrier.d[0])&(data[:,1]<=carrier.d[-1])) and
               np.all(abs(data[:,2])<=carrier.v[-1]) and
               np.all((data[:,3]>=30-1e-8)&(data[:,3]<=300+1e-8)) and
               np.all(data[:,4]>=.001) and np.all(abs(data[:,6])<=500+1e-8) and
               np.all((force>0)&(force<2*FD)))
    return data,dict(cost=float(cost+terminal),force_rmse_N=float(np.sqrt(np.mean((force-FD)**2))),
                     force_peak_N=float(max(force)),final_gain=float(x[2]),minimum_energy=float(min(data[:,4])),
                     sampled_constraints_pass=valid,terminal_cost=float(terminal))
