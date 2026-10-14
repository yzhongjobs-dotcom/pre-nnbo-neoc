"""One R15 spatial soft-to-hard scan; NEOC direct vs warm Newton.

Reduced equations: current manuscript (5),(8)-(11),(28)-(39),(41)-(42).
No online Newton is called by the NEOC path. This is an empirical probe,
not a full PublishOK implementation. Run: python -X utf8 run_experiment.py
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import time
import traceback
import numpy as np
import scipy
from scipy.linalg import qr, lu_factor, lu_solve
from scipy.optimize import least_squares
from scipy.stats import qmc

OUT = Path(__file__).resolve().parent
CFG = dict(R=.03, L=.02, lref=.02, E0=30000., Ef=45000., boundary=0.,
           a0=.003, speed=.008, s0=-.006, T=1.5, dt=.001,
           mass=.2, Br=.2, k0=150., damping0=8., energy0=.02,
           energy_min=.001, chi=.3, kmin=30., kmax=300.,
           dmin=2., dmax=20., krate=500., drate=50., alpha=10.,
           action_weight=.01, gain_weight=.1, seed=20261012,
           node_candidates=512, newton_tolerance=1e-9,
           newton_max_iterations=40, sensitivity_step=1e-5,
           candidate_bellman_residual_allow=.05,
           history_samples=12, history_interval=.005,
           timing_repeats=9,
           thresholds=dict(force_rmse_fraction=.1, cost_relative_extra=.05,
                           force_rmse_extra_fraction=.02, latency_reduction=.3))

def psi(a):
    h = np.sqrt(CFG['lref']**2-a*a)
    return (CFG['lref']**2-CFG['lref']*h+a*a*np.log((CFG['lref']+h)/a))/(2*CFG['R'])

DREF = float(psi(CFG['a0']))
FD = np.pi*CFG['L']*CFG['E0']*CFG['a0']**2/(4*CFG['R'])
SCALE = np.array([DREF, CFG['speed'], CFG['k0'], CFG['damping0'], CFG['energy0']])
CENTER = np.array([DREF, 0., CFG['k0'], CFG['damping0'], CFG['energy0']])
RATE = np.array([CFG['krate'], CFG['drate']])

def halfwidth(delta):
    """Vector bracketed Newton for the monotone cylindrical indentation law."""
    delta = np.asarray(delta)
    if np.any(delta <= 0) or np.any(delta >= psi(CFG['lref']*.999)):
        raise ValueError('Indentation outside the calibrated contact branch')
    lo = np.full_like(delta, 1e-12, dtype=float)
    hi = np.full_like(delta, CFG['lref']*.999, dtype=float)
    a = np.minimum(CFG['a0']*np.sqrt(delta/DREF), hi*.95)
    for _ in range(15):
        h = np.sqrt(CFG['lref']**2-a*a)
        log = np.log((CFG['lref']+h)/a)
        val = (CFG['lref']**2-CFG['lref']*h+a*a*log)/(2*CFG['R'])
        err = val-delta
        if np.max(np.abs(err)) < 2e-15:
            break
        lo = np.where(err < 0, a, lo)
        hi = np.where(err >= 0, a, hi)
        nxt = a-err/(a/CFG['R']*log)
        update = np.where((nxt > lo)&(nxt < hi), nxt, .5*(lo+hi))
        a = np.where(np.abs(err)<2e-15, a, update)
    if np.max(np.abs(psi(a)-delta)) > 1e-12:
        raise ValueError('Geometry inversion did not converge')
    return a

def contact(t, delta, eta, moments=False):
    """eta=(Ef/E0,b/a0); raw My includes known tangential traction."""
    a = halfwidth(delta)
    s = CFG['s0']+CFG['speed']*np.asarray(t)
    r = eta[1]*CFG['a0']-s
    z = np.clip(r/a, -1., 1.)
    theta = np.arccos(z)
    h = np.sqrt(np.maximum(a*a-r*r, 0.))
    area = .5*(a*a*theta-np.clip(r,-a,a)*h)
    efront = eta[0]*CFG['E0']
    const = CFG['L']/(2*CFG['R'])
    fn = const*(CFG['E0']*np.pi*a*a/2+(efront-CFG['E0'])*area)
    slope_a = const*a*(CFG['E0']*np.pi+(efront-CFG['E0'])*theta)
    ah = np.sqrt(CFG['lref']**2-a*a)
    da_dd = CFG['R']/(a*np.log((CFG['lref']+ah)/a))
    slope_d = slope_a*da_dd
    if not moments:
        return fn, slope_d
    m1 = const*(efront-CFG['E0'])*h**3/3
    i2 = np.pi*a**4/8
    i2front = a**4/8*(theta+z*(1-2*z*z)*np.sqrt(np.maximum(1-z*z,0.)))
    tangential_moment = -.2/(2*CFG['R'])*const*(CFG['E0']*i2+(efront-CFG['E0'])*i2front)
    tangential_moment -= CFG['speed']*1.2/(2*CFG['R'])*a*a/4
    return fn, -m1+tangential_moment

# Complete degree <=2 monomials in (t/T, normalized five-state coordinates).
EXP = [np.zeros(6,dtype=int)]
EXP += [np.eye(6,dtype=int)[i] for i in range(6)]
EXP += [np.eye(6,dtype=int)[i]+np.eye(6,dtype=int)[j]
        for i in range(6) for j in range(i,6)]
EXP = np.array(EXP)
assert len(EXP)==28

def features(t, y):
    q = np.column_stack([np.asarray(t)/CFG['T'], y])
    val = np.prod(q[:,None,:]**EXP[None,:,:],axis=2)
    dv = np.zeros((len(q),28,6))
    for j in range(6):
        mask=EXP[:,j]>0
        ex=EXP[mask].copy(); ex[:,j]-=1
        dv[:,mask,j]=EXP[mask,j]*np.prod(q[:,None,:]**ex[None,:,:],axis=2)
    fac=1-q[:,0]
    return fac[:,None]*val, (-val+fac[:,None]*dv[:,:,0])/CFG['T'], fac[:,None,None]*dv[:,:,1:]

def state_terms(t,y,eta):
    x=CENTER+y*SCALE
    d,v,k,damp,energy=x.T
    force,_=contact(t,d,eta)
    accel=(FD-force-k*(d-DREF)-(damp+CFG['Br'])*v)/CFG['mass']
    pd=(damp+CFG['Br'])*v*v
    f=np.column_stack([v/DREF,accel/CFG['speed'],np.zeros(len(y)),np.zeros(len(y)),CFG['chi']*pd/CFG['energy0']])
    G=np.zeros((len(y),5,2)); G[:,2,0]=RATE[0]/SCALE[2]; G[:,3,1]=RATE[1]/SCALE[3]
    G[:,4,0]=-.5*(d-DREF)**2*RATE[0]/SCALE[4]
    terminal_force,terminal_slope=contact(np.full(len(y),CFG['T']),d,eta)
    grad_phi=np.zeros_like(y)
    grad_phi[:,0]=2*y[:,0]+2*(terminal_force/FD-1)*terminal_slope*DREF/FD
    ell=y[:,0]**2+y[:,1]**2+(force/FD-1)**2+CFG['gain_weight']*(y[:,2]**2+y[:,3]**2)
    A=np.zeros((len(y),9,2)); b=np.zeros((len(y),9))
    for row,col,sgn in [(0,0,1),(1,0,-1),(2,1,1),(3,1,-1),(4,0,-1),(5,0,1),(6,1,-1),(7,1,1)]:
        A[:,row,col]=sgn
    b[:,:4]=1
    b[:,4]=CFG['alpha']*(k-CFG['kmin'])/RATE[0]
    b[:,5]=CFG['alpha']*(CFG['kmax']-k)/RATE[0]
    b[:,6]=CFG['alpha']*(damp-CFG['dmin'])/RATE[1]
    b[:,7]=CFG['alpha']*(CFG['dmax']-damp)/RATE[1]
    A[:,8,0]=.5*(d-DREF)**2*RATE[0]/SCALE[4]
    b[:,8]=(CFG['chi']*pd+CFG['alpha']*(energy-CFG['energy_min']))/SCALE[4]
    return f,G,grad_phi,ell,A,b

def action(c,A,b):
    """Exact diagonal 2-D QP: all nine rows are coordinate bounds."""
    n=len(c); lower=np.full((n,2),-np.inf); upper=np.full((n,2),np.inf)
    low_row=np.full((n,2),-1,dtype=int); up_row=low_row.copy()
    for row in range(9):
        for j in range(2):
            aa=A[:,row,j]; nz=np.abs(aa)>1e-18
            bound=np.divide(b[:,row],aa,out=np.zeros(n),where=nz)
            up=nz&(aa>0)&(bound<upper[:,j]); low=nz&(aa<0)&(bound>lower[:,j])
            upper[up,j]=bound[up]; up_row[up,j]=row
            lower[low,j]=bound[low]; low_row[low,j]=row
        if np.any((np.linalg.norm(A[:,row],axis=1)<1e-18)&(b[:,row]<-1e-9)):
            raise ValueError('Zero action constraint row infeasible')
    if np.any(lower>upper+1e-8): raise ValueError('Empty action set')
    free=-c/(2*CFG['action_weight'])
    u=np.minimum(np.maximum(free,lower),upper)
    lam=np.zeros((n,9)); active=np.full((n,2),-1,dtype=int)
    for j in range(2):
        sel=np.where(free[:,j]<lower[:,j],low_row[:,j],np.where(free[:,j]>upper[:,j],up_row[:,j],-1))
        active[:,j]=sel
        ids=np.flatnonzero(sel>=0); rows=sel[ids]
        lam[ids,rows]=-(2*CFG['action_weight']*u[ids,j]+c[ids,j])/A[ids,rows,j]
    return u,lam,active

class RootProblem:
    def __init__(self):
        sample=qmc.Sobol(6,scramble=True,seed=CFG['seed']).random_base2(9)
        lo=np.array([.02,-.35,-.045,-.45,-.30,-.4])
        hi=np.array([.94, .35, .045, .65, .50, .4])
        p=lo+(hi-lo)*sample
        a=halfwidth(DREF*(1+p[:,1])); r=-CFG['s0']-CFG['speed']*p[:,0]*CFG['T']
        p=p[(np.abs(r-a)>.0002)&(np.abs(r+a)>.0002)]
        t=p[:,0]*CFG['T']; y=p[:,1:]
        ph,_,_=features(t,y)
        norms=np.linalg.norm(ph,axis=0)
        _,_,order=qr((ph/norms).T,pivoting=True,mode='economic')
        self.t=t[order[:28]]; self.y=y[order[:28]]
        self.ps,self.pt,self.px=features(self.t,self.y)
        self.weight_scale=1/np.maximum(np.linalg.norm(np.column_stack([self.pt,self.px.reshape(28,-1)]),axis=0)[:28],1e-6)
        # Matrix column preconditioning is fixed, and does not change the root.
        f,G,_,_,_,_=state_terms(self.t,self.y,np.array([1.,0.]))
        L=self.pt+np.einsum('nwj,nj->nw',self.px,f)
        self.weight_scale=1/np.maximum(np.linalg.norm(L,axis=0),1e-8)
        self.row_scale=1/np.maximum(np.linalg.norm(L*self.weight_scale,axis=1),1e-8)
        self.cached={}
    def terms(self,eta):
        key=tuple(eta)
        if key not in self.cached:self.cached[key]=state_terms(self.t,self.y,eta)
        return self.cached[key]
    def evaluate(self,w,eta,jac=False):
        f,G,gphi,ell,A,b=self.terms(eta)
        grad=gphi+np.einsum('nwj,w->nj',self.px,w)
        c=np.einsum('nji,nj->ni',G,grad)
        u,lam,active=action(c,A,b)
        flow=f+np.einsum('nji,ni->nj',G,u)
        bell=self.pt@w+ell+CFG['action_weight']*np.sum(u*u,axis=1)+np.sum(grad*flow,axis=1)
        if jac:
            L=self.pt+np.einsum('nwj,nj->nw',self.px,flow)
            return bell,L,u,lam,active
        return bell,u,lam,active
    def solve(self,eta,w0):
        w=w0.copy(); trace=[]
        for it in range(CFG['newton_max_iterations']+1):
            r,L,u,lam,act=self.evaluate(w,eta,True)
            err=float(np.max(np.abs(r))); trace.append(err)
            if err<CFG['newton_tolerance']:
                return w,dict(success=True,iterations=it,residual=err,trace=trace,
                              jacobian_condition=float(np.linalg.cond(L)),active=act.tolist())
            if it==CFG['newton_max_iterations']:break
            rs=r*self.row_scale
            step=np.linalg.solve(L*self.row_scale[:,None]*self.weight_scale[None,:],-rs)*self.weight_scale
            merit=np.linalg.norm(rs)
            accepted=False
            for k in range(18):
                trial=w+step*2.**(-k)
                rt=self.evaluate(trial,eta)[0]*self.row_scale
                if np.linalg.norm(rt)<merit*(1-1e-4*2.**(-k)):
                    w=trial; accepted=True; break
            if not accepted:break
        return w,dict(success=False,iterations=len(trace)-1,residual=trace[-1],trace=trace)
    def lift(self,w,eta):
        _,u,lam,_=self.evaluate(w,eta)
        return np.concatenate([w,u.ravel(),lam.ravel()])

def sensitivity(problem,w,eta):
    r,L,*_=problem.evaluate(w,eta,True)
    fac=lu_factor(L); step=CFG['sensitivity_step']
    partial=np.column_stack([(problem.evaluate(w,eta+np.eye(2)[j]*step)[0]-problem.evaluate(w,eta-np.eye(2)[j]*step)[0])/(2*step) for j in range(2)])
    cw=lu_solve(fac,-partial)
    cz=np.column_stack([(problem.lift(w+cw[:,j]*step,eta+np.eye(2)[j]*step)-problem.lift(w-cw[:,j]*step,eta-np.eye(2)[j]*step))/(2*step) for j in range(2)])
    return cw,cz,dict(linear_residual=float(np.max(abs(L@cw+partial))),
                      anchor_residual=float(np.max(abs(r))),jacobian_condition=float(np.linalg.cond(L)))

def polynomial(w):
    linear=w[1:7].copy(); H=np.zeros((6,6))
    pos=7
    for i in range(6):
        for j in range(i,6):
            if i==j:H[i,i]=2*w[pos]
            else:H[i,j]=H[j,i]=w[pos]
            pos+=1
    return float(w[0]),linear,H

def fast_width(d):
    if not 0<d<float(psi(CFG['lref']*.999)):raise ValueError('Lost contact / invalid indentation')
    a=CFG['a0']*math.sqrt(d/DREF)
    for _ in range(7):
        h=math.sqrt(CFG['lref']**2-a*a); lg=math.log((CFG['lref']+h)/a)
        err=(CFG['lref']**2-CFG['lref']*h+a*a*lg)/(2*CFG['R'])-d
        if abs(err)<2e-15:return a
        a-=err/(a/CFG['R']*lg)
    return a

def fast_force(t,d,eta):
    a=fast_width(d); r=eta[1]*CFG['a0']-CFG['s0']-CFG['speed']*t
    z=max(-1.,min(1.,r/a)); h=math.sqrt(max(a*a-r*r,0.))
    area=.5*(a*a*math.acos(z)-max(-a,min(a,r))*h)
    return CFG['L']/(2*CFG['R'])*(CFG['E0']*math.pi*a*a/2+(eta[0]-1)*CFG['E0']*area)

def policy(t,x,coeff):
    _,lin,H=coeff
    q=np.r_[t/CFG['T'],(x-CENTER)/SCALE]
    fac=1-t/CFG['T']; e=x[0]-DREF
    gradient=fac*(lin[1:]+H[1:]@q)
    c=np.array([RATE[0]/SCALE[2]*gradient[2]-.5*e*e*RATE[0]/SCALE[4]*gradient[4],RATE[1]/SCALE[3]*gradient[3]])
    free=-c/(2*CFG['action_weight'])
    low=np.maximum(-1.,-CFG['alpha']*(x[2:4]-[CFG['kmin'],CFG['dmin']])/RATE)
    high=np.minimum(1.,CFG['alpha']*([CFG['kmax'],CFG['dmax']]-x[2:4])/RATE)
    pd=(x[3]+CFG['Br'])*x[1]**2
    coef=.5*e*e*RATE[0]
    if coef>1e-18:high[0]=min(high[0],(CFG['chi']*pd+CFG['alpha']*(x[4]-CFG['energy_min']))/coef)
    if np.any(low>high+1e-8):raise ValueError('Action set became empty')
    return np.minimum(np.maximum(free,low),high)

def rhs(t,x,coeff,eta):
    u=policy(t,x,coeff); d,v,k,damp,E=x
    force=fast_force(t,d,eta)
    acc=(FD-force-k*(d-DREF)-(damp+CFG['Br'])*v)/CFG['mass']
    pd=(damp+CFG['Br'])*v*v
    return np.array([v,acc,RATE[0]*u[0],RATE[1]*u[1],CFG['chi']*pd-.5*(d-DREF)**2*RATE[0]*u[0]])

def rollout(start,x0,w,eta,end=None,held_until=None,old_w=None):
    end=CFG['T'] if end is None else end
    coeff=polynomial(w); oldcoeff=polynomial(old_w) if old_w is not None else coeff
    x=x0.copy(); t=start; records=[]
    while True:
        cc=oldcoeff if held_until is not None and t<held_until else coeff
        u=policy(t,x,cc); force=fast_force(t,x[0],eta); y=(x-CENTER)/SCALE
        running=y[0]**2+y[1]**2+(force/FD-1)**2+CFG['gain_weight']*(y[2]**2+y[3]**2)+CFG['action_weight']*sum(u*u)
        records.append([t,*x,force,*u,running])
        if t>=end-1e-12:break
        h=min(CFG['dt'],end-t)
        a=rhs(t,x,cc,eta); b=rhs(t+h/2,x+h*a/2,cc,eta)
        c=rhs(t+h/2,x+h*b/2,cc,eta); d=rhs(t+h,x+h*c,cc,eta)
        x=x+h*(a+2*b+2*c+d)/6; t+=h
    return np.array(records)

def metrics(data,eta):
    t=data[:,0]; x=data[:,1:6]; fn=data[:,6]; u=data[:,7:9]
    cost=float(np.trapz(data[:,-1],t))
    terminal=(x[-1,0]-DREF)**2/DREF**2+(fast_force(CFG['T'],x[-1,0],eta)/FD-1)**2
    margins=dict(indentation_lower=float(x[:,0].min()),indentation_upper=float(psi(.006)-x[:,0].max()),
                 force_lower=float(fn.min()),force_upper=float(2*FD-fn.max()),
                 k_lower=float(x[:,2].min()-CFG['kmin']),k_upper=float(CFG['kmax']-x[:,2].max()),
                 damping_lower=float(x[:,3].min()-CFG['dmin']),damping_upper=float(CFG['dmax']-x[:,3].max()),
                 energy_lower=float(x[:,4].min()-CFG['energy_min']),normalized_rate_margin=float(1-np.abs(u).max()))
    return dict(force_rmse=float(np.sqrt(np.trapz((fn-FD)**2,t)/(t[-1]-t[0]))),
                force_rmse_fraction=float(np.sqrt(np.trapz((fn/FD-1)**2,t)/(t[-1]-t[0]))),
                force_peak=float(fn.max()),cost=cost+float(terminal),
                sampled_constraints_pass=bool(min(margins.values())>=-1e-8),margins=margins,
                contact_loss_samples=int(sum(fn<=0)),k_final=float(x[-1,2]),damping_final=float(x[-1,3]))

def check(problem,w,eta,t,x):
    bell,u,lam,act=problem.evaluate(w,eta)
    _,_,_,_,A,b=problem.terms(eta)
    slack=b-np.einsum('nki,ni->nk',A,u)
    nodal=bool(np.min(slack)>=-1e-8 and np.min(lam)>=-1e-8)
    data=rollout(t,x,w,eta)
    met=metrics(data,eta)
    first=policy(t,x,polynomial(w))
    return dict(accepted=bool(nodal and met['sampled_constraints_pass'] and max(abs(bell))<=CFG['candidate_bellman_residual_allow']),
                max_nodal_bellman_residual=float(max(abs(bell))),nodal_kkt_feasible=nodal,
                rollout_metrics=met,first_normalized_action=first.tolist(),
                label='empirical nodal/rollout screening, not PublishOK')

def save_csv(name,data,header):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f); writer.writerow(header); writer.writerows(data)

def main():
    start=time.perf_counter()
    source=Path(__file__).read_bytes()
    result=dict(config=CFG,dref=DREF,Fd=float(FD),python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
                source_sha256=hashlib.sha256(source).hexdigest(),neoc_online_newton_calls=0,
                neoc_path_executed=False)
    (OUT/'config.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    try:
        problem=RootProblem(); eta0=np.array([1.,0.]); truth=np.array([1.5,0.])
        save_csv('fixed_nodes.csv',np.column_stack([problem.t,problem.y]),['t','normalized_d','normalized_v','normalized_k','normalized_damping','normalized_energy'])
        print('Solving shared 28-weight / 336-variable HJB-KKT anchor...',flush=True)
        ta=time.perf_counter_ns(); w0,anchor=problem.solve(eta0,np.zeros(28)); anchor_ns=time.perf_counter_ns()-ta
        result['anchor']=anchor; result['anchor_solve_ms']=anchor_ns/1e6
        np.savez(OUT/'anchor.npz',w=w0,eta=eta0)
        print(json.dumps(anchor),flush=True)
        if not anchor['success']:
            result['status']='STOP_SHARED_ANCHOR_NOT_SOLVED'; return result
        print('Preparing NEOC sensitivity (no perturbed root solve)...',flush=True)
        ts=time.perf_counter_ns(); cw,cz,sens=sensitivity(problem,w0,eta0); prep_ns=time.perf_counter_ns()-ts
        result['sensitivity']=sens; result['sensitivity_prepare_ms']=prep_ns/1e6
        np.savez(OUT/'sensitivity.npz',cw=cw,cz=cz)
        # Derivative check uses fixed W and the same analytic envelope Jacobian.
        direction=np.random.default_rng(CFG['seed']).normal(size=28)*problem.weight_scale
        eps=1e-5
        fd=(problem.evaluate(w0+eps*direction,eta0)[0]-problem.evaluate(w0-eps*direction,eta0)[0])/(2*eps)
        _,L,*_=problem.evaluate(w0,eta0,True)
        result['analytic_jacobian_direction_relative_error']=float(np.linalg.norm(fd-L@direction)/max(np.linalg.norm(L@direction),1e-12))
        assert result['analytic_jacobian_direction_relative_error']<1e-6
        event=(.006-.003)/CFG['speed']+(CFG['history_samples']-1)*CFG['history_interval']
        prefix=rollout(0.,CENTER.copy(),w0,truth,end=event)
        times=np.linspace(event-(CFG['history_samples']-1)*CFG['history_interval'],event,CFG['history_samples'])
        depths=np.interp(times,prefix[:,0],prefix[:,1])
        force,moment=contact(times,depths,truth,True)
        def residual(eta):
            f,m=contact(times,depths,eta,True)
            return np.r_[(f-force)/FD,(m-moment)/(FD*CFG['a0'])]
        fit=least_squares(residual,[1.2,.2],bounds=([.5,-1.],[3.,1.]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=100)
        eta=fit.x; xevent=prefix[-1,1:6].copy()
        result['estimate']=dict(eta=eta.tolist(),Ef=float(eta[0]*CFG['E0']),boundary=float(eta[1]*CFG['a0']),residual_norm=float(np.linalg.norm(fit.fun)),success=bool(fit.success),nfev=int(fit.nfev),event_time=event,event_state=xevent.tolist())
        print('Shared contact estimate: '+json.dumps(result['estimate']),flush=True)
        save_csv('shared_prefix.csv',prefix,['t','d','v','k','damping','energy','Fn','u_k_normalized','u_d_normalized','running_cost'])
        save_csv('contact_history.csv',np.column_stack([times,depths,force,moment]),['t','d','Fn','raw_My'])
        print('One direct NEOC update and warm Newton reference...',flush=True)
        result['neoc_path_executed']=True
        neoc_w=w0+cw@(eta-eta0)
        neoc_z=problem.lift(w0,eta0)+cz@(eta-eta0)
        newton_w,newton=problem.solve(eta,w0)
        result['newton']=newton
        if not newton['success']:
            result['status']='STOP_WARM_REFERENCE_NOT_SOLVED';return result
        result['node_active_set_changes']=int(np.sum(np.asarray(anchor['active'])!=np.asarray(newton['active'])))
        result['predicted_full_root']=dict(norm_error=float(np.linalg.norm(neoc_z-problem.lift(newton_w,eta))),norm_prediction=float(np.linalg.norm(neoc_z)))
        np.savez(OUT/'updated_weights.npz',neoc=neoc_w,newton=newton_w,eta=eta)
        # Repeated identical query, interleaved order: full screen and first action.
        timings={'NEOC':[],'warm_Newton':[]}; candidates={}
        for repeat in range(CFG['timing_repeats']):
            for method in (['NEOC','warm_Newton'] if repeat%2==0 else ['warm_Newton','NEOC']):
                t0=time.perf_counter_ns()
                if method=='NEOC':w=w0+cw@(eta-eta0)
                else:
                    w,meta=problem.solve(eta,w0)
                    assert meta['success']
                cand=check(problem,w,eta,event,xevent)
                timings[method].append((time.perf_counter_ns()-t0)/1e6)
                candidates[method]=cand
            print('Timing pair %d/%d'%(repeat+1,CFG['timing_repeats']),flush=True)
        result['candidates']=candidates
        result['timing']={m:dict(samples_ms=v,median_ms=float(np.median(v)),p95_ms=float(np.percentile(v,95))) for m,v in timings.items()}
        result['total_update_prepare_cost_ms']={'NEOC':result['timing']['NEOC']['median_ms']+result['sensitivity_prepare_ms'], 'warm_Newton':result['timing']['warm_Newton']['median_ms']}
        result['latency_reduction']=1-result['timing']['NEOC']['median_ms']/result['timing']['warm_Newton']['median_ms']
        met={}
        for method,w in [('NEOC',neoc_w),('warm_Newton',newton_w)]:
            if not candidates[method]['accepted']:
                met[method]=dict(status='REJECTED_NO_EXECUTION');continue
            delay=result['timing'][method]['median_ms']/1000
            physical=rollout(event,xevent,w,truth,held_until=event+delay,old_w=w0)
            full=np.vstack([prefix[:-1],physical]); met[method]=metrics(full,truth)
            save_csv(method+'_trajectory.csv',full,['t','d','v','k','damping','energy','Fn','u_k_normalized','u_d_normalized','running_cost'])
        result['closed_loop']=met
        result['status']='COMPLETED' if all(candidates[m]['accepted'] for m in candidates) else 'CANDIDATE_REJECTED'
        if result['status']=='COMPLETED':
            a,b=met['NEOC'],met['warm_Newton']; th=CFG['thresholds']
            result['judgment']=dict(task_pass=bool(a['sampled_constraints_pass'] and a['force_rmse_fraction']<=th['force_rmse_fraction']),
              reference_task_pass=bool(b['sampled_constraints_pass'] and b['force_rmse_fraction']<=th['force_rmse_fraction']),
              cost_relative_extra=a['cost']/b['cost']-1,
              force_rmse_extra_fraction=(a['force_rmse']-b['force_rmse'])/FD,
              latency_pass=bool(result['latency_reduction']>=th['latency_reduction']),
              prepare_cost_paid_back=bool(result['total_update_prepare_cost_ms']['NEOC']<result['total_update_prepare_cost_ms']['warm_Newton']))
        return result
    except Exception:
        result['status']='FAILED';result['error']=traceback.format_exc();return result
    finally:
        result['duration_seconds']=time.perf_counter()-start
        (OUT/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')

if __name__=='__main__':
    result=main()
    print('FINAL_STATUS: '+result['status'],flush=True)
    if 'error' in result: print(result['error'],flush=True)
