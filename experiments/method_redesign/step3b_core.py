"""Shared causal predictors and unmodified-FINAL integration for step 3b.
Only this directory is new. Units internally: m, s; network uses mm.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_name] = '1'
import sys, math, time, types
from pathlib import Path
import numpy as np
from numba import njit
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
import common, variants, proposed, cohort_eval as CE, scenes_v2 as SV
from imm import IMMTracker
BASE_PLANT = common.P
FAMS = CE.FAMS
DT, L, TAU, W = .008, .012, .019, 63
RES = HERE/'results'
RES.mkdir(exist_ok=True)

# Training/validation contain codes 0..4. Codes 5..7 are strictly OOD.
# Extended p is consumed only by a separate in-memory plant module.
SURFACE = '''@njit(cache=False)
def surface(x,t,p):
    code = int(p[28]) if len(p)>28 else 0
    omega=2*math.pi*p[3]
    ph=omega*t+p[4]
    if code==7:
        # phase-continuous frequency jump; never present in training
        tc=3.25; w2=1.8*omega
        ph=omega*min(t,tc)+w2*max(t-tc,0.)+p[4]
        omega=omega if t<tc else w2
    h=p[2]*math.sin(ph); ht=p[2]*omega*math.cos(ph)
    if code in (1,2,3,4,5,6):
        h+=.3*p[2]*math.sin(1.71*ph+.4)
        ht+=.3*p[2]*1.71*omega*math.cos(1.71*ph+.4)
    w,hr,hx=ridges(x,p); h+=hr
    for k in range(int(p[15])):
        t0=p[16+4*k]; amp=p[17+4*k]; tr=max(p[18+4*k],1e-3); hold=p[19+4*k]
        for direction in range(2):
            u=(t-t0-(tr+hold)*direction)/tr
            z=min(max(u,0.),1.); dz=0.
            if code==1:
                v=.5-.5*math.cos(math.pi*z)
                if 0<u<1: dz=.5*math.pi*math.sin(math.pi*z)
            elif code==2:
                v=z
                if 0<u<1: dz=1.
            elif code==3:
                v=10*z**3-15*z**4+6*z**5
                if 0<u<1: dz=30*z*z*(1-z)**2
            elif code==5:
                # trapezoidal velocity: acceleration, cruise, deceleration
                if z<.25: v=z*z/.375; dz=2*z/.375
                elif z<.75: v=(z-.125)/.75; dz=1/.75
                else: v=1-(1-z)**2/.375; dz=2*(1-z)/.375
                if not 0<u<1: dz=0.
            elif code==6:
                # two-stage movement separated by a pause
                a=min(max(z/.4,0.),1.); b=min(max((z-.6)/.4,0.),1.)
                v=.5*(a*a*(3-2*a)+b*b*(3-2*b))
                if 0<z<.4: dz=3*a*(1-a)/.4
                if .6<z<1: dz=3*b*(1-b)/.4
            else:
                v=z*z*(3-2*z)
                if 0<u<1: dz=6*z*(1-z)
            sign=1-2*direction
            h+=sign*amp*v; ht+=sign*amp*dz/tr
    return h,ht,hx
'''
_PROFILE = None

def profile_plant():
    global _PROFILE
    if _PROFILE is None:
        import re
        src=Path(BASE_PLANT.__file__).read_text(encoding='utf8')
        src=src.replace('@njit(cache=True)', '@njit(cache=False)')
        a=src.index('@njit(cache=False)\ndef surface(')
        b=src.index('@njit(cache=False)\ndef psi(',a)
        src=src[:a]+SURFACE+'\n'+src[b:]
        mod=types.ModuleType('step3b_profile_plant'); mod.__file__=BASE_PLANT.__file__
        sys.modules[mod.__name__]=mod
        exec(compile(src, str(HERE/'step3b_virtual_plant.py'), 'exec'),mod.__dict__)
        _PROFILE=mod
    return _PROFILE

def use_plant(code):
    mod=BASE_PLANT if code==0 else profile_plant()
    common.P=mod; proposed.P=mod
    return mod

def make_scenes(base,n,profiles=(0,)):
    assert base not in range(130000,140000), 'final is forbidden'
    out=[]
    for fi,fam in enumerate(FAMS):
        rng=np.random.default_rng(base+1000*fi)
        for i in range(n):
            p,T=SV.make(fam,rng)
            code=profiles[i%len(profiles)]
            if code:
                p=np.r_[p,float(code)]
                p[2]*=rng.uniform(.75,1.25); p[3]=rng.uniform(.15,.45)
            out.append((fam,base+1000*fi+i,p,T,code))
    return out

def phases(t,p):
    phase=np.full(len(t),'breathing',dtype='<U16')
    for i in range(int(p[15])):
        t0,amp,ramp,hold=p[16+4*i:20+4*i]
        phase[(t>=t0-.05)&(t<t0)]='pre50ms'
        phase[(t>=t0)&(t<t0+ramp)]='ramp'
        phase[(t>=t0+ramp)&(t<t0+ramp+hold)]='hold'
        phase[(t>=t0+ramp+hold)&(t<t0+2*ramp+hold)]='retract'
    return phase

class Recorder(variants.MotionAwareIMM):
    def reset(self,pl):
        super().reset(pl); self.records=[]; self.segment=-1
    def drift_pred(self,o,n):
        prev=self._motion_n
        was_empty=self._fin_prev is None
        q=super().drift_pred(o,n)
        if q is not None and prev!=n:
            reset=(was_empty or prev is None or prev!=n-1)
            if reset: self.segment+=1
            cp=np.array([1.,L+TAU,L*L/2+TAU*L])
            fca=L**5/20+TAU*L**4/4+TAU**2*L**3/3
            qmix=self.tracker.mu@np.array([self.motion_params['q_smooth'],self.motion_params['q_maneuver']])
            self.records.append((o.log[n-1,0],o.log[n-1,1],self.Z,self.tracker.R,
                                 float(q[1]),float(cp@self.tracker.P@cp+qmix*fca),reset,self.segment))
        return q

def data_episode(scene):
    fam,sid,p,T,code=scene
    plant=use_plant(code)
    c=Recorder(**CE.FINAL)
    log,_=common.run(c,p,T,sid)
    arr=np.asarray(c.records,float)
    if not len(arr): raise RuntimeError('empty measurement episode')
    t,x,z,R,p0,v0,reset,segment=arr.T
    target=np.empty(len(z)); surface=np.empty(len(z)); deriv=np.empty(len(z))
    for j in range(len(z)):
        if reset[j]: anchor=plant.surface(x[j],t[j],p)[0]
        h,ht,_=plant.surface(x[j],t[j]+L,p)
        # Plant uses penetration y+surface: eta=surface+constant=-s_eq+constant.
        # This is the positive surface sign returned by proposed.Oracle.
        target[j]=h+TAU*ht-anchor
        surface[j]=h; deriv[j]=ht
    return dict(family=fam,sid=sid,code=code,p=p,t=t,z=z,R=R,target=target,
                p0=p0,v0=v0,reset=reset.astype(bool),phase=phases(t,p),
                full_length=log is not None,segment=segment.astype(int))

@njit(cache=True)
def imm_step(xs,ps,mu,z,R,A,Q,trans):
    m,d=xs.shape
    joint=mu.reshape((m,1))*trans
    prior=joint.sum(axis=0)
    mix=joint/prior.reshape((1,m))
    mx=mix.T@xs
    mp=np.zeros_like(ps)
    for j in range(m):
        for i in range(m):
            delta=xs[i]-mx[j]
            mp[j]+=mix[i,j]*(ps[i]+np.outer(delta,delta))
    ll=np.zeros(m)
    for j in range(m):
        xp=A[j]@mx[j]; pp=A[j]@mp[j]@A[j].T+Q[j]
        inn=z-xp[0]; S=max(pp[0,0]+R,1e-30); K=pp[:,0]/S
        xs[j]=xp+K*inn
        J=np.eye(d); J[:,0]-=K
        ps[j]=J@pp@J.T+R*np.outer(K,K)
        ll[j]=math.log(max(prior[j],1e-30))-.5*(math.log(S)+inn*inn/S)
    weights=np.exp(ll-ll.max()); weights/=weights.sum()
    mu=1e-12+(1-m*1e-12)*weights
    return xs,ps,mu

class MultiIMM:
    def __init__(self,qs=1e-5,qm=.1,qj=100.,stay=.98,third=True):
        m=3 if third else 2; d=4
        Aca=np.eye(d); Aca[0,1]=DT; Aca[0,2]=DT**2/2; Aca[1,2]=DT
        Aca[3,3]=0.
        Acj=Aca.copy(); Acj[3,3]=1.; Acj[0,3]=DT**3/6; Acj[1,3]=DT**2/2; Acj[2,3]=DT
        static=np.zeros((d,d)); static[0,0]=1.
        self.A=np.array([static,Aca,Acj] if third else [Aca,Aca])
        def white(order):
            Q=np.zeros((d,d))
            for i in range(order):
                for j in range(order):
                    power=2*order-1-i-j
                    Q[i,j]=DT**power/(power*math.factorial(order-1-i)*math.factorial(order-1-j))
            return Q
        Qs=np.zeros((d,d)); Qs[0,0]=qs*DT**5/20
        self.Q=np.array([Qs,qm*white(3),qj*white(4)] if third else [qs*white(3),qm*white(3)])
        self.H=np.array([[1,0,0,0],[1,L+TAU,L*L/2+TAU*L,0],
                         [1,L+TAU,L*L/2+TAU*L,L**3/6+TAU*L*L/2]] if third else
                        [[1,L+TAU,L*L/2+TAU*L,0]]*2,dtype=float)
        self.trans=np.full((m,m),(1-stay)/(m-1)); np.fill_diagonal(self.trans,stay)
        fca=L**5/20+TAU*L**4/4+TAU**2*L**3/3
        fcj=L**7/252+TAU*L**6/36+TAU**2*L**5/20
        self.future_var=np.array([Qs[0,0]*L/DT,qm*fca,qj*fcj] if third else [qs*fca,qm*fca])
        self.m=m
    def reset(self,z,R):
        self.xs=np.zeros((self.m,4)); self.xs[:,0]=z
        self.ps=np.tile(np.diag([R,.03**2,.3**2,3.**2]),(self.m,1,1))
        if self.m==2: self.ps[:,3,3]=0.
        self.mu=np.ones(self.m)/self.m
    def update(self,z,R):
        self.xs,self.ps,self.mu=imm_step(self.xs,self.ps,self.mu,z,R,self.A,self.Q,self.trans)
    def output(self):
        means=np.sum(self.xs*self.H,axis=1)
        mean=float(self.mu@means)
        var=float(sum(self.mu[i]*(self.H[i]@self.ps[i]@self.H[i]+(means[i]-mean)**2+self.future_var[i]) for i in range(self.m)))
        return mean,max(var,1e-18)
    def seq(self,e):
        out=[]; var=[]
        for z,R,reset in zip(e['z'],e['R'],e['reset']):
            if reset: self.reset(z,R)
            else: self.update(z,R)
            y,v=self.output(); out.append(y);var.append(v)
        return np.array(out),np.array(var)

@njit(cache=True)
def imm_sequence(z,R,reset,A,Q,H,trans,future_var):
    m=A.shape[0]; d=A.shape[1]
    xs=np.zeros((m,d)); ps=np.zeros((m,d,d)); mu=np.ones(m)/m
    out=np.zeros(len(z)); var=np.zeros(len(z))
    for n in range(len(z)):
        if reset[n]:
            xs[:]=0.; xs[:,0]=z[n]; ps[:]=0.; mu[:]=1/m
            for i in range(m):
                ps[i,0,0]=R[n]; ps[i,1,1]=.03**2; ps[i,2,2]=.3**2
                if m==3: ps[i,3,3]=3.**2
        else: xs,ps,mu=imm_step(xs,ps,mu,z[n],R[n],A,Q,trans)
        means=np.sum(xs*H,axis=1); y=mu@means
        out[n]=y
        for i in range(m):
            var[n]+=mu[i]*(H[i]@ps[i]@H[i]+(means[i]-y)**2+future_var[i])
    return out,np.maximum(var,1e-18)

def seq_imm(model,e):
    return imm_sequence(e['z'],e['R'],e['reset'],model.A,model.Q,model.H,model.trans,model.future_var)

def windows(e):
    z=e['z']; X=np.zeros((len(z),W)); age=np.zeros(len(z),int)
    start=0
    for i in range(len(z)):
        if e['reset'][i]: start=i
        idx=np.maximum(np.arange(i-W+1,i+1),start)
        X[i]=z[idx]-z[i]; age[i]=i-start
    return X,age

class GP:
    """Fixed-window GP, SE + periodic covariance; unknown constant mean profiled out.
    Hyperparameters fitted on training windows, frozen thereafter. Linear functional
    is f(L)+tau f'(L); full posterior variance includes mean uncertainty.
    """
    def __init__(self,params):
        from scipy.linalg import cho_factor,cho_solve
        amp,ell,ap,period,lp,noise=np.exp(params)
        self.params=np.array(params)
        t=np.arange(-W+1,1)*DT
        def ker(d):
            return amp**2*np.exp(-d*d/(2*ell**2))+ap**2*np.exp(-2*np.sin(np.pi*d/period)**2/lp**2)
        K=ker(t[:,None]-t[None,:])+noise**2*np.eye(W)
        d=L-t
        se=amp**2*np.exp(-d*d/(2*ell**2)); per=ap**2*np.exp(-2*np.sin(np.pi*d/period)**2/lp**2)
        cross=se*(1-TAU*d/ell**2)+per*(1-TAU*2*np.pi/period*np.sin(2*np.pi*d/period)/lp**2)
        cf=cho_factor(K); ones=np.ones(W); invone=cho_solve(cf,ones); w=cho_solve(cf,cross)
        self.weights=w+invone*(1-ones@w)/(ones@invone)
        kss=amp**2+ap**2+TAU**2*(amp**2/ell**2+ap**2*(2*np.pi/period)**2/lp**2)
        self.var=max(kss-cross@w+(1-ones@w)**2/(ones@invone),1e-18)
    def batch(self,X,z):
        return X@self.weights+z,np.full(len(z),self.var)

class Net:
    """A 63-sample (0.496 s span) causal temporal convolution (63 -> h), tanh, 1x1 head.
    Applied only at the current tick. Gaussian mean/variance, no future padding.
    Residual level representation enforces translation equivariance.
    """
    def __init__(self,h=48,seed=204901):
        rng=np.random.default_rng(seed)
        self.w=[rng.normal(0,.04,(W,h)),np.zeros(h),rng.normal(0,.04,(h,2)),np.array([0.,-5.])]
    def forward(self,X):
        a=np.tanh(X@self.w[0]+self.w[1]); out=a@self.w[2]+self.w[3]
        return a,out
    def batch(self,X,z):
        _,out=self.forward(X*1000)
        return z+out[:,0]/1000,np.exp(np.clip(out[:,1],-14,8))/1e6
    def step(self,X,Y,opt,lr):
        a,out=self.forward(X); lv=np.clip(out[:,1],-14,8); err=out[:,0]-Y; inv=np.exp(-lv)
        g=np.column_stack([err*inv,.5*(1-err*err*inv)])/len(Y)
        g[:,1]*=(out[:,1]>-14)&(out[:,1]<8)
        ga=(g@self.w[2].T)*(1-a*a)
        grad=[X.T@ga,ga.sum(0),a.T@g,g.sum(0)]
        opt['t']+=1
        norm=np.sqrt(sum(np.sum(v*v) for v in grad)); scale=min(1.,10/max(norm,1e-20))
        for i,(w,gw) in enumerate(zip(self.w,grad)):
            gw=gw*scale; opt['m'][i]=.9*opt['m'][i]+.1*gw; opt['v'][i]=.999*opt['v'][i]+.001*gw*gw
            w-=lr*(opt['m'][i]/(1-.9**opt['t']))/(np.sqrt(opt['v'][i]/(1-.999**opt['t']))+1e-8)
        return float(np.mean(.5*(lv+err*err*inv)))

def load_models(path=None):
    import json
    cfg=json.loads((path or RES/'step3b_predictors.json').read_text(encoding='utf8'))
    m={'P1':MultiIMM(**cfg['selected']['P1']), 'P2':GP(cfg['selected']['P2'])}
    nn=Net(cfg['selected']['P3']['hidden']); nn.w=[np.array(x) for x in cfg['selected']['P3']['weights']]
    m['P3']=nn
    return cfg,m

class Adapter:
    """Proxy for original tracker: original drift_pred performs all measurement work.
    predict() returns the replacement eta_p, with sd=0 so the parent does not
    compensate servo lag twice. Legacy tracker continues for recovery velocity.
    """
    def __init__(self,legacy,model,gain=1.):
        self.legacy=legacy; self.model=model; self.gain=gain; self.buf=[]; self.last=None
    @property
    def initialized(self): return self.legacy.initialized
    @initialized.setter
    def initialized(self,v): self.legacy.initialized=v
    @property
    def R(self): return self.legacy.R
    @R.setter
    def R(self,v): self.legacy.R=v
    @property
    def x(self): return self.legacy.x
    def reset(self,z):
        self.buf=[z]; self.legacy.reset(z)
        if isinstance(self.model,MultiIMM): self.model.reset(z,self.R)
    def update(self,z):
        zz=z*self.gain
        self.buf.append(zz); self.buf=self.buf[-W:]
        self.legacy.update(z)
        if isinstance(self.model,MultiIMM): self.model.update(zz,self.R*self.gain**2)
    def predict(self,L):
        if isinstance(self.model,MultiIMM): y,v=self.model.output()
        else:
            X=np.pad(np.array(self.buf),(max(0,W-len(self.buf)),0),mode='edge')[-W:]
            y,v=self.model.batch((X-X[-1])[None,:],np.array([X[-1]])); y=float(y[0]);v=float(v[0])
        self.last=(y,v)
        return y,0.

class Replacement(variants.MotionAwareIMM):
    def __init__(self,model,gain=1.,**kw):
        self.new_model=model;self.measurement_gain=gain
        super().__init__(**kw)
    def reset(self,pl):
        super().reset(pl)
        self.tracker=Adapter(self.tracker,self.new_model,1.)

    def drift_pred(self,o,n):
        # Stress force/k_ref only. Scaling all eta would not exercise feedback.
        kg,kb=self.k_good,self.kbar
        self.k_good=kg/self.measurement_gain; self.kbar=kb/self.measurement_gain
        try: return super().drift_pred(o,n)
        finally: self.k_good,self.kbar=kg,kb

class GainIMM(variants.MotionAwareIMM):
    def __init__(self,gain=1.,**kw):
        self.measurement_gain=gain;super().__init__(**kw)
    def drift_pred(self,o,n):
        kg,kb=self.k_good,self.kbar
        self.k_good=kg/self.measurement_gain;self.kbar=kb/self.measurement_gain
        try: return super().drift_pred(o,n)
        finally: self.k_good,self.kbar=kg,kb

def build(name,p,models,gain=1.):
    if name=='Proposed':
        if gain==1: return variants.build('Proposed-IMM',CE.FINAL)
        return GainIMM(gain,**CE.FINAL)
    if name=='TrueK': return variants.build('Proposed-IMM',dict(CE.FINAL,oracle_fb=True,oracle_meas=True))
    if name=='Oracle+exc':
        return proposed.Oracle(p,excitation=True,a_force=CE.FINAL.get('a_force',.17),f_exc=CE.FINAL.get('f_exc',20.8),**{k:CE.FINAL[k] for k in ('M1','D1','KPF1','ki')})
    return Replacement(models[name],gain,**CE.FINAL)

def file_pool(module,func,jobs,cache,workers=8):
    """Process parallelism without Windows named pipes (sandbox denies mp.Pipe).
    All transport files are temporary in the caller's cache, never deliverables.
    """
    import subprocess,pickle,importlib
    cache=Path(cache).resolve();cache.mkdir(exist_ok=True,parents=True)
    processes=[];parts=[]
    for i in range(min(workers,len(jobs))):
        part=cache/f'{module}_{func}_{i}.task'
        output=cache/f'{module}_{func}_{i}.out'
        part.write_bytes(pickle.dumps((module,func,list(jobs[i::workers])),protocol=5))
        if output.exists(): output.unlink()
        log=(cache/f'{module}_{func}_{i}.log').open('w',encoding='utf8')
        proc=subprocess.Popen([sys.executable,'-X','utf8',str(Path(__file__).resolve()),'--worker',str(part),str(output)],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        processes.append((proc,log));parts.append(output)
    while any(p.poll() is None for p,_ in processes): time.sleep(.5)
    for p,log in processes:
        log.close()
        if p.returncode: raise RuntimeError(f'file worker failed: see {log.name}')
    shards=[pickle.loads(p.read_bytes()) for p in parts]
    # restore original order, so identical seeds are byte-reproducible across workers
    result=[]
    for i in range(len(jobs)): result.append(shards[i%len(shards)][i//len(shards)])
    return result

if __name__=='__main__' and len(sys.argv)>1 and sys.argv[1]=='--worker':
    import pickle,importlib
    module,func,jobs=pickle.loads(Path(sys.argv[2]).read_bytes())
    fn=getattr(importlib.import_module(module),func)
    rows=[]
    for i,job in enumerate(jobs):
        rows.append(fn(job))
        if i%10==0: print('worker',i,'/',len(jobs),flush=True)
    Path(sys.argv[3]).write_bytes(pickle.dumps(rows,protocol=5))
