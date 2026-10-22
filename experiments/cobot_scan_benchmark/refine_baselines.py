"""Train/dev-only baseline repair and retuning; frozen-candidate physical comparison.

Existing sources, weights and result files remain unchanged. ``tune`` selects
parameters using seeds 1000/2000 only; ``test`` uses the old shared 3000/4000
cohorts after selection has been saved. Adaptation is not exact reproduction.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[_k]='1'
import argparse, csv, hashlib, json, math, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from numba import njit

ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0,str(ROOT))
import plant as P
import scenarios as S
import baselines_mpc as B
import baselines_adaptive as A
import baselines_learned as L
import report_physical_metrics as R
OUT=ROOT/'refined_baselines'
CELLS=R.CELLS

def read(name): return json.loads((ROOT/name).read_text(encoding='utf-8'))
def clean(value):
    if isinstance(value,dict): return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [clean(v) for v in value]
    if isinstance(value,(float,np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value,np.integer): return int(value)
    return value
def save(name,obj):
    OUT.mkdir(exist_ok=True)
    (OUT/name).write_text(json.dumps(clean(obj),indent=1,allow_nan=False),encoding='utf-8')

class XueCorrected(B.XueMPVIC):
    """HC + mode-switch/tank + MPC, adapted to the measured position servo.

    Fixes: active K branch; constrained HC parameter vector (not just a local
    copy); tank/stiffness initialization as stated before Eq 28; anti-overflow
    stiffness projection; measured-state reset for a position interface whose
    commanded acceleration is not instantaneously realized. Optional causal
    contact reanchoring is used only if the original HC coordinate becomes
    nonpositive. Initial surface prior remains the original A1 prior.
    Servo lead is a bounded, tunable position-interface compensation, not a
    passivity guarantee. No true stiffness or future surface position is read.
    """
    def __init__(self,state_reset=True,servo_lead=0.,reanchor=True,est_gain=1.,ck_max=.01,**kw):
        super().__init__(**kw)
        self.state_reset=state_reset; self.servo_lead=servo_lead
        self.reanchor=reanchor; self.est_gain=est_gain; self.ck_max=ck_max
    def _reset(self):
        super()._reset()
        self.Kbar=self.kmax; self.kun_prev=self.kmax
        self._yold=None; self._vel=0.; self._nanchor=0
    def _decide(self,o):
        if o.n==0: return self.BUILTIN+(self.KPF,0.)
        y=float(o.y[-1]); F=float(o.force[-1]); Tc=self.Tc
        if self.yref is None: self.yref=float(o.y[0])
        if self._yold is not None:
            self._vel += (1-math.exp(-Tc/.016))*((y-self._yold)/Tc-self._vel)
        self._yold=y
        if self.state_reset:
            self.e=y-self.yref; self.ed=float(np.clip(self._vel,-self.VMAX,self.VMAX)); self.edd=0.
        if self.th is not None:
            self.th[0]=np.clip(self.th[0],-12.,18.)
            self.th[1]=np.clip(self.th[1],.5,2.5)
            self.th[2]=np.clip(self.th[2],0.,self.ck_max)
        if self.reanchor and self.ys0 is not None and (y-self.ys0)*1e3<=self.zmin_mm and F>.5:
            # Retain the estimated HC stiffness while translating the contact
            # frame. A force-conditioned coordinate is causal; no scene prior
            # is refreshed from the hidden moving surface.
            k=math.exp(float(self.th[0])) if self.th is not None else 1.
            nn=float(self.th[1]) if self.th is not None else self.n0
            z=float(np.clip((max(F,.1)/max(k,1e-6))**(1/max(nn,.5)),self.zmin_mm,15.))
            self.ys0=y-z*1e-3; self.phi_prev=None; self._nanchor+=1
        if self.ys0 is None:
            self.ys0=self.ys0_scene if self.ys0_mode=='A1' and self.ys0_scene is not None else float(o.y[0])-B._znom_mm()*1e-3
        z=(y-self.ys0)*1e3
        zd=(y-float(o.y[-2]))*1e3/Tc if o.n>1 else 0.
        valid=z>self.zmin_mm; zc=max(z,self.zmin_mm)
        phi=np.array([1.,math.log(zc),zd])
        if self.th is None: self.th=np.array([math.log(max(F,.1))-self.n0*math.log(zc),self.n0,0.])
        if self.phi_prev is not None and F>.05 and self.phi_valid:
            ey=math.log(F)-float(self.th@self.phi_prev)
            if abs(ey)>2*self.Delta:
                self.th+=self.est_gain*self.phi_prev*ey/(1+float(self.phi_prev@self.phi_prev)); self.nupd+=1
            self.nest+=1
        self.nfreeze+=int(not valid); self.phi_prev=phi; self.phi_valid=valid
        self.th[0]=np.clip(self.th[0],-12.,18.); self.th[1]=np.clip(self.th[1],.5,2.5); self.th[2]=np.clip(self.th[2],0.,self.ck_max)
        kk=math.exp(self.th[0]); nn=self.th[1]; cc=self.th[2]*kk
        ef=P.FD-F
        if self.efp is not None: self.efd+=(1-math.exp(-Tc/self.fdot_tau))*((ef-self.efp)/Tc-self.efd)
        self.efp=ef
        mu=(B._sig(ef,-self.w1,0.)*(1-B._sig(ef,P.FD-self.w1,P.FD))*
            B._sig(self.efd,-self.efb,-self.efb+self.w2)*(1-B._sig(self.efd,self.efb-self.w2,self.efb)))
        self.mu=mu; d=self.dmax+mu*(self.dmin-self.dmax); kun=mu*self.kmax
        self.Ie=float(np.clip(self.Ie+ef*Tc,-20.,20.)); fr=P.FD+self.ki*self.Ie
        ez=self.e; ezd=self.ed; ek=kun-self.Kbar; kund=(kun-self.kun_prev)/Tc; self.kun_prev=kun
        if abs(ez)<=self.eps: kbd=kund+self.lam*ek
        else:
            p1=self.lam*ek*ez*ez/2 if self.Vtk>0 else 0.
            kbd=min(kund+self.lam*ek,2*(d*ezd*ezd-fr*ezd+p1)/(ez*ez))
        # Enforce the actual discrete change in the energy ledger, not an
        # unclipped derivative followed by a different applied stiffness.
        knew=float(np.clip(self.Kbar+kbd*Tc,0.,self.kmax)); kbd=(knew-self.Kbar)/Tc; self.Kbar=knew
        self.Vtk=max(self.Vtk+((d-self.p2)*ezd*ezd-fr*ezd-.5*kbd*ez*ez)*Tc,0.)
        kc=max(kk+cc*zd,0.); Fhat=kc*zc**nn
        F0=Fhat if self.f0src=='model' else F
        kz=nn*kc*zc**(nn-1)*1e3; cz=cc*zc**nn*1e3
        eact=y-self.yref
        a1=B._xue_mpc(self.M,d+cz,self.Kbar+kz,F0-fr-kz*eact-cz*zd*1e-3,
                      self.e,self.ed,self.edd,Tc,self.hp,self.hc,self.amax)
        if self.state_reset:
            # A current measured-state reset must execute the new acceleration;
            # resetting a jerk chain every tick would otherwise discard it.
            self.e+=Tc*self.ed+.5*Tc*Tc*a1
            self.ed=float(np.clip(self.ed+Tc*a1,-self.VMAX,self.VMAX)); self.edd=a1
        else:
            self.e+=Tc*self.ed; self.ed=float(np.clip(self.ed+Tc*self.edd,-self.VMAX,self.VMAX)); self.edd=a1
        Mb,Db,Kb,kib=self.BUILTIN
        ee=math.exp(-Db*Tc/Mb); vbn=self.vb*ee+(1-ee)*ef/Db
        ycb_next=self.ycb+.5*(self.vb+vbn)*Tc
        lead=float(np.clip(self.servo_lead*self.tau*self.ed,-.003,.003))
        uff=self.yref+self.e+lead-ycb_next
        if o.n%4==0: self.hist.append((o.t[-1],d,self.Kbar,mu,kk,nn,cc/kk,Fhat,z))
        return (Mb,Db,Kb,kib,self.KPF,uff)

class JinCorrected(B.JinMPC):
    """Retain Jin OS/MS optimization; include an action-visible OS endpoint.

    At per=1 the previous ndel+per endpoint precedes the command reaching the
    servo. Extra ticks keep a single held decision but make its force effect
    visible. MS is unchanged. This is an interface correction, not future truth.
    """
    def __init__(self,lookahead=1,**kw):
        super().__init__(**kw)
        self.lookahead=lookahead

    def _decide(self, o):
        n = o.n
        if n < 3 or self.g is None:
            if self.per is None and n > 0: pass
            self._last_n = n
            return (self.M, math.exp(self.ld), self.K, self.ki, self.kp)
        if self.per is None: self.per = max(n - self._last_n, 1)
        per = self.per
        if self.pr is None:
            self.pr = np.array([self.M, self.K, self.ki, self.yref, self.beta, self.af, self.Tc, self.VMAX,
                                1.0, self.qv, 0.01, self.flim, self.alpha_b, self.delta_b])
        f = o.force
        Fbar = f[n - 2] + (f[n - 1] - f[n - 2]) / self.af
        kh = self.est.k * 1e3 * self.kgain          # train-calibrated estimator gain (closed-loop bias)
        st, q = self._pred_state(o, Fbar, kh, self.vs if self.dist else 0.0)
        ldp, kpp = self.ld, self.kp
        if self.mode == 'os':
            x = np.array([ldp, kpp])
            B._os_solve(x, self.ndel + per + self.lookahead, st, q, self.pr, self.rd, self.rc, ldp, kpp, self.use_kpf, kpp, 2,
                      self.rateD, self.rateK, self.ldr, self.kpr)
            ld, kp = float(x[0]), float(x[1])
        else:
            H = self.H
            if self.z is None:
                zd = math.atanh(min(max(2 * (ldp - B.DLO) / (B.DHI - B.DLO) - 1, -0.999), 0.999))
                zk = math.atanh(min(max(2 * kpp / B.KPMAX - 1, -0.999), 0.999))
                self.z = np.array([zd] * H + [zk] * H)
            else:
                self.z[:H - 1] = self.z[1:H].copy(); self.z[H:2 * H - 1] = self.z[H + 1:].copy()
            B._ms_solve(self.z, self.iters, H, per, st, q, self.pr, self.rd, self.rc, ldp, kpp,
                      self.use_kpf, kpp, self.ldr, self.kpr)
            ld = B.DLO + (B.DHI - B.DLO) * 0.5 * (math.tanh(self.z[0]) + 1)
            kp = B.KPMAX * 0.5 * (math.tanh(self.z[H]) + 1) if self.use_kpf else 0.0
            # rate box of Omega_t on the applied move (the planned sequence is penalised by rd)
            ld = min(max(ld, ldp - self.rateD), ldp + self.rateD)
            kp = min(max(kp, kpp - self.rateK), kpp + self.rateK)
        # tank gating of energy-injecting KPF increases (E41 eq. 41-42)
        if self.tank and kp > kpp:
            ef = B.FD - f[n - 1]; inj = 0.5 * ef * ef * (kp - kpp)
            self.nkup += 1
            if self.G - inj < self.Gmin:
                kp = kpp; self.ngate += 1
            else:
                self.G -= inj
        self.ld, self.kp = ld, kp
        self.hist.append((o.t[n - 1], math.exp(ld), kp, self.est.k))
        return (self.M, math.exp(ld), self.K, self.ki, kp)


@njit(cache=True)
def learned_step(log,n,th,base,sc,zprev,g):
    # Same causal feature channels as the previous E27-style adaptation.
    nf=12; nh=8; no=3; w=6
    x=np.zeros(nf)
    if n:
        for k in range(4):
            hi=max(n-k*w,1); lo=min(max(n-(k+1)*w,0),hi-1)
            x[k]=P.FD-np.mean(log[lo:hi,5])
        j=max(n-1-w,0); j4=max(n-1-4*w,0)
        x[4]=(log[n-1,5]-log[j,5])/(w*.008)/20.
        x[5]=(log[n-1,2]-log[j,2])/(w*.008)/.01
        x[6]=(log[n-1,11]-log[n-1,2])/.001
        x[7]=np.std(log[max(n-2*w,0):n,5])/.3
        x[8]=(math.log(max(log[n-1,9],1e-9))-math.log(base[1]))/2.
        x[9]=(math.log(max(log[n-1,10],1e-9))-math.log(max(base[4],1e-6)))/2.
        x[10]=(log[n-1,2]-log[j4,2])/.001
        x[11]=log[n-1,5]-log[j4,5]
    h=np.zeros(nh)
    for a in range(nh):
        v=th[nf*nh+a]
        for k in range(nf): v+=th[a*nf+k]*min(max(x[k],-5.),5.)
        h[a]=math.tanh(v)
    off=nf*nh+nh
    for a in range(no):
        v=th[off+no*nh+a]
        for k in range(nh): v+=th[off+a*nh+k]*h[k]
        zz=math.tanh(v)
        zprev[a] += sc[0]*(zz-zprev[a])
    const=off+no*nh+no
    g[0]=min(max(base[0]*math.exp(.7*min(max(th[const],-4.),4.)),.25),3.)
    g[1]=min(max(base[1]*math.exp(2*zprev[0]+min(max(th[const+1],-4.),4.)),20.),2500.)
    g[2]=base[2]
    g[3]=min(max(base[3]*math.exp(.7*min(max(th[const+2],-4.),4.)),0.),30.)
    g[4]=min(max(base[4]*math.exp(2*zprev[1]+min(max(th[const+3],-4.),4.)),1e-6),1e-3)
    # Reference correction is restored, filtered and rate-bounded. It is not
    # labelled from hidden stiffness or the proposed method.
    du=min(max(sc[1]*zprev[2]-g[5],-sc[2]),sc[2])
    g[5]+=du

class LearnedCorrected:
    def __init__(self,theta,base,period=.008,tau=.04,ref_span=.0005):
        self.theta=np.asarray(theta,float); self.base=np.asarray(base,float)
        self.period=period; self.tau=tau; self.ref_span=ref_span
    def reset(self,pl):
        self.sc=np.array([1-math.exp(-self.period/self.tau),self.ref_span,.0001*self.period/.008])
        self.z=np.zeros(3); self.g=np.r_[self.base,0.]
    def __call__(self,o):
        learned_step(o.log,o.n,self.theta,self.base,self.sc,self.z,self.g)
        return self.g

@njit(cache=True)
def learned_rollout(s,p,theta,base,sc,noise,nt,NC,ndel,tau,fc,vmax,dt):
    """Same tick order/plant as P.episode; eliminate Python callbacks in training."""
    lg=np.zeros((nt+2,P.LOGC)); zz=np.zeros(3); gg=np.zeros(6); gg[:5]=base
    n=0
    while n<nt:
        learned_step(lg,n,theta,base,sc,zz,gg)
        for k in range(6): gg[k]=min(max(gg[k],P.GMIN[k]),P.GMAX[k])
        n=P.advance(s,p,gg,1,NC,ndel,tau,fc,vmax,dt,noise,lg,n)
        if n<0: return lg,False
    return lg[:n],True

def fast_learned_episode(p,pol,T,seed):
    pl=P.PLATFORMS['PC125']; dt=.0005; nt=int(round(T/pl['Tc']))
    pol.reset(pl)
    noise=np.random.default_rng(seed).standard_normal((int(T*1000)+10,2))
    lg,ok=learned_rollout(P.init_state(p,dt),p,pol.theta,pol.base,pol.sc,noise,nt,
                         int(round(pl['Tc']/dt)),int(pl['ndel']),pl['tau_a'],pl['fc'],pl['VMAX'],dt)
    return lg if ok else None


NTH=12*8+8+8*3+3+4
def init_theta(seed):
    z=np.zeros(NTH); z[:96]=np.random.default_rng(seed).standard_normal(96)/math.sqrt(12)
    return z

def make(spec):
    kind,cfg=spec
    if kind=='xue': return XueCorrected(**cfg),None
    if kind=='jin':
        cfg=dict(cfg); update=cfg.pop('update',.016)
        return JinCorrected(**cfg),update
    if kind=='duan': return A.duan_gains(cfg['m'],cfg['b'],cfg['sigma']),None
    if kind=='song': return A.Song(**cfg),None
    if kind=='learn': return LearnedCorrected(**cfg),cfg.get('period',.008)
    if kind=='old_mpc': return B.CLASSES[cfg[0]](**cfg[1]),B.UPDATE[cfg[0]]
    if kind=='old_learn': return L.LearnedVIC(np.asarray(cfg['theta'])),L.UPDATE
    raise ValueError(kind)

def evaluate(job):
    task,spec,seed0,n,physical=job
    ctrl,update=make(spec)
    if physical or spec[0]!='learn':
        pp,_=S.cohort(task,1,1000,CELLS[task])[0]
        if hasattr(ctrl,'bind'): ctrl.bind(pp)
        P.episode(pp,ctrl,T=.2,seed=1000,update=update)
    rows=[]
    for scene,(p,T) in enumerate(S.cohort(task,n,seed0,CELLS[task])):
        ctrl,update=make(spec)
        if hasattr(ctrl,'bind'): ctrl.bind(p)
        timed=R.TimedPolicy(ctrl) if callable(ctrl) and physical else ctrl
        lg=fast_learned_episode(p,ctrl,T,seed0+scene) if spec[0]=='learn' and not physical else P.episode(p,timed,T=T,seed=seed0+scene,update=update)
        m=S.metrics(lg)
        row={k:float(v) for k,v in m.items()}
        if physical:
            row.update(R.physical(lg,T)); row.update(scene=scene,noise_seed=seed0+scene)
            tt=np.array(timed.times)*1000 if callable(ctrl) else np.array([])
            row['policy_median_ms']=float(np.median(tt)) if len(tt) else None
            row['policy_p95_ms']=float(np.percentile(tt,95)) if len(tt) else None
            row['policy_max_ms']=float(np.max(tt)) if len(tt) else None
        rows.append(row)
    js=[r['J'] if np.isfinite(r['J']) else 1e12 for r in rows]
    pk=[r['peak'] if np.isfinite(r['peak']) else 1e9 for r in rows]
    return dict(J=float(np.mean(js)),feasible=float(np.mean(np.array(pk)<=7))>=.95,
                rmse=float(np.nanmean([r['rmse'] for r in rows])),
                peak=float(np.nanmean(pk)),rows=rows if physical else None)

def rank(r): return (not r['feasible'],r['J'])
def run_many(pool,jobs): return list(pool.map(evaluate,jobs,chunksize=1))

def search_space(task,kind):
    rng=np.random.default_rng(2261022+(task=='T1')+10*('xue'==kind))
    mp=next(c for c in read('mpc_results.json')['cells'] if c['task']==task)
    ad=next(c for c in read('adaptive_results.json')['cells'] if c['task']==task)['methods']
    if kind in ('Jin_OS','Jin_MS'):
        old=mp[kind]['params']; configs=[dict(old,update=.016)]
        # Explicit boundary extensions plus joint retuning of the fixed part.
        for rc in (.01,.1,.3,1.,3.,10.,30.,100.,300.):
            for qv in (0.,.1): configs.append(dict(old,rc=rc,qv=qv,update=.008))
        for _ in range(65):
            c=dict(old,rc=float(10**rng.uniform(-2,2.5)),rd=float(10**rng.uniform(-3,1)),
                   qv=float(rng.choice([0.,.01,.1,1.])),kgain=float(rng.choice([.7,1.,1.7,2.5])),
                   tau_k=float(rng.choice([.03,.1,.3,.6])),update=float(rng.choice([.008,.016])),
                   M=float(rng.choice([.3,1.,old['M'],3.])),K=float(rng.choice([0.,old['K']])),
                   ki=float(rng.choice([3.,10.,30.])),H=int(rng.choice([5,10,15])),
                   iters=int(rng.choice([4,8,12])),rateD=float(rng.choice([.15,.4,1.])),
                   rateK=float(rng.choice([1e-5,5e-5,1e-4])))
            configs.append(c)
        return [('jin',dict(c,lookahead=int(rng.choice([1,3,6]))) if kind=='Jin_OS' else c) for c in configs]
    if kind=='Xue':
        old=mp['Xue_A1']['params']
        configs=[]
        for reset in (False,True):
            for lead in (0.,.5,1.):
                for k in (30.,300.,1000.):
                    configs.append(dict(old,kmax=k,state_reset=reset,servo_lead=lead))
        for _ in range(78):
            mass=float(rng.choice([.3,1.,2.,3.]))
            dmin=float(rng.choice([20.,50.,100.,200.])); dmax=float(rng.choice([100.,200.,500.,1000.,2000.]))
            configs.append(dict(old,M=mass,dmin=min(dmin,dmax),dmax=dmax,kmax=float(rng.choice([10.,30.,100.,300.,1000.])),
                ki=float(rng.choice([0.,3.,10.,30.])),KPF=float(rng.choice([0.,2e-5,1e-4,3e-4])),
                Delta=float(rng.choice([.01,.03,.1])),efb=float(rng.choice([2.,10.,40.,100.])),
                n0=float(rng.choice([1.,1.3,1.5])),hp=int(rng.choice([5,10])),
                state_reset=True,servo_lead=float(rng.choice([0.,.5,1.])),reanchor=bool(rng.integers(0,2)),
                est_gain=float(rng.choice([.1,.3,1.])),f0src=str(rng.choice(['model','measured']))))
        return [('xue',c) for c in configs]
    if kind=='Duan':
        old=ad['Duan2018']['tuning']; configs=[dict(m=old['m'],b=old['b'],sigma=old['sigma'])]
        for _ in range(95):
            m=float(np.exp(rng.uniform(np.log(.25),np.log(3.)))); b=float(np.exp(rng.uniform(np.log(20),np.log(2500))))
            ki=min(float(rng.uniform(0,30)),.999*A.duan_sigma_max(m,b)/.008)
            configs.append(dict(m=m,b=b,sigma=ki*.008))
        return [('duan',c) for c in configs]
    if kind=='Song':
        q=ad['Song2025']['tuning']['mapped']; configs=[{k:q[k] for k in ('M','D','ki','KPF')}]
        for _ in range(95):
            m=float(np.exp(rng.uniform(np.log(.25),np.log(3.)))); d=float(np.exp(rng.uniform(np.log(20),np.log(2500))))
            kp=float(np.exp(rng.uniform(np.log(1e-6),np.log(1e-3))))
            ki=min(float(rng.uniform(0,30)),.249/max(kp*d,1e-12),.999*d/m)
            configs.append(dict(M=m,D=d,ki=ki,KPF=kp))
        return [('song',c) for c in configs]
    raise ValueError(kind)

def tune(pool):
    data=read_refined('tuning.json',{})
    for task in CELLS:
        for family in ('Xue','Jin_OS','Jin_MS','Duan','Song'):
            key=task+'_'+family
            if key in data and (family!='Jin_OS' or data[key].get('adaptation_version')==2): continue
            specs=search_space(task,family)
            screened=run_many(pool,[(task,s,1000,8,False) for s in specs])
            ix=sorted(range(len(specs)),key=lambda i:rank(screened[i]))[:10]
            if 0 not in ix: ix.append(0)
            full=run_many(pool,[(task,specs[i],1000,20,False) for i in ix])
            tops=sorted(range(len(ix)),key=lambda i:rank(full[i]))[:3]
            dev=run_many(pool,[(task,specs[ix[i]],2000,20,False) for i in tops])
            # Feasibility is checked on train; among its top three, dev picks
            # the final controller. No test cohort is touched here.
            order=sorted(range(len(tops)),key=lambda j:(not full[tops[j]]['feasible'],not dev[j]['feasible'],dev[j]['J']))
            chosen=order[0]; best=specs[ix[tops[chosen]]]
            data[key]=dict(adaptation_version=2 if family=='Jin_OS' else 1,spec=best,train=full[tops[chosen]],dev=dev[chosen],
                           n_screen=len(specs),screen_scenes=8,full_scenes=20,n_full=len(ix),n_dev=len(tops),
                           screen=screened,full_indices=ix,full=full,dev_indices=[ix[i] for i in tops],dev_all=dev)
            save('tuning.json',data)
            print('SELECTED',key,'train',data[key]['train']['J'],'dev',dev[chosen]['J'],'feasible',data[key]['train']['feasible'],best,flush=True)
    return data

def read_refined(name,default):
    p=OUT/name
    return json.loads(p.read_text()) if p.exists() else default

def train_learning(pool,gens=100):
    records=read_refined('learning.json',{})
    for task in CELLS:
        base=next(c for c in read('gate_PC125.json') if c['task']==task and c['dose']==CELLS[task])['B1']['gains']
        for seed in range(3):
            key=f'{task}_{seed}'
            if key in records: continue
            es=L.SepCMA(init_theta(seed),.025,16,seed)
            def spec(th): return ('learn',dict(theta=th.tolist(),base=base,period=.008,tau=.04,ref_span=.0005))
            start=evaluate((task,spec(es.m),1000,20,False))
            dstart=evaluate((task,spec(es.m),2000,20,False))
            best=((not start['feasible'],not dstart['feasible'],dstart['J']),es.m.copy(),0,start,dstart)
            hist=[]
            for g in range(gens):
                X=es.ask()
                # Training seeds remain in [1000,2000); every generation has
                # common scenes/noise for all population members.
                seed0=1100+g
                scores=run_many(pool,[(task,spec(x),seed0,6,False) for x in X])
                es.tell(X,np.array([v['J'] for v in scores]))
                if (g+1)%10==0:
                    cands=[es.m.copy(),X[int(np.argmin([v['J'] for v in scores]))].copy()]
                    checks=run_many(pool,[(task,spec(x),1000,20,False) for x in cands])
                    dchecks=run_many(pool,[(task,spec(x),2000,20,False) for x in cands])
                    for x,v,dv in zip(cands,checks,dchecks):
                        rr=(not v['feasible'],not dv['feasible'],dv['J'])
                        if rr<best[0]: best=(rr,x.copy(),g+1,v,dv)
                    hist.append(dict(gen=g+1,train=best[3],dev=best[4],sigma=float(es.sigma)))
                    print('LEARN',task,seed,g+1,'train',best[3]['J'],'dev',best[4]['J'],'feas',best[3]['feasible'],flush=True)
            dev=best[4]
            np.savez(OUT/f'learned_{task}_seed{seed}.npz',theta=best[1],base=base)
            records[key]=dict(spec=spec(best[1]),train=best[3],dev=dev,best_gen=best[2],gens=gens,
                              population=16,scenes_per_generation=6,history=hist,
                              episodes=gens*16*6+40+(gens//10)*80,
                              status='finite declared budget; learning curve retained, no global convergence claim')
            save('learning.json',records)
    return records

FAMILIES=('Xue','Xue_nominal','Jin_OS','Jin_MS','Duan','Song')

def choose(pool,task,specs):
    tr=run_many(pool,[(task,sp,1000,20,False) for sp in specs])
    ix=sorted(range(len(specs)),key=lambda i:rank(tr[i]))[:3]
    if 0 not in ix: ix.append(0)
    dv=run_many(pool,[(task,specs[i],2000,20,False) for i in ix])
    j=min(range(len(ix)),key=lambda j:(not tr[ix[j]]['feasible'],not dv[j]['feasible'],dv[j]['J']))
    return dict(spec=specs[ix[j]],train=tr[ix[j]],dev=dv[j],n_full=len(specs),
                full=tr,dev_indices=ix,dev_all=dv)

def local_space(spec):
    kind,cfg=spec; cfg=dict(cfg); out=[(kind,cfg)]
    if kind=='xue':
        axes={'M':(.25,3.),'dmin':(20.,2000.),'dmax':(20.,2500.),'kmax':(1.,5000.),
              'ki':(.1,30.),'KPF':(1e-6,1e-3),'est_gain':(.03,1.)}
    elif kind=='jin':
        axes={'M':(.25,3.),'K':(0.,1000.),'ki':(.1,30.),'D0':(20.,2500.),
              'kpf0':(1e-6,1e-3),'rc':(.001,3000.),'rd':(.0001,100.)}
    elif kind=='duan': axes={'m':(.25,3.),'b':(20.,2500.),'sigma':(.0001,.24)}
    elif kind=='song': axes={'M':(.25,3.),'D':(20.,2500.),'ki':(.1,30.),'KPF':(1e-6,1e-3)}
    else: raise ValueError(kind)
    for a,(lo,hi) in axes.items():
        for factor in (.5,2.):
            c=dict(cfg); c[a]=float(np.clip(c.get(a,1.)*factor,lo,hi))
            if kind=='xue' and c['dmin']>c['dmax']: continue
            if kind=='duan': c['sigma']=min(c['sigma'],.999*A.duan_sigma_max(c['m'],c['b']),30*.008)
            if kind=='song': c['ki']=min(c['ki'],.249/max(c['KPF']*c['D'],1e-12),.999*c['D']/c['M'])
            out.append((kind,c))
    if kind=='xue':
        # Both coordinate strategies are eligible; no post-test switching.
        for a in ('state_reset','reanchor'):
            out.append((kind,dict(cfg,**{a:not cfg.get(a,True)})))
    if kind=='jin':
        if cfg['mode']=='os':
            for h in (1,3,6): out.append((kind,dict(cfg,lookahead=h)))
        else:
            out.append((kind,dict(cfg,H=15,iters=12)))
    return list({json.dumps(x,sort_keys=True):x for x in out}.values())

def polish(pool):
    coarse=read_refined('tuning.json',{}); done=read_refined('selection.json',{})
    for task in CELLS:
        for fam in FAMILIES:
            key=task+'_'+fam
            if key in done: continue
            if fam=='Xue_nominal':
                specs=[(kind,dict(c,ys0_mode='nominal')) for kind,c in search_space(task,'Xue')]
                scr=run_many(pool,[(task,sp,1000,8,False) for sp in specs])
                ix=sorted(range(len(specs)),key=lambda i:rank(scr[i]))[:10]
                best=choose(pool,task,[specs[i] for i in ix])
                start=dict(best,screen=scr,n_screen=len(specs))
            else:
                best=coarse[key]; start=best
            hist=[]
            for rnd in range(2):
                nxt=choose(pool,task,local_space(best['spec']))
                hist.append(nxt); best=nxt
            done[key]=dict(spec=best['spec'],train=best['train'],dev=best['dev'],
                           coarse=start,local_rounds=hist)
            save('selection.json',done)
            print('POLISHED',key,'train',best['train']['J'],'dev',best['dev']['J'],best['spec'],flush=True)
    return done


def calibrate_learning(pool):
    records=read_refined('learning.json',{})
    assert len(records)==6
    for key,r in list(records.items()):
        if 'output_calibration' in r: continue
        task=key.split('_')[0]; kind,cfg=r['spec']
        specs=[r['spec']]
        specs += [(kind,dict(cfg,tau=tau,ref_span=span)) for tau in (.016,.04,.08)
                  for span in (.00025,.0005,.001,.003,.01)]
        specs=list({json.dumps(x,sort_keys=True):x for x in specs}.values())
        selected=choose(pool,task,specs)
        records[key]=dict(r,spec=selected['spec'],train=selected['train'],dev=selected['dev'],
                         before_calibration=dict(spec=r['spec'],train=r['train'],dev=r['dev']),
                         output_calibration=selected,
                         episodes=r['episodes']+20*(len(specs)+len(selected['dev_indices'])))
        save('learning.json',records)
        print('CALIBRATED',key,'dev',selected['dev']['J'],'tau',selected['spec'][1]['tau'],
              'span',selected['spec'][1]['ref_span'],flush=True)
    return records


def selected_checks():
    results={}
    for key,rec in read_refined('selection.json',{}).items():
        kind,c=rec['spec']; r={}
        if kind=='duan':
            cap=A.duan_sigma_max(c['m'],c['b']); assert 0<c['sigma']<cap and c['sigma']/.008<=30.+1e-8
            r=dict(sigma=c['sigma'],sigma_stability_upper_bound=cap,ki=c['sigma']/.008)
        elif kind=='song':
            prod=c['KPF']*c['D']*c['ki']; margin=c['D']-c['M']*c['ki']
            assert prod<=.25+1e-12 and margin>0
            alpha=.5*(1+math.sqrt(1-4*prod)); beta=c['ki']/alpha; gamma=c['KPF']*c['D']
            r=dict(native_md=c['M'],native_bd=c['D'],native_alpha=alpha,native_beta=beta,native_gamma=gamma,
                   realizability_product=prod,stability_margin=margin)
            assert abs(alpha+gamma*beta-1)<1e-10
        elif kind=='xue':
            assert c['kmax']>0 and c['dmax']>=c['dmin']>0
            task=key.split('_')[0]; p,_=S.cohort(task,1,1000,CELLS[task])[0]
            pol,_=make(rec['spec']); pol.bind(p)
            log=P.episode(p,pol,T=2.,seed=1000)
            hist=np.asarray(pol.hist)
            assert log is not None and np.isfinite(log).all() and np.isfinite(pol.th).all()
            assert hist[:,2].max()>0 and hist[:,2].min()>=0 and hist[:,2].max()<=c['kmax']+1e-8
            r=dict(stiffness_min=float(hist[:,2].min()),stiffness_max=float(hist[:,2].max()),
                   replica_max_error=pol.maxrep,projection_parameter_vector=pol.th.tolist())
        elif kind=='jin':
            assert 0<=c.get('K',0)<=P.GMAX[2]
            assert c['mode']!='os' or c.get('lookahead',1)>=1
            r=dict(update_s=c['update'],extra_endpoint_ticks=c.get('lookahead') if c['mode']=='os' else None)
        results[key]=r
    save('selected_checks.json',dict(passed=True,controllers=results))
    return results


def source_files():
    names=['plant.py','scenarios.py','refine_baselines.py','report_physical_metrics.py',
           'baselines_mpc.py','baselines_adaptive.py','baselines_learned.py',
           'candidate1_model.py','candidate1_v2.py','candidate1_v3.py','candidate2_learned.py',
           'candidate1_results.json','candidate1_v2_results.json','gate_PC125.json',
           'physical_results/per_scene.csv','refined_baselines/selection.json','refined_baselines/learning.json']
    names += [f'candidate2_full_seed{i}.npz' for i in range(3)]
    names += [f'refined_baselines/learned_{t}_seed{i}.npz' for t in CELLS for i in range(3)]
    return names

def freeze():
    import datetime,subprocess
    assert len(read_refined('selection.json',{}))==12
    learning=read_refined('learning.json',{})
    assert len(learning)==6 and all('output_calibration' in r for r in learning.values())
    assert read_refined('checks.json',{})['passed']
    selected_checks()
    hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in source_files()}
    protocol=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  train_seeds=1000,dev_seeds=2000,test_cohorts=[3000,4000],n_scenes=40,tasks=CELLS,
                  neural_seeds=[0,1,2],scene_noise_pairing=True,candidates_retuned=False,
                  test_status='historically inspected cohorts; no untouched-test claim',
                  force_metrics='true force, original 1 s warm-up, 8 ms samples',
                  environment=dict(python=sys.version,numpy=np.__version__,command=sys.argv),hashes=hashes)
    previous=read_refined('protocol.json',None)
    if previous:
        assert previous['hashes']==hashes,'Frozen source/configuration changed; inspect before any new test'
        return previous
    save('protocol.json',protocol)
    return protocol

def verify_frozen(protocol):
    for n,h in protocol['hashes'].items():
        assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n


def final_test(pool):
    protocol=freeze(); verify_frozen(protocol)
    tune_data=read_refined('selection.json',{}); learn=read_refined('learning.json',{})
    rows=[]
    for cohort in (3000,4000):
        for task in CELLS:
            specs=[(family,tune_data[task+'_'+family]['spec'],-1) for family in FAMILIES]
            specs += [('Learned_VIC',learn[f'{task}_{seed}']['spec'],seed) for seed in range(3)]
            res=run_many(pool,[(task,spec,cohort,40,True) for _,spec,_ in specs])
            for (name,_,seed),v in zip(specs,res):
                for r in v['rows']:
                    r.update(method=name,training_seed=seed,task=task,cohort=cohort)
                    rows.append(r)
                print('TEST',cohort,task,name,seed,'J',v['J'],'RMSE',v['rmse'],flush=True)
    verify_frozen(protocol)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    rows=[{k:r.get(k) for k in fields} for r in rows]
    R.write_csv(OUT/'per_scene.csv',rows)
    save('test_records.json',rows)
    return rows

def checks():
    # Paper Eq 28 initialization and both stiffness modes must be reachable.
    p,T=S.cohort('T1',1,1000,CELLS['T1'])[0]
    pol=XueCorrected(kmax=300.,M=1.,dmin=50.,dmax=200.,ki=0.,KPF=0.,f0src='measured')
    pol.bind(p); lg=P.episode(p,pol,T=2.,seed=1000)
    hist=np.asarray(pol.hist)
    assert np.isfinite(lg).all()
    assert np.any(hist[:,2]>0), 'Stiffness branch still inactive'
    assert (hist[:,2]>=0).all() and (hist[:,2]<=300.+1e-8).all()
    # Zero learned residual exactly reproduces B1, including K, at same seeds.
    base=next(c for c in read('gate_PC125.json') if c['task']=='T1' and c['dose']==CELLS['T1'])['B1']['gains']
    z=init_theta(0)
    l0=P.episode(p,tuple(base),T=2.,seed=1000)
    l1=P.episode(p,LearnedCorrected(z,base),T=2.,seed=1000)
    error=float(np.max(abs(l0[:,3]-l1[:,3])))
    z[104:128]=np.random.default_rng(31).normal(0,.05,24)
    py=P.episode(p,LearnedCorrected(z,base),T=2.,seed=1000)
    fast=fast_learned_episode(p,LearnedCorrected(z,base),2.,1000)
    fast_error=float(np.max(abs(py-fast)))
    assert fast_error<1e-10,fast_error
    assert error<1e-9,error
    jm0=P.episode(p,B.JinMPC(mode='ms'),T=1.,seed=1000,update=.016)
    jm1=P.episode(p,JinCorrected(mode='ms'),T=1.,seed=1000,update=.016)
    ms_error=float(np.max(abs(jm0-jm1)))
    assert ms_error<1e-10,ms_error
    # The old 8 ms OS endpoint is action-blind through the queue. The
    # corrected endpoint must respond to a finite admissible gain move.
    st=np.array([.002,0.,0.,.002,.002,3.,3.,10000.,0.])
    q=np.array([.002]); pr=np.array([1.,0.,10.,.002,.65,.78,.008,.25,1.,0.,.01,6.5,50.,.5])
    response=[]
    for nt in (2,3):
        a=np.zeros(3); bb=np.zeros(3)
        B._predict(nt,1,np.array([math.log(20.)]),np.array([0.]),st,q,pr,True,a)
        B._predict(nt,1,np.array([math.log(2000.)]),np.array([.001]),st,q,pr,True,bb)
        response.append(float(bb[0]-a[0]))
    assert abs(response[0])<1e-12 and abs(response[1])>1e-3,response
    out=dict(xue_stiffness_min=float(hist[:,2].min()),xue_stiffness_max=float(hist[:,2].max()),
             learned_zero_residual_max_force_difference_N=error,compiled_learned_max_log_difference=fast_error,
             xue_replica_error=pol.maxrep,jin_ms_log_difference=ms_error,jin_old_endpoint_action_response=response[0],jin_corrected_endpoint_action_response=response[1],passed=True)
    save('checks.json',out); print('CHECKS',out,flush=True)

def main():
    ap=argparse.ArgumentParser(__doc__); ap.add_argument('mode',choices=['checks','tune','polish','learn','calibrate','test','all']); ap.add_argument('--generations',type=int,default=100); ap.add_argument('--workers',type=int,default=6)
    args=ap.parse_args(); OUT.mkdir(exist_ok=True)
    if args.mode in ('checks','all'): checks()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        if args.mode in ('tune','all'): tune(pool)
        if args.mode in ('polish','all'): polish(pool)
        if args.mode in ('learn','all'): train_learning(pool,args.generations)
        if args.mode in ('calibrate','all'): calibrate_learning(pool)
        if args.mode in ('test','all'): final_test(pool)

if __name__=='__main__': main()
