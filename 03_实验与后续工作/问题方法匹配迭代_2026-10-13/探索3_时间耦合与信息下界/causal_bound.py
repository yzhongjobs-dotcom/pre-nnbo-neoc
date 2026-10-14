"""Optimal information bound for gain realization with slew/tank constraints.

Convex gain-planning proxy. Two slowly changing material futures have identical
history until 2s; the strongest causal relaxation knows the complete future
immediately after reveal. Common pre-reveal actions encode nonanticipativity.
No claims about whole robot, imaging, HJB critic or NNBO/NEOC necessity.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import minimize,LinearConstraint

P=Path(__file__).resolve().parent
T=10.;N=80;H=T/N;TIMES=np.arange(N)*H
SCALE=1000.;EMIN=.001;RESERVE=.00015;EREF=.0005;ALPHA=10.;CHI=.3
KINIT=float(np.sqrt(2000.**2+1e6)-2000.);DUP=30.;SCAN=.006;BR=1.

def environment(t,soft):
    phase=np.clip((t-2)/6,0,1)
    return 2000.-(1900.*.5*(1-np.cos(np.pi*phase)) if soft else np.zeros_like(phase))

def motion(t):
    phase=np.clip((t-2.8)/.5,0,1)
    e=.00005+.00195*.5*(1-np.cos(np.pi*phase))
    edot=np.where((t>2.8)&(t<3.3),.00195*np.pi*np.sin(np.pi*phase),0.)
    return e,edot

def integrals():
    points,weights=np.polynomial.legendre.leggauss(8)
    average=[];power=[];first_half=[];power_half=[];midp=[];midd=[]
    for t in TIMES:
        s=t+H*(points+1)/2;e,edot=motion(s)
        average.append(float(np.sum(weights*.5*e*e)/2))
        power.append(float(H*np.sum(weights*CHI*(DUP*edot**2+BR*SCAN**2))/2))
        s=t+H*(points+1)/4;e,edot=motion(s)
        first_half.append(float(np.sum(weights*.5*e*e)/4))
        power_half.append(float(H*np.sum(weights*CHI*(DUP*edot**2+BR*SCAN**2))/4))
        e,edot=motion(t+H/2);midp.append(float(.5*e*e));midd.append(float(CHI*(DUP*edot**2+BR*SCAN**2)))
    return tuple(np.array(a) for a in [average,power,first_half,power_half,midp,midd])

AVG,POWER,HALF,PHALF,MIDP,MIDD=integrals()
L=np.tril(np.ones((N,N)));LPREV=L-np.eye(N)
CE=L*AVG[None,:]*SCALE
CEHALF=(LPREV*AVG[None,:]+np.diag(HALF))*SCALE
RESOURCE=RESERVE+np.cumsum(POWER)
RESOURCEHALF=RESERVE+LPREV@POWER+PHALF
TARGET=np.array([np.sqrt(environment(TIMES+H/2,s)**2+1e6)-environment(TIMES+H/2,s) for s in [False,True]])

def solve(rate,lead,full=False,build_only=False):
    # Full means independent plans: branch known from t=0. Causal relaxation
    # permits complete future knowledge after reveal-minus-lead.
    reveal=0. if full else max(0.,2.-lead)
    prefix=int(np.ceil(reveal/H-1e-12))
    cols=2*N-prefix;MAP=np.zeros((2*N,cols))
    for j in range(N):MAP[j,j]=1
    for j in range(N):MAP[N+j,j if j<prefix else N+j-prefix]=1
    LG=np.zeros((2*N,2*N));LG[:N,:N]=L;LG[N:,N:]=L
    G=LG@MAP
    C=np.zeros((2*N,2*N));C[:N,:N]=CE;C[N:,N:]=CE
    CP=np.zeros((2*N,2*N))
    one=ALPHA*CEHALF+np.diag(MIDP*SCALE/H)
    CP[:N,:N]=one;CP[N:,N:]=one
    matrix=np.vstack((G,C@MAP/EREF,CP@MAP/(ALPHA*EREF)))
    lower=np.r_[np.full(2*N,(5-KINIT)/SCALE),np.full(4*N,-np.inf)]
    upper=np.r_[np.full(2*N,(1500-KINIT)/SCALE),np.tile(RESOURCE/EREF,2),
                  np.tile((MIDD+ALPHA*RESOURCEHALF)/(ALPHA*EREF),2)]
    desired=(TARGET.reshape(-1)-KINIT)/SCALE
    Q=.5*H*(G.T@G);grad=-.5*H*G.T@desired
    # Fixed physical effort regularizer, independent of branch and information.
    Q+=.5*H*.0001*(SCALE/(H*rate))**2*(MAP.T@MAP)
    def obj(x):return float(x@Q@x+2*grad@x),2*Q@x+2*grad
    bound=H*rate/SCALE
    if build_only:
        return dict(Q=Q,gradient=grad,matrix=matrix,lower=lower,upper=upper,
                    variable_bound=bound,MAP=MAP,G=G,C=C,prefix=prefix,target=desired)
    tic=time.perf_counter()
    res=minimize(obj,np.zeros(cols),jac=True,method='SLSQP',bounds=[(-bound,bound)]*cols,
                 constraints=[LinearConstraint(matrix,lower,upper)],
                 options=dict(maxiter=250,ftol=1e-10))
    rates=MAP@res.x*SCALE/H;gains=(KINIT+G@res.x*SCALE).reshape(2,N)
    energy=(EMIN+np.tile(RESOURCE,2)-C@MAP@res.x).reshape(2,N)
    actual_cost=.5*H*np.sum(((gains-TARGET)/SCALE)**2)+.5*H*.0001*np.sum((rates/rate)**2)
    violation=float(max(np.max(matrix@res.x-upper),np.max(lower-matrix@res.x),0))
    out=dict(rate=rate,lead=lead,full_information=full,reveal_time=reveal,shared_steps=prefix,
        success=bool(res.success),message=str(res.message),proxy_cost=float(actual_cost),
        gain_rmse_N_m=float(np.sqrt(np.mean((gains-TARGET)**2))),
        soft_final_gain=float(gains[1,-1]),hard_final_gain=float(gains[0,-1]),
        min_energy=float(energy.min()),constraint_max_violation=violation,
        seconds=time.perf_counter()-tic,optimizer_iterations=int(res.nit))
    np.savetxt(P/(f'rate{rate}_lead'+('full' if full else str(lead))+'.csv'),
        np.column_stack((TIMES+H,TARGET[0],gains[0],energy[0],TARGET[1],gains[1],energy[1])),
        delimiter=',',header='t,hard_target,hard_gain,hard_energy,soft_target,soft_gain,soft_energy',comments='')
    assert res.success and violation<1e-7
    return out

def main():
    config=dict(T=T,N=N,h=H,kinit=KINIT,energy_min=EMIN,reserve=RESERVE,
        rates=[150,200,250],lead_times=[0,.25,.5,.8,1.],uncertainty='equiprobable hard/soft futures, identical before 2s',
        motion='planned relative error grows 0.05mm to2mm between2.8and3.3s',
        scope='gain-realization convex proxy, not full robot or clinical calibration',
        dissipation='causal and preview methods both receive maximal D=30 recharge allowance')
    (P/'frozen_config.json').write_text(json.dumps(config,indent=2),encoding='utf8')
    result={}
    for rate in config['rates']:
        result[str(rate)]=dict(full=solve(rate,0,full=True),leads=[])
        for lead in config['lead_times']:
            r=solve(rate,lead);r['extra_cost_vs_full']=r['proxy_cost']/result[str(rate)]['full']['proxy_cost']-1
            result[str(rate)]['leads'].append(r)
            (P/'results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
            print(rate,lead,r['proxy_cost'],r['extra_cost_vs_full'],r['gain_rmse_N_m'],flush=True)
    print('causal information bound complete',flush=True)

if __name__=='__main__':main()
