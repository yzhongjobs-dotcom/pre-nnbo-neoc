"""Independent DOP853 plant with feedback re-evaluated, plus mesh refinement.

Controllers are intentionally shared; plant pressure, actuator equation and costs
are independently written. This checks numerics, not the truth of the tissue model.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
import model as m
import reference_controller as rc
from verify import pressure_discrete,gain_independent
from run import noises
HERE=Path(__file__).resolve().parent


def independent_episode(p,noise,method,cfg0,policy,geometry=True):
    cfg=cfg0.copy();cfg[7]=p[9]
    y=m.initial(p,cfg[:4]);zr=y[1]
    initial_k=1000.
    if geometry:
        initial_k=1.5*m.FD/zr if zr<m.A*m.A/(2*m.R) else m.FD/(zr-m.A*m.A/(6*m.R))
        initial_k=min(7000.,max(200.,initial_k))
    old=np.array([initial_k,0.,0.,3.,.6,zr,0.,0.])
    sensed=np.array([3.,.6,0.]);hist=np.zeros((len(noise),8));metrics=np.zeros(4)
    target=cfg[:4].copy();ff=0.;peak=0.
    for i in range(len(noise)-1):
        t=i*.01
        force=pressure_discrete(t,y,p,33)
        hist[i,:4]=y[:4];hist[i,4:7]=force+noise[i];hist[i,7]=p[9]*t+y[0]
        resp,temporal,j=m.response(hist,i,old);old=resp.copy()
        ob=hist[max(0,i-2),4:7].copy()
        if method!=0:ob[:2]=resp[3:5]
        sensed=cfg[6]*ob+(1-cfg[6])*sensed;resp[3:5]=sensed[:2]
        if i%4==0:
            if method==13:choice=rc.decide(y,zr,cfg,resp,temporal)
            else:choice,_=m.decide(y,zr,cfg,resp,temporal,method,policy)
            target=choice[:4];ff=choice[4]
        def rhs(tt,yy):
            state=yy[:10]
            force=pressure_discrete(tt,state,p,33)
            rate,edot=gain_independent(state,target,zr)
            ef=3.-sensed[0]
            ut=sensed[1]-state[4]*state[0]-state[6]*state[2]+.5*p[9]
            un=3.+cfg[4]*ef+state[9]-state[5]*(state[1]-zr)-state[7]*state[3]+ff
            uta=np.clip(ut,-12.,12.);una=np.clip(un,-12.,12.)
            iy=cfg[5]*ef
            if (state[9]>=5 and iy>0) or (state[9]<=-5 and iy<0):iy=0.
            if (un>12 and iy>0) or (un< -12 and iy<0):iy=0.
            cost=((force[0]-3)/.3)**2+(state[0]/.0015)**2+cfg[8]*((state[1]-zr)/.0015)**2
            cost+=.04*(state[3]/.03)**2+.015*((uta/3)**2+((una-3)/3)**2)
            return np.r_[state[2],state[3],(uta-force[1]-.5*(state[2]+p[9]))/.8,
                (una-force[0]-.5*state[3])/.5,rate,edot,iy,cost,(force[0]-3)**2,state[0]**2,(state[1]-zr)**2]
        init=np.r_[y,np.zeros(4)]
        sol=solve_ivp(rhs,(t,t+.01),init,method='DOP853',rtol=1e-8,atol=1e-11,max_step=.002)
        y=sol.y[:10,-1];metrics+=sol.y[10:,-1]
        peak=max(peak,force[0])
    return metrics,y,peak


def main():
    cfg=np.array(json.loads((HERE/'joint_validation.json').read_text())['cfgs']['fixed_PI'])
    policy=np.zeros(8)
    path=HERE/'final_policies.json'
    if path.exists():policy=np.array(json.loads(path.read_text())['policies']['learned_feedback'])
    ps=m.sample(18,8201);ns=noises(ps,8202)
    output=[]
    for case in [0,2]:
      for method in [0,6,12,13]:
        cc=cfg.copy()
        if method==13:
            cc=np.array(json.loads((HERE/'reference_tuning.json').read_text())['cfg'])
            fast=rc.episode(ps[case],ns[case],cc)[0]
            fine=rc.episode(ps[case],ns[case],cc,dt=.00025,nq=129)[0]
        else:
            fast=m.episode(ps[case],ns[case],method,cc,policy=policy,geometry_prior=True)[0]
            fine=m.episode(ps[case],ns[case],method,cc,policy=policy,geometry_prior=True,dt=.00025,nq=129)[0]
        independent,_,_=independent_episode(ps[case],ns[case],method,cc,policy)
        row={'validation_case':case,'method':method,'baseline_J':float(fast[0]),'fine_J':float(fine[0]),
           'DOP853_J':float(independent[0]),'relative_mesh_difference':float(fine[0]/fast[0]-1),
           'relative_DOP853_difference':float(independent[0]/fast[0]-1)}
        output.append(row);print(row,flush=True)
    result={'scope':'Recomputed causal feedback on independently written DOP853 plant, no latency; validation cases, not test selection',
      'rows':output,'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['model.py','audit_closed_loop.py','verify.py','reference_controller.py']}}
    (HERE/'closed_loop_checks.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
