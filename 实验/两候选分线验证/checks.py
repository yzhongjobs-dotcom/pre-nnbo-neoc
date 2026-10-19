"""Independent derivative, execution and adapted-reference checks."""
import json
import numpy as np
from scipy.optimize import brentq
from optimization import Problem,context
from core import m,plant,local_contact,rollout,encode
from baselines import NNBO,calibrate
from experiment import ROOT,BASE
from rate_core import rollout_rate,advance_rate,target_for_rate


def main():
    cfg=np.array(json.loads((BASE/'tuned_fixed.json').read_text())['cfg']);p=m.sample(1,9011)[0];cfg[7]=p[9]
    y=plant.initial(p,cfg[:4],65);zr=y[1];q=np.array([np.log(p[0]),0,0,0,0,np.log(p[12]),p[10]])
    anchor=np.zeros(2);sensed=local_contact(0.,y,q,anchor,p[9])[0];a=np.tile(np.array([.3,-.3,.3,-.3,.1]),2)
    pr=Problem();records=[]
    for active in (False,True):
        yy=y.copy()
        if active:yy[0]=.004;yy[1]+=.008;yy[8]=m.E_MIN+1e-6
        c=context(0,yy,zr,q,anchor,sensed,cfg)
        value,g=pr.vg(a,c);g=np.array(g).ravel();v2=rollout(0.,yy,zr,np.r_[q,sensed],anchor,p[9],cfg,a,p,False,dt=.004)
        assert abs(float(value)-v2)<1e-6,(value,v2)
        direction=np.random.default_rng(27).normal(size=10);direction/=np.linalg.norm(direction);eps=1e-5
        fd=(float(pr.value(a+eps*direction,c))-float(pr.value(a-eps*direction,c)))/(2*eps)
        err=abs(fd-g@direction)/max(.01,abs(fd));assert err<.01,(err,fd,g@direction)
        records.append({'energy_case':active,'cost_implementation_difference':float(value)-v2,'directional_relative_error':err})
    hc,cal=calibrate(p,zr,1907);nn=NNBO(hc,zr,cfg,omega=2.,learning=0.)
    nn.decide(0,y,np.array([3.01,.8,0.]))
    def equation(z):return .01+(2*nn.w[0]*(z-zr)+nn.w[1]*zr)/(2*nn.b*(z-nn.surface)**nn.n)
    root=brentq(equation,max(1e-7,.5*zr),zr)
    difference=abs(root-nn.reference);assert difference<1e-7,difference
    rate_problem=Problem(rate_actions=True,gains_only=True);rate_records=[]
    ar=np.tile(np.array([.1,.3,-.1,.15,0.]),2)
    for active in (False,True):
        yy=y.copy()
        if active:yy[0]=.004;yy[1]+=.008;yy[8]=m.E_MIN+1e-6
        c=context(0,yy,zr,q,anchor,sensed,cfg);v,g=rate_problem.vg(ar,c)
        v2=rollout_rate(0.,yy,zr,np.r_[q,sensed],anchor,p[9],cfg,ar,p,False,dt=.002)
        assert abs(float(v)-v2)<1e-6,(v,v2)
        direction=np.random.default_rng(14).normal(size=10);direction[[4,9]]=0;direction/=np.linalg.norm(direction);eps=1e-5
        fd=(float(rate_problem.value(ar+eps*direction,c))-float(rate_problem.value(ar-eps*direction,c)))/(2*eps)
        error=abs(fd-np.array(g).ravel()@direction)/max(.01,abs(fd));assert error<.01,error
        target=target_for_rate(yy,m.GRATE*ar[:4]);assert np.all(target>=m.GMIN) and np.all(target<=m.GMAX)
        after,_=advance_rate(0.,yy,p,m.GRATE*ar[:4],0.,zr,sensed,cfg,.1)
        assert after[8]>=m.E_MIN-1e-8;assert np.all(after[4:8]>=m.GMIN-1e-6) and np.all(after[4:8]<=m.GMAX+1e-6)
        rate_records.append({'energy_case':active,'cost_difference':float(v)-v2,'directional_relative_error':error,'energy_margin':after[8]-m.E_MIN})
    out={'derivatives':records,'rate_derivatives_and_execution':rate_records,'nnbo_reference_independent_root_difference':difference,'calibration':cal}
    (ROOT/'checks.json').write_text(json.dumps(out,indent=2));print(out)


if __name__=='__main__':main()
