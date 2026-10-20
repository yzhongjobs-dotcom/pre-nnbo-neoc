"""Numerical/causal checks that matter for interpreting the experiments."""
import json
import numpy as np
from core import m,plant,local_contact,cost_gradient,rollout,encode
from frontend import Frontend
from run import CFG,ROOT


def main():
    p=m.sample(1,9011)[0];p[[1,4,7,11]]=0.;p[13]=3.5
    y=plant.initial(p,CFG[:4],257);q=np.array([np.log(p[0]),0,0,0,0,np.log(p[12]),p[10]])
    f=plant.contact(.1,y,p,257)[0];fl=local_contact(.1,y,q,np.zeros(2),p[9],257)[0]
    np.testing.assert_allclose(f,fl,rtol=1e-12,atol=1e-12)
    a=np.tile(encode(np.r_[CFG[:4],.2]),2);qc=np.r_[q,f];anchor=np.zeros(2)
    y[1]+=.0005;y[3]=.001
    value,g=cost_gradient(.1,y,y[1]-.0005,qc,anchor,p[9],CFG,a,p,True)
    direction=np.random.default_rng(2).normal(size=10);direction/=np.linalg.norm(direction)
    eps=1e-4
    fd=(rollout(.1,y,y[1]-.0005,qc,anchor,p[9],CFG,a+eps*direction,p,True)-rollout(.1,y,y[1]-.0005,qc,anchor,p[9],CFG,a-eps*direction,p,True))/(2*eps)
    rel=abs(fd-g@direction)/max(.01,abs(fd));assert rel<.03,(fd,g@direction,rel)
    # Future/unreleased wrench changes cannot affect the current estimate.
    fa=Frontend(y[1],p[9]);fb=Frontend(y[1],p[9])
    for i in range(11):
        fa.observe(i*.01,y,f);fb.observe(i*.01,y,f+(100 if i>=9 else 0))
    qa,_=fa.update(10);qb,_=fb.update(10);np.testing.assert_array_equal(qa,qb)
    r={"uniform_field_force_error":float(np.max(abs(f-fl))),"directional_gradient_relative_error":rel,"gradient_cost":value,"causal_latest_used_index":8,"current_index":10}
    (ROOT/"checks.json").write_text(json.dumps(r,indent=2));print(r)


if __name__=="__main__":main()
