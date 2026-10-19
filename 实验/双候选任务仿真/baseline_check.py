"""Independent source-equation rollout and gradient checks for Jin MS VIC."""
import json
from pathlib import Path
import numpy as np
from baseline import JinMS
import physics as p

def independent(mu,state,solver):
    x=state.copy();total=0.;dt=solver.T
    for a in mu.reshape(-1,4):
        g=p.GMIN+(np.tanh(a)+1)/2*(p.GMAX-p.GMIN)
        for axis in range(2):
            e,v,f=x[3*axis:3*axis+3];mass=p.MASS[axis];acc=(f-g[axis]*e-g[axis+2]*v)/mass
            ne=e+dt*v+.5*dt*dt*acc;nv=v+dt*acc
            ell=(ne/(.002 if axis==0 else .001))**2+.04*(nv/.015)**2+.1*f*f
            ell+=solver.reg[axis]*((g[axis]/mass)**2+(g[axis+2]/mass)**2+(1/mass)**2)
            distance=np.sqrt(max(abs(ne)/.006-1,0)**2+max(abs(nv)/.06-1,0)**2+1e-20)
            total+=dt*(ell+50*(1-np.exp(-distance/2)));x[3*axis:3*axis+3]=[ne,nv,f]
    return total+.1*((x[0]/.002)**2+(x[3]/.001)**2)

if __name__=='__main__':
    rng=np.random.default_rng(16801);solver=JinMS();mu=rng.normal(0,.2,40);state=np.array([.0008,.001,-.1,-.0001,.001,.02])
    cost,g=solver.objective(mu,state);direction=rng.normal(size=40);h=1e-5
    fd=(independent(mu+h*direction,state,solver)-independent(mu-h*direction,state,solver))/(2*h)
    data={'source_equations':'E41 13/14/21/28/29/32/36--42 + Algorithm1',
          'independent_cost_difference':cost-independent(mu,state,solver),'directional_relative_error':abs(g@direction-fd)/max(abs(fd),1e-10),
          'scope':'method-equation reconstruction plus declared task adaptation; not original Panda experimental reproduction'}
    assert abs(data['independent_cost_difference'])<1e-8 and data['directional_relative_error']<1e-5
    print(json.dumps(data,indent=2));Path(__file__).with_name('baseline_check.json').write_text(json.dumps(data,indent=2))
