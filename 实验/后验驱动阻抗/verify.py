"""Independent ODE integration of physical means, covariance and running cost."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
import mechanism as m


def independent(k, b, K, D, duration, r0, speed, sigma, x0):
    # x0 is physical here. No spectral formula or analytical kernel reused.
    def rhs(t, state):
        q, v = state[:2]; cov = state[2:6].reshape(2, 2)
        r = r0+speed*t
        u = -K*(q-r)-D*(v-speed)
        A = np.array([[0., 1.], [-(k+K)/.5, -(b+D)/.5]])
        W = np.diag([0., (sigma/.5)**2])
        dc = A@cov+cov@A.T+W
        fv = np.array([k, b]); uv = np.array([-K, -D])
        cost = ((k*q+b*v)**2+fv@cov@fv)/.3**2
        cost += ((q-r)**2+cov[0, 0])/.0015**2
        cost += .04*((v-speed)**2+cov[1, 1])/.03**2
        cost += .015*(u*u+uv@cov@uv)/9
        return np.r_[v, u/.5-(k*q+b*v)/.5, dc.ravel(), cost]
    initial = np.r_[x0, np.zeros(5)]
    sol = solve_ivp(rhs, [0, duration], initial, method='DOP853', rtol=1e-11, atol=1e-14)
    assert sol.success
    return sol.y[:, -1]


def main():
    rows=[]; rng=np.random.default_rng(11401)
    for i in range(12):
        k, b=m.MODES[i%4]; K, D=m.GAINS[(i*7)%len(m.GAINS)]
        h=.04 if i%2==0 else .6
        r0=.0002; speed=.006; sigma=.03
        x0=np.array([rng.normal(0,.0001),rng.normal(0,.001)])
        S,c,P,beta,cov=m.kernel(k,b,K,D,h,r0,speed,sigma)
        actual=independent(k,b,K,D,h,r0,speed,sigma,x0)
        pred=(P@(x0/m.SCALE)+beta)*m.SCALE
        y=np.r_[x0/m.SCALE,1.]
        J=y@S@y+c
        row=dict(mean_abs=float(np.max(np.abs(pred-actual[:2]))),
                 covariance_abs=float(np.max(np.abs(cov*np.outer(m.SCALE,m.SCALE)-actual[2:6].reshape(2,2)))),
                 J_relative=float(abs(J-actual[-1])/max(abs(J),1e-8)))
        rows.append(row)
    assert max(r['J_relative'] for r in rows)<1e-7,rows
    ctx=m.sample(1,11402)[0];ctx['prior']=[1.,0.,0.,0.]
    S,c,P,beta,cov=m.prepare(ctx);x=np.array(ctx['x0']);pr=np.array(ctx['prior'])
    means,chol,inv,logdet=m.predictions(x,P,beta,cov);nodes,w=m.sobol_nodes(8)
    ss=[m.scores(x,pr,S,c,means,chol,inv,logdet,ctx['force_noise'],nodes,w,i) for i in range(3)]
    known=max(float(np.max(np.abs(ss[i]-ss[0]))) for i in [1,2])
    assert known<1e-10,known
    ctx['prior']=[.25]*4;pr=np.array(ctx['prior'])
    ss=[m.scores(x,pr,S,c,means,chol,inv,logdet,1e9,nodes,w,i) for i in [1,2]]
    uninformative=float(np.max(np.abs(ss[0]-ss[1])))
    assert uninformative<1e-8,uninformative
    result=dict(independent_ODE=rows,known_mode_max_abs=known,
                uninformative_force_max_abs=uninformative,
                scope='Numerical and information-use checks, not physical validation or an optimality certificate')
    Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__': main()
