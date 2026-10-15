"""Formula-level audit of the strongest prepared value-root chord baseline.

No critic root is solved here. This checks a structural identity for arbitrary
fixed critic weights, so it must not be reported as control/HJB validation.
"""
import json,numpy as np
from contact_stages import Problem,ROOT,environment,dump


def basis(p,t,x):
    center=np.array([0.,0.,200.,6.,.00105]);scale=np.array([.003,.03,500.,15.,.0001])
    coord=(x-center)/scale
    phi=np.r_[1.,t/4.,coord];factor=1-t/4.
    psi=factor*phi;pt=-phi/4.;pt[1]+=factor/4.
    jac=np.zeros((5,7));jac[:,2:]=factor*np.diag(1/scale)
    return psi,pt,jac


def residual(p,t,x,W,theta):
    _,pt,jac=basis(p,t,x);c=p.cfg;e,v,K,D,E=x
    _,gamma=environment(t,p.p,p.pulses)
    qterm=c['terminal']*(1/c['error_scale']**2+c['w_force']*(theta/c['force_scale'])**2)
    grad=jac@W;grad[0]+=2*qterm*e
    linear=np.array([p.rate[0]*(grad[2]-.5*e*e*grad[4]),p.rate[1]*grad[3]])
    lo=np.maximum(-1.,-20*(x[2:4]-p.lo)/p.rate)
    hi=np.minimum(1.,20*(p.hi-x[2:4])/p.rate)
    if e*e>1e-25:
        hi[0]=min(hi[0],(c['recharge']*(c['robot_damping']+D)*v*v+
                        c['energy_alpha']*(E-c['energy_min']))/(.5*e*e*p.rate[0]))
    u=np.clip(-linear/(2*c['w_rate']),lo,hi)
    f0=np.array([v,(-(theta+K)*e-(c['robot_damping']+D)*v+theta*gamma)/c['mass'],
                      0.,0.,c['recharge']*(c['robot_damping']+D)*v*v])
    f=f0+np.array([0.,0.,u[0]*p.rate[0],u[1]*p.rate[1],-.5*e*e*u[0]*p.rate[0]])
    command=-K*e-D*v
    ell=c['w_position']*(e/c['error_scale'])**2+c['w_velocity']*(v/c['velocity_scale'])**2
    ell+=c['w_force']*(theta*(gamma-e)/c['force_scale'])**2+c['w_command']*(command/c['command_scale'])**2
    ell+=c['w_compliance']*np.sum(((x[2:4]-p.pref)/np.array(c['gain_scale']))**2)
    B=pt@W+ell+c['w_rate']*(u@u)+grad@f
    J=pt+jac.T@f  # Envelope derivative: no need to differentiate clipped actor.
    A1=(jac@W)[1]*(gamma-e)/c['mass']
    A2=c['w_force']*((gamma-e)/c['force_scale'])**2+2*c['terminal']*c['w_force']*e*v/c['force_scale']**2
    return B,J,u,A1,A2


def main():
    p=Problem();rng=np.random.default_rng(31);nodes=[]
    for _ in range(7):
        t=rng.uniform(.1,3.7);_,gamma=environment(t,p.p,p.pulses)
        x=np.array([gamma+rng.uniform(-.0003,.0003),rng.uniform(-.005,.005),
             rng.uniform(40.,400.),rng.uniform(2.,10.),.001+rng.uniform(.000003,.0001)])
        nodes.append((t,x))
    W=np.array([.1,.3,-.02,.003,-.0006,.0008,1e-7]);theta_a=400.
    initial=[residual(p,t,x,W,theta_a) for t,x in nodes]
    Fa=np.array([r[0] for r in initial]);Ja=np.array([r[1] for r in initial]);A1=np.array([r[3] for r in initial]);A2=np.array([r[4] for r in initial])
    C1=np.linalg.solve(Ja,-(A1+2*theta_a*A2));C2=np.linalg.solve(Ja,-A2)
    correction=np.linalg.solve(Ja,-Fa);records=[]
    for delta in [-50.,-5.,1.,5.,50.]:
        new=[residual(p,t,x,W,theta_a+delta) for t,x in nodes]
        Fb=np.array([r[0] for r in new]);Jb=np.array([r[1] for r in new])
        polynomial=Fa+(A1+2*theta_a*A2)*delta+A2*delta*delta
        chord=W+np.linalg.solve(Ja,-Fb)
        prepared=W+correction+C1*delta+C2*delta*delta
        records.append(dict(delta=delta,residual_identity_max_abs=float(np.max(abs(Fb-polynomial))),
          actor_difference_max_abs=float(np.max(abs(np.array([r[2] for r in new])-np.array([r[2] for r in initial])))),
          prepared_chord_difference_max_abs=float(np.max(abs(chord-prepared))),
          first_order_only_omitted_residual_norm=float(np.linalg.norm(A2*delta*delta))))
    report=dict(scope='Structural identity at arbitrary fixed critic weights; NOT a solved/validated critic, root accuracy, strategy result or novel method.',
       exact_structure='F(Wa,theta)=A0(Wa)+A1(Wa)*theta+A2*theta^2. Actor and permitted action set are independent of material theta at fixed state and weights.',
       prepared_candidate='Wa - Ja^-1*Fa + C1*delta + C2*delta^2; C1=-Ja^-1*Ftheta; C2=-Ja^-1*A2. For an accepted root anchor, the first correction is zero.',
       baseline_implication='A fair cached-J chord can prepare two response columns and query at the same linear-in-weight-count scale as NEOC. Generic speed claims must confront this structure.',
       limitations='Applies to fixed basis/reference and material-independent action set, not arbitrary eta, varying contact boundaries or state-constrained HJB solutions.',
       Jacobian_condition_number=float(np.linalg.cond(Ja)),anchor_residual_norm=float(np.linalg.norm(Fa)),records=records)
    dump(ROOT/'prepared_residual_structure.json',report);print(json.dumps(report,indent=2),flush=True)
    assert max(r['residual_identity_max_abs'] for r in records)<1e-10
    assert max(r['actor_difference_max_abs'] for r in records)==0
    assert max(r['prepared_chord_difference_max_abs'] for r in records)<1e-10


if __name__=='__main__':main()
