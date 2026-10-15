"""Minimal source-equation NNBO reproduction, paper III-A--C and IV.

Integrates the force-as-input augmented system (24), policy (37), and critic
law (41). This is not a reproduction of R15 or of the robot position loop.
Unreported uncertainty bound/noise duration are explicit, not fitted to plots.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
PDF=HERE.parents[1]/'文献'/'原文'/'NNBO.pdf'
CFG=dict(Ke=1149.,Be=567.,n=1.2,xe=.01,x0=.011,r0=.08,
         command_A=-1.,command_alpha=.03,command_output=1.,R=1.,
         initial_weights=[20.,20.,20.],learning_rate=5.,sample_dt=.001,
         end_time=7.,cost_Q=[500000.,1000000.],
         basis=['e^2','e*r','r^2'],uncertainty='epsilon=1000*e*sin(r)',
         uncertainty_bounds=['tight_known_disturbance','uniform_sine_bound'],
         exploration_duration='entire 7 s; end time not specified in paper',
         integration='continuous ODE, RK45 adaptive internal steps, 1 ms output grid',
         main_rtol=1e-8,main_atol=1e-10,max_internal_step=.001)

def command(t):
    r=.03+.05*np.exp(-t)
    return r,.03-r

def exploration(t):
    return (np.cos(2.4*t)*np.sin(2.4*t)**3+np.sin(1.12*t)**2+
            np.sin(t)**5+np.sin(t)**2*np.cos(t)+
            np.sin(1.2*t)**2*np.cos(.5*t)+np.sin(2*t)**2*np.cos(.1*t))

def quantities(t,y,Q,bound):
    x,w1,w2,w3=y
    r,rd=command(t); e=x-r; delta=x-CFG['xe']
    if delta<=0:
        raise FloatingPointError(f'HC domain violation: t={t}, delta={delta}')
    g=1/(CFG['Be']*delta**CFG['n'])
    grad_e=2*w1*e+w2*r
    fhat=-.5*g*grad_e/CFG['R']
    eps=1000*e*np.sin(r)
    zeta=abs(eps) if bound=='tight_known_disturbance' else 1000*abs(e)
    ed_aux=-CFG['Ke']/CFG['Be']-rd+g*fhat
    nu=np.array([2*e*ed_aux,r*ed_aux+e*rd,2*r*rd])
    stage=Q*e*e+CFG['R']*fhat*fhat
    hamiltonian=np.dot(y[1:],nu)+stage+zeta*zeta
    wdot=-CFG['learning_rate']*nu*hamiltonian/(1+np.dot(nu,nu))**2
    applied=fhat+exploration(t)
    xd=-CFG['Ke']/CFG['Be']+g*applied+eps
    return dict(r=r,e=e,delta=delta,g=g,fhat=fhat,applied=applied,
                epsilon=eps,zeta=zeta,nu=nu,hamiltonian=hamiltonian,
                stage=stage,wdot=wdot,xd=xd,
                value=w1*e*e+w2*e*r+w3*r*r)

def rhs(Q,bound):
    def function(t,y):
        q=quantities(t,y,Q,bound)
        return np.r_[q['xd'],q['wdot']]
    return function

def integrate(Q,bound,method='RK45',rtol=None,atol=None):
    ts=np.arange(7001,dtype=float)*CFG['sample_dt']
    return solve_ivp(rhs(Q,bound),(0.,CFG['end_time']),
                     [CFG['x0'],*CFG['initial_weights']],t_eval=ts,
                     max_step=CFG['max_internal_step'],method=method,
                     rtol=rtol or CFG['main_rtol'],atol=atol or CFG['main_atol'])

def newton_reference(target_force,r,w,start):
    """Offline check of Eq.52/Algorithm1, not connected to ideal position plant.

    Holds measured target force, r and W fixed, uses analytic derivative and
    downhill merit reduction. Domain is respected, never clip force/weights.
    """
    def equation(x):
        delta=x-CFG['xe']
        if delta<=0:return np.inf,np.nan
        g=1/(CFG['Be']*delta**CFG['n'])
        grad=2*w[0]*(x-r)+w[1]*r
        return target_force+.5*g*grad, w[0]*g-.5*CFG['n']*g*grad/delta
    x=start
    for i in range(60):
        v,dv=equation(x)
        if abs(v)<1e-10:return x,i,abs(v)
        step=-v/dv
        for j in range(50):
            trial=x+step*2.**(-j)
            tv,_=equation(trial)
            if np.isfinite(tv) and abs(tv)<abs(v):
                x=trial;break
        else:raise RuntimeError('Reference downhill search failed')
    raise RuntimeError('Reference inversion did not converge')

def checks_and_save(sol,Q,bound):
    assert sol.success and len(sol.t)==7001, sol.message
    assert np.all(np.isfinite(sol.y))
    qs=[quantities(t,y,Q,bound) for t,y in zip(sol.t,sol.y.T)]
    arr=lambda key:np.array([q[key] for q in qs])
    xs=sol.y[0]; stages=arr('stage'); hc_nom=(CFG['Be']*arr('xd')+CFG['Ke'])*arr('delta')**CFG['n']
    # From Eq.23 epsilon=-Delta/(Be*delta^n); match Eq.24 with additive Delta.
    force_unc=-CFG['Be']*arr('delta')**CFG['n']*arr('epsilon')
    hc_measured=hc_nom+force_unc
    identity_error=float(np.max(abs(hc_measured-arr('applied'))))
    assert identity_error<1e-10
    assert np.all(arr('zeta')+1e-12>=abs(arr('epsilon')))
    assert np.min(xs)>CFG['xe'] and CFG['xe']<xs[-1]<arr('r')[-1]
    assert np.min(arr('applied'))>=0
    value=arr('value'); integral=cumulative_trapezoid(stages,sol.t,initial=0.)
    xi=np.array([q['nu']/(1+q['nu']@q['nu']) for q in qs])
    gram=np.trapz(xi[:,:,None]*xi[:,None,:],sol.t,axis=0)
    # Algebraic inversion consistency of recorded force, paper Eq.52.
    inverse=[]
    for i in range(0,len(sol.t),10):
        q=qs[i]; target=q['applied']; r=q['r']; w=sol.y[1:,i]
        xr,it,res=newton_reference(target,r,w,xs[i])
        def eq(x):return target+(2*w[0]*(x-r)+w[1]*r)/(2*CFG['Be']*(x-CFG['xe'])**CFG['n'])
        hi=max(r+1.,xs[i]+1.)
        oracle=brentq(eq,CFG['xe']+1e-10,hi,xtol=1e-14)
        assert abs(xr-oracle)<1e-8
        inverse.append((sol.t[i],xs[i],xr,target,it,res))
    stem=f'Q{int(Q)}_{bound}'
    cols=['time_s','x_m','desired_m','error_m','policy_force_N','applied_force_N',
          'exploration_N','epsilon_m_s','zeta_m_s','W1','W2','W3',
          'stage_cost','finite_horizon_cost_integral','value_estimate','hamiltonian_residual']
    data=np.column_stack([sol.t,xs,arr('r'),arr('e'),arr('fhat'),arr('applied'),
                          exploration(sol.t),arr('epsilon'),arr('zeta'),sol.y[1:].T,
                          stages,integral,value,arr('hamiltonian')])
    with (HERE/f'{stem}.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(cols);writer.writerows(data)
    with (HERE/f'{stem}_reference_check.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(['time_s','current_x_m','inverted_reference_m','fixed_force_N','iterations','residual_N']);writer.writerows(inverse)
    np.savez(HERE/f'{stem}.npz',t=sol.t,y=sol.y,data=data)
    tail=sol.t>=6.
    return dict(case=stem,Q=Q,bound=bound,status='PASS_MINIMAL_MECHANISM',
                output_samples=len(sol.t),rhs_evaluations=sol.nfev,
                minimum_penetration_m=float(min(xs-CFG['xe'])),
                final_position_m=float(xs[-1]),final_desired_m=float(arr('r')[-1]),
                final_weights=sol.y[1:,-1].tolist(),
                final_policy_force_N=float(arr('fhat')[-1]),
                final_applied_force_N=float(arr('applied')[-1]),
                peak_applied_force_N=float(max(arr('applied'))),
                tail_position_minmax_m=[float(min(xs[tail])),float(max(xs[tail]))],
                tail_max_abs_weight_derivative=float(max(np.linalg.norm(arr('wdot')[tail],axis=1))),
                tail_mean_stage_cost=float(np.mean(stages[tail])),
                finite_7s_cost=float(integral[-1]),
                final_hamiltonian_residual=float(arr('hamiltonian')[-1]),
                observed_finite_interval_PE_gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                hc_force_identity_max_abs_error=identity_error,
                reference_inversion_samples=len(inverse),
                reference_inversion_max_residual_N=float(max(row[-1] for row in inverse)),
                reference_inversion_max_iterations=int(max(row[-2] for row in inverse)))

def independent_formula_check():
    # Independent matrix expressions for Eq.37 and 41; not a copy of scalar rhs.
    t=.4;y=np.array([.023,100.,-70.,40.]);Q=500000.
    q=quantities(t,y,Q,'tight_known_disturbance')
    X=np.array([q['e'],q['r']]);dphi=np.array([[2*X[0],0],[X[1],X[0]],[0,2*X[1]]])
    F=np.array([-CFG['Ke']/CFG['Be']-(-q['r']+.03),-q['r']+.03])
    G=np.array([q['g'],0.]);f=-.5*(G@(dphi.T@y[1:]))
    nu=dphi@(F+G*f)
    H=y[1:]@nu+Q*X[0]**2+f*f+q['zeta']**2
    wd=-CFG['learning_rate']*nu*H/(1+nu@nu)**2
    assert abs(f-q['fhat'])<1e-12
    assert np.max(abs(wd-q['wdot']))<1e-10
    return dict(matrix_scalar_policy_agree=True,matrix_scalar_learning_agree=True)

def main():
    start=time.perf_counter()
    (HERE/'config.json').write_text(json.dumps(CFG,indent=2),encoding='utf-8')
    result=dict(status='RUNNING',config=CFG,python=platform.python_version(),
                numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
                source_pdf=str(PDF),source_pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(),
                source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                formula_checks=independent_formula_check(),cases=[])
    for bound in CFG['uncertainty_bounds']:
        for Q in CFG['cost_Q']:
            print(f'Running Q={Q:.0f}, bound={bound}',flush=True)
            sol=integrate(Q,bound)
            metrics=checks_and_save(sol,Q,bound)
            result['cases'].append(metrics)
            print(json.dumps(metrics),flush=True)
    print('Checking first case with independent implicit Radau integration...',flush=True)
    base=np.load(HERE/'Q500000_tight_known_disturbance.npz')
    independent=integrate(500000.,'tight_known_disturbance',method='Radau',rtol=1e-9,atol=1e-11)
    assert independent.success
    dx=float(np.max(abs(base['y'][0]-independent.y[0])))
    dw=float(np.max(abs(base['y'][1:]-independent.y[1:])))
    assert dx<1e-6 and dw<.01
    result['integration_check']=dict(independent_method='Radau',max_position_difference_m=dx,
                                     max_weight_difference=dw,passed=True)
    result['duration_seconds']=time.perf_counter()-start
    result['status']='PASS_MINIMAL_MECHANISM_NOT_EXACT_FIGURE_REPLICATION'
    (HERE/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=result['status'],duration=result['duration_seconds'],integration_check=result['integration_check'])),flush=True)

if __name__=='__main__':main()
