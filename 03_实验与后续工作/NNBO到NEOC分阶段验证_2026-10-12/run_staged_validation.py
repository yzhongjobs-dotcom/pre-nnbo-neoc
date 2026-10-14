"""Stage 1: a fixed root and independent policy-quality checks on source NNBO.

Do not silently turn an online learned weight into a root, redefine an optimum,
or run NEOC when the shared anchor has failed the quality screen.
Sources: NNBO equations 24--41 and local source reproduction config/results.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.linalg import qr
from scipy.optimize import root
from scipy.integrate import solve_ivp, cumulative_trapezoid

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'NNBO原文最小复现_2026-10-12'
CFG=dict(Q=500000.,R=1.,Ke=1149.,Be=567.,n=1.2,xe=.01,
         command_alpha=.03,command_A=-1.,uncertainty_bound='abs(1000*e*sin(r))',
         basis=['e^2','e*r','r^2'],domain_penetration_m=[.006,.02],
         domain_reference_m=[.03,.08],selection_grid=[9,9],
         holdout_midpoint_grid=[16,16],root_tolerance_abs=1e-8,
         holdout_relative_residual_allow=.1,
         root_condition_allow=1e8,
         diagnostic_rollout_seconds=3.,diagnostic_rollout_x0=.022,
         diagnostic_rollout_r0=.04,
         gate='all root, physical-force and independent holdout conditions required',
         threshold_scope='local engineering screen, not a theorem or journal threshold')

def grid(midpoints=False):
    if midpoints:
        nd,nr=CFG['holdout_midpoint_grid']
        d=.006+(np.arange(nd)+.5)*(.02-.006)/nd
        r=.03+(np.arange(nr)+.5)*(.08-.03)/nr
    else:
        nd,nr=CFG['selection_grid']
        d=np.linspace(.006,.02,nd);r=np.linspace(.03,.08,nr)
    return np.array([(di,ri) for di in d for ri in r])

def evaluate(w,nodes):
    d,r=nodes.T;e=d+CFG['xe']-r
    g=1/(CFG['Be']*d**CFG['n']);rd=.03-r
    grad_basis=np.stack([np.c_[2*e,np.zeros(len(e))],
                         np.c_[r,e],np.c_[np.zeros(len(e)),2*r]],axis=1)
    grad=np.einsum('nwi,w->ni',grad_basis,w)
    force=-.5*g*grad[:,0]
    flow=np.c_[-CFG['Ke']/CFG['Be']-rd+g*force,rd]
    nu=np.einsum('nwi,ni->nw',grad_basis,flow)
    stage=CFG['Q']*e*e+force*force+(1000*e*np.sin(r))**2
    h=nu@w+stage
    return h,nu,force,stage

def independent_jacobian_check(w,nodes):
    _,analytic,_,_=evaluate(w,nodes)
    cols=[]
    for j in range(3):
        # Scale finite differences to the very different physical coefficients.
        step=1e-5*max(1.,abs(w[j]))
        direction=np.eye(3)[j]*step
        cols.append((evaluate(w+direction,nodes)[0]-evaluate(w-direction,nodes)[0])/(2*step))
    numerical=np.column_stack(cols)
    return float(np.linalg.norm(numerical-analytic)/np.linalg.norm(analytic))

def rollout(w):
    # Diagnostic only: a frozen shared critic, not NEOC or a new optimal solve.
    def rhs(t,y):
        x=y[0];r=.03+(.04-.03)*np.exp(-t);d=x-.01
        if d<=0:raise FloatingPointError('Lost HC contact domain')
        e=x-r;g=1/(567*d**1.2);f=-.5*g*(2*w[0]*e+w[1]*r)
        return [-1149/567+g*f+1000*e*np.sin(r)]
    ts=np.linspace(0,3,3001)
    sol=solve_ivp(rhs,(0.,3.),[.022],t_eval=ts,rtol=1e-9,atol=1e-11,max_step=.001)
    assert sol.success
    x=sol.y[0];r=.03+.01*np.exp(-ts);d=x-.01;e=x-r
    f=-(2*w[0]*e+w[1]*r)/(2*567*d**1.2)
    stage=CFG['Q']*e*e+f*f
    cost=cumulative_trapezoid(stage,ts,initial=0.)
    return np.c_[ts,x,r,e,f,stage,cost],dict(
        final_position_m=float(x[-1]),final_force_N=float(f[-1]),
        tracking_rmse_m=float(np.sqrt(np.mean(e*e))),finite_3s_cost=float(cost[-1]),
        minimum_penetration_m=float(d.min()),
        trajectory_within_root_domain=bool(np.all((d>=.006)&(d<=.02))),
        exploration=False,critic_frozen=True,not_an_optimal_reference=True)

def save_csv(name,head,data):
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
        wr=csv.writer(f);wr.writerow(head);wr.writerows(data)

def main():
    start=time.perf_counter()
    src=json.loads((SOURCE/'results.json').read_text(encoding='utf-8'))
    learned=np.array(next(c['final_weights'] for c in src['cases']
                          if c['Q']==500000 and c['bound']=='tight_known_disturbance'))
    candidates=grid()
    d,r=candidates.T;e=d+.01-r
    phi=np.c_[e*e,e*r,r*r]
    order=qr((phi/np.linalg.norm(phi,axis=0)).T,pivoting=True,mode='economic')[2]
    nodes=candidates[order[:3]]
    holdout=grid(True)
    save_csv('fixed_nodes.csv',['penetration_m','r_m'],nodes)
    save_csv('holdout_nodes.csv',['penetration_m','r_m'],holdout)
    (HERE/'config.json').write_text(json.dumps(CFG,indent=2),encoding='utf-8')
    results=dict(status='RUNNING_STAGE1',config=CFG,python=platform.python_version(),
                 numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
                 source_reproduction_sha256=hashlib.sha256((SOURCE/'results.json').read_bytes()).hexdigest(),
                 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 selected_nodes=nodes.tolist(),basis_rank=int(np.linalg.matrix_rank(phi)),
                 neoc_path_executed=False,parameter_update_comparison_executed=False,
                 spatial_preview_path_executed=False,root_attempts=[])
    # Different starts diagnose branch dependence, with the same fixed equation.
    roots={}
    for name,w0 in [('paper_initial_weights',np.array([20.,20.,20.])),
                    ('source_learned_weights',learned),('zero',np.zeros(3))]:
        fit=root(lambda w:evaluate(w,nodes)[0],w0,
                 jac=lambda w:evaluate(w,nodes)[1],method='hybr',tol=1e-10)
        h,L,f,_=evaluate(fit.x,nodes)
        hh,_,ff,ss=evaluate(fit.x,holdout)
        rr=np.abs(hh)/(1+ss)
        record=dict(initialization=name,solver_success=bool(fit.success),weights=fit.x.tolist(),
                    max_node_residual=float(max(abs(h))),jacobian_condition=float(np.linalg.cond(L)),
                    jacobian_fd_relative_error=independent_jacobian_check(fit.x,nodes),
                    holdout_count=len(holdout),max_holdout_abs_residual=float(max(abs(hh))),
                    max_holdout_relative_residual=float(max(rr)),
                    rms_holdout_relative_residual=float(np.sqrt(np.mean(rr*rr))),
                    holdout_min_force_N=float(min(ff)),holdout_max_force_N=float(max(ff)),
                    value_quadratic_eigenvalues=np.linalg.eigvalsh(
                        [[fit.x[0],fit.x[1]/2],[fit.x[1]/2,fit.x[2]]]).tolist())
        record['root_precision_pass']=record['max_node_residual']<CFG['root_tolerance_abs']
        record['regularity_point_check_pass']=record['jacobian_condition']<CFG['root_condition_allow']
        record['nonnegative_force_check_pass']=record['holdout_min_force_N']>=-1e-10
        record['holdout_quality_pass']=record['max_holdout_relative_residual']<=CFG['holdout_relative_residual_allow']
        record['stage1_screen_pass']=all(record[k] for k in ['root_precision_pass','regularity_point_check_pass',
                                                            'nonnegative_force_check_pass','holdout_quality_pass'])
        assert record['jacobian_fd_relative_error']<1e-6
        results['root_attempts'].append(record);roots[name]=fit.x
        save_csv(f'holdout_{name}.csv',['penetration_m','r_m','H','force_N','stage','relative_abs_H'],
                 np.c_[holdout,hh,ff,ss,rr])
        print(json.dumps(record),flush=True)
    # Audit source learned W separately; successful online integration was not a root.
    h,_,_,_=evaluate(learned,nodes)
    results['learned_weights_node_residual']=float(max(abs(h)))
    # Follow physics-admissible candidate from the original learned start only.
    selected=roots['source_learned_weights']
    results['diagnostic_rollouts']={}
    for name,w in [('force_admissible_root',selected),('frozen_source_learned',learned)]:
        data,metrics=rollout(w)
        save_csv(f'rollout_{name}.csv',['time_s','x_m','r_m','error_m','force_N','stage','finite_cost'],data)
        results['diagnostic_rollouts'][name]=metrics
    # A physically maintained contact at the final command needs nonzero force,
    # while zero force and zero error are both required for zero running cost.
    equilibrium_x=.03
    f_eq=CFG['Ke']*(equilibrium_x-CFG['xe'])**CFG['n']
    results['zero_error_contact_equilibrium']=dict(x=equilibrium_x,r=.03,
                                                  required_force_N=float(f_eq),
                                                  running_cost=float(f_eq*f_eq),
                                                  scope='nominal static contact, not all possible trajectories')
    rc=results['diagnostic_rollouts']['force_admissible_root']['finite_3s_cost']
    lc=results['diagnostic_rollouts']['frozen_source_learned']['finite_3s_cost']
    results['diagnostic_root_cost_increase_over_frozen_learned_fraction']=rc/lc-1
    passed=any(x['stage1_screen_pass'] for x in results['root_attempts'])
    results['status']='STAGE1_PASS' if passed else 'STOP_STAGE1_ANCHOR_POLICY_QUALITY'
    results['stage2_status']='NOT_EXECUTED_STAGE1_NOT_ACCEPTED'
    results['stage3_status']='NOT_EXECUTED_STAGE2_NOT_ACCEPTED'
    results['duration_seconds']=time.perf_counter()-start
    (HERE/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print(results['status'],flush=True)

if __name__=='__main__':main()
