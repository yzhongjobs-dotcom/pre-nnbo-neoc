"""Exact general-control counterexample: no uniformly near-optimal fixed PD.

Local contact perturbation: m e'' + (b+D)e' + (theta+K)e = 0,
force increment u=-K e-D e'. J=int(qe e^2+R u^2)dt, initial
e=0 and v=1 mm/s. This is an exact stationary LQR family diagnostic,
not the original gain-rate HJB or a time-varying global-optimality claim.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import differential_evolution

P=Path(__file__).resolve().parent
M=.2;B=3.;R=1.;QE=1e6;V0=.001;EPS=.05
THETA=np.linspace(100.,2000.,801)

def optimal(th):
    k=np.sqrt(th*th+QE/R)-th
    d=np.sqrt(B*B+2*M*k)-B
    return k,d

def cost(k,d,th):
    # Exact Lyapunov value P22, independently checked below.
    return V0**2*M/(2*(B+d))*(M*(QE+R*k*k)/(th+k)+R*d*d)

def admissible_stiffness(th):
    k,d=optimal(th)
    # Minimize over ALL nonnegative damping, independently at each theta.
    # J_min(K;theta)<= (1+EPS) J_star(theta) implies this quadratic interval.
    a=((B+(1+EPS)*d)**2-B*B)/M
    discriminant=a*a-4*(QE/R-a*th)
    return max(0.,float((a-np.sqrt(discriminant))/2)),float((a+np.sqrt(discriminant))/2)

def main():
    contract=dict(model='local contact perturbation PD/LQR family',mass=M,passive_damping=B,
        force_cost_weight=R,position_cost_weight=QE,initial_velocity=V0,
        theta_range=[100.,2000.],relative_quality_budget=EPS,
        scope='stationary family regret, not time-varying optimality, gain-rate HJB or ultrasound proof')
    (P/'frozen_config.json').write_text(json.dumps(contract,indent=2),encoding='utf8')
    k,d=optimal(THETA);jstar=cost(k,d,THETA)
    def objective(g):return float(np.max(cost(g[0],g[1],THETA)/jstar-1))
    fit=[]
    for seed in [15,16,17]:
        res=differential_evolution(objective,[(0,2000),(.001,40)],seed=seed,
            popsize=20,tol=1e-10,maxiter=400,polish=True)
        fit.append(dict(gains=res.x.tolist(),worst_regret=float(res.fun),success=bool(res.success)))
    best=min(fit,key=lambda a:a['worst_regret']);fixed=np.array(best['gains'])
    # Endpoint intervals disjoint => impossible even with freely retuned damping.
    soft=admissible_stiffness(THETA[0]);hard=admissible_stiffness(THETA[-1])
    assert soft[0]>hard[1]
    lyapunov_errors=[]
    for th in [100.,500.,2000.]:
        for kk,dd in [optimal(th),fixed]:
            A=np.array([[0.,1.],[-(th+kk)/M,-(B+dd)/M]])
            gain=np.array([kk,dd]);Q=np.diag([QE,0.])+R*np.outer(gain,gain)
            value=solve_continuous_lyapunov(A.T,-Q)
            lyapunov_errors.append(abs(value[1,1]*V0**2-cost(kk,dd,th)))
    # Fixed local sensitivity / predeclared-radius maintenance; no candidate checks.
    smin=np.sqrt(B*B+2*M*optimal(THETA[-1])[0])
    mk=(QE/R)/(THETA[0]**2+QE/R)**1.5
    maxkp=abs(THETA[0]/np.sqrt(THETA[0]**2+QE/R)-1)
    md=M*mk/smin+M*M*maxkp**2/smin**3
    radius=min(np.sqrt(2*2./mk),np.sqrt(2*.03/md))
    anchor=THETA[0];ak,ad=optimal(anchor);roots=1;predictions=[]
    for th in THETA:
        if abs(th-anchor)>radius:
            anchor=th;ak,ad=optimal(anchor);roots+=1
        kp=anchor/np.sqrt(anchor*anchor+QE/R)-1
        dp=M*kp/np.sqrt(B*B+2*M*ak)
        predictions.append([ak+kp*(th-anchor),ad+dp*(th-anchor)])
    predictions=np.array(predictions)
    neoc_regret=cost(predictions[:,0],predictions[:,1],THETA)/jstar-1
    result=dict(endpoint_certificate=dict(soft_stiffness_interval=soft,hard_stiffness_interval=hard,
        disjoint_gap_N_m=soft[0]-hard[1],conclusion='No fixed stiffness can meet 5% regret at both endpoints, even allowing different damping.'),
        fixed_numerical_minimax=best,repeated_minimax_fits=fit,
        optimal_endpoint_gains=[list(optimal(100.)),list(optimal(2000.))],
        neoc=dict(a_priori_radius=radius,query_points=len(THETA),anchors=roots,
            maximum_gain_errors=np.max(abs(predictions-np.column_stack((k,d))),axis=0).tolist(),
            worst_regret=float(np.max(neoc_regret)),candidate_result_checks=0),
        closed_form_adaptive=dict(worst_regret=0.,requires_iterative_solve=False,
            conclusion='This strong adaptive baseline matches the optimal family exactly; numerical optimization/NEOC necessity is not established by this example.'),
        lyapunov_max_absolute_cost_error=float(max(lyapunov_errors)),
        scope=contract['scope'])
    np.savetxt(P/'family.csv',np.column_stack((THETA,k,d,predictions,jstar,cost(*fixed,THETA),neoc_regret)),
               delimiter=',',header='theta,optimal_K,optimal_D,NEOC_K,NEOC_D,optimal_cost,fixed_cost,NEOC_regret',comments='')
    (P/'results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2),flush=True)
    assert max(lyapunov_errors)<1e-10 and np.max(neoc_regret)<EPS
    assert np.max(abs(predictions[:,0]-k))<=2. and np.max(abs(predictions[:,1]-d))<=.03

if __name__=='__main__':main()
