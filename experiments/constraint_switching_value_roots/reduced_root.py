"""Critical-constraint failure and reduced Hamiltonian-root diagnosis.

F(w,theta)=1-2theta*w+min_|u|<=.4 (u^2+2wu) at collocation x=1.
Full complementarity KKT root is singular at weak activity; reduced value root
is C1,1 and regular. This is a one-node critic proxy, not a complete R15 proof.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad

P=Path(__file__).resolve().parent;UB=.4;LO=.9;HI=1.2;BUDGET=.002
SWITCH=(1-UB**2)/(2*UB)

def root(th):
    w=np.sqrt(th*th+1)-th
    return w if w<=UB else (1+UB**2)/(2*(th+UB))

def sensitivity(th):
    if th>=SWITCH:return th/np.sqrt(th*th+1)-1
    return -(1+UB**2)/(2*(th+UB)**2)

def F(w,th):
    u=np.clip(-w,-UB,UB)
    return 1-2*th*w+u*u+2*w*u

def jacobian(w,th):
    u=np.clip(-w,-UB,UB);lam=max(0.,2*(w-UB))
    return np.array([[-2*th+2*u,2*u+2*w,0.],[2.,2.,-1.],[0.,-lam,-u-UB]])

def policy_cost(w,th,x0):
    xs=min(x0,UB/w)
    unsat=(1+w*w)*xs*xs/(2*(th+w))
    tail=quad(lambda x:(x*x+UB*UB)/(th*x+UB),xs,x0,epsabs=1e-12,epsrel=1e-12)[0]
    return unsat+tail

def main():
    m=max((1+UB**2)/(LO+UB)**3,1/(SWITCH*SWITCH+1)**1.5)
    radius=np.sqrt(2*BUDGET/m)
    paths=dict(release=np.linspace(LO,HI,301),activate=np.linspace(HI,LO,301))
    results={};allrows=[]
    for name,path in paths.items():
        anchor=path[0];wa=root(anchor);ca=sensitivity(anchor);ua=np.clip(-wa,-UB,UB)
        du=0. if wa>UB else -ca;anchors=1;rows=[];maxgap=0.;errors=[]
        for th in path:
            if abs(th-anchor)>radius:
                anchor=th;wa=root(anchor);ca=sensitivity(anchor);ua=np.clip(-wa,-UB,UB)
                du=0. if wa>UB else -ca;anchors+=1
            delta=th-anchor;wp=wa+ca*delta
            true=root(th);up=np.clip(-wp,-UB,UB)
            action_raw=ua+du*delta;action_projected=np.clip(action_raw,-UB,UB)
            exact_gain=np.sqrt(th*th+1)-th
            gap=max(policy_cost(wp,th,x)/policy_cost(exact_gain,th,x)-1 for x in [.1,.5,1.,2.])
            maxgap=max(maxgap,gap);err=abs(wp-true);errors.append(err)
            rows.append([th,true,wp,err,up,np.clip(-true,-UB,UB),action_raw,action_projected,gap])
        a=np.array(rows);allrows.append(a)
        results[name]=dict(anchors=anchors,queries=len(path),maximum_root_error=max(errors),
            maximum_root_bound=m*radius*radius/2,raw_action_prediction_violation_count=int(np.sum(abs(a[:,6])>UB+1e-12)),
            projected_action_prediction_max_error=float(np.max(abs(a[:,7]-a[:,5]))),
            value_recovery_max_action_error=float(np.max(abs(a[:,4]-a[:,5]))),
            worst_true_policy_cost_gap=maxgap,value_recovery_constraint_violations=int(np.sum(abs(a[:,4])>UB+1e-12)))
        np.savetxt(P/(name+'.csv'),a,delimiter=',',header='theta,root,predicted_root,root_error,recovered_action,root_action,raw_predicted_action,projected_predicted_action,true_cost_gap',comments='')
    matrix=jacobian(root(SWITCH),SWITCH)
    critical=dict(theta=SWITCH,full_KKT_smallest_singular_value=float(np.linalg.svd(matrix,compute_uv=False)[-1]),
        full_KKT_rank=int(np.linalg.matrix_rank(matrix)),reduced_root_derivative=-2*(SWITCH+UB),
        minimum_reduced_derivative_magnitude=2*(LO+UB),
        left_second_derivative=(1+UB**2)/(SWITCH+UB)**3,
        right_second_derivative=1/(SWITCH*SWITCH+1)**1.5,
        root_first_derivative_left=sensitivity(np.nextafter(SWITCH,0.)),
        root_first_derivative_right=sensitivity(np.nextafter(SWITCH,np.inf)))
    fd=[]
    for th in [LO,1.01,SWITCH,1.10,HI]:
        ep=1e-6;fd.append(abs((root(th+ep)-root(th-ep))/(2*ep)-sensitivity(th)))
    result=dict(critical=critical,budget=BUDGET,curvature_bound=m,a_priori_radius=radius,
        maximum_sensitivity_difference=max(fd),paths=results,
        scope='finite one-node Hamiltonian root; exact constrained scalar policy used only for performance audit',
        novelty='candidate mechanism, not certified novel; projected/multisegment NEOC and reduced-value implicit differentiation are strong comparators')
    (P/'results.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2),flush=True)
    assert critical['full_KKT_rank']<3 and critical['minimum_reduced_derivative_magnitude']>0
    assert max(fd)<1e-6
    for r in results.values():
        assert r['maximum_root_error']<=BUDGET and r['value_recovery_constraint_violations']==0
        assert r['worst_true_policy_cost_gap']<.001

if __name__=='__main__':main()
