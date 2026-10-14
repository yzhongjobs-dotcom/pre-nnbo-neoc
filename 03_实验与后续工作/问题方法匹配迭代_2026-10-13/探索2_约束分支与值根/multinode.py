"""Multi-node smooth critic cross-check with independent constrained policy cost.

Derivative-cardinal cubic basis: V_x at collocation nodes equals w_i. Integrate
the spline to get V(0)=0. Root equations are the minimized scalar HJB at nodes.
No claim that spline matches the exact value between nodes or that this is R15.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

P=Path(__file__).resolve().parent;UB=.4;NODES=np.linspace(0,2,17);XS=NODES[1:]
LO=.9;HI=1.2

def exact_derivatives(theta):
    k=np.sqrt(theta*theta+1)-theta
    p=2*k*XS
    active=p>2*UB
    p[active]=(XS[active]**2+UB**2)/(theta*XS[active]+UB)
    c=2*XS*(theta/np.sqrt(theta*theta+1)-1)
    c[active]=-XS[active]*(XS[active]**2+UB**2)/(theta*XS[active]+UB)**2
    return p,c

def policy_cost(p,theta,x0):
    curve=CubicSpline(NODES,np.r_[0,p],bc_type='natural')
    def integrand(x):
        if x==0:return 0.
        u=np.clip(-curve(x)/2,-UB,UB)
        return (x*x+u*u)/(theta*x-u)
    return quad(integrand,0,x0,epsabs=1e-11,epsrel=2e-7,limit=150,points=list(NODES[(NODES>0)&(NODES<x0)]))[0]

def exact_cost(theta,x0):
    k=np.sqrt(theta*theta+1)-theta;xs=min(x0,UB/k)
    return k*xs*xs+quad(lambda x:(x*x+UB*UB)/(theta*x+UB),xs,x0,epsabs=1e-12,epsrel=1e-12)[0]

def main():
    # Conservative uniform node-wise curvature bounds, including all switches.
    mp=np.maximum(2*XS/(LO*LO+1)**1.5,2*XS**2*(XS**2+UB**2)/(LO*XS+UB)**3)
    budget=.002;radius=float(np.sqrt(2*budget/max(mp)))
    output={};grid=np.linspace(0,2,1001)
    for name,path in dict(release=np.linspace(LO,HI,121),activate=np.linspace(HI,LO,121)).items():
        anchor=path[0];wa,ca=exact_derivatives(anchor);ua=np.clip(-wa/2,-UB,UB)
        action_c=np.where(wa<2*UB,-ca/2,0.);count=1;maxerr=0.;maxgap=0.;maxconfig=0.;violations=0
        min_derivative=1e9;cond=0.;actionerror=0.;rawerror=0.;rows=[]
        for i,theta in enumerate(path):
            if abs(theta-anchor)>radius:
                anchor=theta;wa,ca=exact_derivatives(anchor);ua=np.clip(-wa/2,-UB,UB)
                action_c=np.where(wa<2*UB,-ca/2,0.);count+=1
            truth,c=exact_derivatives(theta);pred=wa+ca*(theta-anchor)
            curve=CubicSpline(NODES,np.r_[0,pred],bc_type='natural')
            min_derivative=min(min_derivative,float(curve(grid).min()))
            u=np.clip(-pred/2,-UB,UB);up=np.clip(ua+action_c*(theta-anchor),-UB,UB)
            viol=int(np.sum(abs(ua+action_c*(theta-anchor))>UB+1e-12));violations+=viol
            err=float(np.max(abs(pred-truth)));maxerr=max(maxerr,err)
            actionerror=max(actionerror,float(np.max(abs(u-np.clip(-truth/2,-UB,UB)))))
            rawerror=max(rawerror,float(np.max(abs(up-np.clip(-truth/2,-UB,UB)))))
            diagonal=-theta*XS+np.clip(-truth/2,-UB,UB)
            cond=max(cond,float(np.max(abs(diagonal))/np.min(abs(diagonal))))
            if i%6==0:
                for x0 in [.1,.5,1.,2.]:
                    exact=exact_cost(theta,x0);base=policy_cost(truth,theta,x0);val=policy_cost(pred,theta,x0)
                    maxconfig=max(maxconfig,base/exact-1);maxgap=max(maxgap,val/exact-1)
            rows.append([theta,err,viol,float(np.max(abs(u-np.clip(-truth/2,-UB,UB))))])
        output[name]=dict(nodes=len(XS),queries=len(path),anchors=count,critic_derivative_error=maxerr,
            derivative_budget=budget,worst_policy_cost_gap=maxgap,worst_representation_cost_gap=maxconfig,
            raw_action_prediction_node_violations=violations,value_recovery_action_error=actionerror,
            projected_action_prediction_error=rawerror,minimum_interpolated_derivative=min_derivative,
            reduced_root_condition_number_max=cond,scope='scalar constrained force-action HJB; not gain-rate/energy R15')
        np.savetxt(P/(name+'_multinode.csv'),rows,delimiter=',',header='theta,gradient_error,raw_node_violation_count,recovered_action_error',comments='')
        assert maxerr<=budget and min_derivative>=-1e-12 and maxgap<.005
    switches=XS/(2*UB)-UB/(2*XS)
    output['critical_nodes']=[dict(x=float(x),theta=float(t),full_complementarity_root_singular=True,
        reduced_diagonal=-float(t*x+UB)) for x,t in zip(XS,switches) if LO<=t<=HI]
    output['radius']=radius
    (P/'multinode_results.json').write_text(json.dumps(output,indent=2),encoding='utf8')
    print(json.dumps(output,indent=2),flush=True)

if __name__=='__main__':main()
