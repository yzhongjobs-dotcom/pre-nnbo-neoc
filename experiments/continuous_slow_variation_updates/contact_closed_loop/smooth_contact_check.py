"""Independent smooth-critic control experiment, not a replacement R15 design.

Linear local soft-contact force control gives a quadratic value critic with an
exact smooth Riccati root. Its optimal force feedback is a normal PD impedance.
This isolates smooth NEOC from the nodal carrier's interpolation branches.
Actual impedance gains are interpolated causally, not jumped or clipped.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import time
import numpy as np
from scipy.linalg import solve_continuous_are,solve_continuous_lyapunov

HERE=Path(__file__).resolve().parent/'02_smooth_quadratic_critic_comparator'
M=.2;BE=3.;FD=.14;R=1/.15**2;Q=np.diag([1/.0005**2,1/.008**2])
B=np.array([[0.],[1/M]]);G=B@B.T/R
LOW=300.;HIGH=450.;DT=.001


def plant_matrix(theta):return np.array([[0.,1.],[-theta/M,-BE/M]])


def exact(theta):
    k=math.sqrt(theta*theta+Q[0,0]/R)-theta
    d=math.sqrt(BE*BE+2*M*k+Q[1,1]/R)-BE
    return np.array([[R*(BE*k+(theta+k)*d),R*M*k],[R*M*k,R*M*d]])


def implicit_sensitivity(P,theta):
    a=plant_matrix(theta)-G@P
    da=np.array([[0.,0.],[-1/M,0.]])
    C=solve_continuous_lyapunov(a.T,-(da.T@P+P@da))
    return (C+C.T)/2


def newton(P,theta,jac_theta=None,jac_P=None):
    a=plant_matrix(theta)
    F=a.T@P+P@a-P@G@P+Q
    j=plant_matrix(theta if jac_theta is None else jac_theta)-G@(P if jac_P is None else jac_P)
    delta=solve_continuous_lyapunov(j.T,-F)
    return (P+delta+(P+delta).T)/2


def theta_at(t,duration,kind):
    phase=min(max(t/duration,0.),1.)
    return LOW+(HIGH-LOW)*(phase if kind=='ramp' else (1-math.cos(2*math.pi*phase))/2)


def surface(t):
    # Prescribed surface motion, equally measured by all controllers.
    amplitude=.00008;omega=2*math.pi
    return amplitude*math.sin(omega*t),amplitude*omega*math.cos(omega*t)


def curvature_budget():
    # Analytic upper bounds throughout [LOW,HIGH], no sampled curvature claim.
    c=Q[0,0]/R
    mk=c/(LOW*LOW+c)**1.5
    lk=1-LOW/math.sqrt(LOW*LOW+c)
    kmin=math.sqrt(HIGH*HIGH+c)-HIGH
    umin=math.sqrt(BE*BE+2*M*kmin+Q[1,1]/R)
    md=M*mk/umin+M*M*lk*lk/umin**3
    radius=min(math.sqrt(2*.1/mk),math.sqrt(2*.002/md))
    return dict(stiffness_second_derivative_bound=mk,damping_second_derivative_bound=md,
                stiffness_prediction_error_allow=.1,damping_prediction_error_allow=.002,
                radius_N_per_m=radius,
                scope='frozen-root feedback gain error only, NOT time-varying global optimality')


def stability_budget(kind,duration,period,budget):
    # V=x^T P*(theta(t)) x. Q scales physical coordinates consistently.
    # Derivative <= -(1 - drift - feedback_perturbation) x^T Q x
    # for the homogeneous error dynamics, including causal gain interpolation.
    c=Q[0,0]/R;lk=1-LOW/math.sqrt(LOW*LOW+c)
    kmax=math.sqrt(LOW*LOW+c)-LOW;kmin=math.sqrt(HIGH*HIGH+c)-HIGH
    dmax=math.sqrt(BE*BE+2*M*kmax+Q[1,1]/R)-BE
    umin=math.sqrt(BE*BE+2*M*kmin+Q[1,1]/R)
    ld=M*lk/umin;smax=math.sqrt(HIGH*HIGH+c)
    p11prime=R*(BE*lk+dmax+smax*ld)
    p12prime=R*M*lk;p22prime=R*M*ld
    weighted_derivative=math.sqrt((p11prime/Q[0,0])**2+
        2*p12prime**2/(Q[0,0]*Q[1,1])+(p22prime/Q[1,1])**2)
    rate=(HIGH-LOW)/duration*(1 if kind=='ramp' else math.pi)
    drift=weighted_derivative*rate
    # Midpoint gain samples are held over DT; include their half-step offset.
    lag=period+DT/2
    delta_k=.1+lk*rate*lag;delta_d=.002+ld*rate*lag
    pb=R*math.sqrt(kmax*kmax/Q[0,0]+dmax*dmax/Q[1,1])
    perturb=2*pb*math.sqrt(delta_k*delta_k/Q[0,0]+delta_d*delta_d/Q[1,1])
    margin=1-drift-perturb
    assert margin>0
    return dict(maximum_theta_rate=rate,weighted_P_derivative_bound=weighted_derivative,
        drift_fraction=drift,feedback_perturbation_fraction=perturb,
        homogeneous_decay_margin=margin,
        scope='local_NEOC homogeneous linear error stability; bounded forcing remains, force/energy verified by sampled rollout')


def simulate(method,kind,duration,period,budget):
    start=time.perf_counter()
    if method=='fixed_midrange':P=exact((LOW+HIGH)/2)
    elif method in ('analytic_full','fixed_gain'):P=exact(LOW)
    else:P=solve_continuous_are(plant_matrix(LOW),B,Q,[[R]])
    initial=time.perf_counter()-start;prep=0.;updates_time=0.;policy_time=0.;anchors=1;derivatives=0
    anchor=P.copy();anchor_theta=LOW;last_theta=LOW;lastP=P.copy();C=None
    if method=='local_NEOC':
        start=time.perf_counter();C=implicit_sensitivity(P,LOW);prep=time.perf_counter()-start;derivatives=1
    gains=(B.T@P/R).ravel();state=np.array([FD/LOW,0.,.02]);rows=[];cost=0.;events=[];errors=[]
    for j in range(round(duration/period)):
        now=j*period;theta=theta_at(now,duration,kind);reset=False
        before=time.perf_counter()
        if j:
            if method=='local_NEOC':
                if abs(theta-anchor_theta)>budget['radius_N_per_m']:
                    P=solve_continuous_are(plant_matrix(theta),B,Q,[[R]])
                    C=implicit_sensitivity(P,theta);anchor=P.copy();anchor_theta=theta
                    anchors+=1;derivatives+=1;reset=True
                else:P=anchor+C*(theta-anchor_theta)
            elif method=='continuous_NEOC':
                C=implicit_sensitivity(P,last_theta);derivatives+=1
                P=P+C*(theta-last_theta)
            elif method=='one_Newton':P=newton(P,theta)
            elif method=='cached_chord':P=newton(P,theta,last_theta,lastP)
            elif method=='full_CARE':P=solve_continuous_are(plant_matrix(theta),B,Q,[[R]]);anchors+=1
            elif method=='analytic_full':P=exact(theta);anchors+=1
        target=(B.T@P/R).ravel();oldgains=gains.copy();lastP=P.copy()
        elapsed=time.perf_counter()-before;updates_time+=elapsed
        events.append([now,theta,elapsed,reset,target[0],target[1]])
        ref=exact(theta);refg=(B.T@ref/R).ravel()
        errors.append([now,float(np.linalg.norm(P-ref)/np.linalg.norm(ref)),
                       abs(float(target[0]-refg[0])),abs(float(target[1]-refg[1]))])
        # Audit is retrospective; error/eigenvalue checks never decide publication.
        if method=='local_NEOC':
            assert errors[-1][2] <= .1+1e-8 and errors[-1][3] <= .002+1e-8
        for i in range(round(period/DT)):
            t=now+i*DT;beta=(i+.5)*DT/period
            k,d=(1-beta)*oldgains+beta*target;kd=(target[0]-oldgains[0])/period
            th=theta_at(t,duration,kind);s,sv=surface(t)
            start=time.perf_counter()
            commanded_error=state[0]-s-FD/theta
            relv=state[1]-sv;u=-k*commanded_error-d*relv
            policy_time+=time.perf_counter()-start
            fn=th*(state[0]-s)+BE*relv
            true_error=state[0]-s-FD/th
            stage=Q[0,0]*true_error**2+Q[1,1]*relv**2+R*u*u
            rows.append([t,th,theta,*state,fn,k,d,u,stage,kd])
            def flow(tt,x):
                ss,vs=surface(tt);ttheta=theta_at(tt,duration,kind)
                # Causal sampled material descriptor; continuous common surface.
                uu=-k*(x[0]-ss-FD/theta)-d*(x[1]-vs)
                ff=ttheta*(x[0]-ss)+BE*(x[1]-vs)
                error=x[0]-ss-FD/theta
                return np.array([x[1],(FD+uu-ff)/M,
                                 .3*(BE+d)*(x[1]-vs)**2-.5*error**2*kd])
            aa=flow(t,state);bb=flow(t+DT/2,state+DT*aa/2)
            cc=flow(t+DT/2,state+DT*bb/2);dd=flow(t+DT,state+DT*cc)
            state+=DT*(aa+2*bb+2*cc+dd)/6;cost+=stage*DT
        gains=target;last_theta=theta
    trace=np.array(rows);err=np.array(errors)
    ok=bool(np.all(trace[:,6]>0) and np.all(trace[:,6]<2*FD) and
            np.all(trace[:,7]>0) and np.all(trace[:,8]>0) and
            np.all(abs(trace[:,11])<=500) and np.all(trace[:,5]>=.001))
    name=f'{kind}_{duration:g}s_{period:g}s_{method}'
    for filename,data,head in [('trace_'+name+'.csv',trace,['time_s','true_theta','sampled_theta','d_m','v_m_s','energy','force_N','gain_K','gain_D','force_command_delta_N','stage','K_rate']),
                               ('events_'+name+'.csv',events,['time_s','theta','update_seconds','anchor_reset','gain_K','gain_D'])]:
        with (HERE/filename).open('w',newline='',encoding='utf-8') as f:
            writer=csv.writer(f);writer.writerow(head);writer.writerows(data)
    return dict(method=method,cost=float(cost),force_rmse_N=float(np.sqrt(np.mean((trace[:,6]-FD)**2))),
        sampled_constraints_pass=ok,maximum_P_relative_error=float(np.max(err[:,1])),
        maximum_K_error=float(np.max(err[:,2])),maximum_D_error=float(np.max(err[:,3])),
        anchor_solves=anchors,sensitivity_evaluations=derivatives,candidate_acceptance_checks=0,
        initial_seconds=initial,prepare_seconds=prep,update_seconds=updates_time,
        policy_seconds=policy_time,total_controller_seconds=initial+prep+updates_time+policy_time,
        maximum_query_seconds=max(row[2] for row in events))


def main():
    HERE.mkdir(exist_ok=True);budget=curvature_budget();audit=[]
    for theta in [300.,375.,450.]:
        P=solve_continuous_are(plant_matrix(theta),B,Q,[[R]]);truth=exact(theta)
        assert np.linalg.norm(P-truth)/np.linalg.norm(truth)<1e-9
        eps=.001;fd=(exact(theta+eps)-exact(theta-eps))/(2*eps)
        C=implicit_sensitivity(P,theta)
        relative=float(np.linalg.norm(C-fd)/np.linalg.norm(fd));assert relative<1e-7
        audit.append(dict(theta=theta,implicit_to_independent_derivative_relative_error=relative))
    methods=['local_NEOC','continuous_NEOC','cached_chord','one_Newton','full_CARE','analytic_full','fixed_gain','fixed_midrange']
    result=dict(model='smooth linear local contact, exact quadratic single critic; force feedback mapped to PD gains',
                full_R15=False,future_descriptor_used=False,budget=budget,derivative_audit=audit,scenarios=[])
    for kind,duration,period in [('ramp',3.,.25),('ramp',6.,.25),('ramp',6.,.5),('cycle',6.,.25)]:
        scores=[simulate(m,kind,duration,period,budget) for m in methods]
        ref=next(row for row in scores if row['method']=='analytic_full')
        for row in scores:
            row['cost_extra_fraction']=row['cost']/ref['cost']-1
            row['force_rmse_extra_fraction']=row['force_rmse_N']/ref['force_rmse_N']-1
        result['scenarios'].append(dict(path=kind,duration_s=duration,period_s=period,
            offline_local_NEOC_stability_budget=stability_budget(kind,duration,period,budget),scores=scores))
        print(kind,duration,period,[(r['method'],r['cost_extra_fraction'],r['maximum_K_error']) for r in scores],flush=True)
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
