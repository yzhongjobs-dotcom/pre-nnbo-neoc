"""Control-facing inverse of the existing locally time-varying SLS model.

F=ke(t)*delta+km*r, rdot=delta_dot-r/tau_m, km locally constant.
For a force difference E and displacement difference c between two trajectories
under the SAME contact operator:
 c_dot + (ke_dot+ke/tau_m)/(ke+km)*c
       = E_dot/(ke+km)+E/((ke+km)*tau_m).
The state c is in metres; changing fitted parameters does not redivide its
history. Unknown initial inverse state remains a transient, not measured depth.
"""
import math
import numpy as np


def inverse_step(c, previous_error, error, theta, tau_m, tw=.256, dt=.008):
    """Exact FOH-error update for midpoint-frozen coefficients of the ODE.

ke_end and ke_slope are the existing endpoint/slope parametrization. The
midpoint freezing is an approximation for varying ke; constant coefficients
and linearly interpolated E are exact. No force derivative is formed.
"""
    ke, slope, km = np.asarray(theta, float)
    ke_dot = slope/tw
    mid_ke = ke-.5*dt*ke_dot
    K = mid_ke+km
    a = (ke_dot+mid_ke/tau_m)/K
    if K <= 0 or mid_ke <= 0 or a <= 0:
        raise ValueError('inverse requires positive gain and locally decaying homogeneous dynamics')
    z = a*dt
    decay = math.exp(-z)
    phi = -math.expm1(-z)/z
    if z < 1e-4:
        # (phi-decay)/a and (1-phi)/a, cancellation-free.
        w0 = dt*(.5-z/3+z*z/8-z**3/30)
        w1 = dt*(.5-z/6+z*z/24-z**3/120)
    else:
        w0, w1 = (phi-decay)/a, (1-phi)/a
    return float(decay*c+phi*(error-previous_error)/K+
                 (w0*previous_error+w1*error)/(K*tau_m))


def step_compliance(theta, tau_m, horizon, tw=.256, subdivisions=12):
    """Unit force-step inverse response, including locally linear ke(t).

RK4 is used only for a fixed-design functional / information calculation;
it is neither a controller trajectory optimizer nor a safety proof.
"""
    ke, slope, km = np.asarray(theta, float)
    derivative = slope/tw
    c = 1/(ke+km)
    if horizon == 0:
        return float(c)
    dt = horizon/subdivisions
    def rhs(t, value):
        now_ke = ke+derivative*t
        K = now_ke+km
        return (1/tau_m-(derivative+now_ke/tau_m)*value)/K
    for j in range(subdivisions):
        t = j*dt
        a = rhs(t,c); b = rhs(t+.5*dt,c+.5*dt*a)
        q = rhs(t+.5*dt,c+.5*dt*b); d = rhs(t+dt,c+dt*q)
        c += dt/6*(a+2*b+2*q+d)
    return float(c)


def compliance_information(r, horizon, age=0.):
    """Propagate the regression design covariance to actual inverse functions.

Conditional on the chosen SLS time constant/model and IID innovation noise.
Not a calibrated constrained posterior or a closed-loop confidence guarantee.
"""
    theta = np.asarray(r['coefficients'],float).copy()
    transform = np.eye(3); transform[0,1]=age/r['tw']
    theta = transform@theta
    cov = transform@np.asarray(r['cov'])@transform.T
    ke,slope,km = theta
    endpoint = ke+slope/r['tw']*horizon
    physical = bool(r['physical'] and min(ke,endpoint)>0 and km>=0 and
                    200<=ke+km<=3e5 and
                    min((slope/r['tw']+ke/r['tau_m']),
                        (slope/r['tw']+endpoint/r['tau_m']))>0)
    result=dict(physical=physical, theta=theta.tolist(), age=age, functions=[])
    if not physical:
        return result|dict(risk=float('inf'))
    inflation = 1+math.sqrt(r.get('dof',max(r['tw']/.008-9,1))*max(r['chi2']-1,0))
    for h in (0.,.008,horizon):
        value=step_compliance(theta,r['tau_m'],h,r['tw'])
        gradient=[]
        for j in range(3):
            step=max(abs(theta[j])*1e-5,.001)
            plus=theta.copy();minus=theta.copy();plus[j]+=step;minus[j]-=step
            gradient.append((step_compliance(plus,r['tau_m'],h,r['tw'])-
                             step_compliance(minus,r['tau_m'],h,r['tw']))/(2*step))
        gradient=np.asarray(gradient)
        sd=math.sqrt(max(float(gradient@cov@gradient),0.))
        result['functions'].append(dict(horizon=h,value=value,sd=sd,relative_sd=sd/max(value,1e-15)))
    result['risk']=inflation*max(f['relative_sd'] for f in result['functions'])
    return result


def reconstruct_increment(r,measured,fc=30.,uncertainty=False,noise=.04,target_force=4.):
    """Warm the SAME inverse with the complete available fitted input window.

Both final displacements use the same model/history. No model-coordinate jump
is included in their difference. The initial inverse state remains assumed;
homogeneous_memory reports its remaining sensitivity, not a certified bound.
"""
    n=len(measured);start=max(1,r['fit_n']-1-round(r['tw']/.008))
    index=np.arange(start,n)
    decay=math.exp(-2*math.pi*fc*.008)
    error=(measured[index,3]-decay*measured[index-1,3])/(1-decay)-target_force
    theta=np.asarray(r['coefficients'],float).copy()
    theta[0]+=theta[1]*(measured[start,0]-r['te'])/r['tw']
    if theta[0]<=0:raise ValueError('nonphysical window initial gain')
    c=error[0]/(theta[0]+theta[2]);values=[c];weight=1.
    force_weights=np.zeros(len(index));force_weights[0]=1/(theta[0]+theta[2])
    for j in range(1,len(index)):
        theta=np.asarray(r['coefficients'],float).copy()
        theta[0]+=theta[1]*(measured[index[j],0]-r['te'])/r['tw']
        c=inverse_step(c,error[j-1],error[j],theta,r['tau_m'],r['tw'])
        if uncertainty:
            A=inverse_step(1.,0.,0.,theta,r['tau_m'],r['tw'])
            B0=inverse_step(0.,1.,0.,theta,r['tau_m'],r['tw'])
            B1=inverse_step(0.,0.,1.,theta,r['tau_m'],r['tw'])
            force_weights=A*force_weights;force_weights[j-1]+=B0;force_weights[j]+=B1
        mid=theta[0]-.5*.008*theta[1]/r['tw']
        a=(theta[1]/r['tw']+mid/r['tau_m'])/(mid+theta[2])
        weight*=math.exp(-a*.008);values.append(c)
    result=dict(increment=float(values[-1]-values[-2]),current=float(values[-1]),
                previous=float(values[-2]),homogeneous_memory=weight,history_ticks=len(index))
    if uncertainty:
        gradient=[]
        for j in range(3):
            step=max(abs(r['coefficients'][j])*1e-5,.001)
            plus=dict(r);minus=dict(r)
            plus['coefficients']=np.asarray(r['coefficients'],float).copy()
            minus['coefficients']=plus['coefficients'].copy()
            plus['coefficients'][j]+=step;minus['coefficients'][j]-=step
            p=reconstruct_increment(plus,measured,fc,target_force=target_force)
            m=reconstruct_increment(minus,measured,fc,target_force=target_force)
            gradient.append((p['current']-m['current'])/(2*step))
        gradient=np.asarray(gradient)
        sigma_model=math.sqrt(max(float(gradient@r['cov']@gradient),0.))
        # Undo the force innovation filter algebra in the influence row. This
        # includes adjacent-sample correlations rather than assuming IID fin.
        influence=np.zeros(len(index)+1)
        influence[1:]+=force_weights/(1-decay)
        influence[:-1]-=decay*force_weights/(1-decay)
        sigma_force=noise*float(np.linalg.norm(influence))
        # No independence between fitted parameters and measured force is
        # claimed: Cauchy gives the conservative sum of marginal SDs within
        # the local linearized model. Unknown initial-state bias is NOT covered.
        result.update(sd_model=sigma_model,sd_force=sigma_force,
                      variance_reference=(sigma_model+sigma_force)**2)
    return result
