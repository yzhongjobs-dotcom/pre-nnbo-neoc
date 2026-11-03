"""Control-relevant dynamic contact response, with delayed causal publication.

125-Hz original controller; contact fitting every four ticks, released four
ticks later. The delay is simulated; this is not real-time certification.
No hidden plant columns or altered excitation/feedback gains.
"""
import math
import time
import numpy as np
from joint_policy import JointPolicy
from common import FD
import memory_contact as MC
import structured_memory as SM
from dynamic_compliance import inverse_step, compliance_information, reconstruct_increment


def free_forecast(r,measured,tau):
    j=np.arange(r['fit_n']-1,len(measured)-1)
    if not len(j):return dict(score=float('inf'),count=0)
    B,U,_=MC.features(measured[j,0],measured[j,2],measured[j+1,2],r['te'],r['ye'],r['tw'],
                      tau,r['tau_m'],r['origin'],r['qend'])
    J=np.column_stack([B,U]);pred=J@np.r_[r['beta'],r['coefficients']]
    variance=r['variance']+np.einsum('ij,jk,ik->i',J,r['full_cov'],J)
    observed=measured[j+1,3]-math.exp(-2*math.pi*30*.008)*measured[j,3]
    return dict(score=float(np.mean(((pred-observed)**2/variance)[-4:])),count=min(len(j),4))


class ResponsePolicy(JointPolicy):
    def __init__(self,response_model='shared',conversion='dynamic',history=False,uncertainty=False,**kw):
        super().__init__(rate=False,**kw)
        self.response_model=response_model;self.conversion=conversion
        self.history=history
        self.uncertainty=uncertainty

    def reset(self,p):
        super().reset(p)
        self.pending=None;self.active_response=None;self.response_fits=[]
        self.response_trace=[];self.response_used=0;self._compliance=None

    def _stiffness_est(self,o,n):
        feedback=super()._stiffness_est(o,n)
        if self.pending is not None and n>=self.pending[0]:
            _,r=self.pending;self.pending=None
            age=float(o.log[n-1,0])-r['te']
            info=compliance_information(r,(self.ndel+.5)*self.Tc+self.tau,age)
            forecast=(SM.forecast_reference if self.response_model=='shared' else free_forecast)(
                r,o.log[:n,[0,1,2,5,11,12,13]],self.tau)
            accepted=bool(info['risk']<=math.log(1.2) and forecast['score']<=self.params['chi2_max'] and
                          r['chi2']<=self.params['chi2_max'] and r.get('optimizer_success',True))
            self.response_fits.append(dict(time=float(o.log[n-1,0]),fit_time=r['te'],accepted=accepted,
                inverse_risk=info['risk'],original_probe_risk=r['risk'],chi2=r['chi2'],
                forecast_score=forecast['score'],coefficients=np.asarray(r['coefficients']).tolist(),
                fit_seconds=r['fit_seconds'],optimizer_success=r.get('optimizer_success',True),
                conditional_covariance=np.asarray(r['cov']).tolist()))
            if accepted:self.active_response=r
        if o.log[n-1,0]>=self.warm+.064+self.information_window*.008 and n%4==0:
            measured=o.log[:n,[0,1,2,5,11,12,13]]
            tic=time.perf_counter()
            if self.response_model=='shared':r=SM.fit(measured,tau=self.tau,noise=self.params['noise'])
            else:r=MC.fit(measured,model='memoryc_040',tau=self.tau,noise=self.params['noise'])
            elapsed=time.perf_counter()-tic
            if r is not None:
                r.update(fit_n=n,fit_seconds=elapsed)
                self.pending=(n+4,r)
        return feedback

    def drift_pred(self,o,n):
        if self._motion_n==n:return self._motion_prediction
        previous_n=self._motion_n;self._motion_n=n;self._motion_prediction=None
        if n<2:return None
        j=n-1;F,Fprevious=o.log[j,5],o.log[j-1,5]
        y0,y1=o.log[j-1,2],o.log[j,2]
        y=float(y0+self._encoder_fraction*(y1-y0))
        decay=(1-self.afc)**self.ns;fin=(F-decay*Fprevious)/(1-decay)
        k=self.k_good if np.isfinite(self.k_good) and self.k_good>0 else self.kbar
        R=(self.params['noise']/k)**2*self.motion_params['r_scale']+self.r_floor**2
        if self.params.get('gain_noise') and self._y_prev is not None and np.isfinite(self.kg_rel):
            R+=(self.kg_rel*(y-self._y_prev))**2
        if self.gate and self.fit is not None and np.isfinite(self.fit.get('chi2',np.nan)):
            R*=max(1.,self.fit['chi2']/self.params['chi2_max'])
        self.tracker.R=max(R,np.finfo(float).tiny)
        used=False;inverse_risk=float('inf');age=float('nan');delta_c=0.;initial_memory=float('nan');extra_R=0.
        if self._fin_prev is None or previous_n!=n-1:
            self.Z=0.;self.tracker.reset(0.);self.q_prev=None;self._compliance=None
        else:
            delta_c=(fin-self._fin_prev)/k
            r=self.active_response
            if r is not None:
                age=float(o.log[j,0])-r['te']
                info=compliance_information(r,(self.ndel+.5)*self.Tc+self.tau,age)
                inverse_risk=info['risk']
                physical=info['physical'] and age<=r['tw']
                if physical:
                    theta=np.asarray(info['theta'])
                    forecast=(SM.forecast_reference if self.response_model=='shared' else free_forecast)(
                        r,o.log[:n,[0,1,2,5,11,12,13]],self.tau)
                    old_c=self._compliance
                    if old_c is None:old_c=(self._fin_prev-FD)/(theta[0]+theta[2])
                    new_c=inverse_step(old_c,self._fin_prev-FD,fin-FD,theta,r['tau_m'],r['tw'])
                    self._compliance=new_c
                    used=bool(info['risk']<=math.log(1.2) and forecast['score']<=self.params['chi2_max'])
                    if used:
                        if self.conversion=='dynamic':
                            if self.history:
                                try:
                                    history=reconstruct_increment(r,o.log[:n,[0,1,2,5,11,12,13]],self.fc,
                                        uncertainty=self.uncertainty,noise=self.params['noise'],target_force=FD)
                                    delta_c=history['increment'];initial_memory=history['homogeneous_memory']
                                    extra_R=history.get('variance_reference',0.)
                                except ValueError:used=False
                            else:delta_c=new_c-old_c
                        else:
                            w=2*math.pi*20.8*r['tau_m']
                            scalar=theta[0]+theta[2]*w*w/(1+w*w)
                            delta_c=(fin-self._fin_prev)/scalar
                        if used:self.response_used+=1
                else:self._compliance=None
            self.Z+=delta_c-(y-self._y_prev)
            if used:self.tracker.R=max(self.tracker.R,extra_R)
            self.tracker.update(self.Z)
        self.response_trace.append([float(o.log[j,0]),float(used),inverse_risk,age,delta_c,float(self.Z),
                                    float(self._compliance) if self._compliance is not None else float('nan'),initial_memory,extra_R])
        self._fin_prev,self._y_prev=float(fin),y
        L=(self.ndel+.5)*self.Tc;s,sd=self.tracker.predict(L)
        self.v_cmd_motion=-float(self.tracker.x[1]) if self.tracker.initialized else 0.
        self._motion_prediction=(None,s+self.tau*sd)
        if not np.isfinite(self.k_good) or self.k_good<=0:self.q_prev=None
        return self._motion_prediction
