"""Contact-model-only closed-loop pilot; inherited feedback/IMM/supervisor unchanged.

This is a synchronous software mechanism test. Slow-fit wall time is measured;
zero computational delay is NOT a deployment claim. Optional update latency
tests hold models until the declared time, with the same original robot plant.
"""
import math
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'method_redesign'))
sys.path.insert(0,str(HERE.parent/'motion_aware_scan'))
from v1_servo_id import ServoID
from dynamic_contact import ContactFit,Memory,Signals


class ContactPolicy(ServoID):
    def __init__(self, spatial=False, memory=Memory(.5,.2), interface='gain', contact='foundation',
                 latency=0.,update_ticks=16,fallback='baseline',data_confidence=True,blocks=4,
                 validate_next=False,stable_jacobian=False,propagate_contact=False,
                 uncertainty_domain='increment',**kw):
        super().__init__(**kw)
        self.contact_args=dict(spatial=spatial,memory=memory,contact=contact,data_confidence=data_confidence,
                               blocks=blocks,stable_jacobian=stable_jacobian)
        self.interface,self.latency,self.update_ticks=interface,float(latency),int(update_ticks)
        self.fallback=fallback
        self.validate_next=bool(validate_next)
        self.propagate_contact=bool(propagate_contact)
        if uncertainty_domain not in ('increment','level'):raise ValueError('uncertainty_domain')
        self.uncertainty_domain=uncertainty_domain

    def reset(self,platform):
        super().reset(platform)
        self.model_fit=ContactFit(**self.contact_args)
        self.signals=Signals(self.contact_args['memory'],fc=self.fc,tau_servo=self.tau)
        self.accepted_fit=None;self.pending=None;self.model_ready=False
        self.fit_records=[];self.model_ticks=0;self.contact_ticks=0;self.model_current=False
        self.signal_n=None
        self.baseline_good=float('nan')
        self.contact_rejections=[];self.contact_trace=[]
        self.observation_variances=[];self.unobserved_corrections=0;self.current_fit=None
        self.last_contact_motion_n=None;self.observation_fit_time=None;self.observation_anchor=None
        self.anchor_gradient=None;self.cumulative_gradient=None

    def step(self,o):
        n=o.n
        if self.signal_n!=n:
            self.signals.tau_servo=self.tau
            self.signals.append(o.log[n-1,[0,1,2,5,11,12,13]])
            self.signal_n=n
        if o.log[n-1,5]<.3:
            self.accepted_fit=None;self.pending=None;self.model_ready=False
            self.model_fit.previous=None
            self.model_fit.history=[];self.model_fit.retired_information=None;self.model_fit.retired_vector=None
        return super().step(o)

    def _stiffness_est(self,o,n):
        previous_good=self.k_good
        if self.fallback=='baseline':self.k_good=self.baseline_good
        feedback=super()._stiffness_est(o,n)
        if self.fallback=='baseline':self.baseline_good=self.k_good
        # Original projection updates feedback only. Model measurement gain has
        # its own residual, variance and acceptance score, not the old gate.
        if self.model_ready and self.fallback=='hold':self.k_good=previous_good
        now=float(o.log[n-1,0]);self.model_current=False;self.current_fit=None
        if self.pending is not None and now>=self.pending[0]:
            self.accepted_fit=self.pending[1];self.pending=None
        if n>=64 and n%self.update_ticks==0:
            obs=self.signals.window()
            if obs is not None:
                hint=self.k_hat if np.isfinite(self.k_hat) and self.k_hat>0 else self.kbar
                try:
                    validation=None
                    if self.validate_next and self.accepted_fit is not None:
                        validation=self.model_fit.validate_contact(self.accepted_fit,obs)
                        if not validation['valid']:
                            self.contact_rejections.append(dict(time=now,**validation))
                            self.accepted_fit=None;self.pending=None;self.model_ready=False
                            self.model_fit.previous=None;self.model_fit.history=[]
                            self.model_fit.retired_information=None;self.model_fit.retired_vector=None
                    fit=self.model_fit.fit(obs,hint)
                    fit['time']=now
                    record={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in fit.items()}
                    acceptable=fit['success'] and fit['risk']<=math.log(1.2) and 200.<=fit['k']<=3e5
                    record.update(time=now,accepted_at_fit=bool(acceptable))
                    if validation is not None:record['previous_contact_validation']=validation
                    self.fit_records.append(record)
                    if acceptable:
                        if self.latency>0:self.pending=(now+self.latency,fit)
                        else:self.accepted_fit=fit
                except (ValueError,np.linalg.LinAlgError,FloatingPointError) as ex:
                    self.fit_records.append(dict(time=now,error=repr(ex)))
        if self.accepted_fit is not None and self.signals.rows[-1][4]>.02:
            fit=dict(self.accepted_fit)
            age=max(now-fit['time'],0.)/.128
            Q=np.diag([.12**2,.25**2,.04**2,.04**2] if self.model_fit.spatial else [.12**2])
            fit['cov']=fit['cov']+Q*age
            fit['data_cov']=fit['data_cov']+Q*age
            signal=self.signals.rows[-1]
            pred=self.model_fit.predict(fit,signal[1],signal[4])
            if pred['risk']<=math.log(1.2) and 200.<=pred['k']<=3e5:
                self.k_good=pred['k'];self.kg_rel=pred['relative_sd']
                self.model_current=True;self.model_ready=True;self.model_ticks+=1
                self.current_fit=fit;self.current_prediction=pred
        self.contact_ticks+=1
        self.contact_trace.append([now,self.baseline_good,self.k_good,int(self.model_current)])
        return feedback

    def drift_pred(self,o,n):
        if not self.model_current or len(self.signals.rows)<2 or (self.interface!='target' and not self.propagate_contact):
            return super().drift_pred(o,n)
        if self._motion_n==n:return self._motion_prediction
        current,previous=self.signals.rows[-1],self.signals.rows[-2]
        if self.propagate_contact and self.uncertainty_domain=='level':
            if self.observation_fit_time!=self.current_fit['time'] or self.last_contact_motion_n!=n-1:
                self.observation_fit_time=self.current_fit['time'];self.observation_anchor=previous.copy()
                self.cumulative_gradient=np.zeros(len(self.current_fit['theta']))
                if self.interface=='target':
                    self.anchor_gradient=self.model_fit.correction_gradient(self.current_fit,previous[1],previous[3],previous[4])[1]
        parameter_variance=0.
        if self.interface=='target':
            if self.propagate_contact:
                increment=self.model_fit.correction_increment(self.current_fit,current,previous)
                if self.uncertainty_domain=='level':
                    level=self.model_fit.functional_uncertainty(self.current_fit,
                                                               increment['gradient_current']-self.anchor_gradient)
                    increment.update(level)
                if increment['unobserved']:
                    self.unobserved_corrections+=1;self.model_ticks-=1;self.model_current=False
                    self.k_good=self.baseline_good
                    return super().drift_pred(o,n)
                dc=increment['increment'];parameter_variance=increment['variance']+increment['bias']**2
            else:
                c=self.model_fit.target_correction(self.accepted_fit,current[1],current[3],current[4])
                p=self.model_fit.target_correction(self.accepted_fit,previous[1],previous[3],previous[4])
                dc=c-p
        else:
            dc=(current[3]-previous[3])/self.k_good
            if self.uncertainty_domain=='level':
                self.cumulative_gradient-=dc*self.current_prediction['gain_gradient']
                stats=self.model_fit.functional_uncertainty(self.current_fit,self.cumulative_gradient)
                if stats['unobserved']:
                    self.unobserved_corrections+=1;self.model_ticks-=1;self.model_current=False
                    self.k_good=self.baseline_good
                    return super().drift_pred(o,n)
                parameter_variance=stats['variance']+stats['bias']**2
            else:
                parameter_variance=dc**2*(self.current_prediction['relative_sd']**2+
                                         self.current_prediction['regularization_log_shift']**2)
        self.last_contact_motion_n=n
        oldn=self._motion_n;self._motion_n=n
        if self._fin_prev is None or oldn!=n-1:
            self.Z=0.;self.tracker.reset(0.);self.q_prev=None
        else:
            # Both corrections use the SAME current model, avoiding artificial
            # level jumps when parameters are updated.
            self.Z+=dc-(current[2]-previous[2])
            force_variance=(self.params['noise']/self.k_good)**2*self.motion_params['r_scale']
            self.tracker.R=max(force_variance+parameter_variance,1e-20)
            self.observation_variances.append([parameter_variance,force_variance])
            self.tracker.update(self.Z)
        self._fin_prev,self._y_prev=current[3],current[2]
        L=(self.ndel+.5)*self.Tc;s,sd=self.tracker.predict(L)
        self.v_cmd_motion=-float(self.tracker.x[1])
        self._motion_prediction=(None,s+self.tau*sd)
        return self._motion_prediction
