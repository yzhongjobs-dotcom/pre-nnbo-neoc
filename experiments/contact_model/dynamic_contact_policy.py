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
                 latency=0.,update_ticks=16,**kw):
        super().__init__(**kw)
        self.contact_args=dict(spatial=spatial,memory=memory,contact=contact)
        self.interface,self.latency,self.update_ticks=interface,float(latency),int(update_ticks)

    def reset(self,platform):
        super().reset(platform)
        self.model_fit=ContactFit(**self.contact_args)
        self.signals=Signals(self.contact_args['memory'],fc=self.fc,tau_servo=self.tau)
        self.accepted_fit=None;self.pending=None;self.model_ready=False
        self.fit_records=[];self.model_ticks=0;self.contact_ticks=0;self.model_current=False
        self.signal_n=None

    def step(self,o):
        n=o.n
        if self.signal_n!=n:
            self.signals.tau_servo=self.tau
            self.signals.append(o.log[n-1,[0,1,2,5,11,12,13]])
            self.signal_n=n
        if o.log[n-1,5]<.3:
            self.accepted_fit=None;self.pending=None;self.model_ready=False
            self.model_fit.previous=None
        return super().step(o)

    def _stiffness_est(self,o,n):
        previous_good=self.k_good
        feedback=super()._stiffness_est(o,n)
        # Original projection updates feedback only. Model measurement gain has
        # its own residual, variance and acceptance score, not the old gate.
        if self.model_ready:self.k_good=previous_good
        now=float(o.log[n-1,0]);self.model_current=False
        if self.pending is not None and now>=self.pending[0]:
            self.accepted_fit=self.pending[1];self.pending=None
        if n>=64 and n%self.update_ticks==0:
            obs=self.signals.window()
            if obs is not None:
                hint=self.k_hat if np.isfinite(self.k_hat) and self.k_hat>0 else self.kbar
                try:
                    fit=self.model_fit.fit(obs,hint)
                    fit['time']=now
                    record={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in fit.items()}
                    record.update(time=now,accepted_at_fit=fit['risk']<=math.log(1.2))
                    self.fit_records.append(record)
                    if fit['risk']<=math.log(1.2) and 200.<=fit['k']<=3e5:
                        if self.latency>0:self.pending=(now+self.latency,fit)
                        else:self.accepted_fit=fit
                except (ValueError,np.linalg.LinAlgError,FloatingPointError) as ex:
                    self.fit_records.append(dict(time=now,error=repr(ex)))
        if self.accepted_fit is not None and self.signals.rows[-1][4]>.02:
            fit=dict(self.accepted_fit)
            age=max(now-fit['time'],0.)/.128
            Q=np.diag([.12**2,.25**2,.04**2,.04**2] if self.model_fit.spatial else [.12**2])
            fit['cov']=fit['cov']+Q*age
            signal=self.signals.rows[-1]
            pred=self.model_fit.predict(fit,signal[1],signal[4])
            if pred['risk']<=math.log(1.2) and 200.<=pred['k']<=3e5:
                self.k_good=pred['k'];self.kg_rel=pred['relative_sd']
                self.model_current=True;self.model_ready=True;self.model_ticks+=1
        self.contact_ticks+=1
        return feedback

    def drift_pred(self,o,n):
        if self.interface!='target' or not self.model_current or len(self.signals.rows)<2:
            return super().drift_pred(o,n)
        if self._motion_n==n:return self._motion_prediction
        oldn=self._motion_n;self._motion_n=n
        current,previous=self.signals.rows[-1],self.signals.rows[-2]
        if self._fin_prev is None or oldn!=n-1:
            self.Z=0.;self.tracker.reset(0.);self.q_prev=None
        else:
            c=self.model_fit.target_correction(self.accepted_fit,current[1],current[3],current[4])
            p=self.model_fit.target_correction(self.accepted_fit,previous[1],previous[3],previous[4])
            # Both corrections use the SAME current model, avoiding artificial
            # level jumps when parameters are updated.
            self.Z+=(c-p)-(current[2]-previous[2])
            self.tracker.R=max((self.params['noise']/self.k_good)**2*self.motion_params['r_scale'],1e-20)
            self.tracker.update(self.Z)
        self._fin_prev,self._y_prev=current[3],current[2]
        L=(self.ndel+.5)*self.Tc;s,sd=self.tracker.predict(L)
        self.v_cmd_motion=-float(self.tracker.x[1])
        self._motion_prediction=(None,s+self.tau*sd)
        return self._motion_prediction
