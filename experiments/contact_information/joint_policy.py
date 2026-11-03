"""First-module candidate. Original feedback/IMM/recovery retained."""
import math
import sys
import time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'method_redesign'))
sys.path.insert(0,str(HERE.parent/'motion_aware_scan'))
from v1_servo_id import ServoID
from local_contact import fit


class JointPolicy(ServoID):
    def __init__(self,wave=None,new_fit=True,rate=True,constant=False,information_window=32,amplitude_scale=1.,warm=1.,replace_value='new',**kw):
        super().__init__(**kw);self.wave=None if wave is None else np.asarray(wave)
        self.new_fit=bool(new_fit);self.rate=bool(rate)
        self.constant=bool(constant)
        if replace_value not in ('new','short'):raise ValueError('replace_value')
        self.replace_value=replace_value
        self.information_window=int(information_window);self.amplitude_scale=float(amplitude_scale);self.warm=float(warm)

    def reset(self,platform):
        super().reset(platform);self.probe_previous=0.;self.baseline_gain=float('nan')
        self.contact_fits=[];self.probe_records=[];self.used=0

    def _stiffness_est(self,o,n):
        if not self.new_fit:return super()._stiffness_est(o,n)
        self.k_good=self.baseline_gain
        feedback=super()._stiffness_est(o,n);self.baseline_gain=self.k_good
        if o.log[n-1,0]<self.warm+.064+self.information_window*.008:return feedback
        tic=time.perf_counter()
        measured=o.log[:n,[0,1,2,5,11,12,13]]
        r=fit(measured,self.information_window,self.tau,self.fc,self.params['noise'],self.rate,self.constant)
        if r is None:return feedback
        accepted=bool(200<=r['k']<=3e5 and r['start_gain']>0 and r['risk']<=math.log(1.2)
                      and r['chi2']<=self.params['chi2_max'])
        self.contact_fits.append(dict(time=float(o.log[n-1,0]),accepted=accepted,
            elapsed=time.perf_counter()-tic,**{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in r.items()}))
        if accepted:
            self.k_good=r['k'] if self.replace_value=='new' else self.k_hat
            self.kg_rel=r['relative_sd'];self.used+=1
        return feedback

    def step(self,o):
        if self.wave is None and self.amplitude_scale==1:return super().step(o)
        p=self.params;enabled=p['excitation'];p['excitation']=False
        try:command=super().step(o)
        finally:p['excitation']=enabled
        F=float(o.log[o.n-1,5]);lost=p['recovery']=='gate' and F<.3
        active=enabled and self.mode==0 and not lost
        t=float(o.log[o.n-1,0]);old=math.sin(2*math.pi*p['f_exc']*(t+self.Tc))
        new=old if self.wave is None else self.wave[o.n%len(self.wave)]
        blend=float(np.clip((t-self.warm)/.064,0.,1.))
        value=p['a_force']*self.exc_gain*((1-blend)*old+blend*self.amplitude_scale*new) if active else 0.
        delta=value-self.probe_previous;self.probe_previous=value
        k=self.k_good if np.isfinite(self.k_good) else self.k_seek
        self.probe_records.append([float(o.log[o.n-1,0]),value,delta/k])
        return command+delta/k if self.mode==0 and not lost else command
