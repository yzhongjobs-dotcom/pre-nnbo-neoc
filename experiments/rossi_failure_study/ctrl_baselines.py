"""Independent full controllers; no Method FF, tracker, H or uncertainty gate."""
import math
import numpy as np

import sweeps  # existing import-path setup
import common
import method3 as M
from comparators import RossiType, DOMPAC


class RossiController(RossiType):
    """Reuse audited Rossi estimator and Eq.(11); scene force replaces old 4 N."""
    def step(self,o):
        row=o.log[o.n-1];t,F,sent=(float(row[i]) for i in (0,5,11))
        self.k_hat=self.estimator.update(F,self.displacement.update(o))
        exc=self.A_E*math.sin(2*math.pi*self.f_exc*t)
        dy=self.K_I/self.k_hat*(self.force_set-F)*self.Tc
        self.exc_prev_delta=exc-self.exc_prev;self.exc_prev=exc
        self.trace.append((t,self.k_hat))
        return sent+dy+self.exc_prev_delta


class DOMPACController(DOMPAC):
    """Shen (17)-(29), outward x=-y. Same exact ESO/MPC reset as existing code."""
    def step(self,o):
        row=o.log[o.n-1];t,F,sent=(float(row[i]) for i in (0,5,11))
        e=F-self.force_set
        if self.sent_prev is not None:
            realized=sent-self.sent_prev
            if abs(realized-self.dy_prev)>1e-12:
                self.velocity=realized/self.Tc
        if self.z is None:
            self.z=np.array([e,0.,0.])
        else:
            self.z=self.eso_A@self.z+self.eso_B@np.array([e,self.held_input])
        # Paper says X=[ef,ef_dot], supplies ef directly and ef_dot from ESO.
        u0=-float(self.mpc_gain@np.array([e,self.z[1]]))
        b0=-self.k_h0/self.m; hn=b0*(e+self.b*self.velocity)
        u=(u0-hn-self.z[2])/b0
        vinf=(-e-u)/self.b
        dy=vinf*self.Tc+(self.velocity-vinf)*self.m/self.b*(1-self.decay)
        vnext=vinf+(self.velocity-vinf)*self.decay
        self.held_input=self.k_h0*(vnext-self.velocity)/self.Tc
        self.velocity=vnext;self.sent_prev=sent;self.dy_prev=dy
        self.trace.append((t,self.k_h0))
        return sent+dy


class XingController(common.DirectPolicy):
    """TASE2026 p6 (21)-(24), single force subspace, exact held-input m,b.

    Use the DISCRETE printed (22). The sign in continuous approximation (32)
    disagrees with (22); it is not used. No extra compensation clipping.
    """
    def __init__(self,m=1.,b=100.,alpha=.99,k1=1.,k2=1.,beta_bar=1.,diff_tau=.04):
        self.params=dict(m=m,b=b,alpha=alpha,k1=k1,k2=k2,beta_bar=beta_bar,diff_tau=diff_tau)
        self.__dict__.update(self.params)
        if not (0<alpha<=1 and min(m,b,k1,k2,beta_bar)>0):raise ValueError(self.params)

    def reset(self,platform):
        super().reset(platform)
        self.comp=0.;self.e_prev=None;self.velocity=0.;self.edot=0.;self.sent_prev=None
        self.decay=math.exp(-self.b*self.Tc/self.m)

    def step(self,o):
        row=o.log[o.n-1];e=float(row[5])-self.force_set;sent=float(row[11])
        if self.e_prev is not None:
            self.edot+=(1-math.exp(-self.Tc/self.diff_tau))*((e-self.e_prev)/self.Tc-self.edot)
            beta=1/(math.exp(min(self.k1*abs(e),700))+math.exp(min(self.k2*abs(self.edot),700))+self.beta_bar)
            self.comp=self.alpha*self.comp-beta*self.e_prev
        if self.sent_prev is not None and abs(sent-self.sent_prev-self.dy_prev)>1e-12:
            self.velocity=(sent-self.sent_prev)/self.Tc
        vinf=(-e+self.comp)/self.b
        dy=vinf*self.Tc+(self.velocity-vinf)*self.m/self.b*(1-self.decay)
        self.velocity=vinf+(self.velocity-vinf)*self.decay
        self.e_prev=e;self.sent_prev=sent;self.dy_prev=dy
        return sent+dy


def build(method,cfg,sc):
    cls={'Rossi-control':RossiController,'DO-MPAC':DOMPACController,'Xing-DAHIC':XingController}[method]
    return M.attach(cls(**cfg),sc)
