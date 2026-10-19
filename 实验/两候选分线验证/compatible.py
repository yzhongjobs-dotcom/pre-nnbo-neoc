"""Candidate 1 iteration: observationally compatible contact-model decisions.

Empirical local model ensemble, NOT a certified posterior or dual controller.
Only already released histories are fit; no future information planning.
Risk-neutral weighted task cost, actual gain-rate dynamics unchanged.
"""
import time
import numpy as np
from scipy.optimize import least_squares,minimize
from optimization import context
from core import residuals
from rate_core import decode_rate


class CompatibleMPC:
    def __init__(self,problem,best_only=False):self.problem=problem;self.a=None;self.log=[];self.best_only=best_only
    def decide(self,t,y,zr,sensed,cfg,front,i):
        tic=time.perf_counter();j=i-2;q=front.q.copy();an=front.anchor.copy();models=[q];errors=[]
        if j<3:
            weights=np.ones(1);fit_ms=0.
        else:
            start=max(0,j-30);states=np.array(front.states[start:j+1]);times=np.array(front.times[start:j+1]);forces=np.array(front.forces[start:j+1])
            def data_error(qq):return residuals(qq,states,times,forces,an,cfg[7],q,np.ones(7))[:3*len(times)]
            errors.append(float(data_error(q)@data_error(q)))
            low=np.array([-2.,-20.,-3.,-10.,np.log(.001),.02]);high=np.array([2.,20.,3.,10.,np.log(.04),.8])
            # Broad physical alternatives, declared before validation; likelihood
            # fit can reject them. These are not newly observed true parameters.
            sd=np.array([.5,2.,1.,2.,1.,.2])
            for factor in (.5,2.):
                fixed=np.clip(q[0]+np.log(factor),np.log(.8e6),np.log(16.5e6))
                def fun(v):
                    qq=np.r_[fixed,v];return np.r_[data_error(qq),(v-q[1:])/sd]
                fit=least_squares(fun,np.clip(q[1:],low+1e-7,high-1e-7),bounds=(low,high),max_nfev=10,loss='soft_l1',f_scale=2.,diff_step=1e-4)
                qq=np.r_[fixed,fit.x];models.append(qq);rr=data_error(qq);errors.append(float(rr@rr))
            variance=max(1.,min(errors)/max(1,3*len(times)-7));logw=-.5*(np.array(errors)-min(errors))/variance
            weights=np.exp(np.maximum(-50.,logw));weights/=sum(weights);fit_ms=1000*(time.perf_counter()-tic)
            if self.best_only:
                weights=np.zeros(len(errors));weights[int(np.argmin(errors))]=1.
        contexts=[context(t,y,zr,qq,an,sensed,cfg) for qq in models]
        def fun(a):
            vv=0.;gg=np.zeros(10)
            for w,c in zip(weights,contexts):
                if w<1e-5:continue
                v,g=self.problem.fun(a,c);vv+=w*v;gg+=w*g
            return vv,gg
        init=np.zeros(10) if self.a is None else np.r_[self.a[5:],self.a[5:]]
        opt=minimize(fun,init,jac=True,bounds=[(-1,1)]*10,method='L-BFGS-B',options={'maxiter':12,'maxls':8,'ftol':1e-6})
        self.a=opt.x.copy();self.log.append({'seconds':time.perf_counter()-tic,'fit_ms':fit_ms,'weights':weights.tolist(),'fit_squared_residuals':errors,'predicted_cost':float(opt.fun),'iterations':int(opt.nit)})
        return decode_rate(opt.x[:5])
