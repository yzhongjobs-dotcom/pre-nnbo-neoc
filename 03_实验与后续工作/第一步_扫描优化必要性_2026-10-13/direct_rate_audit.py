"""Offline diagnostic using the original gain-rate action, without target servo.

Piecewise constant rates connect bounded optimized gain endpoints. Full future
material is supplied only here to look for benefits missed by causal short MPC.
"""
import json,time
import numpy as np
from scipy.optimize import minimize
import scan_validation as v

def rate_step(y,t,h,speed,u,kind):
    def flow(y,t):
        # Invert the target interface so dynamics integrates exactly the proposed
        # constant gain rate. Common tank protection remains active.
        return v.dynamics(y,t,speed,y[:,4:8]+.04*u,None,kind)
    a,l,p,_=flow(y,t);b,lb,pb,_=flow(y+h*a/2,t+h/2)
    c,lc,pc,_=flow(y+h*b/2,t+h/2);d,ld,pd,_=flow(y+h*c,t+h)
    return y+h*(a+2*b+2*c+d)/6,h*(l+2*lb+2*lc+ld)/6,h*(p+2*pb+2*pc+pd)/6

def evaluate(q,kind,speed,gains,dt=.008,record=False):
    q=np.atleast_2d(q);blocks=4;duration=v.CFG['scan_length']/speed
    per=round(duration/blocks/dt);dt=duration/(blocks*per)
    y=np.array([v.initial(speed,gains) for _ in q]);total=np.zeros(len(q));rows=[]
    endpoints=v.LOW+q.reshape(len(q),blocks,4)*v.SPAN
    for block in range(blocks):
        u=(endpoints[:,block]-y[:,4:8])/(duration/blocks)
        for j in range(per):
            t=(block*per+j)*dt
            if record:
                fn,fr,e,ev,n,tan,delta,vt=v.geometry(y,t,speed,None,kind)
                rows.append([t,*y[0],fn[0],*e[0],delta[0]])
            y,l,p=rate_step(y,t,dt,speed,u,kind);total+=l+p
    if not record:return total
    a=np.asarray(rows)
    metric=dict(cost=float(total[0]),force_rmse_N=float(np.sqrt(np.mean((a[:,10]-.14)**2))),
        track_rmse_m=float(np.sqrt(np.mean(a[:,11]**2))),normal_rmse_m=float(np.sqrt(np.mean(a[:,12]**2))),
        force_peak_N=float(a[:,10].max()),min_energy_J=float(a[:,9].min()),
        min_contact_m=float(a[:,13].min()),gain_constraints_pass=bool(np.all((a[:,5:9]>=v.LOW-1e-8)&(a[:,5:9]<=v.HIGH+1e-8))))
    return metric,a

def main():
    gains=np.array(json.loads((v.HERE/'tuning.json').read_text(encoding='utf8'))['selected']['fixed']['gains'])
    result={}
    for kind,speed in [('ramp',.006),('cycle',.003)]:
        q=np.tile((gains-v.LOW)/v.SPAN,4)
        def obj(q):
            ep=2e-5;eye=np.eye(len(q))*ep
            values=evaluate(np.vstack((q,q+eye,q-eye)),kind,speed,gains)
            return values[0],(values[1:1+len(q)]-values[1+len(q):])/(2*ep)
        tic=time.perf_counter()
        res=minimize(obj,q,jac=True,method='L-BFGS-B',bounds=[(0,1)]*len(q),options=dict(maxiter=60,ftol=1e-9,gtol=1e-5))
        metric,a=evaluate(res.x,kind,speed,gains,dt=.002,record=True)
        fine,b=evaluate(res.x,kind,speed,gains,dt=.001,record=True)
        name=kind+('_fast' if speed==.006 else '_slow')
        np.savetxt(v.HERE/(name+'_offline_direct_rate.csv'),a,delimiter=',',header='t,x,z,vx,vz,Kt,Kn,Dt,Dn,tank,Fn,ex,ez,delta',comments='')
        result[name]=dict(success=bool(res.success),message=str(res.message),iterations=int(res.nit),
            metrics=metric,fine_metrics=fine,endpoints=(v.LOW+res.x.reshape(4,4)*v.SPAN).tolist(),
            seconds=time.perf_counter()-tic,information='offline full true future material; original constant gain-rate actions')
        (v.HERE/'direct_rate_oracle_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8')
        print(name,json.dumps(result[name]),flush=True)

if __name__=='__main__':main()
