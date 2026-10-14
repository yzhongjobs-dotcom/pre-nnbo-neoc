"""Independent carrier refinement audit; does not modify the running pilot."""
from pathlib import Path
import json
import time
import numpy as np
from continuous_contact import UniformCarrier,src,DREF,FD,T,force,MASS,DAMP,BR,RU

HERE=Path(__file__).resolve().parent


def refine_gains(c,spacing):
    if spacing is None:return c
    c.k=np.unique(np.r_[np.linspace(30.,300.,19),np.arange(120.,180.+spacing/2,spacing)])
    c.shape=(len(c.d),len(c.v),len(c.k))
    dd,vv,kk=np.meshgrid(c.d,c.v,c.k,indexing='ij')
    c.x=np.c_[dd.ravel(),vv.ravel(),kk.ravel()];c.S=len(c.x)
    c.a=src.halfwidth(c.x[:,0]);c.norm=(c.x[:,2]-150.)/150.
    c.base_force=.02/(2*.03)*30000*np.pi*c.a*c.a/2
    c.model_cache.clear()
    return c


def frozen_rollout(c,V,eta=1.2,dt=.001):
    x=np.array([DREF,0.,150.,.02]);cost=0.;rows=[]
    for i in range(round(T/dt)):
        t=i*dt;d,v,k=x[:3];f=force(d,eta);e=d-DREF
        phase=min(int(t/c.h),c.N-1)
        acc=(FD-f-k*e-(DAMP+BR)*v)/MASS
        geom=c.geometry(np.array([d+c.h*v]),np.array([v+c.h*acc]))
        fraction=min(max(t/c.h-phase,0.),1.)
        line=(1-fraction)*c.lines(V[phase+1],geom)+fraction*c.lines(V[min(phase+2,c.N)],geom)
        u=float(c.minimize(line,np.array([k]))[1][0])
        stage=(e/DREF)**2+(v/.008)**2+(f/FD-1)**2+.1*((k-150)/150)**2+RU*u*u
        rows.append([t,*x,f,u,stage])
        def rhs(y):
            ff=force(y[0],eta);ee=y[0]-DREF
            return np.array([y[1],(FD-ff-y[2]*ee-(DAMP+BR)*y[1])/MASS,
                             u,.3*(DAMP+BR)*y[1]**2-.5*ee**2*u])
        aa=rhs(x);bb=rhs(x+dt*aa/2);cc=rhs(x+dt*bb/2);dd=rhs(x+dt*cc)
        x+=dt*(aa+2*bb+2*cc+dd)/6;cost+=dt*stage
    total=cost+((x[0]-DREF)/DREF)**2+(force(x[0],eta)/FD-1)**2
    return float(total),np.array(rows)


def main():
    # Original pilot's printed independent audit is retained as its source.
    log=(HERE/'pilot_log.txt').read_text(encoding='utf-8-sig')
    ref=json.loads(next(line for line in log.splitlines() if line.startswith('{"shooting"')))
    reference=ref['shooting']['cost'];results=[]
    for nd,nv,spacing in [(33,17,None),(33,17,5.),(33,33,5.),(49,49,5.),(49,49,2.5)]:
        c=refine_gains(UniformCarrier(nd,nv,200),spacing)
        begin=time.perf_counter();V,_=c.solve((1.2,0.));elapsed=time.perf_counter()-begin
        val,trace=frozen_rollout(c,V)
        row=dict(nd=nd,nv=nv,gain_spacing=spacing,shape=c.shape,cost=val,
                 gap_to_shooting=val/reference-1,solve_seconds=elapsed)
        results.append(row);print(json.dumps(row),flush=True)
        (HERE/'resolution_audit.json').write_text(json.dumps(dict(reference=ref,grids=results),indent=2),encoding='utf-8')
        del c,V


if __name__=='__main__':main()
