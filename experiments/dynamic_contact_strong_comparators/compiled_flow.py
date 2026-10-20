"""Same equations/quadrature as the reference; JIT is only implementation work."""
import numpy as np
from numba import njit


@njit(cache=True)
def flow(y,z,th,p):
    e,v,k,d,E=y; uk=z[0]*p[24];ud=z[1]*p[25]
    fr=-k*e-d*v
    dy=np.array([v,(fr-th*e-p[1]*v)/p[0],uk,ud,
                 p[3]*(p[1]+d)*v*v-.5*e*e*uk])
    ell=p[17]*(e/p[9])**2+p[18]*(v/p[10])**2+p[19]*(th*e/p[11])**2
    ell+=p[20]*(fr/p[12])**2+p[21]*(((k-p[13])/p[15])**2+((d-p[14])/p[16])**2)
    ell+=p[22]*(z[0]*z[0]+z[1]*z[1])
    return dy,ell


@njit(cache=True)
def step(y,z,th,p,dt):
    a,la=flow(y,z,th,p);ya=y+dt*a/2
    b,lb=flow(ya,z,th,p);yb=y+dt*b/2
    c,lc=flow(yb,z,th,p);yc=y+dt*c
    d,ld=flow(yc,z,th,p)
    return y+dt*(a+2*b+2*c+d)/6,dt*(la+2*lb+2*lc+ld)/6,ya,yb,yc


@njit(cache=True)
def margins(y,th,p):
    e,v,k,d,E=y;fr=-k*e-d*v;fn=p[2]-th*e
    return np.array([(E-p[4])/p[5],(p[6]**2-fr*fr)/p[6]**2,
                     (fn-p[7])/p[2],(p[8]-fn)/p[2]])


@njit(cache=True)
def rollout(z,th,p,y0,n,sub,dt):
    batch=z.shape[0];steps=n*sub
    costs=np.zeros(batch,dtype=z.dtype)
    cons=np.empty((batch,(4*steps+1)*4),dtype=z.dtype)
    trace=np.empty((steps+1,8),dtype=z.dtype)
    for b in range(batch):
        y=y0.astype(z.dtype);idx=0
        if b==0:
            trace[0,0]=0.;trace[0,1:6]=y;trace[0,6:8]=0.
        for j in range(n):
            u=z[b,2*j:2*j+2]
            for s in range(sub):
                yn,ell,ya,yb,yc=step(y,u,th,p,dt)
                for stage in (y,ya,yb,yc):
                    cons[b,idx:idx+4]=margins(stage,th,p);idx+=4
                y=yn;costs[b]+=ell
                if b==0:
                    t=j*sub+s+1;trace[t,0]=t*dt;trace[t,1:6]=y
                    trace[t,6]=u[0]*p[24];trace[t,7]=u[1]*p[25]
        cons[b,idx:idx+4]=margins(y,th,p)
        costs[b]+=p[23]*((y[0]/p[9])**2+p[19]*(th*y[0]/p[11])**2)
    return costs,cons,trace
