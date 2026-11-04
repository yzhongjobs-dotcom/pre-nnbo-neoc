"""Scored probe-frequency gain against the measured transfer, across tissues. Truth check only."""
import math, numpy as np, json
from check_plant import PP, SC, Hold
rows=[]
for tissue in ('soft','mid','hard'):
    rng=np.random.default_rng(5); sc=SC.make('rest',tissue,rng,T=6.,v=0.,overrides=dict(field_sigma=0.)); sc['h_s'][:]=0.
    for f in (20.8,):
        st={}
        def cmd(o,st=st):
            if 'y' not in st: st['y']=float(o.log[0,11])
            return st['y']+1e-6*math.sin(2*np.pi*f*o.log[o.n-1,0])
        log=PP.episode(sc,Hold(cmd),3,sensor=False); m=log[:,0]>2.; t=log[m,0]
        X=np.c_[np.sin(2*np.pi*f*t),np.cos(2*np.pi*f*t),np.ones(len(t)),t]
        cy=np.linalg.lstsq(X,log[m,2],rcond=None)[0]; cf=np.linalg.lstsq(X,log[m,3],rcond=None)[0]; G=complex(cf[0],cf[1])/complex(cy[0],cy[1])
        rows.append(dict(tissue=tissue,measured=G.real,quad=G.imag,scored=float(np.median(log[m,8])),ratio=G.real/float(np.median(log[m,8]))))
print(json.dumps(rows,indent=1))
