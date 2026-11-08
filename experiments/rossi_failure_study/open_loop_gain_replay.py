import sys;sys.path.insert(0,'.')
import numpy as np, sweeps, scenes3 as S3, physics_plant3 as P3
F0=20.8
def fit(tt,v,f=F0):
    s=(tt-tt.mean());A=np.c_[np.cos(2*np.pi*f*tt),np.sin(2*np.pi*f*tt),np.ones(len(tt)),s,s**2,s**3]
    a=np.linalg.lstsq(A,v,rcond=None)[0];return complex(a[0],-a[1])
for i in (0,1,2):
    sc=S3.cell_scene3('rest','relaxed','nominal',i,base=611700,T=6.,overrides=dict(force=4.,slide='long',dynamics='neural_white',phase=14.))
    saved=np.load(f'results/structured_tvk_audit/rest_{i}.npz',allow_pickle=True)['log']
    cmd=saved[:750,11];nt=len(cmd);tc=np.arange(nt)*.008
    hold=np.convolve(cmd,np.ones(125)/125,'same');hold[:125]=hold[125];hold[-125:]=hold[-126]
    c=hold+60e-6/.37*np.sin(2*np.pi*F0*tc)
    pl=dict(P3.PLATFORM,log_dt=.001,tau_a=sc.get('tau_a_true',P3.PLATFORM['tau_a']))
    L=P3.replay(sc,c,seed=0,sensor=False,platform=pl)
    t,y,fax,k,fm=L[:,0],L[:,2],L[:,18],L[:,21],L[:,5]
    w=(t>1.5)&(t<5.5)
    # continuous gain from 1 ms samples, over one long stretch in 0.5 s pieces
    G=[];Ga=[];Gm=[];Gf=[]
    for a in np.arange(1.5,5,.5):
        m=(t>=a)&(t<a+.5);G.append(fit(t[m],fax[m])/fit(t[m],y[m]))
        idx=np.where(m)[0][::8]
        Ga.append(fit(t[idx],fax[idx])/fit(t[idx],y[idx]))
        Gm.append(fit(t[idx],fm[idx])/fit(t[idx],y[idx]))
        Gf.append(fit(t[m],fax[m],F0+125*0)/fit(t[m],y[m]))
    G=np.median(np.abs(G)),np.median(np.degrees(np.angle(G)));Ga=np.median(np.abs(Ga)),np.median(np.degrees(np.angle(Ga)));Gm=np.median(np.abs(Gm))
    im=[abs(fit(t[(t>=a)&(t<a+.5)],fax[(t>=a)&(t<a+.5)],125-F0))/abs(fit(t[(t>=a)&(t<a+.5)],fax[(t>=a)&(t<a+.5)])) for a in np.arange(1.5,5,.5)]
    f30=abs(1/(1+1j*F0/30.))
    print(f'rest_{i} truth {np.median(k[w]):.0f} | continuous F/y {G[0]:.0f}<{G[1]:.1f} | 8 ms sampled F/y {Ga[0]:.0f}<{Ga[1]:.1f} | sampled measured-force/y {Gm:.0f} (undo 30 Hz: {Gm/f30:.0f}) | 104 Hz image over 20.8 Hz in true force {np.median(im):.2f} | tool {sc["tool_mass"]:.2f} kr {sc["kr"]:.0f} tool*w^2/kr {sc["tool_mass"]*(2*np.pi*F0)**2/sc["kr"]:.3f}',flush=True)
