"""Oracle time-switched schedule at the (known) boundary for T1/T2: does
knowing where the material changes buy anything beyond the best fixed gains?"""
import numpy as np, itertools, json
from multiprocessing import Pool
import plant as P, scenarios as S
KZ=[10.,100.,400.]; DZ=[5.,20.,60.,150.]
def sched(p,before,after,lead):
    tb=(p[3]-0.42)/p[12]; te=(p[4]-0.42)/p[12] if p[4]<5 else 1e9
    def pol(t,o,g):
        return np.array(after if (tb-lead)<=t<te+0.3 else before)
    return pol
def run(args):
    task,i,p,mode=args; best=(1e9,None); base=None
    fixed=[3000.,10.,300.,20.]
    for kz1,dz1,kz2,dz2 in itertools.product(KZ,DZ,KZ,DZ):
        for lead in [0.0,0.3]:
            b=[3000.,kz1,300.,dz1]; a=[3000.,kz2,300.,dz2]
            sc=S.score(P.episode(p,sched(p,b,a,lead),mode=mode,ctrl=(1.5,10.),seed=i)[0])['score']
            if b==a and b==fixed and lead==0: base=sc
            if sc<best[0]: best=(sc,(b,a,lead))
    return task,i,mode,base,best
if __name__=='__main__':
    jobs=[(t,i,p,m) for t in ['T1','T2'] for i,p in enumerate(S.cohort(t,6,2026)) for m in [0,1]]
    with Pool(4) as pool: res=pool.map(run,jobs)
    out={}
    for t in ['T1','T2']:
        for m in [0,1]:
            r=[x for x in res if x[0]==t and x[2]==m]; base=np.array([x[3] for x in r]); best=np.array([x[4][0] for x in r])
            out[f'{t}_mode{m}']={'fixed_ref_mean':float(base.mean()),'switched_oracle_mean':float(best.mean()),'gain_percent':float(100*np.mean((base-best)/base)),'choices':[x[4][1] for x in r]}
            print(t,'PI' if m else 'noPI','fixed %.3f switched-oracle %.3f gain %.1f%%'%(base.mean(),best.mean(),out[f'{t}_mode{m}']['gain_percent']))
    json.dump(out,open('switched_oracle_results.json','w'),indent=1)
