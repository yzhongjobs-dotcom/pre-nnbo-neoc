"""Headroom gate (README section 'Headroom gate'), run before any candidate method.

For each dose cell: B1 best single fixed (tuned on train), B0 fixed extremes, B2 per-scene
grid oracle, B3b truth-stiffness schedule at information delays 0 / 0.1 / 0.2 s (rule
parameters tuned on train per delay). Pass: on the test cohort B1 fails at least one physical
criterion (peak > FD+3 N in >= 20 % of scenes, mean RMSE > 0.5 N, or time in +-0.5 N < 70 %)
while B3b at 0.1 s meets all of them, with the paired J difference CI excluding zero.
"""
import json, sys, time, itertools, numpy as np
import plant as P, scenarios as S, tune as U
from baselines_ref import TruthK

CELLS=[('T3',dict(ramp=0.25,A=0.001)),('T3',dict(ramp=0.25,A=0.004)),
       ('T3',dict(ramp=0.5,A=0.001)),('T3',dict(ramp=0.5,A=0.004)),
       ('T1',dict(ratio=30.,A=0.001)),('T1',dict(ratio=30.,A=0.004))]
RULE=list(itertools.product([20.,40.,80.,160.],[0.5,1.,2.,4.],[3.,10.,30.],[20.,100.]))
NTR,NTE=20,40

def fails(s):
    return dict(peak=s['frac_peak_over']>=0.2,rmse=s['rmse']>0.5,inband=s['inband']<0.7)

def run_cell(task,dose,plat):
    train=S.cohort(task,NTR,S.SPLIT['train'],dose); test=S.cohort(task,NTE,S.SPLIT['test'],dose)
    out={'task':task,'dose':dose,'platform':plat}
    t0=time.time(); b1=U.best_fixed(train,plat,seed0=S.SPLIT['train'])
    g1=b1['gains']; out['B1']={'gains':g1,'train_J':b1['train_J'],'grid_best':b1['grid_best'],'any_feasible_on_grid':b1['any_feasible_on_grid']}
    ext=[(0.3,20.,0.,10.,0.),(0.3,2500.,0.,10.,0.)]
    res=U.evaluate(test,[g1]+ext,plat,S.SPLIT['test'])
    out['B1']['test']=U.summary([r[0] for r in res]); out['B0_low']=U.summary([r[1] for r in res]); out['B0_high']=U.summary([r[2] for r in res])
    jB1=[r[0]['J'] for r in res]
    # B2 per-scene grid oracle (test)
    g=U.evaluate(test,U.FIXED_GRID,plat,S.SPLIT['test']); Jg=U.J(g); a=Jg.argmin(1)
    out['B2']=U.summary([g[i][a[i]] for i in range(len(test))])
    # B3b truth-k schedule
    out['B3b']={}
    for lat in [0.,0.1,0.2]:
        ctr=[TruthK(lat,wI,cP,ki,Dm) for (wI,cP,ki,Dm) in RULE]
        rt=U.evaluate(train,ctr,plat,S.SPLIT['train']); Jt=U.J(rt); feas=U.feasible(rt)
        j=int(np.argmin(np.where(feas,Jt.mean(0),Jt.mean(0)+1e6)))
        rs=U.evaluate(test,[ctr[j]],plat,S.SPLIT['test']); rows=[r[0] for r in rs]
        out['B3b'][str(lat)]={'rule':RULE[j],'test':U.summary(rows),'dJ_vs_B1':U.paired_ci([r['J'] for r in rows],jB1)}
    s1=out['B1']['test']; s3=out['B3b']['0.1']['test']; ci=out['B3b']['0.1']['dJ_vs_B1']
    out['gate']={'B1_fails':fails(s1),'B3b01_fails':fails(s3),'ci_excludes_zero':ci[2]<0,
                 'pass':any(fails(s1).values()) and not any(fails(s3).values()) and ci[2]<0}
    out['seconds']=time.time()-t0
    return out

if __name__=='__main__':
    plat=sys.argv[1] if len(sys.argv)>1 else 'PC125'
    sel=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else range(len(CELLS))
    allres=[]
    for c in sel:
        task,dose=CELLS[c]; r=run_cell(task,dose,plat); allres.append(r)
        b1=r['B1']['test']; b2=r['B2']; b3=r['B3b']['0.1']['test']
        print(f"{task} {dose} {plat}: B1 {tuple(round(x,6) for x in r['B1']['gains'])} rmse {b1['rmse']:.3f} peak {b1['peak']:.2f} over {b1['frac_peak_over']:.2f} inband {b1['inband']:.2f} | "
              f"B2 rmse {b2['rmse']:.3f} peak {b2['peak']:.2f} | B3b@0.1 rmse {b3['rmse']:.3f} peak {b3['peak']:.2f} over {b3['frac_peak_over']:.2f} inband {b3['inband']:.2f} | gate {r['gate']['pass']} ({r['seconds']:.0f}s)",flush=True)
        json.dump(allres,open(f'gate_{plat}.json','w'),indent=1)
