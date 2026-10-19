"""Information causality, implementation equivalence and physical branch ranking."""
from pathlib import Path
import hashlib,json
import numpy as np
import controller as c
import model as m
import coupled_contact as cm
from run import noises
from pilot import context

HERE=Path(__file__).resolve().parent


def true_branch(p,noise,ctx,theta,scale,horizon,latency):
    y,zr,cfg,hist,i,resp,sense,temp,held=ctx
    yy=y.copy();sense=sense.copy();old=resp.copy();n=int(round(horizon/.01))
    history=np.zeros((i+n+1,8));history[:i+1]=hist;total=0.
    for step in range(n):
        t=(i+step)*.01
        if step>0:
            h=i+step;f=cm.contact(t,yy,p)[0]
            history[h,:4]=yy[:4];history[h,4:7]=f+noise[h];history[h,7]=p[9]*t+yy[0]
            rr,temp,_=m.response(history,h,old);old=rr.copy()
            sense[:2]=cfg[6]*rr[3:5]+(1-cfg[6])*sense[:2];rr[3:5]=sense[:2]
        else:rr=resp.copy()
        if step%4==0:target=c.action(yy,zr,cfg,rr,temp,theta,scale)
        # Only the first policy selection is expensive. Subsequent continuation
        # is the cheap base feedback, as in the predictive rollout.
        if step==0 and latency>0:
            first=min(.01,latency)
            yy,a=cm.advance(t,yy,p,held[:4],held[4],zr,sense,cfg,first)
            if first<.01:
                yy,b=cm.advance(t+first,yy,p,target[:4],target[4],zr,sense,cfg,.01-first);a[0]+=b[0]
        else:yy,a=cm.advance(t,yy,p,target[:4],target[4],zr,sense,cfg,.01)
        total+=a[0]
    return total


def main():
    fit=json.loads((HERE/'pilot.json').read_text());cfg=np.array(fit['cfg']);bank=np.array(fit['bank']);scales=np.array(fit['scales'])
    ps=m.sample(6,8201);ns=noises(ps,8202)
    equivalence=[];causality=[];ranking=[]
    for case in range(3):
        p=ps[case];n=ns[case]
        a,ta=c.episode(p,n,cfg,bank,scales,policy_index=0,latency=.000003)
        b,tb=cm.episode(p,n,12,cfg,bank[0],latency=.000003)
        error=float(np.max(np.abs(a-b)));assert error<1e-10
        equivalence.append({'case':case,'max_metric_difference':error})
        ctx=context(p,n,cfg,bank[0],150)
        y,zr,cc,hist,i,resp,sense,temp,held=ctx
        first=c.decide(y,zr,cc,hist,i,resp,sense,temp,bank,scales,.16,True,True,.006,held)
        corrupted=hist.copy();corrupted[max(1,i-1):i+1,4:7]=1e7
        second=c.decide(y,zr,cc,corrupted,i,resp,sense,temp,bank,scales,.16,True,True,.006,held)
        diff=float(np.max(np.abs(first[2]-second[2])));assert diff<1e-12
        causality.append({'case':case,'undelivered_wrench_score_difference':diff})
        nominal=c.decide(y,zr,cc,hist,i,resp,sense,temp,bank,scales,.16,False,True,.006,held)
        actual=np.array([true_branch(p,n,ctx,theta,scale,.16,.006) for theta,scale in zip(bank,scales)])
        ranking.append({'case':case,'family':case%3,'actual_costs':actual.tolist(),
                        'ensemble_predicted_costs':first[2].tolist(),'nominal_predicted_costs':nominal[2].tolist(),
                        'ensemble_choice':int(first[1]),'nominal_choice':int(nominal[1]),
                        'physical_best_policy':int(np.argmin(actual)),
                        'ensemble_absolute_regret':float(actual[first[1]]-actual.min()),
                        'nominal_absolute_regret':float(actual[nominal[1]]-actual.min())})
        print('audit case',case,ranking[-1],flush=True)
    out={'single_policy_equivalence':equivalence,'causality':causality,'physical_branch_ranking':ranking,
         'scope':'Three validation states; true future used ONLY for offline diagnostic labels, never for online decisions. Ranking regret is within this finite policy bank/horizon, not globally optimal regret.',
         'source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['controller.py','pilot.py','audit.py','pilot.json']}}
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
