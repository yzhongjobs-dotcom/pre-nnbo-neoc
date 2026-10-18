"""Charge measured planning time to physical response on the same 8 cases.

Hold the old actual gains while computing; execute the returned rate plan after
release, without pretending it was available at t=0. This diagnoses stale-plan
deployment, not a full RTI/feedback policy or an unavoidable lower bound.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
import transient as t
from verify import independent

HERE=Path(__file__).resolve().parent


def main():
    source=json.loads((HERE/'transient_results.json').read_text())
    bounded=json.loads((HERE/'bounded_results.json').read_text())['cases']
    raw=np.load(HERE/'transient_paired.npz');records=[]
    for row,policy in zip(source['cases'],bounded):
        label,case=row['family'],row['case'];colored=label!='white';dim=4 if colored else 3;n=10 if colored else 4
        arr=raw[label+f'_{case}_initial'];g0=arr[2*dim:2*dim+4]
        mu0=arr[2*dim+4:2*dim+4+n];p0=arr[-n*n:].reshape(n,n)
        rates=np.array(policy['rate_blocks']).reshape(2,4)*t.RATE
        def simulate(delay,hold=False):
            state=np.r_[mu0,p0.ravel(),0.]
            def rhs(tt,y):
                elapsed=max(0.,tt-delay)
                g=g0 if hold else g0+rates[0]*min(elapsed,t.HORIZON/2)+rates[1]*max(0.,elapsed-t.HORIZON/2)
                _,a,c,q,_=independent(row['new_features'],g,colored,colored,False)
                mu=y[:n];cov=y[n:n+n*n].reshape(n,n)
                return np.r_[a@mu,(a@cov+cov@a.T+q).ravel(),mu@c@mu+np.trace(c@cov)]
            breaks=sorted(set([0.,t.HORIZON,min(delay,t.HORIZON),min(delay+t.HORIZON/2,t.HORIZON)]))
            for left,right in zip(breaks[:-1],breaks[1:]):
                if right<=left:continue
                sol=solve_ivp(rhs,(left,right),state,method='DOP853',rtol=1e-10,atol=1e-13)
                state=sol.y[:,-1]
            return float(state[-1])
        baseline=simulate(0,True);measured=policy['history'][-1]['elapsed_ms']/1000
        results=[]
        for name,delay in [('instantaneous',0.),('5ms',.005),('20ms',.020),('50ms',.050),('100ms',.100),('measured_plan_time',measured)]:
            value=simulate(delay)
            results.append(dict(release=name,delay_s=delay,J=value,improvement_over_hold=1-value/baseline))
        records.append(dict(family=label,case=case,hold_J=baseline,releases=results))
    out=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        bounded_sha256=hashlib.sha256((HERE/'bounded_iterations.py').read_bytes()).hexdigest(),
        protocol=dict(task='Unchanged 8 transient instances, independent DOP853 physical evaluation',
            during_compute='Old actual gains continue feedback; gain-rate command is zero',
            after_release='Same previously computed two-block rate plan, remaining horizon only',
            limitations='Not delay-aware replanning, not a learned state-feedback policy; measured times are prototype observations, not WCET',
            status='Post-screen deployment diagnostic; no model/cost/test-case changes'),cases=records)
    (HERE/'delay_results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for label in ['white','colored_observation_noise']:
        rr=[r for r in records if r['family']==label]
        print(label,{name:float(np.mean([next(x for x in r['releases'] if x['release']==name)['improvement_over_hold'] for r in rr])) for name in ['instantaneous','5ms','20ms','50ms','100ms','measured_plan_time']})


if __name__=='__main__':main()
