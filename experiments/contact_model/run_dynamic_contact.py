"""Dynamic-consistent contact: mechanism checks and measured-input causal replay.

All memory hypotheses are declared before comparison. Values are neither plant
parameters nor estimates of them. The original plant/controllers remain intact.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from dynamic_contact import Memory,ContactFit,Signals,project

HERE=Path(__file__).resolve().parent
CONFIGS={
 'static_linear':dict(memory=Memory(0.,.2),spatial=False),
 'memory_linear':dict(memory=Memory(.5,.2),spatial=False),
 'memory_patch':dict(memory=Memory(.5,.2),spatial=True),
 'memory_patch_low':dict(memory=Memory(.2,.08),spatial=True),
 'memory_patch_high':dict(memory=Memory(.8,.6),spatial=True),
 'memory_line':dict(memory=Memory(.5,.2),spatial=True,contact='line')}


def verify():
    rng=np.random.default_rng(84112);elastic=4.+.1*rng.standard_normal(150)
    checks={}
    for alpha,tau in [(0.,.2),(.2,.08),(.5,.2),(.8,.6)]:
        m=Memory(alpha,tau);F=m.forward(elastic,4.,4.)
        recovered=m.inverse(F,4.,4.)
        err=float(np.max(np.abs(recovered-elastic)));assert err<1e-10
        poles=(m.r,(m.r+alpha)/(1+alpha));assert max(poles)<1
        checks[f'{alpha},{tau}']=dict(inverse_error=err,poles=poles)
    # Independently solve one instantaneous force correction, allowing the
    # relaxation spring to change with the elastic force.
    model=ContactFit(memory=Memory(.5,.2),arrival=False)
    fit=dict(theta=np.array([math.log(12000.)]))
    F,Fe=4.7,4.2;delta=model.target_correction(fit,0.,F,Fe)
    expected=(F-4.)/(12000*1.5)
    assert abs(delta-expected)<1e-12
    checks['instantaneous_target_linear_limit']=dict(correction=delta,expected=expected)
    # Exact static limit and causal prefix.
    n=60;t=np.arange(n)*.008;y=.001+.00002*np.sin(2*np.pi*20.8*t)
    rows=np.c_[t,t*.02,y,4.+10000*(y-.001),y,np.zeros(n),np.zeros(n)]
    s1,s2=Signals(),Signals()
    changed=rows.copy();changed[40:]=1e5
    for i in range(40):s1.append(rows[i]);s2.append(changed[i])
    np.testing.assert_array_equal(s1.window()['fe'],s2.window()['fe'])
    checks['causal_prefix']=True
    return checks


def heldout_prediction(fitter,fit,signals,next_row):
    # Conditional new-encoder input; next force hidden until scoring.
    obs=signals.window();m=fitter.model;theta=fit['theta'];mem=fitter.memory
    beta=fit['surface_beta']
    # New encoder interpolation is available; no next Fn/Ft/My used here.
    prev=signals.previous;mm=np.arange(8);af=1-math.exp(-2*math.pi*signals.fc*.001)
    w=af*(1-af)**(7-mm);fr=-np.expm1(-mm*.001/signals.tau_servo)/(-math.expm1(-mem.dt/signals.tau_servo))
    q=float(w@fr/w.sum());xn=prev[1]+q*(next_row[1]-prev[1]);yn=prev[2]+q*(next_row[2]-prev[2])
    tn=(next_row[0]-prev[0])/(16*mem.dt);sn=np.array([tn**i for i in range(4)])@beta
    fen=float(m.response(xn,max(yn-sn,0.),theta,stiffness=False)[0])
    qprev=signals.last_force-signals.last_fe
    fn=fen+mem.r*qprev+mem.alpha*(fen-signals.last_fe)
    decay=math.exp(-2*math.pi*signals.fc*mem.dt)
    return decay*prev[3]+(1-decay)*fn


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify-only',action='store_true')
    ap.add_argument('--max-windows',type=int,default=8)
    ap.add_argument('--families',nargs='+',default=['S2','T3m','T2m-soft','T2m-stiff','T2f-stiff'])
    ap.add_argument('--models',nargs='+',default=['static_linear','memory_linear','memory_patch'])
    ap.add_argument('--label',default='dynamic_contact_replay_v1');args=ap.parse_args()
    checks=verify()
    if args.verify_only:print(json.dumps(checks,indent=2));return
    rows=[]
    for fam in args.families:
        files=list((HERE/'results').glob(f'finite_patch_input_{fam}_*.npz'))
        if len(files)!=1:raise ValueError(f'expected one frozen diagnosis input for {fam}')
        a=np.load(files[0]);z=a['measured'];sid=int(files[0].stem.split('_')[-1])
        ends=np.arange(128,len(z),16);ends=ends[z[ends-1,0]>=1.]
        if args.max_windows>0 and len(ends)>args.max_windows:ends=ends[np.linspace(0,len(ends)-1,args.max_windows).astype(int)]
        selected=set(ends.tolist())
        signals={name:Signals(CONFIGS[name]['memory']) for name in args.models}
        fits={name:ContactFit(**CONFIGS[name]) for name in args.models}
        for j,row in enumerate(z):
            end=j+1
            for s in signals.values():s.tau_servo=float(a['baseline_tau'][j]);s.append(row)
            if end not in selected:continue
            hint=float(a['baseline_kraw'][j]);hint=hint if np.isfinite(hint) and hint>0 else 3e4
            for name in args.models:
                obs=signals[name].window()
                if obs is None:
                    rows.append(dict(family=fam,sid=sid,end=end,model=name,unavailable=True));continue
                r=fits[name].fit(obs,hint)
                out={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in r.items()}
                if end<len(z):out['conditional_force_error']=heldout_prediction(fits[name],r,signals[name],z[end])-z[end,3]
                rows.append(dict(family=fam,sid=sid,end=end,model=name,scoring_k=float(a['scoring_truth_k'][j]),
                                 baseline_kgood=float(a['baseline_kgood'][j]),**out))
            print(f'{fam} end={end}',flush=True)
    out=dict(protocol=dict(window=16,models=args.models,max_windows=args.max_windows,
                           memory_hypotheses={n:vars(CONFIGS[n]['memory']) for n in args.models},
                           selected_parameters_from_truth=False,controller_changed=False,final_used=False),
             verification=checks,rows=rows)
    target=HERE/'results'/(args.label+'.json');target.write_text(json.dumps(out,indent=2),encoding='utf-8')
    for fam in args.families:
        for name in args.models:
            rr=[r for r in rows if r['family']==fam and r['model']==name and 'k'in r and r['scoring_k']>0]
            if not rr:continue
            e=[abs(math.log(r['k']/r['scoring_k'])) for r in rr]
            ef=[r['conditional_force_error'] for r in rr if 'conditional_force_error'in r]
            print(fam,name,'n',len(rr),'log median/P90',np.round(np.percentile(e,[50,90]),3),
                  'conditional force RMSE',round(float(np.sqrt(np.mean(np.square(ef)))),4),
                  'accepted',sum(r['risk']<=math.log(1.2) for r in rr))
    print(target)


if __name__=='__main__':main()
