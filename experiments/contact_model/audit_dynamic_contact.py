"""Measurement-only model diagnostics; reference stiffness enters scoring only.

The exact scalar linear profile is an independent check on the nonlinear
implementation and the common-scale assumption over multiple short windows.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from dynamic_contact import Memory, ContactFit, Signals, project

HERE = Path(__file__).resolve().parent


def linear_profile(obs, memory):
    T = memory.operator(len(obs['f']))
    u, _ = project(T@obs['B'], T@obs['y'])
    v, _ = project(T@obs['B'], T@obs['fe'])
    return u, v


def check_retirement():
    """Independent scalar information recursion verifies no double counting."""
    mem=Memory(); fitter=ContactFit(memory=mem,blocks=4);k=10000.;pieces=[]
    expected_arrival=0.;largest=0.
    for b in range(9):
        t=(np.arange(16)+16*b)*.008;tt=(t-t[-1])/.128
        B=np.column_stack([tt**j for j in range(4)])
        y=.002+.00008*np.sin(2*np.pi*20.8*t);surface=B@np.array([.001,.0001,-.00002,.00001])
        fe=k*(y-surface);x=.01*t
        obs=dict(x=x,y=y,fe=fe,f=mem.forward(fe,fe[0],fe[0]),B=B,
                 ft=np.zeros(16),my=np.zeros(16),xend=float(x[-1]),tend=float(t[-1]))
        u,_=linear_profile(obs,mem);pieces.append(k*k*(u@u)/.04**2)
        fit=fitter.fit(obs,k*1.5)
        assert abs(np.log(fit['k']/(1.5*k)))<2e-5
        if b>=4:
            expected_arrival=expected_arrival/(1+.12**2*expected_arrival)+pieces[b-4]
        actual=float(fitter.retired_information[0,0]);err=abs(actual-expected_arrival)/max(expected_arrival,1.)
        largest=max(largest,err);assert err<1e-4,(actual,expected_arrival)
    valid=fitter.validate_contact(fit,obs)
    wrong=dict(fit,theta=np.array([np.log(k*2.)]))
    invalid=fitter.validate_contact(wrong,obs)
    assert valid['valid'] and not invalid['valid'],(valid,invalid)
    # Weak excitation does NOT become evidence of a correct contact law.
    weak=dict(obs,y=surface+fe.mean()/k,fe=np.full(16,fe.mean()))
    undetectable=fitter.validate_contact(wrong,weak)
    assert undetectable['valid']
    return dict(max_arrival_relative_error=largest,correct_contact=valid,wrong_contact=invalid,
                wrong_contact_with_no_excitation=undetectable)


def main():
    checks=check_retirement();mem=Memory(.5,.2); rows=[]; check=[]
    for path in sorted((HERE/'results').glob('finite_patch_input_*.npz')):
        fam=path.stem.removeprefix('finite_patch_input_').rsplit('_',1)[0]
        a=np.load(path); z=a['measured']; signals=Signals(mem); history=[]; old=None
        numerical=ContactFit(memory=mem,spatial=False,arrival=False,blocks=4)
        for j,row in enumerate(z):
            signals.tau_servo=float(a['baseline_tau'][j]); signals.append(row)
            if (j+1)%16 or j<63: continue
            obs=signals.window()
            if obs is None: history=[];old=None;continue
            u,v=linear_profile(obs,mem)
            heldout=None if old is None else float(np.sqrt(np.mean((old*u-v)**2)))
            history.append((u,v));history=history[-4:]
            uu=np.concatenate([b[0] for b in history]);vv=np.concatenate([b[1] for b in history])
            k=float(uu@vv/max(uu@uu,1e-30));ks=float(u@v/max(u@u,1e-30))
            rr=vv-k*uu; current=v-k*u
            # Different surface polynomials are removed for each short block.
            out=dict(family=fam,time=float(row[0]),k_shared=1.5*k,k_short=1.5*ks,
                     baseline=float(a['baseline_kgood'][j]),reference=float(a['scoring_truth_k'][j]),
                     training_rms=float(np.sqrt(np.mean(rr**2))),current_rms=float(np.sqrt(np.mean(current**2))),
                     heldout_rms=heldout,relative_noise=.04/max(k*np.linalg.norm(uu),1e-30))
            if old is not None and old>0:
                out['previous_contact_validation']=numerical.validate_contact(dict(theta=np.array([np.log(old)])),obs)
            rows.append(out);old=k
            # Check exact linear limit on disjoint windows, including after
            # history retirement. Arrival disabled in BOTH formulations.
            if fam=='S2' and j<160:
                fit=numerical.fit(obs,max(float(a['baseline_kraw'][j]),200.))
                err=abs(np.log(fit['k']/(1.5*k)));assert err<2e-5,(err,fit['k'],k)
                check.append(err)
        rr=[r for r in rows if r['family']==fam]
        good=[r for r in rr if r['reference']>0 and r['k_shared']>0]
        print(fam,'n',len(good),'shared/short log errors',
              np.round([np.median([abs(np.log(r[k]/r['reference'])) for r in good if r[k]>0])
                        for k in ('k_shared','k_short')],3),flush=True)
    out=dict(protocol=dict(input='frozen measured baseline trajectories',memory=dict(alpha=.5,tau=.2),
                           window=16,blocks=4,arrival=False,truth_for_scoring_only=True),
             verification=checks,exact_linear_limit_max_log_error=max(check),rows=rows)
    (HERE/'results'/'dynamic_contact_linear_audit.json').write_text(json.dumps(out,indent=2),encoding='utf8')


if __name__=='__main__':main()
