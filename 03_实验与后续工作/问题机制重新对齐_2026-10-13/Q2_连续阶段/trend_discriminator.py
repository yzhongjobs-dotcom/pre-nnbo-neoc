"""Exclude inaccurate frozen-current forecasts as the cause of update losses."""
import json,numpy as np
from contact_stages import Problem,ROOT,dump
from projected_rival import policy
from independent_audit import independent


def main():
    source=json.loads((ROOT/'strong_rivals'/'results.json').read_text())['best_preview']
    pars=[source['zeta'],source['lead']];out=ROOT/'trend_discriminator';out.mkdir(exist_ok=True)
    results={}
    for span in [250.,700.]:
        p=Problem(theta_span=span);targets=p.targets('preview',pars);results[str(span)]={}
        for fast in [False,True]:
            m,t,tr,z=policy(p,targets,True,causal=True,tail_parameters=pars,trend=True,
                    maxiter=8 if fast else 30,tolerance=1e-5 if fast else 1e-8,
                    use_grid=not fast,warm_start=fast)
            pi=independent(p,z,True);pi.pop('segment_states')
            label='trend_fast' if fast else 'trend'
            results[str(span)][label]=dict(metrics=m,timing=t,independent_actual_plant=pi,
                 forecast='Affine slope from only current and previous observed material; full future material withheld.')
            dump(out/'results.json',results)
            np.save(out/f'span{span}_{label}_desired.npy',z)
            np.savetxt(out/f'span{span}_{label}_trace.csv',tr,delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
            print(span,label,'cost',m['cost'],'worst seconds',t['worst_seconds'],flush=True)


if __name__=='__main__':main()
