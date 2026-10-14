"""Causal material-updating ablation and fair cheap-solver precision repair."""
import json,numpy as np
from contact_stages import Problem,ROOT,dump
from projected_rival import policy
from independent_audit import independent


def main():
    source=json.loads((ROOT/'strong_rivals'/'results.json').read_text())['best_preview']
    pars=[source['zeta'],source['lead']];out=ROOT/'material_update_discriminator';out.mkdir(exist_ok=True)
    all_results={}
    for span in [250.,700.]:
        p=Problem(theta_span=span);targets=p.targets('preview',pars);all_results[str(span)]={}
        for mode in ['frozen','current','current_fast']:
            fast=mode=='current_fast'
            m,t,tr,z=policy(p,targets,True,causal=mode!='frozen',frozen=mode=='frozen',tail_parameters=pars,
                       maxiter=8 if fast else 30,tolerance=1e-5 if fast else 1e-8,
                       use_grid=not fast,warm_start=fast)
            data=dict(metrics=m,timing=t,material_mode=mode,
               solver=dict(maxiter=8 if fast else 30,ftol=1e-5 if fast else 1e-8,grid=not fast,warm=fast))
            # Evaluator integrates actual material, never the frozen forecast.
            pi=independent(p,z,True);pi.pop('segment_states')
            data['independent_actual_plant']=pi
            all_results[str(span)][mode]=data
            dump(out/'results.json',all_results)
            np.save(out/f'span{span}_{mode}_desired.npy',z)
            np.savetxt(out/f'span{span}_{mode}_trace.csv',tr,delimiter=',',header='t,e,v,K,D,E,uK,uD,theta,surface',comments='')
            print(span,mode,'cost',m['cost'],'worst seconds',t['worst_seconds'],
                   'physical margin',m['minimum_margin'],flush=True)
        r=all_results[str(span)]
        r['comparison']=dict(material_update_improvement=1-r['current']['metrics']['cost']/r['frozen']['metrics']['cost'],
          fast_vs_current_excess=r['current_fast']['metrics']['cost']/r['current']['metrics']['cost']-1,
          fast_decisions_over_period=sum(t>p.h for t in r['current_fast']['timing']['decision_seconds']),
          true_material_future_available_to_controller=False,surface_program_shared=True,
          note='Virtual-time quality simulation; computation delay not injected. No claim of WCET or deployed realtime control.')
        dump(out/'results.json',all_results)


if __name__=='__main__':main()
