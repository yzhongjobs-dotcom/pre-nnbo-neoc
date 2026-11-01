"""Finalize serial timing and interpretation after all jobs stop.
No training, selection or evaluation parameters are modified here.
"""
import argparse,json,pickle,copy,platform,hashlib
from step3b_core import *
from step3b_predictors import benchmark
from step3b_closed_loop import write_md as closed_md

def serial_controller(names,models):
    out={}
    for n in names:
        calls=[]
        for scene in make_scenes(110000,1):
            fam,sid,p,T,code=scene;use_plant(code)
            c=build(n,p,copy.deepcopy(models));log,c=common.run(c,p,T,sid)
            calls.extend(np.asarray(c.t_call[100:])*1e3)
        out[n]=dict(p50_ms=float(np.percentile(calls,50)),p99_ms=float(np.percentile(calls,99)),n_calls=len(calls))
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);a=ap.parse_args()
    op,models=load_models();cl=json.loads((RES/'step3b_closed_loop.json').read_text(encoding='utf8'))
    diag=pickle.loads((a.cache/'diag.pkl').read_bytes())
    bm={'P0':None,**models,'P0r':MultiIMM(**op['selected']['P0r'],third=False)}
    from step3b_predictors import evaluate
    ood=[]
    for code in (5,6,7):ood+=pickle.loads((a.cache/f'ood{code}.pkl').read_bytes())
    op['evaluation']={'id':evaluate(bm,diag),'ood':evaluate(bm,ood)}
    op['ood_by_profile']={str(c):evaluate(bm,[e for e in ood if e['code']==c]) for c in (5,6,7)}
    op['covariance_note']='IMM target variance includes posterior state covariance, between-mode covariance and the frozen-mode white-jerk/snap forecast variance integrated to12ms for f+tau*f_dot. Stationary-mode position random-walk variance is scaled by the horizon. It does not include stiffness/contact bias or mode switches during the forecast horizon. Adding process variance changes only reported coverage, never means, fitted parameters or feed-forward.'
    op['timing']=benchmark(bm,diag)
    op['machine']['processor']=platform.processor();op['machine']['logical_cpus']=os.cpu_count()
    op['protocol_correction']='An initial implementation used50 native+50 augmented validation scenes. After discovering the mismatch, selection was replayed with ONLY50 native190000+ scenes. Architecture, training inputs, all61 prespecified IMM candidates, GP optimization and network update trajectory were unchanged; no pilot ID/OOD/closed-loop result was used to add or alter candidates. Initial labels also had an implementation sign error, corrected and independently checked against Oracle before fitting the delivered models. Final numbers use corrected labels and the exact50-scene validation set. This is exploratory diagnostic work, not a preregistered confirmatory evaluation. The correction repeated computation (same candidate sets), and that cost is not hidden as a single executed pass.'
    (RES/'step3b_predictors.json').write_text(json.dumps(op,indent=2,allow_nan=False),encoding='utf8')
    # Refresh tables, keeping the post-selection diagnostics rendered by analyse.
    from step3b_predictors import write_md as open_md
    path=RES/'step3b_predictors.md';old=path.read_text(encoding='utf8')
    tail=old[old.index('\n## Offset-free diagnostics'):]
    open_md(op);path.write_text(path.read_text(encoding='utf8')+tail+'\n## Protocol correction\n\n'+op['protocol_correction']+'\n',encoding='utf8')
    cl['pooled_serial_timing']=serial_controller(cl['models'],models)
    cl['integration_note']='The legacy FINAL IMM remains only to preserve recovery velocity; new FF mean comes entirely from the replacement. Whole-controller timings include this legacy update and the new predictor. Standalone predictor timing excludes the legacy recovery tracker. P0 is NumPy, P1/P0r filtering is Numba; speed comparisons concern these implementations, not intrinsic algorithm complexity.'
    (RES/'step3b_closed_loop.json').write_text(json.dumps(cl,indent=2,allow_nan=False),encoding='utf8');closed_md(cl)
    path=RES/'step3b_closed_loop.md';s=path.read_text(encoding='utf8')
    s+='\nPooled individual-call hot serial timing, no concurrent batch jobs (ms):\n\n| Model | P50 | P99 | n calls |\n|---|---:|---:|---:|\n'
    for n,r in cl['pooled_serial_timing'].items():s+=f"| {n} | {r['p50_ms']:.3f} | {r['p99_ms']:.3f} | {r['n_calls']} |\n"
    s+='\n'+cl['integration_note']+'\n';path.write_text(s,encoding='utf8')
    print('final timing completed')

if __name__=='__main__':main()
