"""Small paired mechanism experiment for shared-motion dynamic compliance."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from run_joint import MeasuredOnly,common,CE
from response_policy import ResponsePolicy

HERE=Path(__file__).resolve().parent


def episode(job):
    scene,model,label=job
    fam,sid,p,T=scene
    policy=ResponsePolicy(response_model='free' if model=='free_dynamic' else 'shared',
        conversion='scalar' if model=='shared_scalar' else 'dynamic',history=model.startswith('shared_history'),
        uncertainty=model=='shared_history_uncertain',single_model=False,**CE.FINAL)
    log,_=common.run(MeasuredOnly(policy),p,T,sid,platform=dict(common.P.PLATFORMS['PC125']))
    frozen=np.load(next((HERE/'results').glob(f'lean_input_{fam}_*.npz')))
    early=log[:,0]<1.32
    np.testing.assert_array_equal(log[early,11],frozen['measured'][early,4])
    file=HERE/'results'/f'{label}_{fam}_{model}.npz'
    if file.exists():raise ValueError('preserve raw observations')
    np.savez_compressed(file,measured=log[:,[0,1,2,5,11,12,13]],
        response_trace=np.asarray(policy.response_trace),controller_trace=np.asarray(policy.trace))
    fits=policy.response_fits
    return dict(family=fam,sid=sid,model=model,metrics=common.task_metrics(log,20.8),
        baseline_metrics=json.loads(str(frozen['metadata']))['metrics'],baseline_sha256=hashlib.sha256(
            next((HERE/'results').glob(f'lean_input_{fam}_*.npz')).read_bytes()).hexdigest(),
        measurements=file.name,measurements_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        used=policy.response_used,published_models=sum(f['accepted'] for f in fits),fits=fits,
        fit_ms_P50_P90_P99=np.percentile([f['fit_seconds']*1000 for f in fits],[50,90,99]).tolist(),
        exact_common_startup=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--label',default='shared_response_v1')
    ap.add_argument('--families',nargs='+',default=['T3m','T2m-stiff'])
    ap.add_argument('--models',nargs='+',default=['free_dynamic','shared_dynamic','shared_scalar'],
                    choices=['free_dynamic','shared_dynamic','shared_scalar','shared_history','shared_history_uncertain'])
    args=ap.parse_args();out=HERE/'results'/(args.label+'.json')
    if out.exists():raise ValueError('preserve primary output')
    scenes=[s for s in CE.cohort('diag',1) if s[0] in args.families]
    result=dict(protocol=dict(question='can shared physical motion and a control-functional dynamic inverse improve the original method?',
        models=args.models,expected_rows=len(scenes)*len(args.models),environment_unchanged=True,
        original_excitation_and_upper_control=True,delayed_publication_seconds=.032,fit_interval_seconds=.032,
        inverse_initial_state='instantaneous at first use; shared_history instead warms the inverse over the full fitted window, with initial sensitivity reported',
        uncertainty='same-model conditional local design reference only',relaxation_time_hypothesis=.04,
        final_used=False,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),HERE/'response_policy.py',HERE/'structured_memory.py',HERE/'dynamic_compliance.py',
             HERE/'joint_policy.py',HERE/'local_contact.py',HERE/'memory_contact.py')}),rows=[])
    with ProcessPoolExecutor(max_workers=2) as pool:
        for r in pool.map(episode,[(s,m,args.label) for s in scenes for m in args.models]):
            result['rows'].append(r);out.write_text(json.dumps(result,indent=2),encoding='utf8')
            print(r['family'],r['model'],'RMSE',r['metrics']['rmse'],'baseline',r['baseline_metrics']['rmse'],
                  'used',r['used'],'published',r['published_models'],'fit ms',r['fit_ms_P50_P90_P99'],flush=True)


if __name__=='__main__':main()
