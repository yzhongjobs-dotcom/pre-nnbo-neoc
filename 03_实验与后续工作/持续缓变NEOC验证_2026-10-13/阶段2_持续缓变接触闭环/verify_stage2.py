from pathlib import Path
import hashlib
import json
import subprocess
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
NUMERIC=['cost','force_rmse_N','sampled_constraints_pass','final_gain','minimum_energy']


def hashes(folder):
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.glob('trace_*.csv'))}


def main():
    original=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    optimized=json.loads((HERE/'01_复用动作计算灵敏度/results.json').read_text(encoding='utf-8'))
    common=[];prediction_differences=[]
    for old,new in zip(original['scenarios'],optimized['scenarios']):
        assert (old['path'],old['duration_s'],old['period_s'])==(new['path'],new['duration_s'],new['period_s'])
        for a,b in zip(old['scores'],new['scores']):
            assert a['method']==b['method']
            if a['method'] not in ('local_NEOC','continuous_NEOC'):
                assert all(a[key]==b[key] for key in NUMERIC)
                stem=f"{old['path']}_{old['duration_s']:g}s_{old['period_s']:g}s_{a['method']}"
                filename='trace_'+stem+'.csv'
                assert (HERE/filename).read_bytes()==(HERE/'01_复用动作计算灵敏度'/filename).read_bytes()
                common.append(stem)
            else:
                gap=abs(a['cost']-b['cost'])/max(abs(a['cost']),1e-15)
                assert gap<1e-6
                prediction_differences.append(dict(path=old['path'],duration=old['duration_s'],
                    period=old['period_s'],method=a['method'],relative_cost_difference=gap))
    smooth=HERE/'02_平滑二次critic对照'
    first=json.loads((smooth/'results.json').read_text(encoding='utf-8'));before=hashes(smooth)
    (smooth/'results_run1.json').write_text(json.dumps(first,indent=2),encoding='utf-8')
    with (HERE/'smooth_reproduction_log.txt').open('w',encoding='utf-8') as log:
        subprocess.run([sys.executable,'-X','utf8',str(HERE/'smooth_contact_check.py')],
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    second=json.loads((smooth/'results.json').read_text(encoding='utf-8'));after=hashes(smooth)
    assert before==after
    smooth_fields=['cost','force_rmse_N','sampled_constraints_pass','maximum_P_relative_error',
                   'maximum_K_error','maximum_D_error','anchor_solves','sensitivity_evaluations']
    for a,b in zip(first['scenarios'],second['scenarios']):
        for old,new in zip(a['scores'],b['scores']):
            assert old['method']==new['method']
            assert all(old[key]==new[key] for key in smooth_fields)
    # Separate repeat of optimized NEOC in the dense physical scenario.
    import run_contact_validation as trial
    from continuous_contact import UniformCarrier
    from sensitivity_fast import fast_sensitivity
    dense=HERE/'03_小邻域与密集更新_探索性'
    repeat=dense/'复现抽查';repeat.mkdir(exist_ok=True)
    reference=json.loads((dense/'results.json').read_text(encoding='utf-8'))
    target=next(r for r in reference['scores'] if r['method']=='local_NEOC')
    UniformCarrier.sensitivity=fast_sensitivity;trial.RADIUS=.0125;trial.HERE=repeat
    row=trial.simulate(UniformCarrier(),'local_NEOC','ramp',6.,.1)
    for key in NUMERIC:assert row[key]==target[key]
    filename='trace_ramp_6s_0.1s_local_NEOC.csv'
    assert (dense/filename).read_bytes()==(repeat/filename).read_bytes()
    (repeat/'results.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
    record=dict(status='VERIFIED_STATED_SUBSETS',common_nonlinear_baseline_trajectories_byte_identical=common,
        optimized_vs_original_NEOC_cost_differences=prediction_differences,
        smooth_all_32_trajectories_two_runs_byte_identical=after,
        dense_NEOC_repeat_trace_sha256=hashlib.sha256((dense/filename).read_bytes()).hexdigest(),
        limits='Timing is variable; not every nonlinear NEOC scenario repeated byte-for-byte; not interval certification.')
    (HERE/'reproducibility.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print('Verified 16 nonlinear baseline trajectories, 32 smooth trajectories, dense NEOC repeat.',flush=True)


if __name__=='__main__':main()
