"""Analyze the existing three-task C1 repair evaluation; runs no simulation."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'contact_update_results'


def read(name):
    return json.loads((OUT/name).read_text(encoding='utf-8'))


def main():
    frozen=read('frozen.json')
    for name,digest in frozen['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    rows=read('test_rows.json')
    assert len(rows)==200
    assert sum(r['method']=='C1_updated' for r in rows)==120
    assert all(r['failed']==0 for r in rows)
    old=list(csv.DictReader((ROOT/'completed_tasks/per_scene.csv').open(encoding='utf-8')))

    def values(method,task,field,cell=None):
        source=rows if method.startswith('C1_') else old
        selected=[r for r in source if r['method']==method and r['task']==task
                  and (cell is None or (int(r['scene'])%4>=2)==(cell=='hard'))]
        # Average neural training seeds within each paired scene.
        groups={}
        for r in selected:
            groups.setdefault(int(r['scene']),[]).append(float(r[field]))
        return np.array([np.mean(groups[i]) for i in sorted(groups)])

    methods=['Jin_MS','Xue_nominal','Learned_VIC','Candidate1_v2','C1_updated','Candidate2']
    table={}
    for method in methods:
        table[method]=[float(values(method,t,'rmse_N',cell).mean()) for t,cell in
                       [('T1_static',None),('T2_motion','soft'),('T2_motion','hard'),('T3_combined',None)]]
    # Five reported paired comparisons, selected after evaluation rather than
    # preregistered. Intervals use Bonferroni coverage for this reported family.
    # Resampling changes neither controllers nor task scenes.
    specifications=[('T2_motion','hard','Candidate1_v2'),
                    ('T2_motion','hard','EtR-A-Haninger'),
                    ('T1_static',None,'Candidate1_v2'),
                    ('T3_combined',None,'Candidate1_v2'),
                    ('T2_motion','hard','C1_old_model_new_recovery')]
    comparisons=[]
    rng=np.random.default_rng(8220)
    for task,cell,reference in specifications:
        d=values('C1_updated',task,'rmse_N',cell)-values(reference,task,'rmse_N',cell)
        boot=d[rng.integers(0,len(d),(50000,len(d)))].mean(axis=1)
        interval=np.quantile(boot,[.005,.995])
        comparisons.append(dict(task=task,cell=cell,reference=reference,n=len(d),
                                mean_difference_N=float(d.mean()),
                                family95_interval_N=interval.tolist()))

    timing={}
    checked=0
    for task in ['T1_static','T2_motion','T3_combined']:
        times=[]
        scenes=np.load(ROOT/'completed_tasks'/f'scenes_{task}.npz')
        for i in range(40):
            data=np.load(OUT/'traces'/f'{task}_C1_updated_{i}.npz')
            assert np.array_equal(data['p'],scenes['p'][i])
            assert np.isfinite(data['log']).all()
            assert np.isfinite(data['commands']).all()
            times.extend(data['update_ms'][10:])
            checked+=1
        timing[task]=dict(p99_ms=float(np.percentile(times,99)),max_ms=float(max(times)))

    recovery={}
    for method in ['Candidate1_v2','C1_old_model_new_recovery','C1_updated']:
        r={}
        for cell in ['soft','hard']:
            recovered=np.concatenate([values(method,'T2_motion',f,cell) for f in ['out_recovered','return_recovered']])
            seconds=np.concatenate([values(method,'T2_motion',f,cell) for f in
                                    ['out_restricted_recovery_s','return_restricted_recovery_s']])
            peak=values(method,'T2_motion','peak_force_N',cell)
            loss=values(method,'T2_motion','low_contact_s',cell)
            r[cell]=dict(recovery_rate=float(recovered.mean()),restricted_mean_recovery_s=float(seconds.mean()),
                         worst_sampled_peak_N=float(peak.max()),mean_sampled_peak_N=float(peak.mean()),
                         mean_low_contact_s=float(loss.mean()),scenes_with_low_contact=int((loss>0).sum()))
        recovery[method]=r
    result=dict(rmse_table_N=table,paired_comparisons=comparisons,
                recovery=recovery,timing=timing,
                verified_original_scene_vectors=checked,
                analysis='Paired mean RMSE differences; 50000 bootstrap resamples; five-comparison Bonferroni family')
    (OUT/'analysis.json').write_text(json.dumps(result,indent=2),encoding='utf-8')

    summary=read('summary.json')
    lines=['# C1 repair on the original three scanning tasks','',
           'The user stopped additional simple-controller tuning and scenario expansion. '
           'Only the selected C1 repair and one estimator ablation were evaluated on the '
           'original seed-7200 tasks (40 paired scenes each). Historical controllers, '
           'weights and results remain unchanged. No neural training was run.','',
           '## Method and capability tested','',
           'The repair profiles a cubic force-drift polynomial out of the stiffness '
           'regression before updating stiffness. It retains the original 128 ms predictive '
           'controller, uses three local solve iterations, and replaces the task-tuned '
           'fixed recovery with stiffness-based recovery and command-state compensation. '
           'The estimated drift is an equivalent contact response, not established anatomical motion.','',
           'The tested capability is maintaining contact and limiting impact during sustained '
           'push/retract on stiff tissue, while retaining useful combined-scan performance. '
           'This is not a first-joint-estimation or near-optimal-control claim.','',
           '## Three-task outcomes','',
           '| Method | T1 stiffness | T2 soft motion | T2 stiff motion | T3 combined |',
           '|---|---:|---:|---:|---:|']
    for method,a in table.items():
        lines.append('| '+method+' | '+' | '.join(f'{v:.3f}' for v in a)+' |')
    lines+=['','RMSE in N. Literature rows are the previously evaluated task adaptations. '
            'T2 was used to develop this repair, whereas those historical rows retain their '
            'T3-selected settings; this is not an equally retuned ranking or a claim about '
            'the original papers\' performance limits.','',
            '## Stiff push/retract and estimator ablation','',
            '| Configuration | RMSE N | Worst sampled peak N | Scenes below 0.5 N | 2 s event recovery |',
            '|---|---:|---:|---:|---:|']
    for method in ['Candidate1_v2','C1_old_model_new_recovery','C1_updated']:
        a=recovery[method]['hard']
        error=float(values(method,'T2_motion','rmse_N','hard').mean())
        lines.append(f"| {method} | {error:.3f} | {a['worst_sampled_peak_N']:.3f} | {a['scenes_with_low_contact']}/20 | {100*a['recovery_rate']:.1f}% |")
    lines+=['','The ablation retains the new recovery and three solve iterations but restores '
            'the original joint nonlinear estimator. Its remaining stiff-motion error '
            'and high-force tail support an additional role for the estimator repair. '
            'The experiment does not separately identify every contribution of iteration '
            'count and recovery.','',
            'The stiff-motion repair costs static accuracy: T1 RMSE increases from 0.060 '
            'to 0.087 N; the worst sampled T1 peak increases from 5.327 to 7.584 N. '
            'T3 mean RMSE remains similar (0.147 to 0.142 N). The new method is not '
            'uniformly better and has no global 7 N guarantee.','',
            '## Paired evidence','',
            '| Task/cell | Reference | Mean difference N | Family 95% interval N |',
            '|---|---|---:|---|']
    for a in comparisons:
        lo,hi=a['family95_interval_N']
        lines.append(f"| {a['task']}/{a['cell'] or 'all'} | {a['reference']} | {a['mean_difference_N']:+.4f} | [{lo:+.4f}, {hi:+.4f}] |")
    lines+=['','Intervals refer to new minus reference, using paired scene resampling and '
            'Bonferroni coverage for the five comparisons above. They do not remove the '
            'different development exposure noted above.','',
            '## Update time and reproducibility','',
            '| Task | P99 ms | Maximum measured ms |','|---|---:|---:|']
    for task,a in timing.items():
        lines.append(f"| {task} | {a['p99_ms']:.3f} | {a['max_ms']:.3f} |")
    lines+=['','Single-process policy call including observation assembly and output '
            'projection, excluding separate JIT warmup, the first 10 calls of each scene, '
            'the simulated plant and physical communications. '
            'The control period is 8 ms; this is not robot real-time certification.','',
            'Reproduce checks and the frozen cohort from the repository root:',
            '```powershell',
            'python -X utf8 experiments/cobot_scan_benchmark/evaluate_contact_update.py check',
            'python -X utf8 experiments/cobot_scan_benchmark/evaluate_contact_update.py evaluate',
            'python -X utf8 experiments/cobot_scan_benchmark/analyze_contact_update.py',
            '```','',
            'Frozen configuration and source hashes: [frozen.json](frozen.json). '
            'Per-scene outcomes: [test_rows.json](test_rows.json). '
            'Raw trajectories, commands and update times: [traces](traces/). '
            'Estimator and truth-input checks: [checks.json](checks.json). '
            'Development/validation data refer only to C1 and predate test evaluation.','',
            'The closest DOSMAC, Cao and Wang comparisons remain unfinished; exploratory '
            'unverified adaptations are not included as publication evidence. '
            'Closed-loop stability, physical motion recovery and image-quality benefits '
            'are not established by this repair.']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
