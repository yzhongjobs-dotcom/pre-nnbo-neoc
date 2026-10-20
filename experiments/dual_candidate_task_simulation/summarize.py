"""Report shared physical errors separately from gain-preference OCP costs."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent

def load(name):
    data=json.loads((ROOT/(name+'.json')).read_text());arrays=np.load(ROOT/(name+'_traces.npz'))
    rows=data['results'];scores=[]
    for r in rows:
        a=arrays[str(r['case'])];dt=float(a[1,0]-a[0,0])
        stage=(a[:,1]/.002)**2+(a[:,2]/.001)**2+.04*np.sum((a[:,3:5]/.015)**2,axis=1)+((a[:,10]-.6)/.1)**2
        terminal=.1*((a[-1,1]/.002)**2+(a[-1,2]/.001)**2+((a[-1,10]-.6)/.1)**2)
        scores.append(float(dt*stage.sum()+terminal))
    return data,rows,np.asarray(scores)

def stats(name):
    data,rows,score=load(name)
    keys=['force_rmse_N','tangent_rmse_mm','normal_rmse_mm','cost','decision_median_ms','decision_p95_ms','mhe_total_seconds']
    out={k:float(np.mean([a[k] for a in rows])) for k in keys}
    out['common_sampled_physical_score']=float(score.mean());out['n_cases']=len(rows)
    out['violating_samples']=sum(a['constraint_violating_samples'] for a in rows)
    out['skipped_updates']=sum(a.get('skipped_busy_updates',0) for a in rows)
    out['expired_policy_queries']=sum(a.get('expired_policy_queries',0) for a in rows)
    out['case_ids']=[a['case'] for a in rows]
    out['strata']={str(j):{k:float(np.mean([r[k] for r in rows if r['kind']==j])) for k in keys[:4]} for j in range(4)}
    return out

def paired(a,b,metric='cost'):
    da,ra,sa=load(a);db,rb,sb=load(b)
    assert [x['case'] for x in ra]==[x['case'] for x in rb]
    assert da['cases']==db['cases']
    va=sa if metric=='physical' else np.array([r[metric] for r in ra]);vb=sb if metric=='physical' else np.array([r[metric] for r in rb])
    d=100*(vb-va)/np.maximum(vb,1e-10);rng=np.random.default_rng(16901)
    boot=np.mean(d[rng.integers(0,len(d),(10000,len(d)))],axis=1)
    return {'metric':metric,'mean_paired_improvement_percent':float(np.mean(d)),'paired_bootstrap_95_percent':np.percentile(boot,[2.5,97.5]).tolist(),
            'ratio_of_means_improvement_percent':float(100*(vb.mean()-va.mean())/vb.mean()),'n':len(d)}

def main():
    names=['candidate1_test','candidate1_without_neighbor_test','candidate1_full_value_reference_test',
           'candidate2_test_seed16111','candidate2_test_seed16112','candidate2_test_seed16113','jin_ms_test',
           'candidate1_preview_audit_full','candidate1_preview_audit_removed','candidate2_preview_audit_full','candidate2_preview_audit_removed',
           'candidate1_shifted_test','candidate2_shifted_test','candidate1_without_moment_test',
           'compiled_timed_candidate1','compiled_timed_candidate1_without_neighbor','compiled_timed_full_value','compiled_timed_candidate2','compiled_timed_jin_ms']
    out={'scope':'Synthetic planar model. Primary nominal results use zero computation wait. Timing rows are a separate isolated 8-case publication experiment. No real robot or journal acceptance claim.',
         'shared_metric':'sampled physical tracking/velocity/force plus endpoint term; no gain preference or gain-rate penalty',
         'confidence_scope':'paired case bootstrap for fixed fitted methods; no simultaneous-inference correction; neural seed summaries are reported separately',
         'methods':{name:stats(name) for name in names},'comparisons':{}}
    for name,a,b,metric in [
        ('neighbor_effect','candidate1_test','candidate1_without_neighbor_test','cost'),
        ('value_quality_reference','candidate1_test','candidate1_full_value_reference_test','cost'),
        ('preview_candidate1','candidate1_preview_audit_full','candidate1_preview_audit_removed','cost'),
        ('preview_candidate2_seed16111','candidate2_preview_audit_full','candidate2_preview_audit_removed','cost'),
        ('candidate1_vs_jin_shared_physical','candidate1_test','jin_ms_test','physical'),
        ('candidate2_vs_jin_shared_physical','candidate2_test_seed16111','jin_ms_test','physical'),
        ('timed_candidate1_vs_jin','compiled_timed_candidate1','compiled_timed_jin_ms','physical'),
        ('timed_candidate2_vs_jin','compiled_timed_candidate2','compiled_timed_jin_ms','physical')]:out['comparisons'][name]=paired(a,b,metric)
    seed_methods=[out['methods'][f'candidate2_test_seed{s}'] for s in [16111,16112,16113]]
    out['withdrawn_ablation']={'files':['candidate1_without_preview_test.json','candidate2_without_preview_test.json'],
        'reason':'old-region-only masking also removed updated current material after crossing; not a valid spatial-preview ablation',
        'correction':'uniform current-contact equivalent retains estimated current Fn/Ft; separate no-preview actor retrained; new seed16501 cohort n24',
        'old_improvement_numbers_not_used':[26.86276488049843,6.172869256143549]}
    out['candidate2_seed_mean_and_range']={k:{'mean':float(np.mean([s[k] for s in seed_methods])),'range':[min(s[k] for s in seed_methods),max(s[k] for s in seed_methods)]}
        for k in ['force_rmse_N','tangent_rmse_mm','normal_rmse_mm','cost','common_sampled_physical_score']}
    (ROOT/'summary.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
    print(json.dumps({'comparisons':out['comparisons'],'neural_seeds':out['candidate2_seed_mean_and_range']},indent=2))
    lines=['# Separate simulations and comparison results for the two candidates','','24 physical scenarios, 6 per category. Three learning seeds reuse the same scenarios; these are not 72 independent physical trials. The main table assumes zero computation wait.',
           '','| Method | Force RMSE/N | Tangential RMSE/mm | Normal RMSE/mm | Shared physical score (lower is better) | Sampled constraint violations |','|---|---:|---:|---:|---:|---:|']
    for label,key in [('Candidate 1: local value feedback + neighboring updates','candidate1_test'),('Candidate 1 without neighboring updates','candidate1_without_neighbor_test'),('More thorough same-family value solve','candidate1_full_value_reference_test'),('Candidate 2: seed16111','candidate2_test_seed16111'),('Candidate 2: seed16112','candidate2_test_seed16112'),('Candidate 2: seed16113','candidate2_test_seed16113'),('Jin multistep scanning adaptation','jin_ms_test')]:
        s=out['methods'][key];lines.append(f"| {label} | {s['force_rmse_N']:.6f} | {s['tangent_rmse_mm']:.4f} | {s['normal_rmse_mm']:.4f} | {s['common_sampled_physical_score']:.5f} | {s['violating_samples']} |")
    lines+=['','The shared physical score contains only motion, velocity and normal force, integrated with the same sampling; it excludes controller-family gain-preference/rate regularization. A large Jin cost under this candidate\'s complete OCP objective does not establish a win. Component errors show the tradeoff between normal position and force/tangential motion.',
            '','| Measured-wait release: 8 scenarios | Force RMSE/N | Tangential RMSE/mm | Normal RMSE/mm | Mean of per-scenario median decision times/ms | Skipped updates |','|---|---:|---:|---:|---:|---:|']
    for label,key in [('Candidate 1','compiled_timed_candidate1'),('Without neighboring prediction','compiled_timed_candidate1_without_neighbor'),('More thorough value solve','compiled_timed_full_value'),('Candidate 2','compiled_timed_candidate2'),('Jin scanning adaptation','compiled_timed_jin_ms')]:
        s=out['methods'][key];lines.append(f"| {label} | {s['force_rmse_N']:.6f} | {s['tangent_rmse_mm']:.4f} | {s['normal_rmse_mm']:.4f} | {s['decision_median_ms']:.3f} | {s['skipped_updates']} |")
    lines+=['','Timing uses CPU prototypes with a shared compiled scalar contact kernel and a 2ms release grid; it is not robot real-time certification. Old timed_* files diagnose the pointwise NumPy version and do not establish current speed. Model updates and individual network inference have different granularity, so their two median times are not an end-to-end speedup ratio. Estimation overhead is reported separately and included in release delay.',
            '','See [summary.json](summary.json) for mechanisms, stratified data, confidence intervals, distribution shifts and training-seed scope. Narrow intervals reflect scenario variation for fixed implementations, not attainment of leading-journal standards.']
    lines+=['','Per-case components and example trajectories: [PNG figure](results.png), [PDF figure](results.pdf).']
    def preview_text(key):
        d=out['comparisons'][key];ci=d['paired_bootstrap_95_percent']
        return f"{d['mean_paired_improvement_percent']:.3f}%, paired 95% interval [{ci[0]:.3f}%,{ci[1]:.3f}%]"
    lines+=['','## What these results establish','',
        'The preview ablation was checked and redone: the old version also removed new-material information after it entered the current region, so the original 26.86%/6.17% is not contribution evidence. The correction retains a uniform model matched to current-contact Fn/Ft and removes spatial variation only; the learning version is independently retrained. Comparison uses 24 new seed16501 scenarios without selecting parameters from the original test.',
        'This comparison matches current Fn/Ft, not every moment or local derivative. It tests the overall role of spatially parameterized information; benefits cannot be attributed solely to temporal anticipation.',
        '',
        'Candidate 1 corrected preview ablation: mean paired improvement in complete task cost is '+preview_text('preview_candidate1')+'. NEOC on/off still gives only about 0.007% improvement, with interval [-0.260%,0.207%]; neighboring prediction is currently slower than disabling it under the same budget. Information and neighboring-update value must be judged separately. More thorough same-family solving is not ground-truth globally optimal control; tiny measured wins cannot mean surpassing the optimum.',
        '',
        'Candidate 2 corrected preview ablation (same training seed16111, separately trained): paired complete-cost improvement is '+preview_text('preview_candidate2_seed16111')+'; this is conditional on a fixed seed. The three main seeds have mean shared physical score 0.06029 versus 0.05646 for the Jin adaptation, about 6.78% higher. Fast individual inference or a one-seed ablation does not establish overall superiority to in-field methods.',
        '',
        'Both candidates have slightly smaller normal-displacement errors but larger force/tangential-motion errors; this is an actual tradeoff. No leading-journal-level method advantage is established. Across 8 softer/faster-surface scenarios, candidate 1/2 mean force RMSE is about 0.01792/0.02482N. This extra cohort has no external baseline and does not establish superiority to prior work. Constraint checks concern sampled trajectories; energy constraints are inactive in the main scenarios, so these results do not establish active-constraint behavior or real-contact safety.',
        '',
        'Next explain and improve why candidate 1\'s extra prediction/correction does not yield better cost than ordinary feedback updates, and why candidate 2 uses preview yet trails strong force/tangential performance; check training objectives against actual task tradeoffs. Existing strong impedance structures or learned corrections may be improved, with explicit changes verified on new scenarios. Changing scores, inventing computational bottlenecks or counting added modules cannot establish success. The actual platform/calibration still needs verification.']
    (ROOT/'results.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
