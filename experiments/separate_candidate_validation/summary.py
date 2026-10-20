"""Paired, per-weight summaries; no selection or new experiments."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def paired(candidate,baseline):
    v=100*(baseline-candidate)/np.maximum(baseline,1e-12);rng=np.random.default_rng(4141)
    means=v[rng.integers(len(v),size=(10000,len(v)))].mean(axis=1)
    return {'n':len(v),'mean_percent_improvement':float(v.mean()),'ci95_percent':np.quantile(means,[.025,.975]).tolist(),'wins':int(np.sum(v>0)),'ratio_of_mean_costs_percent':float(100*(baseline.mean()-candidate.mean())/baseline.mean())}


def main():
    all_results={};lines=['# Independent numerical results for the separate candidates','',
    'This is a reduced test model, not complete B0, an original-paper reproduction or a judgment of leading-journal readiness. There are 24 independent physical scenarios reused at each weight, not 72 independent samples. Percentages are scenario-paired means with bootstrap 95% intervals. Positive values mean lower aggregate candidate cost, not improvement in every force/motion component.','',
    'Main tables use zero release wait and each prototype\'s design frequency (model 12.5Hz, inexpensive feedback/network 100Hz); they do not establish that ordinary MPC cannot run faster. Actual computational release and matched-frequency checks are separate.','']
    for token,w in [('w1',1.),('w12_5',12.5),('w25',25.)]:
        methods={};files=[(f'candidate1_test_{token}.json',None),(f'candidate2_test_{token}.json','candidate2'),(f'candidate2_no_front_test_{token}.json','candidate2_no_front')]
        for name,rename in files:
            path=ROOT/name
            if not path.exists():continue
            d=json.loads(path.read_text());assert d['config']['seed']==9401 and len(d['parameters'])==24
            for r in d['results']:
                key=rename or r['method'];methods.setdefault(key,[]).append(r)
        rows={}
        for method,rs in methods.items():
            rs.sort(key=lambda r:r['case']);cost=np.array([r['metrics'][0] for r in rs]);dt=np.concatenate([r['decision_ms'] for r in rs])
            rows[method]={'n':len(rs),'mean_cost':float(cost.mean()),'mean_force_rms_N':float(np.mean([r['metrics'][1] for r in rs])),'mean_tangent_rms_mm':float(np.mean([r['metrics'][2] for r in rs])*1000),'worst_peak_N':float(max(r['metrics'][3] for r in rs)),'decision_median_ms':float(np.median(dt)),'decision_p95_ms':float(np.percentile(dt,95)),'min_energy_margin':float(min(r['energy_margin_min'] for r in rs)),'max_gain_violation':float(max(r['gain_violation'] for r in rs)),'energy_active_fraction':float(np.mean([r['metrics'][5] for r in rs]))}
            if method=='nnbo_adapted':
                failures=[r['reference_failures'] for r in rs]
                rows[method]['reference_failures_min_median_max']=np.quantile(failures,[0,.5,1]).tolist()
                rows[method]['validated_main_comparator']=False
        comparisons={}
        for candidate in ('rate_neighbor','candidate2'):
            if candidate not in methods:continue
            cc=np.array([r['metrics'][0] for r in methods[candidate]])
            for baseline in ('nnbo_adapted','mpvic_os_adapted','rate_mpc','candidate2_no_front'):
                if baseline in methods:
                    bb=np.array([r['metrics'][0] for r in methods[baseline]]);comparisons[f'{candidate}_vs_{baseline}']=paired(cc,bb)
                    if baseline=='nnbo_adapted':
                        comparisons[f'{candidate}_vs_{baseline}']['valid_for_contribution']=False
                        comparisons[f'{candidate}_vs_{baseline}']['reason']='Nonzero-force task adaptation has substantial reference inversion failures; cannot claim superiority to the original NNBO.'
        all_results[token]={'weight':w,'methods':rows,'paired':comparisons}
        lines+= [f'## Normal-pose weight {w:g}','', '| Method | Mean aggregate cost | Mean force RMS/N | Mean tangential RMS/mm | Peak over all scenarios/N |','|---|---:|---:|---:|---:|']
        for method,row in rows.items():lines.append(f"| {method} | {row['mean_cost']:.6g} | {row['mean_force_rms_N']:.5g} | {row['mean_tangent_rms_mm']:.5g} | {row['worst_peak_N']:.5g} |")
        lines+=['','The NNBO adaptation has substantial inversion failures: this row is diagnostic only, and related improvement percentages are not contribution evidence. Original-paper/nonzero-target-force adaptation still needs repair and verification. The E41 row is limited to this directory\'s declared OS scanning adaptation, not a conclusion about the whole original paper.']
        lines+=['','| Paired comparison | Mean improvement/% | 95% interval/% | Improved scenarios |','|---|---:|---|---:|']
        for name,row in comparisons.items():
            if row.get('valid_for_contribution') is False:continue
            lo,hi=row['ci95_percent'];lines.append(f"| {name} | {row['mean_percent_improvement']:.4f} | [{lo:.4f}, {hi:.4f}] | {row['wins']}/{row['n']} |")
        lines+=['']
    native={};lines+=['## Matched-frequency and computation-wait checks','', '6 seed9201 validation scenarios use 100Hz requests and release according to measured local computation time; they are separate from the 24 independent zero-wait scenarios above. Front-end costs, preparation/tail time and unreleased updates must be included. This is not robot real-time certification. NNBO remains an adaptation whose inversion has not passed verification.','', '| Weight | Method | Median decision/ms | p95/ms | Requests skipped before release/3600 |','|---|---|---:|---:|---:|']
    for token in ('w1','w12_5','w25'):
        path=ROOT/f'native_timing_{token}.json'
        if not path.exists():continue
        data=json.loads(path.read_text());native[token]={}
        for method in ('rate_neighbor','rate_mpc','candidate2'):
            rs=[r for r in data['results'] if r['method']==method];times=np.concatenate([r['decision_ms'] for r in rs]);skipped=sum(r['skipped_updates'] for r in rs)
            row={'n':len(rs),'mean_cost':float(np.mean([r['metrics'][0] for r in rs])),'median_decision_ms':float(np.median(times)),'p95_decision_ms':float(np.percentile(times,95)),'skipped_requests':skipped,'requests':3600}
            native[token][method]=row;lines.append(f"| {data['config']['weight']:g} | {method} | {row['median_decision_ms']:.3f} | {row['p95_decision_ms']:.3f} | {skipped} |")
    (ROOT/'native_summary.json').write_text(json.dumps(native,indent=2))
    (ROOT/'independent_summary.json').write_text(json.dumps(all_results,indent=2));(ROOT/'independent_results.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({w:{k:v for k,v in r['paired'].items() if v.get('valid_for_contribution',True)} for w,r in all_results.items()},indent=2))


if __name__=='__main__':main()
