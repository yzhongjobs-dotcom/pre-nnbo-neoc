"""Reconcile the complete five-family contact comparison without cherry-picking.

One reused diagnosis scenario per family; results are development evidence, not
independent confirmation. Baseline source is unchanged across these runs and its
two repeated episodes must match exactly. Requires all input batches complete.
"""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    data=json.loads((HERE/'results'/name).read_text(encoding='utf8'))
    expected=data['protocol']['expected_rows']
    if len(data['rows'])!=expected:raise ValueError(f'incomplete batch {name}: {len(data["rows"])}/{expected}')
    return data


def main():
    names=['dynamic_contact_five_family_v5.json','dynamic_contact_level_smoke_v7.json',
           'dynamic_contact_level_remaining_v8.json','dynamic_contact_level_static_v8.json']
    base,pair,remaining,static=map(read,names)
    for name in ('dynamic_contact.py','dynamic_contact_policy.py','finite_patch_model.py','line_patch_model.py'):
        hashes=[x['protocol']['source_sha256'][name] for x in (pair,remaining,static)]
        assert len(set(hashes))==1,(name,hashes)
    baseline={r['family']:r for r in base['rows'] if r['model']=='baseline'}
    for r in pair['rows']:
        if r['model']=='baseline':assert r['metrics']==baseline[r['family']]['metrics']
    candidates={}
    for data in (pair,remaining,static):
        for r in data['rows']:
            key=r['family'],r['model']
            if key in candidates:raise ValueError(f'duplicate {key}')
            candidates[key]=r
    families=['S2','T3m','T2m-soft','T2m-stiff','T2f-stiff']
    models=['baseline','patch_stable_target_level','static_patch_stable_target_level','line_stable_target_level']
    rows=[]
    for family in families:
        for model in models:
            r=baseline[family] if model=='baseline' else candidates[family,model]
            rows.append(dict(family=family,model=model,sid=r['sid'],metrics=r['metrics'],
                             active_ticks=r['model_ticks'],contact_ticks=r['contact_ticks'],
                             unobserved_corrections=r.get('unobserved_corrections',0),
                             fit_ms_P50_P90=r['fit_ms_P50_P90'],
                             call_ms_P50_P90_P99=r.get('call_ms_P50_P90_P99'),
                             parameter_to_force_variance_P50_P90=r.get('parameter_to_force_variance_P50_P90')))
        print(family,*[f'{r["metrics"]["rmse"]:.5f}' for r in rows[-4:]],flush=True)
    out=dict(protocol=dict(input_files=names,diagnosis_scenarios_per_family=1,
                           repeated_baseline_exact_match=True,model_source_hashes_match=True,
                           tuning=False,plant_changed=False,final_used=False,
                           actual_computational_delay_simulated=False),
             limitations=['reused development scenes, not confirmation',
                          'memory and geometry hypotheses not covered by parameter covariance',
                          'cross-update parameter correlation and contact model bias unresolved',
                          'not real-time or deployment evidence'],rows=rows,
             memory_sensitivity=[dict(family=r['family'],model=r['model'],sid=r['sid'],metrics=r['metrics'],
                                      active_ticks=r['model_ticks'],
                                      memory=static['protocol']['configs'][r['model']]['memory'])
                                 for r in static['rows'] if r['model'].endswith(('_low','_high'))])
    (HERE/'results'/'dynamic_contact_five_family_summary.json').write_text(json.dumps(out,indent=2),encoding='utf8')


if __name__=='__main__':main()
