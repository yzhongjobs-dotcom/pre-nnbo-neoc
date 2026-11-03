"""Small, purpose-linked evidence chain; preserve every negative pilot."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(label):
    a=json.loads((HERE/'results'/(label+'.json')).read_text(encoding='utf8'))
    assert len(a['rows'])==a['protocol']['expected_rows'],f'incomplete {label}'
    assert not a['protocol']['plant_changed'] and not a['protocol']['final_used']
    return a


def find(data,family,model):
    rows=[r for r in data['rows'] if r['family']==family and r['model']==model]
    assert len(rows)==1,(family,model)
    return rows[0]


def compact(r):
    return {k:r.get(k) for k in ('family','model','sid','noise_sid','metrics','used','call_ms_P50_P90_P99','servo_id_last')}


def main():
    labels=['joint_pilot_v1','joint_response_pilot_v2','joint_abrupt_motion_v3','joint_value_attribution_v4','joint_noise_pair_v5','joint_lag_pair_v6']
    data={n:read(n) for n in labels};v2=data[labels[1]];v4=data[labels[3]];rows=[]
    assert v4['protocol']['probe']['selected_kind']=='original_sine'
    for family in ('S2','T3m'):
        a=find(v2,family,'varying_no_rate');b=find(v4,family,'joint_lean')
        assert a['metrics']==b['metrics'],'the selected original probe must reproduce the lean estimator'
        rows.append(dict(question='material transition' if family=='S2' else 'movement/material confounding',
                         baseline=compact(find(v2,family,'baseline')),candidate=compact(b)))
    stress=data[labels[2]]
    rows.append(dict(question='abrupt motion boundary',baseline=compact(find(stress,'T2m-stiff','baseline')),
                     candidate=compact(find(stress,'T2m-stiff','varying_no_rate'))))
    for label,question in [(labels[4],'same physical S2, changed sensor noise realization'),
                           (labels[5],'same physical S2, 50 percent unknown servo-lag increase')]:
        rows.append(dict(question=question,baseline=compact(find(data[label],'S2','baseline')),
                         candidate=compact(find(data[label],'S2','joint_lean'))))
    result=dict(protocol=dict(input_files=labels,diagnosis_only=True,uniform_large_cohort_requirement=False,
                               one_example_per_scientific_condition=True,
                               selected_probe_matches_original_exactly=True,complete_batches=True),
                paired_evidence=rows,
                ablations=[compact(r) for r in v2['rows']],
                selection_value_diagnostic=[compact(r) for r in v4['rows']],
                negative_command_budget_pilot=[compact(r) for r in data[labels[0]]['rows']],
                limitations=['small reused development cases, not independent confirmation',
                             'conditional covariance does not cover structural bias or closed-loop endogeneity',
                             'no full constitutive identification or hard real-time guarantee',
                             'novelty and T-RO/IJRR contribution not established'])
    (HERE/'results'/'joint_evidence_summary.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    for r in rows:print(r['question'],r['baseline']['metrics']['rmse'],r['candidate']['metrics']['rmse'])


if __name__=='__main__':main()
