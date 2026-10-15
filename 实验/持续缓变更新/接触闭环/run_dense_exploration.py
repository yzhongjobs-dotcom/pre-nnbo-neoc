"""Exploratory smaller-neighborhood/dense-query test, declared before execution.

Motivation: original frozen-radius tests had ~1.6% loss. Halve radius and query
more often, without altering physics/objective/1% quality/30% timing thresholds.
This tests the amortization mechanism; radius remains empirical, not certified.
"""
from pathlib import Path
import json
import run_contact_validation as trial
from continuous_contact import UniformCarrier
from sensitivity_fast import fast_sensitivity

HERE=Path(__file__).resolve().parent/'03_小邻域与密集更新_探索性'


def main():
    HERE.mkdir(exist_ok=True)
    trial.HERE=HERE;trial.RADIUS=.0125
    UniformCarrier.sensitivity=fast_sensitivity
    config=dict(stage='EXPLORATORY_AFTER_PREDECLARED_FOUR_CASES',path='ramp',
        duration_s=6.,period_s=.1,local_radius=.0125,changed_from_original=['radius halved','period .25 -> .1'],
        physics_objective_and_quality_gates_changed=False,
        radius_certified=False,publication_latency_injected=False,
        comparators=['one_Newton','full_solve','fixed_gain'],
        cached_chord_retested=False,continuous_NEOC_retested=False)
    (HERE/'config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    scores=[]
    for method in ['local_NEOC','one_Newton','full_solve','fixed_gain']:
        scores.append(trial.simulate(UniformCarrier(),method,'ramp',6.,.1))
        (HERE/'partial_results.json').write_text(json.dumps(scores,indent=2),encoding='utf-8')
    ref=next(row for row in scores if row['method']=='full_solve')
    for row in scores:
        row['cost_extra_fraction']=row['cost']/ref['cost']-1
        row['force_rmse_extra_fraction']=row['force_rmse_N']/ref['force_rmse_N']-1
        row['quality_gate_pass']=bool(row['cost_extra_fraction']<=.01 and row['force_rmse_extra_fraction']<=.05 and row['sampled_constraints_pass'])
    eligible=[row for row in scores if row['method'] in ('one_Newton','full_solve') and row['quality_gate_pass']]
    fastest=min(eligible,key=lambda row:row['total_controller_seconds'])
    for row in scores:
        row['fastest_measured_optimizer']=fastest['method']
        row['total_time_reduction']=1-row['total_controller_seconds']/fastest['total_controller_seconds']
        row['screen_pass_vs_measured_optimizers']=bool(row['quality_gate_pass'] and row['total_time_reduction']>=.3)
    (HERE/'results.json').write_text(json.dumps(dict(config=config,scores=scores),indent=2),encoding='utf-8')
    print(json.dumps(scores),flush=True)


if __name__=='__main__':main()
