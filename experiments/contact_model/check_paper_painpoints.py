"""Check whether a paper's sliding-friction pain point is represented by this benchmark.

Five existing diagnosis scenarios, unchanged V1+Proposed controller. Change only
the declared friction coefficient to 0.5x/1.5x, use identical sensor noise seeds.
This is a structural sensitivity diagnostic, NOT a candidate performance test.
No changed contact law, actuator, controller tuning or confirmation cohort.
"""
import hashlib
import json
import sys
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent/'method_redesign'))
import common
import cohort_eval as CE
from v1_servo_id import ServoID


def main():
    rows = []
    for fam, sid, params, duration in CE.cohort('diag', 1):
        ref, _ = common.run(ServoID(single_model=False, **CE.FINAL), params, duration, sid)
        ref_metrics = common.task_metrics(ref, 20.8)
        for factor in (.5, 1.5):
            p = params.copy(); p[1] *= factor
            log, _ = common.run(ServoID(single_model=False, **CE.FINAL), p, duration, sid)
            metrics = common.task_metrics(log, 20.8)
            equal = log.shape == ref.shape
            deltas = {}
            for name, cols in {'position':[1,2], 'normal_force':[3,5], 'command':[11],
                               'moment':[13,14,15]}.items():
                same = equal and np.array_equal(log[:, cols].view(np.uint64), ref[:, cols].view(np.uint64))
                deltas[name] = dict(bit_identical=bool(same), max_abs_delta=float(np.max(np.abs(log[:,cols]-ref[:,cols]))) if equal else None)
            row = dict(family=fam, sid=sid, friction_factor=factor,
                       original_mu=float(params[1]), changed_mu=float(p[1]),
                       unchanged_controller='ServoID(single_model=False, **CE.FINAL)',
                       channels=deltas, nominal_force_RMSE=ref_metrics['rmse'],
                       changed_force_RMSE=metrics['rmse'],
                       measured_tangential_RMS_difference=float(np.sqrt(np.mean((log[:,12]-ref[:,12])**2))))
            assert all(v['bit_identical'] for v in deltas.values()), row
            assert row['measured_tangential_RMS_difference'] > 0, row
            rows.append(row)
            print('%s mu x%.1f: Fn/y/command/My bit-identical; Ft RMS change %.4f N' %
                  (fam, factor, row['measured_tangential_RMS_difference']), flush=True)
    output = dict(protocol=dict(cohort='diag', scenarios_per_family=1, controller_retuned=False,
                                plant_law_changed=False, final_used=False,
                                purpose='structural friction-to-normal/task coupling check'),
                  source_hashes={str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in [Path(__file__), HERE.parent/'cobot_scan_benchmark/plant.py',
                                           HERE.parent/'motion_aware_scan/proposed.py',
                                           HERE.parent/'method_redesign/v1_servo_id.py']},
                  rows=rows)
    dest = HERE/'results/paper_painpoint_sensitivity.json'
    dest.write_text(json.dumps(output, indent=2), encoding='utf-8')
    print(dest)


def normal_state_diagnostic():
    """Scoring-only access to hidden contact components, never estimator inputs.

    This distinguishes elastic force/depth from the total measured force and
    log column 8's instantaneous SLS stiffness reference. The plant's explicit
    1.4 factor is used ONLY to interpret this existing scoring definition.
    These are diagnostic decompositions, not a deployable oracle method.
    """
    rows = []
    for fam, sid, params, duration in CE.cohort('diag', 1):
        log, _ = common.run(ServoID(single_model=False, **CE.FINAL), params, duration, sid)
        use = (log[:,0] >= 1.) & (log[:,3] > .5) & (log[:,6] > 0.)
        selected = log[use]
        fe = np.array([common.P.footprint(float(r[1]),float(r[6]),params)[0] for r in selected])
        ratio = selected[:,3]/fe
        # q_local=(dFe/dd)*d/Fe; replace dFe/dd using the LOGGED reference /1.4.
        # Combining total force and the instantaneous reference gives a different
        # apparent depth even if this local elastic exponent is supplied perfectly.
        apparent_depth_ratio = ratio/1.4
        row = dict(family=fam, sid=sid, scoring_samples=int(use.sum()),
                   total_over_elastic_force_P10_P50_P90=np.percentile(ratio,[10,50,90]).tolist(),
                   apparent_elastic_depth_ratio_P10_P50_P90=np.percentile(apparent_depth_ratio,[10,50,90]).tolist(),
                   explanation='scoring-only: static power-law inversion with total force and logged instantaneous reference; not a fitted candidate')
        rows.append(row)
        print('%s total/elastic F P10/P50/P90 %s; mixed static-depth ratio %s' %
              (fam,np.round(row['total_over_elastic_force_P10_P50_P90'],3).tolist(),
               np.round(row['apparent_elastic_depth_ratio_P10_P50_P90'],3).tolist()),flush=True)
    out = dict(protocol=dict(cohort='diag', scenarios_per_family=1, hidden_states_for_scoring_only=True,
                             controller_changed=False, plant_changed=False, final_used=False), rows=rows,
               source_hashes={str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__), HERE.parent/'cobot_scan_benchmark/plant.py']})
    dest = HERE/'results/normal_contact_state_diagnostic.json'
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(dest)


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--normal-state-diagnostic',action='store_true')
    args=ap.parse_args()
    normal_state_diagnostic() if args.normal_state_diagnostic else main()
