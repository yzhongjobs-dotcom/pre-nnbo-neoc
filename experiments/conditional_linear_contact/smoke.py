"""默认参数试通：每族两个诊断场景，不计入 61 次调参。"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent / 'contact_information'))
from cohort_eval import FINAL, cohort  # noqa: E402
import common  # noqa: E402
from joint_policy import JointPolicy  # noqa: E402
from conditional_filter import ConditionalLinearContact, MeasuredOnly  # noqa: E402
from v1_servo_id import ServoID  # noqa: E402


def build(name, params=None):
    p = dict(FINAL)
    if params:
        p.update(params)
    if name == 'v1_proposed':
        return ServoID(single_model=False, **p)
    if name == 'joint_lean':
        return JointPolicy(wave=None, new_fit=True, rate=False,
                            amplitude_scale=1., **p)
    if name == 'conditional_v0':
        return ConditionalLinearContact(single_model=False, **p)
    raise ValueError(name)


def run_one(name, scene, params=None):
    fam, sid, plant_params, T = scene
    pol = build(name, params)
    log, _ = common.run(MeasuredOnly(pol), plant_params, T, sid)
    m = common.task_metrics(log, 20.8)
    return dict(method=name, family=fam, sid=sid, metrics=m,
                startup_equal=None, accepted=getattr(pol, 'accepted_count', None),
                valid=getattr(pol, 'valid_count', None),
                n_modes=len(getattr(pol, 'mode_history', [])),
                timing_ms_P50_P99=(np.percentile(pol.t_call[126:], [50, 99]) * 1000).tolist()
                if len(pol.t_call) > 130 else None)


def startup_check(scene):
    fam, sid, p, T = scene
    b = build('v1_proposed'); c = build('conditional_v0')
    lb, _ = common.run(MeasuredOnly(b), p, T, sid)
    lc, _ = common.run(MeasuredOnly(c), p, T, sid)
    n = min(len(lb), len(lc))
    mask = lb[:n, 0] < 1.0
    diffs = np.max(np.abs(lb[:n][mask, 11] - lc[:n][mask, 11])) if np.any(mask) else np.nan
    return dict(family=fam, sid=sid, rows=int(mask.sum()), max_command_diff=float(diffs),
                exact=bool(np.isfinite(diffs) and diffs == 0.))


def main():
    out = HERE / 'results'
    out.mkdir(exist_ok=True)
    scenes = cohort('diag', 2)
    rows = []
    for s in scenes:
        rows.append(startup_check(s))
        for name in ('v1_proposed', 'joint_lean', 'conditional_v0'):
            r = run_one(name, s)
            rows.append(r)
            print(r['family'], name, r['metrics'].get('rmse'), flush=True)
    payload = dict(protocol=dict(stage='default smoke', families=sorted({s[0] for s in scenes}),
                                 scenarios_per_family=2, seeds='diagnostic 110000+',
                                 confirmation_accessed=False, params=FINAL), rows=rows)
    (out / 'iteration_1_default_smoke.json').write_text(json.dumps(payload, indent=2, default=float), encoding='utf8')
    print(out / 'iteration_1_default_smoke.json')


if __name__ == '__main__':
    main()
