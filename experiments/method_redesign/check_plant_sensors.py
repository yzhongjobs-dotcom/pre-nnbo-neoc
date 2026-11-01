"""Bitwise regression against the pre-sensor plant, plus sensor error statistics.

python -X utf8 experiments/method_redesign/check_plant_sensors.py
Optional --baseline /path/to/plant_before.py; otherwise recover the pinned original
from Git in a temporary directory. No working-tree changes or Git writes.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'experiments/motion_aware_scan'))
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
import numpy as np
import common
import cohort_eval as C

BASE_REV = '727fa20a90230816641954e1edfaa423e5cb8928'


def check(baseline):
    spec = importlib.util.spec_from_file_location('plant_sensor_baseline', baseline)
    old = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = old
    spec.loader.exec_module(old)
    rows = []
    # Entire episodes, identical controller, parameters and normal-noise seeds.
    for fam, sid, p, duration in C.cohort('diag', 1):
        if fam not in ('S2', 'T3m', 'T2f-stiff'):
            continue
        before = old.episode(p, C.build('Proposed'), T=duration, seed=sid)
        after = common.P.episode(p, C.build('Proposed'), T=duration, seed=sid)
        assert before is not None and after is not None
        equal = np.array_equal(before.view(np.uint64), after[:, :12].copy().view(np.uint64))
        assert equal, (fam, np.max(abs(before-after[:, :12]), axis=0))
        r = dict(family=fam, seed=sid, ticks=len(after), bitwise_0_11=equal,
                 max_abs_0_11=np.max(abs(before-after[:, :12]), axis=0).tolist())
        for name, measured, truth in [('Ft', 12, 4), ('My', 13, 15)]:
            error = after[after[:, 0] >= 1., measured]-after[after[:, 0] >= 1., truth]
            r[name] = dict(bias=float(error.mean()), sd=float(error.std()),
                           rmse=float(np.sqrt(np.mean(error**2))),
                           abs_p95=float(np.quantile(abs(error), .95)))
        rows.append(r)
    return dict(baseline_revision=BASE_REV, sensor_errors='measured minus instantaneous truth; includes filter lag',
                rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as td:
        path = args.baseline
        if path is None:
            path = Path(td)/'plant_sensor_original.py'
            path.write_bytes(subprocess.check_output(
                ['git', 'show', BASE_REV+':experiments/cobot_scan_benchmark/plant.py'], cwd=ROOT))
        print(json.dumps(check(path), indent=2))
