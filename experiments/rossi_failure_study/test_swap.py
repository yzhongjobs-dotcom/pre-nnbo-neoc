"""Regression: Projection(16,3,1,1) must be command-identical to ServoID."""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for s in ('motion_aware_scan', 'method_redesign'):
    sys.path.insert(0, str(HERE.parent/s))
sys.path.insert(0, str(HERE))

import physics_plant as PP
from cohort_eval import FINAL
from v1_servo_id import ServoID
from estimator_swap import SwapPolicy
from estimators import ProjectionEstimator
from failure_replay import Measured, cell_scene


def test_projection_is_servoid():
    sc = cell_scene('shift_fast', 'mid', 2, 341234)
    seed = 341234 + 2
    a = Measured(ServoID(single_model=False, **FINAL))
    b = Measured(SwapPolicy(ProjectionEstimator(nw=16, degree=3, sigma_scale=1., chi_scale=1.),
                            single_model=False, **FINAL))
    la = PP.episode(sc, a, seed)
    lb = PP.episode(sc, b, seed)
    assert la is not None and lb is not None
    n = min(len(la), len(lb))
    diff = np.max(np.abs(la[:n, 11]-lb[:n, 11]))
    assert diff == 0.0, f'max command difference = {diff:.17g}'


if __name__ == '__main__':
    test_projection_is_servoid()
    print('PASS test_projection_is_servoid max_command_difference=0')
