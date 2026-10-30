"""Run from experiments/motion_aware_scan: python -X utf8 <codex_out>/smoke.py.

Eight default-controller episodes plus two static Rossi before/after episodes.
No tuning or file output. True force/stiffness columns are read only here for
evaluation (Wij's privileged probe is explicitly identified in comparators).
"""
import math
import sys
from pathlib import Path
from unittest.mock import patch

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(Path.cwd()))
sys.path.insert(0, str(OUT))

import numpy as np
# Keep real Numba JIT execution, but prevent cache creation beside repository
# sources. In restricted Windows sessions Numba's cache writeability probe may
# repeatedly retry blocked tempfile creation. This changes only disk caching.
with patch('numba.core.dispatcher.Dispatcher.enable_caching'):
    import common
import comparators as C
import space_patch as SP


class _RawEncoder:
    def update(self, o):
        return float(o.log[o.n-1, 2])


class BeforeRossi(C.RossiType):
    """Original raw-displacement regression; all other default behavior equal.

    The pre-patch default A_E was positive, so its unused a_force fallback had
    no effect. Replacing only this input reproduces that default controller.
    """
    def reset(self, platform):
        super().reset(platform)
        self.displacement = _RawEncoder()


def check_space():
    """Small parameter/pole checks, not additional plant experiments."""
    r = np.random.default_rng(71)
    for name in ('RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC'):
        modes = set()
        for _ in range(16):
            params = SP.space(name, r)
            c = getattr(C, name)(**params)
            c.reset(common.P.PLATFORMS['PC125'])
            if 'DOMPAC' in name:
                assert SP.nominal_ok(c.q_over_r, c.Np)
                assert c.R == 1. and np.max(np.abs(c.nominal_poles)) < 1.
            if name == 'WijType':
                modes.add(c.forget_mode)
        if name == 'WijType':
            assert modes == {'tanh', 'constant'}
    for n in SP.HORIZONS:
        admissible = []
        for q in np.geomspace(*SP.Q_OVER_R_BOUNDS, 81):
            if SP.nominal_ok(q, n):
                _, poles = C.nominal_mpc(SP.TC, n, q)
                rho = float(np.max(np.abs(poles)))
                assert rho < 1.
                admissible.append((q, rho, -SP.TC/math.log(rho)))
        assert admissible
        print(f'Nominal Np={n}: accepted grid q/r=[{admissible[0][0]:.6g},'
              f'{admissible[-1][0]:.6g}], rho=[{min(x[1] for x in admissible):.9f},'
              f'{max(x[1] for x in admissible):.9f}], '
              f'tau=[{min(x[2] for x in admissible):.6f},'
              f'{max(x[2] for x in admissible):.6f}] s', flush=True)
    c = C.DOMPAC()
    c.reset(common.P.PLATFORMS['PC125'])
    rho = float(np.max(np.abs(c.nominal_poles)))
    assert rho < 1.
    print(f'Legacy default MPC q/r={c.q_over_r:g}, Np={c.Np}, poles={c.nominal_poles}, '
          f'tau={-SP.TC/math.log(rho):.6f} s (outside new speed screen)', flush=True)


def run_one(cls, p, T, seed, label):
    log, controller = common.run(cls(), p.copy(), T, seed)
    if log is None:
        print(f'{label:26s} RMSE=nan N peak=nan N n_loss=nan DIVERGED', flush=True)
        return log, controller
    assert np.isfinite(log).all(), f'{label}: nonfinite episode'
    m = common.task_metrics(log, getattr(controller, 'f_exc', None))
    assert np.isfinite(m['rmse']) and np.isfinite(m['peak'])
    print(f"{label:26s} RMSE={m['rmse']:.6f} N peak={m['peak']:.6f} N "
          f"n_loss={m['n_loss']}", flush=True)
    return log, controller


def main():
    assert Path(C.__file__).resolve() == OUT/'comparators.py'
    print(f'Comparator module: {C.__file__}', flush=True)
    print('Metrics: true contact force, t>=1 s; peak is max force; n_loss counts '
          'entries below 0.5 N. Static ratio uses t>=2 s.', flush=True)
    check_space()
    failures = []
    for task, seed in [('T2', 11), ('T3', 12)]:
        dose = {'stiff': False} if task == 'T2' else None
        p, T = common.S.make(task, np.random.default_rng(seed), dose=dose)
        print(f'{task}: scene seed={seed}, force-noise seed={seed}, '
              f'T={T:g} s, E_soft={p[0]:.6g} Pa', flush=True)
        for cls in (C.RossiType, C.WijType, C.DOMPAC, C.RossiDOMPAC):
            log, _ = run_one(cls, p, T, seed, f'{task} {cls.name}')
            if log is None:
                failures.append(f'{task} {cls.name}')

    # Unmodified T2 seed-11 draw (hard tissue); remove surface motion and scan.
    p, T = common.S.make('T2', np.random.default_rng(11))
    p[2] = 0.; p[15] = 0.; p[5] = 0.
    assert p[6] == 0  # no spatial stiffness ridges
    print('Static homogeneous tissue: T2 seed=11, p[2]=p[15]=p[5]=0; '
          'force noise, SLS relaxation and contact nonlinearity retained.', flush=True)
    ratios = []
    for label, cls in [('before', BeforeRossi), ('after', C.RossiType)]:
        log, c = run_one(cls, p, T, 11, f'Static Rossi {label}')
        trace = np.asarray(c.trace)
        assert len(trace) == len(log) and np.allclose(trace[:, 0], log[:, 0])
        mask = (trace[:, 0] >= 2.) & (log[:, 8] > 0.)
        ratio = float(np.median(trace[mask, 1]/log[mask, 8]))
        ratios.append(ratio)
        print(f'Static {label}: median(k_hat/k_eff)={ratio:.6f}; '
              f'median k_hat={np.median(trace[mask,1]):.6f}, '
              f'median k_eff={np.median(log[mask,8]):.6f} N/m', flush=True)
    assert abs(ratios[1]-1.) < .1, f'Matched static estimate not near 1: {ratios[1]}'
    assert abs(ratios[1]-1.) < abs(ratios[0]-1.), 'Preprocessing did not improve bias'
    print('PASS: both forgetting spaces, screened nominal poles, '
          'and improved static Rossi estimate.', flush=True)
    if failures:
        print('FAIL: divergent default episodes: '+', '.join(failures), flush=True)
        return 1
    print('PASS: all eight default episodes finite.', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
