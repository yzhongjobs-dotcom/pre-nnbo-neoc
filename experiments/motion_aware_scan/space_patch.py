"""Replacement space(method, r) for the four comparator branches in tune_equal.

Keep tune_equal's existing Proposed branch. No tuning is run on import.
Rossi uses fixed displacement A_E, independent of the estimated stiffness.
Wij's measurement noise (5%) is fixed; its privileged probe rate is tunable.

DO-MPAC: search q_over_r=Q/R in [1e3,1e8], R=1, Np in {8,15,25,40}.
Reject candidates unless the ideally compensated A-BK is Schur-stable and
its dominant pole envelope time constant -Tc/log(rho) is in [0.04,2] s.
The upper limit leaves four time constants in the shortest 8 s episode;
the lower limit is five 8 ms ticks (avoids almost instantaneous nominal
dynamics). These are a priori task-scale screens, not test-scene tuning.
Positive Q,R alone do not certify finite-horizon receding-horizon stability.
Servo lag, delay, force LPF, ESO and saturation are absent from this check.
The class default q_over_r=20 preserves the pre-patch controller; it is
Schur-stable but deliberately excluded here because its nominal decay is
too slow. Defaults are diagnostic anchors, not sampled tuning candidates.
The physical model/ESO ranges are retained from tune_equal. In RossiDOMPAC,
k_h0 remains an interface parameter replaced immediately by Rossi k_hat;
it is included for parity but has no effect on that combination's control.

Integration caveat: tune_equal's stage-2 perturb() multiplies all floats
without respecting bounds. To retain these ranges in stage 2, use the same
range/pole checks there and resample invalid proposals. This file only
replaces space(), as requested; it does not alter that tuning procedure.
"""
import math
import numpy as np

from comparators import nominal_mpc


TC = .008
HORIZONS = (8, 15, 25, 40)
Q_OVER_R_BOUNDS = (1e3, 1e8)
TAU_BOUNDS = (.04, 2.)

# Yuen's 133 N/m nominal tissue stiffness through this task's stiff range.
# K_P/K_I are position-interface PI gains, not the paper's torque gains.
# Sensor privilege is fixed in YuenType and never sampled or perturbed.
YUEN_BOUNDS = dict(k_h0=(100., 5e4), K_P=(0., 1.5), K_I=(2., 120.),
                   q_scale=(1e-3, 1e3), refit_threshold=(3e-4, 5e-3))


def yuen_perturb(params, r):
    """Bounded stage-2 search; keep the sensor and integer harmonic count fixed/valid."""
    q = dict(params)
    for key, (lo, hi) in YUEN_BOUNDS.items():
        v = q[key]*math.exp(r.normal(0, .3))
        q[key] = float(np.clip(v, lo, hi))
    if r.random() < .3:
        q['harmonics'] = int(r.choice([1, 2, 3]))
    return q


def logu(r, lo, hi):
    return float(math.exp(r.uniform(math.log(lo), math.log(hi))))


def nominal_ok(q_over_r, Np):
    """Nominal PC125 pole screen; also used by smoke.py for verification."""
    _, poles = nominal_mpc(TC, Np, q_over_r)
    rho = float(np.max(np.abs(poles)))
    if not 0 < rho < 1:
        return False
    tau = -TC/math.log(rho)
    return TAU_BOUNDS[0] <= tau <= TAU_BOUNDS[1]


def _mpc_space(r):
    # Choose horizon first so rejection cannot bias the horizon frequencies.
    Np = int(r.choice(HORIZONS))
    for _ in range(10000):
        q_over_r = logu(r, *Q_OVER_R_BOUNDS)
        if nominal_ok(q_over_r, Np):
            return dict(m=logu(r, .3, 10.), b=logu(r, 50., 5000.),
                        k_h0=logu(r, 500., 5e4), eso_bandwidth=logu(r, 2., 60.),
                        q_over_r=q_over_r, Np=Np)
    raise RuntimeError('No admissible nominal MPC candidate; check bounds')


def space(method, r):
    """Sample a comparator's task-adapted parameters using numpy Generator r."""
    if method == 'YuenType':
        q = {k: (float(r.uniform(lo, hi)) if k == 'K_P' else logu(r, lo, hi))
             for k, (lo, hi) in YUEN_BOUNDS.items()}
        q['harmonics'] = int(r.choice([1, 2, 3]))
        return q
    if method == 'RossiType':
        return dict(K_I=logu(r, 2., 120.), lam=float(r.uniform(.95, .999)),
                    Q=logu(r, 1., 8.), A_E=logu(r, 2e-6, 1e-4),
                    f_exc=logu(r, 2., 30.))
    if method == 'WijType':
        mode = str(r.choice(['tanh', 'constant']))
        params = dict(K_P=float(r.uniform(0., 1.5)), K_I=logu(r, 2., 120.),
                      probe_hz=logu(r, 20., 125.), forget_mode=mode)
        if mode == 'tanh':
            params['innovation_scale'] = logu(r, .2, 5.)
        else:
            params['lam'] = float(r.uniform(.8, .999))
        return params
    if method == 'DOMPAC':
        return _mpc_space(r)
    if method == 'RossiDOMPAC':
        params = _mpc_space(r)
        params.update(lam=float(r.uniform(.95, .999)), bp_Q=logu(r, 1., 8.),
                      A_E=logu(r, 2e-6, 1e-4), f_exc=logu(r, 2., 30.))
        return params
    raise ValueError(method)
