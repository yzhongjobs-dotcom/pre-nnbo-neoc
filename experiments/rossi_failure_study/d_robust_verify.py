"""Independent dense least-squares audit and paired closed-loop timing/equivalence.

Run on Windows: python -X utf8 experiments/rossi_failure_study/d_robust_verify.py
Prints JSON to stdout; writes no result files. Times are microseconds. Paired
calls alternate order and use separate estimator states on the same policy
trajectory. Clean episodes measure method3.metric without paired-call overhead.
Identical-window equivalence does not guarantee identical feedback trajectories;
the independent runs also report their force and command divergence.
The standard Huber convention is h(u)=u^2/2 inside [-c,c]. The implementation
targets residual^2 + 2*rho*s^2*sum h(Dd/s), not residual^2 + rho*sum h(Dd/s).
"""
import copy
import json
import platform
import sys
import time
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import scipy

import sweeps  # establishes the existing experiments' shared import paths
import d_robust as R
import method3 as M
import physics_plant3 as P3
import scenes3 as S3


def percentiles(values):
    return dict(zip(('p50', 'p99'), map(float, np.percentile(values, [50, 99]))))


def compare(a, b):
    """Reference-relative errors, including explicit NaN/Inf pattern checks."""
    a, b = np.asarray(a), np.asarray(b)
    if not (np.array_equal(np.isnan(a), np.isnan(b))
            and np.array_equal(np.isposinf(a), np.isposinf(b))
            and np.array_equal(np.isneginf(a), np.isneginf(b))):
        raise AssertionError(f'Nonfinite pattern mismatch: {a}, {b}')
    valid = np.isfinite(a) & np.isfinite(b)
    rel = np.zeros(a.shape)
    rel[valid] = np.abs(a[valid]-b[valid])/np.maximum(np.abs(b[valid]), 1e-12)
    if np.any(rel > 1e-6):
        raise AssertionError(f'Relative error exceeds 1e-6: {a}, {b}, {rel}')
    return rel


def joint_lstsq(Y, z, D, rho, noise, c, iters):
    """Independent IRLS oracle: solve for (k,d) jointly via augmented dense SVD.

    No inverse, band packing, Schur complement, or production solver is used.
    Final sigma/chi are obtained from augmented residual norms, not I-S.
    """
    m = len(Y); w = np.ones(len(D)); s = noise/np.sqrt(rho)
    for it in range(iters):
        penalty = np.sqrt(rho*w)[:, None]*D
        A = np.block([[Y[:, None], np.eye(m)],
                      [np.zeros((len(D), 1)), penalty]])
        k, *d = np.linalg.lstsq(A, np.r_[z, np.zeros(len(D))], rcond=None)[0]
        if it != iters-1:
            u = np.abs(D@d)/s
            w = np.minimum(1., c/np.maximum(u, 1e-12))
    B = np.vstack((np.eye(m), penalty))
    targets = np.column_stack((Y, z-k*Y))
    targets = np.vstack((targets, np.zeros((len(D), 2))))
    residual = targets-B@np.linalg.lstsq(B, targets, rcond=None)[0]
    info, energy = np.sum(residual**2, axis=0)
    return (R._clip_k(k), noise/np.sqrt(info), energy/noise**2/(m-5)), w


def algebra_audit():
    rng = np.random.default_rng(90210); max_rel = np.zeros(3); max_w = 0.
    max_band_error = 0.; cases = 0
    for fc in (8., 14., 24.):
        for m in range(8, 25):
            policy = SimpleNamespace(Tc=.008, ndel=0, ns=8, eb=.8, afc=.1,
                                     params=dict(noise=.02, enc_y=1e-6))
            est = R.RobustRidge(nw=24, fc=fc); est.reset(policy)
            D, basis = est._difference_bands(m)
            w = rng.uniform(.01, 1., len(D))
            ab = (basis@w).reshape((5, m))
            dense = np.zeros((m, m))
            for q in range(5):
                j = np.arange(m-q)
                dense[j+q, j] = ab[q, :m-q]
                dense[j, j+q] = ab[q, :m-q]
            max_band_error = max(max_band_error, float(np.max(np.abs(dense-est.rho*(D.T*w)@D))))
            np.testing.assert_allclose(dense, est.rho*(D.T*w)@D, rtol=2e-14, atol=1e-10)
            Y = .001*rng.normal(size=m)
            z = 2000.*Y+.02*rng.normal(size=m)
            z[m//2] += .6
            obs = SimpleNamespace(log=np.zeros((m+2, 14)))
            with patch.object(R.common, '_build', return_value=(None, Y, z, np.ones(m, bool))):
                fast = est.estimate(obs, m+2, policy)
            oracle, weights = joint_lstsq(Y, z, D, est.rho, .02, est.c, est.iters)
            max_rel = np.maximum(max_rel, compare(fast, oracle))
            max_w = max(max_w, float(np.max(np.abs(est.weights-weights))))
            cases += 1
    # Startup and a broken contact leaving fewer than MIN_TICKS.
    for n, keep in ((4, np.ones(1, bool)), (14, np.r_[np.zeros(8, bool), np.ones(4, bool)])):
        a = R.RobustRidge(nw=24); b = R.RobustRidge(nw=24)
        a.reset(policy); b.reset(policy)
        obs = SimpleNamespace(log=np.zeros((n, 14)))
        with patch.object(R.common, '_build', return_value=(None, np.ones(len(keep)), np.ones(len(keep)), keep)):
            compare(a.estimate(obs, n, policy), b.estimate_reference(obs, n, policy))
    return dict(cases=cases, max_relative=dict(zip(('k', 'sigma', 'chi'), max_rel.tolist())),
                max_weight_absolute=max_w, max_band_absolute=max_band_error,
                startup_and_broken_contact='passed')


def paired_episode(sc, seed):
    pol = R.build(sc); fast = pol.estimator; ref = copy.deepcopy(fast)
    fast_reset = fast.reset; fast_estimate = fast.estimate
    times = dict(fast=[], reference=[]); errors = np.zeros(3)
    counts = dict(calls=0, finite=0); max_weights = 0.; lengths = set()

    def reset(policy):
        fast_reset(policy); ref.reset(policy)

    def estimate(o, n, policy):
        nonlocal errors, max_weights
        results = {}
        order = ('fast', 'reference') if counts['calls'] % 2 == 0 else ('reference', 'fast')
        for name in order:
            fn = fast_estimate if name == 'fast' else ref.estimate_reference
            t0 = time.perf_counter_ns()
            results[name] = fn(o, n, policy)
            times[name].append((time.perf_counter_ns()-t0)/1000.)
        errors = np.maximum(errors, compare(results['fast'], results['reference']))
        counts['calls'] += 1
        if np.all(np.isfinite(results['fast'])):
            counts['finite'] += 1
            lengths.add(len(fast.weights)+fast.order)
            max_weights = max(max_weights, float(np.max(np.abs(fast.weights-ref.weights))))
            np.testing.assert_allclose(fast.fit['te'], ref.fit['te'], rtol=0, atol=0)
            if fast.dof != ref.dof:
                raise AssertionError('dof mismatch')
        return results['fast']

    fast.reset = reset; fast.estimate = estimate
    log = P3.episode(sc, pol, seed)
    if log is None or not counts['finite']:
        raise AssertionError('Closed loop failed or no finite estimates')
    row = dict(**counts, window_lengths=sorted(lengths),
               max_relative=dict(zip(('k', 'sigma', 'chi'), errors.tolist())),
               max_weight_absolute=max_weights,
               estimator_us={name: percentiles(ts) for name, ts in times.items()})
    return row, times


def clean_episode(sc, seed, reference):
    pol = R.build(sc)
    if reference:
        pol.estimator.estimate = pol.estimator.estimate_reference
    log = P3.episode(sc, pol, seed)
    if log is None:
        raise AssertionError('Clean closed loop failed')
    metric = M.metric(log, pol, sc)
    return {key: metric[key] for key in ('call_p50', 'call_p99', 'J', 'rmse')}, log


def main():
    t0 = time.perf_counter()
    print(json.dumps(dict(environment=dict(python=sys.version, numpy=np.__version__,
                                          scipy=scipy.__version__, platform=platform.platform()),
                          algebra_audit=algebra_audit())), flush=True)
    pooled = dict(fast=[], reference=[]); rows = []
    for i in range(3):
        sc = S3.cell_scene3('combo', 'mixed', 'nominal', i, base=520500,
                            overrides=dict(slide='long', force=8., dynamics='neural_white', phase=28.))
        # Same sensor-noise seed in paired and clean runs, matching d_prior_check.
        row, times = paired_episode(sc, 900+i)
        row['i'] = i; row['seed'] = 900+i
        row['strategy_metrics'] = {}; clean_logs = {}
        for reference in ((False, True) if i % 2 == 0 else (True, False)):
            name = 'reference' if reference else 'fast'
            row['strategy_metrics'][name], clean_logs[name] = clean_episode(sc, 900+i, reference)
        row['independent_trajectory_absolute'] = {}
        for name, col in (('command_m', 11), ('force_N', 18)):
            delta = np.abs(clean_logs['fast'][:, col]-clean_logs['reference'][:, col])
            first = np.flatnonzero(delta > 1e-6)
            row['independent_trajectory_absolute'][name] = dict(
                maximum=float(delta.max()),
                first_above_1e_6_s=float(clean_logs['fast'][first[0], 0]) if len(first) else None)
        for name in pooled:
            pooled[name].extend(times[name])
        rows.append(row)
        print(json.dumps(row), flush=True)
    print(json.dumps(dict(summary=dict(calls=sum(r['calls'] for r in rows),
                                       finite=sum(r['finite'] for r in rows),
                                       max_relative={k: max(r['max_relative'][k] for r in rows)
                                                     for k in ('k', 'sigma', 'chi')},
                                       estimator_us={name: percentiles(ts) for name, ts in pooled.items()},
                                       elapsed_s=time.perf_counter()-t0))), flush=True)


if __name__ == '__main__':
    main()
