"""Six paired closed-loop checks of ordinary least-squares stiffness fits.

Reuses corrected_comparison's frozen scenes, seeds, scoring and all 17 stored
method rows. Only Simple-1 (z=kY) and Simple-2 (z=kY+b) are simulated.
Run: python -X utf8 experiments/rossi_failure_study/simple_ls_validation.py
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import corrected_comparison as CC
import overlay_compare as O
import rib_transition_eval as R
import memory_segment_eval as SE
import common
import id_baselines as I
from estimators import _clip_k, _fit_record
from ridge_estimator import MIN_TICKS

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results/simple_ls_validation.json'
METHODS = ('Simple-1', 'Simple-2')
SCENES = ('A', 'B', 'combo')
SOURCE_METHODS = CC.METHODS


class SimpleLS:
    """OLS on the same compensated, H-preprocessed equation data as Huber-D4."""
    def __init__(self, intercept, nw=24, sigma_scale=1., chi_scale=1.):
        self.intercept = bool(intercept)
        self.nw = int(nw)
        self.sigma_scale = float(sigma_scale)
        self.chi_scale = float(chi_scale)
        self.dof = max(self.nw - int(self.intercept) - 1, 1)

    def reset(self, policy):
        self.fit = _fit_record(degree=0, tw=self.nw * policy.Tc)

    def estimate(self, o, n, policy):
        p = policy.params
        log = o.log
        self.fit = _fit_record(degree=0, tw=self.nw * policy.Tc)
        j0 = max(policy.ndel + 1, n - 1 - self.nw)
        m = n - 2 - j0 + 1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        X, Y, z, keep = common._build(
            log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel,
            te, ye, self.nw * policy.ns * 1e-3, p['enc_y'])
        # Match the selected estimators' contiguous-contact rule.
        broken = np.flatnonzero(~keep)
        if len(broken):
            X, Y, z = X[broken[-1]+1:], Y[broken[-1]+1:], z[broken[-1]+1:]
        m = len(Y)
        npar = 2 if self.intercept else 1
        if m < max(MIN_TICKS, npar + 1):
            return np.nan, np.inf, np.inf
        basis = X[:, 0] if self.intercept else None
        if self.intercept:
            # The filtered constant column is the intercept's basis in z's
            # measurement domain; this retains force-filter/servo compensation.
            A = np.column_stack((Y, basis))
            gram = A.T @ A
            if not np.isfinite(gram).all() or np.linalg.cond(gram) > 1e12:
                return np.nan, np.inf, np.inf
            beta = np.linalg.solve(gram, A.T @ z)
            k_raw, b_raw = float(beta[0]), float(beta[1])
            residual = z - A @ beta
            info = float(Y @ Y - (Y @ basis)**2 / (basis @ basis))
        else:
            info = float(Y @ Y)
            if not np.isfinite(info) or info <= 0:
                return np.nan, np.inf, np.inf
            k_raw = float(Y @ z / info)
            b_raw = 0.
            residual = z - k_raw * Y
        if not np.isfinite(info) or info <= 0:
            return np.nan, np.inf, np.inf
        dof = max(m - npar, 1)
        sigma = self.sigma_scale * p['noise'] / math.sqrt(info)
        chi = self.chi_scale * float(residual @ residual) / (p['noise']**2 * dof)
        k = _clip_k(k_raw)
        self.dof = dof
        self.fit = _fit_record(k, sigma, chi, degree=0, tw=self.nw * policy.Tc)
        self.fit.update(raw_k=k_raw, b=b_raw, info=1./(sigma*sigma), chi2=chi,
                        valid=m, te=float(log[n-1, 0]))
        return k, sigma, chi


def synthetic_checks():
    # Exact known-coefficient checks for both models, including a filtered
    # constant basis rather than a raw-domain force/encoder offset.
    Y = np.linspace(-.01, .02, 40)
    basis = np.linspace(.91, 1., 40)
    k0, b0 = 31750., .42
    z1 = k0 * Y
    z2 = k0 * Y + b0 * basis
    k1 = float(Y @ z1 / (Y @ Y))
    A = np.column_stack((Y, basis))
    k2, b2 = np.linalg.lstsq(A, z2, rcond=None)[0]
    assert abs(k1-k0) < 1e-9 and abs(k2-k0) < 1e-8 and abs(b2-b0) < 1e-10
    return dict(passed=True, method1_k=k1, method2_k=k2,
                method2_measurement_domain_b=b2, tolerance='absolute')


def chosen_cases():
    all_cases = CC.cases()
    out = [c for c in all_cases if c['group'] in SCENES and
           tuple(c['cell'][2:]) == ('neural_white', 28., 0)]
    assert [c['group'] for c in out] == list(SCENES)
    assert [c['index'] for c in out] == [0, 12, 30]
    return out


def run_one(method, case):
    O.setup()
    cell = tuple(case['cell'])
    sc = SE.MC.scene_for(cell, CC.BASE, case['variant'])
    host_cfg = O.chosen()['Huber-D4']
    # Start from the exact existing host policy and retain its method/controller
    # settings. The adapter supplies the common H preprocessing and phase path.
    pol = O.B.build('Huber-D4', host_cfg, sc)
    estimator = SimpleLS(method == 'Simple-2', nw=24,
                         sigma_scale=host_cfg['sigma_scale'],
                         chi_scale=host_cfg['gate'])
    pol.estimator = I.ContactAdapter(estimator)
    assert pol.estimator.nw == 24 and pol.estimator.inner.nw == 24
    assert pol.estimator.inner.sigma_scale == host_cfg['sigma_scale']
    assert pol.estimator.inner.chi_scale == host_cfg['gate']
    assert pol.acquire == O.B.build('Huber-D4', host_cfg, sc).acquire
    assert pol.centre == O.B.build('Huber-D4', host_cfg, sc).centre
    estimates = []
    step = pol.step

    def capture(o):
        command = step(o)
        fit = pol.estimator.fit
        estimates.append((float(o.log[o.n-1, 0]), float(fit.get('raw_k', np.nan)),
                          float(fit.get('b', np.nan)), float(pol.k_good), float(pol.kbar)))
        return command

    pol.step = capture
    log = SE.P3.episode(sc, pol, case['seed'])
    if log is None:
        raise RuntimeError(('plant failed', method, case))
    metric = SE.M.metric(log, pol, sc)
    rec = np.asarray(estimates)
    assert len(rec) == len(log)
    truth = np.asarray(log[:, 26], float)  # instantaneous axial force/encoder gain including arm and body support
    valid = np.isfinite(rec[:, 1]) & np.isfinite(truth) & (truth > 0.)
    err = rec[valid, 1] - truth[valid]
    raw = rec[valid, 1]
    row = dict(method=method, config=method, **case,
        L=float(np.log1p(metric['J'])), force_rmse=float(metric['rmse']),
        peak=float(metric['peak']), J=float(metric['J']),
        tool_mass=sc['tool_mass'], tool_mass_hat=sc['tool_mass_hat'],
        corrected_preprocessing=True, H_preprocessing='ContactAdapter( H )',
        window_ticks=24, phase_estimator='same ContactAdapter phase estimator',
        controller='same Huber-D4 host/controller; identical gates and thresholds',
        k_truth='log[:,26] instantaneous k_eff', k_valid_count=int(valid.sum()),
        k_raw_rmse=float(np.sqrt(np.mean(err**2))) if len(err) else None,
        k_raw_mae=float(np.mean(np.abs(err))) if len(err) else None,
        k_raw_bias=float(np.mean(err)) if len(err) else None,
        k_raw_underestimate_fraction=float(np.mean(raw < truth[valid])) if len(err) else None,
        k_raw_underestimate_mean_magnitude=float(np.mean(np.maximum(-err, 0.))) if len(err) else None,
        k_raw_median=float(np.median(raw)) if len(err) else None,
        offset_estimate_median=float(np.nanmedian(rec[np.isfinite(rec[:, 2]), 2]))
            if np.isfinite(rec[:, 2]).any() else None,
        gate_scales=dict(sigma_scale=host_cfg['sigma_scale'], chi_scale=host_cfg['gate']),
        )
    return row


def main():
    checks = synthetic_checks()
    source = json.loads((HERE/'results/corrected_comparison.json').read_text(encoding='utf-8'))
    assert source['validation']['passed'] and source['validation']['methods'] == 17
    cases = chosen_cases()
    episodes = []
    for method in METHODS:
        for case in cases:
            episodes.append(run_one(method, case))
            print(method, case['group'], 'done', flush=True)
    existing = {(r['method'], r['index']): r for r in source['episodes']}
    assert all((m, c['index']) in existing for m in SOURCE_METHODS for c in cases)
    table = []
    for case in cases:
        for method in SOURCE_METHODS + METHODS:
            row = existing[(method, case['index'])] if method in SOURCE_METHODS else next(
                r for r in episodes if r['method'] == method and r['index'] == case['index'])
            table.append(dict(group=case['group'], index=case['index'], cell=case['cell'],
                              seed=case['seed'], method=method, L=row['L'],
                              force_rmse=row['force_rmse'], peak=row['peak']))
    data = dict(metadata=dict(source='corrected_comparison.json', reused_methods=list(SOURCE_METHODS),
        new_methods=list(METHODS), simulated_episodes=6, selected_cases=cases,
        fixed_before_results=True, selection='A/B/combo; neural_white; 28 degrees; i=0',
        common_settings='same corrected H preprocessing, 24-tick window, phase estimator, gates, controller and parameter settings',
        force_metrics='L=ln(1+J), force RMSE, peak true force; copied from original episodes for existing methods',
        stiffness_diagnostics='Unclipped OLS solution versus instantaneous log[:,26] k_eff, including arm and body support; active controller output is clipped by the existing _clip_k bounds. Raw fit errors are kept explicit so clipping does not conceal estimator excursions.'),
        solver_checks=checks, new_method_episodes=episodes, comparison_table=table,
        validation=dict(passed=True, existing_methods_reused=17, new_methods=2,
                        new_simulations=len(episodes), selected_scenes=3))
    def json_safe(value):
        if isinstance(value, dict):
            return {k: json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(v) for v in value]
        if isinstance(value, (float, np.floating)):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, (np.integer,)):
            return int(value)
        return value
    OUT.write_text(json.dumps(json_safe(data), ensure_ascii=False, allow_nan=False, indent=2), encoding='utf-8')
    print('saved', OUT, flush=True)
    for r in table:
        print(f"{r['group']:5s} {r['method']:14s} L={r['L']:.6f} RMSE={r['force_rmse']:.6f} peak={r['peak']:.6f}")


if __name__ == '__main__':
    main()
