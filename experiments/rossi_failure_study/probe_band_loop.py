"""Closed loop: Huber-D4's window with only the probe-band part of the D4-filtered motion as weight (k = w'z / w'Y, w = P_B M Y).

Same window, preprocessing, controller and ungated definition as Huber-D4_nogate; only the estimate inside the window changes. When the band
gives a relative sd above 0.5 the window returns no value and the controller keeps its last stiffness. No tuning.
Predictions fixed before running: A_slow_thin 7 and combo 35 no longer run away (Oracle .23, .55); normal cases about unchanged; cough over-estimates unchanged.
Writes results/probe_band_loop.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import d_robust as DR
import memory_segment_eval as SE
import overlay_compare as O
from d_robust import common, _clip_k, _fit_record
from ridge_estimator import MIN_TICKS

P3, B, M = SE.P3, O.B, SE.M
HERE = Path(__file__).resolve().parent
CASES = (0, 7, 12, 24, 26, 29, 32, 35, 36, 38)
FREQS = (15.625, 20.833, 26.042)          # the 24-tick DFT bins 3, 4, 5
REL_MAX = .5
_CACHE = {}


def band_projector(m, Tc):
    key = ('P', m)
    if key not in _CACHE:
        t = np.arange(m)*Tc
        P = np.vander(np.arange(m)-(m-1)/2, 4)
        X = np.column_stack([f(2*np.pi*fr*t) for fr in FREQS for f in (np.cos, np.sin)])
        Qp, _ = np.linalg.qr(P); X = X-Qp@(Qp.T@X)
        U, s, _ = np.linalg.svd(X, full_matrices=False); U = U[:, s > 1e-8*s[0]]
        _CACHE[key] = U@U.T
    return _CACHE[key]


def weights(m, rho, Tc):
    if m not in _CACHE:
        t = np.arange(m)*Tc
        P = np.vander(np.arange(m)-(m-1)/2, 4)
        X = np.column_stack([f(2*np.pi*fr*t) for fr in FREQS for f in (np.cos, np.sin)])
        Qp, _ = np.linalg.qr(P); X = X-Qp@(Qp.T@X)
        U, s, _ = np.linalg.svd(X, full_matrices=False); U = U[:, s > 1e-8*s[0]]
        D = np.diff(np.eye(m), 4, axis=0)
        Md = np.eye(m)-np.linalg.inv(np.eye(m)+rho*D.T@D)
        _CACHE[m] = U@U.T@Md
    return _CACHE[m]


class BandEstimate:
    """Replaces RobustRidge.estimate on the same policy object; everything else in the host is unchanged.
    huber=True first runs the host's own Huber IRLS on the window and uses its final weights: w = P_B M_W Y."""
    def __init__(self, huber=False):
        self.huber = huber

    def __call__(self, est, o, n, policy):
        if self.huber:
            est.weights = None
            DR.RobustRidge.estimate(est, o, n, policy)
            hw = getattr(est, 'weights', None)
        p = policy.params; log = est._log(o, n); nw = est.nw
        est.fit = _fit_record(degree=0, tw=nw*policy.Tc)
        j0 = max(policy.ndel+1, n-1-nw); m = n-2-j0+1
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        te, ye = log[n-1, 0], log[n-1, 2]
        _, Y, z, keep = common._build(log, j0, m, 0, policy.ns, policy.eb, policy.afc, policy.ndel, te, ye, nw*policy.ns*1e-3, p['enc_y'])
        broken = np.flatnonzero(~keep)
        if len(broken):
            Y, z = Y[broken[-1]+1:], z[broken[-1]+1:]
        m = len(Y)
        if m < MIN_TICKS:
            return np.nan, np.inf, np.inf
        if self.huber and hw is not None and len(hw) == m-4:
            D = np.diff(np.eye(m), 4, axis=0)
            Mw = np.eye(m)-np.linalg.inv(np.eye(m)+est.rho*(D.T*hw)@D)
            w = band_projector(m, policy.Tc)@Mw@Y
        else:
            w = weights(m, est.rho, policy.Tc)@Y
        a = float(w@Y)
        if a <= 0:
            return np.nan, np.inf, np.inf
        k = float(w@z)/a
        sigma = est.sigma_scale*p['noise']*math.sqrt(float(w@w))/a
        if not np.isfinite(k) or k <= 0 or sigma > REL_MAX*k:
            return np.nan, np.inf, np.inf
        k = _clip_k(k)
        res = z-k*Y
        chi = float(res@res)/p['noise']**2/max(m-6, 1)
        est.fit = _fit_record(k, sigma, chi, degree=0, tw=nw*policy.Tc); est.fit['te'] = float(log[n-1, 0])
        return k, sigma, chi


def job(task):
    arm, index = task
    case = CC.cases()[index]
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    if arm in ('band', 'band_huber'):
        est = pol.estimator; f = BandEstimate(huber=arm == 'band_huber')
        est.estimate = lambda o, n, p: f(est, o, n, p)
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    return dict(arm=arm, index=index, group=case['group'], L=float(np.log1p(m['J'])), peak=float(m['peak']), rmse=float(m['rmse']), call_p50=float(m['call_p50']))


def main():
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, [(a, i) for a in ('huber', 'band', 'band_huber') for i in CASES], chunksize=1))
    stored = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    worst = max(abs(r['L']-stored[('Huber-D4_nogate', r['index'])]['L']) for r in rows if r['arm'] == 'huber')
    print('harness regression to stored Huber-D4_nogate, max |dL|', worst)
    get = {(r['arm'], r['index']): r for r in rows}
    print('idx group       | Huber-D4 ungated L peak | D4 + band L peak | Huber-D4 + band L peak | + IMM (stored) | Oracle (stored)')
    for i in CASES:
        h, b, bh = get[('huber', i)], get[('band', i)], get[('band_huber', i)]
        im, orc = stored[('Huber-D4_imm', i)], stored[('oracle_keff_off', i)]
        print(f"{i:3d} {h['group']:12s} | {h['L']:.3f} {h['peak']:6.1f} | {b['L']:.3f} {b['peak']:6.1f} | {bh['L']:.3f} {bh['peak']:6.1f} | {im['L']:.3f} {im['peak']:6.1f} | {orc['L']:.3f} {orc['peak']:6.1f}")
    print('call ms p50 band', np.median([r['call_p50'] for r in rows if r['arm'] == 'band'])/1e3, 'huber', np.median([r['call_p50'] for r in rows if r['arm'] == 'huber'])/1e3)
    (HERE/'results/probe_band_loop.json').write_text(json.dumps(dict(regression=worst, rows=rows), indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
