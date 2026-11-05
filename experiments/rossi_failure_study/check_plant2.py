"""Regression and mechanism checks for physics_plant2 (no controller ranking)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'cobot_scan_benchmark'))
sys.path.insert(0, str(HERE))
import physics_plant as P1
import physics_plant2 as P2
import scenes as S1
import scenes2 as S2


class Hold:
    direct = True

    def reset(self, platform):
        pass

    def __call__(self, o):
        return float(o.log[o.n - 1, 11])


def zero_scene(T=4.):
    sc = S1.make('rest', 'mid', np.random.default_rng(4201), T=T, v=.02, overrides={'field_sigma': 0.})
    sc['h_s'][:] = 0.
    out = dict(sc)
    n = int(round(T * 1000.)) + 2
    out['xcmd'] = P2.X0 + sc['v'] * np.maximum(np.arange(n) * .001 - P2.PLATFORM['ndel'] * P2.PLATFORM['Tc'], 0.)
    out.update(mu_k=0., mu_s=0., k_d=0., k1=0., k2=0., theta0=0., c_slope=0., xi_c=0., K_t=1e12, m_t=1., zeta_t=.5)
    return sc, out


def held_commands(sc):
    log = P1.episode(sc, Hold(), 4202, sensor=True)
    return log, log[:, 11]


def physical_case(regime='nominal', T=4.):
    sc = S2.cell_scene2('rest', 'mid', regime, 0, base=420000, T=T)
    s, buf = P2._init(sc, P2.PLATFORM)
    return sc, np.full(int(round(T / P2.PLATFORM['Tc'])), buf[-1])


def main():
    out = {}
    old, new = zero_scene()
    old_log, commands = held_commands(old)
    new_log = P2.replay(new, commands, seed=4202, sensor=True)
    diffs = [float(np.max(np.abs(old_log[:, j] - new_log[:, j]))) for j in (3, 5, 8)]
    out['zero_physics_replay_max_diff_columns_3_5_8'] = diffs
    print('zero-new-physics replay max diff: col3 %.9g, col5 %.9g, col8 %.9g' % tuple(diffs))
    if max(diffs) >= 1e-6:
        print('REGRESSION FAIL: at least one requested difference is not below 1e-6')
    else:
        print('REGRESSION PASS: all requested differences < 1e-6')

    sc, cmd = physical_case('nominal')
    sc['mu_s'] = sc['mu_k']; sc['k_d'] = 0.; sc['k1'] = 0.; sc['k2'] = 0.; sc['theta0'] = 0.; sc['c_slope'] = 0.; sc['xi_c'] = 0.
    log = P2.replay(sc, cmd, seed=4203, sensor=False)
    warm = (log[:, 0] > 1.) & (log[:, 3] > .5) & (np.abs(log[:, 16]) > .005)
    ratio = np.abs(log[warm, 19] / np.maximum(log[warm, 3], 1e-9))
    print('uniform sliding |Ft|/Fn median %.6g; set mu_k %.6g; n=%d' % (np.median(ratio), sc['mu_k'], ratio.size))

    rev, rev_cmd = physical_case('reversal')
    rev['mu_s'] = rev['mu_k']; rev['k_d'] = 0.; rev['k1'] = 0.; rev['k2'] = 0.; rev['theta0'] = 0.; rev['c_slope'] = 0.; rev['xi_c'] = 0.
    rlog = P2.replay(rev, rev_cmd, seed=4204, sensor=False)
    sign = np.sign(rlog[:, 16]); flips = np.flatnonzero(sign[1:] * sign[:-1] < 0)
    jumps = [abs(rlog[i + 1, 19] - rlog[i, 19]) for i in flips if rlog[i, 3] > .5]
    print('reversal Fx truth jumps: median %.6g, expected 2 mu F median %.6g, n=%d' % (np.median(jumps) if jumps else np.nan, 2 * rev['mu_k'] * np.median(rlog[:, 3]), len(jumps)))

    stick, stick_cmd = physical_case('nominal')
    stick['K_t'] = 10.
    slog = P2.replay(stick, stick_cmd, seed=4205, sensor=False)
    sv = slog[:, 16]; st = (slog[:, 0] > 1.) & (np.abs(sv) > 1e-5)
    changes = int(np.sum(np.sign(sv[1:]) * np.sign(sv[:-1]) < 0))
    q10, q90 = np.percentile(np.abs(sv[st]), [10, 90]); stick_frac = float(np.mean(np.abs(sv[st]) < .001))
    print('low-Kt stick-slip: vc range [%.6g, %.6g], |vc| P10/P90 %.6g/%.6g, near-stick fraction %.3f, sign changes %d, mu_s/mu_k %.3f' % (np.min(sv[st]), np.max(sv[st]), q10, q90, stick_frac, changes, stick['mu_s'] / stick['mu_k']))

    curved, curved_cmd = physical_case('curved')
    clog = P2.replay(curved, curved_cmd, seed=4206, sensor=False)
    m = (clog[:, 3] > .5) & (clog[:, 0] > 1.)
    err_exact = clog[m, 18] - clog[m, 3] * np.cos(clog[m, 17]) - clog[m, 19] * np.sin(clog[m, 17])
    err_small_angle = clog[m, 18] - clog[m, 3] - clog[m, 19] * np.sin(clog[m, 17])
    print('F_axis exact decomposition residual: median %.6g, max abs %.6g; small-angle residual median %.6g' % (np.median(err_exact), np.max(np.abs(err_exact)), np.median(err_small_angle)))

    asym, asym_cmd = physical_case('nominal')
    asym['xcmd'] = P2.X0 + .03 * np.arange(len(asym['xcmd'])) * .001
    asym['k1'] = .001 * 3.; asym['k2'] = .003 * 3.; asym['mu_k'] = asym['mu_s'] = 0.; asym['k_d'] = 0.; asym['theta0'] = asym['c_slope'] = asym['xi_c'] = 0.
    alog = P2.replay(asym, asym_cmd, seed=4207, sensor=False)
    mm = (alog[:, 0] > 1.) & (alog[:, 3] > .5)
    d = P2._fe_asym(P2.X0, .002, .03, asym['ln_e'], asym['dx'], asym['k1'], asym['k2'])
    d0 = P2._fe(P2.X0, .002, asym['ln_e'], asym['dx'])
    a0 = P2.halfwidth(.002)
    expected_ratio = 1. + (asym['k2'] - asym['k1']) * .03 / (2. * a0)
    observed_ratio = d / d0
    print('v=.03 pressure-asymmetry F(v)/F(0): observed %.9g, eq(14) %.9g, ratio %.9g' % (observed_ratio, expected_ratio, observed_ratio / expected_ratio))
    print('2 kHz finite-state check:', bool(np.all(np.isfinite(alog))), 'rows', len(alog))
    out['physical'] = {'uniform_ratio_median': float(np.median(ratio)), 'reversal_jump_n': len(jumps), 'stick_slip_sign_changes': changes,
                       'projection_error_max': float(np.max(np.abs(err_exact))), 'finite': bool(np.all(np.isfinite(alog)))}
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results' / 'plant2_checks.json').write_text(json.dumps(out, indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
