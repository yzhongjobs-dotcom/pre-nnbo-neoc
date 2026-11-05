"""Candidate regressions and one-scene mechanism checks."""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'cobot_scan_benchmark'))
import contact_candidates as C
import contact_compare as CC
import physics_plant as P1
import physics_plant2 as P2
import scenes as S1
import scenes2 as S2


class FullMeasured(C.Measured2):
    def __call__(self, o):
        a = o.log[:o.n].copy()
        return self.inner(P2.OLD.Obs(a, o.n))


def zero_scene(T=3.):
    sc = S1.make('rest', 'mid', np.random.default_rng(4311), T=T, v=.02, overrides={'field_sigma': 0.})
    sc['h_s'][:] = 0.
    out = dict(sc)
    out['xcmd'] = P2.X0 + sc['v'] * np.maximum(np.arange(int(T * 1000) + 2) * .001 - P2.PLATFORM['ndel'] * P2.PLATFORM['Tc'], 0.)
    out.update(mu_k=0., mu_s=0., k_d=0., k1=0., k2=0., theta0=0., c_slope=0., xi_c=0., K_t=1e12, m_t=1., zeta_t=.5)
    return sc, out


def old_measured(inner):
    class M:
        direct = True
        def reset(self, p): inner.reset(dict(P1.PLATFORM))
        def __call__(self, o):
            a = np.full((o.n, 16), np.nan); a[:, [0, 1, 2, 5, 11, 12, 13]] = o.log[:o.n, [0, 1, 2, 5, 11, 12, 13]]
            return inner(P1.OLD.Obs(a, o.n))
    return M()


def main():
    out = {'regressions': {}, 'mechanisms': {}}
    old, new = zero_scene()
    hp = CC.build('elastic'); cp = C.build('L')
    a = P1.episode(old, old_measured(hp), 4312, sensor=False)
    b = P2.episode(new, C.Measured2(cp), 4312, sensor=False)
    command_gap = float(np.max(np.abs(a[:, 11] - b[:, 11])))
    out['regressions']['L_command_max_gap'] = command_gap
    print('L vs contact_compare elastic command max gap: %.9g' % command_gap)
    print('L regression PASS' if command_gap < 1e-6 else 'L regression FAIL: gap >= 1e-6')

    # Phi=0 must remove the dynamic term in both channels.
    raw = P2.replay(new, b[:, 11], seed=4312, sensor=False)
    obs = P2.OLD.Obs(np.where(np.isfinite(raw), raw, np.nan), len(raw))
    degaps = {}
    for name in ('W', 'S', 'H'):
        lp = C.build('L'); dp = C.build(name)
        lp.reset(P2.PLATFORM); dp.reset(P2.PLATFORM)
        if hasattr(dp.estimator, 'update_model'): dp.estimator.update_model(0.)
        vals = []
        for n in range(24, len(raw)):
            x = lp.estimator.estimate(obs, n, lp); y = dp.estimator.estimate(obs, n, dp)
            vals.extend(np.asarray(x, float) - np.asarray(y, float))
        gap = float(np.nanmax(np.abs(vals))) if vals else np.nan
        degaps[name] = gap
        print('%s phi=0 -> L max numerical gap: %.9g' % (name, gap))
    out['regressions']['phi_zero_gaps'] = degaps
    print('phi=0 degeneration PASS' if max(degaps.values()) < 1e-6 else 'phi=0 degeneration FAIL')

    # Masking truth columns must not alter a causal candidate's trajectory.
    sc = S2.cell_scene2('breath', 'mid', 'memory', 0, base=4320, T=2.)
    p1 = C.build('TR'); p2 = C.build('TR')
    m1 = P2.episode(sc, C.Measured2(p1), 4321, sensor=True)
    m2 = P2.episode(sc, FullMeasured(p2), 4321, sensor=True)
    mask_gap = float(np.max(np.abs(m1[:, 11] - m2[:, 11])))
    out['regressions']['truth_mask_command_gap'] = mask_gap
    print('NaN truth-column mask command max gap: %.9g' % mask_gap)
    print('truth masking PASS' if mask_gap < 1e-12 else 'truth masking FAIL')

    # One mechanism scene per candidate, with only the declared online state printed.
    for name in C.NAMES:
        mem = S2.cell_scene2('breath', 'mid', 'memory', 0, base=4330, T=3.)
        rate = S2.cell_scene2('breath', 'mid', 'rate', 0, base=4330, T=3.)
        curved = S2.cell_scene2('rest', 'mid', 'curved', 0, base=4330, T=3.)
        lm, pm, tm = C.run_one(name, mem, 4331, sensor=True)
        lr, pr, tr = C.run_one(name, rate, 4332, sensor=True)
        lc, pc, tc = C.run_one(name, curved, 4333, sensor=True)
        mu = getattr(pr, 'mu_hat', np.nan); kd = getattr(pr, 'kd_hat', np.nan); ch = getattr(pc, 'c_hat', np.nan)
        k_for_target = pc.k_good if np.isfinite(getattr(pc, 'k_good', np.nan)) else pc.kbar
        target_c = math.tan(curved['theta0']) + curved['xi_c'] * k_for_target / curved['K_t']
        calls = np.asarray(getattr(pm, 't_call', []) + getattr(pr, 't_call', []) + getattr(pc, 't_call', []), float) * 1e6
        row = dict(phi_memory=float(pm.phi), phi_memory_true=float(mem['phase_true']), phi_rate=float(pr.phi), phi_rate_true=float(rate['phase_true']),
                   mu_hat=float(mu), mu_true=float(rate['mu_k']), kd_hat=float(kd), kd_true=float(rate['k_d']), c_hat=float(ch), c_target=float(target_c),
                   call_us_median=float(np.median(calls)) if len(calls) else np.nan, call_us_p99=float(np.percentile(calls, 99)) if len(calls) else np.nan,
                   seconds=float(tm + tr + tc), rows=int(0 if lm is None else len(lm)))
        out['mechanisms'][name] = row
        print('%-3s phi mem %.5f/%.5f rate %.5f/%.5f | mu %.5f/%.5f kd %.5f/%.5f | c %.5f target %.5f | call %.4fs' %
              (name, row['phi_memory'], row['phi_memory_true'], row['phi_rate'], row['phi_rate_true'], row['mu_hat'], row['mu_true'], row['kd_hat'], row['kd_true'], row['c_hat'], row['c_target'], row['seconds']))
        print('    %s single-call time: median %.0f us, P99 %.0f us' % (name, row['call_us_median'], row['call_us_p99']))
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results' / 'candidate_checks.json').write_text(json.dumps(out, indent=2, allow_nan=True), encoding='utf8')


if __name__ == '__main__':
    main()
