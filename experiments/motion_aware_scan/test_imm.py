"""Unit checks, small causal synthetic grid, and one closed-loop T2 smoke scene.

From experiments/motion_aware_scan (no repository writes or JIT cache writes):
  python -B -X utf8 D:/pip_env/temp/claude/E--github-project-scan-impedance-update/87afbb68-b141-417b-8521-15bc71118679/scratchpad/codex_imm/test_imm.py
The script inserts its own directory on sys.path. --unit-only, --synthetic-only
and --smoke-only run individual parts. Output goes to stdout, never result files.

Synthetic selection is in-sample diagnostic tuning, not a held-out comparison.
All methods consume z[0:j+1] only. R equals the known synthetic noise variance.
The control target is s(t+h)+tau*s_dot(t+h), h=12 ms, tau=19 ms (effective
31 ms ahead, first-order lag compensation). We do NOT extrapolate 31 ms and
then add another 19 ms. Scoring starts at 1 s for a common startup exclusion.
Ramp labels include any ramp intersecting [t,t+h+tau]; steady is its complement.
"""
import argparse
import itertools
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
# The requested launch directory contains proposed.py and common.py.
sys.path.insert(1, str(Path.cwd()))

import numpy as np
from imm import IMMTracker


def load_controller():
    """Keep @njit compilation but suppress disk caches in read-only repo modules."""
    import numba
    original = numba.njit

    def no_disk_cache(*args, **kw):
        kw['cache'] = False
        return original(*args, **kw)

    numba.njit = no_disk_cache
    try:
        import variants
    finally:
        numba.njit = original
    return variants


V = load_controller()
import common

TC, H, TAU = .008, .012, .019


class TrackerTests(unittest.TestCase):
    def test_scalar_kf_reference(self):
        f = IMMTracker(q_smooth=.01, single_model=True)
        f.reset(.001)
        Pp = f.A @ f.P @ f.A.T+f.Q[0]
        xp = f.A @ f.x
        K = Pp[:, 0]/(Pp[0, 0]+f.R)
        expected_x = xp+K*(.002-xp[0])
        expected_P = Pp-np.outer(K, Pp[0])
        np.testing.assert_allclose(f.update(.002), expected_x, rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(f.P, expected_P, rtol=1e-9, atol=1e-16)

    def test_equal_models_reduce_to_kf(self):
        imm = IMMTracker(q_smooth=.003, q_maneuver=.003)
        kf = IMMTracker(q_smooth=.003, single_model=True)
        r = np.random.default_rng(9)
        for z in .002*np.sin(np.arange(300)*.02)+r.normal(0, 1e-5, 300):
            np.testing.assert_allclose(imm.update(z), kf.update(z), atol=1e-13)
            np.testing.assert_allclose(imm.P, kf.P, atol=1e-14)
        np.testing.assert_allclose(imm.mu, [.5, .5], atol=1e-12)

    def test_mixture_moments_and_prediction_readonly(self):
        f = IMMTracker()
        f.reset(0.)
        f.x_models[:] = [[.001, .002, .03], [.004, -.001, -.02]]
        f._mu[:] = [.3, .7]
        f._combine()
        x = f.x.copy()
        expected_P = sum(f.mu[i]*(f.P_models[i]+np.outer(v-x, v-x))
                         for i, v in enumerate(f.x_models))
        np.testing.assert_allclose(f.P, expected_P)
        L = .031
        expected = (x[0]+L*x[1]+.5*L*L*x[2], x[1]+L*x[2])
        np.testing.assert_allclose(f.predict(L), expected)
        np.testing.assert_array_equal(f.x, x)

    def test_switching_covariance_and_floor(self):
        f = IMMTracker(q_smooth=1e-8, q_maneuver=1., R=1e-12, p_stay=.99)
        r = np.random.default_rng(7)
        maneuver_weights = []
        for i in range(600):
            t = i*TC
            z = .002*np.sin(t*1.6)+.006*np.clip((t-2)/.15, 0, 1)
            f.update(z+r.normal(0, 1e-6))
            self.assertAlmostEqual(f.mu.sum(), 1.)
            self.assertTrue(np.all(f.mu >= f.prob_floor))
            self.assertTrue(np.isfinite(f.x).all())
            for P in [*f.P_models, f.P]:
                np.testing.assert_allclose(P, P.T, atol=1e-15)
                self.assertGreaterEqual(np.linalg.eigvalsh(P).min(), -1e-14)
            if 250 <= i <= 275:
                maneuver_weights.append(f.mu[1])
        self.assertGreater(max(maneuver_weights), .5)

    def test_validation_and_reset(self):
        with self.assertRaises(ValueError):
            IMMTracker(R=0)
        f = IMMTracker()
        with self.assertRaises(RuntimeError):
            f.predict(.01)
        f.update(.2)
        f.update(.3)
        f.reset(.4)
        np.testing.assert_allclose(f.x, [.4, 0, 0])
        np.testing.assert_allclose(f.mu, [.5, .5])
        with self.assertRaises(ValueError):
            f.update(np.nan)
        with self.assertRaises(ValueError):
            f.predict(-1)

    @staticmethod
    def policy_log():
        pol = V.MotionAwareIMM()
        pol.reset(common.P.PLATFORMS['PC125'])
        pol.k_good = 1000.
        n = 15
        log = np.zeros((n, 12))
        log[:, 0] = np.arange(n)*TC
        log[:, 2] = np.arange(n)*2e-6
        fin = 4.+np.arange(n)*.01
        decay = (1-pol.afc)**pol.ns
        log[0, 5] = 4.
        for j in range(1, n):
            log[j, 5] = decay*log[j-1, 5]+(1-decay)*fin[j]
        return pol, log, fin

    def test_increment_measurement_dedup_and_k_change(self):
        pol, log, fin = self.policy_log()
        obs = SimpleNamespace(n=2, log=log)
        pol.drift_pred(obs, 2)
        self.assertEqual(pol.Z, 0.)
        first = pol.drift_pred(obs, 3)
        self.assertAlmostEqual(pol.Z, (fin[2]-fin[1])/1000-2e-6)
        Z, x, P = pol.Z, pol.tracker.x.copy(), pol.tracker.P.copy()
        self.assertEqual(first, pol.drift_pred(obs, 3))
        self.assertEqual(pol.Z, Z)
        np.testing.assert_array_equal(pol.tracker.x, x)
        np.testing.assert_array_equal(pol.tracker.P, P)
        pol.k_good = 2000.
        pol.drift_pred(obs, 4)
        self.assertAlmostEqual(pol.Z, Z+(fin[3]-fin[2])/2000-2e-6)
        self.assertAlmostEqual(pol.tracker.R, (.04/2000)**2*4)
        # Constant force/position plus a stiffness change produces no Z jump.
        pol, log, _ = self.policy_log()
        log[:, 5], log[:, 2] = 4., 0.
        pol.drift_pred(SimpleNamespace(log=log), 2)
        pol.k_good = 5000.
        pol.drift_pred(SimpleNamespace(log=log), 3)
        self.assertAlmostEqual(pol.Z, 0., places=15)

    def test_future_rows_do_not_affect_prediction(self):
        a, log, _ = self.policy_log()
        b, _, _ = self.policy_log()
        changed = log.copy()
        changed[7:] = 1e9
        for n in range(1, 8):
            self.assertEqual(a.drift_pred(SimpleNamespace(log=log), n),
                             b.drift_pred(SimpleNamespace(log=changed), n))

    def test_startup_fallback_tracks_without_feedforward_amplification(self):
        pol, log, _ = self.policy_log()
        pol.k_good = np.nan
        pol.q_prev = 100.
        prediction = pol.drift_pred(SimpleNamespace(log=log), 2)
        self.assertIsInstance(prediction, tuple)
        self.assertTrue(pol.tracker.initialized)
        self.assertIsNone(pol.q_prev)
        self.assertAlmostEqual(pol.tracker.R, (.04/pol.kbar)**2*4)
        pol.q_prev = prediction[1]       # what base.step stores after this tick
        pol.drift_pred(SimpleNamespace(log=log), 3)
        self.assertIsNone(pol.q_prev)
        pol.k_good = 1000.
        pol.q_prev = 1e-4
        pol.drift_pred(SimpleNamespace(log=log), 4)
        self.assertEqual(pol.q_prev, 1e-4)  # regular -Delta s is enabled

    def test_contact_loss_recovery_reanchors(self):
        pol, log, _ = self.policy_log()
        for n in range(1, 6):
            pol.step(SimpleNamespace(n=n, log=log))
        self.assertTrue(pol.tracker.initialized)
        log[5, 5] = .1
        pol.step(SimpleNamespace(n=6, log=log))
        self.assertEqual(pol.mode, 1)
        self.assertFalse(pol.tracker.initialized)
        self.assertIsNone(pol.q_prev)
        log[6, 5] = .6
        pol.step(SimpleNamespace(n=7, log=log))
        self.assertEqual(pol.mode, 0)
        self.assertTrue(pol.tracker.initialized)
        self.assertEqual(pol.Z, 0.)
        np.testing.assert_allclose(pol.tracker.x, 0.)
        self.assertEqual(pol.q_prev, 0.)

    def test_build_space_and_perturb(self):
        r = np.random.default_rng(1)
        base = V.tuned('Proposed')
        self.assertEqual(V.build('Proposed', base).params['nw'], base['nw'])
        for name in V.DEFAULTS:
            sampled = V.space(name, r)
            pol = V.build(name, lambda m: base, **sampled)
            pol.reset(common.P.PLATFORMS['PC125'])
            self.assertEqual(pol.tracker.n_models, 1 if name.endswith('KF') else 2)
            self.assertEqual(pol.params['D1'], base['D1'])
            new = V.perturb(name, {**base, **sampled}, r)
            for key in base:
                self.assertEqual(new[key], base[key])
            for key in sampled:
                lo, hi = V.BOUNDS[key]
                self.assertTrue(lo <= new[key] <= hi)
        pol = V.build('Proposed-IMM', V.tuned, q_smooth=2e-6)
        self.assertEqual(pol.motion_params['q_smooth'], 2e-6)


CASES = [(2e-3, .20, 3e-3, .50), (3.5e-3, .265, 5.5e-3, .30),
         (5e-3, .33, 8e-3, .15)]


def signal(t, case):
    """Noiseless breathing and two smoothstep push/retract events; analytic speed."""
    A, freq, amp, ramp = case
    omega = 2*np.pi*freq
    s = A*np.sin(omega*t+.7)
    sd = A*omega*np.cos(omega*t+.7)
    for onset, sign, hold in [(2.2, 1., .8), (6.8, -1., 1.1)]:
        for start, direction in [(onset, sign), (onset+ramp+hold, -sign)]:
            a = np.clip((t-start)/ramp, 0., 1.)
            s = s+direction*amp*(3*a*a-2*a*a*a)
            sd = sd+direction*amp*6*a*(1-a)/ramp
    return s, sd


def ramp_mask(t, ramp):
    out = np.zeros(t.shape, dtype=bool)
    for onset, hold in [(2.2, .8), (6.8, 1.1)]:
        for start in [onset, onset+ramp+hold]:
            out |= (t+H+TAU >= start) & (t <= start+ramp)
    return out


def cubic_predict(z, N=24):
    """Last N points, causal cubic evaluated at h plus tau times its derivative."""
    # Scaling time by window width makes the Vandermonde well-conditioned.
    tw = (N-1)*TC
    u = np.arange(1-N, 1)/(N-1)
    X = np.vander(u, 4, increasing=True)
    h = H/tw
    evaluation = np.array([h**d for d in range(4)])
    evaluation[1:] += TAU/tw*np.array([1., 2*h, 3*h*h])
    weights = evaluation @ np.linalg.pinv(X)
    out = np.full(len(z), np.nan)
    for j in range(N-1, len(z)):
        out[j] = weights @ z[j-N+1:j+1]
    return out


def tracker_predict(z, sigma, params, kf=False):
    f = IMMTracker(Tc=TC, R=sigma**2, single_model=kf, **params)
    out = np.empty(len(z))
    for j, value in enumerate(z):
        f.update(value)
        pos, vel = f.predict(H)
        out[j] = pos+TAU*vel
    return out


def synthetic_grid():
    """8 KF q choices and 18 IMM combinations per noise, pooled across 3 motions."""
    t = np.arange(0., 12., TC)
    active = t >= 1.
    # Identical z for every candidate. No physical signal is polluted by measurement noise.
    r = np.random.default_rng(20261030)
    kf_grid = [dict(q_smooth=float(q)) for q in np.logspace(-7, 0, 8)]
    imm_grid = [dict(q_smooth=qs, q_maneuver=qm, p_stay=p)
                for qs, qm, p in itertools.product((1e-7, 1e-5, 1e-3), (.01, .1, 1.), (.95, .99))]
    print('\nSynthetic diagnostic: h=12 ms, tau=19 ms, N=24; score t>=1 s.', flush=True)
    print('Best selected by pooled full RMSE on these same signals (no held-out claim).', flush=True)
    print('| noise [um] | method | full [um] | ramps [um] | steady [um] | best parameters |', flush=True)
    print('|---:|---|---:|---:|---:|---|', flush=True)
    results = []
    for sigma_um in (1., 10., 50.):
        sigma = sigma_um*1e-6
        data = []
        for case in CASES:
            truth, _ = signal(t, case)
            future, speed = signal(t+H, case)
            target = future+TAU*speed
            z = truth+r.normal(0., sigma, len(t))
            data.append((z, target, ramp_mask(t, case[3])))

        def evaluate(method, params):
            sq = [[], [], []]
            for z, target, ramps in data:
                pred = (cubic_predict(z) if method == 'Cubic' else
                        tracker_predict(z, sigma, params, kf=(method == 'KF')))
                err2 = (pred-target)**2
                for arr, mask in zip(sq, (active, active & ramps, active & ~ramps)):
                    arr.extend(err2[mask])
            return tuple(float(np.sqrt(np.mean(a))*1e6) for a in sq)

        for method, grid in [('Cubic', [{}]), ('KF', kf_grid), ('IMM', imm_grid)]:
            scores = [(evaluate(method, p), p) for p in grid]
            rmse, params = min(scores, key=lambda row: row[0][0])
            results.append(dict(sigma_um=sigma_um, method=method, rmse_um=rmse, params=params))
            print(f'| {sigma_um:g} | {method} | {rmse[0]:.3f} | {rmse[1]:.3f} | '
                  f'{rmse[2]:.3f} | {params} |', flush=True)
    return results


def smoke():
    import scenarios
    p, T = scenarios.make('T2', np.random.default_rng(11), dict(stiff=True, ramp=.3))
    pol = V.build('Proposed-IMM', V.tuned)
    log, _ = common.run(pol, p, T, seed=11)
    assert log is not None, 'plant diverged'
    assert len(log) == round(T/TC), 'incomplete trajectory'
    assert np.isfinite(log).all(), 'nonfinite closed-loop values'
    assert np.isfinite(pol.tracker.x).all(), 'nonfinite tracker state'
    m = common.task_metrics(log, f_exc=pol.params['f_exc'])
    assert not m['diverged']
    print('\nT2 stiff smoke (scenario seed=11, noise seed=11, ramp=0.3 s, '
          'platform=PC125, warm=1 s):', flush=True)
    print(f"duration={T:g}s, ticks={len(log)}, force RMSE={m['rmse']:.6f} N, "
          f"peak={m['peak']:.6f} N, loss_fraction={m['loss']:.6f}, "
          f"loss_episodes={m['n_loss']}, diverged={m['diverged']}", flush=True)
    print(f'tracker defaults={pol.motion_params}, final_mu={pol.tracker.mu}', flush=True)
    return m


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--unit-only', action='store_true')
    mode.add_argument('--synthetic-only', action='store_true')
    mode.add_argument('--smoke-only', action='store_true')
    args = parser.parse_args()
    if not (args.synthetic_only or args.smoke_only):
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(TrackerTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        if not result.wasSuccessful():
            raise SystemExit(1)
    if not (args.unit_only or args.smoke_only):
        synthetic_grid()
    if not (args.unit_only or args.synthetic_only):
        smoke()


if __name__ == '__main__':
    main()
