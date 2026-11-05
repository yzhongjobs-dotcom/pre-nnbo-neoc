"""Fixed-distribution scan scenes for the contact-model selection matrix."""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'cobot_scan_benchmark'))
import plant as OLD
import scenes as S
import sweeps as SW

CELLS2 = SW.CELLS + [(c, 'relaxed') for c in ('rest', 'breath', 'shift_slow', 'flinch')]
REGIME_CODE = {'nominal': 0, 'stop_go': 1, 'reversal': 2, 'curved': 3, 'rate': 4, 'memory': 5, 'power': 6, 'rate_sweep': 7}


POWER_TAU = np.array([.003, .03, .3])
POWER_BETA = np.array([.2, .25, .3])
POWER_ALPHA = np.array([[1.463, .795, .277], [2.342, 1.1, .324], [3.581, 1.447, .359]])   # least-squares fit of each slope over .3-100 Hz


def _profile(regime, T, rng):
    n = int(round(T * 1000.)) + 2
    t = np.arange(n) * 1e-3
    if regime in ('nominal', 'curved', 'rate', 'memory', 'power', 'rate_sweep'):
        v = float(rng.uniform(.005, .03))
        return OLD.X0 + v * t, v
    if regime == 'stop_go':
        vmax = float(rng.uniform(.01, .06))
    else:
        vmax = float(rng.uniform(.01, .04))
    accel = float(rng.uniform(.1, 1.))
    target = np.full(n, vmax)
    if regime == 'stop_go':
        for ts in np.linspace(1.5, max(1.6, T-1.5), int(rng.integers(2, 4))):
            hold = float(rng.uniform(.2, .8))
            target[(t >= ts) & (t < ts+hold)] = 0.
    else:
        sign = 1.
        for ts in np.linspace(2., max(2.1, T-2.), int(rng.integers(1, 3))):
            sign = -sign
            target[t >= ts] = sign*vmax
    # the plan follows its target under the acceleration limit, so no reversal or stop has a velocity jump
    vel = np.zeros(n); step = accel*1e-3
    for i in range(1, n):
        vel[i] = vel[i-1]+min(max(target[i]-vel[i-1], -step), step)
    x = np.r_[0., np.cumsum(vel[:-1])*1e-3]
    x = OLD.X0+.002+x-x.min()      # stays on the field, which starts at X0
    return x, float(np.max(np.abs(vel)))


def make2(cond, tissue, regime, rng, T=10., overrides=None, xcmd=None):
    """Call the original disturbance generator and attach the fixed contact law.

    ``rng`` is expected to have been created as
    ``default_rng([base, code, i, regime_code])`` by the caller.  No candidate
    changes these sampled parameters after scene construction.
    """
    if regime not in REGIME_CODE:
        raise ValueError(regime)
    o = dict(overrides or {})
    mu_k = float(rng.uniform(.1, .45))
    mu_s = mu_k * (1. + float(rng.uniform(0., .5)))
    v_s = float(rng.uniform(.002, .01))
    v_e = .0005
    k_d = float(rng.uniform(0., 13.))
    nu = .45
    g = 10. ** float(rng.uniform(0., 1.))
    k1, k2 = .001 * g, .003 * g
    theta0 = float(rng.normal(0., math.radians(2.9)))
    c_slope = float(rng.uniform(.2, .6))
    K_t = 10. ** float(rng.uniform(3.7, 4.7))
    m_t = float(rng.uniform(1., 3.))
    zeta_t = float(rng.uniform(.3, .7))
    xi_c = float(rng.uniform(-.3, .3))
    xtraj, vscale = _profile(regime, T, rng)
    if xcmd is not None:
        xtraj = np.asarray(xcmd, float)
        vscale = float(np.max(np.abs(np.diff(xtraj))) * 1000.)
    if regime == 'curved':
        theta0 = float(rng.choice([-1., 1.]) * rng.uniform(math.radians(5.), math.radians(15.)))
        c_slope = 1.
    alpha = np.array([.35, .35, .35], float)
    lam = float(rng.uniform(.1, .4))
    if regime == 'rate':
        alpha = np.array([.05, .05, .05], float)
        lam = float(rng.uniform(1., 3.))
    elif regime == 'memory':
        alpha = np.array([.6, .6, .6], float)
        lam = .1
    # The original generator allocates the spatial fields using its speed
    # argument.  vscale is the conservative upper bound for all profiles.
    sc = S.make(cond, tissue, rng, T=T, v=max(vscale, .005), overrides=o)
    sc['xcmd'] = np.asarray(xtraj, float)
    sc['v'] = vscale
    sc['mu_k'] = mu_k; sc['mu_s'] = mu_s; sc['v_s'] = v_s; sc['v_e'] = v_e
    sc['k_d'] = k_d; sc['nu'] = nu; sc['k1'] = k1; sc['k2'] = k2
    sc['theta0'] = theta0; sc['c_slope'] = c_slope; sc['K_t'] = K_t
    sc['m_t'] = m_t; sc['zeta_t'] = zeta_t; sc['xi_c'] = xi_c
    sc['regime'] = regime; sc['regime_code'] = REGIME_CODE[regime]
    if regime in ('rate', 'memory'):
        sc['alpha'] = alpha; sc['lam'] = lam
    if regime == 'power':
        # the phase of the rate regime, from the weak power law's own three branches instead of a rate term
        beta = float(rng.uniform(.2, .3))
        sc['alpha'] = np.array([np.interp(beta, POWER_BETA, POWER_ALPHA[:, q]) for q in range(3)]); sc['tau'] = POWER_TAU.copy(); sc['lam'] = .1
    return sc


def cell_scene2(cond, tissue, regime, i, base=420000, T=10., overrides=None):
    code = sum(ord(c) * (j + 1) for j, c in enumerate(cond + tissue)) % 9973
    rng = np.random.default_rng([base, code, i, REGIME_CODE[regime]])
    return make2(cond, tissue, regime, rng, T=T, overrides=overrides)
