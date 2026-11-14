"""Synthetic check of a footprint-structured prior for the stiffness trajectory.

The probe-frequency stiffness is a footprint average of a spatial stiffness density along the scan: k(t) = mean of the density over the
last M ticks of entering tissue (M = footprint length / scan speed / Tc). What enters must leave M ticks later, so a rise is mirrored by a
fall: increments at lag M are negatively correlated. Compare three causal estimators of k now from per-tick noisy observations
(variance R per tick, the information of the 24-tick window spread over its ticks):
  window  : mean of the last 24 observations (the present estimator's noise and lag),
  rw      : random-walk Kalman filter, process variance chosen best for this data (a tuned memory),
  moving  : Kalman filter on the shift-register state of the entering density (the footprint model).
Density model: rib (high) and gap (low) segments of random length, so the truth is not the filter's Gaussian prior.
Run: python -X utf8 footprint_prior_synthetic.py
"""
import numpy as np

TC = .008


def make(n, M, rng, rib=(40, 90), gap=(20, 70), contrast=1.6, jitter=.15):
    """Entering log-stiffness density: alternating gap/rib segments (lengths in ticks), levels 0 and `contrast` with jitter."""
    e = []
    while len(e) < n+M:
        g = rng.integers(*gap); r = rng.integers(*rib)
        e += [rng.normal(0., jitter)]*g+[contrast+rng.normal(0., jitter)]*r
    e = np.array(e[:n+M])
    c = np.cumsum(np.r_[0., e]); k = (c[M:]-c[:-M])/M            # footprint average
    return e, k[:n]


def window(obs, nw=24):
    c = np.cumsum(np.r_[0., obs]); out = np.full(len(obs), np.nan); out[nw-1:] = (c[nw:]-c[:-nw])/nw
    return out


def rw(obs, R, q):
    x, P = obs[0], R; out = np.empty(len(obs))
    for j, z in enumerate(obs):
        P += q; K = P/(P+R); x += K*(z-x); P *= 1-K; out[j] = x
    return out


def moving(obs, R, M, mean, var_cell, cells=None):
    """Shift-register Kalman filter with `cells` states spanning the footprint (each cell = M/cells ticks)."""
    m = cells or M; step = M//m
    x = np.full(m, mean); P = np.eye(m)*var_cell; H = np.ones(m)/m; out = np.empty(len(obs)); frac = 0
    for j, z in enumerate(obs):
        frac += 1
        if frac == step:                      # the oldest cell leaves, a new one enters with the prior
            frac = 0
            x = np.r_[x[1:], mean]; P = np.block([[P[1:, 1:], np.zeros((m-1, 1))], [np.zeros((1, m-1)), np.array([[var_cell]])]])
        PH = P@H; S = H@PH+R; K = PH/S
        x = x+K*(z-H@x); P = P-np.outer(K, PH)
        out[j] = H@x
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    R = 24*.17**2
    print('M ticks (footprint transit) | rms error of ln k: window / best random-walk memory / footprint filter | on rises, on falls (|rate| > 1 /s): same three | sd of ln k')
    for M in (60, 75, 120):
        n = 40000
        e, k = make(n, M, rng)
        obs = k+rng.normal(size=n)*np.sqrt(R)
        dk = np.gradient(k, TC); rise = dk > 1.; fall = dk < -1.
        w = window(obs)
        best = min(((np.nanmean((rw(obs, R, q)[500:]-k[500:])**2), q) for q in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3)))
        r = rw(obs, R, best[1])
        mv = moving(obs, R, M, float(np.mean(e)), float(np.var(e)), cells=M//3)
        s = slice(500, None)
        f = lambda est, msk=None: float(np.sqrt(np.nanmean(((est-k)[s][msk[s]] if msk is not None else (est-k)[s])**2)))
        print(f'{M:4d} ({M*TC:.2f} s) | {f(w):.3f} / {f(r):.3f} (q {best[1]:g}) / {f(mv):.3f} | rise {f(w, rise):.3f} {f(r, rise):.3f} {f(mv, rise):.3f} | fall {f(w, fall):.3f} {f(r, fall):.3f} {f(mv, fall):.3f} | {k.std():.2f}')
        # bias on falls: positive means over-estimate (lag)
        print(f'      mean error on rises / falls: window {np.nanmean((w-k)[s][rise[s]]):+.3f} {np.nanmean((w-k)[s][fall[s]]):+.3f} | rw {np.mean((r-k)[s][rise[s]]):+.3f} {np.mean((r-k)[s][fall[s]]):+.3f} | footprint {np.mean((mv-k)[s][rise[s]]):+.3f} {np.mean((mv-k)[s][fall[s]]):+.3f}')
