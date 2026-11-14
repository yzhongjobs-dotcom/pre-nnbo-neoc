"""Stiffness as a field over scan position, seen through the footprint: a causal estimator with no assumption on levels.

Representation. ln k(t) = (1/L) * integral of e over the footprint [x(t)-L/2, x(t)+L/2], e the tissue density along the scan path.
The field is static; time enters only through the measured scan position x(t). Prior: e is continuous along x (Brownian, variance
qe per metre). The filter is a Kalman filter on the cells of e under the footprint; cells enter at the leading edge with the value of
their neighbour and leave at the trailing edge. Memory length, slope persistence and the fall that follows a rise all come from L and
x(t), which are known, not from tuning per scene.

Test here (representation only): truth = recorded ln k_eff of the 42 stored closed-loop runs (all scene groups), observation = truth +
white noise with the information of the 24-tick window (sd .17 per window). Comparators: the window mean, a random-walk memory in
time, a random-walk memory in position, a local-linear-trend filter in position.
Run: python -X utf8 footprint_field_filter.py
"""
import os, sys
for _k in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'): os.environ[_k] = '1'
sys.dont_write_bytecode = True; sys.path.insert(0, os.getcwd())
from multiprocessing import Pool
import numpy as np

TC = .008; R = 24*.17**2; NSEED = 3
QE = (.3, 1., 3., 10., 30., 100., 300., 1000.); CL = (.35, .5, .75, 1.)
QT = (3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2); QX = (.03, .1, .3, 1., 3., 10., 30., 100.)
QS = (1e2, 1e3, 1e4, 1e5, 1e6, 1e7)


def window(obs, good, nw=24):
    o = np.where(good, obs, 0.); c = np.cumsum(np.r_[0., o]); n = np.cumsum(np.r_[0., good.astype(float)])
    out = np.full(len(obs), np.nan); cnt = n[nw:]-n[:-nw]; out[nw-1:] = (c[nw:]-c[:-nw])/np.maximum(cnt, 1.); return out


def rw(obs, good, q):
    """Random-walk Kalman filter; q is the per-tick process variance (array or scalar)."""
    q = np.broadcast_to(q, obs.shape); j0 = int(np.argmax(good)); x, P = obs[j0], R; out = np.full(len(obs), np.nan)
    for j in range(j0, len(obs)):
        P += q[j]
        if good[j]:
            K = P/(P+R); x += K*(obs[j]-x); P *= 1-K
        out[j] = x
    return out


def llt(obs, good, dx, qs):
    """Local linear trend in position: S += S' dx, S' random walk with variance qs per metre."""
    j0 = int(np.argmax(good)); s = np.array([obs[j0], 0.]); P = np.diag([R, 1e2]); out = np.full(len(obs), np.nan)
    for j in range(j0, len(obs)):
        F = np.array([[1., dx[j]], [0., 1.]]); s = F@s; P = F@P@F.T; P[1, 1] += qs*abs(dx[j])
        if good[j]:
            S = P[0, 0]+R; K = P[:, 0]/S; s = s+K*(obs[j]-s[0]); P = P-np.outer(K, P[0])
        out[j] = s[0]
    return out


def field(obs, good, x, L, qe, ncell=24):
    """Kalman filter on the density cells under the footprint. Returns the causal estimate of the footprint average."""
    dc = L/ncell; j0 = int(np.argmax(good)); out = np.full(len(obs), np.nan)
    lo = int(np.floor((x[j0]-.5*L)/dc)); hi = int(np.floor((x[j0]+.5*L)/dc)); m = hi-lo+1
    idx = np.arange(m); e = np.full(m, obs[j0]); P = 1.+qe*dc*np.minimum.outer(idx, idx)       # Brownian prior, vague level
    for j in range(j0, len(obs)):
        a, b = x[j]-.5*L, x[j]+.5*L
        while (hi+1)*dc < b:                      # a cell enters: its neighbour's value plus a Brownian increment
            e = np.r_[e, e[-1]]; col = P[:, -1]; P = np.block([[P, col[:, None]], [col[None, :], np.array([[P[-1, -1]+qe*dc]])]]); hi += 1
        while (lo+1)*dc <= a and hi > lo:         # a cell has left the footprint
            e = e[1:]; P = P[1:, 1:]; lo += 1
        edges = (lo+np.arange(hi-lo+2))*dc; H = np.clip(np.minimum(edges[1:], b)-np.maximum(edges[:-1], a), 0., None)/L
        if good[j]:
            PH = P@H; S = H@PH+R; K = PH/S; e = e+K*(obs[j]-H@e); P = P-np.outer(K, PH)
        out[j] = H@e
    return out


def load(i):
    import corrected_comparison as CC, overlay_compare as O, memory_segment_eval as SE
    O.setup(); c = CC.cases()[i]; sc = SE.MC.scene_for(tuple(c['cell']), CC.BASE, c['variant']); g = sc['geometry']
    log = np.load(f'results/huber_windows/case_{i:02d}.npz')['log']
    return c['group'], (g['length'] if sc['slide'] == 'long' else g['width']), log[:, 0], log[:, 1], log[:, 3], log[:, 9]


def one(i):
    group, L, t, x, fn, ke = load(i)
    good = (ke > 0) & (fn > .5); lk = np.log(np.where(ke > 0, ke, 1.)); dx = np.r_[0., np.diff(x)]
    ev = good & (t > 1.5); rate = np.gradient(lk, TC); rise = ev & (rate > .5); fall = ev & (rate < -.5)
    res = {}
    for seed in range(NSEED):
        obs = lk+np.random.default_rng(1000*i+seed).normal(size=len(lk))*np.sqrt(R)
        est = {('window',): window(obs, good)}
        for q in QT: est[('rw_t', q)] = rw(obs, good, q)
        for q in QX: est[('rw_x', q)] = rw(obs, good, q*np.abs(dx))
        for q in QS: est[('llt_x', q)] = llt(obs, good, dx, q)
        for q in QE:
            for c in CL: est[('field', q, c)] = field(obs, good, x, c*L, q)
        for key, v in est.items():
            err = v-lk
            res.setdefault(key, []).append((np.nanmean(err[ev]**2), np.nanmean(err[rise]) if rise.sum() > 20 else np.nan,
                                            np.nanmean(err[fall]) if fall.sum() > 20 else np.nan, np.nanpercentile(err[ev], 1)))
    return i, group, {k: np.nanmean(np.array(v), 0) for k, v in res.items()}, float(np.std(np.diff(lk[ev][::24])))


if __name__ == '__main__':
    with Pool(8) as p:
        rows = p.map(one, range(42))
    groups = ('A', 'A_slow_thin', 'B', 'A+B', 'A_large', 'combo', 'cough'); main = groups[:5]
    keys = list(rows[0][2]); fam = lambda name: [k for k in keys if k[0] == name]
    mse = lambda key, sel: float(np.mean([r[2][key][0] for r in rows if sel(r)]))
    dev = lambda r: r[0] % 2 == 0 and r[1] in main; test = lambda r: r[0] % 2 == 1 and r[1] in main
    # one global setting per method: minimise, on the even-numbered cases, the mean over groups of MSE relative to the window
    def score(key, sel):
        return float(np.mean([mse(key, lambda r, g=g: sel(r) and r[1] == g)/mse(('window',), lambda r, g=g: sel(r) and r[1] == g) for g in main]))
    pick = {name: min(fam(name), key=lambda k: score(k, dev)) for name in ('rw_t', 'rw_x', 'llt_x', 'field')}
    print('one global setting per method, chosen on even cases:', pick)
    print('score = mean over the five groups of MSE / window MSE; dev (even cases), test (odd cases):')
    for name, k in pick.items():
        print(f'  {name:6s} dev {score(k, dev):.3f} test {score(k, test):.3f}')
    print('\nrms error of ln k, all cases of the group | window | rw_t global, (best for the group) | rw_x global | llt_x global | field global, (best for the group) || mean error on rises / falls: window, rw_t, field || 1st percentile of error: window, rw_t, field || true change per window, sd')
    for g in groups:
        sel = lambda r, g=g: r[1] == g; rr = [r for r in rows if sel(r)]
        rm = lambda k: np.sqrt(mse(k, sel)); col = lambda k, c: float(np.nanmean([r[2][k][c] for r in rr]))
        bt = min(fam('rw_t'), key=lambda k: mse(k, sel)); bf = min(fam('field'), key=lambda k: mse(k, sel))
        w, a, f = ('window',), pick['rw_t'], pick['field']
        print(f"{g:12s} | {rm(w):.3f} | {rm(a):.3f} ({rm(bt):.3f} q {bt[1]:g}) | {rm(pick['rw_x']):.3f} | {rm(pick['llt_x']):.3f} | {rm(f):.3f} ({rm(bf):.3f} qe {bf[1]:g} c {bf[2]:g}) || "
              f"{col(w, 1):+.3f}/{col(w, 2):+.3f}, {col(a, 1):+.3f}/{col(a, 2):+.3f}, {col(f, 1):+.3f}/{col(f, 2):+.3f} || {col(w, 3):+.2f} {col(a, 3):+.2f} {col(f, 3):+.2f} || {np.mean([r[3] for r in rr]):.3f}")
    print('\nfield filter, rms by group over the grid (rows qe, columns footprint factor c):')
    for g in main:
        sel = lambda r, g=g: r[1] == g
        print(g); [print(f'  qe {q:6g}: '+' '.join(f'{np.sqrt(mse(("field", q, c), sel)):.3f}' for c in CL)) for q in QE]
