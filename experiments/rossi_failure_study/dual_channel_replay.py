"""Two-channel causal replay; no plant imports, simulation, or feedback.

Run from experiments/rossi_failure_study: python -X utf8 dual_channel_replay.py
Dependencies: numpy, numba, matplotlib. Recordings are read-only.
--project-root and --output-root permit running when the project is read-only.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import itertools
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import time
import numpy as np
from numba import njit

TC, LAG, RY = .008, 9, (5e-6)**2
BASE_GROUPS = ('A', 'A_slow_thin', 'A_large', 'B', 'A+B')
SEGMENTS = ('all', 'plateau', 'slow', 'rise', 'fall', 'reversal', 'event')
GRID1 = list(itertools.product((0., 1e-7, 1e-6, 1e-5),
                              (1e-4, 1e-3, 3e-3, 1e-2), (.95, .99)))
GRID2 = [a+b for a in GRID1 for b in itertools.product(
    (1e-6, 1e-5, 1e-4, 1e-3, 1e-2), (1e-8, 1e-6))]


def trailing(a, n):
    """Available-prefix mean at startup; strictly trailing thereafter."""
    sums = np.r_[0., np.cumsum(np.asarray(a, dtype=float))]
    j = np.arange(len(a))
    lo = np.maximum(0, j-n+1)
    return (sums[j+1]-sums[lo])/(j+1-lo)


def centered(a, n):
    """Centred mean with truncated, renormalized edges (scoring only)."""
    sums = np.r_[0., np.cumsum(a)]
    j = np.arange(len(a))
    lo, hi = np.maximum(0, j-n//2), np.minimum(len(a), j+n//2+1)
    return (sums[hi]-sums[lo])/(hi-lo)


def preprocess(t, y, force, raw, Y):
    """Truth-free interface: these five copied arrays are all the inputs."""
    D = np.diff(np.eye(24), n=4, axis=0)
    rho = (2*np.sin(np.pi*14*TC))**-8
    M = np.eye(24)-np.linalg.solve(np.eye(24)+rho*D.T@D, np.eye(24))
    yw = Y[:, -24:]
    info = np.einsum('ij,jk,ik->i', yw, M, yw)
    finite = np.isfinite(raw)
    m = np.full(len(t), np.nan)
    m[finite] = np.log(np.clip(raw[finite], 200., 3e5))
    rel = np.full(len(t), np.inf)
    ok = finite & np.all(np.isfinite(yw), axis=1) & (info > 0) & (raw > 0)
    rel[ok] = .0245/(raw[ok]*np.sqrt(info[ok]))
    valid = ok & (rel <= 1.)
    return (np.array(t, copy=True), trailing(y, 6), trailing(force, 6),
            m, 24*rel**2, valid, finite)


@njit(nogil=True)
def weights(prior, ll):
    lw = np.log(prior)+ll
    w = np.exp(lw-np.max(lw))
    return 1e-12+(1.-2e-12)*w/np.sum(w)


@njit(nogil=True)
def probe_filter(data, pars):
    t, ys, fs, m, R, valid, finite = data
    q0, q1, stay = pars
    qs = np.array([q0, q1])
    trans = np.array([[stay, 1-stay], [1-stay, stay]])
    out = np.full(len(t), np.nan)
    xs, ps, mu = np.zeros(2), np.zeros(2), np.ones(2)*.5
    active, last_raw = False, -1000
    for j in range(len(t)):
        if finite[j]:
            last_raw = j
        if fs[j] < 1. or j-last_raw > 25 or not np.isfinite(fs[j]):
            active = False
            if j > 0:
                out[j] = out[j-1]
            continue
        if not active:
            if not valid[j]:
                if j > 0:
                    out[j] = out[j-1]
                continue
            xs[:] = m[j]
            ps[:] = R[j]
            mu[:] = .5
            active = True
        else:
            prior = mu@trans
            nx, npv, ll = np.zeros(2), np.zeros(2), np.zeros(2)
            for h in range(2):
                mix = mu*trans[:, h]/prior[h]
                xp = np.sum(mix*xs)
                pp = np.sum(mix*(ps+(xs-xp)**2))+qs[h]
                if valid[j]:
                    v = m[j]-xp
                    S = pp+R[j]
                    K = pp/S
                    xp += K*v
                    pp = (1-K)**2*pp+K*K*R[j]
                    ll[h] = -.5*(math.log(2*math.pi*S)+v*v/S)
                nx[h], npv[h] = xp, pp
            xs, ps, mu = nx, npv, weights(prior, ll)
        out[j] = np.sum(mu*xs)
    return out


@njit(nogil=True)
def scalar_update(x, P, H, residual, variance):
    PH = P@H
    S = H@PH+variance
    if not np.isfinite(S) or S <= 0:
        raise ValueError('nonpositive innovation variance')
    K = PH/S
    J = np.eye(5)-np.outer(K, H)
    Pnew = J@P@J.T+variance*np.outer(K, K)
    return x+K*residual, (Pnew+Pnew.T)*.5, -.5*(math.log(2*math.pi*S)+residual**2/S)


@njit(nogil=True)
def probe_model(x, dy, flag):
    """Analytic h and dh/dx for the requested lagged indentation identity."""
    tau = LAG*TC
    p = math.exp(x[3])
    plag = p-dy+x[1]*tau-x[2]*tau*tau/2
    H = np.zeros(5)
    if plag <= .05*p:
        return 0., H, False
    H[1], H[2], H[3], H[4] = -tau/plag, tau*tau/(2*plag), -p/plag, 1.
    return x[4]+math.log(flag)-math.log(plag), H, True


@njit(nogil=True)
def fused_filter(data, pars):
    t, ys, fs, m, R, valid, finite = data
    q0, q1, stay, qb, qc = pars
    trans = np.array([[stay, 1-stay], [1-stay, stay]])
    A = np.eye(5)
    A[0, 1], A[0, 2], A[1, 2] = TC, TC*TC/2, TC
    dt = TC
    jerk = np.array([[dt**5/20, dt**4/8, dt**3/6],
                     [dt**4/8, dt**3/3, dt**2/2],
                     [dt**3/6, dt**2/2, dt]])
    Q = np.zeros((2, 5, 5))
    for h in range(2):
        Q[h, :3, :3] = qb*jerk
        Q[h, 3, 3], Q[h, 4, 4] = (q0 if h == 0 else q1), qc
    out, sd, cout = np.full(len(t), np.nan), np.full(len(t), np.nan), np.full(len(t), np.nan)
    xs, Ps, mu = np.zeros((2, 5)), np.zeros((2, 5, 5)), np.ones(2)*.5
    active, last_raw, skips, resets = False, -1000, 0, 0
    for j in range(len(t)):
        if finite[j]:
            last_raw = j
        if fs[j] < 1. or j-last_raw > 25 or not np.isfinite(fs[j]+ys[j]):
            active = False
            if j > 0:
                out[j], sd[j], cout[j] = out[j-1], sd[j-1], cout[j-1]
            continue
        if not active:
            if not valid[j]:
                if j > 0:
                    out[j], sd[j], cout[j] = out[j-1], sd[j-1], cout[j-1]
                continue
            c = math.log(3.)
            ell = c+math.log(fs[j])-m[j]
            p = math.exp(ell)
            x = np.array([ys[j]-p, 0., 0., ell, c])
            # Independent errors dy, dm, dc, dbd, dbdd; differential definitions:
            # dell=dc-dm, db=dy-p*dell. This preserves all shared c errors.
            J = np.zeros((5, 5))
            J[0, 0], J[0, 1], J[0, 2] = 1., p, -p
            J[1, 3], J[2, 4] = 1., 1.
            J[3, 1], J[3, 2], J[4, 2] = -1., 1., 1.
            V = np.diag(np.array([RY, R[j], math.log(3.)**2, .03**2, .3**2]))
            P = J@V@J.T
            for h in range(2):
                xs[h], Ps[h] = x, P
            mu[:] = .5
            active = True
            resets += 1
        else:
            prior = mu@trans
            nx, nP, ll = np.empty_like(xs), np.empty_like(Ps), np.zeros(2)
            for h in range(2):
                mix = mu*trans[:, h]/prior[h]
                xm = mix@xs
                Pm = np.zeros((5, 5))
                for i in range(2):
                    d = xs[i]-xm
                    Pm += mix[i]*(Ps[i]+np.outer(d, d))
                x = A@xm
                P = A@Pm@A.T+Q[h]
                p = math.exp(x[3])
                H = np.array([1., 0., 0., p, 0.])
                x, P, ly = scalar_update(x, P, H, ys[j]-x[0]-p, RY)
                ll[h] = ly
                if valid[j] and j >= LAG and fs[j-LAG] >= 1.:
                    pred, H, usable = probe_model(x, ys[j]-ys[j-LAG], fs[j-LAG])
                    if usable:
                        x, P, lm = scalar_update(x, P, H, m[j]-pred, R[j])
                        ll[h] += lm
                    else:
                        skips += 1
                nx[h], nP[h] = x, P
            xs, Ps, mu = nx, nP, weights(prior, ll)
        xbar = mu@xs
        Pbar = np.zeros((5, 5))
        for h in range(2):
            d = xs[h]-xbar
            Pbar += mu[h]*(Ps[h]+np.outer(d, d))
        out[j] = xbar[4]+math.log(fs[j])-xbar[3]
        sd[j] = math.sqrt(max(0., Pbar[4, 4]+Pbar[3, 3]-2*Pbar[3, 4]))
        cout[j] = xbar[4]
        if not np.isfinite(out[j]+sd[j]+cout[j]):
            raise ValueError('nonfinite fused filter state')
    return out, sd, cout, skips, resets


def synthetic_check():
    """Prespecified signal-level test, before opening any recordings."""
    t = np.arange(0., 12., TC)
    phase = np.mod(t, 2.4)
    lam = 7.6+.7*np.clip((phase-.8)/.4, 0, 1)-.7*np.clip((phase-2.)/.4, 0, 1)
    ramp = ((phase >= .8) & (phase < 1.2)) | (phase >= 2.)
    b = .004*np.sin(2*np.pi*.25*t)
    y = b+4*8/np.exp(lam)
    rng = np.random.default_rng(20261112)
    R = np.full(len(t), 24*.17**2)
    m = lam[np.maximum(0, np.arange(len(t))-LAG)]+rng.normal(size=len(t))*np.sqrt(R)
    valid = np.ones(len(t), dtype=np.bool_)
    # y is the trajectory observation for this test: no extra smoothing delay
    # added to the explicitly nine-tick-delayed synthetic probe observations.
    data = (t, y, np.full(len(t), 8.), m, R, valid, valid)
    p1, p2 = (1e-6, 3e-3, .99), (1e-6, 3e-3, .99, 1e-5, 1e-8)
    a1 = probe_filter(data, p1)
    a2, sd, c, skips, resets = fused_filter(data, p2)
    result = {'seed': 20261112, 'duration_s': 12., 'score_after_s': 1.5,
              'parameters_A1': p1, 'parameters_A2': p2,
              'observation': 'y directly, probe delayed 9 ticks; no extra smoothing',
              'learned_g': float(np.exp(c[-1])), 'true_g': 4., 'skipped_mode_probe_updates': skips}
    for name, a in [('A1', a1), ('A2', a2)]:
        result[name] = {s: float(np.sqrt(np.mean((a[mask & (t > 1.5)]-lam[mask & (t > 1.5)])**2)))
                        for s, mask in [('ramp', ramp), ('flat', ~ramp)]}
    x = np.array([.001, .003, -.02, math.log(.01), math.log(4.)])
    _, H, _ = probe_model(x, .0003, 8.)
    numerical = np.zeros(5)
    for i in range(5):
        xp, xm = x.copy(), x.copy()
        xp[i] += 1e-7
        xm[i] -= 1e-7
        numerical[i] = (probe_model(xp, .0003, 8.)[0]-probe_model(xm, .0003, 8.)[0])/2e-7
    result['jacobian_max_absolute_error'] = float(np.max(np.abs(H-numerical)))
    assert np.allclose(H, numerical, rtol=1e-6, atol=1e-7)
    assert np.all(np.isfinite(a2)) and np.all(sd >= 0)
    print('Synthetic check (before recordings): '+json.dumps(result), flush=True)
    return result


def raw_arm(raw):
    held, k = np.nan, np.full(len(raw), np.nan)
    for j, value in enumerate(raw):
        if np.isfinite(value):
            held = float(np.clip(value, 200, 3e5))
        k[j] = held
    return np.log(k), k


def masks_for(log, meta):
    t, kt, fax = log[:, 0], log[:, 9], log[:, 7]
    lam = np.log(np.maximum(kt, np.finfo(float).tiny))
    dl = trailing(np.gradient(trailing(lam, 6), TC), 5)
    event = np.zeros(len(t), dtype=bool)
    for ev in meta['event_info']:
        start, end = ev['start'], ev['end']
        event |= (t >= start-.05) & (t <= end+.6)
    rev = np.zeros(len(t), dtype=bool)
    signs = np.sign(dl)
    for j in range(1, len(signs)):
        if signs[j] == 0:
            signs[j] = signs[j-1]
    for j in np.flatnonzero(signs[1:]*signs[:-1] < 0)+1:
        left = (t >= t[j]-.15) & (t < t[j]) & (np.abs(dl) > 1)
        right = (t >= t[j]) & (t <= t[j]+.15) & (np.abs(dl) > 1)
        if np.any(left) and np.any(right):
            rev |= np.abs(t-t[j]) <= .1
    score = (t > 1.5) & (kt > 0) & (fax > .5)
    masks = {'all': score, 'plateau': np.abs(dl) < .2,
             'slow': (np.abs(dl) > .4) & (np.abs(dl) <= 1),
             'rise': dl > 1, 'fall': dl < -1, 'reversal': rev, 'event': event}
    for s in masks:
        if s != 'all':
            masks[s] &= score & (event if s == 'event' else ~event)
    masks['dynamic'] = masks['slow'] | masks['rise'] | masks['fall'] | masks['reversal']
    masks['fast'] = masks['rise'] | masks['fall'] | masks['reversal']
    masks['outside_event'] = score & ~event
    return lam, masks


def rms(a):
    return float(np.sqrt(np.mean(a*a))) if len(a) else 0.


def objective(pred, case):
    e = pred-case['truth']
    return sum(rms(e[case['masks'][s]]) for s in ('plateau', 'dynamic'))


def stats(e):
    if not len(e):
        return {'n': 0, 'rms': None, 'median_abs': None, 'mean': None, 'p05': None}
    if not np.all(np.isfinite(e)):
        raise ValueError('missing estimate on scored ticks')
    return {'n': len(e), 'rms': rms(e), 'median_abs': float(np.median(np.abs(e))),
            'mean': float(np.mean(e)), 'p05': float(np.percentile(e, 5))}


def sustained(e, mask):
    n = 13  # ceil(0.1/TC), 0.104 s; complete contiguous valid windows only.
    ok = mask & np.isfinite(e)
    count = np.convolve(ok.astype(float), np.ones(n), mode='valid')
    avg = np.convolve(np.where(ok, e, 0.), np.ones(n)/n, mode='valid')
    return max(0., -float(np.min(avg[count == n]))) if np.any(count == n) else None


def load_cases(root):
    cases, checks = [], []
    for idx in range(42):
        path = root/'experiments/rossi_failure_study/results/huber_windows'/f'case_{idx:02d}.npz'
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path, allow_pickle=False) as f:
            log, raw, Y = f['log'], f['raw'], f['Y']
            meta = json.loads(str(f['meta']))
        # A0-A2 boundary: log/meta/truth never passed into an estimator.
        inputs = tuple(np.array(a, copy=True) for a in
                       (log[:, 0], log[:, 2], log[:, 4], raw, Y))
        data = preprocess(*inputs)
        truth, masks = masks_for(log, meta)
        a0, k0 = raw_arm(inputs[3])
        finite = np.isfinite(raw)
        checks.append({'case': idx, 'raw_min': float(np.min(raw[finite])),
                       'raw_max': float(np.max(raw[finite])),
                       'A0_exact_raw': bool(np.array_equal(k0[finite], raw[finite]))})
        cases.append({'index': idx, 'group': meta['group'], 'meta': meta, 'data': data,
                      'truth': truth, 'masks': masks, 'log': log,
                      'pred': {'A0': a0}, 'source_sha256': digest})
    return cases, checks


def select(cases, fused, workers):
    dev = [c for c in cases if c['group'] in BASE_GROUPS and c['index'] % 2 == 0]
    grid = GRID2 if fused else GRID1
    func = fused_filter if fused else probe_filter
    def one(pars):
        values = []
        for c in dev:
            out = func(c['data'], pars)
            pred = out[0] if fused else out
            value = objective(pred, c)
            values.append(value if np.isfinite(value) else float('inf'))
        return float(np.mean(values))
    scores = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, score in enumerate(pool.map(one, grid)):
            scores.append(score)
            if (i+1) % 16 == 0:
                print(f'{"A2" if fused else "A1"} grid {i+1}/{len(grid)}; best {min(scores):.6f}', flush=True)
    best = int(np.argmin(scores))
    if not np.isfinite(scores[best]):
        raise ValueError('all grid candidates failed')
    return grid[best], {'dev_cases': [c['index'] for c in dev], 'best_objective': scores[best],
                        'grid': [{'parameters': list(p), 'objective': s if np.isfinite(s) else None}
                                 for p, s in zip(grid, scores)], 'tie_rule': 'first in listed grid'}


def group_table(cases):
    table = {}
    for group in dict.fromkeys(c['group'] for c in cases):
        rows = [c for c in cases if c['group'] == group]
        gr = {'cases': [c['index'] for c in rows], 'segments': {}, 'sustained': {}}
        for segment in SEGMENTS:
            gr['segments'][segment] = {}
            for arm in ('A0', 'A1', 'A2', 'A3'):
                errors = np.concatenate([(c['pred'][arm]-c['truth'])[c['masks'][segment]] for c in rows])
                gr['segments'][segment][arm] = stats(errors)
        for arm in ('A0', 'A1', 'A2', 'A3'):
            values = [c['sustained'][arm] for c in rows if c['sustained'][arm] is not None]
            gr['sustained'][arm] = {'median': float(np.median(values)) if values else None,
                                    'worst': max(values) if values else None}
        table[group] = gr
    return table


def p4_fit(case):
    log = case['log']
    physical = (log[:, 9] > 0) & (log[:, 7] > .5)
    ratio = np.full(len(log), np.nan)
    ratio[physical] = log[physical, 7]/log[physical, 9]
    xx = np.where(physical, log[:, 2]+log[:, 6], np.nan)
    def nanmean_window(a, n, center=False):
        good = np.isfinite(a)
        sums = np.r_[0., np.cumsum(np.where(good, a, 0.))]
        counts = np.r_[0, np.cumsum(good)]
        j = np.arange(len(a))
        lo = np.maximum(0, j-(n//2 if center else n-1))
        hi = np.minimum(len(a), j+(n//2+1 if center else 1))
        nn = counts[hi]-counts[lo]
        result = np.full(len(a), np.nan)
        np.divide(sums[hi]-sums[lo], nn, out=result, where=nn > 0)
        return result
    x, z = nanmean_window(xx, 6), nanmean_window(ratio, 6)
    x, z = x-nanmean_window(x, 251, True), z-nanmean_window(z, 251, True)
    ok = physical & np.isfinite(x+z)
    xc, zc = x[ok]-np.mean(x[ok]), z[ok]-np.mean(z[ok])
    slope = float(xc@zc/(xc@xc))
    fit = 1/slope
    c = case['c2']
    learned = float(np.exp(c[np.flatnonzero(np.isfinite(c))[-1]]))
    return {'g_final': learned, 'g_fit': fit, 'slope': slope, 'n_fit': int(np.sum(ok)),
            'factor': max(learned/fit, fit/learned) if fit > 0 else None,
            'within_factor_1_5': bool(fit > 0 and max(learned/fit, fit/learned) <= 1.5)}


def predictions(cases):
    byidx = {c['index']: c for c in cases}
    def crossing(idx):
        c = byidx[idx]
        mask = c['masks']['fast']
        r1, r2 = [rms((c['pred'][a]-c['truth'])[mask]) for a in ('A1', 'A2')]
        return {'case': idx, 'test': bool(idx % 2), 'A1_rms': r1, 'A2_rms': r2,
                'ratio': r2/r1 if r1 else None, 'A2_sustained': c['sustained']['A2'],
                'met': bool(r2 <= .5*r1 and c['sustained']['A2'] <= .12)}
    p1 = [crossing(i) for i in (36, 37, 39, 40)]
    test = [c for c in cases if c['group'] in BASE_GROUPS and c['index'] % 2]
    p2 = []
    for group in BASE_GROUPS:
        rows = [c for c in test if c['group'] == group]
        rr = [rms(np.concatenate([(c['pred'][a]-c['truth'])[c['masks']['plateau']] for c in rows]))
              for a in ('A1', 'A2')]
        p2.append({'group': group, 'A1': rr[0], 'A2': rr[1], 'delta': rr[1]-rr[0], 'met': rr[1] <= rr[0]+.01})
    p3 = []
    subsets = [('A', [c for c in test if c['group'] == 'A']),
               ('A_slow_thin', [c for c in test if c['group'] == 'A_slow_thin']),
               ('A_large_38', [byidx[38]]), ('A_large_41', [byidx[41]])]
    for name, rows in subsets:
        rr = [rms(np.concatenate([(c['pred'][a]-c['truth'])[c['masks']['all']] for c in rows]))
              for a in ('A1', 'A2')]
        p3.append({'subset': name, 'A1': rr[0], 'A2': rr[1], 'delta': rr[1]-rr[0], 'met': abs(rr[1]-rr[0]) <= .02})
    p4 = [{'case': c['index'], 'group': c['group'], 'test': bool(c['index'] % 2), **c['gfit']} for c in cases]
    p5 = []
    for group in ('B', 'A+B'):
        rows = [c for c in test if c['group'] == group]
        rr = [rms(np.concatenate([(c['pred'][a]-c['truth'])[c['masks']['event']] for c in rows]))
              for a in ('A1', 'A2')]
        p5.append({'group': group, 'A1_event': rr[0], 'A2_event': rr[1], 'ratio': rr[1]/rr[0], 'met': rr[1] > rr[0]})
    return {'P1': {'met': all(r['met'] for r in p1), 'cases': p1},
            'P2': {'met': all(r['met'] for r in p2), 'groups': p2},
            'P3': {'met': all(r['met'] for r in p3), 'subsets': p3, 'segment': 'all scored ticks'},
            'P4': {'met_test_base_groups': all(r['within_factor_1_5'] for r in p4
                                              if r['test'] and r['group'] in BASE_GROUPS),
                   'cases': p4, 'evaluation': 'final episode estimate after repeated crossings'},
            'P5': {'met': all(r['met'] for r in p5), 'groups': p5}}


def make_figures(cases, root, qa_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                         'axes.spines.right': False, 'axes.spines.top': False,
                         'legend.frameon': False})
    out = root/'research/figures/dual_channel'
    out.mkdir(parents=True, exist_ok=True)
    for idx in (36, 39, 0, 6, 12):
        c = cases[idx]
        t = c['data'][0]
        fig, ax = plt.subplots(figsize=(7.1, 2.8))
        fig.subplots_adjust(left=.09, right=.985, bottom=.19, top=.78)
        # ln(k_eff) is undefined at zero stiffness: show gaps, not a log floor.
        truth_display = np.where(c['log'][:, 9] > 0, c['truth'], np.nan)
        ax.plot(t, truth_display, color='#222222', lw=1.2, label='Truth')
        for arm, color, width in [('A0', '#AAAAAA', .55), ('A1', '#3977A8', .85), ('A2', '#CE7435', 1.)]:
            ax.plot(t, c['pred'][arm], color=color, lw=width, label=arm)
        ax.fill_between(t, c['pred']['A2']-c['sd2'], c['pred']['A2']+c['sd2'],
                        color='#CE7435', alpha=.18, linewidth=0, label='A2 ± 1 SD')
        ax.set(xlabel='Time (s)', ylabel='ln k (k in N/m)', xlim=(t[0], t[-1]))
        ax.set_title(f'Case {idx}: {c["group"]} ({"test" if idx % 2 else "development"})', loc='left', pad=28)
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.015), ncol=5, fontsize=7, borderaxespad=0)
        fig.savefig(out/f'case_{idx}.png', dpi=300)
        fig.savefig(qa_dir/f'case_{idx}.pdf')
        fig.savefig(qa_dir/f'case_{idx}.svg')
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--synthetic-only', action='store_true')
    parser.add_argument('--workers', type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--qa-dir', type=Path)
    args = parser.parse_args()
    root, output = args.project_root, args.output_root or args.project_root
    started = time.perf_counter()
    synthetic = synthetic_check()
    if args.synthetic_only:
        return
    cases, raw_checks = load_cases(root)
    print('Loaded 42 read-only recordings; groups: '+str([(c['index'], c['group']) for c in cases]), flush=True)
    p1, selection1 = select(cases, False, args.workers)
    p2, selection2 = select(cases, True, args.workers)
    print(f'Selected A1={p1}, A2={p2}', flush=True)
    for c in cases:
        c['pred']['A1'] = probe_filter(c['data'], p1)
        a2, sd, cc, skips, resets = fused_filter(c['data'], p2)
        c['pred']['A2'], c['sd2'], c['c2'] = a2, sd, cc
        data3 = list(c['data'])
        data3[1] = trailing(c['log'][:, 2]+c['log'][:, 6], 6)
        a3, _, _, skips3, resets3 = fused_filter(tuple(data3), p2)
        c['pred']['A3'] = a3
        c['diagnostics'] = {'A2_mode_probe_skips': skips, 'A3_mode_probe_skips': skips3,
                            'A2_initializations': resets, 'A3_initializations': resets3}
        c['sustained'] = {arm: sustained(pred-c['truth'], c['masks']['outside_event'])
                          for arm, pred in c['pred'].items()}
        c['gfit'] = p4_fit(c)
    test = [c for c in cases if c['group'] in BASE_GROUPS and c['index'] % 2]
    report_only = [c for c in cases if c['group'] not in BASE_GROUPS]
    result = {'selected': {'A1': dict(zip(('q0', 'q1', 'p_stay'), p1)),
                           'A2_A3': dict(zip(('q0', 'q1', 'p_stay', 'q_b', 'q_c'), p2))},
              'selection': {'A1': selection1, 'A2': selection2},
              'test_tables': group_table(test), 'report_only_tables': group_table(report_only),
              'A_large_per_case': {str(c['index']): {'test': bool(c['index'] % 2),
                   **group_table([c])[c['group']]} for c in cases if c['group'] == 'A_large'},
              'predictions': predictions(cases),
              'checks': {'synthetic_before_recordings': synthetic, 'A0_raw_per_case': raw_checks,
                         'A0_all_exact': all(r['A0_exact_raw'] for r in raw_checks),
                         'truth_isolation': 'A0 accepts raw only; A1/A2 accept tuple(t,ys,Fs,m,R,valid,finite_raw), made solely from copied t,y,F_sensor,raw,Y. No log/meta/truth argument or plant import.'},
              'conventions': {'means': 'causal trailing, available prefix at startup; truth segmentation same',
                 'contact_loss_output': 'hold last estimate while not updating; reinitialize at next usable tick',
                 'sustained_window_ticks': 13, 'sustained_window_s': 13*TC,
                 'sustained': 'complete contiguous scored windows outside events only',
                 'group_metrics': 'pooled scored ticks, not average of per-case metrics',
                 'plot_truth': 'all ticks, with gaps where k_eff<=0 (log undefined)',
                 'reversal': 'sign change with abs(dl)>1 within 0.15s on each side, then ±0.1s',
                 'P4': 'OLS with intercept; six-tick mean true fax/k_eff on six-tick mean(y+u_true); subtract centred 251-tick (2.008s) mean; k_eff>0 and fax>.5 only, finite-sample normalized windows and edges',
                 'initialization': 'first usable probe tick, no double update; initial log-k variance R; c variance (ln3)^2',
                 'probe_likelihood': 'sequential conditional Gaussian log likelihoods summed; ineligible measurement contributes zero',
                 'A2_uncertainty': 'mixture covariance projected by H=(0,0,0,-1,1); conditional on observed Fs',
                 'synthetic': 'fixed grid parameters independent of recordings; 12s; alternating ramps every 1.2s, first at 0.8s'},
              'cases': [{'index': c['index'], 'group': c['group'], 'source_sha256': c['source_sha256'],
                          'diagnostics': c['diagnostics'], 'sustained': c['sustained'], 'g_fit': c['gfit'],
                          'n_scored': {s: int(np.sum(c['masks'][s])) for s in SEGMENTS}}
                         for c in cases]}
    if args.qa_dir:
        args.qa_dir.mkdir(parents=True, exist_ok=True)
        make_figures(cases, output, args.qa_dir)
    else:
        with tempfile.TemporaryDirectory(prefix='dual_channel_qa_') as temp:
            make_figures(cases, output, Path(temp))
    result['runtime_s'] = time.perf_counter()-started
    dest = output/'experiments/rossi_failure_study/results/dual_channel_replay.json'
    dest.parent.mkdir(parents=True, exist_ok=True)
    assert all(hashlib.sha256((root/'experiments/rossi_failure_study/results/huber_windows'/f'case_{c["index"]:02d}.npz').read_bytes()).hexdigest()
               == c['source_sha256'] for c in cases)
    result['checks']['source_hashes_unchanged'] = True
    dest.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Written '+str(dest)+f'; runtime {result["runtime_s"]:.1f}s', flush=True)
    print('Predictions: '+json.dumps(result['predictions']), flush=True)


if __name__ == '__main__':
    main()
