"""Offline replay of the stiffness estimators on recorded closed-loop data: projection against a partially linear model of d.

Same equation-error regression as common.window_fit, z = k Y + d + noise, with d given a prior instead of a polynomial:

  min_{k,d} ||z - kY - d||^2 + rho sum_i omega_i (Dn d)_i^2
  fixed omega: d = S(z - kY), S = (I + rho Dn' W Dn)^-1, k = Y'Mz / Y'MY, M = I - S, info = Y'MY/sigma^2
  Huber: omega_i = min(1, delta/|(Dn d)_i/q|), a few reweighting passes

M has frequency response rho t^(2n)/(1 + rho t^(2n)), t = 2 sin(theta/2); rho is set by its corner fc. As rho grows a degree n-1
polynomial d is recovered, so n = 4 contains the cubic projection of the same window. Optional memory passes y through the population
weak power-law response normalised at the probe frequency, so k stays the gain at that frequency.
Selection on seeds 390000, evaluation on 386000 (the seeds of the closed-loop rows), truth column used only for scoring.
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import signal

import sweeps as S      # sets the import paths used below
import common
import physics_plant as PP
from ridge_estimator import POP_ALPHA, POP_TAU
from failure_replay import Measured, cell_scene

HERE = Path(__file__).resolve().parent
LN12 = math.log(1.2)
REL_MAX = .2565
DELTA = 3.


def memory_response(y, Tc=.008, f0=20.8):
    """y through G(s) = (1 + sum alpha tau s/(1 + tau s))/G(j w0) magnitude, phase kept: the force the probe motion produces, per unit k."""
    z = np.exp(-1j*2*math.pi*f0*Tc); e_ = np.exp(-Tc/POP_TAU); b_ = POP_TAU/Tc*(1-e_)
    mag = abs(1+np.sum(POP_ALPHA*b_*(1-z)/(1-e_*z)))      # response of the discrete filter used below
    out = y.astype(float).copy()-y[0]
    total = out.copy()
    for a, t in zip(POP_ALPHA, POP_TAU):
        e = math.exp(-Tc/t); r = t/Tc*(1-e); b = [a*r, -a*r]; den = [1, -e]
        total = total+signal.lfilter(b, den, out)
    return y[0]+total/mag


def smoother_fit(Y, z, sigma, rho, huber, q, order):
    n = len(z); D = np.diff(np.eye(n), order, axis=0)
    w = np.ones(n-order); passes = 4 if huber else 1
    for it in range(passes):
        A = np.eye(n)+rho*(D.T*w)@D
        Sy = np.linalg.solve(A, Y); Sz = np.linalg.solve(A, z)
        info = float(Y@(Y-Sy))
        if info <= 0:
            return np.nan, np.inf, np.inf
        k = float(Y@(z-Sz))/info
        if huber and it < passes-1:
            w = np.minimum(1., DELTA/np.maximum(np.abs(D@(Sz-k*Sy))/q, 1e-9))
    res = z-k*Y
    chi = float(res@(res-np.linalg.solve(A, res)))/sigma**2/max(n-order-1, 1)
    return k, sigma/math.sqrt(info), chi


def window_arrays(pol, log, n, nw):
    j0 = max(pol.ndel+1, n-1-nw); m = n-2-j0+1
    if m < 5:
        return None
    te, ye = log[n-1, 0], log[n-1, 2]
    X, Y, Z, keep = common._build(log, j0, m, 0, pol.ns, pol.eb, pol.afc, pol.ndel, te, ye, nw*pol.ns*1e-3, pol.params['enc_y'])
    if keep.sum() < m or m < nw:
        return None
    return Y, Z


def run_estimators(log, pol, spec):
    n_tick = len(log); out = np.full((n_tick, 3), np.nan)
    sigma = pol.params['noise']
    if spec['kind'] == 'projection':
        src = log
        if spec.get('memory'):
            src = log.copy(); src[:, 2] = memory_response(log[:, 2])
        for n in range(spec['nw']+pol.ndel+3, n_tick):
            f = common.window_fit(src, n, spec['nw'], spec['degree'], sigma, pol.ns, pol.eb, pol.afc, pol.ndel, pol.params['enc_y'])
            if np.isfinite(f['k']) and f['info'] > 0:
                out[n-1] = (f['k'], 1./math.sqrt(f['info']), f['chi2'])   # the window ends at row n-1
        return out
    use = log
    if spec.get('memory'):
        use = log.copy(); use[:, 2] = memory_response(log[:, 2])
    for n in range(spec['nw']+pol.ndel+3, n_tick):
        a = window_arrays(pol, use, n, spec['nw'])
        if a is None:
            continue
        k, s, chi = smoother_fit(a[0], a[1], sigma, spec['rho'], spec['huber'], spec['q'], spec['order'])
        out[n-1] = (k, s, chi)
    return out


def policy(controls, est):
    return S.policy('Projection', controls, est)


def score(log, est, sc, f_exc=20.8):
    t = log[:, 0]; kt = log[:, 8]; F = log[:, 3]
    ok = (t >= 1.) & (F > .5) & (kt > 0) & np.isfinite(est[:, 0]) & (est[:, 0] > 0)
    if ok.sum() < 20:
        return None
    e = np.log(est[ok, 0]/kt[ok]); rel = est[ok, 1]/np.maximum(est[ok, 0], 1.)
    acc = rel <= REL_MAX
    ev = np.zeros(len(t), bool)
    for _, a0, a1 in sc.get('events', []):
        ev |= (t >= a0) & (t <= a1)
    out = dict(all_abs=np.abs(e).tolist(), acc_rate=float(acc.mean()), acc_abs=np.abs(e[acc]).tolist(), acc_over=float(np.mean(np.abs(e[acc]) > LN12)) if acc.any() else None,
               bias=float(np.median(e)), cover=float(np.mean(np.abs(e) <= 2*est[ok, 1]/np.maximum(est[ok, 0], 1.))),
               ev_abs=np.abs(e[ev[ok]]).tolist(), ev_acc=float(acc[ev[ok]].mean()) if ev[ok].any() else None)
    rec = []
    for _, a0, _ in sc.get('events', []):
        after = np.flatnonzero(ok & (t > a0))
        if not len(after):
            continue
        idx = np.flatnonzero(ok)
        good = (t[idx] > a0+.05) & (rel <= REL_MAX) & (np.abs(e) <= LN12)
        bad = (t[idx] > a0) & ((rel > REL_MAX) | (np.abs(e) > LN12))
        fb = np.flatnonzero(bad)
        if len(fb):
            later = np.flatnonzero(good & (t[idx] > t[idx][fb[0]]))
            rec.append(float(t[idx][later[0]]-a0) if len(later) else float('nan'))
    out['recover'] = rec
    return out


def job(a):
    c, t, i, base, controls, est, specs = a
    sc = cell_scene(c, t, i, base); seed = base+sum(ord(x) for x in c+t)*100+i
    pol = policy(controls, est)
    log = PP.episode(sc, Measured(pol), seed)
    if log is None:
        return []
    out = []
    for name, spec in specs.items():
        r = score(log, run_estimators(log, pol, spec), sc)
        if r:
            out.append(dict(name=name, cond=c, tissue=t, i=i, **r))
    return out


def pool(rows, name, pick=lambda r: True):
    rr = [r for r in rows if r['name'] == name and pick(r)]
    if not rr:
        return None
    a = np.concatenate([np.array(r['all_abs']) for r in rr]); acc = np.concatenate([np.array(r['acc_abs']) for r in rr if r['acc_abs']]) if any(r['acc_abs'] for r in rr) else np.array([np.nan])
    cell_p90 = {}
    for r in rr:
        cell_p90.setdefault((r['cond'], r['tissue']), []).append(r['all_abs'])
    p90 = [np.percentile(np.concatenate(v), 90) for v in cell_p90.values()]
    ev = [x for r in rr for x in r['ev_abs']]; rec = [x for r in rr for x in r['recover']]
    return dict(mean_cell_p90=float(np.mean(p90)), p50=float(np.percentile(a, 50)), p90=float(np.percentile(a, 90)), acc_rate=float(np.mean([r['acc_rate'] for r in rr])),
                acc_p90=float(np.nanpercentile(acc, 90)), acc_over=float(np.nanmean([r['acc_over'] for r in rr if r['acc_over'] is not None])),
                bias=float(np.median([r['bias'] for r in rr])), cover=float(np.mean([r['cover'] for r in rr])),
                ev_p90=float(np.percentile(ev, 90)) if ev else None, recover_med=float(np.nanmedian(rec)) if rec else None, recover_never=float(np.mean(np.isnan(rec))) if rec else None)


def run(base, n_scene, specs, controls, est):
    args = [(c, t, i, base, controls, est, specs) for c, t in S.CELLS for i in range(n_scene)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        return [r for g in ex.map(job, args, chunksize=1) for r in g]


def show(rows, names):
    groups = (('all', lambda r: True), ('soft', lambda r: r['tissue'] == 'soft'), ('mid', lambda r: r['tissue'] == 'mid'), ('hard', lambda r: r['tissue'] == 'hard'),
              ('event', lambda r: r['cond'] in ('flinch', 'cough')), ('fast', lambda r: r['cond'] == 'shift_fast'), ('smooth', lambda r: r['cond'] in ('breath', 'shift_slow')),
              ('rib', lambda r: r['cond'] == 'rib'), ('combo', lambda r: r['cond'] == 'combo'))
    for g, pick in groups:
        print(f'[{g}]')
        for n in names:
            s = pool(rows, n, pick)
            if s:
                print(f"  {n:18s} cellP90 {s['mean_cell_p90']:.3f} p50 {s['p50']:.3f} p90 {s['p90']:.3f} acc {s['acc_rate']:.2f} accP90 {s['acc_p90']:.3f} over {s['acc_over']:.3f} bias {s['bias']:+.3f} cover {s['cover']:.2f}"
                      + (f" evP90 {s['ev_p90']:.3f} rec {s['recover_med']:.2f}/never {s['recover_never']:.2f}" if s['recover_med'] is not None else ''))


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']['Projection']['best']['params']
    q = {2: .018, 3: .021, 4: .037}   # measured core of Dn d in newtons, from results/d_structure.json, independent of these seeds
    def ridge(nw, fc, order=4, hub=False, memory=False):
        return dict(kind='smooth', nw=nw, order=order, rho=(2*math.sin(math.pi*fc*.008))**(-2*order), huber=hub, q=q[order], memory=memory)
    specs = {'P3_16': dict(kind='projection', nw=16, degree=3), 'P5_16': dict(kind='projection', nw=16, degree=5), 'P5_24': dict(kind='projection', nw=24, degree=5),
             'D2_24_f14': ridge(24, 14., 2), 'R24': ridge(24, 14.), 'R24_huber': ridge(24, 14., hub=True), 'R24+mem': ridge(24, 14., memory=True), 'R48+mem': ridge(48, 18., memory=True),
             'P3_16+mem': dict(kind='projection', nw=16, degree=3, memory=True)}
    rows = run(386000, 5, specs, controls, est)
    score_sel, best_l2, best_h = {}, 'R24', 'R24_huber'
    (HERE/'results'/'d_estimator_replay.json').write_text(json.dumps(dict(selection=score_sel, chosen=[best_l2, best_h],
        summary={n: pool(rows, n) for n in specs}), indent=1, allow_nan=True), encoding='utf-8')
    show(rows, list(specs))


if __name__ == '__main__':
    main()
