"""Step 1 of the V2 review: re-verify V2 with the exact h1, the note's own solver, a small calibration on half the cases, and per-segment scoring.

Inputs: results/huber_windows (re-captured with H1, H10 = the time-weighted regressor of eq. (1), noisy and noise-free). No new plant runs.
Estimators on identical windows:
  const   solve_v2_window with sv2 = 0, i.e. Huber-D4 on 24 ticks (the note's own limit)
  free    sv2 = inf (slope unpenalized)
  oracle  sv2 = v_true^2, the best any rule for S_v^2 could do
  V2      eqs. (11)-(13), parameters chosen on the development half (even case index within each group), reported on the other half
Segments from k_eff (0.2 s boxcar rate r): plateau as in table 2; rising r > .5/s; falling r < -.5/s; reversal within 0.1 s of a sign change between
runs of |d ln k| >= .15; event = scene event window. Events are excluded from the k-based segments.
For the rule grid, the Huber weights of the sv2 = 0 solve are reused; the chosen configuration is re-run with the full solver as a check.
Writes results/v2_verify.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import huber_prior_check as H

HERE = Path(__file__).resolve().parent
TC, NW = .008, 24
RHO = H.RHO
SIG = H.SIGMA
SD = .04/math.sqrt(RHO)        # the note's initial s_d
UHI = 8.
D = np.diff(np.eye(NW), 4, axis=0)


def solve_v2_window(h0, h1, z, *, rho, sd, sigma_z, sv2, huber_c=1.345, max_iter=50, tol=1e-6, k_bounds=(200.0, 300000.0)):
    """Appendix A of the note, unchanged except that it returns the final weights' smoother for the trend test."""
    n = len(z)
    Dn = np.diff(np.eye(n), 4, axis=0)
    w = np.ones(n-4)
    rhs = np.column_stack([h0, h1, z])
    lam = np.inf if sv2 == 0 else sigma_z**2/sv2
    for iteration in range(max_iter):
        w_used = w.copy()
        S_rhs = np.linalg.solve(np.eye(n)+rho*(Dn.T*w_used)@Dn, rhs)
        M_rhs = rhs-S_rhs
        mh0, mh1, mz = M_rhs.T
        a, b, c1 = float(h0@mh0), float(h0@mh1), float(h1@mh1)
        p0, p1 = float(h0@mz), float(h1@mz)
        if sv2 == 0:
            k, slope = p0/a, 0.0
        else:
            cp = c1+lam
            info = a-b*b/cp
            q = (mh0-(b/cp)*mh1)/info
            k = float(q@z); slope = (p1-b*k)/cp
        k = float(np.clip(k, *k_bounds)) if np.isfinite(k) else np.nan
        d = S_rhs[:, 2]-k*S_rhs[:, 0]-slope*S_rhs[:, 1]
        u = Dn@d/sd
        w_new = np.minimum(1.0, huber_c/np.maximum(abs(u), 1e-15))
        if np.max(abs(w_new-w_used)) <= tol:
            break
        w = w_new
    return k, slope, w_used, dict(a=a, b=b, c1=c1, p0=p0, p1=p1)


def m_of(w):
    return np.eye(NW)-np.linalg.inv(np.eye(NW)+RHO*(D.T*w)@D)


def closed(st, sv2):
    """k for a given sv2 with fixed weights, eqs. (6)-(7)."""
    a, b, c1, p0, p1 = st
    if sv2 <= 0:
        return p0/a
    if not np.isfinite(sv2):
        return (p0*c1-b*p1)/(a*c1-b*b)
    cp = c1+SIG**2/sv2
    return (p0-(b/cp)*p1)/(a-b*b/cp)


def segments(t, k, events, contact):
    ell = np.log(np.maximum(k, 1.))
    s = uniform_filter1d(ell, 25, mode='nearest'); r = np.gradient(s, t)
    ev = np.zeros(len(t), bool)
    for e in events:
        ev |= (t >= e['start']) & (t <= e['end']+.3)
    ok = contact & (t >= 1.)
    sign = np.where(r > .2, 1, np.where(r < -.2, -1, 0))
    rev = np.zeros(len(t), bool)
    # Runs of one sign, merged across gaps under 0.1 s; a reversal where a run of |d ln k| >= .15 meets one of the other sign.
    runs = []; cur = 0; a = 0
    for j in range(len(t)):
        if sign[j] != cur:
            if cur:
                runs.append([a, j, cur])
            a, cur = j, sign[j]
    if cur:
        runs.append([a, len(t), cur])
    merged = []
    for q in runs:
        if merged and merged[-1][2] == q[2] and t[q[0]]-t[merged[-1][1]-1] < .1:
            merged[-1][1] = q[1]
        else:
            merged.append(q)
    big = [q for q in merged if abs(s[q[1]-1]-s[q[0]]) >= .15]
    for p, q in zip(big[:-1], big[1:]):
        if p[2] == -q[2] and t[q[0]]-t[p[1]-1] < .3:
            mid = .5*(t[p[1]-1]+t[q[0]])
            rev |= np.abs(t-mid) <= .1
    out = dict(rising=ok & ~ev & ~rev & (r > .5), falling=ok & ~ev & ~rev & (r < -.5), reversal=ok & ~ev & rev, event=ok & ev)
    out['plateau'] = ok & ~ev & ~rev & (np.abs(r) < .05)
    return out


def case(path):
    data = H.load(path); meta = data['meta']; log, _ = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    seg = segments(t, k, meta['event_info'], k > 200.)
    rec = {q: np.full(n, np.nan) for q in ('k', 'vtrue', 'c', 'f', 'o', 'c0', 'f0', 'o0', 'vt', 'svt', 'vt0', 'kref', 'h1err')}
    st = np.full((n, 5), np.nan); st0 = np.full((n, 5), np.nan)
    for j in range(NW+1, n):
        Y, Y1, z = data['Y'][j], data['H1'][j], data['z'][j]
        Y0, Y10, z0 = data['Y0'][j], data['H10'][j], data['z0'][j]
        if not (np.isfinite(Y).all() and np.isfinite(Y1).all() and np.isfinite(Y0).all() and np.isfinite(Y10).all()) or k[j] <= 200.:
            continue
        h0, h1, zz = Y[-NW:], Y1[-NW:], z[-NW:]
        r = np.arange(j-NW, j); kr = .5*(k[r]+k[r+1])
        if np.any(kr <= 200.):
            continue
        tt = .5*(t[r]+t[r+1])-t[j-1]
        vtrue = float(np.polyfit(tt, kr, 1)[0])
        kc, _, w0, s = solve_v2_window(h0, h1, zz, rho=RHO, sd=SD, sigma_z=SIG, sv2=0.)
        st[j] = [s['a'], s['b'], s['c1'], s['p0'], s['p1']]
        rec['c'][j] = kc
        rec['f'][j] = closed(st[j], np.inf); rec['o'][j] = solve_v2_window(h0, h1, zz, rho=RHO, sd=SD, sigma_z=SIG, sv2=vtrue**2)[0]
        # Trend test, eq. (11), on the constant solve's weights.
        M0 = m_of(w0)
        a0, b0 = h0@M0@h0, h0@M0@h1
        g = h1-(b0/a0)*h0; Jv = g@M0@g; lv = M0@g/Jv
        rec['vt'][j] = lv@zz; rec['svt'][j] = SIG*math.sqrt(lv@lv)
        k0c, _, w00, s0 = solve_v2_window(Y0[-NW:], Y10[-NW:], z0[-NW:], rho=RHO, sd=SD, sigma_z=SIG, sv2=0.)
        st0[j] = [s0['a'], s0['b'], s0['c1'], s0['p0'], s0['p1']]
        rec['c0'][j] = k0c; rec['f0'][j] = closed(st0[j], np.inf)
        rec['o0'][j] = solve_v2_window(Y0[-NW:], Y10[-NW:], z0[-NW:], rho=RHO, sd=SD, sigma_z=SIG, sv2=vtrue**2)[0]
        M00 = m_of(w00); g0 = Y10[-NW:]-((Y0[-NW:]@M00@Y10[-NW:])/(Y0[-NW:]@M00@Y0[-NW:]))*Y0[-NW:]
        rec['vt0'][j] = (M00@g0/(g0@M00@g0))@z0[-NW:]
        rec['k'][j] = k[j]; rec['vtrue'][j] = vtrue
        # How far the interval-centre approximation of the previous check was from the exact h1.
        tc = (np.arange(NW)-(NW-1))*TC-TC+.0035
        rec['h1err'][j] = np.linalg.norm(tc*h0-h1)/np.linalg.norm(h1)
    return dict(index=meta['index'], group=meta['group'], seg=seg, rec=rec, st=st, st0=st0)


def run_rule(res, kappa, tau, beta, use_kref=True):
    """Eq. (13) over each case, every 6 ticks, from the next tick; K_ref = previous accepted (here: previous V2 output)."""
    eta_down = 1-math.exp(-6*TC/tau)
    for r in res:
        rec, st = r['rec'], r['st']
        n = len(rec['k']); out = np.full(n, np.nan); sv2 = 0.; kref = np.nan; sv_trace = np.full(n, np.nan)
        for j in range(n):
            if not np.isfinite(st[j, 0]):
                sv2 = 0.; continue
            kj = closed(st[j], sv2)
            kj = float(np.clip(kj, 200., 3e5))
            out[j] = kj; sv_trace[j] = sv2
            kr = kref if np.isfinite(kref) else kj
            if j % 6 == 0:
                sv = math.sqrt(rec['svt'][j]**2+(kr*beta)**2)
                R = max(rec['vt'][j]**2-kappa*sv*sv, 0.)
                if sv > kr*UHI:
                    pass
                else:
                    eta = 1. if R > sv2 else eta_down
                    sv2 = min((1-eta)*sv2+eta*R, (kr*UHI)**2)
            kref = kj
        r['rec']['v2'] = out; r['rec']['sv2'] = sv_trace


SEGS = ('plateau', 'rising', 'falling', 'reversal', 'event')


def score(res, key, groups=None, cases=None):
    out = {}
    for s in SEGS:
        e = []
        for r in res:
            if (groups and r['group'] not in groups) or (cases is not None and r['index'] not in cases):
                continue
            m = r['seg'][s] & np.isfinite(r['rec'][key]) & np.isfinite(r['rec']['k'])
            e.append(np.log(np.maximum(r['rec'][key][m], 200.)/r['rec']['k'][m]))
        e = np.concatenate(e) if e else np.array([])
        if len(e) < 20:
            out[s] = None; continue
        out[s] = dict(n=len(e), med=float(np.median(np.abs(e))), p90=float(np.percentile(np.abs(e), 90)), rms=float(np.sqrt(np.mean(e*e))), bias=float(np.median(e)), low=float(np.mean(e < math.log(.85))), high=float(np.mean(e > math.log(1.2))))
    return out


def within_sd(res, key, groups):
    vals = []
    for r in res:
        if r['group'] not in groups:
            continue
        m = r['seg']['plateau'] & np.isfinite(r['rec'][key])
        e = np.log(np.maximum(r['rec'][key][m], 200.)/r['rec']['k'][m])
        if len(e) > 30:
            q = np.percentile(e, [25, 75]); vals.append((q[1]-q[0])/1.349)
    return float(np.median(vals)) if vals else None


def objective(sc):
    return float(np.mean([sc[s]['rms'] for s in SEGS if sc[s]]))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rep = {}
    h1err = np.concatenate([r['rec']['h1err'] for r in res]); h1err = h1err[np.isfinite(h1err)]
    rep['h1_centre_approx_rel_error'] = dict(median=float(np.median(h1err)), p90=float(np.percentile(h1err, 90)))
    print('previous interval-centre h1 vs exact, relative error median/P90', rep['h1_centre_approx_rel_error'])
    groups = sorted({r['group'] for r in res})
    dev = {r['index'] for r in res if r['index'] % 2 == 0}; test = {r['index'] for r in res} - dev
    # Calibration grid on the development half.
    grid = list(product((1., 2., 4., 9., 16.), (.048, .192, .768), (0., .5, 1., 2.)))
    best = None
    for kappa, tau, beta in grid:
        run_rule(res, kappa, tau, beta)
        sc = score(res, 'v2', cases=dev); obj = objective(sc)
        if best is None or obj < best[0]:
            best = (obj, kappa, tau, beta)
    rep['chosen'] = dict(objective=best[0], kappa=best[1], tau_release=best[2], beta_v=best[3])
    base_dev = objective(score(res, 'c', cases=dev))
    print('development: mean segment rms, constant', round(base_dev, 4), 'best V2', round(best[0], 4), 'at kappa, tau, beta', best[1:])
    run_rule(res, *best[1:])
    rep['test'] = {}
    print('test half, ln(k/k_eff) by segment: |e| med/P90, rms, bias, >1.2 share, <.85 share')
    for key, name in (('c', 'Huber-D4 (sv2 0)'), ('f', 'free slope'), ('o', 'oracle slope scale'), ('v2', 'V2 calibrated'), ('c0', 'noise-free const'), ('f0', 'noise-free free'), ('o0', 'noise-free oracle')):
        sc = score(res, key, cases=test); rep['test'][name] = sc
        print(f'  {name:20s}', ' | '.join(f"{s} {sc[s]['med']:.3f}/{sc[s]['p90']:.3f} {sc[s]['rms']:.3f} {sc[s]['bias']:+.3f} {sc[s]['high']*100:.0f}/{sc[s]['low']*100:.0f}%" if sc[s] else f'{s} -' for s in SEGS))
    rep['per_group'] = {}
    print('test half by group, rms by segment: const / V2 / oracle')
    for g in groups:
        a, b, o = (score(res, key, groups=[g], cases=test) for key in ('c', 'v2', 'o'))
        rep['per_group'][g] = dict(const=a, v2=b, oracle=o)
        print(f'  {g:12s}', ' | '.join(f"{s} {a[s]['rms']:.3f}/{b[s]['rms']:.3f}/{o[s]['rms']:.3f}" for s in SEGS if a[s] and b[s] and o[s]))
    print('plateau within-case sd, test half: const, free, V2')
    for g in groups:
        print(f'  {g:12s}', within_sd([r for r in res if r['index'] in test], 'c', [g]), within_sd([r for r in res if r['index'] in test], 'f', [g]), within_sd([r for r in res if r['index'] in test], 'v2', [g]))
    # Trend information, exact h1.
    rep['trend'] = {}
    print('trend statistic relative to k (1/s): noise sd, |true| P90, noise-free corr with true, S_v^2 share above its noise level on plateau')
    for g in groups:
        rr = [r for r in res if r['group'] == g]
        svt = np.concatenate([r['rec']['svt']/r['rec']['k'] for r in rr]); vt = np.concatenate([r['rec']['vtrue']/r['rec']['k'] for r in rr]); vt0 = np.concatenate([r['rec']['vt0']/r['rec']['k'] for r in rr])
        pl = np.concatenate([r['seg']['plateau'] for r in rr]); sv = np.concatenate([np.sqrt(r['rec']['sv2'])/r['rec']['k'] for r in rr]); sn = np.concatenate([r['rec']['svt']/r['rec']['k'] for r in rr])
        ok = np.isfinite(vt0) & np.isfinite(vt)
        rep['trend'][g] = dict(noise_sd=float(np.nanmedian(svt)), true_p90=float(np.nanpercentile(np.abs(vt), 90)), corr=float(np.corrcoef(vt0[ok], vt[ok])[0, 1]), plateau_sv_over_noise=float(np.nanmean((sv > .5*sn)[pl])))
        print(f'  {g:12s}', {q: round(v, 3) for q, v in rep['trend'][g].items()})
    # Check the fixed-weight shortcut on the chosen configuration with the full solver on a few cases.
    worst = 0.
    for r in res[:42:7]:
        data = H.load(H.WIN/f"case_{r['index']:02d}.npz")
        for j in range(NW+1, len(r['rec']['k']), 25):
            if np.isfinite(r['rec']['v2'][j]) and np.isfinite(r['rec']['sv2'][j]) and r['rec']['sv2'][j] > 0:
                kf = solve_v2_window(data['Y'][j][-NW:], data['H1'][j][-NW:], data['z'][j][-NW:], rho=RHO, sd=SD, sigma_z=SIG, sv2=r['rec']['sv2'][j])[0]
                worst = max(worst, abs(math.log(kf/r['rec']['v2'][j])))
    rep['fixed_weight_check_max_abs_dln'] = worst
    print('fixed-weight shortcut vs full solver, max |d ln k|', round(worst, 4))
    (HERE/'results/v2_verify.json').write_text(json.dumps(rep, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
