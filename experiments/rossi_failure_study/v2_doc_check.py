"""Data check of active_excitation_joint_theory_v2_unified.md on recorded windows (results/huber_windows), no new simulation.

Approximation: h1_j = (t_j + 3.5 ms - T) h0_j, the sub-sample time weighting replaced by the interval centre; D4 with W = I
(Huber weights are inactive in 86-100% of windows). Truth k_eff is for scoring and for the oracle slope only.
Predictions written before running (in the research note): P1 free slope doubles plateau sd; P2 noise-free slope removes the linear part
of the within-window change but not cough onset or high-stiffness content; P3 even an oracle slope scale gains little in A, some in A_large;
P4 sd of v~/k about 3-4 /s; P5 v~ reacts to non-stiffness content in cough and combo; P6 a 20% target needs at least today's amplitude.
Writes results/v2_doc_check.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H

HERE = Path(__file__).resolve().parent
TC, NW = .008, 24
SIG = H.SIGMA
D = np.diff(np.eye(NW), 4, axis=0)
M = np.eye(NW)-np.linalg.inv(np.eye(NW)+H.RHO*D.T@D)
TIME = (np.arange(NW)-(NW-1))*TC-TC+.0035       # interval centre relative to the window end T
KAPPA, TAU_REL, UHI = 4., .192, 8.


def solve(h0, h1, z, sv2):
    a, b, c1 = h0@M@h0, h0@M@h1, h1@M@h1
    p0, p1 = h0@M@z, h1@M@z
    if sv2 <= 0:
        q = M@h0/a; k = q@z; return k, 0., q
    lam = SIG**2/sv2 if np.isfinite(sv2) else 0.
    cp = c1+lam; info = a-b*b/cp
    q = (M@h0-(b/cp)*(M@h1))/info
    k = q@z
    return k, (p1-b*k)/cp, q


def trend(h0, h1, z):
    a0, b0 = h0@M@h0, h0@M@h1
    g = h1-(b0/a0)*h0; Jv = g@M@g
    lv = M@g/Jv
    return lv@z, SIG*math.sqrt(lv@lv), Jv


def case(path):
    data = H.load(path); meta = data['meta']; log, masks = H.masks_for(data)
    n = len(log); k = log[:, 26]; t = log[:, 0]
    ev = np.zeros(n, bool); onset = np.zeros(n, bool)
    for e in meta['event_info']:
        ev |= (t >= e['start']) & (t <= e['end']+.3)
        onset |= (t >= e['start']) & (t <= e['start']+.15)
    keys = ('c', 'f', 'o', 'r', 'c0', 'f0', 'o0', 'vt', 'vt0', 'svt', 'vtrue', 'sv2', 'sd1', 'k')
    out = {q: np.full(n, np.nan) for q in keys}
    sv2 = 0.
    for j in range(n):
        Y, z = data['Y'][j], data['z'][j]; Y0, z0 = data['Y0'][j], data['z0'][j]
        if not (np.isfinite(Y[-NW:]).all() and np.isfinite(Y0[-NW:]).all()) or k[j] <= 200. or j < NW+1:
            sv2 = 0.; continue
        r = np.arange(j-NW, j); kr = .5*(k[r]+k[r+1])
        if np.any(kr <= 200.):
            sv2 = 0.; continue
        h0, h1 = Y[-NW:], TIME*Y[-NW:]
        g0, g1 = Y0[-NW:], TIME*Y0[-NW:]
        zz, zz0 = z[-NW:], z0[-NW:]
        vtrue = float(np.polyfit(TIME, kr, 1)[0])
        kc, _, q = solve(h0, h1, zz, 0.)
        kf, _, _ = solve(h0, h1, zz, np.inf)
        ko, _, _ = solve(h0, h1, zz, vtrue**2)
        kr_, _, _ = solve(h0, h1, zz, sv2)
        out['c'][j], out['f'][j], out['o'][j], out['r'][j] = kc, kf, ko, kr_
        out['c0'][j] = solve(g0, g1, zz0, 0.)[0]; out['f0'][j] = solve(g0, g1, zz0, np.inf)[0]; out['o0'][j] = solve(g0, g1, zz0, vtrue**2)[0]
        vt, svt, _ = trend(h0, h1, zz); vt0 = trend(g0, g1, zz0)[0]
        out['vt'][j], out['vt0'][j], out['svt'][j], out['vtrue'][j], out['k'][j] = vt, vt0, svt, vtrue, k[j]
        out['sd1'][j] = SIG*math.sqrt(q@q)/kc if kc > 0 else np.nan
        # V2 rule (13), applied every 6 ticks, uncalibrated: kappa 4, release 0.192 s, cap UHI k.
        if j % 6 == 0:
            R = max(vt*vt-KAPPA*svt*svt, 0.)
            eta = 1. if R > sv2 else 1-math.exp(-6*TC/TAU_REL)
            sv2 = min((1-eta)*sv2+eta*R, (k[j]*UHI)**2)
        out['sv2'][j] = sv2
    return dict(meta=dict(index=meta['index'], group=meta['group']), out=out, masks={s: masks[s] for s in ('plateau', 'transition', 'event')}, onset=onset, ev=ev)


def le(x, k):
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.log(np.where(x > 200., x, 200.)/k)


def stat(e):
    e = e[np.isfinite(e)]
    if len(e) < 20:
        return None
    q = np.percentile(e, [25, 75])
    return dict(n=len(e), med=float(np.median(np.abs(e))), p90=float(np.percentile(np.abs(e), 90)), bias=float(np.median(e)), rsd=float((q[1]-q[0])/1.349), rms=float(np.sqrt(np.mean(e*e))))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rep = {}
    sel = [('A', 'plateau'), ('A_slow_thin', 'plateau'), ('A', 'transition'), ('A_large', 'transition'), ('cough', 'event'), ('combo', 'event'), ('combo', 'transition')]
    print('estimator error ln(k/k_eff): |e| med/P90, rms | c constant, f free slope, o oracle slope scale, r V2 rule; noise-free c0 f0 o0')
    for g, s in sel:
        row = {}
        for est in ('c', 'f', 'o', 'r', 'c0', 'f0', 'o0'):
            e = np.concatenate([le(r['out'][est], r['out']['k'])[r['masks'][s]] for r in res if r['meta']['group'] == g])
            row[est] = stat(e)
        rep[f'{g}|{s}'] = row
        print(f'{g:12s}{s:11s}', ' '.join(f"{q} {row[q]['med']:.3f}/{row[q]['p90']:.3f} {row[q]['rms']:.3f}" for q in row if row[q]))
    print('within-case plateau robust sd: constant, free slope, ratio (P1)')
    for g in ('A', 'A_slow_thin', 'B', 'A+B', 'cough'):
        c, f = [], []
        for r in res:
            if r['meta']['group'] != g:
                continue
            m = r['masks']['plateau']
            for est, acc in (('c', c), ('f', f)):
                e = le(r['out'][est], r['out']['k'])[m]; e = e[np.isfinite(e)]
                if len(e) > 30:
                    q = np.percentile(e, [25, 75]); acc.append((q[1]-q[0])/1.349)
        if c:
            rep[f'{g}|plateau_sd'] = dict(constant=float(np.median(c)), free=float(np.median(f)), ratio=float(np.median(f)/np.median(c)))
            print(f'  {g:12s} {np.median(c):.3f} {np.median(f):.3f} ratio {np.median(f)/np.median(c):.2f}')
    print('trend statistic, relative to k (1/s): noise sd of v~/k, noise-free v~/k against the true slope; share of 6-tick checks with S_v^2 > 0 (V2 rule)')
    for g in H.GROUPS:
        rr = [r for r in res if r['meta']['group'] == g]
        svt = np.concatenate([r['out']['svt']/r['out']['k'] for r in rr])
        vt0 = np.concatenate([r['out']['vt0']/r['out']['k'] for r in rr]); vtr = np.concatenate([r['out']['vtrue']/r['out']['k'] for r in rr])
        on = np.concatenate([r['onset'] for r in rr]); pl = np.concatenate([r['masks']['plateau'] for r in rr])
        sv = np.concatenate([r['out']['sv2'] for r in rr])
        ok = np.isfinite(vt0) & np.isfinite(vtr)
        cc = float(np.corrcoef(vt0[ok], vtr[ok])[0, 1]) if ok.sum() > 30 else None
        gain = float(vt0[ok]@vtr[ok]/(vtr[ok]@vtr[ok])) if ok.sum() > 30 else None
        resid = vt0-vtr
        rep[f'{g}|trend'] = dict(noise_sd=float(np.nanmedian(svt)), corr=cc, slope=gain, true_p90=float(np.nanpercentile(np.abs(vtr), 90)),
                                 resid_rms=float(np.sqrt(np.nanmean(resid[ok]**2))), onset_resid=float(np.sqrt(np.nanmean(resid[ok & on]**2))) if (ok & on).sum() > 5 else None,
                                 active_all=float(np.nanmean(sv[np.isfinite(sv)] > 0)), active_plateau=float(np.nanmean(sv[pl & np.isfinite(sv)] > 0)) if (pl & np.isfinite(sv)).sum() else None)
        s = rep[f'{g}|trend']
        print(f"  {g:12s}", {q: (round(v, 3) if isinstance(v, float) else v) for q, v in s.items()})
    print('excitation (P6): amplitude factor needed for R_k <= (0.2 k)^2 with sd ~ 1/A and the model floor beta from noise-free error rms; share of windows needing <=1, <=2, infeasible')
    for g in H.GROUPS:
        rr = [r for r in res if r['meta']['group'] == g]
        sd1 = np.concatenate([r['out']['sd1'] for r in rr])
        e0 = np.concatenate([le(r['out']['c0'], r['out']['k'])[r['masks']['plateau']] for r in rr]); e0 = e0[np.isfinite(e0)]
        beta = float(np.sqrt(np.mean((e0-np.median(e0))**2))) if len(e0) > 20 else .1
        sd1 = sd1[np.isfinite(sd1)]
        need = np.where(beta < .2, sd1/np.sqrt(max(.04-beta*beta, 1e-9)), np.inf)
        rep[f'{g}|excitation'] = dict(beta=beta, sd1_med=float(np.median(sd1)), need_med=float(np.median(need)), le1=float(np.mean(need <= 1)), le2=float(np.mean(need <= 2)), infeasible=float(np.mean(need > 2)))
        s = rep[f'{g}|excitation']
        print(f"  {g:12s} beta {beta:.3f} sd at x1 {s['sd1_med']:.3f} | needed factor median {s['need_med']:.2f} | <=1 {s['le1']:.2f} <=2 {s['le2']:.2f} >2 {s['infeasible']:.2f}")
    (HERE/'results/v2_doc_check.json').write_text(json.dumps(rep, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
