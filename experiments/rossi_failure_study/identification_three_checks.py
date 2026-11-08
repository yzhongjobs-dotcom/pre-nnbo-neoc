"""Offline numerical checks of the derivations for the three identification questions. No closed-loop runs.

A. Passive contact: for K(s)=1+sum a_i s tau_i/(1+s tau_i) with a_i>0 (a complete Bernstein function), Schwarz-Pick on the strip
   gives |d ln|K|/d ln w + i d phi/d ln w| <= sin phi, and K(i w) lies in the hyperbolic disc of radius |ln(w/w0)| about K(i w0) in the upper half plane.
   Check on the plant's dynamic responses.
B. Material coordinate: at fixed material a power law gives k proportional to F^gamma, gamma=(n-1)/n. Replace Y by Y*(F/F_end)^gamma in the window.
C. Memory indexed by distance: log-domain Kalman on non-overlapping window estimates with process noise q_t + q_s*|v|, Huber innovation.
   Parameters chosen on replica 0, scored on replicas 1 and 2. Errors are after removing each episode's median offset, since the logged truth has a known offset.
Usage: python -X utf8 identification_three_checks.py
"""
import glob
import json
from pathlib import Path

import numpy as np

import sweeps  # noqa: F401  (installs shared import paths)
import scenes3 as S3

HERE = Path(__file__).resolve().parent
AUD = HERE/'results'/'structured_tvk_audit'
Tc = .008


def part_a():
    f = np.logspace(np.log10(.3), np.log10(60.), 400); w = 2*np.pi*f; f0 = 20.8
    rows = []
    rng = np.random.default_rng(0)
    for kind in ('neural_white', 'neural_gray', 'viscous'):
        for ph in (8., 14., 21., 28.):
            try:
                alpha, tau, _ = S3._dynamic_response(kind, ph, rng)
            except Exception as e:
                print('skip', kind, ph, e); continue
            alpha, tau = np.asarray(alpha), np.asarray(tau)
            K = 1.+np.sum(alpha[None, :]*(1j*w[:, None]*tau[None, :])/(1.+1j*w[:, None]*tau[None, :]), axis=1)
            lnK, phi = np.log(np.abs(K)), np.unwrap(np.angle(K))
            lw = np.log(w)
            g, p = np.gradient(lnK, lw), np.gradient(phi, lw)
            ratio = np.sqrt(g**2+p**2)/np.sin(phi)
            K0 = 1.+np.sum(alpha*(1j*2*np.pi*f0*tau)/(1.+1j*2*np.pi*f0*tau)); A = abs(K0); ph0 = np.angle(K0)
            row = dict(kind=kind, phase=float(np.degrees(ph0)), slope_ratio_max=float(ratio.max()), monotone=bool(np.all(g > -1e-9)), all_positive=bool(np.all(alpha >= 0)))
            for ft in (5., 10., 35.):
                R = abs(np.log(ft/f0)); x0, y0 = K0.real, K0.imag
                c = complex(x0, y0*np.cosh(R)); r = y0*np.sinh(R)
                true = 1.+np.sum(alpha*(1j*2*np.pi*ft*tau)/(1.+1j*2*np.pi*ft*tau))
                inside = abs(true-c) <= r+1e-9
                row[f'{ft:g}Hz'] = dict(true=float(abs(true)/A), low=float((abs(c)-r)/A), high=float((abs(c)+r)/A), bode=float((ft/f0)**(2*ph0/np.pi)), inside=bool(inside))
            rows.append(row)
            print(kind, f"phi0 {row['phase']:.1f} slope/sin max {row['slope_ratio_max']:.3f} monotone {row['monotone']}",
                  {k: {a: round(b, 3) if isinstance(b, float) else b for a, b in row[k].items()} for k in ('5Hz', '35Hz')}, flush=True)
    return rows


def mat(fc=14., n=24):
    D = np.zeros((n-4, n))
    for i in range(n-4):
        D[i, i:i+5] = (1, -4, 6, -4, 1)
    rho = (2.*np.sin(np.pi*fc*Tc))**-8
    return np.eye(n)-np.linalg.inv(np.eye(n)+rho*D.T@D)


def part_b(gamma=.206):
    M = mat(); res = {}
    for cond in ('rest', 'breath', 'shift_fast', 'rib', 'cough', 'combo'):
        e0, e1, sl = [], [], []
        for path in sorted(glob.glob(str(AUD/f'{cond}_[0-9].npz'))):
            z_ = np.load(path, allow_pickle=True); log = z_['log']
            F, K = log[:, 18], log[:, 21]
            w = (log[:, 0] > 1.) & (F > .5) & (K > 200.)
            if w.sum() > 50 and np.ptp(np.log(F[w])) > .05:
                sl.append(np.polyfit(np.log(F[w]), np.log(K[w]), 1)[0])
            a, b = [], []
            for q, (idx, j0) in enumerate(zip(z_['idx'], z_['j0'])):
                Y, z = z_['Y'][q], z_['z'][q]; k = K[idx]
                if not k > 200. or log[j0:j0+24, 5].min() < .5:
                    continue
                Fw = log[j0:j0+24, 5]; Y2 = Y*(np.maximum(Fw, .3)/max(Fw[-1], .3))**gamma
                a.append(np.log(max((Y@M@z)/(Y@M@Y), 1.)/k)); b.append(np.log(max((Y2@M@z)/(Y2@M@Y2), 1.)/k))
            a, b = np.array(a), np.array(b)
            e0.append(a-np.median(a)); e1.append(b-np.median(b))
        e0, e1 = np.concatenate(e0), np.concatenate(e1)
        res[cond] = dict(base_med=float(np.median(np.abs(e0))), mat_med=float(np.median(np.abs(e1))), base_p90=float(np.percentile(np.abs(e0), 90)), mat_p90=float(np.percentile(np.abs(e1), 90)), slope_lnk_lnF=float(np.median(sl)) if sl else None)
        print('B', cond, {k: round(v, 3) if v is not None else None for k, v in res[cond].items()}, flush=True)
    return res


def kalman(l, R, dv, qt, qs, huber=2.):
    m, P = l[0], R[0]; out = [m]
    for i in range(1, len(l)):
        P = P+qt*.192+qs*dv[i]
        S = P+R[i]; nu = l[i]-m
        r = R[i]
        if huber and abs(nu)/np.sqrt(S) > huber:
            r = (abs(nu)/huber)**2-P
            r = max(r, R[i]); S = P+r
        G = P/S; m = m+G*nu; P = (1.-G)*P
        out.append(m)
    return np.array(out)


def load_series(cond, rep):
    z_ = np.load(AUD/f'{cond}_{rep}.npz', allow_pickle=True); log = z_['log']
    idx = z_['idx']; base = z_['baseline']
    keep = [q for q in range(0, len(idx), 6) if base[q][0] > 200. and log[idx[q], 21] > 200.]
    idx = idx[keep]; base = base[keep]
    x = log[idx, 1]
    ds = np.abs(np.diff(x, prepend=x[0]))*1e3   # mm between windows
    return np.log(base[:, 0]), (base[:, 1]/base[:, 0])**2, ds, np.log(log[idx, 21])


def part_c():
    conds = ('rest', 'breath', 'shift_fast', 'rib', 'cough', 'combo')
    def score(cfg, reps):
        out = {}
        for cond in conds:
            E = []
            for rep in reps:
                l, R, ds, truth = load_series(cond, rep)
                if cfg is None:
                    m = l
                else:
                    m = kalman(l, R, ds, cfg[0], cfg[1])
                e = m-truth; E.append(e-np.median(e))
            E = np.concatenate(E); out[cond] = (float(np.sqrt(np.mean(E**2))), float(np.percentile(np.abs(E), 90)))
        return out
    grid_t = (0., .002, .01, .05, .2)
    grid_s = (0., 2e-4, 1e-3, 5e-3, 2e-2)
    def obj(cfg):
        return np.mean([v[0] for v in score(cfg, (0,)).values()])
    fixed = min(((t, 0.) for t in grid_t if t > 0), key=obj)
    both = min(((t, s) for t in grid_t for s in grid_s if t+s > 0), key=obj)
    res = dict(fixed=fixed, both=both, window=score(None, (1, 2)), kalman_fixed=score(fixed, (1, 2)), kalman_distance=score(both, (1, 2)))
    print('C chosen fixed', fixed, 'distance', both)
    for c in conds:
        print('C', c, 'window rmse/p90', np.round(res['window'][c], 3), 'fixed q', np.round(res['kalman_fixed'][c], 3), 'distance q', np.round(res['kalman_distance'][c], 3), flush=True)
    return res


if __name__ == '__main__':
    r = dict(a=part_a(), b=part_b(), c=part_c())
    (HERE/'results'/'identification_three_checks.json').write_text(json.dumps(r, indent=1))
