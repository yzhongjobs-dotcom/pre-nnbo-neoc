"""How much information about the contact gain lies where, on recorded physical windows. No new simulation.

For each recorded window of results/structured_tvk_audit (z = kY + d + eps, 24 ticks, unchanged Huber-D4 trajectories) compare:
  D4 fc        quadratic D4 nuisance with cutoff fc: is d as narrow as breathing would allow?
  oracle u     regressor Y-U with the true exogenous motion subtracted, then the same nuisance
  oracle depth regressor from the true indentation (everything but stiffness known), constant offset only
Information is x'Mx relative to the fc=14 Hz value. Errors are ln(khat/k) against the true |K(20.8 Hz)| at the window end.
Truth enters scoring and the oracle regressors only; these are diagnostic ceilings, not estimators anyone could run.
Second part: at constant force a power-law contact has k = nF/depth; check ln k against ln depth on the 8 N overlay logs.
Usage: python -X utf8 identification_information_audit.py
"""
import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AUD = HERE/'results'/'structured_tvk_audit'
Tc = .008
afc = 1.-np.exp(-2.*np.pi*30.*.001); a = 1.-afc; NS = 8
W = np.array([afc*a**(NS-m) for m in range(1, NS+1)])
FCS = (14., 7., 3.5, 1.75)


def mat(fc, n=24):
    D = np.zeros((n-4, n))
    for i in range(n-4):
        D[i, i:i+5] = (1, -4, 6, -4, 1)
    rho = (2.*np.sin(np.pi*fc*Tc))**-8
    return np.eye(n)-np.linalg.inv(np.eye(n)+rho*D.T@D)


def const(n=24):
    return np.eye(n)-np.ones((n, n))/n


def filt_1k(u, t):
    """Sum_m w_m u(t + (m-1) ms), u sampled at 1 kHz."""
    i = int(round(t*1000.))
    return float(W@u[i:i+NS])


def h_terms(phi):
    """tau_h and magnitude of the phase-lead preprocessing that Y went through (contact_candidates)."""
    diff = (1.-np.exp(-2j*np.pi*20.8*Tc))/Tc
    tan = np.tan(phi); tau = tan/max(diff.imag-tan*diff.real, 1.)
    return tau, abs(1.+tau*diff)


def filt_8ms(v, j):
    fr = np.arange(NS)/NS
    return float(W@(v[j]+fr*(v[j+1]-v[j])))


def run():
    M = {fc: mat(fc) for fc in FCS}; C = const()
    out = {}
    for cond in ('rest', 'breath', 'shift_fast', 'rib', 'cough', 'combo'):
        acc = {}
        info = {}
        lowres = {}
        for path in sorted(glob.glob(str(AUD/f'{cond}_[0-9].npz'))):
            z_ = np.load(path, allow_pickle=True)
            log = z_['log']; u = z_['sc_u']
            tot = log[:, 2]-log[:, 6]
            for w, (idx, j0) in enumerate(zip(z_['idx'], z_['j0'])):
                Y, z = z_['Y'][w], z_['z'][w]
                k = log[idx, 21]
                if not k > 200.:
                    continue
                # the plant's depth base is y + body + h, and sc_u drives body positively, so in Y-U the motion column is -u;
                # U also goes through the same H preprocessing as Y (corrected after the 2026-11-08 Codex audit)
                tau, mag = h_terms(float(z_['correctedPhi'][w]))
                U = -np.array([filt_1k((u+tau*np.gradient(u, .001))/mag, log[j0+q, 0]) for q in range(24)]); U = U-U[-1]
                T = np.array([filt_8ms(tot, j0+q) for q in range(24)]); T = T-T[-1]
                base = Y@M[14.]@Y
                def rec(name, x, m):
                    den = x@m@x
                    if den <= 0:
                        return
                    acc.setdefault(name, []).append(np.log(max((x@m@z)/den, 1.)/k))
                    if name.startswith('oracle u') or name == 'D4 14 Hz':
                        lowres.setdefault(name, []).append(np.sqrt(np.mean(((z-k*x)-np.mean(z-k*x))**2)))
                    info.setdefault(name, []).append(den/base)
                rec('oracle u, const', Y-U, C)
                for fc in FCS:
                    rec(f'D4 {fc:g} Hz', Y, M[fc])
                    rec(f'oracle u, D4 {fc:g} Hz', Y-U, M[fc])
                rec('oracle depth, const', Y-T, C)
                rec('oracle depth, D4 14 Hz', Y-T, M[14.])
                acc.setdefault('Huber-D4 recorded', []).append(np.log(max(z_['baseline'][w][0], 1.)/k))
                info.setdefault('Huber-D4 recorded', []).append(1.)
        row = {}
        for name, e in acc.items():
            e = np.array(e)
            row[name] = dict(n=len(e), median_abs=float(np.median(np.abs(e))), bias=float(np.median(e)), iqr_sd=float((np.percentile(e, 75)-np.percentile(e, 25))/1.349), p90_abs=float(np.percentile(np.abs(e), 90)), info=float(np.median(info[name])), residual_rms=float(np.median(lowres[name])) if name in lowres else None)
        out[cond] = row
        print(cond)
        for name, r in row.items():
            print(f"  {name:26s} |e| med {r['median_abs']:.3f} bias {r['bias']:+.3f} sd {r['iqr_sd']:.3f} p90 {r['p90_abs']:.3f} info {r['info']:.2f}", flush=True)
    return out


def depth_law():
    res = {}
    for scene, pat in (('breath_rib', 'breath_rib_rib'), ('rib', 'rib_rib'), ('breath', 'breath_relaxed')):
        slopes, r2 = [], []
        for p in sorted(glob.glob(str(HERE/'results'/'baseline_logs'/f'Huber-D4_{pat}_*.npz'))):
            z = np.load(p); L = z['log']; c = list(z['cols']); g = lambda i: L[:, c.index(i)]
            t, f, k, d = g(0), g(18), g(21), g(6)
            w = (t >= 1.) & (f > 5.) & (f < 11.) & (k > 200.) & (d > 1e-4)
            if w.sum() < 100 or np.ptp(np.log(d[w])) < .02:
                continue
            # one-second smoothing keeps the slow relation
            n = 62; sm = lambda v: np.convolve(v, np.ones(n)/n, 'valid')
            x, y = sm(np.log(d[w])), sm(np.log(k[w]/f[w]))
            b = np.polyfit(x, y, 1); r = np.corrcoef(x, y)[0, 1]
            slopes.append(b[0]); r2.append(r*r)
        res[scene] = dict(n=len(slopes), slope_median=float(np.median(slopes)), slope_p10=float(np.percentile(slopes, 10)), slope_p90=float(np.percentile(slopes, 90)), r2_median=float(np.median(r2)))
        print('ln(k/F) vs ln depth', scene, res[scene], flush=True)
    return res


if __name__ == '__main__':
    r = dict(windows=run(), depth_law=depth_law())
    (HERE/'results'/'identification_information_audit.json').write_text(json.dumps(r, indent=1))
