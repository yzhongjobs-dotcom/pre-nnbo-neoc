"""How non-stiffness content enters ungated Huber-D4 and becomes an under- or over-estimate.

For a window with frozen Huber weights the estimate is exactly k_hat = q' z, q = M_W h0 / (h0' M_W h0), q' h0 = 1, and q' kills cubic trends.
So k_hat - k_T = q' r, r = z - k_T h0. Closed loop: Huber-D4 ungated (CC.own_mechanism), unchanged. The same command sequence is replayed
through the plant with one component removed at a time, cumulatively:
  R1 sensor noise off; R2 + true calibration in the sensing stage; R3 + scene-event motion removed; R4 + muscle activation removed;
  R5 + breathing and cardiac motion removed (body held at its initial position).
Each replay is scored against its own true gain, because removing motion or activation moves the operating depth and so the gain itself:
r_s = z_s - k_T^(s) h0_s. Contribution of step s: q' (r_{s-1} - r_s), all with the actual window's q, i.e. how much the estimation error changes
when that component is removed. The within-window stiffness change is removed from every replay first (scored on the actual run), then the steps above; what remains,
q' r_5, is split into the stiffness change inside
the window (cumulative sum of (k_i - k_T) dY_i) and the rest: the plant's own frequency-dependent response and heterogeneity. Each part is
also split by frequency band through the window DFT (q'x = sum over bins). Telescoping is exact by construction; its check is that q' r_0
reproduces the estimator's raw output.
Writes results/huber_error_mechanism.json.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import corrected_comparison as CC
import huber_prior_check as H
import huber_windows as W

P3, B, O, SE = W.P3, W.B, W.O, W.SE
HERE = Path(__file__).resolve().parent
NW = 24
CASES = (0, 7, 24, 26, 29, 36, 32, 35)
STEPS = ('noise', 'calibration', 'event_motion', 'activation', 'breath_cardiac')
BANDS = dict(low=(0, 2), probe=(3, 5), high=(6, 12))      # DFT bins of 5.21 Hz: <=10.4, 15.6-26.0, >=31.2 Hz


def closed_loop(case):
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    CC.own_mechanism(pol)
    est = pol.estimator; inner = est.estimate
    sent, seen, ebs, tau, mag, raw, policy = [], [], {}, {}, {}, {}, {}

    def record(o, n, p):
        policy.update(ndel=p.ndel, ns=p.ns, eb=p.eb, afc=p.afc, enc=bool(p.params['enc_y']), Tc=p.Tc)
        ebs[n-1] = float(p.eb)
        out = inner(o, n, p); raw[n-1] = out[0]
        return out
    est.estimate = record
    inner_log = est._log

    def record_log(o, n):
        own = est.owner
        for j in range(est.count, n):
            seen.append(o.log[j].copy()); tau[j] = float(own.tau_h); mag[j] = float(own.mag_h)
        return inner_log(o, n)
    est._log = record_log

    class Rec:
        def reset(self, platform):
            pol.reset(platform)

        def __call__(self, o):
            u = pol(o); sent.append(float(u)); return u
    log = P3.episode(sc, Rec(), case['seed'])
    rows = np.asarray(seen); n = len(rows)
    tau_h = np.nan_to_num(np.array([tau.get(j, np.nan) for j in range(n)]), nan=0.)
    mag_h = np.array([mag.get(j, np.nan) for j in range(n)]); mag_h = np.where(np.isfinite(mag_h), mag_h, 1.)
    eb = np.array([ebs.get(j, np.nan) for j in range(n)]); eb = np.where(np.isfinite(eb), eb, policy['eb'])
    return sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, np.array([raw.get(j, np.nan) for j in range(n)])


def replay(sc_mod, pol, sent, seed, n, tau_h, mag_h, eb, policy):
    rp = W.Replay(sent, pol, sc_mod)
    log = P3.episode(sc_mod, W.ReplayPolicy(rp), seed)
    rows = np.asarray(rp.seen)[:n]
    Y, z, _ = W.windows(W.hc_log(rows, tau_h[:n], mag_h[:n], policy['Tc']), policy, policy['enc'], eb)
    return Y, z, log[:n, 26]


def band_split(q, x):
    Q, X = np.fft.rfft(q), np.fft.rfft(x)
    part = np.real(np.conj(Q)*X)*2/len(q)
    part[0] /= 2
    if len(q) % 2 == 0:
        part[-1] /= 2
    return {b: float(part[lo:hi+1].sum()) for b, (lo, hi) in BANDS.items()}


def job(index):
    case = CC.cases()[index]
    sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, raw = closed_loop(case)
    n = len(rows)
    Y, z, _ = W.windows(W.hc_log(rows, tau_h, mag_h, policy['Tc']), policy, policy['enc'], eb)
    mods = []
    s = dict(sc, noise_scale=0.); mods.append(dict(s))
    CC.perfect_calibration(s); mods.append(dict(s))
    comp = np.sum(sc['motion_components'], axis=0) if len(sc['motion_components']) else np.zeros_like(sc['u'])
    s['u'] = sc['u']-comp; mods.append(dict(s))
    s['act'] = np.zeros_like(sc['act']); mods.append(dict(s))
    s['u'] = np.full_like(sc['u'], sc['u'][0]); mods.append(dict(s))
    reps = [replay(m, pol, sent, case['seed'], n, tau_h, mag_h, eb, policy) for m in mods]
    k = log[:n, 26]; t = log[:n, 0]
    S = H.Solver()
    keys = ['err', 'check'] + list(STEPS) + ['k_change', 'plant_rest']
    rec = {q: np.full(n, np.nan) for q in keys}
    bands = {q: np.full((n, 3), np.nan) for q in list(STEPS)+['k_change', 'plant_rest']}
    huber_active = np.full(n, np.nan)
    for j in range(NW+1, n):
        y0 = Y[j][-NW:]
        if not np.isfinite(y0).all() or k[j] <= 200. or not np.isfinite(raw[j]):
            continue
        Ys = [y0]+[r[0][j][-NW:] for r in reps]; zs = [z[j][-NW:]]+[r[1][j][-NW:] for r in reps]
        if not all(np.isfinite(a).all() for a in Ys+zs):
            continue
        f = S.fit(y0, zs[0])
        q = f['my']/f['info']
        kT = k[j]
        ks = [k]+[r[2] for r in reps]
        if any(kk[j] <= 200. for kk in ks):
            continue
        # Remove each replay's own within-window stiffness change first, so that a motion or activation step is scored
        # only for what it adds beyond changing the gain itself.
        r_ = np.arange(j-NW, j)
        kparts = []
        for kk, yy in zip(ks, Ys):
            kr = .5*(kk[r_]+kk[r_+1])
            kparts.append(np.r_[0., np.cumsum((kr[1:]-kk[j])*np.diff(yy))])
        rs = [zz-kk[j]*yy-kp for zz, yy, kk, kp in zip(zs, Ys, ks, kparts)]
        rec['err'][j] = (raw[j]-kT)/kT
        rec['check'][j] = (q@(rs[0]+kparts[0]))/kT
        rec['k_change'][j] = (q@kparts[0])/kT
        bands['k_change'][j] = [v/kT for v in band_split(q, kparts[0]).values()]
        for i, name in enumerate(STEPS):
            x = rs[i]-rs[i+1]
            rec[name][j] = (q@x)/kT
            bands[name][j] = [v/kT for v in band_split(q, x).values()]
        rec['plant_rest'][j] = (q@rs[-1])/kT
        bands['plant_rest'][j] = [v/kT for v in band_split(q, rs[-1]).values()]
        huber_active[j] = float(np.mean(f['w'] < 1.))
    # Failure moments: worst 0.1 s mean of the raw relative error, low and high.
    e = rec['err'].copy()
    sm = uniform_filter1d(np.nan_to_num(e), 13, mode='nearest'); sm[~np.isfinite(e) | (t < 1.)] = np.nan
    # How far each replay's own gain moved from the actual one, relative.
    shift = {name: float(np.nanmedian(np.abs(np.log(np.maximum(r[2], 1.)/np.maximum(k, 1.))[t >= 1.]))) for name, r in zip(STEPS, reps)}
    out = dict(index=index, group=case['group'], events=[(ev['kind'], ev['start'], ev['end']) for ev in sc['event_info']], gain_shift=shift,
               check_max=float(np.nanmax(np.abs(rec['check']-np.clip(rec['err'], -1, 10)))) if np.isfinite(rec['check']).any() else None, moments={})
    for name, idx in (('worst_low', int(np.nanargmin(sm))), ('worst_high', int(np.nanargmax(sm)))):
        span = (np.abs(np.arange(n)-idx) <= 6) & np.isfinite(rec['err'])
        m = {q: float(np.nanmean(rec[q][span])) for q in keys}
        m['t'] = float(t[idx]); m['k_eff'] = float(k[idx]); m['huber_active'] = float(np.nanmean(huber_active[span]))
        m['bands'] = {q: dict(zip(BANDS, np.nanmean(bands[q][span], axis=0).tolist())) for q in bands}
        out['moments'][name] = m
    # Whole-episode rms of each part, after 1 s.
    ok = np.isfinite(rec['err']) & (t >= 1.)
    out['rms'] = {q: float(np.sqrt(np.nanmean(rec[q][ok]**2))) for q in keys}
    return out


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(job, CASES, chunksize=1))
    print('telescoping check (q r_0 against the estimator output, relative), max:', [r['check_max'] for r in res])
    parts = list(STEPS)+['k_change', 'plant_rest']
    print('relative error and its parts at the worst 0.1 s low and high moment of each case (q r / k_eff); then band split of the two largest parts')
    for r in res:
        print(f"== case {r['index']} {r['group']} events {[(a, round(b, 2)) for a, b, _ in r['events']]} | replay gain shift |ln| median", {a: round(b, 3) for a, b in r['gain_shift'].items()})
        for name, m in r['moments'].items():
            print(f"  {name:10s} t {m['t']:.2f} k {m['k_eff']:.0f} err {m['err']:+.3f} huber rows {m['huber_active']:.2f} | "+' '.join(f"{p} {m[p]:+.3f}" for p in parts))
            top = sorted(parts, key=lambda p: -abs(m[p]))[:2]
            for p in top:
                print(f"      {p:14s} bands", {b: round(v, 3) for b, v in m['bands'][p].items()})
        print('  rms over the episode', {p: round(v, 3) for p, v in r['rms'].items()})
    (HERE/'results/huber_error_mechanism.json').write_text(json.dumps(res, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
