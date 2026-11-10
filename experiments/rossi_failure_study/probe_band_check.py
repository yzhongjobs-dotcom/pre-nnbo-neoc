"""Can the known excitation frequency separate stiffness information from noise, body motion and mechanical oscillation, and does V2 then work?

Same ungated Huber-D4 closed loops as huber_error_mechanism.py (unchanged estimator in the loop); every estimator below sees the same windows.
Band-limited regression: P_B projects a 24-tick window onto cos/sin at the DFT bins of the band, after removing a cubic (the D4 null space),
so k = Y' P_B z / Y' P_B Y. Bands: bin 4 only (20.8 Hz, a lock-in), bins 3-5 (15.6-26.0 Hz, keeps the +-5.2 Hz sidebands that carry a slope).
V2 in the band: the same 2x2 system of the note with a, b, c1, p0, p1 computed through P_B; h1 is the exact time-weighted regressor.
Hold: when the band information gives a relative sd above 0.5, keep the last accepted value (the note's "information insufficient" status).
Predictions fixed before running: P1 no collapse in band; P2 normal-window sd x1.1-1.25; P3 fall over-estimates unchanged and V2 slope still
below noise; P4 fewer false trend alarms at cough onset.
Writes results/probe_band_check.json.
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
import huber_error_mechanism as E
import huber_prior_check as H
import v2_verify as V

HERE = Path(__file__).resolve().parent
NW, TC, SIG = 24, .008, H.SIGMA
CASES = (0, 7, 12, 24, 26, 29, 32, 35, 36, 38)


def projector(bins):
    t = np.arange(NW)
    P = np.vander(t-t.mean(), 4)
    cols = [f(2*np.pi*b*t/NW) for b in bins for f in (np.cos, np.sin)]
    X = np.column_stack(cols)
    Qp, _ = np.linalg.qr(P)
    X = X-Qp@(Qp.T@X)
    Qx, _ = np.linalg.qr(X)
    return Qx@Qx.T


PB = {'lock-in': projector((4,)), 'band': projector((3, 4, 5))}


def band_fit(P, h0, h1, z, sv2):
    a, b, c1 = h0@P@h0, h0@P@h1, h1@P@h1
    p0, p1 = h0@P@z, h1@P@z
    if a <= 0:
        return np.nan, np.inf
    if sv2 == 0:
        k = p0/a; var = SIG**2/a
    else:
        cp = c1+(SIG**2/sv2 if np.isfinite(sv2) else 0.)
        info = a-b*b/cp
        if info <= 0:
            return np.nan, np.inf
        q = P@h0-(b/cp)*(P@h1); q = q/info
        k = q@z; var = SIG**2*(q@q)
    return k, math.sqrt(max(var, 0.))


DM = np.diff(np.eye(NW), 4, axis=0)
MD4 = np.eye(NW)-np.linalg.inv(np.eye(NW)+H.RHO*DM.T@DM)


def iv_fit(h0, h1, z, sv2, P=None):
    """D4 removes the smooth disturbance; only the probe-band part of the filtered motion is used as weight (instrument).
    Constant: k = w0'z / w0'h0 with w0 = P M h0. V2: instruments w0, w1 = P M h1, normal-equation analogue with the slope prior."""
    P = PB['band'] if P is None else P
    w0, w1 = P@MD4@h0, P@MD4@h1
    if sv2 == 0:
        a = w0@h0
        if a <= 0:
            return np.nan, np.inf, np.nan, np.inf
        q = w0/a
        return q@z, SIG*math.sqrt(q@q), np.nan, np.inf
    lam = SIG**2/sv2 if np.isfinite(sv2) else 0.
    A = np.array([[w0@h0, w0@h1], [w1@h0, w1@h1+lam]])
    try:
        Ai = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return np.nan, np.inf, np.nan, np.inf
    Q = Ai@np.vstack([w0, w1])
    return Q[0]@z, SIG*math.sqrt(Q[0]@Q[0]), Q[1]@z, SIG*math.sqrt(Q[1]@Q[1])


def trend_band(P, h0, h1, z):
    a0, b0 = h0@P@h0, h0@P@h1
    g = h1-(b0/a0)*h0; Jv = g@P@g
    if Jv <= 0:
        return np.nan, np.nan
    lv = P@g/Jv
    return lv@z, SIG*math.sqrt(lv@lv)


def job(index):
    case = CC.cases()[index]
    sc, pol, log, rows, sent, tau_h, mag_h, eb, policy, raw = E.closed_loop(case)
    n = len(rows)
    Y, z, H1 = E.W.windows(E.W.hc_log(rows, tau_h, mag_h, policy['Tc']), policy, policy['enc'], eb)
    k = log[:n, 26]; t = log[:n, 0]
    keys = ['iv_const', 'iv_sd', 'iv_free', 'iv_oracle', 'iv_v', 'iv_vsd', 'huber', 'band_const', 'lock_const', 'band_free', 'band_oracle', 'band_sd', 'huber_sd', 'band_v', 'band_vsd', 'huber_v', 'huber_vsd', 'vtrue']
    out = {q: np.full(n, np.nan) for q in keys}
    S = H.Solver()
    for j in range(NW+1, n):
        h0, h1, zz = Y[j][-NW:], H1[j][-NW:], z[j][-NW:]
        if not (np.isfinite(h0).all() and np.isfinite(h1).all()) or k[j] <= 200.:
            continue
        r_ = np.arange(j-NW, j); kr = .5*(k[r_]+k[r_+1])
        tt = .5*(t[r_]+t[r_+1])-t[j-1]
        vtrue = float(np.polyfit(tt, kr, 1)[0]) if np.all(kr > 200.) else np.nan
        f = S.fit(h0, zz)
        out['huber'][j] = f['k']; out['huber_sd'][j] = SIG/math.sqrt(f['info']) if f['info'] > 0 else np.inf
        out['band_const'][j], out['band_sd'][j] = band_fit(PB['band'], h0, h1, zz, 0.)
        out['lock_const'][j] = band_fit(PB['lock-in'], h0, h1, zz, 0.)[0]
        out['band_free'][j] = band_fit(PB['band'], h0, h1, zz, np.inf)[0]
        if np.isfinite(vtrue):
            out['band_oracle'][j] = band_fit(PB['band'], h0, h1, zz, vtrue**2)[0]
        out['band_v'][j], out['band_vsd'][j] = trend_band(PB['band'], h0, h1, zz)
        M = np.eye(NW)-np.linalg.inv(np.eye(NW)+H.RHO*np.diff(np.eye(NW), 4, axis=0).T@np.diff(np.eye(NW), 4, axis=0))
        out['huber_v'][j], out['huber_vsd'][j] = trend_band(M, h0, h1, zz)
        out['vtrue'][j] = vtrue
        out['iv_const'][j], out['iv_sd'][j], _, _ = iv_fit(h0, h1, zz, 0.)
        out['iv_free'][j], _, out['iv_v'][j], out['iv_vsd'][j] = iv_fit(h0, h1, zz, np.inf)
        if np.isfinite(vtrue):
            out['iv_oracle'][j] = iv_fit(h0, h1, zz, vtrue**2)[0]
    # Hold when the band says the information is insufficient.
    for src, sd in (('band_const', 'band_sd'), ('huber', 'huber_sd'), ('iv_const', 'iv_sd')):
        held = np.full(n, np.nan); last = np.nan
        for j in range(n):
            v, s = out[src][j], out[sd][j]
            if np.isfinite(v) and v > 0 and s/max(v, 1.) <= .5:
                last = v
            held[j] = last if np.isfinite(out[src][j]) else np.nan
        out[src+'_hold'] = held
    ev = np.zeros(n, bool); onset = np.zeros(n, bool)
    for e in sc['event_info']:
        ev |= (t >= e['start']) & (t <= e['end']+.3); onset |= (t >= e['start']) & (t <= e['start']+.2)
    seg = V.segments(t, k, sc['event_info'], k > 200.)
    # Motion spectrum share in the window: what fraction of the window's band-passed motion energy is at the probe bins.
    return dict(index=index, group=case['group'], k=k, t=t, out=out, seg=seg, onset=onset, raw=raw)


def summarize(res):
    names = ('huber', 'huber_hold', 'band_const', 'band_const_hold', 'lock_const', 'band_free', 'band_oracle', 'iv_const', 'iv_const_hold', 'iv_free', 'iv_oracle')
    rep = {}
    print('per case: sustained 0.1 s ln error min / max, share of windows below .5 / above 1.5 of truth, plateau within-case sd')
    for r in res:
        line = {}
        for nm in names:
            v = r['out'][nm]
            ok = np.isfinite(v) & (r['k'] > 200.) & (r['t'] >= 1.)
            e = np.full(len(v), np.nan); e[ok] = np.log(np.maximum(v[ok], 200.)/r['k'][ok])
            sm = uniform_filter1d(np.nan_to_num(e), 13, mode='nearest'); sm[~ok] = np.nan
            pe = e[r['seg']['plateau'] & ok]
            q = np.percentile(pe, [25, 75]) if len(pe) > 30 else [np.nan, np.nan]
            line[nm] = dict(min=float(np.nanmin(sm)), max=float(np.nanmax(sm)), low=float(np.nanmean(e[ok] < math.log(.5))), high=float(np.nanmean(e[ok] > math.log(1.5))), sd=float((q[1]-q[0])/1.349))
        rep[r['index']] = line
        print(f"== {r['index']:2d} {r['group']:12s}", ' | '.join(f"{nm} {line[nm]['min']:+.2f}/{line[nm]['max']:+.2f} {line[nm]['low']*100:.1f}/{line[nm]['high']*100:.1f}% sd {line[nm]['sd']:.3f}" for nm in ('huber', 'band_const', 'iv_const', 'iv_const_hold', 'iv_free')))
    print('by segment over all cases, |e| med/P90 and rms: huber | band const | band hold | band free slope | band oracle slope')
    seg_rep = {}
    for s in V.SEGS:
        row = {}
        for nm in names:
            e = []
            for r in res:
                v = r['out'][nm]; m = r['seg'][s] & np.isfinite(v)
                e.append(np.log(np.maximum(v[m], 200.)/r['k'][m]))
            e = np.concatenate(e)
            row[nm] = None if len(e) < 20 else dict(n=len(e), med=float(np.median(np.abs(e))), p90=float(np.percentile(np.abs(e), 90)), rms=float(np.sqrt(np.mean(e*e))), bias=float(np.median(e)))
        seg_rep[s] = row
        print(f'  {s:9s}', ' | '.join(f"{nm} {row[nm]['med']:.3f}/{row[nm]['p90']:.3f} {row[nm]['rms']:.3f}" for nm in ('huber', 'band_const', 'iv_const', 'iv_const_hold', 'iv_free', 'iv_oracle') if row[nm]))
    print('trend statistic relative to k (1/s): noise sd and error at cough/event onset, D4 versus band')
    tr = {}
    for nm in ('huber', 'band', 'iv'):
        vs = np.concatenate([r['out'][f'{nm}_vsd']/r['k'] for r in res]); vv = np.concatenate([(r['out'][f'{nm}_v']-r['out']['vtrue'])/r['k'] for r in res])
        on = np.concatenate([r['onset'] for r in res]); ok = np.isfinite(vv)
        tr[nm] = dict(noise_sd=float(np.nanmedian(vs)), err_rms=float(np.sqrt(np.mean(vv[ok]**2))), onset_err_rms=float(np.sqrt(np.mean(vv[ok & on]**2))), onset_z_over_2=float(np.mean(np.abs(vv[ok & on]) > 2*np.nanmedian(vs))))
        print(f'  {nm:6s}', {q: round(v, 3) for q, v in tr[nm].items()})
    sd_ratio = np.nanmedian(np.concatenate([r['out']['band_sd']/r['out']['huber_sd'] for r in res]))
    iv_ratio = np.nanmedian(np.concatenate([r['out']['iv_sd']/r['out']['huber_sd'] for r in res]))
    print('band sd / Huber sd, median over windows (P2):', round(float(sd_ratio), 3), '; D4 + band instrument sd / Huber sd:', round(float(iv_ratio), 3))
    return dict(cases=rep, segments=seg_rep, trend=tr, sd_ratio=float(sd_ratio), iv_sd_ratio=float(iv_ratio))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(job, CASES, chunksize=1))
    rep = summarize(res)
    (HERE/'results/probe_band_check.json').write_text(json.dumps(rep, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
