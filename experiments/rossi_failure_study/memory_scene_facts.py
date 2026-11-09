"""Scene facts for the adaptive-memory study: how the true gain changes, what the window estimate's error looks like, and what other measured channel responds to a change.

Reads results/huber_windows/*.npz only; no plant, no estimator choice. Truth is read for analysis, never as an input to any candidate.
Writes results/memory_scene_facts.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import huber_prior_check as H

HERE = Path(__file__).resolve().parent
COL = {c: i for i, c in enumerate([0, 1, 2, 3, 5, 11, 14, 18, 21, 26, 27])}
TC = .008
GROUPS = H.GROUPS
RATE_ON = .5       # /s: "changing"
SIGMA = H.SIGMA
EPS_LOW = math.log(.85)


def monotone_runs(ell, t, contact):
    """Rise and fall runs of the 0.04 s smoothed ln k: size, duration, peak rate; and how each run ends."""
    s = uniform_filter1d(ell, 5, mode='nearest')
    r = np.gradient(s, t)
    runs = []; sign = 0; a = None
    for j in range(len(t)):
        cur = 0 if not contact[j] else (1 if r[j] > RATE_ON else -1 if r[j] < -RATE_ON else 0)
        if cur != sign:
            if sign != 0:
                runs.append((a, j, sign))
            a = j; sign = cur
    if sign != 0:
        runs.append((a, len(t), sign))
    # Merge runs of the same sign separated by under 0.04 s.
    merged = []
    for q in runs:
        if merged and merged[-1][2] == q[2] and t[q[0]]-t[merged[-1][1]-1] < .04:
            merged[-1] = (merged[-1][0], q[1], q[2])
        else:
            merged.append(q)
    out = []
    for i, (a, b, sg) in enumerate(merged):
        d = float(s[b-1]-s[a])
        if abs(d) < .05:
            continue
        gap = float(t[merged[i+1][0]]-t[b-1]) if i+1 < len(merged) else None
        flip = merged[i+1][2] == -sg if i+1 < len(merged) else None
        out.append(dict(start=float(t[a]), dur=float(t[b-1]-t[a]), delta=d, peak_rate=float(np.max(np.abs(r[a:b]))), gap_next=gap, next_opposite=flip))
    return out, r


def case(path):
    data = H.load(path); meta = data['meta']; log = data['log']
    t = log[:, COL[0]]; x = log[:, COL[1]]; k = log[:, COL[26]]; F = log[:, COL[5]]; y = log[:, COL[2]]
    contact = k > 200.
    ell = np.log(np.maximum(k, 1.))
    runs, rate = monotone_runs(ell, t, contact)
    speed = np.abs(np.gradient(x, t))
    S = H.Solver()
    n = len(t)
    keys = ('lnD', 'lnD0', 'pred', 'force', 'k_d', 'noise', 'a_f0', 'a_band', 'res_rate')
    w = {q: np.full(n, np.nan) for q in keys}
    prev_d = None
    for j in range(n):
        Y, z = H.tail(data['Y'][j], 24), H.tail(data['z'][j], 24)
        Y0, z0 = H.tail(data['Y0'][j], 24), H.tail(data['z0'][j], 24)
        if len(Y) < 24 or len(Y0) != 24 or k[j] <= 200.:
            prev_d = None
            continue
        a = S.fit(Y, z, iters=1); b = S.fit(Y0, z0, iters=1)
        if not (a['k'] > 0 and b['k'] > 0):
            prev_d = None
            continue
        w['lnD'][j] = math.log(H.clip(a['k'])); w['lnD0'][j] = math.log(H.clip(b['k']))
        w['pred'][j] = SIGMA/(a['k']*math.sqrt(a['info'])) if a['info'] > 0 else np.nan
        my = a['my']; tt = np.arange(24)
        w['a_band'][j] = k[j]*math.sqrt(2.)*float(np.std(my))
        w['a_f0'][j] = k[j]*2*abs(np.mean(my*np.exp(-2j*np.pi*H.F0*TC*tt)))
    # Force channel: slow force change not explained by motion at the current gain, in N/s over 0.1 s. It contains k du/dt and depth dk/dt.
    L = 12
    good = np.isfinite(w['lnD'])
    kh = np.exp(np.where(good, w['lnD'], np.nan))
    kh = uniform_filter1d(np.nan_to_num(kh, nan=np.nanmedian(kh)), 12, mode='nearest')
    res = np.full(n, np.nan)
    res[L:] = ((F[L:]-F[:-L])-kh[L:]*(y[L:]-y[:-L]))/(L*TC)
    # Remove the trend a slow motion makes: rate now minus rate 0.2 s earlier.
    hp = np.full(n, np.nan)
    hp[25:] = res[25:]-res[:-25]
    e_n = w['lnD']-w['lnD0']; e_0 = w['lnD0']-ell; e = w['lnD']-ell
    return dict(meta={q: meta[q] for q in ('index', 'group', 'event_info')}, t=t, ell=ell, rate=rate, speed=speed, x=x, F=F, contact=contact, runs=runs, e=e, e_n=e_n, e_0=e_0, pred=w['pred'], a_f0=w['a_f0'], a_band=w['a_band'], res=res, hp=hp, lnD=w['lnD'])


def acf(x, lags):
    out = []
    for q in lags:
        vals = []
        for s in x:
            a, b = s[:-q] if q else s, s[q:]
            m = np.isfinite(a) & np.isfinite(b)
            if m.sum() > 100:
                vals.append(np.corrcoef(a[m], b[m])[0, 1])
        out.append(float(np.mean(vals)))
    return out


def auc(pos, neg):
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not len(pos) or not len(neg):
        return None
    r = np.argsort(np.argsort(np.r_[pos, neg]))+1
    return float((r[:len(pos)].sum()-len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg)))


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    report = {}
    print('--- true gain: runs with |d ln k| >= .05 per scene')
    print('group runs/case | |delta| med P90 max | dur med P90 s | peak rate med P90 /s | changing time share | runs followed by opposite within .1 s')
    runs_out = {}
    for g in GROUPS:
        rr = [r for r in res if r['meta']['group'] == g]
        all_runs = [q for r in rr for q in r['runs']]
        d = np.abs([q['delta'] for q in all_runs]); du = np.array([q['dur'] for q in all_runs]); pk = np.array([q['peak_rate'] for q in all_runs])
        gaps = [q for q in all_runs if q['gap_next'] is not None]
        quick = np.mean([q['next_opposite'] and q['gap_next'] < .1 for q in gaps]) if gaps else float('nan')
        share = np.mean([np.mean(np.abs(r['rate'][r['contact']]) > RATE_ON) for r in rr])
        runs_out[g] = dict(n=len(all_runs), delta_med=float(np.median(d)), delta_p90=float(np.percentile(d, 90)), delta_max=float(d.max()), dur_med=float(np.median(du)), dur_p90=float(np.percentile(du, 90)), rate_med=float(np.median(pk)), rate_p90=float(np.percentile(pk, 90)), share=float(share), quick_reversal=float(quick))
        s = runs_out[g]
        print(f"{g:12s} {len(all_runs)/len(rr):5.1f} | {s['delta_med']:.3f} {s['delta_p90']:.3f} {s['delta_max']:.3f} | {s['dur_med']:.3f} {s['dur_p90']:.3f} | {s['rate_med']:.2f} {s['rate_p90']:.2f} | {s['share']:.2f} | {s['quick_reversal']:.2f}")
    report['runs'] = runs_out
    print('--- position versus time: |d ln k/dt| when the probe moves and when it is nearly still (|v|<1 mm/s), no event, after 1 s')
    pos = {}
    for g in GROUPS:
        mv, st, per_mm = [], [], []
        for r in res:
            if r['meta']['group'] != g:
                continue
            ok = r['contact'] & (r['t'] > 1.)
            for e in r['meta']['event_info']:
                ok &= ~((r['t'] >= e['start']) & (r['t'] <= e['end']+.3))
            mv.append(np.abs(r['rate'][ok & (r['speed'] > .001)])); st.append(np.abs(r['rate'][ok & (r['speed'] <= .001)]))
        mv, st = np.concatenate(mv), np.concatenate(st)
        pos[g] = dict(moving=int(len(mv)), still=int(len(st)), moving_p90=float(np.percentile(mv, 90)) if len(mv) else None, still_p90=float(np.percentile(st, 90)) if len(st) else None, moving_share_on=float(np.mean(mv > RATE_ON)) if len(mv) else None, still_share_on=float(np.mean(st > RATE_ON)) if len(st) else None)
        print(g, pos[g])
    report['position'] = pos
    print('--- window error: sensor-noise part and noise-free part; amplitude actually in the window')
    print('group | rms e_n  z rms  z |>3| | e_0 med abs, P90, low<.85 | probe force a_f0 med N, a_band med N | per-tick info I, implied sd(24)')
    noise = {}
    for g in GROUPS+['all']:
        rr = [r for r in res if g == 'all' or r['meta']['group'] == g]
        en = np.concatenate([r['e_n'] for r in rr]); pr = np.concatenate([r['pred'] for r in rr]); e0 = np.concatenate([r['e_0'] for r in rr]); af = np.concatenate([r['a_f0'] for r in rr]); ab = np.concatenate([r['a_band'] for r in rr])
        m = np.isfinite(en) & np.isfinite(pr)
        zz = en[m]/pr[m]
        am = float(np.nanmedian(af)); info = (am/SIGMA)**2/2.*.77
        noise[g] = dict(e_n_rms=float(np.sqrt(np.mean(en[m]**2))), z_rms=float(np.sqrt(np.mean(zz**2))), z_robust=float((np.percentile(zz, 75)-np.percentile(zz, 25))/1.349), z_tail3=float(np.mean(np.abs(zz) > 3)), skew_low=float(np.percentile(zz, 1)), skew_high=float(np.percentile(zz, 99)), e0_med=float(np.nanmedian(np.abs(e0))), e0_p90=float(np.nanpercentile(np.abs(e0), 90)), e0_low=float(np.nanmean(e0[np.isfinite(e0)] < EPS_LOW)), a_f0=am, a_band=float(np.nanmedian(ab)), tick_info=info, sd24=float(1./math.sqrt(info*24.)))
        s = noise[g]
        print(f"{g:12s} | {s['e_n_rms']:.3f} {s['z_rms']:.2f} rob {s['z_robust']:.2f} tail {s['z_tail3']:.3f} p1 {s['skew_low']:.1f} p99 {s['skew_high']:.1f} | {s['e0_med']:.3f} {s['e0_p90']:.3f} {s['e0_low']:.3f} | {s['a_f0']:.3f} {s['a_band']:.3f} | {s['tick_info']:.2f} {s['sd24']:.3f}")
    report['noise'] = noise
    lags = [0, 1, 2, 3, 6, 12, 18, 24, 30]
    report['acf_noise'] = dict(lags=lags, values=acf([r['e_n'] for r in res], lags))
    print('--- autocorrelation of the noise part by lag in ticks', dict(zip(lags, np.round(report['acf_noise']['values'], 3))))
    print('--- what signals a true change (|d ln k/dt|>.5/s) against no change (<.1/s), no event, after 1 s; AUC of each candidate signal')
    sig = {}
    for g in GROUPS+['all']:
        P = {q: [] for q in ('res', 'hp', 'innov')}; N = {q: [] for q in P}
        for r in res:
            if g != 'all' and r['meta']['group'] != g:
                continue
            ok = r['contact'] & (r['t'] > 1.)
            for e in r['meta']['event_info']:
                ok &= ~((r['t'] >= e['start']) & (r['t'] <= e['end']+.3))
            on = ok & (np.abs(r['rate']) > RATE_ON); off = ok & (np.abs(r['rate']) < .1)
            # Innovation of the k channel: this window against the mean of the windows 24 to 48 ticks ago, which share no data.
            old = np.full(len(r['t']), np.nan); old[48:] = np.nanmean(np.lib.stride_tricks.sliding_window_view(r['lnD'], 25)[:len(r['t'])-48], axis=1)
            innov = np.abs(r['lnD']-old)
            for q, v in (('res', np.abs(r['res'])), ('hp', np.abs(r['hp'])), ('innov', innov)):
                P[q].append(v[on]); N[q].append(v[off])
        sig[g] = {q: auc(np.concatenate(P[q]), np.concatenate(N[q])) for q in P}
        sig[g]['on_ticks'] = int(sum(len(v) for v in P['res']))
        print(g, {q: (round(v, 3) if isinstance(v, float) else v) for q, v in sig[g].items()})
    report['signal_auc'] = sig
    (HERE/'results/memory_scene_facts.json').write_text(json.dumps(report, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
