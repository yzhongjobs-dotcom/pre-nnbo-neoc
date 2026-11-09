"""Second look at the scene: separate the slow change of the true gain from its ripple, find where severe under-estimates come from, draw real curves.

Reads results/huber_windows/*.npz and the replay ticks of huber_unified_compare.py (argv[1]). Truth is for analysis only.
Writes results/memory_scene_curves.json and research/figures/memory_scene/curves.png.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d

import huber_prior_check as H
import memory_scene_facts as F

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COL = F.COL
TC = .008
TICKS = Path(sys.argv[1])/'24'
BOX = 25          # 0.2 s, the boxcar of the table's transition detector


def runs_smooth(ell, t, contact, on=.2):
    s = uniform_filter1d(ell, BOX, mode='nearest'); r = np.gradient(s, t)
    sign = 0; a = 0; runs = []
    for j in range(len(t)):
        c = 0 if not contact[j] else (1 if r[j] > on else -1 if r[j] < -on else 0)
        if c != sign:
            if sign:
                runs.append((a, j, sign))
            a, sign = j, c
    if sign:
        runs.append((a, len(t), sign))
    merged = []
    for q in runs:
        if merged and merged[-1][2] == q[2] and t[q[0]]-t[merged[-1][1]-1] < .1:
            merged[-1] = (merged[-1][0], q[1], q[2])
        else:
            merged.append(q)
    return [q for q in merged if abs(s[q[1]-1]-s[q[0]]) >= .15], s, r


def main():
    files = sorted(H.WIN.glob('case_*.npz'))
    res = [F.case(p) for p in files]
    rep = {}
    print('--- ripple: ln k_eff minus its 0.2 s mean; amplitude at the probe frequency; slow change runs at 0.2 s scale')
    print('group | ripple rms | at 20.8 Hz amp | runs/case | |delta| med P90 | dur med P90 | peak rate med P90 | dwell between rise and fall med min | still flat share')
    out = {}
    for g in H.GROUPS:
        rr = [r for r in res if r['meta']['group'] == g]
        rip, amp, deltas, durs, rates, dwell, flat = [], [], [], [], [], [], []
        for r in rr:
            ell, t, c = r['ell'], r['t'], r['contact']
            runs, s, rate = runs_smooth(ell, t, c)
            d = (ell-s)[c & (t > 1.)]
            rip.append(float(np.std(d)))
            tt = t[c & (t > 1.)]
            amp.append(float(2*abs(np.mean(d*np.exp(-2j*np.pi*20.8*tt)))))
            flat.append(float(np.mean(np.abs(rate[c & (t > 1.)]) < .05)))
            for i, (a, b, sg) in enumerate(runs):
                deltas.append(abs(s[b-1]-s[a])); durs.append(t[b-1]-t[a]); rates.append(np.max(np.abs(rate[a:b])))
                if i+1 < len(runs) and runs[i+1][2] == -sg:
                    dwell.append(t[runs[i+1][0]]-t[b-1])
        q = lambda v, p: float(np.percentile(v, p)) if len(v) else float('nan')
        out[g] = dict(ripple=float(np.mean(rip)), amp208=float(np.mean(amp)), runs=len(deltas)/len(rr), delta_med=q(deltas, 50), delta_p90=q(deltas, 90), dur_med=q(durs, 50), dur_p90=q(durs, 90), rate_med=q(rates, 50), rate_p90=q(rates, 90), dwell_med=q(dwell, 50), dwell_min=float(min(dwell)) if dwell else float('nan'), flat=float(np.mean(flat)))
        s = out[g]
        print(f"{g:12s} | {s['ripple']:.3f} | {s['amp208']:.3f} | {s['runs']:.1f} | {s['delta_med']:.3f} {s['delta_p90']:.3f} | {s['dur_med']:.2f} {s['dur_p90']:.2f} | {s['rate_med']:.2f} {s['rate_p90']:.2f} | {s['dwell_med']:.2f} {s['dwell_min']:.2f} | {s['flat']:.2f}")
    rep['gain'] = out

    print('--- severe under-estimate of the raw D4 24 window (e < ln .6): how often, and which part made it')
    sev = {}
    for g in H.GROUPS+['all']:
        rr = [r for r in res if g == 'all' or r['meta']['group'] == g]
        e = np.concatenate([r['e'] for r in rr]); en = np.concatenate([r['e_n'] for r in rr]); e0 = np.concatenate([r['e_0'] for r in rr])
        m = np.isfinite(e) & np.isfinite(en) & np.isfinite(e0)
        low = m & (e < math.log(.6))
        n_low = int(low.sum())
        sev[g] = dict(windows=int(m.sum()), severe=n_low, share=float(n_low/m.sum()), noise_only_low=float(np.mean(e0[low] > math.log(.8))) if n_low else None, model_low=float(np.mean(e0[low] < math.log(.6))) if n_low else None, en_med=float(np.median(en[low])) if n_low else None, e0_med=float(np.median(e0[low])) if n_low else None)
        print(g, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in sev[g].items()})
    rep['severe'] = sev

    print('--- what signals a sustained change (0.2 s scale, |rate|>.5/s) against flat (<.05/s), no event; AUC')
    aucs = {}
    for g in H.GROUPS+['all']:
        P = {q: [] for q in ('res', 'hp', 'innov')}; N = {q: [] for q in P}
        for r in res:
            if g != 'all' and r['meta']['group'] != g:
                continue
            ell, t, c = r['ell'], r['t'], r['contact']
            s = uniform_filter1d(ell, BOX, mode='nearest'); rate = np.gradient(s, t)
            ok = c & (t > 1.)
            for ev in r['meta']['event_info']:
                ok &= ~((t >= ev['start']) & (t <= ev['end']+.3))
            on = ok & (np.abs(rate) > .5); off = ok & (np.abs(rate) < .05)
            old = np.full(len(t), np.nan)
            old[48:] = np.nanmean(np.lib.stride_tricks.sliding_window_view(r['lnD'], 25)[:len(t)-48], axis=1)
            for q, v in (('res', np.abs(r['res'])), ('hp', np.abs(r['hp'])), ('innov', np.abs(r['lnD']-old))):
                P[q].append(v[on]); N[q].append(v[off])
        aucs[g] = {q: F.auc(np.concatenate(P[q]), np.concatenate(N[q])) for q in P}
        aucs[g]['on'] = int(sum(len(x) for x in P['res']))
        print(g, {q: (round(v, 3) if isinstance(v, float) else v) for q, v in aucs[g].items()})
    rep['signal_auc'] = aucs
    print('--- window lag and force channel, no event, t>1, contact: e_0 = -tau * d ln k/dt (noise-free window); res = g * F * d ln k/dt + slow motion')
    lag = {}
    for g in H.GROUPS:
        rate, e0, res_hp, frate = [], [], [], []
        for r in res:
            if r['meta']['group'] != g:
                continue
            ell, t, c = r['ell'], r['t'], r['contact']
            s = uniform_filter1d(ell, BOX, mode='nearest'); rr_ = np.gradient(s, t)
            ok = c & (t > 1.) & np.isfinite(r['e_0']) & np.isfinite(r['res'])
            for ev in r['meta']['event_info']:
                ok &= ~((t >= ev['start']) & (t <= ev['end']+.3))
            # Remove each case's own offset and the slow (breathing) part of both signals, over 1 s.
            hp = r['res']-uniform_filter1d(np.nan_to_num(r['res']), 125, mode='nearest')
            rh = rr_-uniform_filter1d(rr_, 125, mode='nearest')
            rate.append(rr_[ok]); e0.append(r['e_0'][ok]-np.median(r['e_0'][ok])); res_hp.append(hp[ok]); frate.append((r['F']*rh)[ok])
        rate, e0, res_hp, frate = map(np.concatenate, (rate, e0, res_hp, frate))
        tau = -float(rate@e0/(rate@rate))
        cc = float(np.corrcoef(rate, e0)[0, 1])
        gain = float(frate@res_hp/(frate@frate))
        rc = float(np.corrcoef(frate, res_hp)[0, 1])
        lag[g] = dict(tau_s=tau, corr=cc, rate_rms=float(np.sqrt(np.mean(rate**2))), force_gain=gain, force_corr=rc, res_rms=float(np.std(res_hp)), frate_rms=float(np.std(frate)))
        print(g, {k: round(v, 3) for k, v in lag[g].items()})
    rep['lag_force'] = lag
    (HERE/'results/memory_scene_curves.json').write_text(json.dumps(rep, indent=1, default=float), encoding='utf-8')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    picks = (0, 7, 36, 26)
    fig, ax = plt.subplots(3, len(picks), figsize=(17, 8.5), sharex='col', gridspec_kw=dict(height_ratios=[2.2, 1, 1]))
    INK, RAW, MEM = '#222222', '#B8B8B8', '#1F6FB2'
    for col, idx in enumerate(picks):
        r = next(x for x in res if x['meta']['index'] == idx)
        with np.load(TICKS/f'Huber-D4_{idx:02d}.npz') as z:
            raw, imm = z['D4_nogate__raw'].astype(float), z['D4_imm__raw'].astype(float)
        t = r['t']; n = min(len(t), len(raw))
        a = ax[0, col]
        a.plot(t[:n], np.log(np.where(raw[:n] > 0, raw[:n], np.nan)), color=RAW, lw=.8, label='D⁴ window, no memory')
        a.plot(t[:n], np.log(np.where(imm[:n] > 0, imm[:n], np.nan)), color=MEM, lw=1.4, label='D⁴ + IMM')
        a.plot(t[:n], r['ell'][:n], color=INK, lw=1.6, label='true k_eff')
        a.set_title(f"{r['meta']['group']}, case {idx}", fontsize=11, loc='left'); a.grid(alpha=.25); a.spines[['top', 'right']].set_visible(False)
        a.set_ylim(np.nanpercentile(r['ell'][:n], .5)-.9, np.nanpercentile(r['ell'][:n], 99.5)+.6)
        for ev in r['meta']['event_info']:
            for row in ax[:, col]:
                row.axvspan(ev['start'], ev['end'], color='#E9D8A6', alpha=.45, lw=0)
        b = ax[1, col]
        b.plot(t[:n], (np.log(np.where(raw[:n] > 0, raw[:n], np.nan))-r['ell'][:n]), color=RAW, lw=.8)
        b.plot(t[:n], (np.log(np.where(imm[:n] > 0, imm[:n], np.nan))-r['ell'][:n]), color=MEM, lw=1.2)
        b.axhline(0, color=INK, lw=.6); b.set_ylim(-1.3, 1.3); b.grid(alpha=.25); b.spines[['top', 'right']].set_visible(False)
        c = ax[2, col]
        c.plot(t[:n], r['res'][:n], color='#555555', lw=.8); c.grid(alpha=.25); c.spines[['top', 'right']].set_visible(False)
        c.set_xlabel('time (s)')
    ax[0, 0].set_ylabel('ln k (N/m)'); ax[1, 0].set_ylabel('ln(estimate / true)'); ax[2, 0].set_ylabel('dF − k̂ dy (N/s)')
    ax[0, 0].legend(frameon=False, fontsize=9, loc='upper left')
    fig.tight_layout()
    out_dir = ROOT/'research/figures/memory_scene'; out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir/'curves.png', dpi=130)


if __name__ == '__main__':
    main()
