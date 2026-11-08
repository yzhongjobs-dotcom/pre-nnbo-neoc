"""Rescore all estimators against the effective stiffness k_eff (log column 26) on identical observations.

Plants: `old` (quasi-static arm, the current default) and `dyn` (dynamic arm zeta 0.1, sensor free of tool ringing), a sensitivity comparison, since the arm's high-frequency parameters have no measurement.
Usage: python -X utf8 keff_rescore.py
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import id_layer_root_cause as I

L, SE, R, C = I.L, I.SE, I.R, I.C
HERE = Path(__file__).resolve().parent
FIG = HERE.parents[1]/'research/figures/id_layer_keff'
PLANTS = dict(old={}, dyn=dict(arm_zeta=.1, sensor_inertia=0., tool_mass_hat=0.))


def job(a):
    case, plant = a
    original = L.MC.scene_for
    def scene_for(*x, **y):
        sc = original(*x, **y); sc.update(PLANTS[plant]); return sc
    L.MC.scene_for = scene_for
    label = case['group']+'_'+case['cell'][2]+'_i'+str(case['cell'][-1])
    cfg = dict(label=label, cell=case['cell'], base=C.BASE, variant=case['variant'], seed=case['seed'])
    _, sc, log, windows, commands, steps, contexts = L.episode(cfg, full=True)
    count = len(log); t = log[:, 0]
    raw = np.full(count, np.nan)
    for w in windows:
        raw[w['n']-1] = w['k_hat']
    hats = dict(raw=np.log(raw), gated=np.log(np.array([s['kgood'] for s in steps])))
    for kind in I.MEM:
        hats[kind] = I.filter_replay(kind, windows, count)[0]
    out = dict(label=label, plant=plant, plateau={}, plateau_old={}, events=[], events_old=[], keff_over_old=float(np.median((log[:, 26]/np.maximum(log[:, 21], 1e-9))[log[:, 21] > 200.])))
    for tag, col in (('', 26), ('_old', 21)):
        lg = log.copy(); lg[:, 21] = log[:, col]
        truth = np.log(np.maximum(lg[:, 21], 1e-30))
        _, masks = SE.segment_masks(lg, sc)
        _, _, events, _ = R.transitions(lg)
        out['plateau'+tag] = {k: I.stats((h-truth)[masks['plateau']]) for k, h in hats.items()}
        rows = I.event_scores(lg, hats, events)
        out['events'+tag] = rows
        if not tag:
            out['n_plateau'] = int(masks['plateau'].sum())
            mt = masks['transition']
            out['trans_frac'] = {k: dict(below1=float(np.nanmean((h-truth)[mt] < 0)), below85=float(np.nanmean((h-truth)[mt] < np.log(.85))), mean=float(np.nanmean((h-truth)[mt])), sd=float(np.nanstd((h-truth)[mt]))) if mt.any() else None for k, h in hats.items()}
    if plant == 'old' and 'slow' in case['group'] and case['cell'][2] == 'neural_white':
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        FIG.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        ax[0].plot(t, log[:, 26]/1e3, color='#AA0000', lw=2.2, label='k_eff')
        ax[0].plot(t, log[:, 21]/1e3, color='#AA0000', lw=1, ls=':', label='old |K|')
        for k, h in hats.items():
            ax[0].plot(t, np.exp(h)/1e3, color=I.COLORS[k], lw=.9, label=k); ax[1].plot(t, h-np.log(np.maximum(log[:, 26], 1e-30)), color=I.COLORS[k], lw=.9)
        ax[0].set_ylabel('Stiffness (kN/m)'); ax[1].set_ylabel('ln(estimate / k_eff)'); ax[1].set_xlabel('Time (s)'); ax[1].axhline(0, color='k', lw=.6); ax[1].set_ylim(-1.2, .8); ax[0].legend(ncol=6, fontsize=8)
        fig.tight_layout(); fig.savefig(FIG/f'{label}.png', dpi=120); plt.close(fig)
    return out


if __name__ == '__main__':
    start = time.perf_counter()
    cases = [c for c in C.cases() if c['group'] in ('A', 'A_slow_thin') and c['cell'][-1] in (0, 2, 6, 1, 4, 5)]
    jobs = [(c, p) for p in PLANTS for c in cases]
    print('jobs', len(jobs), flush=True)
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, jobs, chunksize=1))
    (HERE/'results'/'keff_rescore.json').write_text(json.dumps(rows, indent=1, default=float))
    print('seconds', round(time.perf_counter()-start))
