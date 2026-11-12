"""Diagnostic ceiling on the A_large runs: the window output corrected by the encoder trajectory, with a per-case gain chosen from truth.

C_hat(t) = mean over the last N ticks of [C_probe - beta*ybar24] + beta*y(t). N = 1 only removes the window's lag; larger N also averages
the level. beta = 1/g with g picked per case from a grid by total rms (privileged, so this is a ceiling for the simple form, not a method).
Run: python -X utf8 dual_channel_delag_check.py
"""
from pathlib import Path

import numpy as np

import dual_channel_replay as D
import dual_channel_v3_replay as V


def comp(c, N, beta):
    t, ys6, fs, m, R, valid, finite = c['data']; y = c['log'][:, 2]; yw = D.trailing(y, 24)
    Cp = np.where(np.isfinite(m), D.trailing(fs, 24)/np.exp(np.nan_to_num(m)), np.nan); r = Cp-beta*yw
    for j in range(1, len(r)):
        if not np.isfinite(r[j]):
            r[j] = r[j-1]
    lvl = D.trailing(np.nan_to_num(r, nan=np.nanmedian(r)), N)
    return np.log(np.maximum(fs, 1.))-np.log(np.maximum(lvl+beta*y, 1e-7))


if __name__ == '__main__':
    cases, _ = D.load_cases(Path(__file__).resolve().parents[2])
    cs = [c for c in cases if c['group'] == 'A_large']
    for c in cs:
        best = (9, 0)
        for g in (1, 1.5, 2, 3, 4, 6, 8, 12, 20, 1e9):
            s = D.rms((comp(c, 24, 1/g)-c['truth'])[c['masks']['outside_event']])
            if s < best[0]:
                best = (s, g)
        c['y'] = np.array(c['log'][:, 2], copy=True); c['yw'] = D.trailing(c['y'], 24)
        c['pred']['A1'] = D.probe_filter(c['data'], (0., .01, .99))
        c['pred']['delag'] = comp(c, 1, 1/best[1]); c['pred']['lvl24'] = comp(c, 24, 1/best[1])
        c['pred']['KFloose'] = V.run(c, (1., 1., 1e-7), V.peak_freq(c['y']))[0]
        print(c['index'], 'gain picked', best[1])
    segs = ('plateau', 'slow', 'rise', 'fall', 'reversal')
    for name, idx in (('fast 36,37,39,40', (36, 37, 39, 40)), ('slow 38,41', (38, 41))):
        sub = [c for c in cs if c['index'] in idx]
        print(name, '| median|e| / P90|e|: plateau slow rise fall reversal | mean e on rise, fall')
        for a in ('A0', 'A1', 'delag', 'lvl24', 'KFloose'):
            row = []
            for seg in segs:
                e = np.concatenate([(c['pred'][a]-c['truth'])[c['masks'][seg]] for c in sub])
                row.append(f'{np.median(np.abs(e)):.3f}/{np.percentile(np.abs(e), 90):.2f}' if len(e) else '  -  ')
            b = [np.mean(np.concatenate([(c['pred'][a]-c['truth'])[c['masks'][s]] for c in sub])) for s in ('rise', 'fall')]
            print(f'   {a:8s}', ' '.join(row), f'| {b[0]:+.3f} {b[1]:+.3f}')
