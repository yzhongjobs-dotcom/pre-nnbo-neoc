"""Does the excitation's exogeneity already help on recorded data? A test of the new idea that needs no new simulation.

On every recorded 24-tick window, D4 uses all of M Y as the regressor. Y also holds motion the robot made in reaction to force
(feedback, feed-forward, the closed-loop response to d and to sensor noise), which is endogenous. Two-stage least squares with the
commanded excitation as instrument, Yhat = P_A M Y with A = {sin, cos}(2 pi f0 t), k = Yhat' M z / Yhat' M Y, keeps only the part of Y
explained by the excitation. It cannot remove a disturbance at the probe frequency (that needs a coded probe), but it removes
endogenous regressor content at every other frequency.

Predictions fixed before running:
 - plateau of ordinary scenes: D4 and 2SLS nearly equal, both noise-limited; 2SLS variance a little larger.
 - where D4's noise-free error is large (A_large, combo, cough events) and comes from endogenous Y content, 2SLS reduces the noise-free error.
 - what 2SLS leaves in noise-free windows bounds what a coded probe could still remove (in-band d) together with in-window k change.
Also reports the bias of the sensor-noise part (closed-loop endogeneity of noise). Truth only for scoring.
Writes results/iv_separation_check.json.
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
F0 = 20.8
NW = 24


def fits(Y, z, t):
    D, _ = H.Solver().bands(len(Y))
    Mm = np.eye(len(Y))-np.linalg.inv(np.eye(len(Y))+H.RHO*D.T@D)
    my = Mm@Y; mz = Mm@z
    kd = float(my@z)/float(my@Y)
    A = np.c_[np.sin(2*np.pi*F0*t), np.cos(2*np.pi*F0*t)]
    MA = Mm@A
    yh = MA@np.linalg.lstsq(MA, my, rcond=None)[0]
    kiv = float(yh@mz)/float(yh@my)
    share = float(yh@yh)/float(my@my)      # fraction of M Y explained by the excitation
    return kd, kiv, share


def case(path):
    data = H.load(path); meta = data['meta']; log, masks = H.masks_for(data)
    t = log[:, 0]; k = log[:, 26]; n = len(t)
    out = {q: np.full(n, np.nan) for q in ('D', 'IV', 'D0', 'IV0', 'share')}
    for j in range(NW+1, n):
        Y, z = H.tail(data['Y'][j], NW), H.tail(data['z'][j], NW)
        Y0, z0 = H.tail(data['Y0'][j], NW), H.tail(data['z0'][j], NW)
        if len(Y) != NW or len(Y0) != NW or k[j] <= 200.:
            continue
        tw = t[j-NW:j]
        a = fits(Y, z, tw); b = fits(Y0, z0, tw)
        if min(a[0], a[1], b[0], b[1]) <= 0:
            continue
        out['D'][j], out['IV'][j], out['share'][j] = math.log(a[0]/k[j]), math.log(a[1]/k[j]), a[2]
        out['D0'][j], out['IV0'][j] = math.log(b[0]/k[j]), math.log(b[1]/k[j])
    return dict(group=meta['group'], out=out, masks={s: masks[s] for s in ('plateau', 'transition', 'event')})


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(case, sorted(H.WIN.glob('case_*.npz'))))
    rep = {}
    print('group seg | windows | share of M Y from excitation med | noise-free |e| med/P90: D4, 2SLS | noisy |e| med/P90: D4, 2SLS | noise part bias D4, 2SLS | noisy below .5: D4, 2SLS')
    for g in H.GROUPS+['all']:
        for s in ('plateau', 'transition', 'event'):
            c = {q: [] for q in ('D', 'IV', 'D0', 'IV0', 'share')}
            for r in res:
                if g != 'all' and r['group'] != g:
                    continue
                m = r['masks'][s] & np.isfinite(r['out']['D']) & np.isfinite(r['out']['D0'])
                for q in c:
                    c[q].append(r['out'][q][m])
            c = {q: np.concatenate(v) for q, v in c.items()}
            if len(c['D']) < 50:
                continue
            f = lambda x: (float(np.median(np.abs(x))), float(np.percentile(np.abs(x), 90)))
            row = dict(n=len(c['D']), share=float(np.median(c['share'])), D0=f(c['D0']), IV0=f(c['IV0']), D=f(c['D']), IV=f(c['IV']),
                       noise_bias_D=float(np.median(c['D']-c['D0'])), noise_bias_IV=float(np.median(c['IV']-c['IV0'])),
                       low_D=float(np.mean(c['D'] < math.log(.5))), low_IV=float(np.mean(c['IV'] < math.log(.5))))
            rep[f'{g}|{s}'] = row
            print(f"{g:11s} {s:10s} | {row['n']:5d} | {row['share']:.2f} | {row['D0'][0]:.3f}/{row['D0'][1]:.3f} {row['IV0'][0]:.3f}/{row['IV0'][1]:.3f} | {row['D'][0]:.3f}/{row['D'][1]:.3f} {row['IV'][0]:.3f}/{row['IV'][1]:.3f} | {row['noise_bias_D']:+.3f} {row['noise_bias_IV']:+.3f} | {row['low_D']*100:.1f}% {row['low_IV']*100:.1f}%")
    (HERE/'results/iv_separation_check.json').write_text(json.dumps(rep, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
