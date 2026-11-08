"""Why scene A's worst episode (breath_rib, rib, neural_white, 28 deg, i=6) runs away under the current Method.

Reruns the frozen overlay case and records the policy's stiffness state each tick. Variants replace only one thing:
  base      the adopted default, must reproduce the saved log
  true_k    window estimate = true probe-frequency gain, same gates (diagnostic=True supplies the truth)
  no_ff     feed-forward increment forced to zero
Usage: python -X utf8 scene_a_worst_diagnosis.py
"""
import json
import sys
from pathlib import Path

import numpy as np

import overlay_compare as O
import method3 as M
import physics_plant3 as P3

B = O.B
HERE = Path(__file__).resolve().parent
CELL = ('breath_rib', 'rib', 'neural_white', 28., int(sys.argv[1]) if len(sys.argv) > 1 else 6)


class Rec:
    def __init__(self, pol):
        self.pol = pol; self.rows = []

    def reset(self, platform):
        self.pol.reset(platform)

    def __call__(self, o):
        u = self.pol(o)
        p = self.pol
        self.rows.append((float(o.log[o.n-1, 0]), float(p.kbar), float(p.k_good), float(p.k_hat), float(p.k_check), float(getattr(p, 'dq_last', np.nan))))
        return u

    def __getattr__(self, name):
        return getattr(self.pol, name)


def run(mode):
    O.setup()
    cfg = O.chosen()['Huber-D4']
    sc = B.scene(CELL, O.BASE, 0.)
    pol = B.build('Huber-D4', cfg, sc)
    if mode == 'true_k':
        est = pol.estimator; original = est.estimate
        def estimate(o, n, policy, original=original):
            k, sig, chi = original(o, n, policy)
            truth = getattr(o, 'truth_k', 0.)
            if np.isfinite(k) and truth > 200.:
                sig = sig*truth/k if k > 0 else sig
                k = truth
            return k, sig, chi
        est.estimate = estimate
    if mode == 'no_ff':
        pol.p = dict(pol.p, ff_gain=0.) if 'ff_gain' in pol.p else pol.p
    rec = Rec(pol)
    log = P3.episode(sc, rec, O.sensor_seed(CELL), diagnostic=mode != 'base')
    return sc, pol, np.array(rec.rows), log


def main():
    out = {}
    for mode in ('base', 'true_k'):
        sc, pol, rows, log = run(mode)
        m = M.metric(log, pol, sc)
        t = log[:, 0]; f = log[:, 18]; k = log[:, 21]
        out[mode] = dict(J=float(m['J']), peak=float(np.nanmax(f[t >= 1])), rmse=float(m['rmse']))
        print(mode, out[mode], flush=True)
        if mode == 'base':
            saved = np.load(HERE/'results'/'baseline_logs'/f'Huber-D4_{CELL[0]}_{CELL[1]}_{CELL[2]}_{int(CELL[3])}_{CELL[4]}.npz')
            cols = list(saved['cols'])
            print('reproduces saved log:', np.allclose(log[:, cols], saved['log'], equal_nan=True), flush=True)
            kt = np.interp(rows[:, 0], t, k)
            for lo, hi in ((1, 3), (3, 5), (5, 6), (6, 7), (7, 8), (8, 9.3)):
                w = (rows[:, 0] >= lo) & (rows[:, 0] < hi)
                ww = (t >= lo) & (t < hi)
                print(f' t {lo}-{hi}: k_true med {np.median(kt[w]):.0f} kbar med {np.median(rows[w,1]):.0f} k_good med {np.nanmedian(rows[w,2]):.0f} '
                      f'kappa=k/kbar med {np.median(kt[w]/rows[w,1]):.2f} kappa_g=k/k_good med {np.nanmedian(kt[w]/rows[w,2]):.2f} max {np.nanmax(kt[w]/rows[w,2]):.2f} '
                      f'F min {f[ww].min():.1f} max {f[ww].max():.1f}', flush=True)
            w = (t >= 5.5) & (t < 9.2)
            e = f[w]-f[w].mean(); sp = np.abs(np.fft.rfft(e*np.hanning(len(e))))**2; fr = np.fft.rfftfreq(len(e), np.median(np.diff(t)))
            top = np.argsort(sp)[::-1][:5]
            print(' force spectrum peaks Hz', np.round(fr[top], 1), flush=True)
            np.savez_compressed(HERE/'results'/'scene_a_worst_trace.npz', rows=rows, log=log[:, [0, 18, 21, 2, 11, 6]])
    (HERE/'results'/'scene_a_worst_diagnosis.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
