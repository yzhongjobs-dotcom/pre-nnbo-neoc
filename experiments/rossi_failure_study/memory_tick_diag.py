"""Why does the tick Kalman still hit the floor and do worse on cough? Replay on recorded observations.

Two explanations are separated: restart after contact loss (implementation) and sharp d in a Gaussian D4 state (model).
If the second holds, a third mode with a large disturbance innovation, same k noise as the quiet mode, should remove the collapse in events and leave plateaus alone:
the model then has three explanations for a large innovation, namely noise, a change of k, and a sudden d. Prediction fixed before running.
Writes results/memory_tick_diag.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import memory_tick_kalman as T

HERE = Path(__file__).resolve().parent
SEVERE = math.log(.5)
Q_RATIO = 100.     # variance ratio of the disturbance innovation in the third mode


class ThreeMode(T.TickKalman):
    """Quiet, k changing, d jumping. The third mode shares the quiet k noise; Q_d is multiplied."""
    def __init__(self, q_rel=(.002, .5), ratio=Q_RATIO, **kw):
        super().__init__(q_rel=(q_rel[0], q_rel[1], q_rel[0]), **kw)
        self.ratio = ratio

    def step(self, dy, z):
        if not np.isfinite(z):
            return np.nan
        dy = 0. if not np.isfinite(dy) else dy*T.UNIT
        n = self.modes
        Tm = np.full((n, n), (1-self.stay)/(n-1)); np.fill_diagonal(Tm, self.stay)
        joint = self.w[:, None]*Tm; pw = joint.sum(axis=0); mix = joint/pw[None, :]
        mx = [sum(mix[a, j]*self.x[a] for a in range(n)) for j in range(n)]
        mP = [sum(mix[a, j]*(self.P[a]+np.outer(self.x[a]-mx[j], self.x[a]-mx[j])) for a in range(n)) for j in range(n)]
        ll = np.zeros(n); xs, Ps = [], []
        for j in range(n):
            A = self.A.copy(); A[1, 0] = dy
            kp = max(mx[j][0], T.KLO/T.UNIT)
            Q = self.Qd*(self.ratio if j == 2 else 1.)
            Q = Q.copy(); Q[0, 0] += self.rel[j]*kp*kp*T.TC
            x = A@mx[j]; P = A@mP[j]@A.T+Q
            if not self.started:
                x = mx[j].copy(); P = mP[j].copy()
            h = np.zeros(self.n); h[self.id] = 1.; h[1] = 1.
            e = z-h@x; Ph = P@h; S = h@Ph+self.sigma**2
            K = Ph/S; B = np.eye(self.n)-np.outer(K, h)
            P2 = B@P@B.T+self.sigma**2*np.outer(K, K)
            xs.append(x+K*e); Ps.append(.5*(P2+P2.T)); ll[j] = -.5*(math.log(2*math.pi*S)+e*e/S)
        self.x, self.P = xs, Ps; self.started = True
        lw = np.log(pw)+ll; lw -= lw.max(); self.w = np.exp(lw); self.w /= self.w.sum()
        return float(np.clip(T.UNIT*sum(self.w[j]*self.x[j][0] for j in range(n)), T.KLO, T.KHI))


CFG = dict(form='incr', d='d4', sigma=T.SIG_TRUE)


def run(path):
    data = H.load(path); meta = data['meta']; log, masks = H.masks_for(data)
    dY, z = T.ticks(data); n = len(z); k = log[:n, 26]; t = log[:n, 0]
    ev = np.zeros(n, bool)
    for e in meta['event_info']:
        ev |= (t >= e['start']) & (t <= e['end']+.3)
    out = {}
    since = np.full(n, 10**6)
    for name, f in (('T2', T.TickKalman(q_rel=(.002, .5), **CFG)), ('T3', ThreeMode(**CFG))):
        est = np.full(n, np.nan); last_reset = -10**6; first = True
        for j in range(n):
            if not np.isfinite(z[j]):
                if f.started:
                    f.reset(f.x[0][0]); last_reset = j
                continue
            if first:
                last_reset = j; first = False
            est[j] = f.step(dY[j], z[j]); since[j] = j-last_reset
        out[name] = est
    win = np.full(n, np.nan); S = H.Solver()
    for j in range(n):
        Y, zz = H.tail(data['Y'][j], 24), H.tail(data['z'][j], 24)
        if len(Y) >= 8:
            win[j] = H.clip(S.fit(Y, zz, iters=1)['k'])
    out['window'] = win
    good = k > 200.
    res = dict(index=meta['index'], group=meta['group'], names={})
    for name, est in out.items():
        e = np.log(np.where(est > 0, est, np.nan))-np.log(np.where(good, k, np.nan))
        ok = np.isfinite(e) & (t >= 1.)
        low = ok & (e < SEVERE)
        res['names'][name] = dict(ticks=int(ok.sum()), severe=int(low.sum()), severe_event=int((low & ev).sum()), severe_after_reset=int((low & (since < 40)).sum()), floor=int((ok & (est <= 200.5)).sum()),
                                  plateau=e[masks['plateau'][:n] & ok].tolist(), transition=e[masks['transition'][:n] & ok].tolist(), event=e[masks['event'][:n] & ok].tolist())
    return res


def main():
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(run, sorted(H.WIN.glob('case_*.npz'))))
    summary = {}
    print('severe windows (estimate below half of truth, t>1): total, in events, within 40 ticks after a restart, at floor; by group')
    for name in ('window', 'T2', 'T3'):
        print('---', name)
        summary[name] = {}
        for g in H.GROUPS+['all']:
            rr = [r['names'][name] for r in res if g == 'all' or r['group'] == g]
            s = {q: sum(x[q] for x in rr) for q in ('ticks', 'severe', 'severe_event', 'severe_after_reset', 'floor')}
            summary[name][g] = s
            print(f"  {g:12s} severe {s['severe']:5d} ({100*s['severe']/s['ticks']:.2f}%) event {s['severe_event']:4d} after restart {s['severe_after_reset']:4d} floor {s['floor']:4d}")
        for seg in ('plateau', 'transition', 'event'):
            e = np.concatenate([np.array(r['names'][name][seg]) for r in res])
            summary[name][seg] = dict(med=float(np.median(np.abs(e))), p90=float(np.percentile(np.abs(e), 90)), bias=float(np.median(e)), below=float(np.mean(e < math.log(.85))))
            print(f"  {seg:10s} {summary[name][seg]['med']:.3f}/{summary[name][seg]['p90']:.3f} bias {summary[name][seg]['bias']:+.3f} below {summary[name][seg]['below']*100:.0f}%")
    (HERE/'results/memory_tick_diag.json').write_text(json.dumps(summary, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
