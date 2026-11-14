"""How much stiffness error does the motion feed-forward loop tolerate when the tracker is fed an average of the motion measurement?

Same frozen linearised platform model as ../motion_aware_scan/ff_loop.py (property 4 of the paper draft). The tracker receives the mean of
the last n motion measurements Z (n = 6 ticks is one 20.83 Hz probe cycle) and predicts `extra` ticks further than before; the tracker
gain is unchanged. Reports the stable interval of k/k_ref around 1 and the frequency of the dominant mode at the baseline boundaries.
Run from this folder: python -X utf8 ff_band_separation.py
"""
import json, math, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'motion_aware_scan'))
import ff_loop as FL

Tc, tau_s, NS, DEC, FRAC, At, Q1 = FL.Tc, FL.tau_s, FL.NS, FL.DEC, FL.FRAC, FL.At, FL.Q1
M1, D1, KPF1 = FL.M1, FL.D1, FL.KPF1
NX = 18
GRID = np.round(np.concatenate([np.arange(.2, 1., .05), np.arange(1., 4.01, .1)]), 3)


def gain(q, k_phys, r_scale=1.):
    R = (FL.NOISE/k_phys)**2*r_scale; Pm = np.eye(3)*1e-6
    for _ in range(4000):
        Pp = At@Pm@At.T+q*Q1; K = Pp[:, 0]/(Pp[0, 0]+R); Pn = (np.eye(3)-np.outer(K, [1, 0, 0]))@Pp
        if np.max(np.abs(Pn-Pm)) < 1e-22*max(1, np.max(np.abs(Pn))):
            break
        Pm = Pn
    return K


def tick(x, kappa, kg, K, n, extra):
    y, cp, Ff, I, wd, cmd, mp = x[:7]; Fp, yp = x[7:9]; xh = x[9:12]; spp = x[12]; zh = x[13:18]
    e, ep = -Ff, -mp; eD = math.exp(-D1*Tc/M1); wd2 = wd*eD+(1-eD)*e/D1
    fin = (Ff-DEC*Fp)/(1-DEC); ybar = yp+FRAC*(y-yp); Z = kg/kappa*fin-ybar
    L = 1.5*Tc+extra*Tc; Zin = (Z+zh[:n-1].sum())/n
    xm = At@xh; xh2 = xm+K*(Zin-xm[0]); sp = np.array([1, L+tau_s, L*L/2+tau_s*L])@xh2
    cmd2 = cmd+0.5*Tc*(wd+wd2)+KPF1*(e-ep)-(sp-spp); cp2 = cmd
    for _ in range(NS):
        Ff = Ff+FL.ALPHA*(kappa*y-Ff); y = y+FL.BETA*(cp2-y)
    return np.concatenate([[y, cp2, Ff, I, wd2, cmd2, -e], [x[2], x[0]], xh2, [sp], [Z], zh[:-1]])


def dominant(kappa, kg, K, n, extra):
    A = np.zeros((NX, NX))
    for i in range(NX):
        x = np.zeros(NX); x[i] = 1.; A[:, i] = tick(x, kappa, kg, K, n, extra)
    keep = [i for i in range(NX) if i != 3]
    ev = np.linalg.eigvals(A[np.ix_(keep, keep)]); ev = ev[np.abs(ev-1) > 1e-7]
    j = int(np.argmax(np.abs(ev)))
    return float(abs(ev[j])), float(abs(np.angle(ev[j]))/(2*np.pi*Tc))


def interval(k, q, kappa, n, extra):
    r = [dominant(kappa, kg, gain(q, k/kg), n, extra)[0] for kg in GRID]; i1 = int(np.argmin(np.abs(GRID-1.)))
    if r[i1] >= 1:
        return None
    j = i1
    while j > 0 and r[j-1] < 1: j -= 1
    h = i1
    while h < len(r)-1 and r[h+1] < 1: h += 1
    return float(GRID[j]), float(GRID[h])


if __name__ == '__main__':
    qs, qm = FL.P['q_smooth'], 3.335          # tracker settings in use (method3.TRACKER), r_scale 1
    out = dict(modes={}, intervals={})
    print('dominant mode at the baseline boundaries (kappa = 1, maneuver mode): |eig|, Hz')
    for k in (1e3, 1e4, 1e5):
        for kg in (.6, .7, 1.3, 1.5):
            r, f = dominant(1., kg, gain(qm, k/kg), 1, 0.)
            out['modes'][f'{k:g}|{kg}'] = (r, f); print(f'   k {k:g} k/k_ref {kg}: {r:.3f} at {f:.1f} Hz')
    cfgs = [(1, 0.), (6, 0.), (6, 1.), (6, 1.5), (6, 2.5), (3, 0.), (3, 1.)]
    print('stable k/k_ref interval; columns (n ticks averaged, extra horizon in ticks): '+' | '.join(f'{n},+{e}' for n, e in cfgs))
    for kappa in (1., 2., .5):
        for k in (1e3, 1e4, 1e5):
            for mode, q in (('maneuver', qm), ('smooth', qs)):
                row = [interval(k, q, kappa, n, e) for n, e in cfgs]
                out['intervals'][f'{kappa}|{k:g}|{mode}'] = row
                print(f'kappa {kappa} k {k:6g} {mode:8s} '+' | '.join(str(v) for v in row))
    json.dump(out, open(Path(__file__).resolve().parent/'results'/'ff_band_separation.json', 'w'), indent=1)
