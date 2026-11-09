"""Tick-level Kalman for the stiffness, on recorded observations from results/huber_windows, no plant.

State (k, p, d, d-1, d-2, d-3): z_i = p_i + d_i + e_i, p_i = p_{i-1} + k_i dY_i, D4 prior on d (or the order-2 model of JointKalman).
Linear in the state given dY, so each mode is an exact Kalman filter. Truth is for scoring only.

Stage 1, ablation ladder on ten cases: why is the existing tick Kalman low, one change at a time. Prediction fixed before running:
the D4 state removes most of the negative bias, because motion below 14 Hz stops informing k; the increment form matters mainly when k changes.
Writes results/memory_tick_kalman.json.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import huber_prior_check as H
import memory_segment_eval as SE
import rib_transition_eval as R

HERE = Path(__file__).resolve().parent
TC, RHO = .008, H.RHO
SIG_TRUE, SIG_NOM = H.SIGMA, .04
UNIT = 1000.
KLO, KHI = 200., 3e5
LADDER_CASES = (0, 3, 7, 10, 12, 16, 24, 26, 33, 36, 38)
SEGS = ('plateau', 'transition', 'event')


def ticks(data):
    """Per row: dY (between consecutive rows, from one window so the origin cancels) and z, with validity; row r is available at log index r+1."""
    n = len(data['Y'])
    dY = np.full(n, np.nan); z = np.full(n, np.nan)
    y = data['Y']; zz = data['z']
    for j in range(n):
        if np.isfinite(y[j, -1]):
            z[j] = zz[j, -1]
            if np.isfinite(y[j, -2]):
                dY[j] = y[j, -1]-y[j, -2]
    return dY, z


class TickKalman:
    def __init__(self, form='incr', d='d4', sigma=SIG_TRUE, q_abs=None, q_rel=None, q_d=300., stay=.98**(1/6)):
        self.form, self.d, self.sigma, self.q_abs, self.q_rel, self.q_d, self.stay = form, d, sigma, q_abs, q_rel, q_d, stay
        self.modes = 1 if q_rel is None or np.isscalar(q_rel) else len(q_rel)
        self.rel = None if q_rel is None else np.atleast_1d(q_rel).astype(float)
        nd = 4 if d == 'd4' else 2
        self.n = 1+(1 if form == 'incr' else 0)+nd
        self.id = 1+(1 if form == 'incr' else 0)     # index of d
        A = np.eye(self.n)
        Qd = np.zeros((self.n, self.n))
        i = self.id
        if d == 'd4':
            A[i:i+4, i:i+4] = 0.
            A[i, i:i+4] = [4., -6., 4., -1.]
            A[i+1:i+4, i:i+3] = np.eye(3)
            Qd[i, i] = sigma*sigma/RHO
        else:
            A[i, i+1] = TC
            Qd[i:i+2, i:i+2] = q_d*q_d*np.array([[TC**3/3, TC*TC/2], [TC*TC/2, TC]])
        self.A = A; self.Qd = Qd
        self.reset()

    def reset(self, k0=30.):
        P0 = np.eye(self.n)*100.**2
        P0[0, 0] = 300.**2
        if self.form == 'incr':
            P0[1, 1] = 0.
        self.x = [np.r_[k0, np.zeros(self.n-1)] for _ in range(self.modes)]
        self.P = [P0.copy() for _ in range(self.modes)]
        self.w = np.ones(self.modes)/self.modes
        self.Yabs = 0.
        self.started = False

    def step(self, dy, z):
        if not np.isfinite(z):
            return np.nan
        dy = 0. if not np.isfinite(dy) else dy*UNIT
        if self.modes > 1:
            T = np.full((self.modes, self.modes), (1-self.stay)/(self.modes-1)); np.fill_diagonal(T, self.stay)
            joint = self.w[:, None]*T; pw = joint.sum(axis=0); mix = joint/pw[None, :]
            mx = [sum(mix[a, j]*self.x[a] for a in range(self.modes)) for j in range(self.modes)]
            mP = [sum(mix[a, j]*(self.P[a]+np.outer(self.x[a]-mx[j], self.x[a]-mx[j])) for a in range(self.modes)) for j in range(self.modes)]
        else:
            pw = np.ones(1); mx, mP = self.x, self.P
        self.Yabs += dy
        ll = np.zeros(self.modes); xs, Ps = [], []
        for j in range(self.modes):
            A = self.A.copy()
            if self.form == 'incr':
                A[1, 0] = dy
            kp = max(mx[j][0], KLO/UNIT)
            q = self.q_abs*self.q_abs*TC/1e6 if self.q_abs is not None else self.rel[j]*kp*kp*TC
            Q = self.Qd.copy(); Q[0, 0] += q
            if self.started:
                x = A@mx[j]; P = A@mP[j]@A.T+Q
            else:
                x = mx[j].copy(); P = mP[j].copy()
            h = np.zeros(self.n); h[self.id] = 1.
            if self.form == 'incr':
                h[1] = 1.
            else:
                h[0] = self.Yabs
            e = z-h@x; Ph = P@h; S = h@Ph+self.sigma*self.sigma
            K = Ph/S
            B = np.eye(self.n)-np.outer(K, h)
            P2 = B@P@B.T+self.sigma*self.sigma*np.outer(K, K)
            xs.append(x+K*e); Ps.append(.5*(P2+P2.T))
            ll[j] = -.5*(math.log(2*math.pi*S)+e*e/S)
        self.x, self.P = xs, Ps
        self.started = True
        if self.modes > 1:
            lw = np.log(pw)+ll; lw -= lw.max(); self.w = np.exp(lw); self.w /= self.w.sum()
        return float(np.clip(UNIT*sum(self.w[j]*self.x[j][0] for j in range(self.modes)), KLO, KHI))


LADDER = dict(
    K1_row_o2=dict(form='row', d='o2', sigma=SIG_NOM, q_abs=300.),
    K2_row_d4=dict(form='row', d='d4', sigma=SIG_NOM, q_abs=300.),
    K3_incr_o2=dict(form='incr', d='o2', sigma=SIG_NOM, q_abs=300.),
    K4_incr_d4=dict(form='incr', d='d4', sigma=SIG_NOM, q_abs=300.),
    K5_true_sigma=dict(form='incr', d='d4', sigma=SIG_TRUE, q_abs=300.),
    K6_rel_q=dict(form='incr', d='d4', sigma=SIG_TRUE, q_rel=.02),
    T2_lo=dict(form='incr', d='d4', sigma=SIG_TRUE, q_rel=(.002, .5)),
    T2_mid=dict(form='incr', d='d4', sigma=SIG_TRUE, q_rel=(.002, 2.)),
    T2_hi=dict(form='incr', d='d4', sigma=SIG_TRUE, q_rel=(.02, 2.)),
)


def run_case(args):
    path, ticks_dir, configs = args
    data = H.load(path); meta = data['meta']; log, masks = H.masks_for(data)
    dY, z = ticks(data); n = len(z)
    out = {}; cost = {}
    for name, cfg in configs.items():
        f = TickKalman(**cfg); est = np.full(n, np.nan); last = 30.; tic = time.perf_counter()
        for j in range(n):
            if not np.isfinite(z[j]):
                if f.started:
                    last = f.x[0][0]; f.reset(last)
                continue
            est[j] = f.step(dY[j], z[j])
        cost[name] = (time.perf_counter()-tic)/max(np.isfinite(est).sum(), 1)*1e6
        out[name] = est
    with np.load(Path(ticks_dir)/f"Huber-D4_{meta['index']:02d}.npz") as t:
        out['raw window'] = t['D4_nogate__raw'].astype(float)[:n]; out['window + IMM'] = t['D4_imm__raw'].astype(float)[:n]
        out['Kalman baseline'] = t['Kalman__raw'].astype(float)[:n] if 'Kalman__raw' in t.files else None
    if out['Kalman baseline'] is None:
        with np.load(HERE/f"results/id_replay_ticks/Huber-D4_{meta['index']:02d}.npz") as t:
            out['Kalman baseline'] = t['Kalman__raw'].astype(float)[:n]
    # Same transition scoring as table 3.
    sc = dict(event_info=meta['event_info'])
    t = log[:, 0]
    tr = {}
    for name, est in out.items():
        e = np.where(np.isfinite(est) & (est > 0), est, np.nan)
        scored = R.score(log[:n], np.c_[t[:n], e, e], sc, names=('raw', 'k_good'), keep_after_loss=True)
        tr[name] = [dict(direction=q['direction'], size=q['size_class'], overlaps=q['overlaps_event'], change=q['channels']['raw']['change_mae'], under=q['channels']['raw']['raw_under_magnitude']) for q in scored['transitions']]
    ell = np.log(np.where(log[:n, 26] > 0, log[:n, 26], np.nan))
    seg = {s: {name: np.log(np.where(est > 0, est, np.nan))[masks[s][:n]]-ell[masks[s][:n]] for name, est in out.items()} for s in SEGS}
    return dict(index=meta['index'], group=meta['group'], seg=seg, tr=tr, cost=cost, names=list(out))


def rsd(x):
    q = np.percentile(x, [25, 75]); return (q[1]-q[0])/1.349


def main():
    ticks_dir = Path(sys.argv[1])/'24'
    only = LADDER_CASES if len(sys.argv) < 3 or sys.argv[2] != 'all' else None
    files = [p for p in sorted(H.WIN.glob('case_*.npz')) if only is None or int(p.stem[-2:]) in only]
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(run_case, [(p, str(ticks_dir), LADDER) for p in files], chunksize=1))
    names = res[0]['names']
    report = dict(cases=[r['index'] for r in res], segments={}, transitions={}, cost={k: float(np.median([r['cost'][k] for r in res])) for k in LADDER})
    print('us per tick', {k: round(v) for k, v in report['cost'].items()})
    for s in SEGS:
        print(f'--- {s}: |e| med/P90, bias, within-case sd median, below .85')
        report['segments'][s] = {}
        for nm in names:
            arrs = [r['seg'][s][nm] for r in res if len(r['seg'][s][nm])]
            e = np.concatenate(arrs); e = e[np.isfinite(e)]
            within = np.median([rsd(a[np.isfinite(a)]) for a in arrs if np.isfinite(a).sum() > 10])
            row = dict(med=float(np.median(np.abs(e))), p90=float(np.percentile(np.abs(e), 90)), bias=float(np.median(e)), within=float(within), below=float(np.mean(e < math.log(.85))), n=len(e))
            report['segments'][s][nm] = row
            print(f"  {nm:16s} {row['med']:.3f}/{row['p90']:.3f}  {row['bias']:+.3f}  {row['within']:.3f}  {row['below']*100:4.0f}%")
    print('--- transitions, no event, |dlnk|>=.3 then .15-.3: change error median, under median, under worst; up then down')
    for size in ('large', 'small'):
        report['transitions'][size] = {}
        for nm in names:
            line = []
            cell = {}
            for d in ('rising', 'falling'):
                ev = [q for r in res for q in r['tr'][nm] if q['size'] == size and q['direction'] == d and not q['overlaps']]
                ch = [q['change'] for q in ev if q['change'] is not None]; un = [q['under'] for q in ev if q['under'] is not None]
                cell[d] = dict(n=len(ev), change=float(np.median(ch)) if ch else None, under=float(np.median(un)) if un else None, worst=float(np.max(un)) if un else None)
                line.append(f"n{len(ev)} {cell[d]['change']:.3f} {cell[d]['under']:.3f} {cell[d]['worst']:.3f}" if ch and un else f'n{len(ev)} -')
            report['transitions'][size][nm] = cell
            print(f'  {size:5s} {nm:16s} '+' | '.join(line))
    tag = 'all' if only is None else 'ladder'
    (HERE/f'results/memory_tick_kalman_{tag}.json').write_text(json.dumps(report, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
