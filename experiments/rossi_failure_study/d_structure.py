"""What the unknown disturbance d looks like in closed-loop contact, measured from truth columns, before any model of it.

Definitions on the logged 8 ms ticks, inside contact and after 1 s, with k the true gain at the probe frequency:
  s_i = sum_j (dF_j/k_j - dy_j)   motion the true gain cannot explain, in metres; d_i = k_i s_i in newtons
Checks:
  in-band   rms of d in 16-26 Hz relative to the probe force, and the signed error it causes, <d_bp, y_bp>/(k <y_bp, y_bp>)
  scaling   does that relative contamination grow with k, as the body-inertia estimate k f/(m w^2) predicts?
  innovation  w = D^n s for n = 1..4: heavy tails, how much energy the outliers carry, how long outlier bursts last,
            whiteness of the core, and core size against the force-noise floor
The plant generates motion by body dynamics, so a good fit here supports the model only within this simulation.
"""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.special import comb

import sweeps as S

HERE = Path(__file__).resolve().parent
BASE = 410000
N_SCENE = 5
TC = .008
NOISE = .04      # design force-noise scale used by the estimators [N]
BAND = signal.butter(4, [16., 26.], btype='band', fs=1./TC, output='sos')


def segments(mask, min_len=64):
    d = np.diff(np.r_[0, mask.astype(int), 0]); a = np.flatnonzero(d == 1); b = np.flatnonzero(d == -1)
    return [(i, j) for i, j in zip(a, b) if j-i >= min_len]


def burst_lengths(flag, gap=2):
    idx = np.flatnonzero(flag)
    if not len(idx):
        return []
    return [int(g[-1]-g[0]+1) for g in np.split(idx, np.flatnonzero(np.diff(idx) > gap)+1)]


def innovation(s, k_med, n):
    w = np.diff(s, n=n)
    core = 1.4826*np.median(np.abs(w-np.median(w)))
    if core <= 0:
        return None
    out = np.abs(w) > 5*core
    inner = w[~out]; inner = inner-inner.mean()
    ac = float(np.sum(inner[1:]*inner[:-1])/max(np.sum(inner*inner), 1e-30))
    floor = NOISE/k_med*math.sqrt(comb(2*n, n))
    return dict(core=float(core), kurt=float(np.mean((w-w.mean())**4)/max(np.var(w)**2, 1e-60)), frac_out=float(out.mean()),
                energy_out=float(np.sum(w[out]**2)/max(np.sum(w*w), 1e-60)), bursts=burst_lengths(out),
                lag1=ac, core_over_noise=float(core/floor), n_samples=int(len(w)))


def job(a):
    c, t, i, controls, est = a
    sc = S.cell_scene(c, t, i, BASE); seed = BASE+sum(ord(x) for x in c+t)*100+i
    log, _ = S._run(sc, 'Projection', controls, est, seed)
    if log is None:
        return dict(cond=c, tissue=t, i=i, failed=True)
    tt, y, F, k = log[:, 0], log[:, 2], log[:, 3], log[:, 8]
    keep = (tt >= 1.) & (F > .5) & (k > 0)
    ev = np.zeros(len(tt), bool)
    for _, a0, a1 in sc.get('events', []):
        ev |= (tt >= a0) & (tt <= a1)
    out = dict(cond=c, tissue=t, i=i, failed=False, k_med=float(np.median(k[keep])) if keep.any() else np.nan,
               m=float(sc['m']), inband=[], innov={str(n): [] for n in (1, 2, 3, 4)})
    for j0, j1 in segments(keep):
        dF, dy = np.diff(F[j0:j1]), np.diff(y[j0:j1]); kk = k[j0:j1-1]
        s = np.r_[0., np.cumsum(dF/kk-dy)]
        d = s*k[j0:j1]; yy = y[j0:j1]; ks = k[j0:j1]
        dbp = signal.sosfiltfilt(BAND, d); ybp = signal.sosfiltfilt(BAND, yy)
        evseg = ev[j0:j1]
        for name, m in (('quiet', ~evseg), ('event', evseg)):
            if m.sum() < 32:
                continue
            kmean = float(np.mean(ks[m])); py = float(np.sum(ybp[m]**2))
            basis = np.c_[ybp, np.imag(signal.hilbert(ybp))][m]
            fit = np.linalg.lstsq(basis, dbp[m], rcond=None)[0]
            out['inband'].append(dict(state=name, k=kmean, ticks=int(m.sum()),
                                      coherent=float(np.sum((basis@fit)**2)/max(np.sum(dbp[m]**2), 1e-30)),
                                      rel_free=float(np.sqrt(np.mean((dbp[m]-basis@fit)**2))/max(kmean*np.sqrt(np.mean(ybp[m]**2)), 1e-12)),
                                      in_phase=float(fit[0]/kmean), quadrature=float(fit[1]/kmean),
                                      rel=float(np.sqrt(np.mean(dbp[m]**2))/max(kmean*np.sqrt(np.mean(ybp[m]**2)), 1e-12)),
                                      bias=float(np.sum(dbp[m]*ybp[m])/max(kmean*py, 1e-30)),
                                      d_rms=float(np.std(d[m]))))
        for n in (1, 2, 3, 4):
            r = innovation(s, float(np.median(ks)), n)
            if r is not None:
                out['innov'][str(n)].append(r)
    return out


def table(rows):
    res = {}
    for group, pick in (('quiet conds', lambda r: r['cond'] in ('rest', 'breath')), ('slow shift', lambda r: r['cond'] == 'shift_slow'),
                        ('fast shift', lambda r: r['cond'] == 'shift_fast'), ('flinch+cough', lambda r: r['cond'] in ('flinch', 'cough')),
                        ('tremor', lambda r: r['cond'] == 'tremor'), ('rib', lambda r: r['cond'] == 'rib'), ('combo', lambda r: r['cond'] == 'combo')):
        rr = [r for r in rows if pick(r) and not r['failed']]
        res[group] = {}
        for n in ('1', '2', '3', '4'):
            x = [q for r in rr for q in r['innov'][n]]
            if not x:
                continue
            b = [v for q in x for v in q['bursts']]
            res[group][n] = dict(kurt=float(np.median([q['kurt'] for q in x])), frac_out=float(np.mean([q['frac_out'] for q in x])),
                                 energy_out=float(np.median([q['energy_out'] for q in x])), lag1=float(np.median([q['lag1'] for q in x])),
                                 core_over_noise=float(np.median([q['core_over_noise'] for q in x])),
                                 burst_p50=float(np.median(b)) if b else 0., burst_p90=float(np.percentile(b, 90)) if b else 0.)
    return res


def scaling(rows):
    res = {}
    pts = [(r['cond'], e) for r in rows if not r['failed'] for e in r['inband']]
    for state in ('quiet', 'event'):
        for t in (None,):
            x = [(e['k'], e['rel'], e['bias']) for _, e in pts if e['state'] == state and e['k'] > 0 and e['rel'] > 0]
            if len(x) < 5:
                continue
            k, rel, bias = np.array(x).T
            slope, icept = np.polyfit(np.log(k), np.log(rel), 1)
            bins = [(0, 4e3), (4e3, 1e4), (1e4, 2e4), (2e4, 1e9)]
            res[state] = dict(slope=float(slope), n=len(k), by_k=[dict(k_lo=a, k_hi=b, n=int(((k >= a) & (k < b)).sum()),
                              rel_p50=float(np.median(rel[(k >= a) & (k < b)])) if ((k >= a) & (k < b)).any() else None,
                              rel_p90=float(np.percentile(rel[(k >= a) & (k < b)], 90)) if ((k >= a) & (k < b)).any() else None,
                              bias_p50=float(np.median(bias[(k >= a) & (k < b)])) if ((k >= a) & (k < b)).any() else None) for a, b in bins])
    return res


def main():
    t = json.loads((HERE/'results'/'tuning.json').read_text(encoding='utf-8'))
    controls = t['public']['C_neutral']['best']['params']; est = t['estimators']['C_neutral']['Projection']['best']['params']
    args = [(c, ti, i, controls, est) for c, ti in S.CELLS for i in range(N_SCENE)]
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, args, chunksize=1))
    out = dict(innovation=table(rows), scaling=scaling(rows), n=len(rows), failed=sum(r['failed'] for r in rows))
    (HERE/'results'/'d_structure.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
    print('scenes', out['n'], 'failed', out['failed'])
    print('group n kurt frac_out energy_out lag1 core/noise burst_p50/p90')
    for g, v in out['innovation'].items():
        for n, q in v.items():
            print(f"{g:13s} D{n} {q['kurt']:9.1f} {q['frac_out']:.3f} {q['energy_out']:.3f} {q['lag1']:+.2f} {q['core_over_noise']:7.2f} {q['burst_p50']:.0f}/{q['burst_p90']:.0f}")
    pts = [e for r in rows if not r['failed'] for e in r['inband']]
    out['probe_locked'] = {}
    for state in ('quiet', 'event'):
        for a, b in ((0, 4e3), (4e3, 1e4), (1e4, 2e4), (2e4, 1e9)):
            x = [e for e in pts if e['state'] == state and a <= e['k'] < b]
            if x:
                q = {f: float(np.median([e[f] for e in x])) for f in ('coherent', 'rel', 'rel_free', 'in_phase', 'quadrature')}
                out['probe_locked'][f'{state} {a:.0f}-{b:.0f}'] = dict(q, n=len(x))
                print(f"locked {state} k {a:.0f}-{b:.0f} n{len(x)} coherent {q['coherent']:.2f} rel {q['rel']:.3f} free {q['rel_free']:.3f} in-phase {q['in_phase']:+.3f} quad {q['quadrature']:+.3f}")
    (HERE/'results'/'d_structure.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
    for s, v in out['scaling'].items():
        print(s, 'slope ln(rel)~ln(k)', round(v['slope'], 2), 'n', v['n'])
        for b in v['by_k']:
            if b['n']:
                print(f"   k {b['k_lo']:.0f}-{b['k_hi']:.0f} n{b['n']} rel p50 {b['rel_p50']:.3f} p90 {b['rel_p90']:.3f} signed bias p50 {b['bias_p50']:+.3f}")


if __name__ == '__main__':
    main()
