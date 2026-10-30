"""Measurement-noise model of the force-derived motion measurement (paper draft §三.5, §五.6).

The earlier tracker set R = r (sigma_F/k)^2, i.e. force noise only; on stiff tissue R shrinks to ~1 um^2 and every
IMM mode becomes fast, which tracker_dilemma.py showed (IMM ~ high-bandwidth Kalman in quiet stiff segments) and
which ff_loop.py identifies as the stiffness range where the feed-forward loop tolerates the least gain. The
refinement adds a displacement-domain floor: R = r (sigma_F/k)^2 + sigma_y^2 (encoder interpolation, the (y-s) dln k
term and contact-model error do not shrink with stiffness).
Protocol (declared before running): tracker parameters (q_smooth, q_maneuver, p_stay, r_scale, sigma_y) are tuned
on the development split with the equal 61-evaluation budget, all other parameters fixed at the selected values; the
default candidate is the current IMM (sigma_y = 0). On the selection split (seeds 80000+) the refinement is adopted
only if its mean log(1+J) is lower than the current IMM's by more than 5% (it has one more parameter). Only then is
it evaluated on the confirmation cohort.
Run: python -X utf8 tracker_refine.py -> results_v2/tracker_refine.json, tuned_v2/Proposed-IMMf.json
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import scenes_v2 as SV
import variants

H = common.HERE
RES = os.path.join(H, 'results_v2')
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
KEYS = ('q_smooth', 'q_maneuver', 'p_stay', 'r_scale', 'r_floor')
B = dict(variants.BOUNDS, r_floor=(1e-7, 1e-4))


def build(params):
    return variants.build('Proposed-IMM', dict(params))


def _ep(a):
    params, (fam, sid, p, T) = a
    try:
        log, _ = common.run(build(params), p, T, sid)
        m = common.task_metrics(log, 20.8)
        J = m['J'] if np.isfinite(m['J']) else 1e3
    except Exception as ex:
        return dict(family=fam, sid=sid, L=math.log1p(1e3), rmse=np.nan, peak=np.nan, loss=np.nan, ok=False, err=str(ex))
    return dict(family=fam, sid=sid, L=math.log1p(min(J, 1e3)), rmse=m['rmse'], peak=m['peak'], loss=m['loss'],
                ok=not m['diverged'])


def score(pool, cands, sc):
    res = pool.map(_ep, [(c, s) for c in cands for s in sc], chunksize=1)
    return [(float(np.mean([x['L'] for x in res[i*len(sc):(i+1)*len(sc)]])), c) for i, c in enumerate(cands)]


def draw(r):
    q = dict(BASE)
    for k in KEYS:
        lo, hi = B[k]
        q[k] = float(math.exp(r.uniform(math.log(lo), math.log(hi))))
    return q


def perturb(c, r):
    q = dict(c)
    for k in KEYS:
        lo, hi = B[k]
        v = q[k] if q[k] > 0 else 1e-6
        q[k] = float(np.clip(v*math.exp(r.normal(0, .35)), lo, hi))
    return q


def boot(d, n=4000, seed=0):
    r = np.random.default_rng(seed); d = np.asarray(d)
    b = [d[r.integers(0, len(d), len(d))].mean() for _ in range(n)]
    return float(d.mean()), [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


if __name__ == '__main__':
    cur = dict(BASE, r_floor=0.)
    out = {}
    with Pool(16) as pool:
        dev = SV.scenes('dev'); r = np.random.default_rng(778)
        rows = score(pool, [cur]+[draw(r) for _ in range(40)], dev)
        best = min(rows, key=lambda x: x[0])
        for _ in range(4):
            rows += score(pool, [perturb(best[1], r) for _ in range(5)], dev)
            best = min(rows, key=lambda x: x[0])
        new = best[1]
        out['dev'] = dict(current=rows[0][0], refined=best[0], n_eval=len(rows), params={k: new[k] for k in KEYS})
        print('dev mean log(1+J): current %.3f refined %.3f' % (rows[0][0], best[0]), {k: new[k] for k in KEYS}, flush=True)
        json.dump(dict(best=dict(params=new, L=best[0]), n_eval=len(rows)),
                  open(os.path.join(H, 'tuned_v2', 'Proposed-IMMf.json'), 'w'), indent=1)
        sel = SV.scenes('select')
        Ls = score(pool, [cur, new], sel)
        adopt = Ls[1][0] < 0.95*Ls[0][0]
        out['select'] = dict(current=Ls[0][0], refined=Ls[1][0], adopt=bool(adopt))
        print('select mean log(1+J): current %.3f refined %.3f -> adopt %s' % (Ls[0][0], Ls[1][0], adopt), flush=True)
        if adopt:
            conf = SV.scenes('confirm')
            res = pool.map(_ep, [(new, s) for s in conf], chunksize=1)
            out['confirm_rows'] = res
            ref = json.load(open(os.path.join(RES, 'confirm', 'main.json')))
            summ = {}
            for fam in SV.FAMILIES:
                R = {x['sid']: x for x in res if x['family'] == fam}
                ok = [x for x in R.values() if x['ok']]
                e = dict(n_fail=len(R)-len(ok), rmse=float(np.mean([x['rmse'] for x in ok])),
                         peak=float(max(x['peak'] for x in ok)), loss=float(np.mean([x['loss'] for x in ok])))
                line = '%-9s fail %d rmse %.3f peak %.1f loss %.4f |' % (fam, e['n_fail'], e['rmse'], e['peak'], e['loss'])
                for m in ('Proposed', 'RossiType', 'WijType', 'DOMPAC', 'RossiDOMPAC'):
                    C = {x['sid']: x for x in ref if x['method'] == m and x['family'] == fam and x['ok']}
                    d = [C[s]['rmse']-R[s]['rmse'] for s in C if s in R and R[s]['ok']]
                    mu, ci = boot(d)
                    e['vs_'+m] = dict(mean=mu, ci=ci, wins=int(sum(x > 0 for x in d)), n=len(d))
                    line += ' %s %+.3f[%+.3f,%+.3f] %d/%d' % (m[:5], mu, *ci, e['vs_'+m]['wins'], len(d))
                summ[fam] = e
                print(line, flush=True)
            out['confirm'] = summ
    json.dump(out, open(os.path.join(RES, 'tracker_refine.json'), 'w'), indent=1)
