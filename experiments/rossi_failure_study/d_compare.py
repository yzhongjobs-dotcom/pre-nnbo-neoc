"""D⁴ prior against polynomial projection, both with the selected contact model H, in the full Method. Same tuning budget for each host.

Each host gets three settings of its own structure, two gate scales and four tracker settings (24), chosen on development scenes;
then both are tested on new seeds across the kinds of d the scan produces: slow motion, sudden events, stiffness change while sliding.
Logs of the test runs are kept for module replay (diagnostic truth only for scoring).
"""
import json
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

import contact_candidates as C
import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
BASE = 520000
GATE = dict(sigma_scale=.796)
STRUCT = {'D4': [dict(nw=24, order=4, fc=14.), dict(nw=32, order=4, fc=14.), dict(nw=24, order=4, fc=10.)],
          'P': [dict(nw=16, degree=5), dict(nw=24, degree=5), dict(nw=16, degree=3)]}
CHI = (.551, 1.)
TRACK = [dict(r_scale=r, q_maneuver=q) for r in (1., 4.) for q in (3.335, .8)]
CONFIGS = list(product(range(3), range(2), range(4)))
DEV = [(c, t) for c, t in (('breath', 'relaxed'), ('shift_fast', 'relaxed'), ('cough', 'active'), ('rib', 'rib'))]
TEST = [('breath', 'relaxed'), ('shift_fast', 'relaxed'), ('flinch', 'active'), ('cough', 'relaxed'), ('rib', 'rib'), ('combo', 'mixed')]


def make(host, cfg, sc):
    s, g, t = cfg
    if host == 'D4':
        C.RIDGE = dict(STRUCT['D4'][s], chi_scale=CHI[g], **GATE); name = 'H'
    else:
        C.P5 = dict(STRUCT['P'][s], chi_scale=CHI[g], **GATE); name = 'P5:H'
    return M.build(name, sc, TRACK[t])


def _job(a):
    host, cfg, cond, tissue, dyn, phase, i, base, keep = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=base, overrides=dict(slide='long', force=8., dynamics=dyn, phase=phase))
    pol = make(host, cfg, sc)
    log = P3.episode(sc, pol, 900+i)
    m = M.metric(log, pol, sc)
    m.update(host=host, cfg=list(cfg), cond=cond, tissue=tissue, dynamics=dyn, phase=phase, i=i)
    if keep and log is not None:
        m['log'] = log[:, [0, 1, 2, 3, 5, 6, 11, 18, 20, 21, 22]].astype(np.float32)
    return m


def run(jobs):
    print(len(jobs), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        return list(ex.map(_job, jobs, chunksize=1))


L = lambda q: float(np.mean([np.log1p(r['J']) for r in q]))


def line(q):
    ok = [r for r in q if not r['failed']]
    return (f"L {L(q):.3f} rmse {np.median([r['rmse'] for r in ok]):.3f}/{np.mean([r['rmse'] for r in ok]):.3f} event {np.nanmedian([r['event_rmse'] for r in ok]):.3f} "
            f">15 {sum(r['peak'] > 15 for r in ok)} >36.5 {sum(r['peak'] > 36.5 for r in ok)} long {sum(r['long_loss'] for r in q)} fail {sum(r['failed'] for r in q)} "
            f"accept {np.mean([r['accept'] for r in ok]):.2f} gain err {np.nanmedian([abs(r['est_abs']) for r in ok]):.3f} call {np.median([r['call_p50'] for r in ok]):.0f}/{np.median([r['call_p99'] for r in ok]):.0f} us")


def main():
    dev = [(h, cfg, c, t, 'neural_white', 21., i, BASE, False) for h in ('D4', 'P') for cfg in CONFIGS for c, t in DEV for i in range(2)]
    rows = run(dev)
    chosen = {}
    for h in ('D4', 'P'):
        scores = {cfg: L([r for r in rows if r['host'] == h and tuple(r['cfg']) == cfg]) for cfg in CONFIGS}
        chosen[h] = min(scores, key=scores.get)
        default = (0, 0, 0)
        print(f"{h}: default {scores[default]:.3f}, best {scores[chosen[h]]:.3f} with structure {STRUCT[h][chosen[h][0]]}, chi {CHI[chosen[h][1]]}, tracker {TRACK[chosen[h][2]]}")
        print('   by structure:', [round(min(v for k, v in scores.items() if k[0] == s), 3) for s in range(3)], 'by chi:', [round(min(v for k, v in scores.items() if k[1] == g), 3) for g in range(2)])
    test = [(h, chosen[h], c, t, dyn, ph, i, BASE+500, True) for h in ('D4', 'P') for c, t in TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    rows2 = run(test)
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    A = {key(r): r for r in rows2 if r['host'] == 'D4'}; B = {key(r): r for r in rows2 if r['host'] == 'P'}
    def diff(sel):
        d = np.array([np.log1p(A[x]['J'])-np.log1p(B[x]['J']) for x in A if sel(A[x])])
        ci = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
        return f"{d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}] n={len(d)}"
    print('\ntest, new seeds')
    for h in ('D4', 'P'):
        print(' ', h, line([r for r in rows2 if r['host'] == h]))
    print('  D4 - P, all:', diff(lambda r: True))
    for c, t in TEST:
        print(f"  {c:10s} {t:7s} D4 {np.median([r['rmse'] for r in rows2 if r['host'] == 'D4' and r['cond'] == c]):.3f} P {np.median([r['rmse'] for r in rows2 if r['host'] == 'P' and r['cond'] == c]):.3f}  D4-P {diff(lambda r, c=c: r['cond'] == c)}")
    for ph in (14., 28.):
        print(f"  phase {ph:.0f}: D4-P {diff(lambda r, ph=ph: r['phase'] == ph)}")
    (HERE/'results'/'d_compare.json').write_text(json.dumps(dict(chosen={h: dict(structure=STRUCT[h][c[0]], chi=CHI[c[1]], tracker=TRACK[c[2]]) for h, c in chosen.items()},
                                                               dev=rows, test=[{k: v for k, v in r.items() if k != 'log'} for r in rows2]), indent=1, default=float), encoding='utf-8')
    Path(sys.argv[1] if len(sys.argv) > 1 else HERE/'results'/'d_compare_logs.pkl').write_bytes(pickle.dumps(dict(chosen=chosen, rows=rows2)))


if __name__ == '__main__':
    main()
