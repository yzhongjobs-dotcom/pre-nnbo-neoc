"""Contact-model selection on the identified dynamic truth, in the full Method. Each question only where the models can differ.

dynamic  L, H, W, S: tuned with one budget on development seeds (white and viscous shapes, 21 deg); tested on new seeds, all three shapes, 14/21/28 deg.
         neural_grey is never seen in tuning.
static   L, Hz; sliding L, LF, AS, TR: small confirmations, the physics these address did not change with the new dynamic truth.
host     the projection host, only for the conclusion that decides the choice, default settings.
Usage: python select4.py dynamic | static | sliding | host NAME [NAME]
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
BASE = 470000
QUIET = dict(slide='long', pose_error=0., rho=0., force=8.)
CONFIGS = [dict(r_scale=r, q_maneuver=q) for r in (.2264, 1., 4., 16.) for q in (3.335, .8)]      # index 0 is the default
DYNAMIC = ('L', 'H', 'W', 'S')


def _job(a):
    q, name, c, cond, tissue, regime, i, base, o = a
    sc = S3.cell_scene3(cond, tissue, regime, i, base=base, overrides=o)
    pol = M.build(name, sc, CONFIGS[c] if c is not None else None)
    log = P3.episode(sc, pol, 800+i)
    m = M.metric(log, pol, sc)
    m.update(question=q, name=name, c=c, cond=cond, tissue=tissue, regime=regime, i=i, **{k: o.get(k) for k in ('dynamics', 'phase', 'probe', 'contact_law', 'force', 'mu')})
    return m


def run(jobs):
    print(len(jobs), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        return list(ex.map(_job, jobs, chunksize=1))


L = lambda q: float(np.mean([np.log1p(r['J']) for r in q]))


def summary(q):
    ok = [r for r in q if not r['failed']]
    return (f"L {L(q):.3f} rmse {np.median([r['rmse'] for r in ok]):.3f}/{np.mean([r['rmse'] for r in ok]):.3f} event {np.nanmedian([r['event_rmse'] for r in ok]):.3f} "
            f">15 {sum(r['peak'] > 15 for r in ok):2d} >36.5 {sum(r['peak'] > 36.5 for r in ok)} long {sum(r['long_loss'] for r in q)} fail {sum(r['failed'] for r in q)} "
            f"accept {np.mean([r['accept'] for r in ok]):.2f} gain err {np.nanmedian([abs(r['est_abs']) for r in ok]):.3f} call {np.median([r['call_p50'] for r in ok]):.0f}/{np.median([r['call_p99'] for r in ok]):.0f} us")


def paired(rows, a, b, sel):
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    A = {key(r): r for r in rows if r['name'] == a[0] and r['c'] == a[1] and sel(r)}
    B = {key(r): r for r in rows if r['name'] == b[0] and r['c'] == b[1] and sel(r)}
    d = np.array([np.log1p(A[x]['J'])-np.log1p(B[x]['J']) for x in A if x in B])
    ci = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
    return f"{d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}] n={len(d)}"


def dynamic():
    dev = [('dev', n, c, cond, 'relaxed', 'nominal', i, BASE, dict(QUIET, dynamics=dyn, phase=21.))
           for n in DYNAMIC for c in range(len(CONFIGS)) for cond in ('breath', 'shift_fast') for dyn in ('neural_white', 'viscous') for i in range(2)]
    rows = run(dev)
    chosen = {}
    print('development (21 deg, white and viscous): mean L per configuration')
    for n in DYNAMIC:
        scores = [L([r for r in rows if r['name'] == n and r['c'] == c]) for c in range(len(CONFIGS))]
        chosen[n] = int(np.argmin(scores))
        print(f"  {n}: "+' '.join(f"{s:.3f}" for s in scores)+f"  -> {CONFIGS[chosen[n]]}")
    test = [('test', n, c, cond, t, 'nominal', i, BASE+500, dict(QUIET, dynamics=dyn, phase=ph))
            for n, c in [(n, chosen[n]) for n in DYNAMIC]+[('L', 0)] for cond in ('breath', 'shift_fast') for t in ('relaxed', 'active')
            for dyn in S3.SHAPES for ph in (14., 21., 28.) for i in range(3)]
    rows2 = run(test)
    print('\nnew seeds, per phase (all shapes):')
    for ph in (14., 21., 28.):
        print(f' phase {ph:.0f}')
        for n, c in [('L', 0)]+[(n, chosen[n]) for n in DYNAMIC]:
            print(f"   {n} {'default' if c == 0 else 'tuned  '} "+summary([r for r in rows2 if r['name'] == n and r['c'] == c and r['phase'] == ph]))
    print('\npaired against tuned L, mean difference in L [95% CI]; negative is better')
    for n in ('H', 'W', 'S'):
        for ph in (14., 21., 28.):
            line = [paired(rows2, (n, chosen[n]), ('L', chosen['L']), lambda r, ph=ph: r['phase'] == ph)]
            line += [dyn[:7]+' '+paired(rows2, (n, chosen[n]), ('L', chosen['L']), lambda r, ph=ph, dyn=dyn: r['phase'] == ph and r['dynamics'] == dyn) for dyn in S3.SHAPES]
            print(f"  {n} {ph:.0f}: all {line[0]} | "+' | '.join(line[1:]))
    print('\nbetween the dynamic forms, all phases and shapes, and on the shape never used in tuning')
    for a, b in (('H', 'W'), ('H', 'S'), ('W', 'S')):
        print(f"  {a}-{b}: all {paired(rows2, (a, chosen[a]), (b, chosen[b]), lambda r: True)}  grey only {paired(rows2, (a, chosen[a]), (b, chosen[b]), lambda r: r['dynamics'] == 'neural_grey')}  28 deg {paired(rows2, (a, chosen[a]), (b, chosen[b]), lambda r: r['phase'] == 28.)}")
    (HERE/'results'/'select4_dynamic.json').write_text(json.dumps(dict(chosen={n: CONFIGS[c] for n, c in chosen.items()}, dev=rows, test=rows2), indent=1), encoding='utf-8')


def static():
    jobs = [('static', n, None, cond, 'relaxed', 'nominal', i, BASE+900, dict(QUIET, dynamics=dyn, phase=14., probe=probe, contact_law=law, force=force))
            for n in ('L', 'Hz') for probe, law in (('convex', 'foundation'), ('linear', 'halfspace')) for force in (4., 12.) for cond in ('shift_fast',)
            for dyn in ('neural_white',) for i in range(3)]
    rows = run(jobs)
    for n in ('L', 'Hz'):
        print(' ', n, summary([r for r in rows if r['name'] == n]))
    for probe in ('convex', 'linear'):
        for f in (4., 12.):
            q = lambda n: [r for r in rows if r['name'] == n and r['probe'] == probe and r['force'] == f]
            print(f"  {probe} {f:.0f} N: rmse L {np.median([r['rmse'] for r in q('L')]):.3f} Hz {np.median([r['rmse'] for r in q('Hz')]):.3f}")
    (HERE/'results'/'select4_static.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')


def sliding():
    jobs = [('sliding', n, None, 'rest', 'relaxed', regime, i, BASE+950, dict(dynamics='neural_white', phase=14., probe='linear', contact_law='halfspace', force=8., slide='long', rho=0., mu=.3, k_d=3., K_t=1.e4))
            for n in ('L', 'LF', 'AS', 'TR') for regime in ('stop_go', 'curved_untracked') for i in range(3)]
    rows = run(jobs)
    for regime in ('stop_go', 'curved_untracked'):
        for n in ('L', 'LF', 'AS', 'TR'):
            q = [r for r in rows if r['name'] == n and r['regime'] == regime]
            print(f"  {regime:16s} {n:3s} axial {np.median([r['rmse'] for r in q]):.3f} normal {np.median([r['normal_rmse'] for r in q]):.3f} call {np.median([r['call_p50'] for r in q]):.0f}")
    (HERE/'results'/'select4_sliding.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')


def host(*names):
    jobs = [('host', 'P5:'+n, None, cond, t, 'nominal', i, BASE+500, dict(QUIET, dynamics=dyn, phase=ph))
            for n in ('L',)+names for cond in ('breath', 'shift_fast') for t in ('relaxed', 'active') for dyn in S3.SHAPES for ph in (14., 28.) for i in range(3)]
    rows = run(jobs)
    for r in rows:
        r['c'] = 0
    for ph in (14., 28.):
        for n in ('L',)+names:
            print(f"  P5:{n} {ph:.0f} "+summary([r for r in rows if r['name'] == 'P5:'+n and r['phase'] == ph]))
        for n in names:
            print(f"  P5:{n} - P5:L {ph:.0f}: {paired(rows, ('P5:'+n, 0), ('P5:L', 0), lambda r, ph=ph: r['phase'] == ph)}")
        if len(names) > 1:
            print(f"  P5:{names[0]} - P5:{names[1]} {ph:.0f}: {paired(rows, ('P5:'+names[0], 0), ('P5:'+names[1], 0), lambda r, ph=ph: r['phase'] == ph)}")
    (HERE/'results'/'select4_host.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')


if __name__ == '__main__':
    {'dynamic': dynamic, 'static': static, 'sliding': sliding}.get(sys.argv[1], lambda: host(*sys.argv[2:]))()
