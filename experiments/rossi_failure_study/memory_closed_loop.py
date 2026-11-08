"""Closed-loop comparison of memory mechanisms inside the adopted Method. Development seeds select; fresh seeds confirm.

Stage dev: base 520000 (never the overlay test grid), scenes A (breath_rib), B (breath_fast), A+B (breath_rib_fast), with the probe mass present.
Usage: python -X utf8 memory_closed_loop.py dev
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import overlay_compare as O
import method3 as M
import memory_filters as MF
import physics_plant3 as P3

B = O.B
HERE = Path(__file__).resolve().parent
DEV_BASE = 520000
SEED = 700000

DEV = dict(
    window=('window', {}),
    fix_002=('kalman', dict(q_t=.002, q_s=0.)), fix_01=('kalman', dict(q_t=.01, q_s=0.)), fix_05=('kalman', dict(q_t=.05, q_s=0.)),
    dist_2e4=('kalman', dict(q_t=0., q_s=2e-4)), dist_1e3=('kalman', dict(q_t=0., q_s=1e-3)), dist_5e3=('kalman', dict(q_t=0., q_s=5e-3)),
    both=('kalman', dict(q_t=.002, q_s=1e-3)),
    lepski_2=('lepski', dict(c=2.)), lepski_3=('lepski', dict(c=3.)),
    forget=('forgetting', dict(lam_max=.995, lam_min=.8, gain=1.)),
    sparse=('sparse', dict(gamma0=1e-4, a=.5, b=1e-5)), sparse_d=('sparse', dict(gamma0=1e-4, a=.5, b=1e-5, q_s=5e-4)),
)


DEV3 = dict(
    window=('window', {}),
    fix_002=('kalman', dict(q_t=.002, q_s=0.)), fix_0005=('kalman', dict(q_t=.0005, q_s=0.)), fix_0001=('kalman', dict(q_t=.0001, q_s=0.)),
    dist_2e4=('kalman', dict(q_t=0., q_s=2e-4)), dist_5e5=('kalman', dict(q_t=0., q_s=5e-5)), dist_1e5=('kalman', dict(q_t=0., q_s=1e-5)),
    forget=('forgetting', dict()), forget_998=('forgetting', dict(lam_max=.998)), forget_d=('forgetting', dict(q_s=5e-5)), forget_998d=('forgetting', dict(lam_max=.998, q_s=5e-5)),
    sparse=('sparse', dict(gamma0=1e-4, a=.5, b=1e-5)), sparse_d5=('sparse', dict(gamma0=1e-4, a=.5, b=1e-5, q_s=5e-5)),
)
DEV7 = dict(
    window=('window', {}),
    forget=('forgetting', dict(lam_max=.998)),
    D4=('none', dict(method='D4')),
    D4_forget=('forgetting', dict(lam_max=.998, method='D4')),
    Rossi=('none', dict(method='Rossi-ID')),
    MHE=('none', dict(method='MHE')),
    Xing=('none', dict(method='Xing-DAHIC')),
)
DEV6 = dict(
    window=('window', {}),
    forget=('forgetting', dict(lam_max=.998)),
    forget_keep=('forgetting', dict(lam_max=.998, keep_sigma=True)),
    forget_keep_995=('forgetting', dict(lam_max=.995, keep_sigma=True)),
    kal_002_keep=('kalman', dict(q_t=.002, q_s=0., huber=0., keep_sigma=True)),
    lepski_keep=('lepski', dict(c=2., keep_sigma=True)),
    ard_keep=('huang_ard', dict(learn_rho=False, keep_sigma=True)),
    forget_keep_d=('forgetting', dict(lam_max=.998, q_s=5e-5, keep_sigma=True)),
)
DEV5 = dict(
    window=('window', {}),
    forget=('forgetting', dict(lam_max=.998)),
    forget_nogate=('forgetting', dict(lam_max=.998, raw_gate=False)),
    amp70=('forgetting', dict(lam_max=.998, raw_gate=False, a_scale=.7)),
    amp50=('forgetting', dict(lam_max=.998, raw_gate=False, a_scale=.5)),
    shape35=('forgetting', dict(lam_max=.998, shape_fc=35.)),
    shape50=('forgetting', dict(lam_max=.998, shape_fc=50.)),
    forget_d5=('forgetting', dict(lam_max=.998, q_s=5e-5)),
)
DEV4 = dict(
    window=('window', {}),
    kal_01=('kalman', dict(q_t=.01, q_s=0., huber=0.)), kal_05=('kalman', dict(q_t=.05, q_s=0., huber=0.)), kal_002=('kalman', dict(q_t=.002, q_s=0., huber=0.)),
    lepski_2=('lepski', dict(c=2.)),
    forget=('forgetting', dict(lam_max=.998)), forget_d=('forgetting', dict(lam_max=.998, q_s=1e-5)),
    tard_23=('threshold_ard', dict(q_t=.002, q_s=0., kappa_up=2., kappa_down=3.)), tard_33=('threshold_ard', dict(q_t=.002, q_s=0., kappa_up=3., kappa_down=3.)),
    tard_153=('threshold_ard', dict(q_t=.002, q_s=0., kappa_up=1.5, kappa_down=3.)), tard_23d=('threshold_ard', dict(q_t=0., q_s=1e-4, kappa_up=2., kappa_down=3.)),
    ard_fix=('huang_ard', dict(learn_rho=False)),
)


def scene_for(cell, base, variant):
    if variant == 'nominal':
        return B.scene(cell, base, 0.)
    c, t, dy, ph, i = cell
    rng = np.random.default_rng([base, i, int(ph), dy == 'viscous', 7])
    o = dict(slide='long', force=8., dynamics=dy, phase=ph, speed=float(rng.uniform(.004, .010)), cover=.007)
    return B.S3.cell_scene3(c, t, 'nominal', i, base=base, overrides=o)


DEV2 = dict(
    window=('window', {}),
    huang_ard=('huang_ard', dict()), huang_ard_d=('huang_ard', dict(q_s=1e-3)), huang_ard_fixrho=('huang_ard', dict(learn_rho=False)),
)


def job(a):
    name, kind, params, cell, base, mass, seed = a[:7]
    variant = a[7] if len(a) > 7 else 'nominal'
    O.setup()
    params = dict(params); method = params.pop('method', 'Huber-D4')
    cfg = O.chosen()[method]
    sc = scene_for(cell, base, variant)
    if not mass:
        sc['tool_mass'] = 0.; sc['tool_mass_hat'] = 0.
    pol = B.build(method, cfg, sc)
    if kind == 'dual':
        MF.attach_dual(pol, **params)
    elif kind != 'none':
        MF.attach(pol, kind, **params)
    log = P3.episode(sc, pol, seed)
    m = M.metric(log, pol, sc)
    c, t, dy, ph, i = cell
    return dict(variant=variant, name=name, scene=c, dynamics=dy, phase=ph, i=i, L=float(np.log1p(m['J'])), rmse=float(m['rmse']), peak=float(m['peak']), mass=mass)


def boot(d, groups, n=5000, seed=0):
    rng = np.random.default_rng(seed); keys = sorted(set(groups)); idx = {k: np.array([q for q, g in enumerate(groups) if g == k]) for k in keys}
    out = [np.mean(np.concatenate([d[idx[keys[p]]] for p in rng.choice(len(keys), len(keys))])) for _ in range(n)]
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def report(rows, ref='window'):
    key = lambda r: (r['scene'], r['dynamics'], r['phase'], r['i'], r['mass'])
    names = set(r['name'] for r in rows)
    ref = ref if ref in names else 'Huber-D4'
    base = {key(r): r for r in rows if r['name'] == ref or r['name'] == ref+' no mass'}
    names = list(dict.fromkeys(r['name'] for r in rows))
    for s in sorted(set(r['scene'] for r in rows)):
        print('--', s)
        for n in names:
            q = [r for r in rows if r['name'] == n and r['scene'] == s]
            if not q:
                continue
            d = np.array([r['L']-base[key(r)]['L'] for r in q]); ci = boot(d, [r['i'] for r in q], 2000)
            print(f"{n:14s} L {np.mean([r['L'] for r in q]):.3f}/{np.median([r['L'] for r in q]):.3f} rmse {np.median([r['rmse'] for r in q]):.3f} peak {max(r['peak'] for r in q):5.1f} >15N {sum(r['peak']>15 for r in q):2d} dL {d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]", flush=True)
    print('-- all')
    for n in names:
        q = [r for r in rows if r['name'] == n]
        d = np.array([r['L']-base[key(r)]['L'] for r in q]); ci = boot(d, [(r['scene'], r['i']) for r in q], 2000)
        print(f"{n:14s} L {np.mean([r['L'] for r in q]):.3f}/{np.median([r['L'] for r in q]):.3f} peak {max(r['peak'] for r in q):5.1f} >15N {sum(r['peak']>15 for r in q):2d} dL {d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]", flush=True)


def main():
    stage = sys.argv[1]
    start = time.perf_counter()
    if stage == 'dev':
        cells = O.cells()
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV.items() for q, c in enumerate(cells)]
    if stage == 'dev3':
        cells = O.cells()
        slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV3.items() for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, DEV_BASE, True, SEED+q, 'slow_thin') for n, (k, p) in DEV3.items() for q, c in enumerate(slow)]
    if stage == 'dev4':
        cells = O.cells()
        slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV4.items() for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, DEV_BASE, True, SEED+q, 'slow_thin') for n, (k, p) in DEV4.items() for q, c in enumerate(slow)]
    if stage == 'dev7':
        cells = O.cells()
        slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV7.items() for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, DEV_BASE, True, SEED+q, 'slow_thin') for n, (k, p) in DEV7.items() for q, c in enumerate(slow)]
    if stage == 'dev6':
        cells = O.cells()
        slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV6.items() for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, DEV_BASE, True, SEED+q, 'slow_thin') for n, (k, p) in DEV6.items() for q, c in enumerate(slow)]
    if stage == 'dev5':
        cells = O.cells()
        slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV5.items() for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, DEV_BASE, True, SEED+q, 'slow_thin') for n, (k, p) in DEV5.items() for q, c in enumerate(slow)]
    if stage == 'test':
        # never used for any choice: fresh object and sensor seeds. Selected configuration frozen before this run.
        base = 774500; seed = 880000
        sel = ('forgetting', dict(lam_max=.998))
        win = ('window', {})
        others = [(x, ('none', dict(method=x))) for x in ('D4', 'Rossi-ID', 'Kalman', 'MHE', 'Rossi-control', 'DO-MPAC', 'Xing-DAHIC')]
        methods = [('Huber-D4', win), ('Proposed', sel)]+others
        cells = O.cells(); slow = [c for c in cells if c[0] == 'breath_rib']
        jobs = [(n, k, p, c, base, True, seed+q) for n, (k, p) in methods for q, c in enumerate(cells)]
        jobs += [(n, k, p, c, base, True, seed+q, 'slow_thin') for n, (k, p) in methods for q, c in enumerate(slow)]
        for n, (k, p) in methods[:2]:
            jobs += [(n+' no mass', k, p, c, base, False, seed+q) for q, c in enumerate(cells)]
    if stage == 'test6':
        base = 774500; seed = 890000
        sel = ('forgetting', dict(lam_max=.998)); win = ('window', {})
        cells = B.cells()
        jobs = [(n, k, p, c, base, True, seed+q) for n, (k, p) in (('Huber-D4', win), ('Proposed', sel)) for q, c in enumerate(cells)]
    if stage == 'test6b':
        base = 774500; seed = 890000
        cells = B.cells()
        jobs = [(x, 'none', dict(method=x), c, base, True, seed+q) for x in ('D4', 'Rossi-ID', 'MHE') for q, c in enumerate(cells)]
    if stage == 'true_base':
        # the adopted Method with nothing attached. The earlier 'window' reference went through attach, which held its output for 6 ticks and doubled sigma.
        base = 774500; seed = 880000; cells = O.cells(); slow = [c for c in cells if c[0] == 'breath_rib']; tb = ('none', dict(method='Huber-D4'))
        jobs = [('Huber-true', tb[0], tb[1], c, base, True, seed+q) for q, c in enumerate(cells)]
        jobs += [('Huber-true', tb[0], tb[1], c, base, True, seed+q, 'slow_thin') for q, c in enumerate(slow)]
        jobs += [('Huber-true no mass', tb[0], tb[1], c, base, False, seed+q) for q, c in enumerate(cells)]
        jobs += [('Huber-true dev', tb[0], tb[1], c, DEV_BASE, True, SEED+q) for q, c in enumerate(cells)]
        jobs += [('Huber-true dev', tb[0], tb[1], c, DEV_BASE, True, SEED+q, 'slow_thin') for q, c in enumerate(slow)]
        jobs += [('Huber-true six', tb[0], tb[1], c, base, True, 890000+q) for q, c in enumerate(B.cells())]
    if stage == 'dev2':
        cells = O.cells()
        jobs = [(n, k, p, c, DEV_BASE, True, SEED+q) for n, (k, p) in DEV2.items() for q, c in enumerate(cells)]
    print('jobs', len(jobs), flush=True)
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = list(ex.map(job, jobs, chunksize=2))
    (HERE/'results'/f'memory_closed_loop_{stage}.json').write_text(json.dumps(rows, indent=1))
    if stage == 'true_base':
        print('seconds', round(time.perf_counter()-start), flush=True)
        return
    for v in sorted(set(r.get('variant', 'nominal') for r in rows)):
        print('==== variant', v)
        report([r for r in rows if r.get('variant', 'nominal') == v])
    print('seconds', round(time.perf_counter()-start), flush=True)


if __name__ == '__main__':
    main()
