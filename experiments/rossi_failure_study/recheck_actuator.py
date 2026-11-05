"""After the actuator correction (held setpoint, two poles, inertia from the same motion): the assumption on the fast pole, and the two conclusions it could affect.

A  fast pole 1, 2, 4 ms: ideal system and all system factors, L
B  contact model: tuned L, H, W at 14 and 28 deg, settings as chosen in select4.py, new seeds
C  d handling: D4 against projection with the settings chosen in d_compare.py, its test scenes
"""
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import check_system3 as CS
import d_compare as D
import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
TUNED = {'L': dict(r_scale=1., q_maneuver=3.335), 'H': dict(r_scale=1., q_maneuver=3.335), 'W': dict(r_scale=1., q_maneuver=3.335)}
D_CHOSEN = {'D4': (1, 0, 0), 'P': (0, 0, 0)}     # d_compare.py: D4 32 ticks, chi .551, tracker r 1 q 3.335; P 16 ticks quintic, chi .551, same tracker


def _job(a):
    part = a[0]
    if part == 'A':
        _, t2, change, case = a
        P3.ACTUATOR_T2 = t2
        c, t, g, i = case
        sc = S3.cell_scene3(c, t, g, i, base=492000, overrides=dict(force=8., dynamics='neural_white', phase=14., **CS.CHANGES[change]))
        pol = M.build('L', sc); log = P3.episode(sc, pol, 950+i)
        m = M.metric(log, pol, sc); m.update(part=part, t2=t2, change=change, case=case)
        return m
    if part == 'B':
        _, name, cond, tissue, dyn, ph, i = a
        sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=470000+500, overrides=dict(slide='long', pose_error=0., rho=0., force=8., dynamics=dyn, phase=ph))
        pol = M.build(name, sc, TUNED[name]); log = P3.episode(sc, pol, 800+i)
        m = M.metric(log, pol, sc); m.update(part=part, name=name, cond=cond, tissue=tissue, dynamics=dyn, phase=ph, i=i)
        return m
    _, host, cond, tissue, dyn, ph, i = a
    sc = S3.cell_scene3(cond, tissue, 'nominal', i, base=D.BASE+500, overrides=dict(slide='long', force=8., dynamics=dyn, phase=ph))
    pol = D.make(host, D_CHOSEN[host], sc); log = P3.episode(sc, pol, 900+i)
    m = M.metric(log, pol, sc); m.update(part=part, host=host, cond=cond, tissue=tissue, dynamics=dyn, phase=ph, i=i)
    return m


def paired(rows, sel_a, sel_b, sel):
    key = lambda r: (r['cond'], r['tissue'], r['dynamics'], r['phase'], r['i'])
    A = {key(r): r for r in rows if sel_a(r) and sel(r)}; B = {key(r): r for r in rows if sel_b(r) and sel(r)}
    d = np.array([np.log1p(A[x]['J'])-np.log1p(B[x]['J']) for x in A if x in B])
    ci = np.percentile(d[np.random.default_rng(1).integers(0, len(d), (5000, len(d)))].mean(1), [2.5, 97.5])
    return f"{d.mean():+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}] n={len(d)}"


def main():
    jobs = [('A', t2, ch, c) for t2 in (.001, .002, .004) for ch in ('ideal system, no probe mass', 'all system factors') for c in CS.CASES]
    jobs += [('B', n, cond, t, dyn, ph, i) for n in ('L', 'H', 'W') for cond in ('breath', 'shift_fast') for t in ('relaxed', 'active') for dyn in S3.SHAPES for ph in (14., 28.) for i in range(3)]
    jobs += [('C', h, c, t, dyn, ph, i) for h in ('D4', 'P') for c, t in D.TEST for dyn in ('neural_white', 'viscous') for ph in (14., 28.) for i in range(3)]
    print(len(jobs), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    print('A. fast pole: L, axial RMSE median, accept')
    for t2 in (.001, .002, .004):
        for ch in ('ideal system, no probe mass', 'all system factors'):
            q = [r for r in rows if r['part'] == 'A' and r['t2'] == t2 and r['change'] == ch]
            print(f"  {t2*1e3:.0f} ms {ch:30s} L {np.mean([np.log1p(r['J']) for r in q]):.3f} rmse {np.median([r['rmse'] for r in q]):.3f} accept {np.mean([r['accept'] for r in q]):.2f} >15 {sum(r['peak'] > 15 for r in q)}")
    B = [r for r in rows if r['part'] == 'B']
    print('B. contact model, paired, negative is better')
    for ph in (14., 28.):
        for n in ('L', 'H', 'W'):
            print(f"  {n} {ph:.0f}: "+D.line([r for r in B if r['name'] == n and r['phase'] == ph]))
        for a, b in (('H', 'L'), ('W', 'L'), ('H', 'W')):
            print(f"  {a}-{b} {ph:.0f}: {paired(B, lambda r, a=a: r['name'] == a, lambda r, b=b: r['name'] == b, lambda r, ph=ph: r['phase'] == ph)}")
    C = [r for r in rows if r['part'] == 'C']
    print('C. d handling')
    for h in ('D4', 'P'):
        print(f"  {h}: "+D.line([r for r in C if r['host'] == h]))
    print(f"  D4-P all: {paired(C, lambda r: r['host'] == 'D4', lambda r: r['host'] == 'P', lambda r: True)}")
    for ph in (14., 28.):
        print(f"  D4-P {ph:.0f}: {paired(C, lambda r: r['host'] == 'D4', lambda r: r['host'] == 'P', lambda r, ph=ph: r['phase'] == ph)}")
    (HERE/'results'/'recheck_actuator.json').write_text(json.dumps(rows, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
