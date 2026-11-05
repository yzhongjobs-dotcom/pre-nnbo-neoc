"""Does the Method read the new plant correctly? Signal checks, not a comparison: force regulation per setpoint, gain against truth, phase, real-time cost."""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import scenes3 as S3

HERE = Path(__file__).resolve().parent
QUIET = dict(pose_error=0., mu=.1, rho=0., accel=.1, speed=.015, k_d=2., alpha_scale=1., tau_scale=1., kr=1.e5, K_t=2.e4, surface_amp=.0005)


def _job(a):
    name, tissue, cond, regime, o = a
    sc = S3.cell_scene3(cond, tissue, regime, 0, base=440000, overrides=o)
    log, pol, _ = M.run(name, sc, seed=1)
    m = M.metric(log, pol, sc)
    m.update(name=name, tissue=tissue, cond=cond, regime=regime, force=o.get('force'), law=o.get('contact_law'), probe=o.get('probe'), phi_hc=o.get('phi_hc'))
    return m


def main():
    jobs = []
    for name in ('L', 'P5:L', 'H', 'P5:H', 'W', 'P5:W'):
        for tissue in ('relaxed', 'active'):
            for force in (4., 8., 12.):
                for law in ('foundation', 'halfspace'):
                    for probe in ('convex', 'linear'):
                        for phi in (5., 30.):
                            jobs.append((name, tissue, 'rest', 'nominal', dict(QUIET, force=force, contact_law=law, probe=probe, slide='long', phi_hc=phi, field_sigma=.1)))
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    (HERE/'results'/'method3_checks.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    deg = 180/math.pi
    print('name  tissue  force  n | fails loss | mean Fn err, axial err | RMSE median | est/Re, est/|K| | phi hat, plant | accept | call p50 p99 us')
    for name in ('L', 'P5:L', 'H', 'P5:H', 'W', 'P5:W'):
        for tissue in ('relaxed', 'active'):
            for force in (4., 8., 12.):
                q = [r for r in rows if r['name'] == name and r['tissue'] == tissue and r['force'] == force]
                ok = [r for r in q if not r['failed']]
                med = lambda k: float(np.nanmedian([r[k] for r in ok])) if ok else float('nan')
                print(f"{name:5s} {tissue:7s} {force:3.0f} {len(q):3d} | {sum(r['failed'] or r['peak'] > 50 for r in q)} {sum(r['long_loss'] for r in q)} | {med('mean_err'):+.3f} {med('mean_normal_err'):+.3f} | {med('rmse'):.3f} | {math.exp(med('est_re')):.3f} {math.exp(med('est_abs')):.3f} | {med('phi')*deg:5.1f} {med('phi_true')*deg:5.1f} | {med('accept'):.2f} | {med('call_p50'):.0f} {med('call_p99'):.0f} | k {med('k_true')/1e3:.2f} d {med('depth')*1e3:.1f}")
    print('\nby law, probe, phi: RMSE median and fails, all names and tissues')
    for name in ('L', 'P5:L', 'H', 'P5:H'):
        for law in ('foundation', 'halfspace'):
            for probe in ('convex', 'linear'):
                for phi in (5., 30.):
                    q = [r for r in rows if r['name'] == name and r['law'] == law and r['probe'] == probe and r['phi_hc'] == phi]
                    print(f"{name:5s} {law:10s} {probe:6s} phi {phi:2.0f}: rmse {np.nanmedian([r['rmse'] for r in q]):.3f} max {np.nanmax([r['rmse'] for r in q]):.3f} osc {np.nanmedian([r['osc'] for r in q]):.3f} J {np.nanmedian([r['J'] for r in q]):.3f} fail {sum(r['failed'] or r['peak'] > 50 for r in q)}")


if __name__ == '__main__':
    main()
