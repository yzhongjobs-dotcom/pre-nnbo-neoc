"""Contact-model selection: each question only in the conditions where the models can differ, in the full Method. Other host only for what could change the choice.

dynamic   does the phase of the probe's own force response matter, whichever mechanism causes it?  L, H, W, S
static    does force-dependent stiffness matter?                                                  L, Hz
sliding   does friction reach the axial force enough to need a tangential model and sensor?       L, LF, AS, TR
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import scenes3 as S3

HERE = Path(__file__).resolve().parent
BASE = 460000
QUIET = dict(slide='long', pose_error=0., rho=0.)


def plan():
    jobs = []
    for dyn in ('hc', 'power', 'kb'):
        for phase in (11., 23.):
            for cond in ('rest', 'breath'):
                for probe, law in (('linear', 'halfspace'), ('convex', 'foundation')):
                    for name in ('L', 'H', 'W', 'S'):
                        for i in range(3):
                            jobs.append(('dynamic', name, cond, 'relaxed', 'nominal', i, dict(QUIET, dynamics=dyn, phase=phase, probe=probe, contact_law=law, force=8.)))
    for probe, law in (('convex', 'foundation'), ('linear', 'halfspace')):
        for force in (4., 12.):
            for cond in ('shift_fast', 'cough'):
                for name in ('L', 'Hz'):
                    for i in range(3):
                        jobs.append(('static', name, cond, 'relaxed', 'nominal', i, dict(QUIET, dynamics='power', phase=11., probe=probe, contact_law=law, force=force)))
    for regime in ('stop_go', 'reversal', 'curved_untracked'):
        for mu in (.1, .45):
            for name in ('L', 'LF', 'AS', 'TR'):
                for i in range(3):
                    jobs.append(('sliding', name, 'rest', 'relaxed', regime, i, dict(dynamics='power', phase=11., probe='linear', contact_law='halfspace', force=8., slide='long', rho=0., mu=mu, k_d=3., K_t=1.e4)))
    return jobs


def _job(a):
    q, name, cond, tissue, regime, i, o = a
    sc = S3.cell_scene3(cond, tissue, regime, i, base=BASE, overrides=o)
    log, pol, sec = M.run(name, sc, seed=500+i)
    m = M.metric(log, pol, sc)
    m.update(question=q, name=name, cond=cond, regime=regime, i=i, **{k: o.get(k) for k in ('dynamics', 'phase', 'probe', 'contact_law', 'force', 'mu')})
    return m


def table(rows, q, keys, names):
    print(f'\n== {q}: axial RMSE median, normal RMSE, event RMSE, peak max, time over FD+3, loss time, accept, gain error median, call p50/p99 us')
    groups = sorted({tuple(r[k] for k in keys) for r in rows if r['question'] == q}, key=str)
    for g in groups:
        print('  ', dict(zip(keys, g)))
        for n in names:
            s = [r for r in rows if r['question'] == q and tuple(r[k] for k in keys) == g and r['name'] == n]
            ok = [r for r in s if not r['failed']]
            if not ok:
                print('    ', n, 'numerically failed'); continue
            print(f"    {n:6s} {np.median([r['rmse'] for r in ok]):6.3f} {np.median([r['normal_rmse'] for r in ok]):6.3f} {np.nanmedian([r['event_rmse'] for r in ok]) if any(np.isfinite(r['event_rmse']) for r in ok) else float('nan'):6.3f} {max(r['peak'] for r in ok):6.1f} {np.mean([r['over_band'] for r in ok]):6.3f} {np.mean([r['loss_time'] for r in ok]):6.3f} {np.mean([r['accept'] for r in ok]):5.2f} {np.nanmedian([abs(r['est_abs']) for r in ok]):6.3f} {np.median([r['call_p50'] for r in ok]):5.0f}/{np.median([r['call_p99'] for r in ok]):5.0f} fail {sum(r['failed'] for r in s)}")


def main():
    host = sys.argv[1] if len(sys.argv) > 1 else ''
    jobs = plan()
    if host:
        jobs = [(a[0], host+a[1], *a[2:]) for a in jobs]
    print(len(jobs), 'runs', flush=True)
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    if host:
        for r in rows:
            r['name'] = r['name'].replace(host, '')
    (HERE/'results'/f"select3{'_p5' if host else ''}.json").write_text(json.dumps(rows, indent=1), encoding='utf-8')
    table(rows, 'dynamic', ('dynamics', 'phase', 'cond'), ('L', 'H', 'W', 'S'))
    table(rows, 'static', ('probe', 'force', 'cond'), ('L', 'Hz'))
    table(rows, 'sliding', ('regime', 'mu'), ('L', 'LF', 'AS', 'TR'))


if __name__ == '__main__':
    main()
