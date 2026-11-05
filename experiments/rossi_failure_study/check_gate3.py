"""Is the prior-gain gate reliable where the prior is wrong? Stiffness jumps up and down, first contact, re-contact, for both gates."""
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import method3 as M
import physics_plant3 as P3
import scenes3 as S3

HERE = Path(__file__).resolve().parent
NAMES = ('L', 'P5:L', 'H', 'P5:H')
BASE = dict(pose_error=0., mu=.1, rho=0., accel=.1, speed=.01, k_d=2., kr=1.e5, K_t=2.e4, surface_amp=.0003, field_sigma=0., force=8.,
            contact_law='foundation', probe='linear', slide='short', phase=14., dynamics='neural_white')
JUMP_T = 4.
CASES = {'up x3': dict(ratio=3.), 'up x10': dict(ratio=10.), 'down /3': dict(ratio=1/3), 'down /10': dict(ratio=.1),
         'first contact soft': dict(gap=.004, tissue='relaxed'), 'first contact stiff': dict(gap=.004, tissue='active'), 'recontact': dict(pull=.012, tissue='relaxed')}


def scene(case, i):
    spec = CASES[case]
    sc = S3.cell_scene3('rest', spec.get('tissue', 'relaxed' if spec.get('ratio', 1.) > 1. else 'active'), 'nominal', i, base=445000, overrides=BASE)
    if 'ratio' in spec:
        grid = sc['field_x0']+np.arange(len(sc['q_field']))*sc['dx']
        xj = np.interp(JUMP_T, np.arange(len(sc['xcmd']))*.001, sc['xcmd'])
        step = np.exp(math.log(spec['ratio'])*S3._smooth((grid-xj)/.004))
        sc['q_field'] = sc['q_field']*step; sc['E_field'] = sc['E_field']*step
    if 'gap' in spec:
        sc['start_gap'] = spec['gap']
    if 'pull' in spec:
        t = np.arange(len(sc['u']))*.001
        sc['u'] = sc['u']-spec['pull']*(S3._smooth((t-JUMP_T)/.1)-S3._smooth((t-JUMP_T-.4)/.1))
    return sc


def _job(a):
    case, name, gate, i = a
    sc = scene(case, i)
    pol = M.build(name, sc); pol.gate_on_prior = gate
    log = P3.episode(sc, pol, 100+i)
    if log is None:
        return dict(case=case, name=name, gate=gate, i=i, failed=True)
    t = log[:, 0]; f = log[:, 3]; fd = sc['FD']
    start = 0. if 'gap' in CASES[case] else JUMP_T
    after = t >= start
    # settling: last time outside the band, in the two seconds that follow
    window = after & (t < start+3.)
    contact = t[np.argmax(f > .5)] if (f > .5).any() else start
    if 'gap' in CASES[case]:
        start = float(contact); after = t >= start; window = after & (t < start+3.)
    out = np.where(np.abs(f[window]-fd) > 1.)[0]
    settle = float(t[window][out[-1]]-start) if len(out) else 0.
    acc = np.array(pol.accept_trace)
    ok = acc[:, 4] > 0 if len(acc) else np.zeros(0, bool)
    errs = {}
    for lag in (.25, .5, 1., 2.):
        tk = start+lag
        m = acc[:, 0] <= tk
        j = int(np.searchsorted(t, tk))
        if len(acc) and ok[m].any() and j < len(t):
            errs[lag] = float(acc[m][ok[m]][-1, 1]/max(log[j, 21], 1.))
    return dict(case=case, name=name, gate=gate, i=i, failed=False, peak=float(f[after].max()), settle=settle, err=errs, k_true=float(log[-1, 21]),
                k_good=float(pol.k_good), k_err_end=float(pol.k_good/log[-1, 21]), ever=float(t[np.argmax(f > .5)]) if (f > .5).any() else float('nan'),
                rmse=float(np.sqrt(np.mean((f[window]-fd)**2))))


def main():
    jobs = [(c, n, g, i) for c in CASES for n in NAMES for g in (True, False) for i in range(5)]
    with ProcessPoolExecutor(max_workers=15) as ex:
        rows = list(ex.map(_job, jobs, chunksize=1))
    (HERE/'results'/'gate3_checks.json').write_text(json.dumps(rows, indent=1), encoding='utf-8')
    med = lambda q, k: float(np.nanmedian([r[k] for r in q])) if q else float('nan')
    print('case                 name gate   runaway peak  settle s  rmse | last accepted gain / |K| at .25 .5 1 2 s after | end k_good/|K|')
    for c in CASES:
        for n in NAMES:
            for g in (True, False):
                q = [r for r in rows if r['case'] == c and r['name'] == n and r['gate'] == g]
                ok = [r for r in q if not r['failed']]
                run = sum(r['failed'] or r['peak'] > 50 for r in q)
                e = [float(np.nanmedian([r['err'].get(lag, np.nan) for r in ok])) if ok else float('nan') for lag in (.25, .5, 1., 2.)]
                print(f"{c:20s} {n:5s} {'prior' if g else 'self ':5s} {run}/{len(q)} {med(ok, 'peak'):6.1f} {med(ok, 'settle'):6.2f} {med(ok, 'rmse'):6.2f} | {e[0]:5.2f} {e[1]:5.2f} {e[2]:5.2f} {e[3]:5.2f} | {med(ok, 'k_err_end'):5.2f}")


if __name__ == '__main__':
    main()
