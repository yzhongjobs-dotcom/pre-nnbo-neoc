"""Same-input replay of the spectral window against projection and lean. Truth scores only."""
import argparse, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
HERE = Path(__file__).resolve().parent
for sub in ('conditional_linear_contact', 'contact_information', 'motion_aware_scan', 'method_redesign'):
    sys.path.insert(0, str(HERE.parent/sub))
import common
from common import P
from cohort_eval import FINAL, FAMS, cohort
from v1_servo_id import ServoID
from bottleneck_support import twin
import local_contact as LC
from spectral_window import SpectralWindow, AdaptiveSpectralWindow

AFC = 1-math.exp(-2*math.pi*30*.001); EB = math.exp(-.001/.019); DT = .008
CFG = {'sw32_c16': dict(n=32, corner_hz=16.)}
FILTERS = {k: SpectralWindow(**v) for k, v in CFG.items()}
ADAPT = {'ad32': dict(n=32), 'ad32_s1': dict(n=32, change_prior=1.), 'ad32_s.5': dict(n=32, change_prior=.5), 'ad32_s.2': dict(n=32, change_prior=.2), 'ad32_const': dict(n=32, slope=False), 'ad48_s.5': dict(n=48, change_prior=.5)}
FILTERS.update({k: AdaptiveSpectralWindow(**v) for k, v in ADAPT.items()})
CFG = dict(CFG); CFG.update(ADAPT)
CHI = FINAL['chi2_max']


def job(scene):
    fam, sid, p, T = scene
    pol = ServoID(single_model=False, **FINAL)
    log, _ = common.run(pol, p, T, sid)
    if log is None:
        return dict(family=fam, rows=None)
    clean = twin(p, log, T)
    N = len(log); kt = log[:, 8]; pa = np.asarray(p, float)
    res = {k: [np.full(N, np.nan), np.full(N, np.nan), np.full(N, np.nan), np.full(N, np.nan)] for k in list(CFG)+['win3', 'win5', 'lean32']}
    mn = log[:, [0, 1, 2, 5, 11, 12, 13]]; mc = clean[:, [0, 1, 2, 5, 11, 12, 13]]
    for n in range(int(1.0/DT), N):
        if np.any(log[max(n-52, 0):n, 5] <= .5):
            continue
        for d, name in ((3, 'win3'), (5, 'win5')):
            a = common.window_fit(log, n, 16, d, .04, 8, EB, AFC, 1, True); b = common.window_fit(clean, n, 16, d, .04, 8, EB, AFC, 1, True)
            res[name][0][n-1], res[name][1][n-1] = a['k'], b['k']
            res[name][2][n-1] = 1/math.sqrt(a['info']) if a['info'] > 0 else np.nan; res[name][3][n-1] = a['chi2']
        r = LC.fit(mn[:n], 32, .019, 30., .04, False); s = LC.fit(mc[:n], 32, .019, 30., .04, False)
        if r and s:
            res['lean32'][0][n-1], res['lean32'][1][n-1], res['lean32'][2][n-1], res['lean32'][3][n-1] = r['k'], s['k'], r['sigma'], r['chi2']
        for name, f in FILTERS.items():
            m = f.n
            a = f.fit(log[n-m-1:n, 2], log[n-m-1:n, 5]); b = f.fit(clean[n-m-1:n, 2], clean[n-m-1:n, 5])
            if a and b:
                res[name][0][n-1], res[name][1][n-1], res[name][2][n-1], res[name][3][n-1] = a['k'], b['k'], a['sigma'], a['chi2']
    dln = np.zeros(N); ratio = np.ones(N); speed = np.zeros(N); exc = np.zeros(N)
    for i in range(17, N):
        w = kt[i-16:i+1]
        if w.min() > 1:
            dln[i] = math.log(w.max()/w.min()); ratio[i] = w[-1]/w[0]
        speed[i] = abs(P.surface(float(log[i, 1]), float(log[i, 0]), pa)[1]); exc[i] = np.max(np.abs(log[i-16:i+1, 3]-4.))
    mask0 = (log[:, 0] >= 1.5) & (kt > 0) & (log[:, 5] > .5)
    groups = dict(quiet=(dln < .03) & (speed < .005), motion_const=(dln < .03) & (speed >= .02), motion_small=(dln >= .03) & (dln < .15) & (speed >= .02),
                  rise=(dln >= .3) & (ratio > 1) & (speed < .005), fall=(dln >= .3) & (ratio < 1) & (speed < .005),
                  change_motion=(dln >= .15) & (speed >= .01), large_force=(exc > 1.5) & (speed >= .01))
    rows = {}
    for g, mask in groups.items():
        m = mask & mask0
        if m.sum() < 5:
            continue
        rows[g] = {}
        for name, (a, c, s, chi) in res.items():
            ok = m & np.isfinite(a) & np.isfinite(c) & (a > 0) & (c > 0) & np.isfinite(s)
            if ok.sum() < 3:
                continue
            e = np.log(a[ok]/kt[ok])
            rows[g][name] = dict(total=e.tolist(), struct=np.log(c[ok]/kt[ok]).tolist(), noise=np.log(a[ok]/c[ok]).tolist(),
                                 z=(np.abs(a[ok]-kt[ok])/np.maximum(s[ok], 1e-9)).tolist(), rel=(s[ok]/a[ok]).tolist(), chi=chi[ok].tolist())
    return dict(family=fam, rows=rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--scenes', type=int, default=12); ap.add_argument('--label', default='spectral_replay_v1')
    a = ap.parse_args()
    with ProcessPoolExecutor(max_workers=15) as e:
        out = list(e.map(job, cohort('diag', a.scenes), chunksize=1))
    names = ['win3', 'win5', 'lean32']+list(CFG); summary = {}
    for fam in FAMS:
        for g in ('quiet', 'motion_const', 'motion_small', 'rise', 'fall', 'change_motion', 'large_force'):
            line = {}
            for name in names:
                vs = [o['rows'][g][name] for o in out if o['family'] == fam and o['rows'] and g in o['rows'] and name in o['rows'][g]]
                if not vs:
                    continue
                c = {k: np.concatenate([v[k] for v in vs]) for k in vs[0]}
                line[name] = dict(n=len(c['total']), total=float(np.percentile(abs(c['total']), 90)), struct=float(np.percentile(abs(c['struct']), 90)),
                                  noise=float(np.std(c['noise'])), cov2=float(np.mean(c['z'] <= 2)), rel=float(np.median(c['rel'])))
            if line:
                summary[f'{fam}|{g}'] = line
                print(f'{fam:10s}{g:14s}', ' '.join(f"{k}={v['total']:.2f}/{v['struct']:.2f}/{v['noise']:.2f}" for k, v in line.items()), 'n', next(iter(line.values()))['n'], flush=True)
    (HERE/'results').mkdir(exist_ok=True); (HERE/'results'/f'{a.label}.json').write_text(json.dumps(summary, indent=1), encoding='utf8')


if __name__ == '__main__':
    main()
