"""Reproducible tuning, open-loop replay and paired diagnostic closed loop.

python -X utf8 experiments/method_redesign/row3_gls.py --workers 4
Only writes row3_gls.json/md in this directory's results/. No final cohort.
--resume reuses completed summaries (not trajectories); no re-selection.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import step2_estimators as S  # establishes single-thread native kernels
import numpy as np
from gls_estimator import GLS, GLSController
import sweep

HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results/row3_gls.json'


def candidates():
    # Same 12 evaluations as E1/E2/E3, fixed before diagnostic access.
    return [[ell, alpha] for ell in (.008, .02, .05, .12)
            for alpha in (.03, .3, 3.)]


def replay(log, method, cfg):
    if method != 'E4':
        k, sd, ms, _, _ = S.replay(log, method, cfg)
        return k[:, 0], sd[:, 0], ms, np.zeros(len(log), bool)
    obs = S.sanitize(log)
    est = GLS(cfg)
    k, sd, ms, gate = [], [], [], []
    for i in range(len(log)):
        start = time.perf_counter_ns()
        pred, cov = est.step(obs[:i+1], i+1)
        ms.append((time.perf_counter_ns()-start)*1e-6)
        k.append(pred[0]); sd.append(np.sqrt(cov[0, 0]))
        f = est.last
        gate.append(f.get('gate', False) and f['chi2'] <= S.PARAM['chi2_max'])
    return np.array(k), np.array(sd), np.array(ms), np.array(gate)


def tune_job(job):
    cfg, logs = job
    scores = {}
    for fam, sid, log in logs:
        k, _, _, _ = replay(log, 'E4', cfg)
        mask = S.scoring_mask(log)
        e = abs(np.log(k[mask]/log[mask, 8]))
        scores.setdefault(fam, []).append(float(np.median(e)+.25*np.quantile(e, .95)))
    return dict(params=cfg, objective=float(np.mean([np.mean(x) for x in scores.values()])))


def ramp_mask(log, p):
    t = log[:, 0]
    mask = np.zeros(len(log), bool)
    for i in range(int(p[15])):
        onset, amp, ramp, hold = p[16+4*i:20+4*i]
        mask |= ((t >= onset) & (t < onset+ramp)) | ((t >= onset+ramp+hold) & (t < onset+2*ramp+hold))
    return mask


def open_job(job):
    scene, cfg, e2 = job
    fam, sid, p, T = scene
    _, _, log = S.generate(scene)
    good = S.scoring_mask(log)
    with np.errstate(divide='ignore', invalid='ignore'):
        grad = abs(np.gradient(np.log(np.maximum(log[:, 8], 1.)), log[:, 1]))
    transition = good & (np.where(np.isfinite(grad), grad, 0.) > 10.)
    masks = dict(all=good, transition=transition, ramp=good & ramp_mask(log, p))
    rows = []
    for method, params in [('E0', None), ('E2', e2), ('E4', cfg)]:
        k, sd, ms, accepted = replay(log, method, params)
        ratio = k/np.maximum(log[:, 8], 1e-30)
        values = dict(error=abs(np.log(np.maximum(ratio, 1e-30))), under=ratio < 1/1.2,
                      over=ratio > 1/.75, cover1=abs(k-log[:, 8]) <= sd,
                      cover2=abs(k-log[:, 8]) <= 2*sd, ms=ms)
        row = dict(family=fam, sid=sid, method=method, segments={})
        for segment, mask in masks.items():
            row['segments'][segment] = {key: value[mask] for key, value in values.items()}
        if method == 'E4':
            mask = good & accepted
            row['segments']['accepted'] = {key: value[mask] for key, value in values.items()}
            row['acceptance'] = dict(eligible_ticks=int(good.sum()), accepted_ticks=int(mask.sum()))
        rows.append(row)
    return rows


def dist(values):
    if not len(values):
        return dict(ticks=0, median=None, p95=None)
    return dict(ticks=len(values), median=float(np.median(values)), p95=float(np.quantile(values, .95)))


def stats(values):
    result = dist(values['error'])
    for key in ('under', 'over', 'cover1', 'cover2'):
        result[key] = float(np.mean(values[key])) if len(values[key]) else None
    result['p50_ms'] = float(np.median(values['ms'])) if len(values['ms']) else None
    result['p99_ms'] = float(np.quantile(values['ms'], .99)) if len(values['ms']) else None
    return result


def summarize_open(rows):
    summary = {}
    for fam in (*S.C.FAMS, 'overall'):
        summary[fam] = {}
        for method in ('E0', 'E2', 'E4'):
            rr = [r for r in rows if r['method'] == method and (fam == 'overall' or r['family'] == fam)]
            result = {}
            for segment in rr[0]['segments']:
                values = {key: np.concatenate([r['segments'][segment][key] for r in rr])
                          for key in rr[0]['segments'][segment]}
                result[segment] = stats(values)
            if method == 'E4':
                result['acceptance_rate'] = sum(r['acceptance']['accepted_ticks'] for r in rr)/sum(r['acceptance']['eligible_ticks'] for r in rr)
            summary[fam][method] = result
    episodes = [dict(family=r['family'], sid=r['sid'], method=r['method'],
                     segments={k: stats(v) for k, v in r['segments'].items()},
                     acceptance=r.get('acceptance')) for r in rows]
    return dict(summary=summary, episodes=episodes)


def closed_job(job):
    method, group, dose, sid, p, T, cfg = job
    policy = S.C.build('Proposed') if method == 'Proposed' else GLSController(cfg)
    start = time.perf_counter()
    # Warm JIT once per worker outside policy timing; plant episode itself is
    # Python, its force kernel is compiled on the first episode of this worker.
    log, policy = S.common.run(policy, p, T, sid)
    metric = S.common.task_metrics(log, 20.8)
    row = dict(method=method, group=group, dose=dose, sid=sid,
               ok=not metric['diverged'], wall_s=time.perf_counter()-start)
    row.update({k: float(metric[k]) for k in ('rmse', 'peak', 'loss', 'loss_max')})
    # Exclude first second, including compilation/startup and initial QR cache.
    calls = np.array(policy.t_call)[126:]*1000
    row.update(p50_ms=float(np.median(calls)), p99_ms=float(np.quantile(calls, .99)))
    if method == 'E4':
        row.update(accepted=policy.accepted, attempted=policy.attempted)
    return row


def summarize_closed(rows):
    summary = {}
    groups = sorted(set((r['group'], r['dose']) for r in rows), key=str)
    for group, dose in groups:
        rr = [r for r in rows if r['group'] == group and r['dose'] == dose]
        result = {}
        for method in ('Proposed', 'E4'):
            allrows = [r for r in rr if r['method'] == method]
            valid = [r for r in allrows if r['ok']]
            rmse = [r['rmse'] for r in valid]
            result[method] = dict(n=len(allrows), failed=len(allrows)-len(valid),
                median=float(np.median(rmse)), mean=float(np.mean(rmse)),
                runaway=sum(r['peak'] > 50 for r in valid),
                long_loss=sum(r['loss_max'] > .5 for r in valid),
                loss_mean=float(np.mean([r['loss'] for r in valid])),
                p50_ms=float(np.median([r['p50_ms'] for r in valid])),
                p99_ms=float(np.median([r['p99_ms'] for r in valid])))
        a = {r['sid']: r for r in rr if r['method'] == 'Proposed' and r['ok']}
        b = {r['sid']: r for r in rr if r['method'] == 'E4' and r['ok']}
        diff = np.array([b[i]['rmse']-a[i]['rmse'] for i in sorted(a.keys() & b.keys())])
        rng = np.random.default_rng(291101)
        boot = diff[rng.integers(len(diff), size=(10000, len(diff)))].mean(axis=1)
        result['paired'] = dict(n=len(diff), mean_delta=float(diff.mean()),
                                ci95=np.quantile(boot, [.025, .975]).tolist())
        summary[group if dose is None else f'{group}:{dose}'] = result
    return summary


def verify():
    rng = np.random.default_rng(8123)
    t = np.arange(16)*S.DT
    X = np.vander(t/t[-1], 4, increasing=True)
    Y = rng.normal(size=16)*1e-4
    Z = 4000*Y+X@rng.normal(size=4)+rng.normal(size=16)*S.NOISE
    est = GLS((.02, .3))
    V, eig, _ = est.geometry(X, np.arange(16))
    def calculate(z, alpha):
        yr, zr = V.T@Y, V.T@z
        omega = S.NOISE**2+alpha**2*eig
        return (yr*zr/omega).sum()/(yr*yr/omega).sum(), 1/(yr*yr/omega).sum()
    k, var = calculate(Z, .3)
    # Exact invariance even for very large polynomial nuisance.
    assert np.isclose(k, calculate(Z+X@rng.normal(size=4)*100, .3)[0], rtol=1e-9)
    Q = np.linalg.qr(X, mode='complete')[0][:, 4:]
    assert np.isclose(calculate(Z, 0)[0], (Q.T@Y)@(Q.T@Z)/np.linalg.norm(Q.T@Y)**2, rtol=1e-12)
    K = .3**2*np.exp(-.5*((t[:,None]-t[None,:])/.02)**2)
    C = S.NOISE**2*np.eye(12)+Q.T@K@Q
    y, z = Q.T@Y, Q.T@Z
    assert np.isclose(k, y@np.linalg.solve(C,z)/(y@np.linalg.solve(C,y)), rtol=1e-10)
    # Full joint GLS with flat beta agrees with quotient formula.
    A = np.column_stack((Y, X))
    omega = S.NOISE**2*np.eye(16)+Q@(Q.T@K@Q)@Q.T
    H = A.T@np.linalg.solve(omega,A)
    joint = np.linalg.solve(H, A.T@np.linalg.solve(omega,Z))
    assert np.isclose(joint[0], k, rtol=1e-9)
    assert np.isclose(np.linalg.inv(H)[0,0], var, rtol=1e-9)
    # Model-based variance using 20k independent synthetic draws.
    draws = np.linalg.cholesky(C)@rng.normal(size=(12,20000))
    weights = np.linalg.solve(C,y)/(y@np.linalg.solve(C,y))
    empirical = float(np.var(weights@draws, ddof=1))
    assert abs(empirical/var-1) < .04
    scene = S.make_cohort(True)[0]
    _, _, log = S.generate(scene)
    k0, sd0, _, _ = replay(log[:65], 'E4', (.02,.3))
    poison = log[:90].copy()
    poison[:, [3,4,6,7,8,9,10,14,15]] = 1e90
    poison[65:, 5] += 10
    k1, sd1, _, _ = replay(poison, 'E4', (.02,.3))
    assert np.array_equal(k0,k1[:65]) and np.array_equal(sd0,sd1[:65])
    original = S.replay(log[:90], 'E0', None)[0][:,0]
    zero = replay(log[:90], 'E4', (.02,0.))[0]
    assert np.allclose(original, zero, rtol=1e-6, atol=1e-3)
    return dict(polynomial_invariance=True, zero_prior_limit=True, dense_quotient_equivalence=True,
                joint_GLS_equivalence=True, variance_equivalence=True,
                synthetic_variance_ratio=empirical/var, truth_and_future_poisoning=True,
                baseline_limit_max_abs=float(np.max(abs(original-zero))))


def save(result):
    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if args.verify_only:
        print(json.dumps(verify(), indent=2)); return
    result = json.loads(OUTPUT.read_text(encoding='utf-8')) if args.resume and OUTPUT.exists() else {}
    if not result:
        result = dict(status='in_progress', verification=verify(), metadata=dict(
            python=sys.version, numpy=np.__version__, platform=platform.platform(), workers=args.workers,
            tuning_seeds='190000+1000*family_index, 10/family', diagnostic_seeds='110000+1000*family_index',
            sweep_seeds='160000+ (ramps), 170000+ (ratios)', final_access=False,
            bootstrap='10000 paired scene resamples, seed 291101; per-comparison 95%, unadjusted',
            sources={str(p.relative_to(S.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in [Path(__file__), HERE/'gls_estimator.py', HERE/'step2_estimators.py',
                               S.ROOT/'experiments/motion_aware_scan/common.py',
                               S.ROOT/'experiments/cobot_scan_benchmark/plant.py',
                               S.ROOT/'experiments/motion_aware_scan/tuned_v2/Proposed-IMMb.json']}))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        if 'selected' not in result:
            print('Generate 50 tuning trajectories', flush=True)
            logs = list(pool.map(S.generate, S.make_cohort(True)))
            rows = []
            for row in pool.map(tune_job, [(cfg, logs) for cfg in candidates()]):
                rows.append(row); print('tune', row, flush=True)
            result['tuning'] = rows
            result['selected'] = min(rows, key=lambda r:r['objective'])
            save(result)
        cfg = result['selected']['params']
        print('LOCKED CONFIG', cfg, flush=True)
        if 'open_loop' not in result:
            prior = json.loads((HERE/'results/step2_estimators.json').read_text(encoding='utf-8'))
            e2 = prior['selected']['E2']['params']
            result['baseline_E2'] = prior['selected']['E2']
            rows = []
            for i, batch in enumerate(pool.map(open_job, [(s, cfg, e2) for s in S.make_cohort(False)])):
                rows.extend(batch)
                if (i+1)%20 == 0: print('open loop', i+1, '/100', flush=True)
            result['open_loop'] = summarize_open(rows)
            # Check unchanged trajectories/implementations reproduce stored E0/E2 aggregates.
            for fam in (*S.C.FAMS, 'overall'):
                for method in ('E0','E2'):
                    for key in ('median','p95','under','over','cover1','cover2'):
                        assert np.isclose(result['open_loop']['summary'][fam][method]['all'][key], prior['summary'][fam][method][key], rtol=1e-12), (fam,method,key)
            result['verification']['step2_E0_E2_aggregate_reproduction'] = True
            save(result)
        for kind in ('diagnostic', 'sweep'):
            if kind in result: continue
            scenes = [(fam, None, sid, p, T) for fam, sid, p, T in S.C.cohort('diag', 50)] if kind == 'diagnostic' else sweep.scenes()
            jobs = [(method, group, dose, sid, p, T, cfg) for group, dose, sid, p, T in scenes for method in ('Proposed','E4')]
            rows = []
            for i, row in enumerate(pool.map(closed_job, jobs, chunksize=1)):
                rows.append(row)
                if (i+1)%40 == 0: print(kind, i+1, '/', len(jobs), flush=True)
            result[kind] = dict(summary=summarize_closed(rows), episodes=rows)
            save(result)
    result['status'] = 'complete'
    save(result)
    from row3_gls_report import markdown
    OUTPUT.with_suffix('.md').write_text(markdown(result), encoding='utf-8')
    print('COMPLETE', OUTPUT, flush=True)


if __name__ == '__main__':
    main()
