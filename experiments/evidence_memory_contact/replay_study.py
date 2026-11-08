"""Offline study on frozen V1+Proposed trajectories. Truth is used for scoring only.

Part 1 refines weakness 1: does a higher drift order remove the bias during fast motion?
Part 2 compares contact estimators per condition: the 16-tick projection, the 32-tick lean
estimator, single-memory linear filters, and an evidence-weighted bank of memories.

Each linear filter keeps its own likelihood history, so evidence about which memory fits
accumulates across intervals. An IMM merges states every interval and loses that evidence;
a per-interval probe response is far below the sensor noise.
"""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for sub in ('conditional_linear_contact', 'contact_information', 'motion_aware_scan', 'method_redesign'):
    sys.path.insert(0, str(ROOT/sub))
import common
from common import P, K_MIN, K_MAX
from cohort_eval import FINAL, FAMS, cohort
from v1_servo_id import ServoID
from bottleneck_support import twin
from conditional_filter import interval_measurement
import local_contact as LC
from memory_bank import MemoryBank

AFC = 1-math.exp(-2*math.pi*30*.001); EB = math.exp(-.001/.019); DT = .008
CHI = FINAL['chi2_max']; QM = FINAL['q_maneuver']


def run_bank(src, seed_time, **kw):
    bank = MemoryBank(**kw)
    out = np.full(len(src), np.nan); sig = np.full(len(src), np.nan)
    for j in range(1, len(src)-1):
        item = interval_measurement(src, j, 8, EB, AFC, .019)
        if item is None:
            bank.contact_lost(); continue
        H, z = item
        if not bank.ready:
            if src[j, 0] < seed_time:
                continue
            seed = common.window_fit(src, j+2, 16, 3, .04, 8, EB, AFC, 1, True)['k']
            if not np.isfinite(seed) or seed <= 0:
                continue
            bank.start(float(np.clip(seed, K_MIN, K_MAX)), float(src[j, 5]))
        bank.step(z, H, float(src[j+1, 2]-src[j, 2]))
        out[j+1], sig[j+1] = bank.k, bank.sigma
    return out, sig


ESTIMATORS = {
    'single_.3': dict(memories=(.3,), q_motion=QM),
    'single_1': dict(memories=(1.,), q_motion=QM),
    'bank_RW': dict(memories=(.1, .3, 1., 3.), q_motion=QM),
    'bank_RW_h.5': dict(memories=(.1, .3, 1., 3.), q_motion=QM, horizon=.5),
    'bank_slope': dict(memories=(.1, .3, 1., 3.), q_motion=QM, slope=True),
}


def job(scene):
    fam, sid, p, T = scene
    pol = ServoID(single_model=False, **FINAL)
    log, _ = common.run(pol, p, T, sid)
    if log is None:
        return dict(family=fam, rows=None)
    clean = twin(p, log, T)
    kt = log[:, 8]; pa = np.asarray(p, float); n_all = len(log)
    est = {}
    for name, kw in ESTIMATORS.items():
        est[name] = (run_bank(log, 1.0, **kw), run_bank(clean, 1.0, **kw)[0])
    win = {d: (np.full(n_all, np.nan), np.full(n_all, np.nan)) for d in (3, 5, 7)}
    chi = np.full(n_all, np.nan); lean = (np.full(n_all, np.nan), np.full(n_all, np.nan))
    mn = log[:, [0, 1, 2, 5, 11, 12, 13]]; mc = clean[:, [0, 1, 2, 5, 11, 12, 13]]
    for n in range(int(1.0/DT), n_all):
        if np.any(log[n-34:n, 5] <= .5):
            continue
        for d in win:
            a = common.window_fit(log, n, 16, d, .04, 8, EB, AFC, 1, True)
            b = common.window_fit(clean, n, 16, d, .04, 8, EB, AFC, 1, True)
            win[d][0][n-1], win[d][1][n-1] = a['k'], b['k']
            if d == 3:
                chi[n-1] = a['chi2']
        r = LC.fit(mn[:n], 32, .019, 30., .04, False); s = LC.fit(mc[:n], 32, .019, 30., .04, False)
        if r and s:
            lean[0][n-1], lean[1][n-1] = r['k'], s['k']
    dln = np.zeros(n_all); ratio = np.ones(n_all); speed = np.zeros(n_all); excursion = np.zeros(n_all)
    for i in range(17, n_all):
        w = kt[i-16:i+1]
        if w.min() > 1:
            dln[i] = math.log(w.max()/w.min()); ratio[i] = w[-1]/w[0]
        speed[i] = abs(P.surface(float(log[i, 1]), float(log[i, 0]), pa)[1])
        excursion[i] = np.max(np.abs(log[i-16:i+1, 3]-4.))
    mask0 = (log[:, 0] >= 1.5) & (kt > 0) & (log[:, 5] > .5)
    groups = dict(quiet=(dln < .03) & (speed < .005),
                  motion_const=(dln < .03) & (speed >= .02),
                  motion_small_change=(dln >= .03) & (dln < .15) & (speed >= .02),
                  rise=(dln >= .3) & (ratio > 1) & (speed < .005),
                  fall=(dln >= .3) & (ratio < 1) & (speed < .005),
                  change_motion=(dln >= .15) & (speed >= .01),
                  large_force=(excursion > 1.5) & (speed >= .01))
    rows = {}
    for g, mask in groups.items():
        m = mask & mask0
        if m.sum() < 5:
            continue
        block = {}
        items = dict(win3=win[3], win5=win[5], win7=win[7], lean32=lean)
        items.update({k: (v[0][0], v[1]) for k, v in est.items()})
        for name, (a, c) in items.items():
            ok = m & np.isfinite(a) & np.isfinite(c) & (a > 0) & (c > 0)
            if ok.sum() < 3:
                continue
            block[name] = dict(total=np.log(a[ok]/kt[ok]).tolist(), struct=np.log(c[ok]/kt[ok]).tolist(), noise=np.log(a[ok]/c[ok]).tolist())
        block['chi_flag'] = float(np.mean(chi[m & np.isfinite(chi)] > CHI)) if np.any(m & np.isfinite(chi)) else None
        rows[g] = block
    return dict(family=fam, rows=rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--scenes', type=int, default=12); ap.add_argument('--label', default='replay_v1')
    a = ap.parse_args()
    with ProcessPoolExecutor(max_workers=15) as e:
        out = list(e.map(job, cohort('diag', a.scenes), chunksize=1))
    summary = {}
    names = ['win3', 'win5', 'win7', 'lean32']+list(ESTIMATORS)
    for fam in FAMS:
        gs = {g for o in out if o['family'] == fam and o['rows'] for g in o['rows']}
        for g in sorted(gs):
            line = {}
            for name in names:
                vs = [o['rows'][g][name] for o in out if o['family'] == fam and o['rows'] and g in o['rows'] and name in o['rows'][g]]
                if not vs:
                    continue
                t = np.concatenate([v['total'] for v in vs]); s = np.concatenate([v['struct'] for v in vs]); n = np.concatenate([v['noise'] for v in vs])
                line[name] = dict(n=len(t), total=float(np.percentile(abs(t), 90)), struct=float(np.percentile(abs(s), 90)), noise=float(np.std(n)), med=float(np.median(s)))
            summary[f'{fam}|{g}'] = line
            print(f'{fam:10s}{g:20s}', ' '.join(f"{k}={v['total']:.3f}/{v['struct']:.3f}/{v['noise']:.3f}" for k, v in line.items()), 'n', line['win3']['n'] if 'win3' in line else 0, flush=True)
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/f'{a.label}.json').write_text(json.dumps(summary, indent=1), encoding='utf8')


if __name__ == '__main__':
    main()
