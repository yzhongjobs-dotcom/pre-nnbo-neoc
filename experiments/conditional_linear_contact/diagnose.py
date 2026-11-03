"""Separate implementation, interface and method causes of the V0 result.

Truth is read only for scoring. Diagnostic cohort only.
  tracker   V1+Proposed with only the tuned tracker parameters that leaked into the
            motion IMM (p_stay, r_scale, noise): is the median gain from those?
  isolated  Conditional filter with the motion IMM kept at FINAL and the tuned filter
            parameters given only to the new filter, gate at ln 1.2 or tuned.
  replay    Single-mode filter on frozen V1+Proposed trajectories across q_lnk, with
            and without the IMM: does long memory actually beat the 16-tick window?
  modes     Stiffness-mode identifiability from stored mode probabilities.
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
sys.path.insert(0, str(HERE.parent/'motion_aware_scan'))
sys.path.insert(0, str(HERE.parent/'method_redesign'))
sys.path.insert(0, str(HERE.parent/'contact_information'))
import common
from cohort_eval import FINAL, FAMS, cohort
from v1_servo_id import ServoID
import conditional_filter as CF
from conditional_filter import ConditionalLinearContact, FourModeIMM, MeasuredOnly, interval_measurement

RES = HERE/'results'
TUNED = json.loads((RES/'tuned_params.json').read_text())
NOISE = .04


class Isolated(ConditionalLinearContact):
    """Filter gets its own noise/r_scale/p_stay; tracker and startup stay at FINAL."""
    def __init__(self, filter_noise, filter_r, filter_p_stay, **kw):
        super().__init__(**kw)
        self.filter_args = (filter_noise, filter_r, filter_p_stay)

    def reset(self, platform):
        super().reset(platform)
        n, r, s = self.filter_args
        p = self.params
        self.cl_filter = FourModeIMM(dt=self.Tc, q_smooth=p['q_smooth'], q_maneuver=p['q_maneuver'],
                                     p_stay_motion=s, q_lnk0=p['q_lnk0'], q_lnk1=p['q_lnk1'],
                                     p_stay_stiff=p['p_stay_stiff'], R=n**2, r_scale=r)


def metrics(log, pol):
    m = common.task_metrics(log, 20.8)
    return dict(rmse=float(m['rmse']), peak=float(m['peak']), loss=float(m['loss_max']))


def tracker_job(job):
    scene, name = job; fam, sid, p, T = scene
    q = dict(FINAL)
    if name == 'tracker_tuned':
        q.update(p_stay=TUNED['p_stay'], r_scale=TUNED['r_scale'], noise=NOISE*TUNED['noise_scale'])
    elif name == 'tracker_r':
        q.update(r_scale=TUNED['r_scale'])
    elif name == 'tracker_p':
        q.update(p_stay=TUNED['p_stay'])
    pol = ServoID(single_model=False, **q)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    return dict(config=name, family=fam, sid=sid, **metrics(log, pol)) if log is not None else dict(config=name, family=fam, sid=sid, rmse=None, peak=None, loss=None)


ISO = {
    'iso_tuned_gate': dict(gain_rel_max=TUNED['gain_rel_max']),
    'iso_band_gate': dict(gain_rel_max=math.log(1.2)),
    'iso_band_gate_fb': dict(gain_rel_max=math.log(1.2), c=FINAL['c'], rho_d=FINAL['rho_d']),
}


def iso_job(job):
    scene, name = job; fam, sid, p, T = scene
    q = dict(FINAL)
    q.update(q_lnk0=TUNED['q_lnk0'], q_lnk1=TUNED['q_lnk1'], p_stay_stiff=TUNED['p_stay_stiff'],
             consistency_max=TUNED['consistency_max'], c=TUNED['c'], rho_d=TUNED['rho_d'])
    q.update(ISO[name])
    pol = Isolated(NOISE*TUNED['noise_scale'], TUNED['r_scale'], TUNED['p_stay'], single_model=False, **q)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    if log is None:
        return dict(config=name, family=fam, sid=sid, rmse=None, peak=None, loss=None)
    mu = np.asarray([x['mu'] for x in pol.mode_history]) if pol.mode_history else np.zeros((1, 4))
    return dict(config=name, family=fam, sid=sid, **metrics(log, pol),
                accept=pol.accepted_count/max(len(pol.accept_history), 1),
                stiff_mode_mean=float(mu[:, [1, 3]].sum(axis=1).mean()),
                stiff_mode_sd=float(mu[:, [1, 3]].sum(axis=1).std()))


class SingleMode:
    """One linear filter, same observation and time update as the four-mode IMM."""
    def __init__(self, q, q_lnk, R):
        self.f = FourModeIMM(dt=.008, q_smooth=q, q_maneuver=q, p_stay_motion=1., q_lnk0=q_lnk, q_lnk1=q_lnk,
                             p_stay_stiff=1., R=R, r_scale=1.)


def replay_job(scene):
    fam, sid, p, T = scene
    pol = ServoID(single_model=False, **FINAL)
    log, _ = common.run(pol, p, T, sid)
    if log is None:
        return dict(family=fam, sid=sid, rows=None)
    from bottleneck_support import twin
    clean = twin(p, log, T)
    afc = 1-math.exp(-2*math.pi*30*.001); eb = math.exp(-.001/.019)
    R = NOISE**2
    grid = [.03, .1, .3, 1., 3.]
    qs = (FINAL['q_smooth'], FINAL['q_maneuver'])
    run = {}
    for q in qs:
        for qk in grid:
            for label, src in (('n', log), ('c', clean)):
                f = SingleMode(q, qk, R).f; out = np.full(len(log), np.nan); sig = np.full(len(log), np.nan)
                for j in range(1, len(log)-1):
                    if min(src[j, 5], src[j+1, 5]) <= .5:
                        continue
                    item = interval_measurement(src, j, 8, eb, afc, .019)
                    if item is None:
                        continue
                    H, z = item
                    if not f.initialized:
                        if j < 130:
                            continue     # seed from a projection window at about 1 s, never from truth
                        seed = common.window_fit(src, j+2, 16, 3, .04, 8, eb, afc, 1, True)['k']
                        f.initialize_if_needed(float(np.clip(seed if np.isfinite(seed) else 3e4, 200., 3e5)), float(src[j, 5]))
                    dy = float(src[j+1, 2]-src[j, 2])
                    F = np.array([[1., 0, 0, 0], [dy, 1., .008, .008**2/2], [0, 0, 1., .008], [0, 0, 0, 1.]])
                    f.update_interval(z, H, F, max(float(f.x[0]), 200.))
                    out[j+1] = f.x[0]; sig[j+1] = math.sqrt(max(f.P[0, 0], 0.))
                run[(q, qk, label)] = (out, sig)
    win = np.full(len(log), np.nan); winc = np.full(len(log), np.nan)
    for n in range(int(1.0/.008), len(log)):
        if np.all(log[n-18:n, 5] > .5):
            a = common.window_fit(log, n, 16, 3, .04, 8, eb, afc, 1, True); b = common.window_fit(clean, n, 16, 3, .04, 8, eb, afc, 1, True)
            win[n-1], winc[n-1] = a['k'], b['k']
    t = log[:, 0]; kt = log[:, 8]
    dln = np.zeros(len(log)); ht = np.zeros(len(log))
    P = common.P; pa = np.asarray(p, float)
    for i in range(32, len(log)):
        w = kt[i-32:i]
        if w.min() > 1:
            dln[i] = math.log(w.max()/w.min())
    for i in range(len(log)):
        ht[i] = abs(P.surface(float(log[i, 1]), float(log[i, 0]), pa)[1])
    groups = dict(calm=(ht < .005) & (dln < .05), fast=ht >= .02, change=dln >= .15)
    mask0 = (t >= 1.5) & (kt > 0) & (log[:, 5] > .5)

    def score(e, ec, s=None, g=None):
        m = mask0 & g & np.isfinite(e) & np.isfinite(ec) & (e > 0) & (ec > 0)
        if m.sum() < 5:
            return None
        return dict(n=int(m.sum()), total=np.log(e[m]/kt[m]).tolist(), struct=np.log(ec[m]/kt[m]).tolist(),
                    noise=np.log(e[m]/ec[m]).tolist(),
                    cover=None if s is None else float(np.mean(np.abs(e[m]-kt[m]) <= 2*s[m])))
    rows = {}
    for g, mask in groups.items():
        rows[f'{g}|window'] = score(win, winc, None, mask)
        for q in qs:
            for qk in grid:
                rows[f'{g}|q{q:.3g}|lnk{qk}'] = score(run[(q, qk, 'n')][0], run[(q, qk, 'c')][0], run[(q, qk, 'n')][1], mask)
    return dict(family=fam, sid=sid, rows=rows)


def summarize_closed(rows, configs):
    table = {}
    for c in configs:
        for fam in FAMS:
            r = [x for x in rows if x['config'] == c and x['family'] == fam]
            ok = [x['rmse'] for x in r if x['rmse'] is not None]
            table[f'{c}|{fam}'] = dict(median=float(np.median(ok)), mean=float(np.mean(ok)),
                r=sum(1 for x in r if x['rmse'] is None or x['peak'] > 50),
                L=sum(1 for x in r if x['loss'] is not None and x['loss'] > .5),
                accept=float(np.mean([x['accept'] for x in r if 'accept' in x])) if any('accept' in x for x in r) else None,
                stiff_mean=float(np.mean([x['stiff_mode_mean'] for x in r if 'stiff_mode_mean' in x])) if any('stiff_mode_mean' in x for x in r) else None,
                stiff_sd=float(np.mean([x['stiff_mode_sd'] for x in r if 'stiff_mode_sd' in x])) if any('stiff_mode_sd' in x for x in r) else None)
    return table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('part', choices=['tracker', 'isolated', 'replay'])
    ap.add_argument('--scenes', type=int, default=50)
    ap.add_argument('--workers', type=int, default=15)
    a = ap.parse_args()
    scenes = cohort('diag', a.scenes)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        if a.part == 'tracker':
            configs = ('tracker_tuned', 'tracker_r', 'tracker_p')
            rows = list(pool.map(tracker_job, [(s, c) for c in configs for s in scenes], chunksize=4))
            out = summarize_closed(rows, configs)
        elif a.part == 'isolated':
            configs = tuple(ISO)
            rows = list(pool.map(iso_job, [(s, c) for c in configs for s in scenes], chunksize=4))
            out = summarize_closed(rows, configs)
        else:
            rows = list(pool.map(replay_job, scenes, chunksize=1))
            out = {}
            for fam in FAMS:
                keys = set()
                for r in rows:
                    if r['family'] == fam and r['rows']:
                        keys |= {k for k, v in r['rows'].items() if v}
                for k in sorted(keys):
                    vs = [r['rows'][k] for r in rows if r['family'] == fam and r['rows'] and r['rows'].get(k)]
                    tot = np.concatenate([v['total'] for v in vs]); st = np.concatenate([v['struct'] for v in vs]); no = np.concatenate([v['noise'] for v in vs])
                    cov = [v['cover'] for v in vs if v['cover'] is not None]
                    out[f'{fam}|{k}'] = dict(n=int(len(tot)), total_p90=float(np.percentile(abs(tot), 90)), low_p5=float(np.percentile(tot, 5)),
                                             struct_p90=float(np.percentile(abs(st), 90)), noise_sd=float(np.std(no)),
                                             cover2=float(np.mean(cov)) if cov else None)
    (RES/f'diagnose_{a.part}.json').write_text(json.dumps(dict(summary=out, rows=None if a.part == 'replay' else rows), indent=1), encoding='utf8')
    if a.part != 'replay':
        for k, v in out.items():
            print(k, {x: (round(y, 3) if isinstance(y, float) else y) for x, y in v.items() if y is not None})
    else:
        for k, v in out.items():
            print(k, {x: (round(y, 3) if isinstance(y, float) else y) for x, y in v.items() if y is not None})


if __name__ == '__main__':
    main()
