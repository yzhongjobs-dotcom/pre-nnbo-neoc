"""Three comparison tables on one scoring rule, reference k_eff; no tuning.

Run in order; each stage reads only what the previous one saved:
  python -X utf8 identification_completion.py --stage closed   # results/completion_closed_loop.json
  python -X utf8 identification_completion.py --stage replay   # results/id_replay_ticks/*.npz
  python -X utf8 identification_completion.py --stage score    # results/identification_replay.json
  python -X utf8 identification_completion.py --stage oracle_diagnosis  # results/oracle_diagnosis.json
  python -X utf8 identification_completion.py --stage tables   # comparison_tables.md

closed: corrected_comparison.cases() (36 paired cases). Only methods whose control
changed run again; other rows in completion_closed_loop.json are kept.

replay: every identifier on identical recorded observations of the same 36 cases.
Primary source Huber-D4; Rossi-ID is a second source, only to see whether
conclusions depend on the host. Per-tick data is saved, so scoring needs no plant.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import corrected_comparison as CC
import id_baselines as I
import overlay_compare as O
import rib_transition_eval as R
import memory_segment_eval as SE
import simple_ls_validation as S

P3, B = SE.P3, O.B
HERE = Path(__file__).resolve().parent
TICKS = HERE/'results/id_replay_ticks'
# Xue outputs a tangent stiffness that cannot be converted to k_eff without truth: control table only.
MAIN = ('Huber-D4', 'D4', 'P5', 'Rossi-ID', 'Kalman', 'MHE', 'Simple-1', 'Simple-2')
CONTROL = ('Huber-D4', 'D4', 'P5', 'Rossi-ID', 'Kalman', 'MHE', 'Xue-const', 'Simple-1', 'Simple-2')
INTERNAL = ('off_none',) + CC.MEMORY
SOURCES = {'Huber-D4': MAIN + INTERNAL, 'Rossi-ID': MAIN}
FULL = ('Rossi-control', 'DO-MPAC', 'Xing-DAHIC')
UNGATED = CC.OWN_MECHANISM
RERUN = UNGATED + CC.IDEAL_CAL
KEPT = ('Huber-D4', 'D4', 'P5') + FULL + tuple(CC.ORACLES)
SEGMENTS = ('plateau', 'transition', 'event', 'event_only', 'event_transition')
CHANNELS = ('raw', 'k_good')
LOG12 = math.log(1.2)
LOG085 = math.log(.85)
SOURCE_COLS = [0, 1, 2, 3, 5, 18, 21, 26, 27]


def stored_rows():
    data = json.loads((HERE/'results/corrected_comparison.json').read_text(encoding='utf-8'))
    return {(r['method'], r['index']): r for r in data['episodes']}


def closed_job(task):
    row = CC.job(task)
    for key in ('transitions', 'excluded_prior_contact_loss', 'plateau', 'plateau_ticks'):
        row.pop(key, None)  # stiffness scoring belongs to the identical-trajectory replay
    return row


def closed():
    start = time.perf_counter()
    path = HERE/'results/completion_closed_loop.json'
    previous = json.loads(path.read_text(encoding='utf-8'))
    kept = [r for r in previous['rows'] if r['method'] in KEPT]
    assert len(kept) == 36*len(KEPT)
    tasks = [(m, c) for m in RERUN for c in CC.cases()]
    with ProcessPoolExecutor(8) as ex:
        fresh = [dict(r, method=m) for (m, _), r in zip(tasks, ex.map(closed_job, tasks, chunksize=1))]
    # Sanity on kept rows: gates-on Huber must still match its stored 612-grid values.
    old = stored_rows()
    worst = max(max(abs(r[k]-old[('Huber-D4', r['index'])][k]) for k in ('L', 'force_rmse', 'peak')) for r in kept if r['method'] == 'Huber-D4')
    assert worst == 0.
    rows = kept+fresh
    now = {(r['method'], r['index']): r for r in rows}
    gate = [abs(now[('oracle_keff_on', c['index'])]['L']-now[('oracle_keff_off', c['index'])]['L']) for c in CC.cases()]
    metadata = dict(previous['metadata'])
    metadata.update(
        own_mechanism=', '.join(UNGATED)+': no Method acceptance gates, feed-forward jump rule or CENTRE; own estimator mechanisms, finite and contact checks, gain clipping and k_bar bound with rate limit kept. Huber-D4, D4 and P5 keep the gates they were designed with.',
        xue_chi='Dropped: its only difference from Xue-const was a residual gate this project added; a separately selected configuration is not a different method.',
        ideal_cal='Huber-D4_idealcal, oracle_keff_idealcal: same scenes with exact mass, force offset and servo lag in the sensing stage, a diagnostic, not a plausible calibration.',
        huber_regression_max_abs=worst, oracle_gate_effect_max_abs_dL=max(gate))
    for key in ('rerun', 'rerun_regression_max_abs', 'simple_regression_max_abs'):
        metadata.pop(key, None)
    B.dump('completion_closed_loop.json', dict(metadata=metadata, rows=rows, seconds=time.perf_counter()-start))
    print('closed done', round(time.perf_counter()-start), 's; Huber regression', worst, flush=True)


class Recorder:
    """Run a policy unchanged and keep each pre-tick observation row."""
    def __init__(self, pol):
        self.pol = pol; self.observed = []

    def reset(self, platform):
        self.observed = []; self.pol.reset(platform)

    def __call__(self, o):
        self.observed.append(o.log[o.n-1, B.OBS_COLS].copy())
        return self.pol(o)


def build_replay(name, sc, case):
    host = 'Huber-D4' if name in INTERNAL or name.startswith('Simple') else name
    cfg = O.chosen()[host]
    pol = B.build(host, cfg, sc)
    if name.startswith('Simple'):
        pol.estimator = I.ContactAdapter(S.SimpleLS(name == 'Simple-2', nw=24, sigma_scale=cfg['sigma_scale'], chi_scale=cfg['gate']))
    if name in INTERNAL:
        SE.gates_off(pol)
    if name in UNGATED:
        CC.own_mechanism(pol)
    if name in CC.MEMORY:
        R.make_filter(pol, name)
    return pol


def replay_method(name, sc, case, log, observed, source):
    pol = build_replay(name, sc, case)
    pol.reset(dict(P3.PLATFORM, tau_a=sc['tau_a_hat']))
    count = len(log)
    raw = np.full((count, 3), np.nan); elapsed = []
    inner = pol.estimator.estimate

    def record(o, n, p):
        tic = time.perf_counter()
        k, sig, chi = inner(o, n, p)
        elapsed.append(time.perf_counter()-tic)
        raw[n-1] = k, sig, chi
        return k, sig, chi

    pol.estimator.estimate = record
    good, bar = np.full(count, np.nan), np.full(count, np.nan)
    worst = 0.
    limit = P3.PLATFORM['VMAX']*P3.PLATFORM['Tc']
    for n in range(1, count+1):
        command = pol(B.source_observation(log, observed, n))
        good[n-1] = pol.k_good; bar[n-1] = pol.kbar
        if name == source:
            sent = observed[n-1, B.OBS_COLS.index(11)]
            worst = max(worst, abs(float(np.clip(command, sent-limit, sent+limit)-log[n-1, 11])))
    if name == source:
        assert worst < 1e-10, ('Causal replay differs', name, case['index'], worst)
    ok = np.zeros(count, bool)
    t = log[:, 0]
    for tk, _, _, _, accepted in pol.accept_trace:
        if accepted:
            ok[min(int(np.searchsorted(t, tk-1e-9)), count-1)] = True
    valid = np.isfinite(raw[:, 0]) & (raw[:, 0] > 0)
    return {'raw': raw[:, 0].astype(np.float32), 'sigma': raw[:, 1].astype(np.float32), 'chi': raw[:, 2].astype(np.float32),
            'good': good.astype(np.float32), 'bar': bar.astype(np.float32), 'ok': ok, 'valid': valid,
            'call_us': (np.asarray(elapsed)*1e6).astype(np.float32)}, worst


def case_job(task):
    source, case, names = task
    O.setup()
    cell = tuple(case['cell'])
    sc = SE.MC.scene_for(cell, CC.BASE, case['variant'])
    pol = build_replay(source, sc, case)
    recorder = Recorder(pol)
    log = P3.episode(sc, recorder, case['seed'])
    if log is None:
        raise RuntimeError(('Source failed', source, case))
    metric = SE.M.metric(log, pol, sc)
    observed = np.asarray(recorder.observed)
    arrays = dict(source_log=log[:, SOURCE_COLS], source_cols=np.asarray(SOURCE_COLS))
    errors = {}
    for name in names:
        out, worst = replay_method(name, sc, case, log, observed, source)
        errors[name] = worst
        arrays.update({f'{name}__{k}': v for k, v in out.items()})
    meta = dict(source=source, index=case['index'], group=case['group'], cell=case['cell'], seed=case['seed'],
                variant=case['variant'], event_info=sc['event_info'], methods=list(names),
                source_L=float(np.log1p(metric['J'])), source_rmse=metric['rmse'], source_peak=metric['peak'],
                command_replay_error=errors[source])
    arrays['meta'] = np.asarray(json.dumps(meta, default=float))
    np.savez_compressed(TICKS/f"{source}_{case['index']:02d}.npz", **arrays)
    return meta


def replay():
    start = time.perf_counter()
    TICKS.mkdir(exist_ok=True)
    tasks = [(s, c, names) for s, names in SOURCES.items() for c in CC.cases()]
    with ProcessPoolExecutor(8) as ex:
        metas = list(ex.map(case_job, tasks, chunksize=1))
    closed_rows = {(r['method'], r['index']): r for r in json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))['rows']}
    worst = max(max(abs(m['source_L']-closed_rows[(m['source'], m['index'])]['L']),
                    abs(m['source_rmse']-closed_rows[(m['source'], m['index'])]['force_rmse']),
                    abs(m['source_peak']-closed_rows[(m['source'], m['index'])]['peak'])) for m in metas)
    assert worst == 0., ('Source trajectory differs from the closed-loop row', worst)
    print('replay done', round(time.perf_counter()-start), 's; source regression', worst,
          'command replay', max(m['command_replay_error'] for m in metas), flush=True)


def load_case(path):
    with np.load(path) as z:
        data = {k: z[k] for k in z.files}
    meta = json.loads(str(data.pop('meta')))
    source_log = data.pop('source_log'); cols = data.pop('source_cols')
    log = np.full((len(source_log), P3.LOGC), np.nan)
    log[:, cols] = source_log
    return meta, log, data


def score_case(path):
    meta, log, data = load_case(path)
    t = log[:, 0]
    sc = dict(event_info=meta['event_info'])
    _, masks = SE.segment_masks(log, sc)
    contact = SE.contact_mask(log)
    ell = np.full(len(t), np.nan); ell[contact] = np.log(log[contact, SE.K_REF])
    row = dict({k: meta[k] for k in ('source', 'index', 'group', 'cell', 'seed', 'source_L', 'source_peak')},
               segment_ticks={s: int(masks[s].sum()) for s in SEGMENTS},
               contact_lost_after_1s=float(np.mean(~contact[t >= 1.])), methods={})
    samples = {}
    for name in meta['methods']:
        raw, sigma = data[name+'__raw'].astype(float), data[name+'__sigma'].astype(float)
        valid, ok = data[name+'__valid'], data[name+'__ok']
        if name in UNGATED:
            ok = valid  # k_good follows every valid output; the host trace still marks floor values as imprecise
        good = data[name+'__good'].astype(float)
        # raw: only ticks where the estimator produced a value, no hold. k_good is what
        # the controller uses at each tick, so it holds by construction.
        raw = np.where(valid, raw, np.nan)
        with np.errstate(divide='ignore', invalid='ignore'):
            err = dict(raw=np.log(raw)-ell, k_good=np.log(good)-ell)
            window = err['raw']
        scored = R.score(log, np.c_[t, raw, good], sc, names=CHANNELS, keep_after_loss=True)
        m = dict(call_us_p50=float(np.percentile(data[name+'__call_us'], 50)),
                 call_us_p99=float(np.percentile(data[name+'__call_us'], 99)), segments={},
                 transitions=[{k: e[k] for k in ('t0', 't1', 'delta', 'direction', 'size_class', 'overlaps_event',
                                                 'after_contact_loss', 'contact_loss_in_interval', 'channels')}
                              for e in scored['transitions']])
        samples[name] = {}
        for s in SEGMENTS:
            mask = masks[s]
            samples[name][s] = {c: err[c][mask] for c in CHANNELS}
            if not mask.any():
                continue
            used = mask & ok & valid
            e = window[used]
            sg = sigma[used]/raw[used]
            cell = dict(ticks=int(mask.sum()), valid=float(valid[mask].mean()), accepted=float(ok[mask].mean()),
                        usable=float(used.sum()/mask.sum()),
                        k_good_available=float(np.isfinite(good[mask]).mean()),
                        within_sd={c: SE.distribution(err[c][mask])['robust_sd'] for c in CHANNELS},
                        accepted_ticks=int(used.sum()))
            if used.any():
                cell.update(false_update=float(np.mean(abs(e) > LOG12)), low_update=float(np.mean(e < -LOG12)),
                            coverage=float(np.mean(abs(e) <= 2*sg)))
            m['segments'][s] = cell
        row['methods'][name] = m
    return row, samples


def km_median(times, observed):
    """Kaplan-Meier median with right censoring; None if survival stays above one half.

    Censored events stay in the risk set up to their censoring time, so slow or
    unrecovered events cannot be dropped and cannot pull the median down.
    """
    order = sorted(zip(times, [not o for o in observed]))
    at_risk, survival = len(order), 1.
    for time_, censored in order:
        if not censored:
            survival *= 1.-1./at_risk
            if survival <= .5:
                return float(time_)
        at_risk -= 1
    return None


def pooled(x):
    finite = x[np.isfinite(x)]
    if not len(finite):
        return dict(n=0, missing=len(x))
    a = abs(finite)
    return dict(n=len(finite), missing=int(len(x)-len(finite)), abs_median=float(np.median(a)), abs_p90=float(np.percentile(a, 90)),
                bias=float(np.median(finite)), below_085=float(np.mean(finite < LOG085)), worst_low=float(finite.min()))


def wmean(cells, field):
    v = [(c[field], c['ticks' if field not in ('false_update', 'low_update', 'coverage') else 'accepted_ticks'])
         for c in cells if c.get(field) is not None]
    w = sum(x[1] for x in v)
    return float(sum(x[0]*x[1] for x in v)/w) if w else None


def summarize(rows, samples, source, names):
    idx = [i for i, r in enumerate(rows) if r['source'] == source]
    segments, tracking, cost = [], [], []
    for group in CC.SCENES+('all',):
        use = [i for i in idx if group == 'all' or rows[i]['group'] == group]
        for name in names:
            for s in SEGMENTS:
                cells = [rows[i]['methods'][name]['segments'][s] for i in use if s in rows[i]['methods'][name]['segments']]
                if not cells:
                    continue
                item = dict(source=source, group=group, method=name, segment=s, ticks=sum(c['ticks'] for c in cells), cases=len(cells))
                for c in CHANNELS:
                    item[c] = pooled(np.concatenate([samples[i][name][s][c] for i in use]))
                    sds = [x['within_sd'][c] for x in cells if x['within_sd'][c] is not None]
                    item[c]['within_case_sd'] = float(np.median(sds)) if sds else None
                for field in ('valid', 'accepted', 'usable', 'k_good_available', 'false_update', 'low_update', 'coverage'):
                    item[field] = wmean(cells, field)
                segments.append(item)
    for family, groups, coupled, size in (('no_event', CC.SCENES, False, 'large'), ('core_scenes', ('A', 'A_slow_thin', 'A+B'), False, 'large'),
                                          ('no_event_small', CC.SCENES, False, 'small'), ('event_coupled', CC.SCENES, True, 'large')):
        use = [i for i in idx if rows[i]['group'] in groups]
        for name in names:
            for direction in ('rising', 'falling'):
                events = [e for i in use for e in rows[i]['methods'][name]['transitions']
                          if e['size_class'] == size and e['direction'] == direction and e['overlaps_event'] == coupled]
                item = dict(source=source, family=family, method=name, direction=direction, n=len(events),
                            loss_related=sum(e['after_contact_loss'] or e['contact_loss_in_interval'] for e in events))
                for c in CHANNELS:
                    cc = [e['channels'][c] for e in events]
                    # Lag: reaching the 1.5 s grid bound is a censored value; undefined lag (too little data) is censored at zero.
                    lags = [(1.5, False) if x['lag_at_upper_bound'] else (0., False) if x['lag'] is None else (x['lag'], True) for x in cc]
                    # Recovery: censored at 2 s, or earlier when the episode ends.
                    recov = [(x['recovery_raw']['observed_horizon'], False) if x['recovery_raw']['censored'] else (x['recovery_raw']['time'], True) for x in cc]
                    truncated = sum(x['recovery_raw']['censored'] and x['recovery_raw']['observed_horizon'] < 2.-1e-8 for x in cc)
                    unrecovered = sum(x['recovery_raw']['censored'] for x in cc)-truncated
                    item[c] = dict(lag_median=km_median(*zip(*lags)) if lags else None, lag_censored=sum(x['lag_at_upper_bound'] for x in cc),
                        lag_undefined=sum(x['lag'] is None for x in cc),
                        under_raw=R.dist([x['raw_under_magnitude'] for x in cc]), no_estimate=sum(x['raw_under_magnitude'] is None for x in cc),
                        over_raw=R.dist([x['raw_over_magnitude'] for x in cc]),
                        under_local=R.dist([x['under_magnitude'] for x in cc]),
                        under_fraction=R.dist([x['under_fraction'] for x in cc]),
                        missing_fraction=R.dist([x['missing_fraction'] for x in cc]),
                        recovery_median=km_median(*zip(*recov)) if recov else None,
                        recovered=sum(not x['recovery_raw']['censored'] for x in cc), unrecovered=unrecovered, truncated=truncated)
                tracking.append(item)
    for name in names:
        cost.append(dict(source=source, method=name,
            call_us_p50=float(np.median([rows[i]['methods'][name]['call_us_p50'] for i in idx])),
            call_us_p99=float(np.median([rows[i]['methods'][name]['call_us_p99'] for i in idx]))))
    return segments, tracking, cost


def score():
    start = time.perf_counter()
    paths = sorted(TICKS.glob('*.npz'))
    assert len(paths) == sum(36 for _ in SOURCES), len(paths)
    with ProcessPoolExecutor(8) as ex:
        results = list(ex.map(score_case, paths, chunksize=1))
    rows = [r for r, _ in results]; samples = [s for _, s in results]
    segments, tracking, cost = [], [], []
    for source, names in SOURCES.items():
        a, b, c = summarize(rows, samples, source, names)
        segments += a; tracking += b; cost += c
    B.dump('identification_replay.json', dict(
        metadata=dict(reference='k_eff, log[:,26]; contact from normal force and tissue gain',
            cases='corrected_comparison.cases(); observations recorded once per source; every method sees identical rows and no truth',
            segments='plateau, transition: no scene event; rate of ln k_eff (0.2 s boxcar) below .05/s away from transitions, above .2/s; between them unassigned. event: every tick in [start, end+.3]; event_only and event_transition split it by whether the reference is also changing.',
            channels='raw: estimator output on ticks where it produced a value, nothing held. k_good: the stiffness the controller uses at each tick, held by construction; Huber-D4, D4 and P5 apply their acceptance gates and CENTRE 1.054, the others use the estimate directly. Memory methods: raw is the filter output.',
            rates='valid: raw finite and positive; accepted: per host rule, equal to valid for methods without gates; usable: both; false_update: |ln(k/k_eff)|>ln1.2 among ticks that update k_good; low_update: its low side; coverage: |err|<=2 sigma among those ticks.',
            tracking='Transitions by family. no_event: all scenes, |d ln k|>=.3, not overlapping an event. core_scenes: its A, A_slow_thin, A+B part. no_event_small: .15<=|d ln k|<.3. event_coupled: all scenes, |d ln k|>=.3, overlapping an event. Events after contact loss are kept and flagged. Lag on pre-transition-debiased error, needing data on half the interval; under_raw = max(0, -min err), never negative; over_raw the high side. Recovery within ln1.2 for .2 s from the transition end. Lag and recovery medians are Kaplan-Meier: reaching the bound, not recovering in 2 s and episode end are right-censored, not capped values. None means the median lies beyond the censoring horizon.',
            limits='Open-loop replay: an estimator does not influence the trajectory it sees. Configurations frozen, selected before the H correction. Six cases per scene, descriptive only. Timings are Python prototype under parallel workers, not real-time evidence.'),
        segments=segments, tracking=tracking, cost=cost, cases=rows, seconds=time.perf_counter()-start))
    print('score done', round(time.perf_counter()-start), 's', flush=True)


LABELS = {'Simple-1': 'Simple-1 (F=ky)', 'Simple-2': 'Simple-2 (F=ky+b)', 'Xue-const': 'Xue',
          'oracle_keff_off': 'Oracle (1.0 k_eff)', 'old_K_x1.2_margin': 'Old |K|×1.2 margin',
          'Huber-D4_idealcal': 'Huber-D⁴，理想标定', 'oracle_keff_idealcal': 'Oracle (1.0 k_eff)，理想标定'}


def f(x, d=3):
    return '—' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'


def pct(x, d=0):
    return '—' if x is None else f'{100*x:.{d}f}%'


def table1(rows):
    stored = stored_rows()
    get = lambda m: [r for r in rows if r['method'] == m]
    head = ['方法', 'L 均值', 'L 中位', 'L A', 'L A慢薄', 'L B', 'L A+B', 'L 咳嗽', 'L 综合', 'RMSE 中位 (N)', '最大力 (N)', '>15 N', '>36.5 N', '丢失时间占比', '长失联', '调用 p50/p99 (ms)']
    out = []

    def line(m):
        rr = get(m)
        L = lambda g: np.mean([r['L'] for r in rr if g is None or r['group'] == g])
        return '| ' + ' | '.join([LABELS.get(m, m), f(float(L(None))), f(float(np.median([r['L'] for r in rr])))] + [f(float(L(g))) for g in CC.SCENES] + [
            f(float(np.median([r['force_rmse'] for r in rr]))), f(max(r['peak'] for r in rr), 1),
            pct(np.mean([r['peak'] > 15. for r in rr])), pct(np.mean([r['peak'] > 36.5 for r in rr])),
            pct(float(np.mean([r['loss_time'] for r in rr]))/9., 2), pct(np.mean([r['long_loss'] for r in rr])),
            f"{np.median([r['call_p50'] for r in rr])/1e3:.2f}/{np.median([r['call_p99'] for r in rr])/1e3:.2f}"]) + ' |'
    out += ['| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
    out += ['| **相同控制器：Method 及其前身（自带门控）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in ('Huber-D4', 'D4', 'P5')]
    out += ['| **相同控制器：其他辨识方法（无 Method 门控）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in ('Rossi-ID', 'Kalman', 'MHE', 'Xue-const', 'Simple-1', 'Simple-2')]
    out += ['| **完整文献控制方案（各自架构）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in FULL]
    out += ['| **诊断参考，不参与排名** |'+' |'*(len(head)-1), line('oracle_keff_off'), line('oracle_keff_idealcal'), line('Huber-D4_idealcal')]
    margin = [r for (m, _), r in stored.items() if m == CC.MARGIN]
    out += ['| '+LABELS[CC.MARGIN]+' | '+f(float(np.mean([r['L'] for r in margin])))+' | '+f(float(np.median([r['L'] for r in margin])))+' | '+' | '.join(
        f(float(np.mean([r['L'] for r in margin if r['group'] == g]))) for g in CC.SCENES)+' | '+
        f(float(np.median([r['force_rmse'] for r in margin])))+' | '+f(max(r['peak'] for r in margin), 1)+' | '+
        pct(np.mean([r['peak'] > 15. for r in margin]))+' | '+pct(np.mean([r['peak'] > 36.5 for r in margin]))+' | 未存 | 未存 | — |']
    return out


def table2(seg):
    head = ['方法', '通道', '|误差| 中位', '|误差| P90', '偏差', '组内波动', '低于 0.85 占比', '有效率', '接受率', '更新率', '错误更新', '覆盖率']
    out = []
    for s, title in (('plateau', '平稳段'), ('transition', '刚度变化段'), ('event', '突发扰动段')):
        ticks = next(x['ticks'] for x in seg if x['segment'] == s and x['method'] == MAIN[0] and x['source'] == 'Huber-D4' and x['group'] == 'all')
        out += ['', f'**{title}**（{ticks} 拍）', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')
            for c, label in (('raw', '原始'), ('k_good', '控制用')):
                p = x[c]
                gated = m not in UNGATED
                shared = [pct(x['valid']), pct(x['accepted']) if gated else '无门控', pct(x['usable']), pct(x['false_update'], 1), pct(x['coverage'])] if c == 'raw' else ['—', '—', '—', '—', '—']
                out.append('| '+' | '.join([LABELS.get(m, m) if c == 'raw' else '', label, f(p.get('abs_median')), f(p.get('abs_p90')), f(p.get('bias')),
                                            f(p.get('within_case_sd')), pct(p.get('below_085'))]+shared)+' |')
    return out


def table2_scene(seg):
    out = []
    for s, title, groups in (('transition', '刚度变化段', ('A', 'A_slow_thin', 'A+B', 'cough', 'combo')), ('event', '突发扰动段', ('B', 'A+B', 'cough', 'combo'))):
        head = ['方法']+list(groups)
        out += ['', f'**{title}，原始输出 |误差| 中位 / P90**', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            cells = []
            for g in groups:
                x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == g)
                cells.append(f"{f(x['raw'].get('abs_median'))} / {f(x['raw'].get('abs_p90'))}")
            out.append('| '+' | '.join([LABELS.get(m, m)]+cells)+' |')
    return out


def beyond(x, bound):
    return f'>{bound}' if x is None else f(x, 2)


def table3(track, family):
    head = ['方法', '方向', '事件数', '滞后中位 (s)', '滞后删失/不可算', '最大低估中位', '最大低估最差', '最大高估中位', '低估占比中位', '恢复中位 (s)', '已恢复', '未恢复', '回合结束截断', '失联相关']
    out = []
    for c, title in (('raw', '原始输出'), ('k_good', '控制用 k_good')):
        out += ['', f'**{title}**', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            for d, label in (('rising', '上升'), ('falling', '下降')):
                x = next(i for i in track if i['method'] == m and i['direction'] == d and i['family'] == family and i['source'] == 'Huber-D4')
                p = x[c]
                out.append('| '+' | '.join([LABELS.get(m, m) if d == 'rising' else '', label, str(x['n']), beyond(p['lag_median'], 1.5), f"{p['lag_censored']}/{p['lag_undefined']}",
                    f(p['under_raw']['median']), f(p['under_raw']['max']), f(p['over_raw']['median']), pct(p['under_fraction']['median']), beyond(p['recovery_median'], 2),
                    str(p['recovered']), str(p['unrecovered']), str(p['truncated']), str(x['loss_related'])])+' |')
    return out


def table_memory(seg, track):
    stored = stored_rows()
    head = ['候选', '闭环 L（36 例）', '平稳段 |误差| 中位/P90', '变化段 |误差| 中位/P90', '事件段 |误差| 中位/P90', '上升滞后 (s)', '上升最大低估中位', '上升未恢复/截断']
    out = ['| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
    for m in ('Huber-D4',)+INTERNAL:
        cells = []
        for s in ('plateau', 'transition', 'event'):
            x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')['raw']
            cells.append(f"{f(x['abs_median'])} / {f(x['abs_p90'])}")
        t = next(i for i in track if i['method'] == m and i['direction'] == 'rising' and i['family'] == 'no_event' and i['source'] == 'Huber-D4')
        L = float(np.mean([r['L'] for (n, _), r in stored.items() if n == m]))
        p = t['raw']
        out.append('| '+' | '.join([m, f(L)]+cells+[beyond(p['lag_median'], 1.5), f(p['under_raw']['median']), f"{p['unrecovered']}/{p['truncated']}（共 {t['n']}）"])+' |')
    return out


def sense_ratio(sc, k, phase):
    """Closed-form gain of compensated force over aligned encoder at the probe frequency.

    The sensor reads F+m*a at offset delta. Compensation subtracts m_hat*a_hat, with
    a_hat from sent commands through the calibrated servo. Only an explanation:
    H preprocessing and online servo identification are not in this expression.
    """
    w = 2*math.pi*20.8; t2 = P3.ACTUATOR_T2
    servo = lambda tau: 1/((1+1j*w*(tau-t2))*(1+1j*w*t2))
    offset = np.exp(-1j*w*(sc['sensor_delay_ms']-sc['sensor_delay_hat'])*1e-3)
    return abs(((k*np.exp(1j*phase)-w*w*sc['tool_mass'])*offset+w*w*sc['tool_mass_hat']*servo(sc['tau_a_hat'])/servo(sc['tau_a_true'])))/k


def diagnosis_job(a):
    index, kind, calibration, good_gain, bar_gain, window, noise = a
    O.setup()
    case = CC.cases()[index]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    if calibration == 'mass':
        sc['tool_mass_hat'] = sc['tool_mass']
    elif calibration == 'servo_lag':
        sc['tau_a_hat'] = sc['tau_a_true']
    elif calibration == 'offset':
        sc['sensor_delay_hat'] = sc['sensor_delay_ms']
    elif calibration == 'all':
        CC.perfect_calibration(sc)
    elif calibration == 'no_inertia':
        sc['tool_mass'] = sc['tool_mass_hat'] = 0.
        sc['sensor_delay_ms'] = sc['sensor_delay_hat'] = 0.
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    rec, step = [], pol.step
    if kind == 'oracle':
        SE.gates_off(pol); pol.centre = False
        inner = pol.estimator.estimate; history = []

        def estimate(o, n, policy):
            inner(o, n, policy)
            k = float(o.truth_keff)
            if not k > 0.:
                return np.nan, np.inf, np.inf
            history.append(k)
            pol.estimator.fit = dict(pol.estimator.fit, info=1.)
            return bar_gain*float(np.mean(history[-window:])), 0., 0.

        pol.estimator.estimate = estimate
        original = pol._stiffness_est

        def separate(o, n):
            out = original(o, n)
            if good_gain != bar_gain and np.isfinite(pol.k_good):
                pol.k_good *= good_gain/bar_gain
            return out

        pol._stiffness_est = separate

    def capture(o):
        u = step(o)
        rec.append((pol.k_good, pol.kbar, getattr(pol, 'k_hat', np.nan)))
        return u

    pol.step = capture
    log = R.oracle_run(sc, pol, case['seed']) if kind == 'oracle' else P3.episode(sc, pol, case['seed'], sensor=noise)
    metric = SE.M.metric(log, pol, sc)
    w = log[:, 0] >= 1.
    x = log[w, 18]-8.
    spectrum = abs(np.fft.rfft((x-x.mean())*np.hanning(len(x))))
    freq = np.fft.rfftfreq(len(x), .008)
    spectrum[abs(freq-20.8) < .4] = 0.
    rec = np.asarray(rec); k = log[:, SE.K_REF]
    return dict(index=index, group=case['group'], kind=kind, calibration=calibration, k_good_gain=good_gain, k_bar_gain=bar_gain, window=window, noise=noise,
                L=float(np.log1p(metric['J'])), force_rmse=metric['rmse'], peak=metric['peak'],
                dominant_hz=float(freq[3+np.argmax(spectrum[3:])]), keff_median=float(np.median(k[w])),
                k_good_over_keff=float(np.nanmedian(rec[w, 0]/k[w])), k_bar_over_keff=float(np.nanmedian(rec[w, 1]/k[w])), raw_over_keff=float(np.nanmedian(rec[w, 2]/k[w])))


def oracle_diagnosis():
    """Three oscillating push cases, one that does not oscillate, one A case. Diagnostic only."""
    cases = (13, 16, 17, 14, 0)
    oracle = [('none', 1., 1., 1), ('mass', 1., 1., 1), ('servo_lag', 1., 1., 1), ('offset', 1., 1., 1), ('all', 1., 1., 1), ('no_inertia', 1., 1., 1),
              ('none', 1.2, 1., 1), ('none', 1., 1.2, 1), ('none', 1.2, 1.2, 1), ('none', 1., 1., 24)]
    huber = [('none', True), ('none', False), ('all', False)]
    jobs = [(i, 'oracle', c, g, b, n, True) for i in cases for c, g, b, n in oracle]
    jobs += [(i, 'huber', c, 1., 1., 1, noise) for i in cases for c, noise in huber]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(diagnosis_job, jobs, chunksize=1))
    ratios = []
    for c in CC.cases():
        sc = SE.MC.scene_for(tuple(c['cell']), CC.BASE, c['variant'])
        _, log, data = load_case(TICKS/f"Huber-D4_{c['index']:02d}.npz")
        m = (log[:, 0] >= 1.) & (log[:, 3] > 0.)
        k = log[m, 26]
        ratios.append(dict(index=c['index'], group=c['group'], keff=float(np.median(k)),
                           predicted=float(np.median(sense_ratio(sc, k, log[m, 27]))),
                           huber_raw=float(np.nanmedian(data['Huber-D4__raw'][m]/k)),
                           mass_error=sc['tool_mass_hat']/sc['tool_mass']-1, servo_lag_error=sc['tau_a_hat']/sc['tau_a_true']-1))
    a = np.log([[r['predicted'], r['huber_raw']] for r in ratios])
    B.dump('oracle_diagnosis.json', dict(note='Diagnostic only, on cases chosen after seeing the Oracle row.', rows=rows, sense_ratio=ratios,
                                         log_correlation=float(np.corrcoef(a[:, 0], a[:, 1])[0, 1])))
    print('oracle diagnosis done', flush=True)


def tables():
    closed_data = json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))
    replay_data = json.loads((HERE/'results/identification_replay.json').read_text(encoding='utf-8'))
    diagnosis = json.loads((HERE/'results/oracle_diagnosis.json').read_text(encoding='utf-8'))
    seg, track = replay_data['segments'], replay_data['tracking']
    sens = []
    for s in ('plateau', 'transition', 'event'):
        a = [next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')['raw']['abs_median'] for m in MAIN]
        b = [next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Rossi-ID' and i['group'] == 'all')['raw']['abs_median'] for m in MAIN]
        sens.append(f"{s} {float(np.corrcoef(np.argsort(np.argsort(a)), np.argsort(np.argsort(b)))[0, 1]):.2f}")
    text = ['# 三张主要对比表', '',
        '由 `identification_completion.py --stage tables` 从保存的数据生成，不做仿真。配置冻结自 `baseline_test_r2.json`，是 H 修正前且带 Method 门控时选定的，本轮未重调。36 个配对案例：A、A 慢薄、B、A+B、咳嗽、综合各 6 例（A=呼吸＋跨肋骨，B=呼吸＋快速推拉）。每场景样本少，只作描述，不作统计结论。', '',
        '**门控口径。** Huber-D⁴、D⁴、P5 是本方法及其前身，保留设计时的接受门控和 CENTRE。Rossi-ID、Kalman、MHE、Xue、Simple 不用这些门控和前馈跳变规则，估计直接进入控制；保留各自机制（如 Xue 的死区）、有限值与接触检查、增益上下限，以及共用反馈控制器中 k_bar 的边界与下降速率限制。Xue-chi 只比 Xue-const 多一个本项目加的残差门控，已不单列。', '',
        '## 表 1：总体闭环控制性能', '',
        'L = ln(1+J)。辨识方法用同一控制器；完整文献方案保持各自架构，两类不直接排名。超限比例按回合计；丢失时间占比按 t≥1 s 的 9 s 计。调用时间为 Python 原型，不是实时性证据。Xue 输出切线刚度，闭环中被当作控制增益使用，这是接口适配。', '']
    text += table1(closed_data['rows'])
    rows = diagnosis['rows']

    def cell(i, calibration, g=1., b=1., n=1):
        r = next(x for x in rows if x['kind'] == 'oracle' and x['index'] == i and x['calibration'] == calibration and x['k_good_gain'] == g and x['k_bar_gain'] == b and x['window'] == n)
        return f"{r['L']:.3f}"
    text += ['', f"门控在真值注入下不起作用：`oracle_keff_on` 与 `oracle_keff_off` 的 L 最大差 {closed_data['metadata']['oracle_gate_effect_max_abs_dL']:.2e}，表中只列一行。", '',
        '**Oracle 比 Huber-D⁴ 差的原因。** 真值注入本身正确：注入值与日志差 <1e-8，配对回放实测的探测频率增益与 `log[:,26]` 的中位比在 36 例中为 0.974–1.002。问题在于控制回路用的不是接触力，而是经过探头惯性补偿的力。软组织上惯性项（约 6700 N/m）远大于 k_eff（约 700 N/m），伺服滞后、质量等标定误差留下的残差明显改变回路增益，所以回路实际增益高于 k_eff。Huber-D⁴ 从同一补偿后的力估计，估计值跟随回路实际增益；Oracle 注入物理 k_eff，k_good 反而偏低，回路比值超过前馈容许上界（约 1.23–1.29），出现 17–19 Hz 振荡。L 对比：', '',
        '| 案例 | k_eff (N/m) | 伺服滞后误差 | Huber 原始/k_eff（无噪声） | 同上·理想标定 | Oracle | 只修质量 | 只修伺服 | 只修时间偏移 | 全部修正 | 无惯性 | k_good×1.2 | k_bar×1.2 | 真值 24 拍平均 |', '|'+'---|'*14]
    for i in (13, 16, 17, 14, 0):
        h = lambda c: next(x for x in rows if x['kind'] == 'huber' and x['index'] == i and x['calibration'] == c and not x['noise'])['raw_over_keff']
        q = next(r for r in diagnosis['sense_ratio'] if r['index'] == i)
        text.append(f"| {i} | {q['keff']:.0f} | {100*q['servo_lag_error']:+.1f}% | {h('none'):.2f} | {h('all'):.2f} | {cell(i, 'none')} | {cell(i, 'mass')} | {cell(i, 'servo_lag')} | {cell(i, 'offset')} | {cell(i, 'all')} | {cell(i, 'no_inertia')} | {cell(i, 'none', 1.2, 1.)} | {cell(i, 'none', 1., 1.2)} | {cell(i, 'none', n=24)} |")
    text += ['', f"所以这不是真值注入的实现错误，也不是控制器普遍需要增益裕度：标定准确时 1.0 k_eff 不振荡，提高增益只是抵消了这个偏差。把 k_good 换成 1.2 倍比换 k_bar 更有效，说明关键在前馈回路。惯性残差的闭式式（不含 H 和在线伺服辨识）与 Huber 原始/k_eff 在 36 例上的对数相关为 {diagnosis['log_correlation']:.2f}。", '',
        '## 表 2：刚度辨识精度与可靠性', '',
        '相同观测（Huber-D⁴ 闭环记录），真值 k_eff，误差为 ln(k̂/k_eff)。原始：只统计估计器实际给出值的拍，不保持历史值；有效率单列。控制用：控制器每拍实际使用的 k_good，按定义保持上次更新值，Huber-D⁴、D⁴、P5 含 CENTRE 1.054。组内波动为每例稳健标准差的中位数。接受率、更新率、错误更新（更新拍中 |误差|>ln1.2）和覆盖率（|误差|≤2σ）按原始输出计。', '',
        '**真值含义。** k_eff 是物理接触增益，所有估计器读的都是补偿后的力，因此软组织上的误差含共同的传感偏差（见上）。同一观测下这对各方法相同，比较仍公平，但误差数值不能全归于估计器。', '',
        '**适配边界。** 表中文献方法都是原文机制在本任务中的适配，不是原论文的完整复现，见 `baseline_protocol.md`。Kalman 与 MHE 的冻结配置过程噪声相同，同一观测下原始估计完全相同（线性高斯下 MHE 退化为 Kalman），不能算两份独立证据。Xue 输出切线刚度，不进入表 2、表 3：换算到 k_eff 需要探测频率动态项、串联柔度和身体支撑，还需要不受身体运动影响的绝对压入深度，这些不用真值拿不到；实测每例 Xue/k_eff 中位比在 0.17–5.8 之间，不存在可用的固定换算。']
    text += table2(seg)
    text += table2_scene(seg)
    text += ['', f'换用 Rossi-ID 闭环作记录源时，8 种方法按原始 |误差| 中位数排序的秩相关：{"，".join(sens)}。', '',
        '## 表 3：刚度变化跟踪与恢复能力', '',
        '与表 2 同一次回放。大过渡 |Δln k|≥0.3，全部场景，不与事件重叠。滞后在扣除过渡前偏差后计。最大低估 = max(0, −最小误差)，不会为负；高估单列。恢复：|误差|≤ln1.2 并保持 0.2 s。', '',
        '滞后和恢复的中位数用 Kaplan–Meier 估计：滞后到 1.5 s 上界、2 s 内未恢复、回合结束，都按右删失处理，不当作实际时间，也不删除；“>1.5”“>2”表示中位数超出观察范围。滞后不可算（区间内有效数据不足一半）单独计数。已恢复、未恢复、回合结束截断三项之和等于事件数。', '',
        '**事件数少。** 换成 k_eff 后串联柔度压缩了肋骨处的增益变化，共同轨迹上的大过渡只有这些，主要来自 A 慢薄和综合，A 只有 1 个，A+B 没有。0.15≤|Δln k|<0.3 的小过渡见 `identification_replay.json` 的 `no_event_small`。']
    text += table3(track, 'no_event')
    text += ['', '### 表 3 附：与突发事件重叠的刚度变化', '', '全部场景，包含咳嗽活化和推拉改变工作深度。不计入上表，单独列出。']
    text += table3(track, 'event_coupled')
    text += ['', '## 补充：内部记忆候选（不是文献基线）', '', 'Huber-D⁴ 窗口输出关门控后接各记忆滤波。闭环 L 来自 `corrected_comparison.json`，其余来自同一次回放。', '']
    text += table_memory(seg, track)
    (HERE/'comparison_tables.md').write_text('\n'.join(text)+'\n', encoding='utf-8')
    print('tables written', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=('closed', 'replay', 'score', 'oracle_diagnosis', 'tables'), required=True)
    args = ap.parse_args()
    O.setup()
    globals()[args.stage]()
