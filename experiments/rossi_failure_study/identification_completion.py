"""Three comparison tables on one scoring rule, reference k_eff; no tuning.

Run in order; each stage reads only what the previous one saved:
  python -X utf8 identification_completion.py --stage closed   # results/completion_closed_loop.json
  python -X utf8 identification_completion.py --stage replay   # results/id_replay_ticks/*.npz
  python -X utf8 identification_completion.py --stage score    # results/identification_replay.json
  python -X utf8 identification_completion.py --stage oracle_diagnosis  # results/oracle_diagnosis.json
  python -X utf8 identification_completion.py --stage tables   # comparison_tables.md

closed: corrected_comparison.cases() (42 paired cases). Only methods whose control
changed run again; other rows in completion_closed_loop.json are kept.

replay: every identifier on identical recorded observations of the same 42 cases.
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
MAIN = ('Huber-D4', 'Huber-D4_nogate', 'Huber-D4_imm', 'D4', 'D4_nogate', 'D4_imm', 'P5', 'Rossi-ID', 'Kalman', 'MHE', 'Simple-1', 'Simple-2')
CONTROL = MAIN + ('Xue-const',)
INTERNAL = tuple(m for m in CC.MEMORY if m != 'imm')
SOURCES = {'Huber-D4': MAIN + INTERNAL, 'Rossi-ID': MAIN}
FULL = ('Rossi-control', 'DO-MPAC', 'Xing-DAHIC')
UNGATED = CC.OWN_MECHANISM + CC.ARMS + INTERNAL
RERUN = CC.ARMS
TABLE1_METHODS = CONTROL + FULL + ('oracle_keff_off',) + CC.IDEAL_CAL + (CC.MARGIN,)
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
    kept = previous['rows']
    present = {(r['method'], r['index']) for r in kept}
    methods = TABLE1_METHODS
    tasks = [(m, c) for m in methods for c in CC.cases()
             if (m, c['index']) not in present and (m in RERUN or c['index'] >= 36)]
    print('closed tasks', len(tasks), 'workers 8', flush=True)
    with ProcessPoolExecutor(8) as ex:
        fresh = []
        for (m, _), r in zip(tasks, ex.map(closed_job, tasks, chunksize=1)):
            fresh.append(dict(r, method=m))
            B.dump('completion_closed_loop.json', dict(previous, rows=kept+fresh))
            if len(fresh) % 24 == 0:
                print('closed', len(fresh), '/', len(tasks), flush=True)
    # Sanity on kept rows: gates-on Huber must still match its stored 612-grid values.
    old = stored_rows()
    worst = max(max(abs(r[k]-old[('Huber-D4', r['index'])][k]) for k in ('L', 'force_rmse', 'peak')) for r in kept if r['method'] == 'Huber-D4' and r['index'] < 36)
    assert worst == 0.
    rows = kept+fresh
    now = {(r['method'], r['index']): r for r in rows}
    gate = [abs(now[('oracle_keff_on', c['index'])]['L']-now[('oracle_keff_off', c['index'])]['L']) for c in CC.cases()[:36]]
    metadata = dict(previous['metadata'])
    metadata.update(
        own_mechanism=', '.join(UNGATED)+': CC.own_mechanism; no acceptance gates, phase sigma gate, feed-forward jump rule or CENTRE. _imm adds unchanged R.make_filter(imm). Huber-D4, D4 and P5 keep their gates.',
        cases=42, reused_old_cases=True,
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
    host = name.rsplit('_', 1)[0] if name in CC.ARMS else 'Huber-D4' if name in INTERNAL or name.startswith('Simple') else name
    cfg = O.chosen()[host]
    pol = B.build(host, cfg, sc)
    if name.startswith('Simple'):
        pol.estimator = I.ContactAdapter(S.SimpleLS(name == 'Simple-2', nw=24, sigma_scale=cfg['sigma_scale'], chi_scale=cfg['gate']))
    if name in UNGATED:
        CC.own_mechanism(pol)
    if name in INTERNAL or name.endswith('_imm'):
        R.make_filter(pol, 'imm' if name.endswith('_imm') else name)
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
                                                 'after_contact_loss', 'contact_loss_in_interval', 'next_t0',
                                                 'raw_window_end', 'raw_interval_ticks', 'channels')}
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
    for family, groups, coupled, size in (('no_event', CC.SCENES[:-1], False, 'large'), ('core_scenes', ('A', 'A_slow_thin', 'A+B'), False, 'large'),
                                          ('no_event_small', CC.SCENES[:-1], False, 'small'), ('A_large_no_event', ('A_large',), False, 'large'),
                                          ('event_coupled', CC.SCENES, True, 'large')):
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
                    # Recovery: stop at 2 s, episode end, or the next detected change.
                    recov = [(x['recovery_raw']['observed_horizon'], False) if x['recovery_raw']['censored'] else (x['recovery_raw']['time'], True) for x in cc]
                    reasons = [x['recovery_raw']['censored_by'] for x in cc]
                    recovered = sum(not x['recovery_raw']['censored'] for x in cc)
                    unrecovered = reasons.count('horizon')
                    truncated = reasons.count('episode_end')
                    truncated_next = reasons.count('next_change')
                    assert recovered+unrecovered+truncated+truncated_next == len(cc)
                    item[c] = dict(lag_median=km_median(*zip(*lags)) if lags else None, lag_censored=sum(x['lag_at_upper_bound'] for x in cc),
                        change_mae=R.dist([x['change_mae'] for x in cc]),
                        lag_undefined=sum(x['lag'] is None for x in cc),
                        under_raw=R.dist([x['raw_under_magnitude'] for x in cc]), no_estimate=sum(x['raw_under_magnitude'] is None for x in cc),
                        over_raw=R.dist([x['raw_over_magnitude'] for x in cc]),
                        under_local=R.dist([x['under_magnitude'] for x in cc]),
                        under_fraction=R.dist([x['under_fraction'] for x in cc]),
                        missing_fraction=R.dist([x['missing_fraction'] for x in cc]),
                        recovery_median=km_median(*zip(*recov)) if recov else None,
                        recovered=recovered, unrecovered=unrecovered, truncated=truncated,
                        truncated_next_change=truncated_next)
                tracking.append(item)
    for name in names:
        cost.append(dict(source=source, method=name,
            call_us_p50=float(np.median([rows[i]['methods'][name]['call_us_p50'] for i in idx])),
            call_us_p99=float(np.median([rows[i]['methods'][name]['call_us_p99'] for i in idx]))))
    return segments, tracking, cost


def score():
    start = time.perf_counter()
    paths = sorted(TICKS.glob('*.npz'))
    assert len(paths) == sum(42 for _ in SOURCES), len(paths)
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
            tracking='Transitions by family. no_event: existing 36 scenes, |d ln k|>=.3, not overlapping an event. core_scenes: its A, A_slow_thin, A+B part. A_large_no_event: A_large only, large, no event. no_event_small: existing scenes, .15<=|d ln k|<.3. event_coupled: all scenes, |d ln k|>=.3, overlapping an event. Events after contact loss are kept and flagged. Lag on pre-transition-debiased error, needing data on half the interval; change_mae: mean absolute log error on valid ticks in [t0,t1]; under_raw = max(0, -min err), never negative; over_raw the high side. Raw under/over/fractions use [t0, min(t1+.5, next_t0)), where next_t0 includes kept and excluded transitions. Signed runs: rate>.2 rising and rate<-.2 falling, same-sign gaps<.1s merge only across continuous contact and without an intervening opposite run; abs(delta)>=.15. Recovery within ln1.2 for .2 s from the transition end; hold must complete before next_t0/episode end. censored_by: next_change, episode_end or horizon; zero horizons are retained. truncated denotes episode_end; truncated_next_change is separate. Lag and recovery medians are Kaplan-Meier: reaching the bound, not recovering in 2 s, next change and episode end are right-censored, not capped values. None means the median lies beyond the censoring horizon.',
            limits='Open-loop replay: an estimator does not influence the trajectory it sees. Configurations frozen, selected before the H correction. Six cases per scene, descriptive only. Timings are Python prototype under parallel workers, not real-time evidence.'),
        segments=segments, tracking=tracking, cost=cost, cases=rows, seconds=time.perf_counter()-start))
    print('score done', round(time.perf_counter()-start), 's', flush=True)


LABELS = {'Simple-1': 'Simple-1 (F=ky)', 'Simple-2': 'Simple-2 (F=ky+b)', 'Xue-const': 'Xue',
          'oracle_keff_off': 'Oracle (1.0 k_eff)', 'old_K_x1.2_margin': 'Old |K|×1.2 margin',
          'Huber-D4_idealcal': 'Huber-D⁴，理想标定', 'oracle_keff_idealcal': 'Oracle (1.0 k_eff)，理想标定'}


LABELS.update({m: base+'（'+label+'）' for base in ('Huber-D4', 'D4') for m, label in ((base, '门控'), (base+'_nogate', '无门控'), (base+'_imm', '无门控＋IMM'))})


def f(x, d=3):
    return '—' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'


def pct(x, d=0):
    return '—' if x is None else f'{100*x:.{d}f}%'


def md_row(cells):
    return '| '+' | '.join(str(x).replace('|', r'\|') for x in cells)+' |'


def table1(rows):
    stored = stored_rows()
    get = lambda m: [r for r in rows if r['method'] == m]
    head = ['方法', 'L 均值', 'L 中位', 'L A', 'L A慢薄', 'L B', 'L A+B', 'L 咳嗽', 'L 综合', 'L A_large', 'RMSE 中位 (N)', '最大力 (N)', '>15 N', '>36.5 N', '丢失时间占比', '长失联', '调用 p50/p99 (ms)']
    out = []

    def line(m):
        rr = get(m)
        L = lambda g: np.mean([r['L'] for r in rr if g is None or r['group'] == g])
        return '| ' + ' | '.join([LABELS.get(m, m), f(float(L(None))), f(float(np.median([r['L'] for r in rr])))] + [f(float(L(g))) for g in CC.SCENES] + [
            f(float(np.median([r['force_rmse'] for r in rr]))), f(max(r['peak'] for r in rr), 1),
            pct(np.mean([r['peak'] > 15. for r in rr])), pct(np.mean([r['peak'] > 36.5 for r in rr])),
            pct(float(np.mean([r['loss_time'] for r in rr]))/9., 2), pct(np.mean([r['long_loss'] for r in rr])),
            f"{np.median([r['call_p50'] for r in rr])/1e3:.2f}/{np.median([r['call_p99'] for r in rr])/1e3:.2f}"]) + ' |'
    out += [md_row(head), '|'+'---|'*len(head)]
    out += ['| **相同控制器：Method 三臂及其前身** |'+' |'*(len(head)-1)]
    out += [line(m) for m in MAIN[:7]]
    out += ['| **相同控制器：其他辨识方法（无 Method 门控）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in ('Rossi-ID', 'Kalman', 'MHE', 'Xue-const', 'Simple-1', 'Simple-2')]
    out += ['| **完整文献控制方案（各自架构）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in FULL]
    out += ['| **诊断参考，不参与排名** |'+' |'*(len(head)-1), line('oracle_keff_off'), line('oracle_keff_idealcal'), line('Huber-D4_idealcal')]
    margin = [r for (m, _), r in stored.items() if m == CC.MARGIN]
    margin += get(CC.MARGIN)
    out.append(table1_margin(margin))
    return out


def table1_margin(rows):
    # Historical diagnostic has no stored controller timings/loss metrics on old cases.
    cells = [LABELS[CC.MARGIN], f(float(np.mean([r['L'] for r in rows]))), f(float(np.median([r['L'] for r in rows])))]
    cells += [f(float(np.mean([r['L'] for r in rows if r['group'] == g]))) for g in CC.SCENES]
    cells += [f(float(np.median([r['force_rmse'] for r in rows]))), f(max(r['peak'] for r in rows), 1),
              pct(np.mean([r['peak'] > 15. for r in rows])), pct(np.mean([r['peak'] > 36.5 for r in rows])), '未存', '未存', '—']
    return md_row(cells)


def table2(seg):
    head = ['方法', '通道', '|误差| 中位', '|误差| P90', '偏差', '组内波动', '低于 0.85 占比', '有效率', '接受率', '更新率', '错误更新', '覆盖率']
    out = []
    for s, title in (('plateau', '平稳段'), ('transition', '刚度变化段'), ('event', '突发扰动段')):
        ticks = next(x['ticks'] for x in seg if x['segment'] == s and x['method'] == MAIN[0] and x['source'] == 'Huber-D4' and x['group'] == 'all')
        out += ['', f'**{title}**（{ticks} 拍）', '', md_row(head), '|'+'---|'*len(head)]
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
    for s, title, groups in (('transition', '刚度变化段', ('A', 'A_slow_thin', 'A+B', 'cough', 'combo', 'A_large')), ('event', '突发扰动段', ('B', 'A+B', 'cough', 'combo'))):
        head = ['方法']+list(groups)
        out += ['', f'**{title}，原始输出 |误差| 中位 / P90**', '', md_row(head), '|'+'---|'*len(head)]
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
    head = ['方法', '方向', '事件数', '变化段误差中位', '最大低估中位', '最大低估最差', '恢复时间中位 (s)', '已恢复', '未恢复', '下次变化截断', '回合结束截断']
    out = []
    for c, title in (('raw', '原始输出'), ('k_good', '控制用 k_good')):
        out += ['', f'**{title}**', '', md_row(head), '|'+'---|'*len(head)]
        for m in MAIN:
            for d, label in (('rising', '上升'), ('falling', '下降')):
                x = next(i for i in track if i['method'] == m and i['direction'] == d and i['family'] == family and i['source'] == 'Huber-D4')
                p = x[c]
                recovery = '—' if x['n'] == 0 else '未达到' if p['recovery_median'] is None else f(p['recovery_median'], 2)
                out.append('| '+' | '.join([LABELS.get(m, m) if d == 'rising' else '', label, str(x['n']),
                    f(p['change_mae']['median']), f(p['under_raw']['median']), f(p['under_raw']['max']), recovery,
                    str(p['recovered']), str(p['unrecovered']), str(p['truncated_next_change']), str(p['truncated'])])+' |')
    out += ['', '滞后保留在数据中，不作为主指标，因为曲线平移没有直接的响应含义。']
    out += ['', '下一次刚度变化开始前未恢复的事件按右删失处理，最大低估窗口也在下一次变化开始时结束。']
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
    seg, track = replay_data['segments'], replay_data['tracking']
    text = ['# 三张主要对比表', '',
        '由 `identification_completion.py --stage tables` 从保存数据生成。冻结配置 O.chosen() 未重调。42 个配对案例：A、A 慢薄、B、A+B、咳嗽、综合、A_large 各 6 例。A=呼吸＋跨肋骨，B=呼吸＋快速推拉；A_large = 与 A 相同场景，探头长轴平行于肋骨，覆盖层为 18 mm。每组仅 6 例，只作描述。', '',
        '门控臂保持原样；无门控臂采用 `CC.own_mechanism`，取消接受门控、CENTRE、相位检测 sigma 门槛及前馈跳变规则，保留估计器自身机制、有限值与接触检查、增益边界及 k_bar 下降速率限制。无门控＋IMM 在此基础上加原有 `R.make_filter(pol, "imm")`。P5 仅门控。', '',
        '除 IMM 外的记忆形式仅保留在 Huber-D4 源回放数据中，不进入表格。旧 `off_none` 和旧闭环记忆行仍保存于 corrected_comparison.json，含 CENTRE 和相位检测 sigma 门槛，已被当前口径取代。', '',
        '## 表 1：总体闭环控制性能', '',
        'L=ln(1+J)。辨识方法使用同一控制器；完整文献控制方案各自保留架构。Oracle 为诊断参考。旧 36 例复用保存结果；新增臂跑 42 例，既有方法补 6 例 A_large。超限比例按回合计，丢失时间占比按 t≥1 s 的 9 s 计。Python 调用时间不是实时性证据。', '']
    text += table1(closed_data['rows'])
    text += ['', '## 表 2：刚度辨识精度与可靠性', '',
        '共同观测来自 Huber-D4 闭环，参考为物理 k_eff（log[:,26]），误差为 ln(k̂/k_eff)。原始输出只在估计器实际给值的拍统计；k_good 为控制器每拍实际使用值，含历史保持。门控臂与 P5 保留 CENTRE 1.054。接受率、更新率、错误更新（更新拍 |误差|>ln1.2）和覆盖率（|误差|≤2σ）按原始输出计。组内波动为每例稳健标准差的中位数。', '',
        '估计器使用惯性补偿后的力，因此相对于物理 k_eff 的误差包含共同传感偏差。文献方法是原文机制在本任务中的适配，见 baseline_protocol.md。冻结 Kalman 与 MHE 在线性高斯情形下原始估计相同。Xue 输出切线刚度，转换为 k_eff 需要额外不可得信息，不进入表 2、3。既有 Oracle 标定诊断仍见 results/oracle_diagnosis.json，未追加该类实验。']
    text += table2(seg)
    text += table2_scene(seg)
    text += ['', '## 表 3：刚度变化跟踪与恢复能力', '',
        '与表 2 同次回放，以下三族均不与场景事件重叠，A_large 不进入既有场景两族。上升和下降分别检测，同方向短间隙仅在没有反向变化且连续接触时合并。变化段误差是每个事件在 [t0,t1] 有效拍的平均 |ln(k̂/k_eff)|，再取事件间中位数。最大低估取 [t0,min(t1+0.5 s,next_t0)) 的原始误差 max(0,−min ln(k̂/k_eff))，不扣偏差；next_t0 包括保留和排除的下一过渡。恢复从变化结束算起，|误差|≤ln1.2 连续保持 0.2 s，观察上限 2 s，保持段必须在下一变化前完成。中位数按 Kaplan–Meier 估计，“未达到”表示删失数据未能识别中位数，不等于恢复时间大于 2 s。已恢复、未恢复（完整观察 2 s）、下次变化截断、回合结束截断四项之和等于事件数；零观察窗也计为删失。缺失估计及失联标记保留于 JSON。事件耦合结果仅保留于 JSON。']
    for family, title in (('no_event', '既有场景，|Δln k|≥0.3'), ('no_event_small', '既有场景，0.15≤|Δln k|<0.3'), ('A_large_no_event', 'A_large，|Δln k|≥0.3')):
        text += ['', '### '+title]
        text += table3(track, family)
    (HERE/'comparison_tables.md').write_text('\n'.join(text)+'\n', encoding='utf-8')
    print('tables written', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=('closed', 'replay', 'score', 'oracle_diagnosis', 'tables'), required=True)
    args = ap.parse_args()
    O.setup()
    globals()[args.stage]()
