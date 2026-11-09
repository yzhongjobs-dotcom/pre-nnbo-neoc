"""Three comparison tables on one scoring rule, reference k_eff; no tuning.

Run in order; each stage reads only what the previous one saved:
  python -X utf8 identification_completion.py --stage closed   # results/completion_closed_loop.json
  python -X utf8 identification_completion.py --stage replay   # results/id_replay_ticks/*.npz
  python -X utf8 identification_completion.py --stage score    # results/identification_replay.json
  python -X utf8 identification_completion.py --stage oracle_margin  # results/oracle_margin_check.json
  python -X utf8 identification_completion.py --stage tables   # comparison_tables.md

closed: corrected_comparison.cases() (36 paired cases). corrected_comparison.json
lacks contact-loss fields, so its table methods run again; their L, RMSE and peak
must equal the stored values bitwise. Memory variants and the margin row are not rerun.

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
MAIN = ('Huber-D4', 'D4', 'P5', 'Rossi-ID', 'Kalman', 'MHE', 'Xue-const', 'Xue-chi', 'Simple-1', 'Simple-2')
INTERNAL = ('off_none',) + CC.MEMORY
SOURCES = {'Huber-D4': MAIN + INTERNAL, 'Rossi-ID': MAIN}
FULL = ('Rossi-control', 'DO-MPAC', 'Xing-DAHIC')
RERUN = ('Huber-D4', 'D4', 'Rossi-ID', 'Kalman', 'MHE') + FULL  # already in the 612; need contact-loss fields
NEW = ('P5', 'Xue-const', 'Xue-chi', 'Simple-1', 'Simple-2') + tuple(CC.ORACLES)
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
    tasks = [(m, c) for m in RERUN+NEW for c in CC.cases()]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(closed_job, tasks, chunksize=1))
    old = stored_rows()
    worst = 0.
    for r in rows:
        if r['method'] in RERUN:
            o = old[(r['method'], r['index'])]
            worst = max(worst, *(abs(r[k]-o[k]) for k in ('L', 'force_rmse', 'peak')))
    assert worst == 0., ('Control results moved', worst)
    simple = json.loads((HERE/'results/simple_ls_validation.json').read_text(encoding='utf-8'))['new_method_episodes']
    now = {(r['method'], r['index']): r for r in rows}
    simple_worst = max(abs(now[(s['method'], s['index'])]['L']-s['L']) for s in simple)
    assert simple_worst == 0., ('Simple methods moved', simple_worst)
    gate = [abs(now[('oracle_keff_on', c['index'])]['L']-now[('oracle_keff_off', c['index'])]['L']) for c in CC.cases()]
    B.dump('completion_closed_loop.json', dict(
        metadata=dict(cases='corrected_comparison.cases(), same scenes and sensor seeds', base=CC.BASE,
            oracle='1.0*k_eff, log[:,26], at the current tick from a shadow tick on copied state; sigma=chi=0; CENTRE off; k_bar keeps its downward rate limit (rho_d). Diagnostic, not a ranked algorithm.',
            frozen='Configurations from baseline_test_r2 selected on development scenes before the H correction, not retuned. Xue first-contact position calibrated per scene as in baseline_compare.',
            rerun='Already in corrected_comparison but lacking contact-loss fields: '+', '.join(RERUN),
            rerun_regression_max_abs=worst, simple_regression_max_abs=simple_worst,
            oracle_gate_effect_max_abs_dL=max(gate)),
        rows=rows, seconds=time.perf_counter()-start))
    print('closed done', round(time.perf_counter()-start), 's; regression', worst, simple_worst, 'oracle gate dL', max(gate), flush=True)


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
    if host.startswith('Xue'):
        sc = dict(sc, initial_contact_observed=CC.calibrate_contact(sc, case['seed']))
    cfg = O.chosen()[host]
    pol = B.build(host, cfg, sc)
    if name.startswith('Simple'):
        pol.estimator = I.ContactAdapter(S.SimpleLS(name == 'Simple-2', nw=24, sigma_scale=cfg['sigma_scale'], chi_scale=cfg['gate']))
    if name in INTERNAL:
        SE.gates_off(pol)
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
    pol = B.build(source, O.chosen()[source], sc)
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
    old = stored_rows()
    worst = max(max(abs(m['source_L']-old[(m['source'], m['index'])]['L']),
                    abs(m['source_rmse']-old[(m['source'], m['index'])]['force_rmse']),
                    abs(m['source_peak']-old[(m['source'], m['index'])]['peak'])) for m in metas)
    assert worst == 0., ('Source trajectory differs from corrected_comparison', worst)
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
        good = data[name+'__good'].astype(float)
        # A consumer of the raw window holds its last valid value.
        last = np.where(valid, np.arange(len(t)), 0)
        np.maximum.accumulate(last, out=last)
        held = np.where(valid[last], raw[last], np.nan)
        with np.errstate(divide='ignore', invalid='ignore'):
            err = dict(raw=np.log(held)-ell, k_good=np.log(good)-ell)
            window = np.log(raw)-ell
        scored = R.score(log, np.c_[t, held, good], sc, names=CHANNELS, keep_after_loss=True)
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
                    item[c] = dict(lag=R.dist([x['lag'] for x in cc]), lag_at_bound=sum(x['lag_at_upper_bound'] for x in cc),
                        no_lag=sum(x['lag'] is None for x in cc),
                        under_raw=R.dist([x['raw_under_magnitude'] for x in cc]), no_estimate=sum(x['raw_under_magnitude'] is None for x in cc),
                        under_local=R.dist([x['under_magnitude'] for x in cc]),
                        under_fraction=R.dist([x['under_fraction'] for x in cc]),
                        missing_fraction=R.dist([x['missing_fraction'] for x in cc]),
                        recovery=R.dist([x['recovery_raw']['time'] for x in cc]),
                        unrecovered=sum(x['recovery_raw']['censored'] for x in cc),
                        short_horizon=sum(x['recovery_raw']['censored'] and x['recovery_raw']['observed_horizon'] < 2.-1e-8 for x in cc),
                        settling_local=R.dist([x['settling']['time'] for x in cc]),
                        unsettled_local=sum(x['settling']['censored'] for x in cc))
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
            channels='raw: window output, last valid held. k_good: after that method\'s gates, consumed by the feed-forward; gates-on hosts apply CENTRE 1.054. Memory methods: raw is the filter output.',
            rates='valid: raw finite and positive; accepted: precise and consistent per host gates; usable: both; false_update: |ln(k/k_eff)|>ln1.2 among accepted ticks; low_update: its low side; coverage: |err|<=2 sigma among accepted ticks. Xue sigma is a fixed interface addition, not an uncertainty.',
            tracking='Transitions by family. no_event: all scenes, |d ln k|>=.3, not overlapping an event. core_scenes: its A, A_slow_thin, A+B part. no_event_small: .15<=|d ln k|<.3. event_coupled: all scenes, |d ln k|>=.3, overlapping an event. Events after contact loss are kept and flagged. Lag on pre-transition-debiased error; under_raw and recovery on raw error; recovery: within ln1.2 for .2 s, capped at 2 s, otherwise unrecovered.',
            limits='Open-loop replay: an estimator does not influence the trajectory it sees. Configurations frozen, selected before the H correction. Six cases per scene, descriptive only. Timings are Python prototype under parallel workers, not real-time evidence.'),
        segments=segments, tracking=tracking, cost=cost, cases=rows, seconds=time.perf_counter()-start))
    print('score done', round(time.perf_counter()-start), 's', flush=True)


LABELS = {'Simple-1': 'Simple-1 (F=ky)', 'Simple-2': 'Simple-2 (F=ky+b)', 'oracle_keff_off': 'Oracle (1.0 k_eff, diagnostic)',
          'old_K_x1.2_margin': 'Old |K|×1.2 margin (diagnostic)'}
LOW = math.log(1.2)


def f(x, d=3):
    return '—' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'


def pct(x, d=0):
    return '—' if x is None else f'{100*x:.{d}f}%'


def table1(rows):
    stored = stored_rows()
    get = lambda m: [r for r in rows if r['method'] == m]
    head = ['方法', 'L 全部', 'L 中位', 'L A', 'L A慢薄', 'L B', 'L A+B', 'L 咳嗽', 'L 综合', 'RMSE 中位 (N)', '最大力 (N)', '>15 N', '>36.5 N', '丢失时间占比', '长失联', '调用 p50/p99 (ms)']
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
    out += ['| **相同控制器：辨识方法** |'+' |'*(len(head)-1)]
    out += [line(m) for m in MAIN]
    out += ['| **完整文献控制方案（各自架构）** |'+' |'*(len(head)-1)]
    out += [line(m) for m in FULL]
    out += ['| **诊断参考，不参与排名** |'+' |'*(len(head)-1), line('oracle_keff_off')]
    margin = [r for (m, _), r in stored.items() if m == CC.MARGIN]
    out += ['| '+LABELS[CC.MARGIN]+' | '+f(float(np.mean([r['L'] for r in margin])))+' | '+f(float(np.median([r['L'] for r in margin])))+' | '+' | '.join(
        f(float(np.mean([r['L'] for r in margin if r['group'] == g]))) for g in CC.SCENES)+' | '+
        f(float(np.median([r['force_rmse'] for r in margin])))+' | '+f(max(r['peak'] for r in margin), 1)+' | '+
        pct(np.mean([r['peak'] > 15. for r in margin]))+' | '+pct(np.mean([r['peak'] > 36.5 for r in margin]))+' | 未存 | 未存 | — |']
    return out


def table2(seg):
    head = ['方法', '通道', '|误差| 中位', '|误差| P90', '偏差', '组内波动', '低于 0.85 占比', '有效率', '接受率', '可用率', '错误更新', '覆盖率']
    out = []
    for s, title in (('plateau', '平稳段'), ('transition', '刚度变化段'), ('event', '突发扰动段')):
        ticks = next(x['ticks'] for x in seg if x['segment'] == s and x['method'] == MAIN[0] and x['source'] == 'Huber-D4' and x['group'] == 'all')
        out += ['', f'**{title}**（{ticks} 拍）', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')
            for c, label in (('raw', '原始'), ('k_good', '门控后')):
                p = x[c]
                shared = [pct(x['valid']), pct(x['accepted']), pct(x['usable']), pct(x['false_update'], 1), pct(x['coverage'])] if c == 'raw' else ['—', '—', pct(x['k_good_available']), '—', '—']
                out.append('| '+' | '.join([LABELS.get(m, m) if c == 'raw' else '', label, f(p.get('abs_median')), f(p.get('abs_p90')), f(p.get('bias')),
                                            f(p.get('within_case_sd')), pct(p.get('below_085'))]+shared)+' |')
    return out


def table3(track, family):
    head = ['方法', '方向', '事件数', '滞后中位 (s)', '滞后不可算或到上界', '最大低估中位', '最大低估最差', '低估占比中位', '恢复时间中位 (s)', '未恢复', '其中数据截断', '失联相关']
    out = []
    for c, title in (('raw', '原始输出'), ('k_good', '门控后 k_good')):
        out += ['', f'**{title}**', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            for d, label in (('rising', '上升'), ('falling', '下降')):
                x = next(i for i in track if i['method'] == m and i['direction'] == d and i['family'] == family and i['source'] == 'Huber-D4')
                p = x[c]; n = x['n']
                out.append('| '+' | '.join([LABELS.get(m, m) if d == 'rising' else '', label, str(n), f(p['lag']['median']), str(p['lag_at_bound']+p['no_lag']),
                    f(p['under_raw']['median']), f(p['under_raw']['max']), pct(p['under_fraction']['median']), f(p['recovery']['median'], 2),
                    f"{p['unrecovered']}/{n}", str(p['short_horizon']), str(x['loss_related'])])+' |')
    return out


def table2_scene(seg):
    out = []
    for s, title, groups in (('transition', '刚度变化段', ('A', 'A_slow_thin', 'A+B', 'cough', 'combo')), ('event', '突发扰动段', ('B', 'A+B', 'cough', 'combo'))):
        head = ['方法']+[f'{g}' for g in groups]
        out += ['', f'**{title}，原始输出 |误差| 中位 / P90**', '', '| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
        for m in MAIN:
            cells = []
            for g in groups:
                x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == g)
                cells.append(f"{f(x['raw'].get('abs_median'))} / {f(x['raw'].get('abs_p90'))}")
            out.append('| '+' | '.join([LABELS.get(m, m)]+cells)+' |')
    return out


def table_memory(seg, track):
    stored = stored_rows()
    head = ['候选', '闭环 L（36 例）', '平稳段 |误差| 中位/P90', '变化段 |误差| 中位/P90', '事件段 |误差| 中位/P90', '上升滞后 (s)', '上升最大低估中位', '上升未恢复']
    out = ['| '+' | '.join(head)+' |', '|'+'---|'*len(head)]
    for m in ('Huber-D4',)+INTERNAL:
        cells = []
        for s in ('plateau', 'transition', 'event'):
            x = next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')['raw']
            cells.append(f"{f(x['abs_median'])} / {f(x['abs_p90'])}")
        t = next(i for i in track if i['method'] == m and i['direction'] == 'rising' and i['family'] == 'no_event' and i['source'] == 'Huber-D4')
        L = float(np.mean([r['L'] for (n, _), r in stored.items() if n == m]))
        out.append('| '+' | '.join([m, f(L)]+cells+[f(t['raw']['lag']['median']), f(t['raw']['under_raw']['median']), f"{t['raw']['unrecovered']}/{t['n']}"])+' |')
    return out


def margin_job(a):
    index, good_gain, bar_gain = a
    O.setup()
    case = CC.cases()[index]
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    pol = B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    SE.gates_off(pol); pol.centre = False
    inner = pol.estimator.estimate

    def estimate(o, n, policy):
        inner(o, n, policy)
        k = float(o.truth_keff)
        if not k > 0.:
            return np.nan, np.inf, np.inf
        pol.estimator.fit = dict(pol.estimator.fit, info=1.)
        # The returned value sets k_good; k_bar is value + c*sigma.
        return good_gain*k, (bar_gain-good_gain)*k/pol.params['c'], 0.

    pol.estimator.estimate = estimate
    log = R.oracle_run(sc, pol, case['seed'])
    metric = SE.M.metric(log, pol, sc)
    w = log[:, 0] >= 1.
    x = log[w, 18]-8.
    spectrum = abs(np.fft.rfft((x-x.mean())*np.hanning(len(x))))
    freq = np.fft.rfftfreq(len(x), .008)
    return dict(index=index, group=case['group'], k_good_gain=good_gain, k_bar_gain=bar_gain, L=float(np.log1p(metric['J'])),
                force_rmse=metric['rmse'], peak=metric['peak'], dominant_hz=float(freq[3+np.argmax(spectrum[3:])]),
                keff_median=float(np.median(log[w, SE.K_REF])))


def oracle_margin():
    """Three push cases where 1.0 k_eff oscillates, one A case where it does not."""
    jobs = [(i, g, b) for i in (13, 16, 17, 0) for g, b in ((1., 1.), (1., 1.2), (1.2, 1.2))]
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(margin_job, jobs, chunksize=1))
    B.dump('oracle_margin_check.json', dict(note='Diagnostic only, chosen after seeing the Oracle row. k_good and k_bar set separately as multiples of k_eff, gates and CENTRE off.', rows=rows))
    print('oracle margin done', flush=True)


def tables():
    closed_data = json.loads((HERE/'results/completion_closed_loop.json').read_text(encoding='utf-8'))
    replay_data = json.loads((HERE/'results/identification_replay.json').read_text(encoding='utf-8'))
    seg, track = replay_data['segments'], replay_data['tracking']
    # Sensitivity to the recording host, on the main methods only.
    sens = []
    for s, c in (('plateau', 'raw'), ('transition', 'raw'), ('event', 'raw')):
        a = [next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Huber-D4' and i['group'] == 'all')[c]['abs_median'] for m in MAIN]
        b = [next(i for i in seg if i['method'] == m and i['segment'] == s and i['source'] == 'Rossi-ID' and i['group'] == 'all')[c]['abs_median'] for m in MAIN]
        ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
        sens.append(f"{s} {float(np.corrcoef(ra, rb)[0, 1]):.2f}")
    text = ['# 三张主要对比表', '',
        '由 `identification_completion.py --stage tables` 从保存的数据生成，不做仿真。方法配置冻结自 `baseline_test_r2.json`（H 修正前在开发场景选定，本轮未重调）。36 个配对案例：A、A 慢薄、B、A+B、咳嗽、综合各 6 例（A=呼吸＋跨肋骨，B=呼吸＋快速推拉）。每场景样本少，只作描述，不作统计结论。', '',
        '## 表 1：总体闭环控制性能', '',
        'L = ln(1+J)，各列为案例均值。辨识方法用同一控制器；完整文献方案保持各自架构，两类不直接排名。超限比例按回合计；丢失时间占比按 t≥1 s 的 9 s 计。调用时间为 Python 原型，不是实时性证据。', '']
    text += table1(closed_data['rows'])
    margin = json.loads((HERE/'results/oracle_margin_check.json').read_text(encoding='utf-8'))['rows']
    text += ['', f"**读表说明。** L 均值受综合工况的失控回合支配，应同时看中位数和分工况列。门控在真值注入下不起作用：`oracle_keff_on` 与 `oracle_keff_off` 的 L 最大差 {closed_data['metadata']['oracle_gate_effect_max_abs_dL']:.2e}，表中只列一行。", '',
        '**Oracle 不是性能上限。** 真值注入后，k_good 和 k_bar 都等于 k_eff，控制器其余参数不变。部分软组织回合出现持续 17–19 Hz 振荡（表中 B 的 L 因此变大），峰值不高。下面是事后诊断，只用于解释该现象：', '',
        '| 案例 | 场景 | k_eff 中位 (N/m) | k_good / k_bar = 1.0 / 1.0 | 1.0 / 1.2 | 1.2 / 1.2 |', '|---|---|---|---|---|---|']
    for i in (13, 16, 17, 0):
        rr = [r for r in margin if r['index'] == i]
        cell = lambda g, b: next(f"L={r['L']:.3f}，{r['dominant_hz']:.1f} Hz" for r in rr if r['k_good_gain'] == g and r['k_bar_gain'] == b)
        text.append(f"| {i} | {rr[0]['group']} | {rr[0]['keff_median']:.0f} | {cell(1., 1.)} | {cell(1., 1.2)} | {cell(1.2, 1.2)} |")
    text += ['', '主频 20.8 Hz 是探测激励本身。这说明当前冻结的控制参数依赖高于 k_eff 的增益裕度，精确的 k_eff 不等于最优的控制用增益。', '',
        '## 表 2：刚度辨识精度与可靠性', '',
        '相同观测（Huber-D⁴ 闭环记录），真值 k_eff。误差为 ln(k̂/k_eff)。原始：窗口输出，无效时保持上一有效值；门控后：该方法门控后的 k_good，含 CENTRE 1.054。组内波动为每例稳健标准差的中位数。接受率、错误更新（已接受拍中 |误差|>ln1.2）和覆盖率（|误差|≤2σ）按原始窗口计。门控后行的“可用率”列指 k_good 已存在的比例；接受率为 0 的方法，k_good 停在首次获取值。“突发扰动段”含事件中刚度也在变的拍，分开的结果在 `identification_replay.json`。', '',
        '**适配边界。** 表中 Rossi-ID、Kalman、MHE、Xue 都是原文机制在本任务中的适配，不是原论文的完整复现，见 `baseline_protocol.md`。Xue 的回归量输出切线刚度，不是探测频率增益；其 σ 为固定接口值，覆盖率不能解释为不确定度质量。Kalman 与 MHE 的冻结配置过程噪声相同，同一观测下原始估计完全相同（线性高斯下 MHE 退化为 Kalman），两者差别只在门控和启动，不能算两份独立证据。']
    text += table2(seg)
    text += table2_scene(seg)
    text += ['', f'换用 Rossi-ID 闭环作记录源时，10 种方法按原始 |误差| 中位数排序的秩相关：{"，".join(sens)}。', '',
        '## 表 3：刚度变化跟踪与恢复能力', '',
        '与表 2 同一次回放。大过渡 |Δln k|≥0.3，全部场景，不与事件重叠。滞后在扣除过渡前偏差后计；最大低估和恢复用原始误差，最大低估为负表示全程高于真值。恢复：|误差|≤ln1.2 并保持 0.2 s，上限 2 s，否则计未恢复。未恢复、失联事件都保留，失联相关指事件发生在首次失联之后或区间内失联。“其中数据截断”指回合结束时观察不足 2 s。', '',
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
    ap.add_argument('--stage', choices=('closed', 'replay', 'score', 'oracle_margin', 'tables'), required=True)
    args = ap.parse_args()
    O.setup()
    globals()[args.stage]()
