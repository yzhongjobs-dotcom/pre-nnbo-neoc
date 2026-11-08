"""Simulation-only 2x2 Huber-D4 segment evaluation; baseline must pass first.

Run with python -X utf8; writes only results/memory_segment_eval.json.
No truth is exposed to the policy. Scoring uses the returned plant log only.
"""
import argparse
import inspect
import json
import math
import textwrap
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import MethodType

# Import before numpy: baseline_compare sets BLAS thread limits to one.
import overlay_compare as O
import numpy as np
from scipy.ndimage import uniform_filter1d

import memory_closed_loop as MC
import memory_filters as MF
import method3 as M
import physics_plant3 as P3
import estimator_swap as ES
import contact_compare as CC
import contact_candidates as C

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results/memory_segment_eval.json'
CONFIGS = ('on_none', 'on_memory', 'off_none', 'off_memory')
SEGMENTS = ('transition', 'plateau', 'event')


def bind_source(pol, fn, replacements):
    """Private method copies only; fail if the host source no longer matches."""
    source = textwrap.dedent(inspect.getsource(fn))
    for old, new in replacements:
        assert source.count(old) == 1, old
        source = source.replace(old, new)
    namespace = {}
    exec(compile(source, '<instance segment ablation>', 'exec'),
         fn.__globals__, namespace)
    return MethodType(namespace[fn.__name__], pol)


def gates_off(pol):
    # Keep the original finite/info check, clipping, diagnostics, and the entire
    # first-acquisition branch. Replace only the k_good acceptance block and the
    # post-acquisition k_bar acceptance branches. ModelPolicy's acquisition and
    # CENTRE wrapper and CandidatePolicy's hooks remain identical private copies.
    original = textwrap.dedent(inspect.getsource(ES.SwapPolicy._stiffness_est))
    start = original.index("    if p['kgood_rule'] == 'filter':")
    end = original.index("    if not p['conservative']:", start)
    acceptance = original[start:end]
    direct = """    self.k_good = k_hat
    self.f_set = max(float(o.log[n-1, 5]), .3)
    self.kg_rel = sig/k_hat
"""
    tail_start = original.index('    if precise and consistent:\n        return max(cand,')
    tail = original[tail_start:]
    pol._segment_swap = bind_source(pol, ES.SwapPolicy._stiffness_est, [
        (acceptance, direct),
        (tail, "    return max(cand, (1-p['rho_d'])*self.kbar), True\n")])
    pol._segment_model = bind_source(pol, CC.ModelPolicy._stiffness_est, [
        ('super()._stiffness_est(o, n)', 'self._segment_swap(o, n)')])
    pol._stiffness_est = bind_source(pol, C.CandidatePolicy._stiffness_est, [
        ('super()._stiffness_est(o, n)', 'self._segment_model(o, n)')])
    assert pol.params['ff_band'] is None
    pol.params['ff_jump'] = None


def expand(mask, t, seconds):
    hits = t[mask]
    if not len(hits):
        return np.zeros(len(t), bool)
    j = np.searchsorted(hits, t)
    left = np.maximum(j-1, 0)
    right = np.minimum(j, len(hits)-1)
    return np.minimum(abs(t-hits[left]), abs(t-hits[right])) <= seconds+1e-10


def segment_masks(log, sc):
    t, k = log[:, 0], log[:, 21]
    contact = (k > 0.) & (log[:, 3] > 0.) & np.isfinite(k)
    r = np.full(len(t), np.nan)
    # Centered 0.2 s arithmetic boxcar of k, nearest endpoint extension within
    # each continuous contact run; no smoothing across contact-loss boundaries.
    starts = np.flatnonzero(contact & ~np.r_[False, contact[:-1]])
    ends = np.flatnonzero(contact & ~np.r_[contact[1:], False])+1
    dt = float(np.median(np.diff(t)))
    width = max(1, int(round(.2/dt)))
    for a, b in zip(starts, ends):
        if b-a >= 2:
            smooth = uniform_filter1d(k[a:b], size=width, mode='nearest')
            r[a:b] = abs(np.gradient(np.log(smooth), t[a:b]))
    transition = expand(contact & (r > .2), t, .1)
    plateau = contact & (r < .05) & ~expand(transition, t, .3)
    event = np.zeros(len(t), bool)
    if sc['cond'] == 'breath_fast':
        for _, a, b in sc.get('events', []):
            event |= (t >= a) & (t <= b+.3)
    eligible = (t >= 1.) & contact
    masks = dict(event=eligible & event,
                 transition=eligible & transition & ~event,
                 plateau=eligible & plateau & ~event)
    assert not np.any(masks['transition'] & masks['plateau'])
    return eligible, masks


def distribution(x, fractions=False):
    x = x[np.isfinite(x)]
    if not len(x):
        out = dict(n=0, median=None, robust_sd=None)
        if fractions:
            out.update(below_1=None, below_085=None)
        return out
    q25, med, q75 = np.percentile(x, [25, 50, 75])
    out = dict(n=len(x), median=float(med), robust_sd=float((q75-q25)/1.349))
    if fractions:
        out.update(below_1=float(np.mean(x < 0.)),
                   below_085=float(np.mean(x < math.log(.85))))
    return out


def job(task):
    config, cell, variant, q, off_hop = task
    O.setup()
    sc = MC.scene_for(cell, MC.DEV_BASE, variant)
    pol = O.B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    if config.startswith('off'):
        gates_off(pol)
    state = dict(raw=np.nan)
    estimate = pol.estimator.estimate

    def raw_capture(o, n, policy):
        out = estimate(o, n, policy)
        state['raw'] = out[0]
        return out

    pol.estimator.estimate = raw_capture
    if config == 'on_memory':
        MF.attach(pol, 'forgetting', lam_max=.998)
    elif config == 'off_memory':
        flt = MF.attach(pol, 'forgetting', lam_max=.998, hop=off_hop)
        update = flt.update
        flt.update = lambda ell, R, ds, usable: update(ell, R, ds, True)
    records = []
    step = pol.step

    def capture(o):
        state['raw'] = np.nan
        command = step(o)
        records.append((float(o.log[o.n-1, 0]), pol.k_good, pol.kbar, state['raw']))
        return command

    pol.step = capture
    log = P3.episode(sc, pol, MC.SEED+q)
    if log is None:
        raise RuntimeError(f'Plant failed: {task}')
    rec = np.asarray(records)
    assert len(rec) == len(log)
    assert np.allclose(rec[:, 0], log[:, 0], rtol=0., atol=1e-9)
    metric = M.metric(log, pol, sc)
    eligible, masks = segment_masks(log, sc)
    with np.errstate(divide='ignore', invalid='ignore'):
        errors = np.log(rec[:, 1:4]/log[:, 21, None])
    peak_j = int(np.argmax(np.where(log[:, 0] >= 1., log[:, 18], -np.inf)))
    peak_segment = next((s for s in SEGMENTS if masks[s][peak_j]), 'unscored')
    row = dict(config=config, scene=cell[0], variant=variant,
               dynamics=cell[2], phase=cell[3], i=cell[4], seed=MC.SEED+q,
               L=float(np.log1p(metric['J'])), peak=float(metric['peak']),
               peak_segment=peak_segment, eligible_ticks=int(eligible.sum()), segments={})
    for s in SEGMENTS:
        mask = masks[s]
        force_error = log[mask, 18]-8.
        row['segments'][s] = dict(
            ticks=int(mask.sum()), share=float(mask.sum()/eligible.sum()),
            e_g=distribution(errors[mask, 0], fractions=True), e_bar=distribution(errors[mask, 1]),
            e_raw=distribution(errors[mask, 2]),
            force_rmse=float(np.sqrt(np.mean(force_error**2))) if len(force_error) else None,
            force_p99=float(np.percentile(abs(force_error), 99)) if len(force_error) else None,
            peak_over_15=int(row['peak'] > 15. and peak_segment == s),
            _force_abs=abs(force_error).tolist())
    return row


def mean_present(values):
    values = [v for v in values if v is not None]
    return float(np.mean(values)) if values else None


def aggregate(rows):
    groups = []
    for scene, variant in (('breath_rib', 'nominal'), ('breath_rib', 'slow_thin'),
                           ('breath_fast', 'nominal')):
        for config in CONFIGS:
            rr = [r for r in rows if (r['scene'], r['variant'], r['config']) ==
                  (scene, variant, config)]
            if not rr:
                continue
            group = dict(scene=scene, variant=variant, config=config, episodes=len(rr),
                         L_mean=float(np.mean([r['L'] for r in rr])),
                         L_median=float(np.median([r['L'] for r in rr])),
                         max_peak=max(r['peak'] for r in rr), segments={})
            for s in SEGMENTS:
                parts = [r['segments'][s] for r in rr]
                force = np.asarray([v for p in parts for v in p['_force_abs']])
                summary = dict(episodes_with_segment=sum(p['ticks'] > 0 for p in parts),
                               ticks=sum(p['ticks'] for p in parts),
                               share=float(np.mean([p['share'] for p in parts])),
                               force_rmse=float(np.sqrt(np.mean(force**2))) if len(force) else None,
                               force_p99=float(np.percentile(force, 99)) if len(force) else None,
                               peak_over_15=sum(p['peak_over_15'] for p in parts))
                for channel in ('e_g', 'e_bar', 'e_raw'):
                    fields = ('median', 'robust_sd')
                    if channel == 'e_g':
                        fields += ('below_1', 'below_085')
                    summary[channel] = {key: mean_present([p[channel][key] for p in parts])
                                        for key in fields}
                group['segments'][s] = summary
            group['peak_over_15_unscored'] = sum(
                r['peak'] > 15. and r['peak_segment'] == 'unscored' for r in rr)
            groups.append(group)
    return groups


def check_baseline(rows):
    stored = json.loads((HERE/'results/memory_closed_loop_true_base.json').read_text(encoding='utf-8'))
    key = lambda r: (r['scene'], r['variant'], r['dynamics'], r['phase'], r['i'])
    refs = {key(r): r for r in stored if r['name'] == 'Huber-true dev' and r['mass']}
    checks = []
    for r in rows:
        ref = refs[key(r)]
        checks.append(dict(key=list(key(r)), L_delta=r['L']-ref['L'],
                           peak_delta=r['peak']-ref['peak']))
    passed = all(abs(c['L_delta']) <= 1e-10 and abs(c['peak_delta']) <= 1e-8 for c in checks)
    return dict(passed=passed, episodes=len(checks), L_atol=1e-10, peak_atol=1e-8,
                max_abs_L_delta=max(abs(c['L_delta']) for c in checks),
                max_abs_peak_delta=max(abs(c['peak_delta']) for c in checks), checks=checks)


def run_pool(ex, tasks):
    rows = []
    futures = [ex.submit(job, task) for task in tasks]
    for f in as_completed(futures):
        rows.append(f.result())
        if len(rows) % 12 == 0 or len(rows) == len(tasks):
            print(f'completed {len(rows)}/{len(tasks)}', flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    parser.add_argument('--off-memory-hop', type=int, choices=(1, 6), default=1)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--resume', action='store_true', help='Reuse a passed baseline from this script only.')
    args = parser.parse_args()
    start = time.perf_counter()
    # slow_thin uses the filtered-grid index in dev3. A is the first 36 cells,
    # so that index is exactly the original O.cells index in this selection.
    episodes = [(c, 'nominal', q) for q, c in enumerate(O.cells())
                if c[0] in ('breath_rib', 'breath_fast')]
    episodes += [(c, 'slow_thin', q) for q, c in enumerate(O.cells()) if c[0] == 'breath_rib']
    assert len(episodes) == 84
    metadata = dict(base=MC.DEV_BASE, seed_base=MC.SEED, probe_mass=True,
                    host='unchanged Huber-D4', workers=args.workers,
                    centering_factor=float(CC.CENTRE), window_ticks=24,
                    on_memory="MF.attach(pol, 'forgetting', lam_max=.998), unchanged",
                    off_memory_hop=args.off_memory_hop,
                    off_memory_overlap=f'nw/{args.off_memory_hop}',
                    off_overrides=[
                        'Instance-only private copies of SwapPolicy, ModelPolicy, CandidatePolicy stiffness methods; no class/global mutation.',
                        'Finite estimate and fit.info>0 required; clipping 200..300000 retained.',
                        'Every valid estimate sets k_good=k_hat, f_set=max(F_measured,.3), kg_rel=sigma/k_hat; existing CENTRE wrapper retained.',
                        'Precision, chi and measurement-bound tests do not block k_good or post-acquisition k_bar.',
                        'Original conservative first-acquisition branch and ModelPolicy acquisition wrapper retained.',
                        'Post-acquisition k_bar=max(clip(k_hat+c*sigma,200,300000),(1-rho_d)*previous_k_bar).',
                        'ff_jump=None; ff_band already None. Contact supervisor and servo identification unchanged.',
                        'Off-memory: every filter update has usable=True; finite/positive input guard retained; filter sigma used.',
                    ],
                    scoring=dict(truth='log column 21 |K(20.8)|, offline scoring only',
                                 eligible='t>=1, true normal force>0 and k>0',
                                 smoothing='0.2 s centered arithmetic boxcar of k; nearest ends within each continuous contact run',
                                 transition='abs(gradient(log(smoothed k)))>.2/s, expanded by .1 s each side',
                                 plateau='r<.05/s, farther than .3 s from expanded transition',
                                 event='scene B sc.events through end+.3 s; priority over both other segments',
                                 share='mean of episode segment ticks / episode eligible contact ticks',
                                 distributions='mean of per-episode median, IQR/1.349 and fractions; empty episodes omitted per channel',
                                 raw='current window output before memory, before CENTRE, only ticks actually estimated; missing raw never forward-filled',
                                 force='RMSE and P99 pooled segment ticks across episodes, true axial force log[:,18]-8',
                                 peak='whole-episode peak from t>=1; segment at peak tick; unscored peaks retained separately'))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        if args.resume:
            saved = json.loads(OUT.read_text(encoding='utf-8'))
            assert saved['baseline_check']['passed']
            rows = [r for r in saved['episodes'] if r['config'] == 'on_none']
            assert len(rows) == 84 and all('_force_abs' in r['segments']['transition'] for r in rows)
        else:
            rows = run_pool(ex, [('on_none', c, v, q, args.off_memory_hop) for c, v, q in episodes])
        check = check_baseline(rows)
        print('baseline_check', json.dumps({k: v for k, v in check.items() if k != 'checks'}), flush=True)
        saved = dict(metadata=metadata, baseline_check=check, episodes=rows)
        OUT.write_text(json.dumps(saved, ensure_ascii=False, indent=1, allow_nan=False), encoding='utf-8')
        if not check['passed']:
            raise SystemExit('STOP: Gates ON/no-memory did not reproduce stored L and peak.')
        if args.baseline_only:
            return
        rows += run_pool(ex, [(config, c, v, q, args.off_memory_hop)
                              for config in CONFIGS[1:] for c, v, q in episodes])
    saved['summary'] = aggregate(rows)
    for row in rows:
        for part in row['segments'].values():
            del part['_force_abs']
    rows.sort(key=lambda r: (r['variant'], r['scene'], CONFIGS.index(r['config']), r['seed']))
    saved['seconds'] = time.perf_counter()-start
    OUT.write_text(json.dumps(saved, ensure_ascii=False, indent=1, allow_nan=False), encoding='utf-8')
    print('saved', OUT, 'seconds', round(saved['seconds']), flush=True)


if __name__ == '__main__':
    main()
