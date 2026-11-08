"""Does the arm approximation change identification or the loop? Same 36 cases of corrected_comparison, old equation against the dynamic arm.

The option off must reproduce results/corrected_comparison.json exactly.
Usage: python -X utf8 plant_arm_effect.py
"""
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import overlay_compare as O
import numpy as np

import corrected_comparison as CC
import memory_segment_eval as SE
import method3 as M
import physics_plant3 as P3
import rib_transition_eval as R

HERE = Path(__file__).resolve().parent
ARMS = (('old', {}), ('z.03', dict(arm_zeta=.03)), ('z.1', dict(arm_zeta=.1)), ('z.3', dict(arm_zeta=.3)),
        ('z.1 pole1', dict(arm_zeta=.1, actuator_t2=.001)), ('z.1 pole4', dict(arm_zeta=.1, actuator_t2=.004)))
METHODS = ('Huber-D4', 'off_none', 'kal_short')


def job(a):
    method, case, arm, extra = a
    O.setup()
    sc = SE.MC.scene_for(tuple(case['cell']), CC.BASE, case['variant'])
    sc.update(extra)
    pol = O.B.build('Huber-D4', O.chosen()['Huber-D4'], sc)
    if method != 'Huber-D4':
        SE.gates_off(pol)
    if method == 'kal_short':
        R.make_filter(pol, method)
    raw = []
    est = pol.estimator; orig = est.estimate
    def estimate(o, n, policy):
        out = orig(o, n, policy)
        raw.append((n-1, out[0], pol.k_good))
        return out
    est.estimate = estimate
    log = P3.episode(sc, pol, case['seed'])
    m = M.metric(log, pol, sc)
    if log is None:
        return dict(method=method, arm=arm, group=case['group'], index=case['index'], L=float(np.log1p(m['J'])), rmse=float('nan'), peak=5000., raw_bias=float('nan'), raw_sd=float('nan'), good_bias=float('nan'), good_sd=float('nan'), keff_over_col21=float('nan'), diverged=True)
    r = np.array(raw).reshape(-1, 3); idx = r[:, 0].astype(int); idx = np.minimum(idx, len(log)-1)
    w = (log[idx, 0] >= 1.) & (log[idx, 26] > 200.) & np.isfinite(r[:, 2]) & (r[:, 1] > 0)
    e_raw = np.log(r[w, 1]/log[idx[w], 26]); e_good = np.log(r[w, 2]/log[idx[w], 26])
    iqr = lambda v: float((np.percentile(v, 75)-np.percentile(v, 25))/1.349) if len(v) else float('nan')
    med = lambda v: float(np.median(v)) if len(v) else float('nan')
    c = log[:, 21] > 200.
    return dict(method=method, arm=arm, group=case['group'], index=case['index'], L=float(np.log1p(m['J'])), rmse=float(m['rmse']), peak=float(m['peak']),
                raw_bias=med(e_raw), raw_sd=iqr(e_raw), good_bias=med(e_good), good_sd=iqr(e_good), keff_over_col21=med(log[c, 26]/log[c, 21]), diverged=False)


def main():
    start = time.perf_counter()
    cases = list(CC.cases())
    jobs = [(m, c, a, e) for a, e in ARMS for m in METHODS for c in cases if not (a.startswith('z.1 pole') and m != 'Huber-D4')]
    print('jobs', len(jobs), flush=True)
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(job, jobs, chunksize=2))
    (HERE/'results'/'plant_arm_effect.json').write_text(json.dumps(rows, indent=1))
    stored = json.load(open(HERE/'results'/'corrected_comparison.json'))['episodes']
    st = {(e['method'], e['index']): e for e in stored}
    d = [abs(r['L']-st[(r['method'], r['index'])]['L']) for r in rows if r['arm'] == 'old' and r['method'] in ('Huber-D4', 'off_none', 'kal_short')]
    print('diverged', [(r['method'], r['arm'], r['group'], r['index']) for r in rows if r['diverged']], flush=True)
    print('old arm vs stored: max |dL|', max(d), flush=True)
    for m in METHODS:
        print('--', m)
        for a, _ in ARMS:
            q = [r for r in rows if r['method'] == m and r['arm'] == a]
            if not q:
                continue
            line = ' '.join(f"{g}:{np.mean([r['L'] for r in q if r['group']==g]):.3f}/{max(r['peak'] for r in q if r['group']==g):.0f}" for g in CC.SCENES)
            print(f"{a:10s} {line} | k_raw bias {np.nanmedian([r['raw_bias'] for r in q]):+.3f} sd {np.nanmedian([r['raw_sd'] for r in q]):.3f} k_good bias {np.nanmedian([r['good_bias'] for r in q]):+.3f} sd {np.nanmedian([r['good_sd'] for r in q]):.3f} keff/col21 {np.nanmedian([r['keff_over_col21'] for r in q]):.3f} diverged {sum(r['diverged'] for r in q)}", flush=True)
    print('seconds', round(time.perf_counter()-start))


if __name__ == '__main__':
    main()
