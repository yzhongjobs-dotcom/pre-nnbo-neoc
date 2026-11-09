"""Port of rib_contrast_check.py / rib_contrast_check2.py; no parameter tuning.

Run after identification_completion.py --stage score, with python -X utf8.
All ratios, including individual crossings, are saved for inspection.
"""
import json
import math
import argparse
import hashlib
from pathlib import Path

import identification_completion as IC  # thread limits before numpy
import numpy as np
from scipy.signal import find_peaks

CC, SE, O, P, R = IC.CC, IC.SE, IC.O, IC.P3, IC.R
HERE = IC.HERE
FIG = dict(
    measured_p1=[6.85, 13.2, 6.3, 10.35, 5.25, 7.65, 5.5],
    fea_p1=[3.75, 6.4, 3.4, 7.15, 2.85, 6.8, 2.5],
    fea_p2=[3.75, 4.55, 3.5, 4.1, 3.2, 3.8, 3.3],
    fea_p3=[9.8, 11.1, 10.0, 10.95, 9.4, 10.75, 9.7],
)
GROUPS = ('A', 'A_slow_thin', 'A+B', 'A_large')


def crossings(y):
    """Consecutive local extrema; same log(1.03) prominence as source script."""
    ly = np.log(y)
    pk, _ = find_peaks(ly, prominence=math.log(1.03))
    tr, _ = find_peaks(-ly, prominence=math.log(1.03))
    idx = np.sort(np.r_[pk, tr])
    return np.exp(abs(np.diff(ly[idx]))) if len(idx) >= 2 else np.array([])


def distribution(v):
    a = np.asarray(v, float)
    return dict(n=len(a), median=float(np.median(a)) if len(a) else None,
        p10=float(np.quantile(a, .1)) if len(a) else None,
        p90=float(np.quantile(a, .9)) if len(a) else None,
        min=float(a.min()) if len(a) else None, max=float(a.max()) if len(a) else None)


def figure():
    result = {}
    for name, values in FIG.items():
        d = np.asarray(values)
        adj = np.maximum(d[1:], d[:-1])/np.minimum(d[1:], d[:-1])
        result[name] = dict(depth_mm=values, adjacent_ratios=adj.tolist(),
            adjacent=distribution(adj), overall_max_min=float(d.max()/d.min()))
    return result


def static():
    # Same six neural_white rib cells as rib_contrast_check2.py, same base.
    cases = [c for c in CC.cases() if c['group'] in ('A', 'A_slow_thin') and c['cell'][2] == 'neural_white']
    rows, soft = [], []
    for slide in ('long', 'short'):
        for cover in (.007, .018, .032):
            per_case = []
            for case in cases:
                c, t, dy, ph, i = case['cell']
                sc = O.B.S3.cell_scene3(c, t, 'nominal', i, base=CC.BASE,
                    overrides=dict(slide=slide, force=8., dynamics=dy, phase=ph,
                                   cover=cover, field_sigma=0., surface_amp=0.))
                ribs = sc['ribs']
                xs = np.arange(float(ribs[0, 0]-.02), float(ribs[min(3, len(ribs)-1), 0]+.02), .0005)
                d20 = np.array([P.solve_depth(sc, 20., x) for x in xs])
                d8 = np.array([P.solve_depth(sc, 8., x) for x in xs])
                k8 = np.array([P.contact(sc, d, x)[1] for d, x in zip(d8, xs)])
                gain = np.array([abs(P.complex_gain(sc, k, 8.)) for k in k8])
                supported = 1./(1./gain+1./sc['K_t'])
                per_case.append(dict(index=case['index'], cell=case['cell'],
                    depth20_ratios=crossings(d20).tolist(), tangent8_ratios=crossings(k8).tolist(),
                    supported_gain8_ratios=crossings(supported).tolist()))
                # Remove ribs from exactly this static scene; do not change orientation or moduli.
                no_ribs = dict(sc, ribs=np.empty((0, 4)))
                x = float(sc['xcmd'][0])
                depth = P.solve_depth(no_ribs, 8., x)
                tangent = P.contact(no_ribs, depth, x)[1]
                g = float(abs(P.complex_gain(no_ribs, tangent, 8.)))
                soft.append(dict(index=case['index'], slide=slide, cover_mm=cover*1000,
                    depth8_mm=depth*1000, tangent8=tangent, abs_G8=g,
                    supported_gain8=float(1./(1./g+1./sc['K_t']))))
            summary = {}
            for key in ('depth20_ratios', 'tangent8_ratios', 'supported_gain8_ratios'):
                summary[key] = distribution([v for r in per_case for v in r[key]])
            rows.append(dict(slide=slide, cover_mm=cover*1000, cases=per_case, summary=summary))
    ratios = []
    for short in [s for s in soft if s['slide'] == 'short']:
        long = next(s for s in soft if s['slide'] == 'long' and (s['index'], s['cover_mm']) == (short['index'], short['cover_mm']))
        ratios.append(dict(index=short['index'], cover_mm=short['cover_mm'],
            abs_G8_short_over_long=short['abs_G8']/long['abs_G8'],
            supported_gain8_short_over_long=short['supported_gain8']/long['supported_gain8']))
    return rows, dict(values=soft, ratios=ratios,
        abs_G8_ratio=distribution([r['abs_G8_short_over_long'] for r in ratios]),
        supported_gain8_ratio=distribution([r['supported_gain8_short_over_long'] for r in ratios]))


def logged():
    rows = []
    for case in CC.cases():
        if case['group'] not in GROUPS:
            continue
        _, log, _ = IC.load_case(IC.TICKS/f"Huber-D4_{case['index']:02d}.npz")
        # The transition detector already implements the exact .1 s boxcar,
        # separately inside contiguous contact runs (no stitching across loss).
        smooth, contact, _, _ = R.transitions(log)
        valid = contact & (log[:, 0] >= 1.) & np.isfinite(smooth)
        ratios = [v for a, b in R.spans(valid) for v in crossings(np.exp(smooth[a:b]))]
        rows.append(dict(index=case['index'], group=case['group'], ratios=ratios, summary=distribution(ratios)))
    return dict(cases=rows, groups={g: distribution([v for r in rows if r['group'] == g for v in r['ratios']]) for g in GROUPS})


def transitions():
    data = json.loads((HERE/'results/identification_replay.json').read_text(encoding='utf-8'))
    rows = []
    for group in CC.SCENES:
        events = [dict(index=c['index'], **{k: e[k] for k in ('t0', 't1', 'delta', 'direction', 'size_class', 'after_contact_loss', 'contact_loss_in_interval')})
            for c in data['cases'] if c['source'] == 'Huber-D4' and c['group'] == group
            for e in c['methods']['Huber-D4']['transitions'] if not e['overlaps_event']]
        size = [abs(e['delta']) for e in events]
        rows.append(dict(group=group, events=events, abs_delta_log=distribution(size),
            ratios=distribution(np.exp(size)), large=sum(s >= .3 for s in size),
            small=sum(s < .3 for s in size),
            duration=distribution([e['t1']-e['t0'] for e in events])))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stage = parser.add_mutually_exclusive_group()
    stage.add_argument('--static-only', action='store_true', help='Save the requested static and soft-region calculations while replay runs')
    stage.add_argument('--reuse-static', action='store_true', help='Finish from the saved static calculations without repeating them')
    stage.add_argument('--scored-only', action='store_true', help='Update only scored transitions from identification_replay.json')
    args = parser.parse_args()
    path = HERE/'results/rib_contrast_evidence.json'
    if args.scored_only:
        result = json.loads(path.read_text(encoding='utf-8'))
        result['scored_no_event'] = transitions()
        path.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2), encoding='utf-8')
        print('updated scored_no_event only', path, flush=True)
        return
    O.setup()
    sources = {Path(p).name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        for p in (__file__, CC.__file__, SE.MC.__file__, P.__file__, O.B.S3.__file__)}
    if args.reuse_static:
        prior = json.loads(path.read_text(encoding='utf-8'))
        assert prior['metadata']['static_source_sha256'] == sources, 'Static source changed'
        stat, soft = prior['static_ribs_only'], prior['soft_region']
    else:
        stat, soft = static()
    result = dict(metadata=dict(figure_source='read by eye from Al-Zogbi 2021 Fig. 4C',
        source_scripts=['rib_contrast_check.py', 'rib_contrast_check2.py'],
        static='Same six neural_white cells from A and A_slow_thin; field_sigma=0, surface_amp=0; .5 mm sampling over first four ribs; no body motion in contact-depth solve.',
        logged='Huber-D4 source; log[:,26], t>=1, contiguous true-contact runs, exact .1 s centered log boxcar, log(1.03) extrema prominence. Source scratchpad used a 13-tick approximation; here the existing exact detector boxcar is reused.',
        crossings='Ratios between consecutive interior local extrema; endpoints not counted.',
        soft='Same static scenes, ribs removed. |G| is abs(P.complex_gain) at 8 N; body-support series gain also reported. Ratio short/long; no adjustment.',
        static_source_sha256=sources, complete=not args.static_only),
        figure4C=figure(), static_ribs_only=stat, soft_region=soft)
    if not args.static_only:
        result.update(logged_keff=logged(), scored_no_event=transitions())
    path.write_text(json.dumps(CC.clean(result), ensure_ascii=False, allow_nan=False, indent=2), encoding='utf-8')
    print('saved', path, flush=True)


if __name__ == '__main__':
    main()
