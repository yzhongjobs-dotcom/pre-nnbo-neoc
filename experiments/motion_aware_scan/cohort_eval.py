"""Evaluate controller configurations on large seeded cohorts (robustness / failure-rate study, paper draft §五).

Cohorts are drawn per family with their own seeds (not the 12-scene splits): 'diag' (110000+, used for diagnosing
the stiff-push failures and screening fixes) and 'final' (130000+, untouched until the final evaluation).
Configurations are named in CONFIGS; comparators use their confirmed tuned parameters.
Run: CONFIGS=final,F1 COHORT=diag N=50 python -X utf8 cohort_eval.py -> results_v2/cohort_<cohort>.json (merged)
"""
import json, math, os
import numpy as np
from multiprocessing import Pool
import common
import scenes_v2 as SV
import tune_v2 as TV
import variants

H = common.HERE
RES = os.path.join(H, 'results_v2')
BASE = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMM.json')))['best']['params']
FAMS = ('S2', 'T3m', 'T2m-soft', 'T2m-stiff', 'T2f-stiff')
SEEDS = {'diag': 110000, 'final': 130000}


def tuned(m):
    return json.load(open(os.path.join(H, 'tuned_v2', m+'.json')))['best']['params']


def extra_configs():
    path = os.path.join(H, 'tuned_v2', 'cohort_configs.json')
    return json.load(open(path)) if os.path.exists(path) else {}


FINAL = json.load(open(os.path.join(H, 'tuned_v2', 'Proposed-IMMb.json')))['best']['params']
ABL = {'Abl-noProjection': dict(degree=0), 'Abl-noConservative': dict(conservative=False),
       'Abl-noFeedforward': dict(feedforward=False), 'Abl-noExcitation': dict(excitation=False),
       'Abl-noBound': dict(meas_bound=None), 'Abl-noRecovery': dict(recovery=False)}


def build(name, p_scene=None):
    import proposed as Pm
    if name == 'Proposed':
        return variants.build('Proposed-IMM', dict(FINAL))
    if name in ABL:
        return variants.build('Proposed-IMM', dict(FINAL, **ABL[name]))
    if name == 'Abl-KF':
        kf = json.load(open(os.path.join(RES, 'mechanism_checks.json')))['kf_matched']['best']
        common_ = {k: v for k, v in FINAL.items() if k not in ('q_smooth', 'q_maneuver', 'p_stay', 'r_scale')}
        return variants.build('Proposed-KF', dict(common_, **kf))
    if name == 'Abl-polyFF':
        return Pm.MotionAware(**{k: v for k, v in FINAL.items() if k not in ('q_smooth', 'q_maneuver', 'p_stay', 'r_scale')})
    if name.startswith('NoNorm:'):
        import mechanism_checks as MC
        return MC.FixedFeedback(float(name.split(':')[1]), single_model=False, **FINAL)
    if name == 'Oracle':
        return Pm.Oracle(p_scene, **{k: FINAL[k] for k in ('M1', 'D1', 'KPF1', 'ki')})
    if name.startswith('Fixed:'):
        return Pm.FixedEnd(float(name.split(':')[1]), **{k: FINAL[k] for k in ('M1', 'D1', 'KPF1', 'ki')})
    if name == 'YuenType':
        from comparators import YuenType
        return YuenType(p_scene, **tuned(name))
    fixes = {'final': {}, 'F1': dict(kgood_rule='raise'), 'F2': dict(follow_tau=0.3),
             'F12': dict(kgood_rule='raise', follow_tau=0.3)}
    fixes.update(extra_configs())
    if name in fixes:
        # names starting with 'P:' modify the final method, others the earlier selected parameters
        return variants.build('Proposed-IMM', dict(FINAL if name.startswith('P:') else BASE, **fixes[name]))
    return TV.make(name, tuned(name))


def cohort(name, n):
    out = []
    for fi, fam in enumerate(FAMS):
        rng = np.random.default_rng(SEEDS[name]+1000*fi)
        for i in range(n):
            p, T = SV.make(fam, rng)
            out.append((fam, SEEDS[name]+1000*fi+i, p, T))
    return out


def _ep(a):
    cfg, (fam, sid, p, T) = a
    try:
        log, _ = common.run(build(cfg, p), p, T, sid)
        m = common.task_metrics(log, 20.8)
    except Exception as ex:
        return dict(cfg=cfg, family=fam, sid=sid, ok=False, err=str(ex))
    return dict(cfg=cfg, family=fam, sid=sid, ok=not m['diverged'], rmse=m['rmse'], peak=m['peak'],
                loss=m['loss'], loss_max=m['loss_max'], J=m['J'])


def summarize(rows, cfgs):
    for fam in FAMS:
        for c in cfgs:
            R = [r for r in rows if r['cfg'] == c and r['family'] == fam]
            ok = [r for r in R if r['ok']]
            if not R:
                continue
            pk = np.array([r['peak'] for r in ok]); rm = np.array([r['rmse'] for r in ok]); lm = np.array([r['loss_max'] for r in ok])
            print('%-9s %-12s n=%d fail %d | RMSE median %.3f mean %.3f | peak>15N %d >50N %d | longest loss>0.5s %d >1s %d | mean loss %.4f'
                  % (fam, c, len(R), len(R)-len(ok), np.median(rm), rm.mean(), (pk > 15).sum(), (pk > 50).sum(),
                     (lm > .5).sum(), (lm > 1).sum(), np.mean([r['loss'] for r in ok])), flush=True)


if __name__ == '__main__':
    cfgs = [c for c in os.environ.get('CONFIGS', 'final').split(',') if c]
    name = os.environ.get('COHORT', 'diag'); n = int(os.environ.get('N', '50'))
    fams = os.environ.get('FAMS')
    sc = cohort(name, n)
    if fams:
        sc = [s for s in sc if s[0] in fams.split(',')]
    path = os.path.join(RES, 'cohort_%s.json' % name)
    old = json.load(open(path)) if os.path.exists(path) else []
    old = [r for r in old if r['cfg'] not in cfgs or (fams and r['family'] not in fams.split(','))]
    with Pool(16) as po:
        rows = po.map(_ep, [(c, s) for c in cfgs for s in sc], chunksize=1)
    json.dump(old+rows, open(path, 'w'), indent=1)
    summarize(old+rows, cfgs)
