"""P1b performance attribution of the final method (diagnosis cohort 110000+, 20 scenes per family, all families).
Component replacement (nominal platform): Proposed, TrueMeas, TrueFb, TrueK, TrueK+TrueMotion (excitation on),
Oracle (TrueK+TrueMotion, excitation off), NoRecovery.
Conditions for the final method only (controller keeps the nominal model): servo lag x0.6, 20 Hz force filter,
extra force noise (x3 variance, proxy for a coarser cobot F/T sensor). A 10 ms command period (RM65) needs a re-parametrised controller and is left to a separate check.
Metrics: RMSE median/mean, runaways (>50 N or diverged), long contact loss (>0.5 s), mean loss fraction.
Run: python -X utf8 p1_attribution.py -> results/p1_attribution.json
"""
import json, math, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
MA = os.path.join(HERE, '..', 'motion_aware_scan'); sys.path.insert(0, MA); sys.path.insert(0, HERE)
import common, variants, cohort_eval as CE, robustness as RB
import proposed as Pm
from step3_gap_decomposition import build as build_gap
from common import P

F = dict(CE.FINAL)
def build(cfg, p):
    if cfg == 'Proposed': return variants.build('Proposed-IMM', F)
    if cfg == 'TrueMeas': return variants.build('Proposed-IMM', dict(F, oracle_meas=True))
    if cfg == 'TrueFb': return variants.build('Proposed-IMM', dict(F, oracle_fb=True))
    if cfg == 'TrueK': return variants.build('Proposed-IMM', dict(F, oracle_meas=True, oracle_fb=True))
    if cfg == 'NoRecovery': return variants.build('Proposed-IMM', dict(F, recovery=False))
    if cfg == 'TrueK+TrueMotion': return build_gap('Oracle+exc', p)
    if cfg == 'Oracle': return build_gap('Oracle', p)
    raise KeyError(cfg)

COMP = ('Proposed', 'TrueMeas', 'TrueFb', 'TrueK', 'TrueK+TrueMotion', 'Oracle', 'NoRecovery')
COND = {'lag x0.6': dict(tau_a=0.019*0.6), 'filter 20 Hz': dict(fc=20.), 'noise x3 (RM65 F/T proxy)': 'noise'}

def _ep(a):
    kind, cfg, (fam, sid, p, T) = a
    try:
        if kind == 'comp':
            log, _ = common.run(build(cfg, p), p, T, sid)
        else:
            spec = COND[cfg]; plat = dict(RB.NOMINAL); nom = None; extra = 0.
            if isinstance(spec, dict): plat.update(spec)
            elif spec == 'noise': extra = 0.08*math.sqrt(2.)
            pol = RB.NominalModel(build('Proposed', p), extra_noise=extra, seed=sid)
            log = P.episode(p, pol, T=T, platform=plat, seed=sid)
        m = common.task_metrics(log, 20.8)
        return dict(kind=kind, cfg=cfg, family=fam, sid=sid, ok=not m['diverged'], rmse=m['rmse'], peak=m['peak'],
                    loss=m['loss'], loss_max=m['loss_max'])
    except Exception as ex:
        return dict(kind=kind, cfg=cfg, family=fam, sid=sid, ok=False, err=repr(ex))

if __name__ == '__main__':
    sc = CE.cohort('diag', 20)
    jobs = [('comp', c, s) for c in COMP for s in sc]+[('cond', c, s) for c in COND for s in sc]
    with Pool(16) as po: rows = po.map(_ep, jobs, chunksize=1)
    json.dump(rows, open(os.path.join(HERE, 'results', 'p1_attribution.json'), 'w'), indent=1)
    for kind, names in (('comp', COMP), ('cond', tuple(COND))):
        for fam in CE.FAMS:
            line = '%-10s' % fam
            for c in names:
                R = [r for r in rows if r['kind'] == kind and r['cfg'] == c and r['family'] == fam]
                ok = [r['rmse'] for r in R if r['ok']]
                bad = sum((not r['ok']) or r.get('peak', 0) > 50 for r in R)
                ll = sum(r.get('loss_max', 0) > .5 for r in R)
                line += ' | %s %.3f/%.3f r%d L%d' % (c, np.median(ok) if ok else np.nan, np.mean(ok) if ok else np.nan, bad, ll)
            print(line, flush=True)
