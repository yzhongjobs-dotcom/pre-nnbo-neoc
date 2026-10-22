"""Reproducible, bounded development and paired evaluation of contact updates.

Develop on seed 1000, validate the selected instances on seed 2000, then freeze
configuration/source hashes before the existing seed-7200 task cohort. No test
retuning, network training, or replacement of historical results. Run from any
directory. Commands: check, evaluate, summarize. Parameter selection is complete.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0,str(ROOT))
import plant as P
import scenarios as S
import report_physical_metrics as R
import complete_task_evidence as E
from contact_update import ContactUpdate, profile_contact

OUT = ROOT/'contact_update_results'


def save(name,value):
    OUT.mkdir(exist_ok=True)
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def read(name):
    return json.loads((OUT/name).read_text(encoding='utf-8'))


def cohort(task,seed,n):
    if seed==7200:
        return E.scenes(task)[:n]
    if task=='T1_static':
        return S.cohort('T1',n,seed,dict(ratio=30.,A=0.))
    if task=='T3_combined':
        return S.cohort('T3',n,seed,R.CELLS['T3'])
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        stiff = i%4>=2
        modulus = rng.uniform(5e5,1.5e6) if stiff else rng.uniform(1e4,4e4)
        p = S._base(rng,modulus,0.,.015,8.)
        p[15] = 1
        p[16:20] = [2.,.005*(1 if i%2==0 else -1),.5,2.]
        out.append((p,8.))
    return out


def original_config(task):
    source = 'T1' if task=='T1_static' else 'T3'
    return next(s[2] for s in R.specifications(source) if s[0]=='Candidate1_v2')


def build(method,cfg,task):
    if method=='C1_updated':
        return ContactUpdate(**dict(original_config(task),**cfg))
    if method=='C1_old_model_new_recovery':
        return ContactUpdate(**dict(original_config(task),**dict(cfg,projection=False)))
    if method=='C1_updated_gains_only':
        return ContactUpdate(**dict(original_config(task),**dict(cfg,reference=False)))
    raise ValueError(method)


def run(method,cfg,task,seed,n,keep=False):
    rows = []
    for i,(p,T) in enumerate(cohort(task,seed,n)):
        pol = build(method,cfg,task)
        if hasattr(pol,'bind'):
            pol.bind(p)
        log,tt,commands = E.episode(p,pol,T,seed+i,.008,timed=keep)
        row = dict(task=task,scene=i,method=method,**E.metrics(log,T,task))
        row['selection_J'] = S.metrics(log)['J']
        if tt.size:
            # JIT excluded by warming the full policy separately before test.
            row['update_p99_ms'] = float(np.percentile(tt[10:],99))
        rows.append(row)
        if keep:
            path = OUT/'traces'
            path.mkdir(exist_ok=True)
            np.savez_compressed(path/f'{task}_{method}_{i}.npz',
                                log=log if log is not None else np.empty((0,P.LOGC)),
                                commands=commands,p=p,
                                update_ms=tt,
                                details=np.array(getattr(pol,'details',[])))
    return rows


def score(rows):
    if any(r['failed'] for r in rows):
        return 1e6
    # Equal task/cell weight: the stiff-motion cell must not be diluted by
    # soft scenes. Original scalar J, with its declared force/loss penalties.
    groups = {}
    for r in rows:
        key = r['task']+('_hard' if r['scene']%4>=2 else '_soft') if r['task']=='T2_motion' else r['task']
        j = r['selection_J']
        groups.setdefault(key,[]).append(j)
    return float(np.mean([np.log(.01+np.mean(v)) for v in groups.values()]))


def compact(rows):
    out = {}
    for task in E.TASKS:
        for suffix in (['_soft','_hard'] if task=='T2_motion' else ['']):
            subset = [r for r in rows if r['task']==task and
                      (not suffix or (r['scene']%4>=2)==(suffix=='_hard'))]
            if subset:
                ok = [r for r in subset if not r['failed']]
                out[task+suffix] = dict(n=len(subset),failed=len(subset)-len(ok),
                    rmse_N=float(np.mean([r['rmse_N'] for r in ok])) if ok else None,
                    mean_peak_N=float(np.mean([r['peak_force_N'] for r in ok])) if ok else None,
                    worst_peak_N=float(max(r['peak_force_N'] for r in ok)) if ok else None,
                    mean_low_contact_s=float(np.mean([r['low_contact_s'] for r in ok])) if ok else None)
    return out




def check():
    # Known polynomial drift + independently excited encoder trajectory, using
    # exactly the servo interpolation/filter model of the estimator.
    checks = []
    for stiffness in (1500.,40000.):
        for excitation in (False,True):
            n = 100
            Tc = .008
            log = np.zeros((n,P.LOGC))
            log[:,0] = np.arange(n)*Tc
            log[:,11] = 5e-5*np.sin(np.arange(n)*2*np.pi/6) if excitation else 0.
            afc = 1-math_exp(-2*np.pi*30*.001)
            eb = math_exp(-.001/.019)
            y,F = 0.,4.
            for j in range(n):
                log[j,2],log[j,5] = y,F
                target = log[j-1,11] if j>=1 else 0.
                for mm in range(1,9):
                    y = target+(y-target)*eb
                    t = j*Tc+mm*.001
                    force = 4.+stiffness*y+.1*t+.03*t*t+.01*t**3
                    F += afc*(force-F)
            fit = profile_contact(log,n,24,.95,1,eb,afc,8,30000.,1e7,.04,3)
            if excitation:
                assert abs(fit[0]/stiffness-1)<1e-4,(stiffness,fit)
            else:
                assert abs(fit[0]-30000.)<1e-6,fit
                assert fit[4]<1e-12,fit
            checks.append(dict(stiffness=stiffness,excitation=excitation,
                               estimate=float(fit[0]),data_information=float(fit[4])))
    # Hidden truth columns must not affect a policy output.
    cfg = read('selection.json')['C1_updated'] if (OUT/'selection.json').exists() else {}
    p,T = cohort('T3_combined',2000,1)[0]
    pol = build('C1_updated',cfg,'T3_combined')
    log = P.episode(p,pol,T=.6,seed=2000,update=.008)
    a,b = build('C1_updated',cfg,'T3_combined'),build('C1_updated',cfg,'T3_combined')
    a.reset(P.PLATFORMS['PC125']);b.reset(P.PLATFORMS['PC125'])
    masked = log.copy();masked[:,[1,3,4,6,7,8,9,10]] = np.nan
    maxdiff = 0.
    for n in range(len(log)):
        g1,g2 = np.array(a(P.Obs(log,n))),np.array(b(P.Obs(masked,n)))
        assert np.isfinite(g2).all()
        maxdiff = max(maxdiff,float(np.max(abs(g1-g2))))
    assert maxdiff==0.
    save('checks.json',dict(polynomial_identification=checks,hidden_truth_max_difference=maxdiff))
    print('CHECKS PASS',checks,flush=True)




def math_exp(x):
    return float(np.exp(x))




def hashes():
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in ['contact_update.py','evaluate_contact_update.py',
                         'plant.py','candidate1_model.py','scenarios.py']}


def evaluate():
    E.verify()
    assert (OUT/'checks.json').exists()
    selected = {'C1_updated':read('selection.json')['C1_updated']}
    frozen = dict(configurations=selected,source_sha256=hashes(),seed=7200,n=40,
                  selection='seed-1000 training; seed-2000 validation reported without retuning',
                  task_protocol='same completed_tasks seed-7200 cohort; no test retuning',
                  information={'C1_updated':'force, encoder, sent commands'},
                  scope='Original three tasks only; no further simple-controller tuning or scene expansion',
                  comparisons='Frozen representative literature adaptations from completed_tasks; no new baseline tuning')
    if (OUT/'frozen.json').exists():
        assert read('frozen.json')==frozen
    else:
        save('frozen.json',frozen)
    rows = []
    for method,cfg in selected.items():
        # Full policy/JIT warmup outside timed episodes.
        run(method,cfg,'T3_combined',1000,1)
        for task in E.TASKS:
            rr = run(method,cfg,task,7200,40,keep=True)
            rows += rr
            print('TEST',method,task,compact(rr),flush=True)
            save('test_rows.json',rows)
    for method in ('C1_old_model_new_recovery',):
        for task in ('T2_motion','T3_combined'):
            rr = run(method,selected['C1_updated'],task,7200,40,keep=True)
            rows += rr
            print('ABLATION',method,task,compact(rr),flush=True)
            save('test_rows.json',rows)
    summarize()


def summarize():
    rows = read('test_rows.json')
    summary = {name:compact([r for r in rows if r['method']==name])
               for name in sorted(set(r['method'] for r in rows))}
    save('summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['check','evaluate','summarize'])
    args = parser.parse_args()
    start = time.perf_counter()
    globals()[args.command]()
    print('elapsed_s',round(time.perf_counter()-start,2),flush=True)
