"""Closed-loop paired evaluation of the V2 contact bank on the diagnostic cohort. Truth is scoring only."""
import argparse
import json
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
from joint_policy import JointPolicy
from conditional_filter import MeasuredOnly
from bank_policy import BankContact

CONFIGS = {
    'v1': None, 'lean': None,
    'gain_only': dict(feedback='projection'),
    'gain_fb': dict(feedback='filter'),
    'noprobe': dict(feedback='projection', probe_rescale=False),
    'slow.3': dict(feedback='projection', q_slow=.3),
    'slow.3_fb': dict(feedback='filter', q_slow=.3),
    'gain_only_q1': dict(feedback='projection', q_motion=.6),
    'gain_only_q5': dict(feedback='projection', q_motion=6.),
    'gain_only_slow.3': dict(feedback='projection', q_slow=.3),
    'gain_only_fast3': dict(feedback='projection', q_fast=3.),
}


def job(a):
    name, scene, lag = a; fam, sid, p, T = scene
    if name == 'v1':
        pol = ServoID(single_model=False, **FINAL)
    elif name == 'lean':
        pol = JointPolicy(wave=None, new_fit=True, rate=False, amplitude_scale=1., **FINAL)
    else:
        pol = BankContact(single_model=False, **CONFIGS[name], **FINAL)
    plat = dict(common.P.PLATFORMS['PC125']); plat['tau_a'] *= lag
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=plat)
    if log is None:
        return dict(config=name, family=fam, sid=sid, lag=lag, rmse=None, peak=None, loss=None)
    m = common.task_metrics(log, 20.8)
    calls = np.asarray(pol.t_call[126:])*1000
    return dict(config=name, family=fam, sid=sid, lag=lag, rmse=float(m['rmse']), peak=float(m['peak']), loss=float(m['loss_max']),
                reseeds=getattr(pol, 'reseeds', 0), pub=getattr(pol, 'published', 0)/max(getattr(pol, 'total', 1), 1),
                p50=float(np.percentile(calls, 50)), p99=float(np.percentile(calls, 99)))


def boot(d, B=4000):
    rng = np.random.default_rng(7); d = np.asarray(d)
    m = d[rng.integers(0, len(d), (B, len(d)))].mean(axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--configs', nargs='+', default=['v1', 'lean', 'gain_only', 'gain_fb'])
    ap.add_argument('--scenes', type=int, default=50)
    ap.add_argument('--lag', type=float, default=1.)
    ap.add_argument('--label', required=True)
    a = ap.parse_args()
    scenes = cohort('diag', a.scenes)
    with ProcessPoolExecutor(max_workers=15) as pool:
        rows = list(pool.map(job, [(c, s, a.lag) for c in a.configs for s in scenes], chunksize=3))
    out = {}
    for c in a.configs:
        for fam in FAMS:
            r = [x for x in rows if x['config'] == c and x['family'] == fam]
            ok = [x['rmse'] for x in r if x['rmse'] is not None]
            b = {x['sid']: x for x in rows if x['config'] == 'v1' and x['family'] == fam}
            d = [x['rmse']-b[x['sid']]['rmse'] for x in r if x['rmse'] is not None and b.get(x['sid'], {}).get('rmse') is not None]
            s = dict(med=round(float(np.median(ok)), 3), mean=round(float(np.mean(ok)), 3),
                     r=sum(1 for x in r if x['rmse'] is None or x['peak'] > 50), L=sum(1 for x in r if x['loss'] and x['loss'] > .5),
                     diff=round(float(np.mean(d)), 3) if c != 'v1' else 0., ci=[round(v, 3) for v in boot(d)] if c != 'v1' else None,
                     reseeds=round(float(np.mean([x.get('reseeds', 0) for x in r])), 1), pub=round(float(np.mean([x.get('pub', 0) for x in r])), 2),
                     p99=round(float(np.median([x['p99'] for x in r if x.get('p99')])), 2))
            out[f'{c}|{fam}'] = s
            print(c.ljust(18), fam.ljust(10), s, flush=True)
    (HERE/'results'/f'bank_{a.label}.json').write_text(json.dumps(dict(summary=out, rows=rows), indent=1), encoding='utf8')


if __name__ == '__main__':
    main()
