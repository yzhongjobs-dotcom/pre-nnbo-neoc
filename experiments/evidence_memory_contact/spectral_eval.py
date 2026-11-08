"""Paired closed-loop evaluation on the diagnostic cohort. Truth is scoring only."""
import argparse, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
HERE = Path(__file__).resolve().parent
for sub in ('contact_information', 'motion_aware_scan', 'method_redesign', 'conditional_linear_contact'):
    sys.path.insert(0, str(HERE.parent/sub))
import common
from cohort_eval import FINAL, FAMS, cohort
from v1_servo_id import ServoID
from joint_policy import JointPolicy
from conditional_filter import MeasuredOnly
from spectral_policy import SpectralContact

CONFIGS = {'v1': None, 'lean': None}
for s in (1., .5):
    for nz in (.04, .028):
        CONFIGS[f's{s}_n{nz}'] = dict(change_prior=s, fit_noise=nz)
CONFIGS.update({
    's1_n.028_chi': dict(change_prior=1., fit_noise=.028, chi_max=FINAL['chi2_max']),
    's1_n.028_fb': dict(change_prior=1., fit_noise=.028, feedback=True),
    's.5_n.028_fb': dict(change_prior=.5, fit_noise=.028, feedback=True),
    's1_n.028_a.7': dict(change_prior=1., fit_noise=.028, amplitude_scale=.7),
    's1_n.028_a.5': dict(change_prior=1., fit_noise=.028, amplitude_scale=.5),
    's.5_n.028_a.7': dict(change_prior=.5, fit_noise=.028, amplitude_scale=.7),
    's.5_n.028_a.5': dict(change_prior=.5, fit_noise=.028, amplitude_scale=.5),
    'lean_a.7': 'lean_a.7', 'v1_a.7': 'v1_a.7',
})


def job(a):
    name, scene, lag = a; fam, sid, p, T = scene
    c = CONFIGS[name]
    if name == 'v1':
        pol = ServoID(single_model=False, **FINAL)
    elif name == 'lean':
        pol = JointPolicy(wave=None, new_fit=True, rate=False, amplitude_scale=1., **FINAL)
    elif name == 'lean_a.7':
        pol = JointPolicy(wave=None, new_fit=True, rate=False, amplitude_scale=.7, **FINAL)
    elif name == 'v1_a.7':
        pol = JointPolicy(wave=None, new_fit=False, rate=False, amplitude_scale=.7, **FINAL)
    else:
        pol = SpectralContact(single_model=False, **c, **FINAL)
    plat = dict(common.P.PLATFORMS['PC125']); plat['tau_a'] *= lag
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=plat)
    if log is None:
        return dict(config=name, family=fam, sid=sid, rmse=None, peak=None, loss=None)
    m = common.task_metrics(log, 20.8)
    calls = np.asarray(pol.t_call[126:])*1000
    tr = np.asarray(getattr(pol, 'spectral_trace', []) or [[0, 0, 0, 0, 0, 0]])
    return dict(config=name, family=fam, sid=sid, rmse=float(m['rmse']), peak=float(m['peak']), loss=float(m['loss_max']), ripple=float(m['ripple']),
                accept=float(tr[:, 5].mean()), corner=float(tr[:, 3].mean()), p50=float(np.percentile(calls, 50)), p99=float(np.percentile(calls, 99)))


def boot(d, B=4000):
    rng = np.random.default_rng(7); d = np.asarray(d)
    m = d[rng.integers(0, len(d), (B, len(d)))].mean(axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--configs', nargs='+', required=True); ap.add_argument('--scenes', type=int, default=50)
    ap.add_argument('--lag', type=float, default=1.); ap.add_argument('--label', required=True); ap.add_argument('--ref', default='v1')
    a = ap.parse_args()
    cfgs = list(dict.fromkeys(['v1', 'lean']+a.configs))
    scenes = cohort('diag', a.scenes)
    with ProcessPoolExecutor(max_workers=15) as e:
        rows = list(e.map(job, [(c, s, a.lag) for c in cfgs for s in scenes], chunksize=3))
    out = {}
    for c in cfgs:
        print('---', c)
        for fam in FAMS:
            r = [x for x in rows if x['config'] == c and x['family'] == fam]
            ok = [x['rmse'] for x in r if x['rmse'] is not None]
            def paired(ref):
                b = {x['sid']: x for x in rows if x['config'] == ref and x['family'] == fam}
                d = [x['rmse']-b[x['sid']]['rmse'] for x in r if x['rmse'] is not None and b[x['sid']]['rmse'] is not None]
                return round(float(np.mean(d)), 3), [round(v, 3) for v in boot(d)]
            s = dict(med=round(float(np.median(ok)), 3), mean=round(float(np.mean(ok)), 3), r=sum(1 for x in r if x['rmse'] is None or x['peak'] > 50),
                     L=sum(1 for x in r if x['loss'] and x['loss'] > .5), vs_v1=paired('v1'), vs_lean=paired('lean'),
                     acc=round(float(np.mean([x.get('accept', 0) for x in r])), 2), corner=round(float(np.mean([x.get('corner', 0) for x in r])), 1),
                     ripple=round(float(np.median([x.get('ripple', 0) for x in r])), 3), p99=round(float(np.median([x.get('p99', 0) for x in r])), 2))
            out[f'{c}|{fam}'] = s
            print(f'{fam:10s}', s, flush=True)
    (HERE/'results'/f'closed_{a.label}.json').write_text(json.dumps(dict(summary=out, rows=rows, lag=a.lag), indent=1), encoding='utf8')


if __name__ == '__main__':
    main()
