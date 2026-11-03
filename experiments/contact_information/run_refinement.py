"""Targeted contact-model refinement on unchanged, reused diagnosis conditions."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
import math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from run_joint import MeasuredOnly, common, CE, JointPolicy
from contact_jet import MODELS, fit, predict, SCALE
from curved_policy import CurvedContactPolicy
from memory_policy import MemoryContactPolicy
import memory_contact

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'


class RecordingMeasurements(MeasuredOnly):
    def reset(self, p):
        super().reset(p)
        self.taus = []

    def __call__(self, o):
        value = super().__call__(o)
        if o.n:
            self.taus.append((o.n-1, self.inner.tau))
        return value


def freeze(scene):
    fam, sid, p, T = scene
    path = OUT/f'lean_input_{fam}_{sid}.npz'
    if path.exists():
        return str(path)
    policy = JointPolicy(rate=False, single_model=False, **CE.FINAL)
    wrapper = RecordingMeasurements(policy)
    log, _ = common.run(wrapper, p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    taus = np.full(len(log), np.nan)
    for index, tau in wrapper.taus:
        taus[index] = tau
    # Distinct arrays: estimator/replay are passed ONLY measured and causal_tau.
    np.savez_compressed(path, measured=log[:, [0, 1, 2, 5, 11, 12, 13]], causal_tau=taus,
                        scoring_reference_k=log[:, 8], scoring_true_force=log[:, 3],
                        metadata=json.dumps(dict(family=fam, sid=sid, metrics=common.task_metrics(log, 20.8),
                            policy='JointPolicy(rate=False), original 20.8Hz excitation', plant_changed=False)))
    return str(path)


def replay(label, models=MODELS):
    target = OUT/(label+'.json')
    if target.exists():
        raise ValueError('preserve primary result; choose a new label')
    rows, summaries, inputs = [], [], {}
    for path in sorted(OUT.glob('lean_input_*.npz')):
        inputs[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        data = np.load(path)
        m, tau_trace = data['measured'], data['causal_tau']
        fam = json.loads(str(data['metadata']))['family']
        for end in range(192, len(m)-4, 8):
            tau = float(tau_trace[end-1])
            if not np.isfinite(tau):
                continue
            for model in models:
                module = memory_contact if model in memory_contact.TIMES else None
                r = (module.fit if module else fit)(m[:end], tau=tau, model=model)
                if r is None:
                    continue
                accepted = bool(200 <= r['k'] <= 3e5 and r['min_gain'] > 0 and
                                r['risk'] <= math.log(1.2) and r['chi2'] <= CE.FINAL['chi2_max'])
                rr = dict(family=fam, end=end, time=float(m[end-1, 0]), model=model,
                          k=r['k'], sigma=r['sigma'], risk=r['risk'], accepted=accepted,
                          chi2=r['chi2'], condition=r['condition'], residual_rms=r['residual_rms'],
                          coefficients=r['coefficients'].tolist(), coefficient_sd=np.sqrt(np.diag(r['cov'])).tolist())
                if module:
                    rr.update(physical=r['physical'], low_frequency_gain=r['low_frequency_gain'],
                              low_frequency_sd=r['low_frequency_sd'], quadrature_gain=r['quadrature_gain'])
                if model in ('normal_curvature', 'centered_curvature'):
                    # How strongly can this trajectory constrain curvature over its observed excursion?
                    span = float(np.ptp(m[end-33:end, 2]))
                    if r['center'] is not None:
                        ix = np.arange(end-33, end)
                        displacement = m[ix, 2]-r['ye']-np.polynomial.polynomial.polyval((m[ix, 0]-r['te'])/r['tw'], r['center'])
                        span = float(np.ptp(displacement))
                    rr['curvature_force_effect'] = float(.5*r['coefficients'][2]/SCALE*span**2)
                    rr['curvature_force_sd'] = float(.5*np.sqrt(r['cov'][2, 2])/SCALE*span**2)
                for h in (1, 4):
                    fp, var = (module.predict if module else predict)(r, m, end, h, tau, model)
                    error = fp-float(m[end-1+h, 3])
                    rr[f'h{h}'] = dict(error=error, sd=math.sqrt(var), nll=.5*(math.log(2*math.pi*var)+error**2/var))
                rows.append(rr)
        for model in models:
            rr = [r for r in rows if r['family'] == fam and r['model'] == model]
            summary = dict(family=fam, model=model, windows=len(rr), accepted=sum(r['accepted'] for r in rr),
                           sigma_k_median=float(np.median([r['sigma'] for r in rr])),
                           conditional={f'h{h}':dict(rmse=float(np.sqrt(np.mean([r[f'h{h}']['error']**2 for r in rr]))),
                               nll=float(np.mean([r[f'h{h}']['nll'] for r in rr])),
                               coverage95=float(np.mean([abs(r[f'h{h}']['error']) <= 1.96*r[f'h{h}']['sd'] for r in rr]))) for h in (1, 4)})
            if model in ('normal_curvature', 'centered_curvature'):
                summary['curvature_resolved_fraction'] = float(np.mean([abs(r['coefficients'][2]) > 1.96*r['coefficient_sd'][2] for r in rr]))
                summary['curvature_effect_abs_median'] = float(np.median([abs(r['curvature_force_effect']) for r in rr]))
                summary['curvature_effect_sd_median'] = float(np.median([r['curvature_force_sd'] for r in rr]))
            if model in memory_contact.TIMES:
                summary['physical_fraction'] = float(np.mean([r['physical'] for r in rr]))
                summary['low_frequency_gain_median'] = float(np.median([r['low_frequency_gain'] for r in rr]))
                summary['low_frequency_sd_median'] = float(np.median([r['low_frequency_sd'] for r in rr]))
            summaries.append(summary)
            print(fam, model, json.dumps(summary), flush=True)
    target.write_text(json.dumps(dict(protocol=dict(input_hashes=inputs, source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in (Path(__file__), HERE/'contact_jet.py', HERE/'memory_contact.py')},
                  stride=8, horizons_ms=[8, 32], windows_independent=False, conditioning='future encoder only; no future forces',
                  environment_unchanged=True, baseline='effective lean candidate', tuning=False), summary=summaries, rows=rows), indent=2), encoding='utf8')


def closed_episode(job):
    scene, name = job
    fam, sid, p, T = scene
    if name == 'lean':
        pol = JointPolicy(rate=False, single_model=False, **CE.FINAL)
    elif name in memory_contact.TIMES:
        pol = MemoryContactPolicy(memory_model=name, single_model=False, **CE.FINAL)
    else:
        pol = CurvedContactPolicy(model='normal_curvature' if name == 'normal_raw' else 'centered_curvature',
                                  require_resolved=name == 'centered_resolved', single_model=False, **CE.FINAL)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    return dict(family=fam, sid=sid, model=name, metrics=common.task_metrics(log, 20.8),
                used=getattr(pol, 'curvature_used', getattr(pol, 'memory_used', 0)),
                fits=getattr(pol, 'curvature_fits', getattr(pol, 'memory_fits', [])),
                call_ms_P50_P90_P99=np.percentile(np.asarray(pol.t_call)[126:]*1000, [50, 90, 99]).tolist())


def closed(label, memory=False):
    target = OUT/(label+'.json')
    if target.exists():
        raise ValueError('preserve primary result; choose a new label')
    scenes = [s for s in CE.cohort('diag', 1) if s[0] in ('T3m', 'T2m-stiff')]
    names = ('memoryc_040', 'drift_040') if memory else ('lean', 'normal_raw', 'centered_raw', 'centered_resolved')
    jobs = [(scene, name) for scene in scenes for name in names]
    result = dict(protocol=dict(environment_unchanged=True, parameters=CE.FINAL, seed_pairing=True,
                               baseline='effective lean candidate', new_model_covariance='conditional only',
                               integration=('in-phase 20.8Hz gain only; no dynamic compliance inverse' if memory else
                                            'current endpoint tangent gain, original motion conversion'),
                               resolution_rule=('positive endpoint moduli plus inherited gain-risk gate' if memory else
                                                'positive H / conditional SD > 1.96 for resolved variant only'),
                               expected_rows=len(jobs), final_used=False, computational_delay_simulated=False,
                               source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                   (Path(__file__), HERE/'curved_policy.py', HERE/'contact_jet.py', HERE/'joint_policy.py',
                                    HERE/'memory_policy.py', HERE/'memory_contact.py')}), rows=[])
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(closed_episode, jobs):
            result['rows'].append(row)
            target.write_text(json.dumps(result, indent=2), encoding='utf8')
            print(row['family'], row['model'], 'RMSE', row['metrics']['rmse'], 'peak', row['metrics']['peak'],
                  'candidate_updates', row['used'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['freeze', 'replay', 'closed', 'memory_closed'])
    parser.add_argument('--label', default='contact_jet_replay_v1')
    parser.add_argument('--models', nargs='+', choices=list(MODELS)+list(memory_contact.TIMES), default=list(MODELS))
    args = parser.parse_args()
    if args.action == 'freeze':
        scenes = [s for s in CE.cohort('diag', 1) if s[0] in ('S2', 'T3m', 'T2m-stiff')]
        with ProcessPoolExecutor(max_workers=2) as pool:
            for name in pool.map(freeze, scenes):
                print(name, flush=True)
    elif args.action == 'replay':
        replay(args.label, args.models)
    else:
        closed(args.label, memory=args.action == 'memory_closed')
