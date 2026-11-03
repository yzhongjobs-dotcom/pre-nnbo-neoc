"""Task value of the existing probe, with the preserved lean estimator fixed.

Two interventions, not an amplitude optimizer: half/off after a common one-
second startup. Same environment, noise, gain estimator and upper control.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from run_joint import JointPolicy, MeasuredOnly, common, CE

HERE = Path(__file__).resolve().parent


class MatchedStartupProbe(JointPolicy):
    """Use the original arithmetic path through the entire common startup."""
    def step(self, o):
        if self.amplitude_scale == 1.:
            return super().step(o)
        if o.log[o.n-1, 0] < self.warm:
            scale = self.amplitude_scale
            self.amplitude_scale = 1.
            try:
                command = super().step(o)
            finally:
                self.amplitude_scale = scale
            self.probe_previous = self.x_prev
            return command
        # In the replacement branch the parent probe is disabled; the separate
        # probe_previous register carries the true prior probe level.
        self.x_prev = 0.
        return super().step(o)


def episode(job):
    scene, scale = job
    family, sid, p, T = scene
    pol = MatchedStartupProbe(rate=False, amplitude_scale=scale, single_model=False, **CE.FINAL)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    metrics = common.task_metrics(log, 20.8)
    mask = log[:, 0] >= 1.5
    t = log[mask, 0]; force = log[mask, 3]
    B = np.column_stack([np.ones(len(t)), t-t.mean()])
    tones = np.column_stack([np.sin(2*np.pi*20.8*t), np.cos(2*np.pi*20.8*t)])
    coef = np.linalg.lstsq(np.column_stack([B, tones]), force, rcond=None)[0]
    metrics['probe_coherent_force_rms'] = float(np.sqrt(np.mean((tones@coef[2:])**2)))
    fits = pol.contact_fits
    return dict(family=family, sid=sid, amplitude_scale_after_startup=scale, metrics=metrics,
        accepted_updates=pol.used, estimates=len(fits),
        median_conditional_relative_sd=float(np.median([r['relative_sd'] for r in fits])),
        fits=fits,
        common_startup_commands=log[log[:, 0] < 1., 11].tolist(),
        common_startup_commands_sha256=hashlib.sha256(np.ascontiguousarray(log[log[:, 0] < 1., 11]).tobytes()).hexdigest(),
        call_ms_P50_P90_P99=np.percentile(np.asarray(pol.t_call)[126:]*1000, [50, 90, 99]).tolist())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default='probe_value_intervention_v2')
    args = parser.parse_args()
    target = HERE/'results'/(args.label+'.json')
    if target.exists():
        raise ValueError('preserve original result')
    scenes = [s for s in CE.cohort('diag', 1) if s[0] in ('S2', 'T3m')]
    result = dict(protocol=dict(question='marginal contact-information task benefit versus physical force ripple',
        startup_seconds=1., transition_seconds=.064, families=['S2', 'T3m'],
        scales=[1., .5, 0.], no_gain_or_control_retuning=True, no_plant_changes=True,
        expected_rows=6, truth_for_scoring_only=True, final_used=False,
        computational_delay_simulated=False,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__), HERE/'joint_policy.py', HERE/'local_contact.py', HERE/'run_joint.py')}), rows=[])
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(episode, [(s, scale) for s in scenes for scale in (1., .5, 0.)]):
            result['rows'].append(row)
            target.write_text(json.dumps(result, indent=2), encoding='utf8')
            print(row['family'], row['amplitude_scale_after_startup'], 'RMSE', row['metrics']['rmse'],
                  'probe RMS', row['metrics']['probe_coherent_force_rms'], 'accepted', row['accepted_updates'], flush=True)
    for family in ('S2', 'T3m'):
        rows = [r for r in result['rows'] if r['family'] == family]
        assert len(set(r['common_startup_commands_sha256'] for r in rows)) == 1
        frozen = np.load(next((HERE/'results').glob(f'lean_input_{family}_*.npz')))
        assert rows[0]['metrics']['rmse'] == json.loads(str(frozen['metadata']))['metrics']['rmse']
    result['verified_common_startup_and_exact_baseline'] = True
    target.write_text(json.dumps(result, indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
