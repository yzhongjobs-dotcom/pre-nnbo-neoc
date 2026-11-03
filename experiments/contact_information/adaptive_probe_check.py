"""Can the inherited adaptive probe address the measured task tradeoff?

Compare its original short-window precision input with the current lean
estimator's precision. Same law and bounds; not a new excitation theorem.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import json
import hashlib
import math
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from run_joint import JointPolicy, MeasuredOnly, common, CE

HERE = Path(__file__).resolve().parent


class AdaptiveProbe(JointPolicy):
    def __init__(self, precision_source='legacy', **kw):
        super().__init__(rate=False, **kw)
        self.precision_source = precision_source

    def reset(self, p):
        super().reset(p)
        self.amplitude_trace = []

    def _stiffness_est(self, o, n):
        feedback = super()._stiffness_est(o, n)
        if self.precision_source in ('lean', 'both') and self.contact_fits:
            r = self.contact_fits[-1]
            if r['time'] == float(o.log[n-1, 0]):
                # rel is consumed only by the inherited amplitude law after
                # stiffness/feedback selection. Do not modify returned feedback.
                if self.precision_source == 'lean':
                    self.rel = r['relative_sd'] if r['k'] > 0 else float('inf')
                else:
                    eps_feedback = (1-self.params['sector_lo'])/(2*self.params['c'])
                    risk_motion = r['risk'] if r['k'] > 0 and r['start_gain'] > 0 else float('inf')
                    # The inherited law takes rel/eps_feedback as its demand.
                    # Preserve its shape while requiring BOTH feedback precision
                    # and the current motion-conversion acceptance precision.
                    self.rel = eps_feedback*max(self.rel/eps_feedback, risk_motion/math.log(1.2))
        return feedback

    def step(self, o):
        t = float(o.log[o.n-1, 0])
        self.params['adaptive_exc'] = t >= self.warm+.064+self.information_window*.008
        command = super().step(o)
        self.amplitude_trace.append([t, self.exc_gain, self.rel])
        return command


def episode(job):
    scene, source = job
    family, sid, p, T = scene
    pol = AdaptiveProbe(precision_source=source, single_model=False, **CE.FINAL)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    metrics = common.task_metrics(log, 20.8)
    mask = log[:, 0] >= 1.5; t = log[mask, 0]; force = log[mask, 3]
    B = np.column_stack([np.ones(len(t)), t-t.mean()])
    tones = np.column_stack([np.sin(2*np.pi*20.8*t), np.cos(2*np.pi*20.8*t)])
    beta = np.linalg.lstsq(np.column_stack([B, tones]), force, rcond=None)[0]
    metrics['probe_coherent_force_rms'] = float(np.sqrt(np.mean((tones@beta[2:])**2)))
    frozen = np.load(next((HERE/'results').glob(f'lean_input_{family}_*.npz')))
    warm = log[:, 0] < 1.32
    np.testing.assert_array_equal(log[warm, 11], frozen['measured'][warm, 4])
    trace = np.asarray(pol.amplitude_trace)
    return dict(family=family, sid=sid, precision_source=source, metrics=metrics, accepted_updates=pol.used,
                mean_amplitude_scale_after_startup=float(np.mean(trace[trace[:, 0] >= 1.5, 1])),
                amplitude_trace=pol.amplitude_trace, fits=pol.contact_fits,
                startup_exactly_matches_frozen_baseline=True,
                call_ms_P50_P90_P99=np.percentile(np.asarray(pol.t_call)[126:]*1000, [50, 90, 99]).tolist())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sources', nargs='+', choices=['legacy', 'lean', 'both'], default=['legacy', 'lean'])
    parser.add_argument('--label', default='adaptive_probe_alignment_v1')
    args = parser.parse_args()
    target = HERE/'results'/(args.label+'.json')
    if target.exists():
        raise ValueError('preserve primary result')
    scenes = [s for s in CE.cohort('diag', 1) if s[0] in ('S2', 'T3m')]
    result = dict(protocol=dict(question='does inherited amplitude adaptation need alignment with the active estimator?',
        controller_parameters=CE.FINAL, original_law='proposed.py adaptive_exc; sector-based precision, smoothing .1, scale [0,1]',
        motion_gain_confidence_not_certified_by_sector_target=True,
        both_demand='max(original short-window feedback precision / sector eps, lean risk / log(1.2))',
        no_retuning=True, no_plant_changes=True, no_new_frequency_or_waveform=True,
        fixed_baseline_result='probe_value_intervention_v2.json', expected_rows=2*len(args.sources),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__), HERE/'joint_policy.py', HERE/'local_contact.py')}), rows=[])
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(episode, [(s, source) for s in scenes for source in args.sources]):
            result['rows'].append(row)
            target.write_text(json.dumps(result, indent=2), encoding='utf8')
            print(row['family'], row['precision_source'], 'RMSE', row['metrics']['rmse'],
                  'mean amplitude', row['mean_amplitude_scale_after_startup'], 'accepted', row['accepted_updates'], flush=True)


if __name__ == '__main__':
    main()
