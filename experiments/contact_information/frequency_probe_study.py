"""A physical-frequency hypothesis, not an unconstrained waveform search.

Lower the carrier to the nominal servo corner. Match the predicted physical
force ripple by reducing command amplitude, and compare with the SAME command
amplitude at the original frequency. Only measured channels reach the policy.
"""
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
from run_joint import JointPolicy, MeasuredOnly, common, CE
from local_contact import sine_information, sine_reference_rms

HERE = Path(__file__).resolve().parent
TAU = .019
CORNER = 1/(2*math.pi*TAU)
SCALE = sine_reference_rms(20.8)/sine_reference_rms(CORNER)


class FrequencyProbe(JointPolicy):
    def __init__(self, carrier=20.8, scale=1., **kw):
        self.carrier = float(carrier)
        wave = None if carrier == 20.8 else np.sin(2*np.pi*carrier*.008*np.arange(32768))
        super().__init__(rate=False, wave=wave, amplitude_scale=scale, **kw)

    def step(self, o):
        if self.wave is None and self.amplitude_scale == 1.:
            return super().step(o)
        if o.log[o.n-1, 0] < self.warm:
            wave, scale = self.wave, self.amplitude_scale
            self.wave, self.amplitude_scale = None, 1.
            try:
                command = super().step(o)
            finally:
                self.wave, self.amplitude_scale = wave, scale
            self.probe_previous = self.x_prev
            return command
        self.x_prev = 0.
        return super().step(o)


def design():
    rows = []
    reference = sine_reference_rms()
    # A design calculation, not repeated physical trials. The carrier at the
    # servo corner is selected from hardware, before examining task outcomes.
    for frequency in sorted(set([4., 6., CORNER, 12., 16., 20.8, 24., 28.])):
        rms = sine_reference_rms(frequency)
        scale = min(1., reference/rms,
                    math.sin(math.pi*20.8*.008)/math.sin(math.pi*frequency*.008),
                    (math.sin(math.pi*20.8*.008)/math.sin(math.pi*frequency*.008))**2)
        info = sine_information(frequency=frequency, n=32, rate=False)[0]*scale**2
        rows.append(dict(frequency=frequency, command_scale=scale,
                         nominal_physical_force_rms=rms*scale,
                         worst_phase_endpoint_information=info))
    return dict(rows=rows, servo_corner=CORNER, corner_scale=SCALE,
                basis='nominal 19 ms position servo; same nominal physical ripple, command peak/slew/second-difference limits',
                limitation='linear nominal gain/servo only; no actual-force safety guarantee or contact-memory identification',
                endpoint_model='unchanged 32-tick gain+slope; B4 nuisance; phase-grid minimum')


def episode(job):
    scene, mechanism = job
    family, sid, p, T = scene
    carrier, scale = ((20.8, 1.) if mechanism == 'baseline' else
                      (CORNER, SCALE) if mechanism == 'servo_corner' else (20.8, SCALE))
    pol = FrequencyProbe(carrier=carrier, scale=scale, single_model=False, **CE.FINAL)
    log, _ = common.run(MeasuredOnly(pol), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    frozen = np.load(next((HERE/'results').glob(f'lean_input_{family}_*.npz')))
    warm = log[:, 0] < 1.
    np.testing.assert_array_equal(log[warm, 11], frozen['measured'][warm, 4])
    if mechanism == 'baseline':
        np.testing.assert_array_equal(log[:, [0, 1, 2, 5, 11, 12, 13]], frozen['measured'])
    metrics = common.task_metrics(log, carrier)
    trace = np.asarray(pol.trace)
    row = dict(family=family, sid=sid, mechanism=mechanism, frequency=carrier, scale=scale,
               metrics=metrics, accepted_updates=pol.used, fits=pol.contact_fits,
               exact_common_startup=True, baseline_exactly_reproduced=mechanism == 'baseline',
               median_conditional_relative_sd=float(np.median([f['relative_sd'] for f in pol.contact_fits])))
    # Raw observation trajectories: scoring truth is neither saved here nor
    # passed back into the policy. This permits paired future model diagnosis.
    name = HERE/'results'/f'frequency_input_{family}_{mechanism}_v1.npz'
    if name.exists():
        raise ValueError('preserve raw trajectory')
    np.savez_compressed(name, measured=log[:, [0, 1, 2, 5, 11, 12, 13]], trace=trace,
                        metadata=json.dumps({k:v for k,v in row.items() if k != 'fits'}))
    row['measured_file'] = name.name
    row['measured_sha256'] = hashlib.sha256(name.read_bytes()).hexdigest()
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default='frequency_probe_v1')
    parser.add_argument('--families', nargs='+', default=['S2', 'T3m'])
    parser.add_argument('--mechanisms', nargs='+', default=['baseline', 'servo_corner', 'amplitude_control'],
                        choices=['baseline', 'servo_corner', 'amplitude_control'])
    args = parser.parse_args()
    target = HERE/'results'/(args.label+'.json')
    if target.exists():
        raise ValueError('preserve original result')
    result = dict(protocol=dict(question='does lower-frequency contact information improve force control at a matched nominal ripple budget?',
        design=design(), expected_rows=len(args.families)*len(args.mechanisms),
        same_amplitude_control='separates carrier from command amplitude; its force ripple is intentionally lower',
        unchanged_estimator=True, unchanged_plant=True, no_control_retuning=True,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__), HERE/'joint_policy.py', HERE/'local_contact.py')}), rows=[])
    scenes = [s for s in CE.cohort('diag', 1) if s[0] in args.families]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(episode, [(s, m) for s in scenes for m in args.mechanisms]):
            result['rows'].append(row)
            target.write_text(json.dumps(result, indent=2), encoding='utf8')
            print(row['family'], row['mechanism'], 'RMSE', row['metrics']['rmse'],
                  'peak', row['metrics']['peak'], 'accepted', row['accepted_updates'], flush=True)


if __name__ == '__main__':
    main()
