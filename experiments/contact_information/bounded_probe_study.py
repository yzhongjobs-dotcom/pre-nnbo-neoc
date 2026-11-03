"""Compare two coordinate-consistent implementations of the SAME sine probe.

Original increment: Delta x / k_current. Proposed offset: p=x/k, command
increment Delta p. A phase-held version changes k only at a sampled zero
crossing. All have the original amplitude parameter, carrier and controller.
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

HERE = Path(__file__).resolve().parent


class BoundedProbe(JointPolicy):
    def __init__(self, normalization='offset', **kw):
        if normalization not in ('legacy', 'offset', 'phase'):
            raise ValueError('normalization')
        super().__init__(rate=False, **kw)
        self.normalization = normalization

    def reset(self, p):
        super().reset(p)
        self.offset_previous = 0.
        self.held_gain = self.k_seek
        self.coordinate_trace = []

    def step(self, o):
        old_x = self.x_prev
        old_offset = self.offset_previous
        command = super().step(o)
        t = float(o.log[o.n-1, 0])
        k = self.k_good if np.isfinite(self.k_good) else self.k_seek
        active = self.params['excitation'] and self.mode == 0 and o.log[o.n-1, 5] >= .3
        enabled = t >= self.warm+.064+self.information_window*.008
        if not active:
            self.offset_previous = 0.
            self.held_gain = k
            return command
        if not enabled or self.normalization != 'phase' or old_x*self.x_prev <= 0:
            self.held_gain = k
        offset = self.x_prev/self.held_gain
        legacy_increment = (self.x_prev-old_x)/k
        new_increment = offset-old_offset
        used_increment = legacy_increment
        if enabled and self.normalization != 'legacy':
            # Change only the deliberate excitation contribution; leave the
            # full feedback/feedforward computation and state updates intact.
            command += new_increment-legacy_increment
            used_increment = new_increment
        self.offset_previous = offset
        self.coordinate_trace.append([t, float(k), float(self.held_gain), float(self.x_prev),
                                      offset, legacy_increment, new_increment, used_increment,
                                      float(enabled)])
        return command


def episode(job):
    scene, mechanism = job
    family, sid, p, T = scene
    policy = BoundedProbe(normalization=mechanism, single_model=False, **CE.FINAL)
    log, _ = common.run(MeasuredOnly(policy), p, T, sid, platform=dict(common.P.PLATFORMS['PC125']))
    frozen = np.load(next((HERE/'results').glob(f'lean_input_{family}_*.npz')))
    startup = log[:, 0] < 1.32
    np.testing.assert_array_equal(log[startup, 11], frozen['measured'][startup, 4])
    if mechanism == 'legacy':
        np.testing.assert_array_equal(log[:, [0, 1, 2, 5, 11, 12, 13]], frozen['measured'])
    trace = np.asarray(policy.coordinate_trace)
    rows = trace[:, 8] > 0
    # On uninterrupted active segments the exact offset increments telescope.
    # The trace retains increments as emitted, so this is auditable separately.
    telescope = float(np.max(np.abs(np.cumsum(trace[rows, 7])-np.cumsum(trace[rows, 6]))))
    name = HERE/'results'/f'bounded_input_{family}_{mechanism}_v1.npz'
    if name.exists():
        raise ValueError('preserve raw result')
    np.savez_compressed(name, measured=log[:, [0, 1, 2, 5, 11, 12, 13]], coordinate_trace=trace,
                        controller_trace=np.asarray(policy.trace))
    return dict(family=family, sid=sid, normalization=mechanism, metrics=common.task_metrics(log, 20.8),
                accepted_updates=policy.used, fits=policy.contact_fits,
                exact_common_startup=True, emitted_minus_offset_cumulative_max_m=telescope,
                measured_file=name.name, measured_sha256=hashlib.sha256(name.read_bytes()).hexdigest(),
                call_ms_P50_P90_P99=np.percentile(np.asarray(policy.t_call)[126:]*1000, [50, 90, 99]).tolist())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', default='bounded_probe_v1')
    parser.add_argument('--families', nargs='+', default=['S2', 'T3m'])
    parser.add_argument('--mechanisms', nargs='+', choices=['legacy', 'offset', 'phase'],
                        default=['legacy', 'offset', 'phase'])
    args = parser.parse_args()
    output = HERE/'results'/(args.label+'.json')
    if output.exists():
        raise ValueError('preserve original result')
    result = dict(protocol=dict(question='does removing unintended probe coordinate drift improve the original method?',
        expected_rows=len(args.families)*len(args.mechanisms), original_frequency=20.8,
        original_force_amplitude=CE.FINAL['a_force'], no_retuning=True, unchanged_plant=True,
        unchanged_gain_and_motion_estimators=True, final_used=False,
        common_startup_seconds=1.32, computational_delay_simulated=False,
        physical_force_budget_not_guaranteed=True,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__), HERE/'joint_policy.py', HERE/'local_contact.py')}), rows=[])
    scenes = [s for s in CE.cohort('diag', 1) if s[0] in args.families]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(episode, [(s, m) for s in scenes for m in args.mechanisms]):
            result['rows'].append(row)
            output.write_text(json.dumps(result, indent=2), encoding='utf8')
            print(row['family'], row['normalization'], 'RMSE', row['metrics']['rmse'],
                  'peak', row['metrics']['peak'], 'accepted', row['accepted_updates'], flush=True)


if __name__ == '__main__':
    main()
