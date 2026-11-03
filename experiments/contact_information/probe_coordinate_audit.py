"""Check whether force-normalized increments implement a bounded probe offset.

Only seven measured columns, the policy's own probe register and its own gain.
The policy command replay must be EXACT; no plant state is available here.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from run_joint import JointPolicy, common, CE

HERE = Path(__file__).resolve().parent


def replay(measured):
    policy = JointPolicy(rate=False, single_model=False, **CE.FINAL)
    policy.reset(dict(common.P.PLATFORMS['PC125']))
    log = np.full((len(measured), 16), np.nan)
    log[:, [0, 1, 2, 5, 11, 12, 13]] = measured
    own_command = 0.
    rows = []
    for j in range(len(measured)):
        log[j, 11] = measured[j-1, 4] if j else measured[0, 2]
        old = policy.x_prev
        issued = policy(common.P.Obs(log, j+1))
        np.testing.assert_equal(issued, measured[j, 4])
        log[j, 11] = measured[j, 4]
        gain = policy.k_good if np.isfinite(policy.k_good) else policy.k_seek
        delta = (policy.x_prev-old)/gain
        own_command += delta
        intended = policy.x_prev/gain
        rows.append([measured[j, 0], gain, policy.x_prev, delta, own_command, intended,
                     own_command-intended])
    return np.asarray(rows)


def main():
    output = HERE/'results/probe_coordinate_audit_v1.json'
    if output.exists():
        raise ValueError('preserve result')
    result = dict(question='does integrating force-probe increments with a changing gain create an unintended position component?',
                  algebra='sum Delta(x)/k - x/k = -sum x_previous Delta(1/k), for initial x=0',
                  caveat='a mathematical command decomposition, not a plant-state counterfactual or proof of task harm', rows=[])
    for path in sorted((HERE/'results').glob('lean_input_*.npz')):
        data = np.load(path)
        trace = replay(data['measured'])
        mask = trace[:, 0] >= 1.5
        wandering = trace[mask, 6]
        nominal = trace[mask, 5]
        row = dict(family=json.loads(str(data['metadata']))['family'],
                   input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), exact_command_replay=True,
                   intended_offset_rms_m=float(np.sqrt(np.mean(nominal**2))),
                   extra_command_span_m=float(np.ptp(wandering)),
                   extra_command_end_m=float(wandering[-1]), trace=trace.tolist())
        result['rows'].append(row)
        print(row['family'], 'intended RMS um', row['intended_offset_rms_m']*1e6,
              'extra span um', row['extra_command_span_m']*1e6,
              'end um', row['extra_command_end_m']*1e6)
    output.write_text(json.dumps(result, indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
