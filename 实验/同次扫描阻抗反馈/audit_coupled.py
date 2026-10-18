"""Frozen-policy transfer: paired uncertainty and a targeted mesh refinement.

The changed plant is a sensitivity study. Policies are NOT retuned here; these
results do not compare the best controllers achievable on the corrected plant.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import coupled_contact as cm
from run import noises

HERE = Path(__file__).resolve().parent


def comparison(a, b):
    rng = np.random.default_rng(8691)
    ids = rng.integers(0, len(a), size=(10000, len(a)))
    benefit = 1 - a[ids, 0].mean(1) / b[ids, 0].mean(1)
    return {
        'relative_mean_J_benefit': float(1 - a[:, 0].mean() / b[:, 0].mean()),
        'paired_95pct': np.quantile(benefit, [.025, .975]).tolist(),
        'n': len(a),
    }


def main():
    data = np.load(HERE / 'coupled_paired.npz')
    final = json.loads((HERE / 'final_policies.json').read_text())
    pure = json.loads((HERE / 'pure_gain_fit.json').read_text())
    timing = json.loads((HERE / 'results.json').read_text())['timing']
    comparisons = {}
    for candidate in ['learned_feedback', 'learned_gain_only', 'spatial_joint']:
        for baseline in ['analytic_inverse', 'continuous_gain_only']:
            comparisons[candidate + '_vs_' + baseline] = comparison(data[candidate], data[baseline])
    ps = data['test_params']
    ns = noises(ps, 8602)
    # Both first cases and the largest pure-gain loss; post-hoc diagnostics only.
    worst = int(np.argmax(data['learned_gain_only'][:, 0] - data['continuous_gain_only'][:, 0]))
    rows = []
    for case in sorted(set([0, 2, worst])):
        for name in ['continuous_gain_only', 'learned_gain_only', 'learned_feedback']:
            theta = np.array((final if name == 'learned_feedback' else pure)['policies'][name])
            latency = timing['learned_feedback']['p95_ms'] / 1000
            metric, _ = cm.episode(ps[case], ns[case], 12, np.array(final['cfg']),
                                   theta, latency, dt=.00025, nq=257)
            old = float(data[name][case, 0])
            row = {'case': case, 'family': case % 3, 'method': name,
                   'original_J': old, 'fine_J': float(metric[0]),
                   'relative_difference': float(metric[0] / old - 1)}
            rows.append(row)
            print(row, flush=True)
    out = {'comparisons': comparisons, 'refinement': rows,
           'scope': 'Frozen-controller sensitivity only; 24 fresh cases, not jointly retuned controllers or validated real tissue. Refinement cases are not new independent data.',
           'source_sha256': {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest()
                             for f in ['model.py', 'coupled_contact.py', 'audit_coupled.py',
                                       'coupled_protocol.json', 'coupled_paired.npz']}}
    (HERE / 'coupled_audit.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(comparisons, indent=2), flush=True)


if __name__ == '__main__':
    main()
