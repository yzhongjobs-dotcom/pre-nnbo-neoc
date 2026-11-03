"""Completed checkpoint: numeric regression, paired baseline and result integrity."""
import hashlib
import json
from pathlib import Path
import numpy as np
import contact_jet as J
import memory_contact as M

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'


def main():
    replay_names = ['contact_jet_replay_v1', 'contact_jet_replay_centered_v2', 'contact_memory_replay_v1',
                    'contact_memory_constrained_v2', 'contact_memory_drift_ablation_v3']
    checked = 0
    for name in replay_names:
        result = json.loads((OUT/(name+'.json')).read_text(encoding='utf8'))
        for path, expected in result['protocol']['input_hashes'].items():
            assert hashlib.sha256((OUT/path).read_bytes()).hexdigest() == expected
        for summary in result['summary']:
            rows = [r for r in result['rows'] if (r['family'], r['model']) == (summary['family'], summary['model'])]
            assert len(rows) == summary['windows']
            path = next(OUT.glob(f'lean_input_{summary["family"]}_*.npz'))
            data = np.load(path); measured = data['measured']
            module = M if summary['model'] in M.TIMES else J
            # After local implementation revisions, reproduce representative
            # early/central/late saved points without creating new evidence n.
            for row in [rows[0], rows[len(rows)//2], rows[-1]]:
                end = row['end']; tau = float(data['causal_tau'][end-1])
                r = module.fit(measured[:end], tau=tau, model=summary['model'])
                np.testing.assert_allclose(r['k'], row['k'], rtol=1e-7, atol=1e-5)
                for h in (1, 4):
                    prediction, _ = module.predict(r, measured, end, h, tau, summary['model'])
                    error = prediction-measured[end-1+h, 3]
                    np.testing.assert_allclose(error, row[f'h{h}']['error'], rtol=1e-7, atol=1e-7)
                checked += 1
    curve = json.loads((OUT/'contact_curvature_closed_v1.json').read_text(encoding='utf8'))
    memory = json.loads((OUT/'contact_memory_gain_closed_v1.json').read_text(encoding='utf8'))
    assert len(curve['rows']) == curve['protocol']['expected_rows'] == 8
    assert len(memory['rows']) == memory['protocol']['expected_rows'] == 4
    baseline = {r['family']:r for r in curve['rows'] if r['model'] == 'lean'}
    for family, row in baseline.items():
        metadata = json.loads(str(np.load(next(OUT.glob(f'lean_input_{family}_*.npz')))['metadata']))
        assert metadata['metrics'] == row['metrics']
    gated = next(r for r in curve['rows'] if r['family'] == 'T3m' and r['model'] == 'centered_resolved')
    assert gated['metrics'] == baseline['T3m']['metrics'] and gated['used'] == 0
    result = dict(saved_replay_points_reproduced=checked, input_hashes_verified=True,
                  lean_closed_loop_equals_frozen_measurements=True,
                  empty_curvature_gate_preserves_lean_exactly=True, complete_closed_loop_rows=12,
                  closed_loop_summary=[dict(family=r['family'], model=r['model'], rmse=r['metrics']['rmse'],
                       peak=r['metrics']['peak'], used=r['used'], callback_p99_ms=r['call_ms_P50_P90_P99'][2])
                       for r in curve['rows']+memory['rows']],
                  metadata_note='memory_gain_closed_v1 inherited curvature resolution_rule text incorrectly; actual code used positive SLS coefficients plus inherited gain gate, NOT a curvature test. Raw result and original metadata retained; current runner fixes the description.',
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       [HERE/'contact_jet.py', HERE/'memory_contact.py', HERE/'run_refinement.py']})
    (OUT/'contact_refinement_audit.json').write_text(json.dumps(result, indent=2), encoding='utf8')
    print('refinement audit passed:', checked, 'replayed numeric points; 12 complete closed-loop rows')


if __name__ == '__main__':
    main()
