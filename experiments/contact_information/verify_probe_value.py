"""Algebra and paired-result checks for the probe-value checkpoint."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_probe_information import instrument_fit
from local_contact import features, fit

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'


def main():
    checks = []
    t = np.arange(80)*.008; n = 64
    for extra in (0., .00004):
        y = .002+.00003*np.sin(2*np.pi*20.8*t)+.0004*t+extra*np.cos(2*np.pi*9.1*t)
        measured = np.c_[t, .02*t, y, np.full(len(t), 4.), np.zeros((len(t), 3))]
        j = np.arange(n-33, n-1)
        B, U = features(t[j], y[j], y[j+1], t[n-1], y[n-1], .256)
        theta = np.array([12000., 2000.]); beta = np.array([4., .02, -.01, 0., 0.])
        z = U[:, :2]@theta+B@beta; decay = np.exp(-2*np.pi*30*.008)
        for index, value in zip(j, z):
            measured[index+1, 3] = decay*measured[index, 3]+value
        r = instrument_fit(measured[:n])
        ordinary = fit(measured[:n], nw=32, rate=False)
        np.testing.assert_allclose(r['coefficients'], theta, rtol=1e-8, atol=1e-5)
        np.testing.assert_allclose(ordinary['coefficients'], theta, rtol=1e-8, atol=1e-5)
        assert 0 < r['probe_fraction'] <= 1+1e-9
        if extra == 0:
            np.testing.assert_allclose(r['probe_fraction'], 1., atol=1e-9)
        else:
            assert r['probe_fraction'] < .99
        checks.append(dict(extra_nonprobe_displacement=extra, both_estimators_recover=True,
                           probe_information_fraction=r['probe_fraction']))
    audit = json.loads((OUT/'probe_information_actual_command_v2.json').read_text())
    for path, expected in audit['protocol']['inputs_sha256'].items():
        assert hashlib.sha256((OUT/path).read_bytes()).hexdigest() == expected
    assert all(0 < r['probe_fraction'] <= 1+1e-8 for r in audit['rows'])
    assert all(r['max_command_error'] == 0 for r in audit['protocol']['replayed_traces'].values())
    intervention = json.loads((OUT/'probe_value_intervention_v2.json').read_text())
    assert len(intervention['rows']) == intervention['protocol']['expected_rows']
    assert intervention['verified_common_startup_and_exact_baseline']
    summary = [dict(family=r['family'], mechanism='fixed_'+str(r['amplitude_scale_after_startup']),
                    rmse=r['metrics']['rmse'], peak=r['metrics']['peak'], accepted=r['accepted_updates'])
               for r in intervention['rows']]
    result = dict(algebra=checks, valid_inputs_and_conditional_information_order=True,
                  exact_causal_command_replay=True, complete_closed_loop_rows=len(intervention['rows']), summary=summary,
                  timing_is_python_prototype_not_real_time_certification=True,
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__), HERE/'audit_probe_information.py', HERE/'probe_value_intervention.py')})
    (OUT/'probe_value_checkpoint.json').write_text(json.dumps(result, indent=2), encoding='utf8')
    print('probe-value verification passed: exact replay, information algebra and', len(intervention['rows']), 'paired rows')


if __name__ == '__main__':
    main()
