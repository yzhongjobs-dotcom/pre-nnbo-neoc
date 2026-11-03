"""Independent algebra, causal replay and paired-result checks for this study."""
import hashlib
import json
from pathlib import Path
import numpy as np
from run_joint import common, CE
from bounded_probe_study import BoundedProbe
from local_contact import features, residualize
from frequency_probe_study import SCALE, CORNER
from response_band_information import bands

HERE = Path(__file__).resolve().parent
OUT = HERE/'results'


def main():
    angle = np.linspace(0,2*np.pi,301)
    x = .167*np.sin(angle)
    k = 12000*(1+.3*np.cos(angle)+.15*np.sin(angle))
    legacy = np.cumsum(np.r_[0.,np.diff(x)]/k)
    extra = legacy-x/k
    independent = -np.cumsum(np.r_[0.,x[:-1]*np.diff(1/k)])
    np.testing.assert_allclose(extra, independent, atol=1e-18)
    constant = np.cumsum(np.r_[0.,np.diff(x)]/12000)
    np.testing.assert_allclose(constant,x/12000,atol=1e-18)
    np.testing.assert_allclose(np.cumsum(np.r_[0.,np.diff(x/k)]),x/k,atol=1e-18)
    assert abs(legacy[-1]) > 1e-7
    # Same nominal force-RMS, lower command peak/slew/second difference.
    assert 0 < SCALE < 1
    assert SCALE*np.sin(np.pi*CORNER*.008) < np.sin(np.pi*20.8*.008)
    assert SCALE*np.sin(np.pi*CORNER*.008)**2 < np.sin(np.pi*20.8*.008)**2
    assert sum(a for f,a in bands('dual')) <= 1+1e-12
    tables = {}
    for name in ('frequency_probe_v1','bounded_probe_v1'):
        d = json.loads((OUT/(name+'.json')).read_text())
        assert len(d['rows']) == d['protocol']['expected_rows']
        assert all(r['exact_common_startup'] for r in d['rows'])
        for r in d['rows']:
            assert hashlib.sha256((OUT/r['measured_file']).read_bytes()).hexdigest() == r['measured_sha256']
        tables[name] = [{k:r[k] for k in ('family', 'normalization' if name.startswith('bounded') else 'mechanism')}
                        | dict(rmse=r['metrics']['rmse'], peak=r['metrics']['peak']) for r in d['rows']]
    # Reproduce candidate commands from prefix observations alone. Unseen
    # columns remain NaN, including all physical contact scoring states.
    replay = []
    for mechanism in ('offset','phase'):
        source = OUT/f'bounded_input_T3m_{mechanism}_v1.npz'
        data = np.load(source); m = data['measured']
        policy = BoundedProbe(normalization=mechanism,single_model=False,**CE.FINAL)
        policy.reset(dict(common.P.PLATFORMS['PC125']))
        for j in range(len(m)):
            log = np.full((j+1,16),np.nan)
            log[:,[0,1,2,5,11,12,13]] = m[:j+1]
            log[j,11] = m[j-1,4] if j else m[0,2]
            issued = policy(common.P.Obs(log,j+1))
            np.testing.assert_equal(issued,m[j,4])
        trace = data['coordinate_trace']; index = np.flatnonzero(trace[:,8] > 0)
        assert np.all(np.diff(index) == 1)  # maintained contact in these instances
        cumulative = np.cumsum(trace[index,7])
        levels = trace[index,4]-trace[index[0]-1,4]
        np.testing.assert_allclose(cumulative,levels,atol=1e-16)
        replay.append(dict(mechanism=mechanism, exact_measured_only_prefix_replay=True,
                           bounded_offset_telescoping_max_error=float(np.max(abs(cumulative-levels)))))
    # Check the noise-domain covariance against independent noise realizations
    # on the scalar endpoint design (not the real adaptive closed loop).
    t = np.arange(33)*.008; y=.1*np.sin(2*np.pi*20.8*t)
    B,U=features(t[:-1],y[:-1],y[1:],t[-1],y[-1],t[-1])
    X = residualize(B,U[:,:2]); pseudoinverse=np.linalg.pinv(X)
    rng=np.random.default_rng(71423)
    errors=pseudoinverse@residualize(B,rng.normal(0,.04,(32,10000)))
    predicted=.04**2*np.linalg.inv(X.T@X)
    np.testing.assert_allclose(np.cov(errors),predicted,rtol=.05,atol=.0001)
    result=dict(probe_normalization_algebra=True, constant_gain_equivalence=True,
                nominal_probe_constraints=True, covariance_same_model_reference=True,
                causal_prefix_replay=replay, paired_results=tables,
                no_universally_better_replacement_adopted=True)
    target=OUT/'excitation_mechanism_checkpoint.json'
    target.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
